package sam

import (
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"strconv"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// BAMMagic is the 4-byte signature that introduces the binary header of a
// BAM file: the bytes "BAM" followed by 0x01.
var BAMMagic = []byte{'B', 'A', 'M', 0x01}

// seqLookup maps 4-bit packed nucleotide codes to their character. The codes
// are taken from the SAM/BAM spec, Table "Nibble encoding of nucleotides":
//
//	0 = '=', 1 = A, 2 = C, 3 = M, 4 = G, 5 = R, 6 = S, 7 = V,
//	8 = T, 9 = W, 10 = Y, 11 = H, 12 = K, 13 = D, 14 = B, 15 = N.
var seqLookup = [...]byte{'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}

// seqEncodeTable maps an ASCII nucleotide character (uppercased) to its 4-bit
// BAM code. Used by the BAM writer when packing SEQ.
var seqEncodeTable = func() [256]byte {
	var t [256]byte
	for i := range t {
		t[i] = 0xff
	}
	for i, c := range seqLookup {
		t[c] = byte(i)
		// Map lowercase forms too.
		if c >= 'A' && c <= 'Z' {
			t[c+32] = byte(i)
		}
	}
	return t
}()

// ErrNotBAM indicates the stream did not start with the BAM magic.
var ErrNotBAM = errors.New("sam: input is not a BAM file (missing BAM\\1 magic)")

// BAMReader decodes alignment records from a BGZF-wrapped BAM stream.
type BAMReader struct {
	// src is the BAM byte stream. For the standard BGZF-wrapped input this is
	// a *bgzip.Reader; for already-decompressed input (e.g. one routed
	// through pkg/htsgo/iohelper, which strips the BGZF layer) it is
	// the raw io.Reader directly.
	src   io.Reader
	bgz   *bgzip.Reader // non-nil when src is the BGZF reader
	hdr   *Header
	refs  []Reference // copied from hdr.Refs for fast indexed lookup
	scrat []byte      // reusable buffer for record bodies
	err   error
}

// NewBAMReader constructs a BAMReader that consumes BGZF-encoded BAM bytes
// from r. The header is parsed eagerly so failures surface up front.
func NewBAMReader(r io.Reader) (*BAMReader, error) {
	bgz, err := bgzip.NewReader(r)
	if err != nil {
		return nil, err
	}
	br := &BAMReader{src: bgz, bgz: bgz}
	if err := br.readHeader(); err != nil {
		return nil, err
	}
	return br, nil
}

// newBAMReaderRaw constructs a BAMReader that reads from an already-
// decompressed BAM stream (one that begins directly with the "BAM\1" magic).
// Internal helper used by NewReader when iohelper has already stripped BGZF.
func newBAMReaderRaw(r io.Reader) (*BAMReader, error) {
	br := &BAMReader{src: r}
	if err := br.readHeader(); err != nil {
		return nil, err
	}
	return br, nil
}

// NewBAMBodyReader constructs a BAMReader that decodes records from r
// using the supplied header for reference resolution. r must already be
// positioned at the start of a record (i.e. past the header / @SQ table
// bytes). The reader does not own r — callers are responsible for closing
// the underlying source.
//
// NewBAMBodyReader is the entry point used by region-query seek paths:
// after seeking the BGZF stream to a chunk's compressed offset and
// skipping the in-block uncompressed bytes, the next byte is a record's
// block_size prefix and NewBAMBodyReader can decode from there.
func NewBAMBodyReader(r io.Reader, hdr *Header) *BAMReader {
	refs := make([]Reference, len(hdr.Refs))
	copy(refs, hdr.Refs)
	return &BAMReader{src: r, hdr: hdr, refs: refs}
}

// Header returns the parsed BAM header.
func (br *BAMReader) Header() *Header { return br.hdr }

// readHeader parses the BAM header: magic, text header, then reference table.
func (br *BAMReader) readHeader() error {
	var magic [4]byte
	if _, err := io.ReadFull(br.src, magic[:]); err != nil {
		return err
	}
	if !bytes.Equal(magic[:], BAMMagic) {
		return ErrNotBAM
	}
	var lText int32
	if err := binary.Read(br.src, binary.LittleEndian, &lText); err != nil {
		return err
	}
	if lText < 0 {
		return fmt.Errorf("sam: negative l_text %d", lText)
	}
	text := make([]byte, lText)
	if _, err := io.ReadFull(br.src, text); err != nil {
		return err
	}
	// Strip trailing NULs that htslib likes to pad the header text with.
	text = bytes.TrimRight(text, "\x00")
	hdr, err := ParseHeaderText(string(text))
	if err != nil {
		return err
	}
	var nRef int32
	if err := binary.Read(br.src, binary.LittleEndian, &nRef); err != nil {
		return err
	}
	if nRef < 0 {
		return fmt.Errorf("sam: negative n_ref %d", nRef)
	}
	binRefs := make([]Reference, 0, nRef)
	for i := int32(0); i < nRef; i++ {
		var lName int32
		if err := binary.Read(br.src, binary.LittleEndian, &lName); err != nil {
			return err
		}
		if lName <= 0 {
			return fmt.Errorf("sam: bad l_name %d for ref %d", lName, i)
		}
		nameBuf := make([]byte, lName)
		if _, err := io.ReadFull(br.src, nameBuf); err != nil {
			return err
		}
		// l_name includes a trailing NUL.
		name := string(bytes.TrimRight(nameBuf, "\x00"))
		var lRef int32
		if err := binary.Read(br.src, binary.LittleEndian, &lRef); err != nil {
			return err
		}
		binRefs = append(binRefs, Reference{Name: name, Length: lRef})
	}
	// Prefer text header @SQ entries if they cover every binary ref; otherwise
	// synthesise a header from the binary table to ensure RNames can be
	// resolved during decoding.
	if len(hdr.Refs) == 0 {
		for _, r := range binRefs {
			hl := HeaderLine{
				Tag: "SQ",
				Fields: []HeaderField{
					{Tag: "SN", Value: r.Name},
					{Tag: "LN", Value: strconv.FormatInt(int64(r.Length), 10)},
				},
			}
			hdr.Lines = append(hdr.Lines, hl)
			hdr.Refs = append(hdr.Refs, r)
		}
	}
	br.hdr = hdr
	br.refs = binRefs
	return nil
}

// ParseEncodedBAMHeader parses a Header from the encoded (uncompressed) BAM
// header bytes (magic + l_text + text + n_ref + reference table). It is used by
// raw-block passthrough paths that inflate the header out of band and need the
// parsed @SQ table without consuming the alignment records.
func ParseEncodedBAMHeader(data []byte) (*Header, error) {
	br := &BAMReader{src: bytes.NewReader(data)}
	if err := br.readHeader(); err != nil {
		return nil, err
	}
	return br.hdr, nil
}

// BAMHeaderEncodedLen returns the number of bytes the encoded BAM header
// (magic + l_text + text + n_ref + reference table) occupies at the start of
// data. data must contain the full header; ErrShortHeader-style errors are
// returned as a non-nil error when it does not. It is used by raw-block
// passthrough paths (reheader / cat) to locate the byte boundary between the
// header and the first alignment record after inflating the leading BGZF
// blocks.
func BAMHeaderEncodedLen(data []byte) (int, error) {
	if len(data) < 12 {
		return 0, io.ErrUnexpectedEOF
	}
	if !bytes.Equal(data[:4], BAMMagic) {
		return 0, ErrNotBAM
	}
	off := 4
	lText := int(int32(binary.LittleEndian.Uint32(data[off : off+4])))
	off += 4
	if lText < 0 {
		return 0, fmt.Errorf("sam: negative l_text %d", lText)
	}
	off += lText
	if off+4 > len(data) {
		return 0, io.ErrUnexpectedEOF
	}
	nRef := int(int32(binary.LittleEndian.Uint32(data[off : off+4])))
	off += 4
	if nRef < 0 {
		return 0, fmt.Errorf("sam: negative n_ref %d", nRef)
	}
	for i := 0; i < nRef; i++ {
		if off+4 > len(data) {
			return 0, io.ErrUnexpectedEOF
		}
		lName := int(int32(binary.LittleEndian.Uint32(data[off : off+4])))
		off += 4
		if lName <= 0 {
			return 0, fmt.Errorf("sam: bad l_name %d for ref %d", lName, i)
		}
		off += lName + 4 // name bytes + l_ref
		if off > len(data) {
			return 0, io.ErrUnexpectedEOF
		}
	}
	return off, nil
}

// Read returns the next BAM record, or io.EOF when no more records are
// available.
func (br *BAMReader) Read() (*Record, error) {
	if br.err != nil {
		return nil, br.err
	}
	var blockSize int32
	if err := binary.Read(br.src, binary.LittleEndian, &blockSize); err != nil {
		if err == io.EOF {
			br.err = io.EOF
		}
		return nil, err
	}
	if blockSize < 32 {
		return nil, fmt.Errorf("sam: BAM block too small (%d)", blockSize)
	}
	if cap(br.scrat) < int(blockSize) {
		br.scrat = make([]byte, blockSize)
	} else {
		br.scrat = br.scrat[:blockSize]
	}
	if _, err := io.ReadFull(br.src, br.scrat); err != nil {
		return nil, err
	}
	return br.decodeRecord(br.scrat)
}

// Close releases the underlying BGZF reader, when one is wrapping the input.
// For raw (already-decompressed) BAM streams Close is a no-op — callers own
// the source io.Reader and are responsible for closing it.
func (br *BAMReader) Close() error {
	if br.bgz != nil {
		return br.bgz.Close()
	}
	return nil
}

// VirtualOffset returns the BGZF virtual offset of the next byte that Read
// will consume from the underlying BAM stream. It is only meaningful for
// readers constructed with NewBAMReader (i.e. with a real BGZF layer). For
// raw, already-decompressed streams (newBAMReaderRaw) it always returns 0.
//
// Callers use VirtualOffset to record the start position of each record as
// they iterate the stream — invoke it *before* calling Read to capture the
// current record's offset; the value is the byte just past the previous
// record (or just past the header for the very first record).
func (br *BAMReader) VirtualOffset() uint64 {
	if br.bgz != nil {
		return br.bgz.VirtualOffset()
	}
	return 0
}

// decodeRecord deserialises one BAM record body (everything after block_size).
func (br *BAMReader) decodeRecord(buf []byte) (*Record, error) {
	if len(buf) < 32 {
		return nil, fmt.Errorf("sam: BAM record body too small (%d)", len(buf))
	}
	rec := &Record{}
	refID := int32(binary.LittleEndian.Uint32(buf[0:4]))
	pos := int32(binary.LittleEndian.Uint32(buf[4:8]))
	lReadName := buf[8]
	mapq := buf[9]
	_ = binary.LittleEndian.Uint16(buf[10:12]) // bin — ignored on read
	nCigarOp := binary.LittleEndian.Uint16(buf[12:14])
	flag := binary.LittleEndian.Uint16(buf[14:16])
	lSeq := int32(binary.LittleEndian.Uint32(buf[16:20]))
	nextRefID := int32(binary.LittleEndian.Uint32(buf[20:24]))
	nextPos := int32(binary.LittleEndian.Uint32(buf[24:28]))
	tlen := int32(binary.LittleEndian.Uint32(buf[28:32]))
	off := 32

	// Read name (l_read_name bytes including trailing NUL).
	if off+int(lReadName) > len(buf) {
		return nil, fmt.Errorf("sam: truncated read name")
	}
	nameBytes := buf[off : off+int(lReadName)]
	if len(nameBytes) > 0 && nameBytes[len(nameBytes)-1] == 0 {
		nameBytes = nameBytes[:len(nameBytes)-1]
	}
	rec.QName = string(nameBytes)
	off += int(lReadName)

	// CIGAR ops: nCigarOp uint32s.
	if off+int(nCigarOp)*4 > len(buf) {
		return nil, fmt.Errorf("sam: truncated CIGAR")
	}
	if nCigarOp > 0 {
		rec.Cigar = make(Cigar, nCigarOp)
		for i := 0; i < int(nCigarOp); i++ {
			rec.Cigar[i] = CigarOp(binary.LittleEndian.Uint32(buf[off : off+4]))
			off += 4
		}
	}

	// Packed SEQ: (l_seq+1)/2 bytes; high nibble first.
	seqLen := int((lSeq + 1) / 2)
	if off+seqLen > len(buf) {
		return nil, fmt.Errorf("sam: truncated SEQ")
	}
	if lSeq > 0 {
		seqOut := make([]byte, lSeq)
		for i := int32(0); i < lSeq; i++ {
			b := buf[off+int(i/2)]
			var nibble byte
			if i%2 == 0 {
				nibble = b >> 4
			} else {
				nibble = b & 0x0f
			}
			seqOut[i] = seqLookup[nibble]
		}
		rec.Seq = string(seqOut)
	}
	off += seqLen

	// QUAL: lSeq bytes of Phred.
	if off+int(lSeq) > len(buf) {
		return nil, fmt.Errorf("sam: truncated QUAL")
	}
	if lSeq > 0 {
		qual := make([]byte, lSeq)
		copy(qual, buf[off:off+int(lSeq)])
		rec.Qual = qual
	}
	off += int(lSeq)

	// AUX: parse remaining bytes as a stream of tag/type/value triples.
	if off < len(buf) {
		aux, err := decodeBAMAux(buf[off:])
		if err != nil {
			return nil, err
		}
		rec.Aux = aux
	}

	// Reference resolution.
	if refID >= 0 && int(refID) < len(br.refs) {
		rec.RName = br.refs[refID].Name
	}
	if nextRefID >= 0 && int(nextRefID) < len(br.refs) {
		if nextRefID == refID {
			rec.RNext = "="
		} else {
			rec.RNext = br.refs[nextRefID].Name
		}
	}

	rec.Flag = flag
	rec.MapQ = mapq
	// BAM POS is 0-based; SAM POS is 1-based. -1 → 0.
	if pos >= 0 {
		rec.Pos = pos + 1
	} else {
		rec.Pos = 0
	}
	if nextPos >= 0 {
		rec.PNext = nextPos + 1
	} else {
		rec.PNext = 0
	}
	rec.TLen = tlen
	return rec, nil
}

// decodeBAMAux walks the binary aux stream and returns parsed Aux entries.
func decodeBAMAux(buf []byte) ([]Aux, error) {
	var out []Aux
	for len(buf) > 0 {
		if len(buf) < 3 {
			return nil, fmt.Errorf("sam: truncated aux header")
		}
		tag := string(buf[:2])
		typ := buf[2]
		buf = buf[3:]
		a := Aux{Tag: tag, Type: typ}
		switch typ {
		case 'A':
			if len(buf) < 1 {
				return nil, fmt.Errorf("sam: truncated aux 'A'")
			}
			a.Value = string(buf[:1])
			buf = buf[1:]
		case 'c':
			if len(buf) < 1 {
				return nil, fmt.Errorf("sam: truncated aux 'c'")
			}
			a.Value = int64(int8(buf[0]))
			buf = buf[1:]
		case 'C':
			if len(buf) < 1 {
				return nil, fmt.Errorf("sam: truncated aux 'C'")
			}
			a.Value = int64(buf[0])
			buf = buf[1:]
		case 's':
			if len(buf) < 2 {
				return nil, fmt.Errorf("sam: truncated aux 's'")
			}
			a.Value = int64(int16(binary.LittleEndian.Uint16(buf[:2])))
			buf = buf[2:]
		case 'S':
			if len(buf) < 2 {
				return nil, fmt.Errorf("sam: truncated aux 'S'")
			}
			a.Value = int64(binary.LittleEndian.Uint16(buf[:2]))
			buf = buf[2:]
		case 'i':
			if len(buf) < 4 {
				return nil, fmt.Errorf("sam: truncated aux 'i'")
			}
			a.Value = int64(int32(binary.LittleEndian.Uint32(buf[:4])))
			buf = buf[4:]
		case 'I':
			if len(buf) < 4 {
				return nil, fmt.Errorf("sam: truncated aux 'I'")
			}
			a.Value = int64(binary.LittleEndian.Uint32(buf[:4]))
			buf = buf[4:]
		case 'f':
			if len(buf) < 4 {
				return nil, fmt.Errorf("sam: truncated aux 'f'")
			}
			a.Value = float64(math.Float32frombits(binary.LittleEndian.Uint32(buf[:4])))
			buf = buf[4:]
		case 'Z', 'H':
			end := bytes.IndexByte(buf, 0)
			if end < 0 {
				return nil, fmt.Errorf("sam: unterminated aux 'Z'/'H'")
			}
			a.Value = string(buf[:end])
			buf = buf[end+1:]
		case 'B':
			if len(buf) < 5 {
				return nil, fmt.Errorf("sam: truncated aux 'B' header")
			}
			sub := buf[0]
			count := binary.LittleEndian.Uint32(buf[1:5])
			buf = buf[5:]
			a.ArrayType = sub
			var elemSize int
			switch sub {
			case 'c', 'C':
				elemSize = 1
			case 's', 'S':
				elemSize = 2
			case 'i', 'I', 'f':
				elemSize = 4
			default:
				return nil, fmt.Errorf("sam: unknown aux 'B' subtype %q", sub)
			}
			need := int(count) * elemSize
			if len(buf) < need {
				return nil, fmt.Errorf("sam: truncated aux 'B' body")
			}
			for j := uint32(0); j < count; j++ {
				off := int(j) * elemSize
				switch sub {
				case 'c':
					a.ArrayValues = append(a.ArrayValues, int64(int8(buf[off])))
				case 'C':
					a.ArrayValues = append(a.ArrayValues, int64(buf[off]))
				case 's':
					a.ArrayValues = append(a.ArrayValues, int64(int16(binary.LittleEndian.Uint16(buf[off:off+2]))))
				case 'S':
					a.ArrayValues = append(a.ArrayValues, int64(binary.LittleEndian.Uint16(buf[off:off+2])))
				case 'i':
					a.ArrayValues = append(a.ArrayValues, int64(int32(binary.LittleEndian.Uint32(buf[off:off+4]))))
				case 'I':
					a.ArrayValues = append(a.ArrayValues, int64(binary.LittleEndian.Uint32(buf[off:off+4])))
				case 'f':
					a.ArrayValues = append(a.ArrayValues, float64(math.Float32frombits(binary.LittleEndian.Uint32(buf[off:off+4]))))
				}
			}
			buf = buf[need:]
		default:
			return nil, fmt.Errorf("sam: unknown aux type %q", typ)
		}
		out = append(out, a)
	}
	return out, nil
}
