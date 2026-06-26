package sam

import (
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"strconv"
	"unsafe"

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
	// seqScratch holds the expanded SEQ nibbles between decode and the single
	// string conversion, reused across records to avoid a per-record slice.
	seqScratch []byte
	// numScratch is a small reusable buffer for integer/float formatting on the
	// BAM->SAM text fast path (WriteSAMBody), so serialising a record's numeric
	// fields and aux integers allocates nothing. 32 bytes is ample: the longest
	// base-10 int64 is 20 chars and the longest float32 %g is well under 32.
	numScratch [32]byte
	// textScratch is a reusable buffer the BAM->SAM fast path fills with a
	// record's expanded SEQ or ASCII-33 QUAL so each is written to the bufio
	// writer in a single Write rather than one WriteByte per base. Reused
	// across records to stay allocation-free.
	textScratch []byte
	err         error
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

// ReadRaw returns the next BAM record's body bytes — everything after the
// 4-byte block_size prefix — in a freshly allocated, caller-owned slice,
// WITHOUT decoding them into a *Record. It returns io.EOF when no more records
// are available, exactly like Read.
//
// Unlike Read/ReadInto, the returned slice does not alias the reader's reused
// scratch buffer: every call allocates a new []byte sized to the record body,
// so callers may retain or buffer the bytes across further reads (samtools
// sort buffers and spills them). The bytes are the verbatim on-disk BAM record
// body — the fixed 32-byte prefix followed by read name, CIGAR, packed SEQ,
// QUAL and aux — suitable for WriteRaw to emit unchanged.
func (br *BAMReader) ReadRaw() ([]byte, error) {
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
	body := make([]byte, blockSize)
	if _, err := io.ReadFull(br.src, body); err != nil {
		return nil, err
	}
	return body, nil
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

// ReadInto decodes the next record into dst, reusing dst's Cigar and Qual
// backing arrays (and dst itself) to avoid the per-record allocation that
// Read incurs. It is meant for consume-and-discard scans (flagstat, depth,
// stats): the caller must not retain pointers into dst — including its Seq,
// Qual, Cigar or Aux — across calls, since the next ReadInto overwrites them.
// It returns io.EOF at end of stream, like Read.
func (br *BAMReader) ReadInto(dst *Record) error {
	if br.err != nil {
		return br.err
	}
	var blockSize int32
	if err := binary.Read(br.src, binary.LittleEndian, &blockSize); err != nil {
		if err == io.EOF {
			br.err = io.EOF
		}
		return err
	}
	if blockSize < 32 {
		return fmt.Errorf("sam: BAM block too small (%d)", blockSize)
	}
	if cap(br.scrat) < int(blockSize) {
		br.scrat = make([]byte, blockSize)
	} else {
		br.scrat = br.scrat[:blockSize]
	}
	if _, err := io.ReadFull(br.src, br.scrat); err != nil {
		return err
	}
	// owned=false: ReadInto's contract is that the caller must not retain dst
	// past the next call, so QName/Seq may alias the reused buffers.
	return br.decodeInto(dst, br.scrat, false)
}

// ReadShallowInto decodes only the fixed-prefix fields of the next record
// into dst — RName/Pos/MapQ/Flag/RNext/PNext/TLen — and skips the entire
// variable-length region (read name, CIGAR, SEQ, QUAL, aux) without parsing
// it. The variable fields on dst (QName, Cigar, Seq, Qual, Aux) are reset to
// empty/zero so a reused record never carries stale data from a previous
// decode. The stream is advanced past the full record, exactly as ReadInto
// would advance it.
//
// It is meant for counters that touch only flags, mapping quality and the
// mate reference (e.g. samtools flagstat), where decoding the read name and
// expanding the packed SEQ nibbles is pure waste. Like ReadInto, it is
// allocation-free across calls and the caller must not retain dst past the
// next call. It returns io.EOF at end of stream.
func (br *BAMReader) ReadShallowInto(dst *Record) error {
	if br.err != nil {
		return br.err
	}
	var blockSize int32
	if err := binary.Read(br.src, binary.LittleEndian, &blockSize); err != nil {
		if err == io.EOF {
			br.err = io.EOF
		}
		return err
	}
	if blockSize < 32 {
		return fmt.Errorf("sam: BAM block too small (%d)", blockSize)
	}
	if cap(br.scrat) < int(blockSize) {
		br.scrat = make([]byte, blockSize)
	} else {
		br.scrat = br.scrat[:blockSize]
	}
	if _, err := io.ReadFull(br.src, br.scrat); err != nil {
		return err
	}
	return br.decodeShallowInto(dst, br.scrat)
}

// decodeShallowInto deserialises only the 32-byte fixed prefix of a BAM record
// body into dst, leaving the variable-length region unparsed. The fields it
// does not populate are reset so a reused record never leaks stale data.
func (br *BAMReader) decodeShallowInto(dst *Record, buf []byte) error {
	if len(buf) < 32 {
		return fmt.Errorf("sam: BAM record body too small (%d)", len(buf))
	}
	refID := int32(binary.LittleEndian.Uint32(buf[0:4]))
	pos := int32(binary.LittleEndian.Uint32(buf[4:8]))
	mapq := buf[9]
	flag := binary.LittleEndian.Uint16(buf[14:16])
	nextRefID := int32(binary.LittleEndian.Uint32(buf[20:24]))
	nextPos := int32(binary.LittleEndian.Uint32(buf[24:28]))
	tlen := int32(binary.LittleEndian.Uint32(buf[28:32]))

	// Clear the variable-length fields so a reused record never carries stale
	// data from a previous (possibly full) decode.
	dst.QName, dst.Seq = "", ""
	dst.Cigar = dst.Cigar[:0]
	dst.Qual = dst.Qual[:0]
	dst.Aux = dst.Aux[:0]
	dst.auxIndex = nil
	dst.RName, dst.RNext = "", ""

	// Reference resolution (mirrors decodeInto).
	if refID >= 0 && int(refID) < len(br.refs) {
		dst.RName = br.refs[refID].Name
	}
	if nextRefID >= 0 && int(nextRefID) < len(br.refs) {
		if nextRefID == refID {
			dst.RNext = "="
		} else {
			dst.RNext = br.refs[nextRefID].Name
		}
	}

	dst.Flag = flag
	dst.MapQ = mapq
	// BAM POS is 0-based; SAM POS is 1-based. -1 → 0.
	if pos >= 0 {
		dst.Pos = int64(pos) + 1
	} else {
		dst.Pos = 0
	}
	if nextPos >= 0 {
		dst.PNext = int64(nextPos) + 1
	} else {
		dst.PNext = 0
	}
	dst.TLen = int64(tlen)
	return nil
}

// ReadDepthInto decodes only the fields samtools depth consumes — RName, Pos,
// Flag, MapQ and CIGAR — into dst, plus QUAL when needQual is true. It skips
// the read name, the mate reference, the packed SEQ expansion and the entire
// aux stream (the dominant decode cost on real BAMs, none of which depth uses).
// The variable-length fields it does not populate (QName, Seq, RNext, Aux, and
// Qual when needQual is false) are reset to empty so a reused record never
// carries stale data from a previous decode. Like ReadInto it is allocation-
// free across calls and advances the stream past the full record. It returns
// io.EOF at end of stream.
func (br *BAMReader) ReadDepthInto(dst *Record, needQual bool) error {
	if br.err != nil {
		return br.err
	}
	var blockSize int32
	if err := binary.Read(br.src, binary.LittleEndian, &blockSize); err != nil {
		if err == io.EOF {
			br.err = io.EOF
		}
		return err
	}
	if blockSize < 32 {
		return fmt.Errorf("sam: BAM block too small (%d)", blockSize)
	}
	if cap(br.scrat) < int(blockSize) {
		br.scrat = make([]byte, blockSize)
	} else {
		br.scrat = br.scrat[:blockSize]
	}
	if _, err := io.ReadFull(br.src, br.scrat); err != nil {
		return err
	}
	return br.decodeDepthInto(dst, br.scrat, needQual)
}

// decodeDepthInto deserialises the depth-relevant fields of one BAM record
// body into dst. It mirrors decodeInto for RName/Pos/Flag/MapQ/CIGAR (and QUAL
// when needQual) but skips read-name, mate-reference, SEQ and aux parsing.
func (br *BAMReader) decodeDepthInto(dst *Record, buf []byte, needQual bool) error {
	if len(buf) < 32 {
		return fmt.Errorf("sam: BAM record body too small (%d)", len(buf))
	}
	refID := int32(binary.LittleEndian.Uint32(buf[0:4]))
	pos := int32(binary.LittleEndian.Uint32(buf[4:8]))
	lReadName := buf[8]
	mapq := buf[9]
	nCigarOp := binary.LittleEndian.Uint16(buf[12:14])
	flag := binary.LittleEndian.Uint16(buf[14:16])
	lSeq := int32(binary.LittleEndian.Uint32(buf[16:20]))

	// Clear the fields depth never reads so a reused record carries no stale
	// data from an earlier (possibly full) decode.
	dst.QName, dst.Seq, dst.RNext = "", "", ""
	dst.Aux = dst.Aux[:0]
	dst.auxIndex = nil

	off := 32 + int(lReadName)

	// CIGAR ops: nCigarOp uint32s. Reuse dst.Cigar's backing array.
	if off+int(nCigarOp)*4 > len(buf) {
		return fmt.Errorf("sam: truncated CIGAR")
	}
	if cap(dst.Cigar) >= int(nCigarOp) {
		dst.Cigar = dst.Cigar[:nCigarOp]
	} else {
		dst.Cigar = make(Cigar, nCigarOp)
	}
	for i := 0; i < int(nCigarOp); i++ {
		dst.Cigar[i] = CigarOp(binary.LittleEndian.Uint32(buf[off : off+4]))
		off += 4
	}

	if needQual {
		// Skip the packed SEQ, then copy QUAL (lSeq Phred bytes).
		seqLen := int((lSeq + 1) / 2)
		off += seqLen
		if off+int(lSeq) > len(buf) {
			return fmt.Errorf("sam: truncated QUAL")
		}
		if cap(dst.Qual) >= int(lSeq) {
			dst.Qual = dst.Qual[:lSeq]
		} else {
			dst.Qual = make([]byte, lSeq)
		}
		copy(dst.Qual, buf[off:off+int(lSeq)])
	} else {
		dst.Qual = dst.Qual[:0]
	}

	// Reference resolution (RName only; depth never reads the mate reference).
	dst.RName = ""
	if refID >= 0 && int(refID) < len(br.refs) {
		dst.RName = br.refs[refID].Name
	}

	dst.Flag = flag
	dst.MapQ = mapq
	// BAM POS is 0-based; SAM POS is 1-based. -1 → 0.
	if pos >= 0 {
		dst.Pos = int64(pos) + 1
	} else {
		dst.Pos = 0
	}
	return nil
}

// RawAuxOffset returns the byte offset within a BAM record body at which the
// aux stream begins — i.e. just past the fixed 32-byte prefix, read name, CIGAR,
// packed SEQ and QUAL. It is the offset callers pass to FindRawAuxTag's slice
// (body[off:]) to scan a single record's aux fields without decoding the whole
// record. It returns an error when the body is too short to hold the prefix or
// the fixed-length region it describes overruns the body.
func RawAuxOffset(body []byte) (int, error) {
	if len(body) < 32 {
		return 0, fmt.Errorf("sam: BAM record body too small (%d)", len(body))
	}
	lReadName := int(body[8])
	nCigarOp := int(binary.LittleEndian.Uint16(body[12:14]))
	lSeq := int(int32(binary.LittleEndian.Uint32(body[16:20])))
	off := 32 + lReadName + nCigarOp*4 + (lSeq+1)/2 + lSeq
	if off > len(body) {
		return 0, fmt.Errorf("sam: BAM record fixed region overruns body (%d > %d)", off, len(body))
	}
	return off, nil
}

// DecodeRecordBody decodes one verbatim BAM record body (everything after the
// 4-byte block_size, exactly as ReadRaw returns it) into a freshly allocated,
// caller-owned *Record, resolving reference names against br's header. It is the
// inverse of ReadRaw for callers that buffered raw bodies but need a decoded
// record for a non-BAM output (e.g. samtools sort writing SAM text). The decoded
// record is identical to what Read would have produced for the same bytes.
func (br *BAMReader) DecodeRecordBody(body []byte) (*Record, error) {
	return br.decodeRecord(body)
}

// decodeRecord deserialises one BAM record body (everything after block_size)
// into a freshly allocated Record.
func (br *BAMReader) decodeRecord(buf []byte) (*Record, error) {
	rec := &Record{}
	// owned: Read returns this record to the caller, who may retain it, so its
	// QName/Seq strings must own their memory.
	if err := br.decodeInto(rec, buf, true); err != nil {
		return nil, err
	}
	return rec, nil
}

// decodeInto deserialises one BAM record body into rec, reusing rec's Cigar and
// Qual backing arrays where they are large enough. Every field rec carries is
// reset so a reused record never leaks stale data from a previous decode.
// decodeInto decodes one BAM record body into rec. When owned is true the QName
// and Seq strings are copied so the caller may retain rec (the Read path); when
// false they alias the reader's reused buffers (br.scrat / br.seqScratch) with
// no per-record string allocation — valid only until the next read, which is
// exactly ReadInto's no-retain contract. QName and Seq are the two big
// per-record allocations on the BAM hot path, so eliminating them for the
// streaming consumers (samtools view/stats/depth, bcftools) is the memory win.
func (br *BAMReader) decodeInto(rec *Record, buf []byte, owned bool) error {
	if len(buf) < 32 {
		return fmt.Errorf("sam: BAM record body too small (%d)", len(buf))
	}
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

	// Clear the fields that the rest of the decode sets only conditionally,
	// so a reused record never carries stale data from a previous decode.
	rec.RName, rec.RNext, rec.Seq = "", "", ""
	rec.auxIndex = nil

	// Read name (l_read_name bytes including trailing NUL).
	if off+int(lReadName) > len(buf) {
		return fmt.Errorf("sam: truncated read name")
	}
	nameBytes := buf[off : off+int(lReadName)]
	if len(nameBytes) > 0 && nameBytes[len(nameBytes)-1] == 0 {
		nameBytes = nameBytes[:len(nameBytes)-1]
	}
	if owned || len(nameBytes) == 0 {
		rec.QName = string(nameBytes)
	} else {
		// Alias the block buffer (reused on the next read) — no allocation.
		rec.QName = unsafe.String(&nameBytes[0], len(nameBytes))
	}
	off += int(lReadName)

	// CIGAR ops: nCigarOp uint32s. Reuse rec.Cigar's backing array.
	if off+int(nCigarOp)*4 > len(buf) {
		return fmt.Errorf("sam: truncated CIGAR")
	}
	if cap(rec.Cigar) >= int(nCigarOp) {
		rec.Cigar = rec.Cigar[:nCigarOp]
	} else {
		rec.Cigar = make(Cigar, nCigarOp)
	}
	for i := 0; i < int(nCigarOp); i++ {
		rec.Cigar[i] = CigarOp(binary.LittleEndian.Uint32(buf[off : off+4]))
		off += 4
	}

	// Packed SEQ: (l_seq+1)/2 bytes; high nibble first. The nibbles are
	// expanded into a reused scratch buffer and converted to the Seq string
	// once, avoiding a separate per-record intermediate slice.
	seqLen := int((lSeq + 1) / 2)
	if off+seqLen > len(buf) {
		return fmt.Errorf("sam: truncated SEQ")
	}
	if lSeq > 0 {
		if cap(br.seqScratch) < int(lSeq) {
			br.seqScratch = make([]byte, lSeq)
		} else {
			br.seqScratch = br.seqScratch[:lSeq]
		}
		for i := int32(0); i < lSeq; i++ {
			b := buf[off+int(i/2)]
			var nibble byte
			if i%2 == 0 {
				nibble = b >> 4
			} else {
				nibble = b & 0x0f
			}
			br.seqScratch[i] = seqLookup[nibble]
		}
		if owned {
			rec.Seq = string(br.seqScratch)
		} else {
			// Alias seqScratch (reused on the next read) — no allocation.
			rec.Seq = unsafe.String(&br.seqScratch[0], lSeq)
		}
	}
	off += seqLen

	// QUAL: lSeq bytes of Phred. Reuse rec.Qual's backing array.
	if off+int(lSeq) > len(buf) {
		return fmt.Errorf("sam: truncated QUAL")
	}
	if cap(rec.Qual) >= int(lSeq) {
		rec.Qual = rec.Qual[:lSeq]
	} else {
		rec.Qual = make([]byte, lSeq)
	}
	copy(rec.Qual, buf[off:off+int(lSeq)])
	off += int(lSeq)

	// AUX: parse remaining bytes as a stream of tag/type/value triples,
	// reusing rec.Aux's backing array.
	rec.Aux = rec.Aux[:0]
	if off < len(buf) {
		aux, err := decodeBAMAuxInto(rec.Aux, buf[off:])
		if err != nil {
			return err
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
		rec.Pos = int64(pos) + 1
	} else {
		rec.Pos = 0
	}
	if nextPos >= 0 {
		rec.PNext = int64(nextPos) + 1
	} else {
		rec.PNext = 0
	}
	rec.TLen = int64(tlen)
	return nil
}

// decodeBAMAux walks the binary aux stream and returns parsed Aux entries.
func decodeBAMAux(buf []byte) ([]Aux, error) {
	return decodeBAMAuxInto(nil, buf)
}

// DecodeBAMAux parses a raw on-disk BAM aux byte block (the form carried by
// Record.RawAux) into a fresh []Aux. It is the inverse of serialising a []Aux
// through AppendBAMAux: re-encoding the result reproduces buf byte-for-byte
// (modulo the canonical compact-integer choice the writer would also make). It
// is the exported entry point for callers that hold a raw aux block and need the
// decoded fields.
func DecodeBAMAux(buf []byte) ([]Aux, error) {
	return decodeBAMAuxInto(nil, buf)
}

// rawAuxStep returns the total byte length of one aux entry at the start of buf
// (the 3-byte tag+type header plus the value), so a caller can walk a raw BAM
// aux block tag by tag without decoding any values. It returns an error when the
// entry is truncated or its type byte is unknown.
func rawAuxStep(buf []byte) (int, error) {
	if len(buf) < 3 {
		return 0, fmt.Errorf("sam: truncated aux header")
	}
	typ := buf[2]
	rest := buf[3:]
	switch typ {
	case 'A', 'c', 'C':
		if len(rest) < 1 {
			return 0, fmt.Errorf("sam: truncated aux value")
		}
		return 4, nil
	case 's', 'S':
		if len(rest) < 2 {
			return 0, fmt.Errorf("sam: truncated aux value")
		}
		return 5, nil
	case 'i', 'I', 'f':
		if len(rest) < 4 {
			return 0, fmt.Errorf("sam: truncated aux value")
		}
		return 7, nil
	case 'Z', 'H':
		end := bytes.IndexByte(rest, 0)
		if end < 0 {
			return 0, fmt.Errorf("sam: unterminated aux 'Z'/'H'")
		}
		return 3 + end + 1, nil
	case 'B':
		if len(rest) < 5 {
			return 0, fmt.Errorf("sam: truncated aux 'B' header")
		}
		sub := rest[0]
		count := int(binary.LittleEndian.Uint32(rest[1:5]))
		var elemSize int
		switch sub {
		case 'c', 'C':
			elemSize = 1
		case 's', 'S':
			elemSize = 2
		case 'i', 'I', 'f':
			elemSize = 4
		default:
			return 0, fmt.Errorf("sam: unknown aux 'B' subtype %q", sub)
		}
		need := count * elemSize
		if len(rest) < 5+need {
			return 0, fmt.Errorf("sam: truncated aux 'B' body")
		}
		return 3 + 5 + need, nil
	default:
		return 0, fmt.Errorf("sam: unknown aux type %q", typ)
	}
}

// WalkRawAux scans a raw on-disk BAM aux byte block without decoding any value.
// It reports whether the block carries an MD tag and an NM tag, and the byte
// offset at which a trailing data-series RG tag begins (or -1 when the final tag
// is not RG). It is used by the CRAM→BAM raw-aux MD/NM regeneration to decide
// what to splice and where, mirroring the eager path's struct scan
// (regenerateMDNM) and trailing-RG splice (insertBeforeTrailingRG) but operating
// on raw bytes so no []Aux is materialised. It returns an error only when the
// block is malformed.
func WalkRawAux(buf []byte) (hasMD, hasNM bool, trailingRGStart int, err error) {
	trailingRGStart = -1
	off := 0
	for off < len(buf) {
		if len(buf)-off < 3 {
			return false, false, -1, fmt.Errorf("sam: truncated aux header")
		}
		t0, t1 := buf[off], buf[off+1]
		step, serr := rawAuxStep(buf[off:])
		if serr != nil {
			return false, false, -1, serr
		}
		isRG := t0 == 'R' && t1 == 'G'
		if t0 == 'M' && t1 == 'D' {
			hasMD = true
		}
		if t0 == 'N' && t1 == 'M' {
			hasNM = true
		}
		// Record the start of this tag if it is RG; the loop overwrites it for
		// every RG, but since the data-series RG is always the final tag the
		// last write is the trailing-RG offset. Reset to -1 when a non-RG tag
		// follows so only a genuinely trailing RG is reported.
		if isRG {
			trailingRGStart = off
		} else {
			trailingRGStart = -1
		}
		off += step
	}
	return hasMD, hasNM, trailingRGStart, nil
}

// FindRawAuxTag walks the binary aux stream in buf and decodes ONLY the entry
// whose 2-char tag equals tag, returning it and true. It stops at the first
// match without decoding the rest of the stream, so callers that compare a
// single tag (e.g. samtools sort -t) need not materialise every aux field. The
// decoded Aux — its Type byte, boxed Value, ArrayType and ArrayValues — is
// byte-for-byte identical to what a full decodeBAMAuxInto would produce for the
// same tag. It returns (Aux{}, false) when the tag is absent and an error only
// when the aux stream is malformed before the tag is found.
func FindRawAuxTag(buf []byte, tag string) (Aux, bool, error) {
	if len(tag) != 2 {
		return Aux{}, false, fmt.Errorf("sam: aux tag must be 2 chars, got %q", tag)
	}
	for len(buf) > 0 {
		if len(buf) < 3 {
			return Aux{}, false, fmt.Errorf("sam: truncated aux header")
		}
		t0, t1 := buf[0], buf[1]
		typ := buf[2]
		buf = buf[3:]
		match := t0 == tag[0] && t1 == tag[1]
		var a Aux
		if match {
			a = Aux{Tag: internTag(t0, t1), Type: typ}
		}
		switch typ {
		case 'A':
			if len(buf) < 1 {
				return Aux{}, false, fmt.Errorf("sam: truncated aux 'A'")
			}
			if match {
				a.Value = string(buf[:1])
			}
			buf = buf[1:]
		case 'c':
			if len(buf) < 1 {
				return Aux{}, false, fmt.Errorf("sam: truncated aux 'c'")
			}
			if match {
				a.Value = int64(int8(buf[0]))
			}
			buf = buf[1:]
		case 'C':
			if len(buf) < 1 {
				return Aux{}, false, fmt.Errorf("sam: truncated aux 'C'")
			}
			if match {
				a.Value = int64(buf[0])
			}
			buf = buf[1:]
		case 's':
			if len(buf) < 2 {
				return Aux{}, false, fmt.Errorf("sam: truncated aux 's'")
			}
			if match {
				a.Value = int64(int16(binary.LittleEndian.Uint16(buf[:2])))
			}
			buf = buf[2:]
		case 'S':
			if len(buf) < 2 {
				return Aux{}, false, fmt.Errorf("sam: truncated aux 'S'")
			}
			if match {
				a.Value = int64(binary.LittleEndian.Uint16(buf[:2]))
			}
			buf = buf[2:]
		case 'i':
			if len(buf) < 4 {
				return Aux{}, false, fmt.Errorf("sam: truncated aux 'i'")
			}
			if match {
				a.Value = int64(int32(binary.LittleEndian.Uint32(buf[:4])))
			}
			buf = buf[4:]
		case 'I':
			if len(buf) < 4 {
				return Aux{}, false, fmt.Errorf("sam: truncated aux 'I'")
			}
			if match {
				a.Value = int64(binary.LittleEndian.Uint32(buf[:4]))
			}
			buf = buf[4:]
		case 'f':
			if len(buf) < 4 {
				return Aux{}, false, fmt.Errorf("sam: truncated aux 'f'")
			}
			if match {
				a.Value = float64(math.Float32frombits(binary.LittleEndian.Uint32(buf[:4])))
			}
			buf = buf[4:]
		case 'Z', 'H':
			end := bytes.IndexByte(buf, 0)
			if end < 0 {
				return Aux{}, false, fmt.Errorf("sam: unterminated aux 'Z'/'H'")
			}
			if match {
				a.Value = string(buf[:end])
			}
			buf = buf[end+1:]
		case 'B':
			if len(buf) < 5 {
				return Aux{}, false, fmt.Errorf("sam: truncated aux 'B' header")
			}
			sub := buf[0]
			count := binary.LittleEndian.Uint32(buf[1:5])
			buf = buf[5:]
			var elemSize int
			switch sub {
			case 'c', 'C':
				elemSize = 1
			case 's', 'S':
				elemSize = 2
			case 'i', 'I', 'f':
				elemSize = 4
			default:
				return Aux{}, false, fmt.Errorf("sam: unknown aux 'B' subtype %q", sub)
			}
			need := int(count) * elemSize
			if len(buf) < need {
				return Aux{}, false, fmt.Errorf("sam: truncated aux 'B' body")
			}
			if match {
				a.ArrayType = sub
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
			}
			buf = buf[need:]
		default:
			return Aux{}, false, fmt.Errorf("sam: unknown aux type %q", typ)
		}
		if match {
			return a, true, nil
		}
	}
	return Aux{}, false, nil
}

// tagInternTable holds a pre-allocated string for every 2-char aux tag built
// from printable ASCII bytes (0x20..0x7e). It is built once at package init and
// only read thereafter, so it is safe to share across concurrent BAMReaders. The
// strings it returns are byte-for-byte identical to string(buf[:2]); interning
// them just avoids a fresh 2-byte allocation on the per-record decode hot path.
var tagInternTable = func() *[95][95]string {
	t := new([95][95]string)
	for i := 0; i < 95; i++ {
		for j := 0; j < 95; j++ {
			t[i][j] = string([]byte{byte(0x20 + i), byte(0x20 + j)})
		}
	}
	return t
}()

// internTag returns the interned 2-char tag string for bytes a,b. For the
// printable-ASCII tags that every real BAM uses it returns a shared, pre-built
// string (no allocation); for any out-of-range byte it falls back to a fresh
// string so the decoded tag is always exactly string([]byte{a, b}).
func internTag(a, b byte) string {
	if a >= 0x20 && a <= 0x7e && b >= 0x20 && b <= 0x7e {
		return tagInternTable[a-0x20][b-0x20]
	}
	return string([]byte{a, b})
}

// decodeBAMAuxInto is decodeBAMAux that appends into dst (use dst[:0] to reuse
// a record's existing Aux backing array), returning the grown slice.
//
// FIX 2: it reuses both the dst Aux backing array (callers pass rec.Aux[:0]) and
// each reused entry's ArrayValues backing slice — for the common case of
// re-decoding records into the same *Record, a record with B-array aux no longer
// allocates a fresh []interface{} per read. It also interns the 2-char tag via a
// fixed lookup (internTag) instead of allocating a new 2-byte string each call.
// None of this changes a single decoded value: the Aux struct, its Value boxing,
// and the parsed contents are identical to the allocate-fresh path.
func decodeBAMAuxInto(dst []Aux, buf []byte) ([]Aux, error) {
	// fullCap is dst's backing array length: entries in dst[len:cap] still hold
	// the previous decode's ArrayValues backing slices, which we can reuse.
	out := dst[:0]
	fullCap := dst[:cap(dst)]
	for len(buf) > 0 {
		if len(buf) < 3 {
			return nil, fmt.Errorf("sam: truncated aux header")
		}
		tag := internTag(buf[0], buf[1])
		typ := buf[2]
		buf = buf[3:]
		// Reuse the slot at index len(out) if it exists in the backing array, so
		// its ArrayValues backing slice can be recycled; otherwise start fresh.
		var a Aux
		idx := len(out)
		var reuseArr []interface{}
		if idx < len(fullCap) {
			reuseArr = fullCap[idx].ArrayValues[:0]
		}
		a = Aux{Tag: tag, Type: typ}
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
			// Seed the array with the reused backing slice (len 0, prior cap) so
			// the appends below recycle it instead of allocating from nil.
			a.ArrayValues = reuseArr
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
