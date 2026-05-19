package sam

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"math"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

// BAMWriter emits BAM-encoded records on top of a BGZF stream.
type BAMWriter struct {
	bw     *bgzip.Writer
	hdr    *Header
	refMap map[string]int32
	closed bool
}

// NewBAMWriter constructs a BAMWriter that writes BGZF-compressed BAM bytes
// to w. Close must be called to finalise the BGZF EOF block.
func NewBAMWriter(w io.Writer) *BAMWriter {
	return &BAMWriter{bw: bgzip.NewWriter(w)}
}

// WriteHeader serialises the header to BAM: magic, l_text, header text,
// n_ref, then per-ref name and length.
func (bw *BAMWriter) WriteHeader(h *Header) error {
	if h == nil {
		h = &Header{}
	}
	bw.hdr = h
	bw.refMap = make(map[string]int32, len(h.Refs))
	for i, r := range h.Refs {
		bw.refMap[r.Name] = int32(i)
	}
	var buf bytes.Buffer
	buf.Write(BAMMagic)
	text := h.Text()
	if err := binary.Write(&buf, binary.LittleEndian, int32(len(text))); err != nil {
		return err
	}
	buf.WriteString(text)
	if err := binary.Write(&buf, binary.LittleEndian, int32(len(h.Refs))); err != nil {
		return err
	}
	for _, r := range h.Refs {
		nameBytes := []byte(r.Name)
		nameBytes = append(nameBytes, 0)
		if err := binary.Write(&buf, binary.LittleEndian, int32(len(nameBytes))); err != nil {
			return err
		}
		buf.Write(nameBytes)
		if err := binary.Write(&buf, binary.LittleEndian, r.Length); err != nil {
			return err
		}
	}
	_, err := bw.bw.Write(buf.Bytes())
	return err
}

// Write serialises one record into the BAM stream.
func (bw *BAMWriter) Write(rec *Record) error {
	body, err := bw.encodeRecord(rec)
	if err != nil {
		return err
	}
	var sz [4]byte
	binary.LittleEndian.PutUint32(sz[:], uint32(len(body)))
	if _, err := bw.bw.Write(sz[:]); err != nil {
		return err
	}
	_, err = bw.bw.Write(body)
	return err
}

// Close flushes the BGZF stream and emits the EOF block.
func (bw *BAMWriter) Close() error {
	if bw.closed {
		return nil
	}
	bw.closed = true
	return bw.bw.Close()
}

// encodeRecord serialises one record into the BAM record body format.
func (bw *BAMWriter) encodeRecord(rec *Record) ([]byte, error) {
	var buf bytes.Buffer
	refID := int32(-1)
	if rec.RName != "" && rec.RName != "*" {
		if id, ok := bw.refMap[rec.RName]; ok {
			refID = id
		} else {
			return nil, fmt.Errorf("sam: BAM writer: unknown RNAME %q (not in header)", rec.RName)
		}
	}
	nextRefID := int32(-1)
	switch rec.RNext {
	case "", "*":
		// stays -1
	case "=":
		nextRefID = refID
	default:
		if id, ok := bw.refMap[rec.RNext]; ok {
			nextRefID = id
		} else {
			return nil, fmt.Errorf("sam: BAM writer: unknown RNEXT %q (not in header)", rec.RNext)
		}
	}

	// SAM POS is 1-based, BAM 0-based; 0 → -1.
	var bamPos int32 = -1
	if rec.Pos > 0 {
		bamPos = rec.Pos - 1
	}
	var bamPNext int32 = -1
	if rec.PNext > 0 {
		bamPNext = rec.PNext - 1
	}

	nameBytes := []byte(rec.QName)
	nameBytes = append(nameBytes, 0)
	if len(nameBytes) > 255 {
		return nil, fmt.Errorf("sam: BAM read name too long (%d > 255)", len(nameBytes))
	}

	lSeq := int32(len(rec.Seq))
	bin := reg2bin(int(bamPos), int(bamPos)+rec.Cigar.ReferenceLength())

	// Fixed 32-byte header.
	binary.Write(&buf, binary.LittleEndian, refID)
	binary.Write(&buf, binary.LittleEndian, bamPos)
	buf.WriteByte(byte(len(nameBytes)))
	buf.WriteByte(rec.MapQ)
	binary.Write(&buf, binary.LittleEndian, uint16(bin))
	binary.Write(&buf, binary.LittleEndian, uint16(len(rec.Cigar)))
	binary.Write(&buf, binary.LittleEndian, rec.Flag)
	binary.Write(&buf, binary.LittleEndian, lSeq)
	binary.Write(&buf, binary.LittleEndian, nextRefID)
	binary.Write(&buf, binary.LittleEndian, bamPNext)
	binary.Write(&buf, binary.LittleEndian, rec.TLen)

	// Read name (with trailing NUL).
	buf.Write(nameBytes)

	// CIGAR.
	for _, op := range rec.Cigar {
		binary.Write(&buf, binary.LittleEndian, uint32(op))
	}

	// Packed SEQ.
	seqEncoded := encodeSeq(rec.Seq)
	buf.Write(seqEncoded)

	// QUAL.
	if int32(len(rec.Qual)) == lSeq && lSeq > 0 {
		buf.Write(rec.Qual)
	} else {
		// Either no quality or mismatched length: emit lSeq bytes of 0xff.
		pad := make([]byte, lSeq)
		for i := range pad {
			pad[i] = 0xff
		}
		buf.Write(pad)
	}

	// AUX.
	for _, a := range rec.Aux {
		if err := encodeBAMAux(&buf, a); err != nil {
			return nil, err
		}
	}
	return buf.Bytes(), nil
}

// encodeSeq packs an ASCII nucleotide string into the BAM 4-bit form.
func encodeSeq(s string) []byte {
	n := (len(s) + 1) / 2
	out := make([]byte, n)
	for i := 0; i < len(s); i++ {
		code := seqEncodeTable[s[i]]
		if code == 0xff {
			code = 15 // unknown → N
		}
		if i%2 == 0 {
			out[i/2] = code << 4
		} else {
			out[i/2] |= code
		}
	}
	return out
}

// encodeBAMAux serialises one aux field to BAM binary form.
func encodeBAMAux(buf *bytes.Buffer, a Aux) error {
	if len(a.Tag) != 2 {
		return fmt.Errorf("sam: aux tag must be 2 chars, got %q", a.Tag)
	}
	buf.WriteString(a.Tag)
	switch a.Type {
	case 'A':
		buf.WriteByte('A')
		s, _ := a.Value.(string)
		if len(s) >= 1 {
			buf.WriteByte(s[0])
		} else {
			buf.WriteByte(0)
		}
	case 'c':
		buf.WriteByte('c')
		v, _ := a.Value.(int64)
		buf.WriteByte(byte(int8(v)))
	case 'C':
		buf.WriteByte('C')
		v, _ := a.Value.(int64)
		buf.WriteByte(byte(v))
	case 's':
		buf.WriteByte('s')
		v, _ := a.Value.(int64)
		var b [2]byte
		binary.LittleEndian.PutUint16(b[:], uint16(int16(v)))
		buf.Write(b[:])
	case 'S':
		buf.WriteByte('S')
		v, _ := a.Value.(int64)
		var b [2]byte
		binary.LittleEndian.PutUint16(b[:], uint16(v))
		buf.Write(b[:])
	case 'i':
		v, _ := a.Value.(int64)
		// Choose the most compact integer encoding when writing from a
		// generic 'i' aux to keep BAM files small.
		writeBAMIntCompact(buf, v)
	case 'I':
		buf.WriteByte('I')
		v, _ := a.Value.(int64)
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], uint32(v))
		buf.Write(b[:])
	case 'f':
		buf.WriteByte('f')
		v, _ := a.Value.(float64)
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], math.Float32bits(float32(v)))
		buf.Write(b[:])
	case 'Z':
		buf.WriteByte('Z')
		s, _ := a.Value.(string)
		buf.WriteString(s)
		buf.WriteByte(0)
	case 'H':
		buf.WriteByte('H')
		s, _ := a.Value.(string)
		buf.WriteString(s)
		buf.WriteByte(0)
	case 'B':
		buf.WriteByte('B')
		buf.WriteByte(a.ArrayType)
		var sz [4]byte
		binary.LittleEndian.PutUint32(sz[:], uint32(len(a.ArrayValues)))
		buf.Write(sz[:])
		for _, v := range a.ArrayValues {
			switch a.ArrayType {
			case 'c':
				n, _ := v.(int64)
				buf.WriteByte(byte(int8(n)))
			case 'C':
				n, _ := v.(int64)
				buf.WriteByte(byte(n))
			case 's':
				n, _ := v.(int64)
				var b [2]byte
				binary.LittleEndian.PutUint16(b[:], uint16(int16(n)))
				buf.Write(b[:])
			case 'S':
				n, _ := v.(int64)
				var b [2]byte
				binary.LittleEndian.PutUint16(b[:], uint16(n))
				buf.Write(b[:])
			case 'i':
				n, _ := v.(int64)
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], uint32(int32(n)))
				buf.Write(b[:])
			case 'I':
				n, _ := v.(int64)
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], uint32(n))
				buf.Write(b[:])
			case 'f':
				f, _ := v.(float64)
				var b [4]byte
				binary.LittleEndian.PutUint32(b[:], math.Float32bits(float32(f)))
				buf.Write(b[:])
			default:
				return fmt.Errorf("sam: unsupported B array subtype %q", a.ArrayType)
			}
		}
	default:
		return fmt.Errorf("sam: unknown aux type %q", a.Type)
	}
	return nil
}

// writeBAMIntCompact writes an integer aux value using the smallest BAM
// integer type that can hold it.
func writeBAMIntCompact(buf *bytes.Buffer, v int64) {
	switch {
	case v >= 0 && v <= 0xff:
		buf.WriteByte('C')
		buf.WriteByte(byte(v))
	case v >= -128 && v < 0:
		buf.WriteByte('c')
		buf.WriteByte(byte(int8(v)))
	case v >= 0 && v <= 0xffff:
		buf.WriteByte('S')
		var b [2]byte
		binary.LittleEndian.PutUint16(b[:], uint16(v))
		buf.Write(b[:])
	case v >= -32768 && v < 0:
		buf.WriteByte('s')
		var b [2]byte
		binary.LittleEndian.PutUint16(b[:], uint16(int16(v)))
		buf.Write(b[:])
	case v >= 0 && v <= 0xffffffff:
		buf.WriteByte('I')
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], uint32(v))
		buf.Write(b[:])
	default:
		buf.WriteByte('i')
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], uint32(int32(v)))
		buf.Write(b[:])
	}
}

// reg2bin computes the UCSC bin number for a 0-based half-open [beg, end)
// interval, delegating to the shared implementation in pkg/htsgo/tabix.
// For BAM records the writer guards against a degenerate empty CIGAR (which
// gives end == beg) by bumping end up so Reg2bin's half-open contract is met.
func reg2bin(beg, end int) int {
	if end <= beg {
		end = beg + 1
	}
	return tabix.Reg2bin(beg, end)
}
