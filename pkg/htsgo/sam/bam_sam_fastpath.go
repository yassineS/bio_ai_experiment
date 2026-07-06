package sam

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"strconv"
)

// FastFields holds the small set of decoded fields the samtools view filter
// pipeline consults before deciding whether to emit a record: the flag, mapping
// quality, reference name, 1-based position and the reference-span length of
// the CIGAR. It is populated by BAMReader.ReadSAMInto without building a full
// Record or expanding SEQ/QUAL/aux, so the BAM->SAM text fast path can apply
// flag/MAPQ/region filters cheaply and only serialise records that survive.
//
// RawBody points at the BAMReader's internal scratch buffer holding the current
// record's body (the bytes after block_size). It is only valid until the next
// read on the same reader and must not be retained.
type FastFields struct {
	Flag    uint16
	MapQ    uint8
	RName   string // resolved reference name ("" when refID < 0, formatted as "*")
	Pos     int64  // 1-based; 0 when unmapped
	RefID   int32  // raw refID (for region resolution by caller)
	RefSpan int    // reference bases consumed by the CIGAR (for region overlap)
	RawBody []byte // the record body bytes (after block_size); transient
}

// ReadSAMInto reads the next BAM record, decodes only the cheap fixed-prefix
// fields needed for filtering into ff, and leaves the record body in the
// reader's scratch buffer (exposed via ff.RawBody) for a later WriteSAMBody.
// It returns io.EOF at end of stream.
//
// This is the entry point for the samtools view BAM->SAM text fast path: the
// caller inspects ff to decide whether to keep the record, then — only for
// kept records — calls WriteSAMBody to serialise the body straight from the
// raw bytes, skipping the per-field/per-aux string allocations that a full
// Record decode incurs. ff is valid only until the next read on this reader.
func (br *BAMReader) ReadSAMInto(ff *FastFields) error {
	if br.err != nil {
		return br.err
	}
	// Read the 4-byte little-endian block_size prefix directly rather than via
	// reflection-based binary.Read — this loop runs once per record, so the
	// reflection cost is material at hundreds of thousands of records. Use the
	// reader's reusable sizeBuf field (not a stack-local array) so the slice
	// handed to io.ReadFull does not escape to a fresh heap allocation on every
	// record — that per-record escape was the dominant fast-path allocation.
	if _, err := io.ReadFull(br.src, br.sizeBuf[:]); err != nil {
		if err == io.EOF {
			br.err = io.EOF
		}
		return err
	}
	blockSize := int32(binary.LittleEndian.Uint32(br.sizeBuf[:]))
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
	buf := br.scrat
	if len(buf) < 32 {
		return fmt.Errorf("sam: BAM record body too small (%d)", len(buf))
	}
	refID := int32(binary.LittleEndian.Uint32(buf[0:4]))
	pos := int32(binary.LittleEndian.Uint32(buf[4:8]))
	lReadName := buf[8]
	mapq := buf[9]
	nCigarOp := int(binary.LittleEndian.Uint16(buf[12:14]))
	flag := binary.LittleEndian.Uint16(buf[14:16])

	ff.Flag = flag
	ff.MapQ = mapq
	ff.RefID = refID
	if refID >= 0 && int(refID) < len(br.refs) {
		ff.RName = br.refs[refID].Name
	} else {
		ff.RName = ""
	}
	if pos >= 0 {
		ff.Pos = int64(pos) + 1
	} else {
		ff.Pos = 0
	}
	// Reference span of the CIGAR (used by region-overlap filtering). The CIGAR
	// ops live right after the read name. We sum only ops that consume the
	// reference, reading the packed uint32s directly without materialising a
	// Cigar slice.
	span := 0
	cigOff := 32 + int(lReadName)
	if cigOff+nCigarOp*4 <= len(buf) {
		for i := 0; i < nCigarOp; i++ {
			op := binary.LittleEndian.Uint32(buf[cigOff : cigOff+4])
			cigOff += 4
			if cigarOpConsumesRef(op & 0xf) {
				span += int(op >> 4)
			}
		}
	}
	ff.RefSpan = span
	ff.RawBody = buf
	return nil
}

// WriteSAMBody serialises one BAM record body (the bytes after block_size, as
// captured in FastFields.RawBody) directly to w as a tab-delimited SAM line,
// mirroring htslib's sam_format1. It bypasses the intermediate Record entirely
// — QName, CIGAR, RNEXT, SEQ, QUAL and every aux tag are written straight from
// the raw little-endian fields — so the only allocations are the integer
// formatting scratch reused on the BAMReader.
//
// The output is byte-identical to formatting a fully decoded Record through
// SAMWriter.Write: SEQ "*"/QUAL "*" (all-0xFF) conventions, integer/float aux
// formatting and the '=' RNEXT shorthand all match. It must be called with the
// FastFields produced by the immediately preceding ReadSAMInto on the same
// reader; RawBody is only valid until the next read.
func (br *BAMReader) WriteSAMBody(w *bufio.Writer, ff *FastFields) error {
	buf := ff.RawBody
	if len(buf) < 32 {
		return fmt.Errorf("sam: BAM record body too small (%d)", len(buf))
	}
	refID := int32(binary.LittleEndian.Uint32(buf[0:4]))
	pos := int32(binary.LittleEndian.Uint32(buf[4:8]))
	lReadName := int(buf[8])
	mapq := buf[9]
	nCigarOp := int(binary.LittleEndian.Uint16(buf[12:14]))
	flag := binary.LittleEndian.Uint16(buf[14:16])
	lSeq := int(binary.LittleEndian.Uint32(buf[16:20]))
	nextRefID := int32(binary.LittleEndian.Uint32(buf[20:24]))
	nextPos := int32(binary.LittleEndian.Uint32(buf[24:28]))
	tlen := int32(binary.LittleEndian.Uint32(buf[28:32]))

	off := 32
	// QNAME (read name minus its trailing NUL).
	if off+lReadName > len(buf) {
		return fmt.Errorf("sam: truncated read name")
	}
	name := buf[off : off+lReadName]
	if len(name) > 0 && name[len(name)-1] == 0 {
		name = name[:len(name)-1]
	}
	w.Write(name)
	off += lReadName
	w.WriteByte('\t')

	// FLAG.
	w.Write(br.appendUint(uint64(flag)))
	w.WriteByte('\t')

	// RNAME.
	if refID >= 0 && int(refID) < len(br.refs) {
		w.WriteString(br.refs[refID].Name)
	} else {
		w.WriteByte('*')
	}
	w.WriteByte('\t')

	// POS (1-based; BAM 0-based, -1 -> 0).
	if pos >= 0 {
		w.Write(br.appendUint(uint64(pos) + 1))
	} else {
		w.WriteByte('0')
	}
	w.WriteByte('\t')

	// MAPQ.
	w.Write(br.appendUint(uint64(mapq)))
	w.WriteByte('\t')

	// CIGAR.
	if off+nCigarOp*4 > len(buf) {
		return fmt.Errorf("sam: truncated CIGAR")
	}
	if nCigarOp == 0 {
		w.WriteByte('*')
	} else {
		for i := 0; i < nCigarOp; i++ {
			op := binary.LittleEndian.Uint32(buf[off : off+4])
			off += 4
			w.Write(br.appendUint(uint64(op >> 4)))
			oc := op & 0xf
			if int(oc) < len(cigarOpChars) {
				w.WriteByte(cigarOpChars[oc])
			} else {
				w.WriteByte('?')
			}
		}
	}
	w.WriteByte('\t')

	// RNEXT ('=' when the mate ref equals this record's ref).
	if nextRefID >= 0 && int(nextRefID) < len(br.refs) {
		if nextRefID == refID {
			w.WriteByte('=')
		} else {
			w.WriteString(br.refs[nextRefID].Name)
		}
	} else {
		w.WriteByte('*')
	}
	w.WriteByte('\t')

	// PNEXT.
	if nextPos >= 0 {
		w.Write(br.appendUint(uint64(nextPos) + 1))
	} else {
		w.WriteByte('0')
	}
	w.WriteByte('\t')

	// TLEN (signed).
	w.Write(br.appendInt(int64(tlen)))
	w.WriteByte('\t')

	// SEQ: (lSeq+1)/2 packed bytes, high nibble first. "*" when empty.
	seqBytes := (lSeq + 1) / 2
	if off+seqBytes > len(buf) {
		return fmt.Errorf("sam: truncated SEQ")
	}
	if lSeq == 0 {
		w.WriteByte('*')
	} else {
		// Expand the packed nibbles into the reusable text scratch and write
		// the whole SEQ in one Write — far fewer bufio calls than a WriteByte
		// per base on this, the longest per-record field.
		packed := buf[off : off+seqBytes]
		s := br.growText(lSeq)
		for i := 0; i < lSeq; i++ {
			b := packed[i>>1]
			var nib byte
			if i&1 == 0 {
				nib = b >> 4
			} else {
				nib = b & 0x0f
			}
			s[i] = seqLookup[nib]
		}
		w.Write(s)
	}
	off += seqBytes
	w.WriteByte('\t')

	// QUAL: lSeq Phred bytes. "*" when empty or all 0xFF.
	if off+lSeq > len(buf) {
		return fmt.Errorf("sam: truncated QUAL")
	}
	qual := buf[off : off+lSeq]
	off += lSeq
	if lSeq == 0 || allBytes0xFF(qual) {
		w.WriteByte('*')
	} else {
		s := br.growText(lSeq)
		for i, q := range qual {
			s[i] = q + 33
		}
		w.Write(s)
	}

	// AUX: stream the remaining bytes, writing each tag straight to text.
	if off < len(buf) {
		if err := br.writeAuxSAM(w, buf[off:]); err != nil {
			return err
		}
	}
	return w.WriteByte('\n')
}

// writeAuxSAM walks the binary aux stream and writes each entry as
// "\tTAG:TYPE:VALUE", byte-identical to formatting the decoded Aux via
// Aux.FormatSAM. It handles every BAM aux type (A,c,C,s,S,i,I,f,Z,H,B with all
// B subtypes), mapping the integer widths to the SAM 'i' letter and floats via
// strconv's shortest %g (FormatFloat 'g', -1, 32), exactly as FormatSAM does.
func (br *BAMReader) writeAuxSAM(w *bufio.Writer, buf []byte) error {
	for len(buf) > 0 {
		if len(buf) < 3 {
			return fmt.Errorf("sam: truncated aux header")
		}
		w.WriteByte('\t')
		w.WriteByte(buf[0])
		w.WriteByte(buf[1])
		w.WriteByte(':')
		typ := buf[2]
		buf = buf[3:]
		switch typ {
		case 'A':
			if len(buf) < 1 {
				return fmt.Errorf("sam: truncated aux 'A'")
			}
			w.WriteString("A:")
			w.WriteByte(buf[0])
			buf = buf[1:]
		case 'c':
			if len(buf) < 1 {
				return fmt.Errorf("sam: truncated aux 'c'")
			}
			w.WriteString("i:")
			w.Write(br.appendInt(int64(int8(buf[0]))))
			buf = buf[1:]
		case 'C':
			if len(buf) < 1 {
				return fmt.Errorf("sam: truncated aux 'C'")
			}
			w.WriteString("i:")
			w.Write(br.appendInt(int64(buf[0])))
			buf = buf[1:]
		case 's':
			if len(buf) < 2 {
				return fmt.Errorf("sam: truncated aux 's'")
			}
			w.WriteString("i:")
			w.Write(br.appendInt(int64(int16(binary.LittleEndian.Uint16(buf[:2])))))
			buf = buf[2:]
		case 'S':
			if len(buf) < 2 {
				return fmt.Errorf("sam: truncated aux 'S'")
			}
			w.WriteString("i:")
			w.Write(br.appendInt(int64(binary.LittleEndian.Uint16(buf[:2]))))
			buf = buf[2:]
		case 'i':
			if len(buf) < 4 {
				return fmt.Errorf("sam: truncated aux 'i'")
			}
			w.WriteString("i:")
			w.Write(br.appendInt(int64(int32(binary.LittleEndian.Uint32(buf[:4])))))
			buf = buf[4:]
		case 'I':
			if len(buf) < 4 {
				return fmt.Errorf("sam: truncated aux 'I'")
			}
			w.WriteString("i:")
			w.Write(br.appendInt(int64(binary.LittleEndian.Uint32(buf[:4]))))
			buf = buf[4:]
		case 'f':
			if len(buf) < 4 {
				return fmt.Errorf("sam: truncated aux 'f'")
			}
			w.WriteString("f:")
			f := float64(math.Float32frombits(binary.LittleEndian.Uint32(buf[:4])))
			w.Write(strconv.AppendFloat(br.numScratch[:0], f, 'g', -1, 32))
			buf = buf[4:]
		case 'Z':
			end := indexZeroByte(buf)
			if end < 0 {
				return fmt.Errorf("sam: unterminated aux 'Z'")
			}
			w.WriteString("Z:")
			w.Write(buf[:end])
			buf = buf[end+1:]
		case 'H':
			end := indexZeroByte(buf)
			if end < 0 {
				return fmt.Errorf("sam: unterminated aux 'H'")
			}
			w.WriteString("H:")
			w.Write(buf[:end])
			buf = buf[end+1:]
		case 'B':
			if len(buf) < 5 {
				return fmt.Errorf("sam: truncated aux 'B' header")
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
				return fmt.Errorf("sam: unknown aux 'B' subtype %q", sub)
			}
			need := int(count) * elemSize
			if len(buf) < need {
				return fmt.Errorf("sam: truncated aux 'B' body")
			}
			w.WriteByte('B')
			w.WriteByte(':')
			w.WriteByte(sub)
			for j := 0; j < int(count); j++ {
				eoff := j * elemSize
				w.WriteByte(',')
				switch sub {
				case 'c':
					w.Write(br.appendInt(int64(int8(buf[eoff]))))
				case 'C':
					w.Write(br.appendInt(int64(buf[eoff])))
				case 's':
					w.Write(br.appendInt(int64(int16(binary.LittleEndian.Uint16(buf[eoff : eoff+2])))))
				case 'S':
					w.Write(br.appendInt(int64(binary.LittleEndian.Uint16(buf[eoff : eoff+2]))))
				case 'i':
					w.Write(br.appendInt(int64(int32(binary.LittleEndian.Uint32(buf[eoff : eoff+4])))))
				case 'I':
					w.Write(br.appendInt(int64(binary.LittleEndian.Uint32(buf[eoff : eoff+4]))))
				case 'f':
					f := float64(math.Float32frombits(binary.LittleEndian.Uint32(buf[eoff : eoff+4])))
					w.Write(strconv.AppendFloat(br.numScratch[:0], f, 'g', -1, 32))
				}
			}
			buf = buf[need:]
		default:
			return fmt.Errorf("sam: unknown aux type %q", typ)
		}
	}
	return nil
}

// appendUint formats v into the reader's reusable numeric scratch buffer and
// returns the slice. The bytes are valid until the next appendUint/appendInt
// call on the same reader; WriteSAMBody copies them into the bufio writer
// before the next call, so reuse is safe.
func (br *BAMReader) appendUint(v uint64) []byte {
	return strconv.AppendUint(br.numScratch[:0], v, 10)
}

// appendInt formats a signed v into the reader's reusable numeric scratch
// buffer; see appendUint for the validity contract.
func (br *BAMReader) appendInt(v int64) []byte {
	return strconv.AppendInt(br.numScratch[:0], v, 10)
}

// growText returns the reader's reusable text scratch sliced to exactly n
// bytes, growing the backing array when needed. The bytes are valid until the
// next growText call on the same reader; WriteSAMBody writes them to the bufio
// writer before reusing the buffer, so this is safe.
func (br *BAMReader) growText(n int) []byte {
	if cap(br.textScratch) < n {
		br.textScratch = make([]byte, n)
	} else {
		br.textScratch = br.textScratch[:n]
	}
	return br.textScratch
}

// indexZeroByte returns the index of the first NUL in b, or -1 if none.
func indexZeroByte(b []byte) int {
	for i := 0; i < len(b); i++ {
		if b[i] == 0 {
			return i
		}
	}
	return -1
}

// allBytes0xFF reports whether every byte equals 0xFF (the SAM "no quality"
// sentinel, written as "*").
func allBytes0xFF(q []byte) bool {
	for _, b := range q {
		if b != 0xff {
			return false
		}
	}
	return true
}
