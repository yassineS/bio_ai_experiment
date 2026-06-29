package sam

import (
	"errors"
	"fmt"
	"strconv"
	"strings"
)

// SAM flag bits, mirroring the SAMv1 spec.
const (
	FlagPaired        uint16 = 0x1
	FlagProperPair    uint16 = 0x2
	FlagUnmapped      uint16 = 0x4
	FlagMateUnmapped  uint16 = 0x8
	FlagReverse       uint16 = 0x10
	FlagMateReverse   uint16 = 0x20
	FlagRead1         uint16 = 0x40
	FlagRead2         uint16 = 0x80
	FlagSecondary     uint16 = 0x100
	FlagQCFail        uint16 = 0x200
	FlagDuplicate     uint16 = 0x400
	FlagSupplementary uint16 = 0x800
)

// CIGAR op codes (4-bit values inside packed CIGAR uint32).
const (
	CigarMatch     = 0 // M
	CigarInsertion = 1 // I
	CigarDeletion  = 2 // D
	CigarSkipped   = 3 // N
	CigarSoftClip  = 4 // S
	CigarHardClip  = 5 // H
	CigarPadding   = 6 // P
	CigarEqual     = 7 // =
	CigarMismatch  = 8 // X
	CigarBack      = 9 // B (rarely used, present for completeness)
)

// cigarOpChars maps a 4-bit CIGAR op to its character representation.
var cigarOpChars = [...]byte{'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B'}

// cigarOpConsumesQuery reports whether the given op consumes bases from the
// query (read) sequence — i.e. advances the position within SEQ/QUAL.
func cigarOpConsumesQuery(op uint32) bool {
	switch op {
	case CigarMatch, CigarInsertion, CigarSoftClip, CigarEqual, CigarMismatch:
		return true
	}
	return false
}

// cigarOpConsumesRef reports whether the given op consumes bases from the
// reference — i.e. advances the alignment position on the reference.
func cigarOpConsumesRef(op uint32) bool {
	switch op {
	case CigarMatch, CigarDeletion, CigarSkipped, CigarEqual, CigarMismatch:
		return true
	}
	return false
}

// CigarOp is a single packed CIGAR operation: the bottom 4 bits are the op
// code (one of the Cigar* constants) and the top 28 bits are the length.
type CigarOp uint32

// Op returns the operation code.
func (c CigarOp) Op() uint32 { return uint32(c) & 0xf }

// Length returns the operation length.
func (c CigarOp) Length() uint32 { return uint32(c) >> 4 }

// Char returns the SAM character for the operation (e.g. 'M', 'I').
func (c CigarOp) Char() byte {
	if int(c.Op()) >= len(cigarOpChars) {
		return '?'
	}
	return cigarOpChars[c.Op()]
}

// Cigar is a sequence of CigarOps with helper methods for parsing, encoding,
// and length math.
type Cigar []CigarOp

// String returns the standard SAM textual encoding of the CIGAR (e.g. "10M2I88M").
// An empty CIGAR returns "*", matching the SAM convention.
func (c Cigar) String() string {
	if len(c) == 0 {
		return "*"
	}
	var sb strings.Builder
	for _, op := range c {
		sb.WriteString(strconv.FormatUint(uint64(op.Length()), 10))
		sb.WriteByte(op.Char())
	}
	return sb.String()
}

// QueryLength returns the number of query (SEQ) bases consumed by the CIGAR.
func (c Cigar) QueryLength() int {
	n := 0
	for _, op := range c {
		if cigarOpConsumesQuery(op.Op()) {
			n += int(op.Length())
		}
	}
	return n
}

// ReferenceLength returns the number of reference bases consumed by the CIGAR.
// Used to compute the end coordinate of an alignment.
func (c Cigar) ReferenceLength() int {
	n := 0
	for _, op := range c {
		if cigarOpConsumesRef(op.Op()) {
			n += int(op.Length())
		}
	}
	return n
}

// ParseCigar parses a textual CIGAR like "10M2I88M". The single character "*"
// is recognised as the empty (unavailable) CIGAR and returns a nil slice.
func ParseCigar(s string) (Cigar, error) {
	if s == "" {
		return nil, fmt.Errorf("sam: empty CIGAR")
	}
	if s == "*" {
		return nil, nil
	}
	var ops Cigar
	i := 0
	for i < len(s) {
		j := i
		for j < len(s) && s[j] >= '0' && s[j] <= '9' {
			j++
		}
		if j == i || j == len(s) {
			return nil, fmt.Errorf("sam: malformed CIGAR %q", s)
		}
		length, err := strconv.ParseUint(s[i:j], 10, 28)
		if err != nil {
			return nil, fmt.Errorf("sam: bad CIGAR length in %q: %w", s, err)
		}
		op, ok := cigarCharOp(s[j])
		if !ok {
			return nil, fmt.Errorf("sam: bad CIGAR op %q in %q", s[j], s)
		}
		ops = append(ops, CigarOp(uint32(length)<<4|op))
		i = j + 1
	}
	return ops, nil
}

// cigarCharOp maps a CIGAR character to its 4-bit op code.
func cigarCharOp(c byte) (uint32, bool) {
	switch c {
	case 'M':
		return CigarMatch, true
	case 'I':
		return CigarInsertion, true
	case 'D':
		return CigarDeletion, true
	case 'N':
		return CigarSkipped, true
	case 'S':
		return CigarSoftClip, true
	case 'H':
		return CigarHardClip, true
	case 'P':
		return CigarPadding, true
	case '=':
		return CigarEqual, true
	case 'X':
		return CigarMismatch, true
	case 'B':
		return CigarBack, true
	}
	return 0, false
}

// Aux is one optional auxiliary field on an alignment record.
//
// Type encodes the SAM aux type letter:
//   - 'A' single ASCII character (Value is a string of length 1)
//   - 'c'/'C' signed/unsigned 8-bit integer (Value stored as int64)
//   - 's'/'S' signed/unsigned 16-bit integer (Value stored as int64)
//   - 'i'/'I' signed/unsigned 32-bit integer (Value stored as int64)
//   - 'f'     single-precision float (Value is float32 in float64)
//   - 'Z'     string
//   - 'H'     hex string
//   - 'B'     array; ArrayType is the element subtype and ArrayValues is a
//     slice of the parsed values (int64 for integer subtypes, float64
//     for 'f').
//
// In text SAM, integer aux values use the type letter 'i' regardless of width;
// when a value is read from BAM with type 'c', 's', etc, the original width is
// preserved in Type so SAM→BAM round-trips choose the most compact encoding.
type Aux struct {
	Tag         string
	Type        byte
	Value       interface{}
	ArrayType   byte
	ArrayValues []interface{}
}

// Int returns the aux value as an int64 if it is an integer type, and reports
// whether the conversion was possible.
func (a Aux) Int() (int64, bool) {
	switch a.Type {
	case 'c', 'C', 's', 'S', 'i', 'I':
		if v, ok := a.Value.(int64); ok {
			return v, true
		}
	}
	return 0, false
}

// String returns the aux value as a string if Type is 'Z' or 'H', otherwise
// returns "" and false.
func (a Aux) String() (string, bool) {
	if a.Type == 'Z' || a.Type == 'H' {
		if v, ok := a.Value.(string); ok {
			return v, true
		}
	}
	return "", false
}

// FormatSAM returns the text SAM encoding of this aux field, e.g. "NM:i:3".
func (a Aux) FormatSAM() string {
	var sb strings.Builder
	sb.WriteString(a.Tag)
	sb.WriteByte(':')
	switch a.Type {
	case 'A':
		sb.WriteString("A:")
		if s, ok := a.Value.(string); ok && len(s) > 0 {
			sb.WriteByte(s[0])
		}
	case 'c', 'C', 's', 'S', 'i', 'I':
		sb.WriteString("i:")
		if v, ok := a.Value.(int64); ok {
			sb.WriteString(strconv.FormatInt(v, 10))
		}
	case 'f':
		sb.WriteString("f:")
		if v, ok := a.Value.(float64); ok {
			sb.WriteString(strconv.FormatFloat(v, 'g', -1, 32))
		}
	case 'Z':
		sb.WriteString("Z:")
		if v, ok := a.Value.(string); ok {
			sb.WriteString(v)
		}
	case 'H':
		sb.WriteString("H:")
		if v, ok := a.Value.(string); ok {
			sb.WriteString(v)
		}
	case 'B':
		sb.WriteString("B:")
		sb.WriteByte(a.ArrayType)
		for _, v := range a.ArrayValues {
			sb.WriteByte(',')
			switch a.ArrayType {
			case 'f':
				if f, ok := v.(float64); ok {
					sb.WriteString(strconv.FormatFloat(f, 'g', -1, 32))
				}
			default:
				if i, ok := v.(int64); ok {
					sb.WriteString(strconv.FormatInt(i, 10))
				}
			}
		}
	}
	return sb.String()
}

// ParseAux parses a single "TAG:TYPE:VALUE" SAM aux field.
func ParseAux(field string) (Aux, error) {
	if len(field) < 5 || field[2] != ':' || field[4] != ':' {
		return Aux{}, fmt.Errorf("sam: bad aux field %q", field)
	}
	a := Aux{Tag: field[:2], Type: field[3]}
	val := field[5:]
	switch a.Type {
	case 'A':
		if len(val) != 1 {
			return Aux{}, fmt.Errorf("sam: aux 'A' value must be one char, got %q", val)
		}
		a.Value = val
	case 'i':
		v, err := strconv.ParseInt(val, 10, 64)
		if err != nil {
			return Aux{}, fmt.Errorf("sam: aux 'i' parse: %w", err)
		}
		a.Value = v
	case 'f':
		v, err := strconv.ParseFloat(val, 64)
		if err != nil {
			return Aux{}, fmt.Errorf("sam: aux 'f' parse: %w", err)
		}
		a.Value = v
	case 'Z', 'H':
		a.Value = val
	case 'B':
		if len(val) < 2 || val[1] != ',' && len(val) != 1 {
			return Aux{}, fmt.Errorf("sam: aux 'B' bad header in %q", val)
		}
		a.ArrayType = val[0]
		rest := ""
		if len(val) > 1 {
			rest = val[1:]
		}
		if rest != "" && rest[0] == ',' {
			rest = rest[1:]
		}
		if rest != "" {
			for _, p := range strings.Split(rest, ",") {
				switch a.ArrayType {
				case 'f':
					f, err := strconv.ParseFloat(p, 64)
					if err != nil {
						return Aux{}, fmt.Errorf("sam: aux 'B:f' parse: %w", err)
					}
					a.ArrayValues = append(a.ArrayValues, f)
				default:
					n, err := strconv.ParseInt(p, 10, 64)
					if err != nil {
						return Aux{}, fmt.Errorf("sam: aux 'B:%c' parse: %w", a.ArrayType, err)
					}
					a.ArrayValues = append(a.ArrayValues, n)
				}
			}
		}
	default:
		return Aux{}, fmt.Errorf("sam: unknown aux type %q", a.Type)
	}
	return a, nil
}

// Record is one SAM/BAM alignment record. All fields use Go-friendly types
// (1-based POS as stored in SAM, with -1 meaning "*" for RNAME). When read
// from BAM, fields preserve the binary types where appropriate.
//
// Pos, PNext and TLen are int64 (htslib's hts_pos_t) so SAM and CRAM, which
// support reference coordinates beyond 2^31, round-trip losslessly. BAM stores
// POS/PNEXT/TLEN in 32-bit on-disk fields, so the BAM writer rejects a record
// whose coordinate exceeds that range rather than truncating it — a format
// limit, matching htslib.
type Record struct {
	QName string
	Flag  uint16
	RName string
	Pos   int64 // 1-based position; 0 means unmapped (BAM stores 0-based + 1)
	MapQ  uint8
	Cigar Cigar
	RNext string
	PNext int64
	TLen  int64
	Seq   string
	Qual  []byte // raw Phred scores; nil/length zero or all 0xff means "*"
	Aux   []Aux

	// RawAux is the record's auxiliary fields as a raw on-disk BAM aux byte
	// block (a run of tag[2]+type[1]+value entries in BAM binary layout). When
	// it is non-nil the BAM writer emits it verbatim and Aux is left nil; it is
	// mutually exclusive with Aux. It is a memory-lean passthrough used on the
	// CRAM→BAM view path, where the aux bytes already arrive in BAM layout and
	// never need to be materialised into the heavier []Aux form. Any consumer
	// that reads Aux (GetAux, the SAM/CRAM writers) self-heals by decoding
	// RawAux into Aux on first access, so correctness never depends on the
	// passthrough gate being perfect.
	RawAux []byte

	// RawSeq is the record's SEQ as the on-disk BAM 4-bit-packed nibble block
	// ((SeqLen+1)/2 bytes, high nibble first). When it is non-nil the BAM writer
	// emits it verbatim — skipping the ASCII→nibble pack — and Seq is left "".
	// It is the SEQ analogue of RawAux: a memory-lean passthrough used on the
	// CRAM→BAM view path, where SEQ is the dominant fat per-record field. It is
	// mutually exclusive with a populated Seq; SeqLen carries the base count so
	// the writer's l_seq and the QUAL length are known without unpacking. Any
	// consumer that reads Seq (the SAM/CRAM writers, mdnm) self-heals via
	// MaterialiseSeq or reads the nibbles directly through SeqBaseAt/SeqBytes, so
	// correctness never depends on the passthrough gate being perfect.
	RawSeq []byte

	// SeqLen is the number of SEQ bases when RawSeq carries the packed sequence.
	// It is meaningful only while RawSeq != nil; otherwise len(Seq) is the count.
	SeqLen int

	// AuxIndex maps a 2-char tag to its position in Aux for quick lookups.
	// Populated lazily by GetAux.
	auxIndex map[string]int
}

// materialiseAux decodes a pending RawAux byte block into the Aux slice when
// Aux has not yet been populated, then clears RawAux so the two never disagree.
// It is the defensive lazy guard that lets any Aux-reading code path stay
// correct even if a record reached it with RawAux still set (the gated view
// passthrough should ensure that never happens, but this makes correctness
// independent of the gate). A record with neither RawAux nor Aux is left
// untouched. The decoded Aux is byte-for-byte what decodeBAMAux would produce
// for the same bytes, so the materialised record is identical to one decoded
// eagerly.
func (r *Record) materialiseAux() {
	if r.Aux == nil && r.RawAux != nil {
		r.Aux, _ = decodeBAMAuxInto(nil, r.RawAux)
		r.RawAux = nil
	}
}

// MaterialiseAux decodes a pending RawAux byte block into the Aux slice, then
// clears RawAux. It is the exported form of the defensive lazy guard for
// consumers in other packages (e.g. the CRAM writer) that read Aux directly: a
// record handed to them with RawAux set self-heals into the decoded Aux form.
// It is a no-op when Aux is already populated or no aux is present.
func (r *Record) MaterialiseAux() { r.materialiseAux() }

// materialiseSeq unpacks a pending RawSeq nibble block into the Seq string when
// Seq has not yet been populated, then clears RawSeq so the two never disagree.
// It is the SEQ analogue of materialiseAux: the defensive lazy guard that lets
// any Seq-reading code path stay correct even if a record reached it with
// RawSeq still set (the gated view passthrough should ensure that never happens,
// but this makes correctness independent of the gate). A record with neither
// RawSeq nor an eager Seq is left untouched. The unpacked Seq is byte-for-byte
// what the BAM reader would produce for the same nibbles.
func (r *Record) materialiseSeq() {
	if r.Seq == "" && r.RawSeq != nil {
		r.Seq = string(UnpackSeq(r.RawSeq, r.SeqLen))
		r.RawSeq = nil
		r.SeqLen = 0
	}
}

// MaterialiseSeq unpacks a pending RawSeq nibble block into the Seq string, then
// clears RawSeq. It is the exported form of the defensive lazy guard for
// consumers in other packages (e.g. the CRAM writer) that read Seq directly: a
// record handed to them with RawSeq set self-heals into the unpacked Seq form.
// It is a no-op when Seq is already populated or no SEQ is present.
func (r *Record) MaterialiseSeq() { r.materialiseSeq() }

// SeqLength returns the number of SEQ bases regardless of which representation
// the record currently holds: SeqLen when SEQ is carried packed in RawSeq, else
// len(Seq).
func (r *Record) SeqLength() int {
	if r.RawSeq != nil {
		return r.SeqLen
	}
	return len(r.Seq)
}

// SeqBaseAt returns the ASCII nucleotide of the SEQ base at index i, reading
// directly from the packed RawSeq nibbles when SEQ is carried packed, or from
// the Seq string otherwise. It lets a consumer (e.g. mdnm) walk the sequence
// without forcing RawSeq to be materialised into a string. i must be in
// [0, SeqLength()).
func (r *Record) SeqBaseAt(i int) byte {
	if r.RawSeq != nil {
		b := r.RawSeq[i/2]
		var nibble byte
		if i%2 == 0 {
			nibble = b >> 4
		} else {
			nibble = b & 0x0f
		}
		return seqLookup[nibble]
	}
	return r.Seq[i]
}

// GetAux returns the aux field with the given tag, or false if absent.
func (r *Record) GetAux(tag string) (Aux, bool) {
	r.materialiseAux()
	if r.auxIndex == nil && len(r.Aux) > 0 {
		r.auxIndex = make(map[string]int, len(r.Aux))
		for i, a := range r.Aux {
			r.auxIndex[a.Tag] = i
		}
	}
	if i, ok := r.auxIndex[tag]; ok {
		return r.Aux[i], true
	}
	return Aux{}, false
}

// InvalidateAuxIndex drops the cached tag→position lookup map. Call it after
// directly mutating, appending to, or deleting from r.Aux so that the next
// GetAux rebuilds the index from the current slice contents.
func (r *Record) InvalidateAuxIndex() { r.auxIndex = nil }

// IsUnmapped reports whether the read is flagged as unmapped (flag bit 0x4).
func (r *Record) IsUnmapped() bool { return r.Flag&FlagUnmapped != 0 }

// IsMapped reports whether the read is mapped (i.e. not unmapped).
func (r *Record) IsMapped() bool { return r.Flag&FlagUnmapped == 0 }

// IsPaired reports whether the read is paired (flag bit 0x1).
func (r *Record) IsPaired() bool { return r.Flag&FlagPaired != 0 }

// IsProperPair reports whether the read is in a proper pair (flag bit 0x2).
func (r *Record) IsProperPair() bool { return r.Flag&FlagProperPair != 0 }

// IsMateUnmapped reports whether the mate is unmapped (flag bit 0x8).
func (r *Record) IsMateUnmapped() bool { return r.Flag&FlagMateUnmapped != 0 }

// IsRead1 reports whether the read is the first in pair (flag bit 0x40).
func (r *Record) IsRead1() bool { return r.Flag&FlagRead1 != 0 }

// IsRead2 reports whether the read is the second in pair (flag bit 0x80).
func (r *Record) IsRead2() bool { return r.Flag&FlagRead2 != 0 }

// IsSecondary reports whether the alignment is secondary (flag bit 0x100).
func (r *Record) IsSecondary() bool { return r.Flag&FlagSecondary != 0 }

// IsSupplementary reports whether the alignment is supplementary (flag bit 0x800).
func (r *Record) IsSupplementary() bool { return r.Flag&FlagSupplementary != 0 }

// IsDuplicate reports whether the read is marked a duplicate (flag bit 0x400).
func (r *Record) IsDuplicate() bool { return r.Flag&FlagDuplicate != 0 }

// IsQCFail reports whether the read failed quality controls (flag bit 0x200).
func (r *Record) IsQCFail() bool { return r.Flag&FlagQCFail != 0 }

// IsPrimary reports whether the alignment is the primary alignment, which by
// definition is neither secondary nor supplementary.
func (r *Record) IsPrimary() bool {
	return r.Flag&(FlagSecondary|FlagSupplementary) == 0
}

// EndPosition returns the 1-based inclusive end coordinate of the alignment
// on the reference, computed from Pos + Cigar.ReferenceLength() - 1. Returns
// Pos when the CIGAR consumes no reference bases (e.g. empty / pure insertion).
func (r *Record) EndPosition() int64 {
	n := r.Cigar.ReferenceLength()
	if n == 0 {
		return r.Pos
	}
	return r.Pos + int64(n) - 1
}

// PackedSeq returns the BAM 4-bit-packed encoding of the record's SEQ:
// (len(Seq)+1)/2 bytes with the high nibble holding the first base. It is
// byte-identical to the buffer pointed at by htslib's bam_get_seq, which makes
// it suitable for checksumming. A SEQ of "*" (empty Seq) yields an empty slice.
// When SEQ is carried packed in RawSeq the bytes are already in this layout, so
// they are returned directly without an unpack/repack round-trip.
func (r *Record) PackedSeq() []byte {
	if r.RawSeq != nil {
		return r.RawSeq
	}
	return encodeSeq(r.Seq)
}

// errBadRecord is returned for malformed text SAM body lines.
var errBadRecord = errors.New("sam: malformed record")
