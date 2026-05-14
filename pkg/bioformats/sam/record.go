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
type Record struct {
	QName string
	Flag  uint16
	RName string
	Pos   int32 // 1-based position; 0 means unmapped (BAM stores 0-based + 1)
	MapQ  uint8
	Cigar Cigar
	RNext string
	PNext int32
	TLen  int32
	Seq   string
	Qual  []byte // raw Phred scores; nil/length zero or all 0xff means "*"
	Aux   []Aux

	// AuxIndex maps a 2-char tag to its position in Aux for quick lookups.
	// Populated lazily by GetAux.
	auxIndex map[string]int
}

// GetAux returns the aux field with the given tag, or false if absent.
func (r *Record) GetAux(tag string) (Aux, bool) {
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
func (r *Record) EndPosition() int32 {
	n := r.Cigar.ReferenceLength()
	if n == 0 {
		return r.Pos
	}
	return r.Pos + int32(n) - 1
}

// errBadRecord is returned for malformed text SAM body lines.
var errBadRecord = errors.New("sam: malformed record")
