package cram

import (
	"encoding/binary"
	"math"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// decodeTagValue interprets the raw bytes of one CRAM tag value into a
// sam.Aux. The three-byte tag key carries the two-character SAM tag name
// and the one-byte BAM value-type letter; the raw bytes are the value in
// BAM binary layout (little-endian, no type prefix), exactly as CRAM
// stores a tag's data series.
func decodeTagValue(key tagKey, raw []byte) (sam.Aux, error) {
	aux := sam.Aux{Tag: string(key[:2]), Type: key[2]}
	switch key[2] {
	case 'A':
		if len(raw) != 1 {
			return aux, errFormat("tag 'A' value must be one byte, got %d", len(raw))
		}
		aux.Value = string(raw[:1])
	case 'c':
		v, err := need(raw, 1, "c")
		if err != nil {
			return aux, err
		}
		aux.Value = int64(int8(v[0]))
	case 'C':
		v, err := need(raw, 1, "C")
		if err != nil {
			return aux, err
		}
		aux.Value = int64(v[0])
	case 's':
		v, err := need(raw, 2, "s")
		if err != nil {
			return aux, err
		}
		aux.Value = int64(int16(binary.LittleEndian.Uint16(v)))
	case 'S':
		v, err := need(raw, 2, "S")
		if err != nil {
			return aux, err
		}
		aux.Value = int64(binary.LittleEndian.Uint16(v))
	case 'i':
		v, err := need(raw, 4, "i")
		if err != nil {
			return aux, err
		}
		aux.Value = int64(int32(binary.LittleEndian.Uint32(v)))
	case 'I':
		v, err := need(raw, 4, "I")
		if err != nil {
			return aux, err
		}
		aux.Value = int64(binary.LittleEndian.Uint32(v))
	case 'f':
		v, err := need(raw, 4, "f")
		if err != nil {
			return aux, err
		}
		aux.Value = float64(math.Float32frombits(binary.LittleEndian.Uint32(v)))
	case 'Z', 'H':
		// Z and H values are NUL-terminated in BAM/CRAM; drop a single
		// trailing NUL if present.
		s := raw
		if len(s) > 0 && s[len(s)-1] == 0 {
			s = s[:len(s)-1]
		}
		aux.Value = string(s)
	case 'B':
		return decodeArrayTag(aux, raw)
	default:
		return aux, errFormat("unknown tag value type %#02x (%q)", key[2], string(key[2]))
	}
	return aux, nil
}

// decodeArrayTag interprets a 'B' (array) tag value. The CRAM layout
// stores the value without the leading BAM type letter: a one-byte
// element subtype, a 4-byte little-endian count, then count elements of
// that subtype.
func decodeArrayTag(aux sam.Aux, raw []byte) (sam.Aux, error) {
	if len(raw) < 5 {
		return aux, errFormat("tag 'B' value truncated (%d bytes, need at least 5)", len(raw))
	}
	sub := raw[0]
	count := binary.LittleEndian.Uint32(raw[1:5])
	aux.ArrayType = sub
	body := raw[5:]
	elemSize, ok := arrayElemSize(sub)
	if !ok {
		return aux, errFormat("tag 'B' has unknown element subtype %#02x (%q)", sub, string(sub))
	}
	if uint64(count)*uint64(elemSize) != uint64(len(body)) {
		return aux, errFormat("tag 'B:%c' declares %d elements (%d bytes) but %d bytes follow",
			sub, count, uint64(count)*uint64(elemSize), len(body))
	}
	aux.ArrayValues = make([]interface{}, 0, count)
	for i := uint32(0); i < count; i++ {
		off := int(i) * elemSize
		switch sub {
		case 'c':
			aux.ArrayValues = append(aux.ArrayValues, int64(int8(body[off])))
		case 'C':
			aux.ArrayValues = append(aux.ArrayValues, int64(body[off]))
		case 's':
			aux.ArrayValues = append(aux.ArrayValues, int64(int16(binary.LittleEndian.Uint16(body[off:]))))
		case 'S':
			aux.ArrayValues = append(aux.ArrayValues, int64(binary.LittleEndian.Uint16(body[off:])))
		case 'i':
			aux.ArrayValues = append(aux.ArrayValues, int64(int32(binary.LittleEndian.Uint32(body[off:]))))
		case 'I':
			aux.ArrayValues = append(aux.ArrayValues, int64(binary.LittleEndian.Uint32(body[off:])))
		case 'f':
			aux.ArrayValues = append(aux.ArrayValues, float64(math.Float32frombits(binary.LittleEndian.Uint32(body[off:]))))
		}
	}
	return aux, nil
}

// arrayElemSize returns the byte width of one 'B'-array element of the
// given subtype.
func arrayElemSize(sub byte) (int, bool) {
	switch sub {
	case 'c', 'C':
		return 1, true
	case 's', 'S':
		return 2, true
	case 'i', 'I', 'f':
		return 4, true
	}
	return 0, false
}

// need returns the first n bytes of raw, erroring when fewer are
// present. The label names the tag type in the error message.
func need(raw []byte, n int, label string) ([]byte, error) {
	if len(raw) < n {
		return nil, errFormat("tag %q value truncated (%d bytes, need %d)", label, len(raw), n)
	}
	return raw[:n], nil
}
