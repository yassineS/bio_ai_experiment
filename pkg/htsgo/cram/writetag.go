package cram

import (
	"encoding/binary"
	"fmt"
	"math"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// cramTagType returns the one-byte CRAM value-type letter for an
// auxiliary field. It is the third byte of the tag's three-byte key and
// selects how decodeTagValue interprets the stored bytes. The writer
// uses the aux field's own Type letter so a BAM-sourced narrow integer
// ('c', 'S', …) round-trips at its original width; an aux with no usable
// type is reported as an error by encodeTagValue.
func cramTagType(a sam.Aux) byte {
	switch a.Type {
	case 'A', 'c', 'C', 's', 'S', 'i', 'I', 'f', 'Z', 'H', 'B':
		return a.Type
	default:
		// An aux parsed from text SAM with no width hint is a plain
		// integer; fall back to 'i'.
		if _, ok := a.Value.(int64); ok {
			return 'i'
		}
		return a.Type
	}
}

// encodeTagValue serialises one auxiliary field's value into the raw
// byte form CRAM stores for that tag's data series — BAM binary layout
// with no leading type letter. It is the writer-side inverse of
// decodeTagValue: the produced bytes, fed back through decodeTagValue
// with the same three-byte key, reconstruct an equal sam.Aux.
//
// A value whose Go type does not match its declared aux Type is
// reported as an error so the writer never emits a tag it cannot
// round-trip. The tag series are BYTE_ARRAY_LEN encoded, so a value may
// contain any byte.
func encodeTagValue(a sam.Aux) ([]byte, error) {
	typ := cramTagType(a)
	switch typ {
	case 'A':
		s, ok := a.Value.(string)
		if !ok || len(s) != 1 {
			return nil, fmt.Errorf("'A' tag value must be a one-character string")
		}
		return []byte{s[0]}, nil
	case 'c', 'C':
		v, err := auxInt(a)
		if err != nil {
			return nil, err
		}
		return []byte{byte(v)}, nil
	case 's', 'S':
		v, err := auxInt(a)
		if err != nil {
			return nil, err
		}
		buf := make([]byte, 2)
		binary.LittleEndian.PutUint16(buf, uint16(v))
		return buf, nil
	case 'i', 'I':
		v, err := auxInt(a)
		if err != nil {
			return nil, err
		}
		buf := make([]byte, 4)
		binary.LittleEndian.PutUint32(buf, uint32(v))
		return buf, nil
	case 'f':
		f, ok := a.Value.(float64)
		if !ok {
			return nil, fmt.Errorf("'f' tag value must be a float")
		}
		buf := make([]byte, 4)
		binary.LittleEndian.PutUint32(buf, math.Float32bits(float32(f)))
		return buf, nil
	case 'Z', 'H':
		s, ok := a.Value.(string)
		if !ok {
			return nil, fmt.Errorf("%q tag value must be a string", typ)
		}
		out := append([]byte(s), 0) // NUL-terminated in BAM/CRAM.
		return out, nil
	case 'B':
		return encodeArrayTag(a)
	default:
		return nil, fmt.Errorf("unsupported aux value type %q", typ)
	}
}

// auxInt extracts the int64 value of an integer-typed auxiliary field.
func auxInt(a sam.Aux) (int64, error) {
	v, ok := a.Value.(int64)
	if !ok {
		return 0, fmt.Errorf("integer tag value is not stored as int64")
	}
	return v, nil
}

// encodeArrayTag serialises a 'B' (array) auxiliary value: the one-byte
// element subtype, a 4-byte little-endian element count, then the
// elements in that subtype's binary layout. It mirrors decodeArrayTag.
func encodeArrayTag(a sam.Aux) ([]byte, error) {
	elemSize, ok := arrayElemSize(a.ArrayType)
	if !ok {
		return nil, fmt.Errorf("'B' tag has unknown element subtype %q", a.ArrayType)
	}
	out := make([]byte, 5+len(a.ArrayValues)*elemSize)
	out[0] = a.ArrayType
	binary.LittleEndian.PutUint32(out[1:5], uint32(len(a.ArrayValues)))
	for i, v := range a.ArrayValues {
		off := 5 + i*elemSize
		switch a.ArrayType {
		case 'c', 'C':
			n, err := arrayInt(v)
			if err != nil {
				return nil, err
			}
			out[off] = byte(n)
		case 's', 'S':
			n, err := arrayInt(v)
			if err != nil {
				return nil, err
			}
			binary.LittleEndian.PutUint16(out[off:], uint16(n))
		case 'i', 'I':
			n, err := arrayInt(v)
			if err != nil {
				return nil, err
			}
			binary.LittleEndian.PutUint32(out[off:], uint32(n))
		case 'f':
			f, ok := v.(float64)
			if !ok {
				return nil, fmt.Errorf("'B:f' element is not a float")
			}
			binary.LittleEndian.PutUint32(out[off:], math.Float32bits(float32(f)))
		}
	}
	return out, nil
}

// arrayInt extracts the int64 value of one 'B'-array integer element.
func arrayInt(v interface{}) (int64, error) {
	n, ok := v.(int64)
	if !ok {
		return 0, fmt.Errorf("'B' integer element is not stored as int64")
	}
	return n, nil
}
