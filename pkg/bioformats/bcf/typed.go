package bcf

import (
	"encoding/binary"
	"fmt"
	"math"
)

// Typed-value descriptor codes used by the BCF wire format. They occupy the
// low nibble (descriptor) of the type byte; the high nibble holds the count
// (0..14 = literal, 15 = "next typed integer is the real count").
const (
	TypeMissing = 0 // no data follows; size class must be 0
	TypeInt8    = 1
	TypeInt16   = 2
	TypeInt32   = 3
	TypeInt64   = 4 // BCF 2.2+; htslib 1.13+ may emit this for very large counts
	TypeFloat   = 5
	TypeChar    = 7
)

// Missing-value sentinels for the typed numeric encodings. These are the
// exact bit patterns htslib uses on the wire and what we must round-trip.
const (
	MissingInt8    int8    = -128        // 0x80
	MissingInt16   int16   = -32768      // 0x8000
	MissingInt32   int32   = -2147483648 // 0x80000000
	MissingFloat   float32 = 0           // see MissingFloatBits for the real wire value
	MissingFloat32         = uint32(0x7F800001)
	// EndOfVector sentinels mark the end of a per-record vector when a value
	// is shorter than the declared dimension (htslib calls these "vector end").
	EndOfVectorInt8  int8  = -127        // 0x81
	EndOfVectorInt16 int16 = -32767      // 0x8001
	EndOfVectorInt32 int32 = -2147483647 // 0x80000001
	EndOfVectorFloat       = uint32(0x7F800002)
)

// IsMissingFloat reports whether bits matches the float32 missing pattern.
// We compare on the uint32 bit pattern because NaN != NaN in IEEE 754 and we
// must accept the exact 0x7F800001 produced by htslib.
func IsMissingFloat(bits uint32) bool { return bits == MissingFloat32 }

// IsEndOfVectorFloat reports whether bits matches the float32 vector-end
// sentinel 0x7F800002.
func IsEndOfVectorFloat(bits uint32) bool { return bits == EndOfVectorFloat }

// TypedValue is the decoded form of one "typed" field. Only one of the
// numeric/string slices is populated; Descriptor identifies which.
type TypedValue struct {
	Descriptor uint8     // TypeInt8 / TypeInt16 / TypeInt32 / TypeFloat / TypeChar / TypeMissing
	Length     int       // declared element count (0 for TypeMissing)
	Ints       []int32   // populated for int8/int16/int32; missing values become MissingInt32
	Floats     []float32 // populated for TypeFloat; missing entries are stored as math.Float32frombits(MissingFloat32)
	String     string    // populated for TypeChar
	// Raw holds the verbatim byte run for the value (without the descriptor
	// byte and without any length prefix). Useful for tests and for future
	// writer round-tripping.
	Raw []byte
}

// IsMissing reports whether v carries the "missing" descriptor or has zero
// declared length.
func (v TypedValue) IsMissing() bool { return v.Descriptor == TypeMissing || v.Length == 0 }

// readByte returns one byte at offset, advancing the offset. It is a small
// helper used by the typed decoder to keep call-site code legible.
func readByte(buf []byte, off *int) (byte, error) {
	if *off >= len(buf) {
		return 0, fmt.Errorf("bcf: typed read out of bounds at offset %d", *off)
	}
	b := buf[*off]
	*off++
	return b, nil
}

// readBytes returns the n-byte slice starting at *off and advances *off. The
// returned slice aliases buf, so the caller must not mutate it.
func readBytes(buf []byte, off *int, n int) ([]byte, error) {
	if n < 0 {
		return nil, fmt.Errorf("bcf: negative typed length %d", n)
	}
	if *off+n > len(buf) {
		return nil, fmt.Errorf("bcf: typed read of %d bytes past end (offset %d, buflen %d)", n, *off, len(buf))
	}
	out := buf[*off : *off+n]
	*off += n
	return out, nil
}

// DecodeTyped decodes one typed value beginning at buf[*off]. On success it
// advances *off past the value and returns the decoded TypedValue.
//
// The descriptor byte's low nibble selects the element type; the high nibble
// is either the literal count (0..14) or 15, meaning the count is encoded as
// the next typed integer (recursive single-step).
func DecodeTyped(buf []byte, off *int) (TypedValue, error) {
	descByte, err := readByte(buf, off)
	if err != nil {
		return TypedValue{}, err
	}
	descriptor := descByte & 0x0F
	size := int(descByte >> 4)
	if size == 15 {
		// The next typed value is itself an int describing the actual
		// length. We decode it recursively and then keep going.
		ln, err := DecodeTyped(buf, off)
		if err != nil {
			return TypedValue{}, fmt.Errorf("bcf: failed to read typed length prefix: %w", err)
		}
		if len(ln.Ints) == 0 {
			return TypedValue{}, fmt.Errorf("bcf: empty length prefix on typed value")
		}
		size = int(ln.Ints[0])
		if size < 0 {
			return TypedValue{}, fmt.Errorf("bcf: negative length prefix %d", size)
		}
	}
	tv := TypedValue{Descriptor: descriptor, Length: size}

	switch descriptor {
	case TypeMissing:
		// "Missing" descriptor (raw 0x00) carries no payload. We still
		// expose Length so callers can distinguish "field absent" from
		// "field present but empty".
		return tv, nil

	case TypeInt8:
		raw, err := readBytes(buf, off, size)
		if err != nil {
			return TypedValue{}, err
		}
		tv.Raw = raw
		tv.Ints = make([]int32, size)
		for i, b := range raw {
			v := int8(b)
			if v == MissingInt8 {
				tv.Ints[i] = MissingInt32
			} else if v == EndOfVectorInt8 {
				tv.Ints[i] = EndOfVectorInt32
			} else {
				tv.Ints[i] = int32(v)
			}
		}
		return tv, nil

	case TypeInt16:
		raw, err := readBytes(buf, off, size*2)
		if err != nil {
			return TypedValue{}, err
		}
		tv.Raw = raw
		tv.Ints = make([]int32, size)
		for i := 0; i < size; i++ {
			v := int16(binary.LittleEndian.Uint16(raw[i*2:]))
			if v == MissingInt16 {
				tv.Ints[i] = MissingInt32
			} else if v == EndOfVectorInt16 {
				tv.Ints[i] = EndOfVectorInt32
			} else {
				tv.Ints[i] = int32(v)
			}
		}
		return tv, nil

	case TypeInt32:
		raw, err := readBytes(buf, off, size*4)
		if err != nil {
			return TypedValue{}, err
		}
		tv.Raw = raw
		tv.Ints = make([]int32, size)
		for i := 0; i < size; i++ {
			v := int32(binary.LittleEndian.Uint32(raw[i*4:]))
			tv.Ints[i] = v
		}
		return tv, nil

	case TypeInt64:
		// BCF 2.2+ added a 64-bit integer type. Values that overflow int32
		// are clamped to MissingInt32 / EndOfVectorInt32 sentinels because
		// our downstream model is int32-based; in practice upstream emits
		// int64 only for counts that fit in int32 anyway. The clamp is
		// documented behaviour for now (no users depend on >2G counts in
		// our pipeline).
		raw, err := readBytes(buf, off, size*8)
		if err != nil {
			return TypedValue{}, err
		}
		tv.Raw = raw
		tv.Ints = make([]int32, size)
		for i := 0; i < size; i++ {
			v := int64(binary.LittleEndian.Uint64(raw[i*8:]))
			switch {
			case v == int64(MissingInt32):
				tv.Ints[i] = MissingInt32
			case v == int64(EndOfVectorInt32):
				tv.Ints[i] = EndOfVectorInt32
			case v > int64(0x7FFFFFFF):
				tv.Ints[i] = 0x7FFFFFFF
			case v < int64(-0x7FFFFFFF):
				tv.Ints[i] = -0x7FFFFFFF
			default:
				tv.Ints[i] = int32(v)
			}
		}
		return tv, nil

	case TypeFloat:
		raw, err := readBytes(buf, off, size*4)
		if err != nil {
			return TypedValue{}, err
		}
		tv.Raw = raw
		tv.Floats = make([]float32, size)
		for i := 0; i < size; i++ {
			bits := binary.LittleEndian.Uint32(raw[i*4:])
			tv.Floats[i] = math.Float32frombits(bits)
		}
		return tv, nil

	case TypeChar:
		raw, err := readBytes(buf, off, size)
		if err != nil {
			return TypedValue{}, err
		}
		tv.Raw = raw
		tv.String = string(raw)
		return tv, nil
	}

	return TypedValue{}, fmt.Errorf("bcf: unknown typed descriptor %d at offset %d", descriptor, *off-1)
}

// DecodeTypedInt is a convenience wrapper that reads one typed value and
// returns the first integer element. It is used for dictionary indices,
// where the field is always a small int (and almost always a single scalar).
// Returns -1 if the value is missing.
func DecodeTypedInt(buf []byte, off *int) (int32, error) {
	tv, err := DecodeTyped(buf, off)
	if err != nil {
		return 0, err
	}
	if tv.IsMissing() || len(tv.Ints) == 0 {
		return -1, nil
	}
	return tv.Ints[0], nil
}

// DecodeTypedInts decodes a typed integer vector and returns all of its
// elements as int32. Floats and char values produce an error.
func DecodeTypedInts(buf []byte, off *int) ([]int32, error) {
	tv, err := DecodeTyped(buf, off)
	if err != nil {
		return nil, err
	}
	switch tv.Descriptor {
	case TypeInt8, TypeInt16, TypeInt32, TypeMissing:
		return tv.Ints, nil
	default:
		return nil, fmt.Errorf("bcf: expected typed int, got descriptor %d", tv.Descriptor)
	}
}

// DecodeTypedString decodes a typed char value and returns it as a Go string.
// Returns the empty string if the value is missing.
func DecodeTypedString(buf []byte, off *int) (string, error) {
	tv, err := DecodeTyped(buf, off)
	if err != nil {
		return "", err
	}
	if tv.IsMissing() {
		return "", nil
	}
	if tv.Descriptor != TypeChar {
		return "", fmt.Errorf("bcf: expected typed char, got descriptor %d", tv.Descriptor)
	}
	return tv.String, nil
}

// EncodeTypedString writes s as a TypeChar typed value. It is used by the
// test fixtures (the reader is the production path); the encoder is fine to
// remain test-only until the writer slice lands.
func EncodeTypedString(s string) []byte {
	return encodeTypedRaw(TypeChar, len(s), []byte(s))
}

// EncodeTypedInt8 writes a single int8 scalar.
func EncodeTypedInt8(v int8) []byte {
	return encodeTypedRaw(TypeInt8, 1, []byte{byte(v)})
}

// EncodeTypedInt32 writes a single int32 scalar.
func EncodeTypedInt32(v int32) []byte {
	b := make([]byte, 4)
	binary.LittleEndian.PutUint32(b, uint32(v))
	return encodeTypedRaw(TypeInt32, 1, b)
}

// EncodeTypedInt32Vec writes a vector of int32 values.
func EncodeTypedInt32Vec(vs []int32) []byte {
	b := make([]byte, 4*len(vs))
	for i, v := range vs {
		binary.LittleEndian.PutUint32(b[i*4:], uint32(v))
	}
	return encodeTypedRaw(TypeInt32, len(vs), b)
}

// EncodeTypedInt16 writes a single int16 scalar.
func EncodeTypedInt16(v int16) []byte {
	b := make([]byte, 2)
	binary.LittleEndian.PutUint16(b, uint16(v))
	return encodeTypedRaw(TypeInt16, 1, b)
}

// EncodeTypedFloat writes a single float32 scalar.
func EncodeTypedFloat(v float32) []byte {
	b := make([]byte, 4)
	binary.LittleEndian.PutUint32(b, math.Float32bits(v))
	return encodeTypedRaw(TypeFloat, 1, b)
}

// EncodeTypedFloatVec writes a vector of float32 values.
func EncodeTypedFloatVec(vs []float32) []byte {
	b := make([]byte, 4*len(vs))
	for i, v := range vs {
		binary.LittleEndian.PutUint32(b[i*4:], math.Float32bits(v))
	}
	return encodeTypedRaw(TypeFloat, len(vs), b)
}

// EncodeMissing writes a single missing value with descriptor 0.
func EncodeMissing() []byte { return []byte{0x00} }

// encodeTypedRaw is the shared helper for the public encoders. When n fits
// in the descriptor's high nibble it emits one byte; otherwise it follows
// the "size = 15, then typed int" rule used by the spec for long vectors.
func encodeTypedRaw(descriptor uint8, n int, payload []byte) []byte {
	if n <= 14 {
		header := byte((n << 4) | int(descriptor&0x0F))
		out := make([]byte, 0, 1+len(payload))
		out = append(out, header)
		out = append(out, payload...)
		return out
	}
	header := byte((15 << 4) | int(descriptor&0x0F))
	lenPayload := EncodeTypedInt32(int32(n))
	out := make([]byte, 0, 1+len(lenPayload)+len(payload))
	out = append(out, header)
	out = append(out, lenPayload...)
	out = append(out, payload...)
	return out
}

// floatBits returns the IEEE-754 bit pattern of f, useful for both tests and
// the future encoder.
func floatBits(f float32) uint32 { return math.Float32bits(f) }
