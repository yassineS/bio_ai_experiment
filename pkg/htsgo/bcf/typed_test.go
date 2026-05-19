package bcf

import (
	"encoding/binary"
	"math"
	"reflect"
	"testing"
)

func TestEncodeDecodeTypedInt8(t *testing.T) {
	buf := EncodeTypedInt8(7)
	if buf[0] != 0x11 {
		t.Fatalf("descriptor: got 0x%02x, want 0x11 (size=1, type=int8)", buf[0])
	}
	off := 0
	tv, err := DecodeTyped(buf, &off)
	if err != nil {
		t.Fatal(err)
	}
	if tv.Descriptor != TypeInt8 || tv.Length != 1 || len(tv.Ints) != 1 || tv.Ints[0] != 7 {
		t.Fatalf("got %+v", tv)
	}
	if off != len(buf) {
		t.Fatalf("off=%d, want %d", off, len(buf))
	}
}

func TestEncodeDecodeTypedInt16(t *testing.T) {
	buf := EncodeTypedInt16(-300)
	off := 0
	tv, err := DecodeTyped(buf, &off)
	if err != nil {
		t.Fatal(err)
	}
	if tv.Descriptor != TypeInt16 || tv.Ints[0] != -300 {
		t.Fatalf("got %+v", tv)
	}
}

func TestEncodeDecodeTypedInt32Vec(t *testing.T) {
	vs := []int32{10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160}
	buf := EncodeTypedInt32Vec(vs)
	// Length is 16 (>14) so descriptor must be 0xF3 followed by a typed int.
	if buf[0] != 0xF3 {
		t.Fatalf("descriptor: got 0x%02x, want 0xF3 (size=15, type=int32)", buf[0])
	}
	off := 0
	tv, err := DecodeTyped(buf, &off)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(tv.Ints, vs) {
		t.Fatalf("ints: got %v, want %v", tv.Ints, vs)
	}
}

func TestEncodeDecodeTypedFloat(t *testing.T) {
	buf := EncodeTypedFloat(3.5)
	off := 0
	tv, err := DecodeTyped(buf, &off)
	if err != nil {
		t.Fatal(err)
	}
	if tv.Descriptor != TypeFloat || tv.Floats[0] != 3.5 {
		t.Fatalf("got %+v", tv)
	}
}

func TestEncodeDecodeTypedString(t *testing.T) {
	buf := EncodeTypedString("rs123")
	off := 0
	tv, err := DecodeTyped(buf, &off)
	if err != nil {
		t.Fatal(err)
	}
	if tv.String != "rs123" {
		t.Fatalf("got %q", tv.String)
	}
	if tv.Length != 5 {
		t.Fatalf("length: got %d, want 5", tv.Length)
	}
}

func TestEncodeDecodeTypedStringLong(t *testing.T) {
	long := "abcdefghijklmnopqrst" // 20 > 14, forces the length-prefix variant
	buf := EncodeTypedString(long)
	if buf[0] != 0xF7 {
		t.Fatalf("descriptor: got 0x%02x, want 0xF7", buf[0])
	}
	off := 0
	tv, err := DecodeTyped(buf, &off)
	if err != nil {
		t.Fatal(err)
	}
	if tv.String != long {
		t.Fatalf("got %q", tv.String)
	}
}

func TestDecodeMissingInt(t *testing.T) {
	// One int8 with the missing sentinel.
	buf := EncodeTypedInt8(MissingInt8)
	off := 0
	tv, err := DecodeTyped(buf, &off)
	if err != nil {
		t.Fatal(err)
	}
	if tv.Ints[0] != MissingInt32 {
		t.Fatalf("missing int8 did not promote to MissingInt32: %d", tv.Ints[0])
	}
}

func TestDecodeMissingFloat(t *testing.T) {
	buf := []byte{0x15} // size=1, type=float
	mb := make([]byte, 4)
	binary.LittleEndian.PutUint32(mb, MissingFloat32)
	buf = append(buf, mb...)
	off := 0
	tv, err := DecodeTyped(buf, &off)
	if err != nil {
		t.Fatal(err)
	}
	if !IsMissingFloat(math.Float32bits(tv.Floats[0])) {
		t.Fatalf("expected missing float, got bits 0x%08x", math.Float32bits(tv.Floats[0]))
	}
}

func TestDecodeMissingDescriptor(t *testing.T) {
	// 0x00 means "no value follows" (descriptor 0 / size 0).
	buf := []byte{0x00}
	off := 0
	tv, err := DecodeTyped(buf, &off)
	if err != nil {
		t.Fatal(err)
	}
	if !tv.IsMissing() {
		t.Fatalf("expected missing, got %+v", tv)
	}
}

func TestDecodeTypedUnknownDescriptor(t *testing.T) {
	buf := []byte{0x14} // descriptor 4 (reserved/unknown) with size 1
	off := 0
	if _, err := DecodeTyped(buf, &off); err == nil {
		t.Fatal("expected error for unknown descriptor")
	}
}

func TestDecodeTypedShortInput(t *testing.T) {
	off := 0
	if _, err := DecodeTyped(nil, &off); err == nil {
		t.Fatal("expected error reading from empty buffer")
	}
	buf := []byte{0x21} // size=2, type=int16, but no payload
	off = 0
	if _, err := DecodeTyped(buf, &off); err == nil {
		t.Fatal("expected error for truncated typed value")
	}
}

func TestDecodeTypedIntsString(t *testing.T) {
	buf := EncodeTypedString("abc")
	off := 0
	if _, err := DecodeTypedInts(buf, &off); err == nil {
		t.Fatal("expected error reading string as ints")
	}
}

func TestDecodeTypedStringWrongType(t *testing.T) {
	buf := EncodeTypedInt8(1)
	off := 0
	if _, err := DecodeTypedString(buf, &off); err == nil {
		t.Fatal("expected error reading int as string")
	}
}

func TestDecodeTypedIntMissing(t *testing.T) {
	// EncodeMissing => single byte 0x00 (descriptor 0 / size 0).
	buf := EncodeMissing()
	off := 0
	v, err := DecodeTypedInt(buf, &off)
	if err != nil {
		t.Fatal(err)
	}
	if v != -1 {
		t.Fatalf("missing descriptor should decode to -1, got %d", v)
	}
}

func TestDecodeTypedStringMissing(t *testing.T) {
	buf := EncodeMissing()
	off := 0
	s, err := DecodeTypedString(buf, &off)
	if err != nil {
		t.Fatal(err)
	}
	if s != "" {
		t.Fatalf("missing descriptor should decode to empty string, got %q", s)
	}
}

func TestFloatBits(t *testing.T) {
	if floatBits(0) != 0 {
		t.Fatalf("0 should have bits 0")
	}
}
