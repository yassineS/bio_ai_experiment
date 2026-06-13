package cram

import (
	"testing"
)

// encEncodingV4 builds the on-wire bytes of a CRAM v4 encoding: uint7
// codec id, uint7 parameter length, then the parameter bytes (which the
// caller has already uint7-encoded). It is the v4 analogue of encEncoding.
func encEncodingV4(id EncodingID, params []byte) []byte {
	out := appendUint7(nil, uint64(uint32(int32(id))))
	out = appendUint7(out, uint64(len(params)))
	return append(out, params...)
}

// TestParseEncodingV4VarintUnsigned parses a v4 VARINT_UNSIGNED encoding
// (content id then a signed 64-bit offset) and decodes values through it.
func TestParseEncodingV4VarintUnsigned(t *testing.T) {
	// content id 11, offset -2.
	params := appendUint7(nil, 11)
	params = appendSint7(params, -2)
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingVarintUnsigned, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding VARINT_UNSIGNED: %v", err)
	}
	if enc.ID != EncodingVarintUnsigned || enc.ExternalID != 11 || enc.VarintOffset != -2 {
		t.Fatalf("VARINT_UNSIGNED parsed wrong: %+v", enc)
	}
	if enc.major != 4 {
		t.Errorf("encoding major = %d, want 4", enc.major)
	}

	// The block holds uint7 values 5 and 100; with offset -2 they decode
	// to 3 and 98.
	block := appendUint7(nil, 5)
	block = appendUint7(block, 100)
	s := &seriesSource{
		core:     newBitReader(nil),
		external: map[int32]*byteCursor{},
		blocks:   map[int32][]byte{11: block},
		reader:   newIntReader(4),
	}
	got, err := enc.decodeInts(s, 2)
	if err != nil {
		t.Fatalf("decodeInts: %v", err)
	}
	if got[0] != 3 || got[1] != 98 {
		t.Errorf("decoded %v, want [3 98]", got)
	}
}

// TestParseEncodingV4VarintSigned decodes through a v4 VARINT_SIGNED
// encoding, whose external block stores zig-zag values.
func TestParseEncodingV4VarintSigned(t *testing.T) {
	params := appendUint7(nil, 7) // content id 7
	params = appendSint7(params, 0)
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingVarintSigned, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding VARINT_SIGNED: %v", err)
	}
	if enc.ID != EncodingVarintSigned || enc.ExternalID != 7 {
		t.Fatalf("VARINT_SIGNED parsed wrong: %+v", enc)
	}
	// zig-zag values for -1, 1, -100.
	block := appendSint7(nil, -1)
	block = appendSint7(block, 1)
	block = appendSint7(block, -100)
	s := &seriesSource{
		core:     newBitReader(nil),
		external: map[int32]*byteCursor{},
		blocks:   map[int32][]byte{7: block},
		reader:   newIntReader(4),
	}
	got, err := enc.decodeInts(s, 3)
	if err != nil {
		t.Fatalf("decodeInts: %v", err)
	}
	if got[0] != -1 || got[1] != 1 || got[2] != -100 {
		t.Errorf("decoded %v, want [-1 1 -100]", got)
	}
}

// TestParseEncodingV4Const parses CONST_INT and CONST_BYTE encodings and
// confirms each decodes to its constant without reading any block.
func TestParseEncodingV4Const(t *testing.T) {
	encInt, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingConstInt, appendSint7(nil, 60)), 0)
	if err != nil {
		t.Fatalf("parseEncoding CONST_INT: %v", err)
	}
	if encInt.ID != EncodingConstInt || encInt.ConstValue != 60 {
		t.Fatalf("CONST_INT parsed wrong: %+v", encInt)
	}
	s := &seriesSource{core: newBitReader(nil), external: map[int32]*byteCursor{}, blocks: map[int32][]byte{}, reader: newIntReader(4)}
	got, err := encInt.decodeInts(s, 3)
	if err != nil {
		t.Fatalf("decodeInts CONST_INT: %v", err)
	}
	for i, v := range got {
		if v != 60 {
			t.Errorf("CONST_INT value %d = %d, want 60", i, v)
		}
	}

	// CONST_INT with a negative constant (it is stored signed).
	encNeg, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingConstInt, appendSint7(nil, -1)), 0)
	if err != nil {
		t.Fatalf("parseEncoding CONST_INT negative: %v", err)
	}
	if v, derr := encNeg.decodeInt(s); derr != nil || v != -1 {
		t.Errorf("CONST_INT(-1) decoded %d, %v; want -1, nil", v, derr)
	}
}

// TestParseEncodingV4ByteArrayStop confirms the v4 BYTE_ARRAY_STOP content
// id is read as a uint7 varint, not ITF-8.
func TestParseEncodingV4ByteArrayStop(t *testing.T) {
	// stop byte 0x00, content id 200 (which differs between ITF-8 and
	// uint7, so a wrong reader would mis-parse it).
	params := []byte{0x00}
	params = appendUint7(params, 200)
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingByteArrayStop, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding BYTE_ARRAY_STOP: %v", err)
	}
	if enc.ID != EncodingByteArrayStop || enc.StopByte != 0x00 || enc.ExternalID != 200 {
		t.Fatalf("BYTE_ARRAY_STOP parsed wrong: %+v", enc)
	}
}

// TestParseEncodingV4Transform confirms a v4 XPACK transform codec parses
// its nbits, reverse map and wrapped sub-encoding (the full decode is
// exercised in transform_test.go).
func TestParseEncodingV4Transform(t *testing.T) {
	// nbits=2, nval=4, map P,A,C,K, wrapping EXTERNAL content id 5.
	params := appendUint7(nil, 2)
	params = appendUint7(params, 4)
	for _, v := range []byte{'P', 'A', 'C', 'K'} {
		params = appendUint7(params, uint64(v))
	}
	sub := appendUint7(nil, 5) // EXTERNAL content id
	params = append(params, encEncodingV4(EncodingExternal, sub)...)

	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingXPack, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding XPACK: %v", err)
	}
	if enc.ID != EncodingXPack || enc.PackBits != 2 || len(enc.PackMap) != 4 {
		t.Fatalf("XPACK parsed wrong: %+v", enc)
	}
	if string(enc.PackMap) != "PACK" {
		t.Errorf("XPACK map = %q, want PACK", enc.PackMap)
	}
	if enc.SubEnc == nil || enc.SubEnc.ID != EncodingExternal || enc.SubEnc.ExternalID != 5 {
		t.Errorf("XPACK sub-encoding parsed wrong: %+v", enc.SubEnc)
	}
}
