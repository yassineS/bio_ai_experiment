package cram

import (
	"bytes"
	"testing"
)

// newV4Source builds a v4 seriesSource whose external blocks are the given
// content-id to payload map. It is the harness for the transform-codec unit
// vectors, which feed a hand-built sub-codec block and check the expanded
// output against a hand-computed expectation.
func newV4Source(blocks map[int32][]byte) *seriesSource {
	return &seriesSource{
		core:     newBitReader(nil),
		external: map[int32]*byteCursor{},
		blocks:   blocks,
		reader:   newIntReader(4),
	}
}

// extSub builds the on-wire bytes of an EXTERNAL sub-encoding (id, length,
// content-id param) as a transform codec stores its wrapped sub-codec.
func extSub(contentID uint64) []byte {
	return encEncodingV4(EncodingExternal, appendUint7(nil, contentID))
}

// xpackParams builds the parameter bytes of an XPACK encoding: nbits, nval,
// the reverse map, then a wrapped EXTERNAL sub-encoding.
func xpackParams(nbits int, mapBytes []byte, subID uint64) []byte {
	p := appendUint7(nil, uint64(nbits))
	p = appendUint7(p, uint64(len(mapBytes)))
	for _, v := range mapBytes {
		p = appendUint7(p, uint64(v))
	}
	return append(p, extSub(subID)...)
}

// TestXPackDecode2Bit checks the XPACK reverse of htscodecs hts_unpack for a
// 2-bit (4-symbol) alphabet: each source byte expands LSB-first into four
// symbols indexed through the reverse map "PACK".
func TestXPackDecode2Bit(t *testing.T) {
	// Map index 0->P, 1->A, 2->C, 3->K.
	// Packed byte holds symbols LSB-first: P(0) A(1) C(2) K(3)
	//   = 0 | 1<<2 | 2<<4 | 3<<6 = 0xE4 -> "PACK".
	// Second byte: A A P P = 1 | 1<<2 | 0 | 0 = 0x05 -> "AAPP".
	block := []byte{0xE4, 0x05}
	params := xpackParams(2, []byte("PACK"), 9)
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingXPack, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding XPACK: %v", err)
	}
	s := newV4Source(map[int32][]byte{9: block})
	got, err := enc.readTransform(s, 8)
	if err != nil {
		t.Fatalf("readTransform: %v", err)
	}
	if want := "PACKAAPP"; string(got) != want {
		t.Errorf("XPACK decoded %q, want %q", got, want)
	}
}

// TestXPackDecode4Bit checks a 4-bit (2-symbol-per-byte) XPACK alphabet.
func TestXPackDecode4Bit(t *testing.T) {
	// nbits=4: each byte holds two nibbles LSB-first. map index i -> 'a'+i.
	mapBytes := []byte{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p'}
	// byte 0x10 -> low nibble 0 ('a'), high nibble 1 ('b') -> "ab".
	// byte 0xFE -> low nibble 14 ('o'), high nibble 15 ('p') -> "op".
	block := []byte{0x10, 0xFE}
	params := xpackParams(4, mapBytes, 3)
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingXPack, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding XPACK 4-bit: %v", err)
	}
	s := newV4Source(map[int32][]byte{3: block})
	got, err := enc.readTransform(s, 4)
	if err != nil {
		t.Fatalf("readTransform: %v", err)
	}
	if want := "abop"; string(got) != want {
		t.Errorf("XPACK 4-bit decoded %q, want %q", got, want)
	}
}

// TestXPackDecode1Bit checks the 1-bit (8-symbol-per-byte) XPACK case.
func TestXPackDecode1Bit(t *testing.T) {
	// nbits=1: 8 bits LSB-first. map 0->'.', 1->'#'.
	// byte 0b10110001 = 0xB1 -> bits LSB-first 1,0,0,0,1,1,0,1
	//   -> # . . . # # . # = "#...##.#".
	block := []byte{0xB1}
	params := xpackParams(1, []byte{'.', '#'}, 7)
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingXPack, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding XPACK 1-bit: %v", err)
	}
	s := newV4Source(map[int32][]byte{7: block})
	got, err := enc.readTransform(s, 8)
	if err != nil {
		t.Fatalf("readTransform: %v", err)
	}
	if want := "#...##.#"; string(got) != want {
		t.Errorf("XPACK 1-bit decoded %q, want %q", got, want)
	}
}

// TestXPackDecodeConstAlphabet checks the nbits==0 degenerate case, where
// every output byte is map[0].
func TestXPackDecodeConstAlphabet(t *testing.T) {
	// nbits=0, single-symbol alphabet 'Z'. Source has 5 placeholder bytes;
	// each yields one 'Z'.
	block := []byte{0, 0, 0, 0, 0}
	params := xpackParams(0, []byte{'Z'}, 2)
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingXPack, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding XPACK nbits=0: %v", err)
	}
	s := newV4Source(map[int32][]byte{2: block})
	got, err := enc.readTransform(s, 5)
	if err != nil {
		t.Fatalf("readTransform: %v", err)
	}
	if want := "ZZZZZ"; string(got) != want {
		t.Errorf("XPACK const decoded %q, want %q", got, want)
	}
}

// xrleParams builds the parameter bytes of an XRLE encoding: the run-symbol
// set, then the length and literal sub-encodings (each EXTERNAL).
func xrleParams(rleSyms []byte, lenID, litID uint64) []byte {
	p := appendUint7(nil, uint64(len(rleSyms)))
	for _, v := range rleSyms {
		p = appendUint7(p, uint64(v))
	}
	p = append(p, extSub(lenID)...)
	return append(p, extSub(litID)...)
}

// TestXRLEDecode checks the XRLE reverse of htscodecs hts_rle_decode: a
// literal whose byte is an RLE symbol consumes one run length L and expands
// to L+1 copies; any other literal stands for a single byte.
func TestXRLEDecode(t *testing.T) {
	// RLE symbol set {'A'}. Literals "ABA":
	//   'A' (RLE) consumes run length 2 -> "AAA"
	//   'B' (literal)                   -> "B"
	//   'A' (RLE) consumes run length 0 -> "A"
	// -> "AAABA", 5 bytes.
	lit := []byte("ABA")
	lenBlock := appendUint7(nil, 5) // total output size
	lenBlock = appendUint7(lenBlock, 2)
	lenBlock = appendUint7(lenBlock, 0)

	params := xrleParams([]byte{'A'}, 10, 11)
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingXRLE, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding XRLE: %v", err)
	}
	if len(enc.RLESyms) != 1 || enc.RLESyms[0] != 'A' {
		t.Fatalf("XRLE syms parsed wrong: %v", enc.RLESyms)
	}
	s := newV4Source(map[int32][]byte{10: lenBlock, 11: lit})
	got, err := enc.readTransform(s, 5)
	if err != nil {
		t.Fatalf("readTransform: %v", err)
	}
	if want := "AAABA"; string(got) != want {
		t.Errorf("XRLE decoded %q, want %q", got, want)
	}
}

// TestXRLEDecodeNoRuns checks an XRLE stream where no literal is an RLE
// symbol: every literal is copied verbatim and the length stream holds only
// the output size.
func TestXRLEDecodeNoRuns(t *testing.T) {
	lit := []byte("hello")
	lenBlock := appendUint7(nil, 5)
	params := xrleParams([]byte{'X'}, 1, 2) // 'X' never appears
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingXRLE, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding XRLE: %v", err)
	}
	s := newV4Source(map[int32][]byte{1: lenBlock, 2: lit})
	got, err := enc.readTransform(s, 5)
	if err != nil {
		t.Fatalf("readTransform: %v", err)
	}
	if want := "hello"; string(got) != want {
		t.Errorf("XRLE no-runs decoded %q, want %q", got, want)
	}
}

// xdeltaParams builds the parameter bytes of an XDELTA encoding: the word
// size, then a wrapped EXTERNAL sub-encoding.
func xdeltaParams(wordSize int, subID uint64) []byte {
	p := appendUint7(nil, uint64(wordSize))
	return append(p, extSub(subID)...)
}

// TestXDeltaDecodeInt checks the XDELTA value-at-a-time integer path
// (cram_xdelta_decode_int): each value is a zig-zag delta from its
// predecessor, prefix-summed onto a running last value.
func TestXDeltaDecodeInt(t *testing.T) {
	// Reconstructed values [5, 3, 8] -> deltas [5, -2, 5] -> zig-zag
	// [10, 3, 10], stored as uint7 in the external block. The EXTERNAL
	// sub-codec reads each as an unsigned uint7 value.
	block := appendUint7(nil, 10)
	block = appendUint7(block, 3)
	block = appendUint7(block, 10)
	params := xdeltaParams(4, 9)
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingXDelta, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding XDELTA: %v", err)
	}
	if enc.DeltaWordSize != 4 || enc.SubEnc == nil {
		t.Fatalf("XDELTA parsed wrong: %+v", enc)
	}
	s := newV4Source(map[int32][]byte{9: block})
	got, err := enc.decodeInts(s, 3)
	if err != nil {
		t.Fatalf("decodeInts: %v", err)
	}
	if want := []int32{5, 3, 8}; !equalInts(got, want) {
		t.Errorf("XDELTA int decoded %v, want %v", got, want)
	}
}

// TestXDeltaDecodeBlock checks the XDELTA 2-byte block path
// (cram_xdelta_decode_block): zig-zag deltas reconstruct little-endian
// int16 values.
func TestXDeltaDecodeBlock(t *testing.T) {
	// Values [5, 3, 8] -> deltas [5, -2, 5] -> zig-zag [10, 3, 10].
	block := appendUint7(nil, 10)
	block = appendUint7(block, 3)
	block = appendUint7(block, 10)
	params := xdeltaParams(2, 6)
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingXDelta, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding XDELTA block: %v", err)
	}
	s := newV4Source(map[int32][]byte{6: block})
	got, err := enc.transformBytes(s)
	if err != nil {
		t.Fatalf("transformBytes: %v", err)
	}
	// int16 little-endian 5, 3, 8.
	want := []byte{5, 0, 3, 0, 8, 0}
	if !bytes.Equal(got, want) {
		t.Errorf("XDELTA block decoded % x, want % x", got, want)
	}
}

// TestXDeltaDecodeBlockNegative checks XDELTA block reconstruction across a
// negative running value.
func TestXDeltaDecodeBlockNegative(t *testing.T) {
	// Values [10, 4, 4, -1] -> deltas [10, -6, 0, -5]
	//   -> zig-zag [20, 11, 0, 9].
	block := appendUint7(nil, 20)
	block = appendUint7(block, 11)
	block = appendUint7(block, 0)
	block = appendUint7(block, 9)
	params := xdeltaParams(2, 1)
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingXDelta, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding XDELTA: %v", err)
	}
	s := newV4Source(map[int32][]byte{1: block})
	got, err := enc.transformBytes(s)
	if err != nil {
		t.Fatalf("transformBytes: %v", err)
	}
	// int16 LE: 10, 4, 4, -1 (0xFFFF).
	want := []byte{10, 0, 4, 0, 4, 0, 0xFF, 0xFF}
	if !bytes.Equal(got, want) {
		t.Errorf("XDELTA block (negative) decoded % x, want % x", got, want)
	}
}

// TestTransformNestedSubCodec checks that a transform whose sub-codec is
// itself a transform expands recursively: an XRLE literal stream that is
// itself produced by XPACK.
func TestTransformNestedSubCodec(t *testing.T) {
	// Inner XPACK: 2-bit map "PACK", source byte 0xE4 -> "PACK".
	innerBlock := []byte{0xE4}
	xpack := xpackParams(2, []byte("PACK"), 20)

	// Outer XRLE wrapping the XPACK as its literal sub-codec. No literal is
	// an RLE symbol, so output == "PACK".
	lenBlock := appendUint7(nil, 4) // output size 4
	outer := appendUint7(nil, 0)    // zero RLE symbols
	outer = append(outer, extSub(21)...)
	outer = append(outer, encEncodingV4(EncodingXPack, xpack)...)

	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingXRLE, outer), 0)
	if err != nil {
		t.Fatalf("parseEncoding nested XRLE/XPACK: %v", err)
	}
	if enc.LitSubEnc == nil || enc.LitSubEnc.ID != EncodingXPack {
		t.Fatalf("XRLE literal sub-encoding is not XPACK: %+v", enc.LitSubEnc)
	}
	s := newV4Source(map[int32][]byte{20: innerBlock, 21: lenBlock})
	got, err := enc.readTransform(s, 4)
	if err != nil {
		t.Fatalf("readTransform nested: %v", err)
	}
	if want := "PACK"; string(got) != want {
		t.Errorf("nested transform decoded %q, want %q", got, want)
	}
}

// TestXPackDrain checks that a transform-wrapped byte series drains its full
// expanded output through DrainSeries' raw-byte path.
func TestXPackDrain(t *testing.T) {
	block := []byte{0xE4} // -> "PACK"
	params := xpackParams(2, []byte("PACK"), 9)
	enc, _, err := parseEncoding(newIntReader(4), encEncodingV4(EncodingXPack, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding XPACK: %v", err)
	}
	s := newV4Source(map[int32][]byte{9: block})
	got, err := enc.drainRawBytes(s)
	if err != nil {
		t.Fatalf("drainRawBytes: %v", err)
	}
	if want := "PACK"; string(got) != want {
		t.Errorf("XPACK drained %q, want %q", got, want)
	}
}

// equalInts reports whether two int32 slices are element-wise equal.
func equalInts(a, b []int32) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
