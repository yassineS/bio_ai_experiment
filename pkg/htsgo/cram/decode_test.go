package cram

import (
	"bytes"
	"testing"
)

// newTestSource builds a seriesSource from a CORE bitstream and a map of
// external block payloads.
func newTestSource(core []byte, blocks map[int32][]byte) *seriesSource {
	if blocks == nil {
		blocks = map[int32][]byte{}
	}
	return &seriesSource{
		core:     newBitReader(core),
		external: make(map[int32]*byteCursor),
		blocks:   blocks,
	}
}

// TestByteCursor exercises the forward-only external-block cursor.
func TestByteCursor(t *testing.T) {
	c := &byteCursor{data: []byte{1, 2, 3, 4, 5}}
	if c.remaining() != 5 {
		t.Fatalf("remaining = %d, want 5", c.remaining())
	}
	b, err := c.readByte()
	if err != nil || b != 1 {
		t.Fatalf("readByte = %d, %v; want 1, nil", b, err)
	}
	n, err := c.readN(3)
	if err != nil || !bytes.Equal(n, []byte{2, 3, 4}) {
		t.Fatalf("readN(3) = %v, %v", n, err)
	}
	if c.exhausted() {
		t.Errorf("cursor should not be exhausted with one byte left")
	}
	if _, err := c.readN(2); err == nil {
		t.Errorf("over-long readN should error")
	}
	if _, err := c.readN(-1); err == nil {
		t.Errorf("negative readN should error")
	}
}

// TestDecodeExternalInts decodes an EXTERNAL integer series.
func TestDecodeExternalInts(t *testing.T) {
	var blk bytes.Buffer
	want := []int32{0, 5, 130, 16384, 7}
	for _, v := range want {
		blk.Write(encITF8(v))
	}
	enc := &Encoding{ID: EncodingExternal, ExternalID: 1}
	s := newTestSource(nil, map[int32][]byte{1: blk.Bytes()})
	got, err := enc.decodeInts(s, len(want))
	if err != nil {
		t.Fatalf("decodeInts: %v", err)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("value %d: got %d want %d", i, got[i], want[i])
		}
	}
}

// TestDecodeExternalMissingBlock checks decoding from an absent external
// block errors cleanly.
func TestDecodeExternalMissingBlock(t *testing.T) {
	enc := &Encoding{ID: EncodingExternal, ExternalID: 99}
	if _, err := enc.decodeInts(newTestSource(nil, nil), 1); err == nil {
		t.Errorf("decoding from a missing external block should error")
	}
}

// TestDrainExternalInts drains an EXTERNAL integer series to exhaustion.
func TestDrainExternalInts(t *testing.T) {
	var blk bytes.Buffer
	want := []int32{1, 2, 3, 4}
	for _, v := range want {
		blk.Write(encITF8(v))
	}
	enc := &Encoding{ID: EncodingExternal, ExternalID: 2}
	got, err := enc.drainInts(newTestSource(nil, map[int32][]byte{2: blk.Bytes()}))
	if err != nil {
		t.Fatalf("drainInts: %v", err)
	}
	if len(got) != 4 {
		t.Fatalf("drained %d values, want 4", len(got))
	}
}

// TestDrainNonExternal checks draining a non-EXTERNAL integer encoding
// is rejected.
func TestDrainNonExternal(t *testing.T) {
	enc := &Encoding{ID: EncodingBeta, NumBits: 4}
	if _, err := enc.drainInts(newTestSource(nil, nil)); err == nil {
		t.Errorf("draining a BETA series should error")
	}
}

// TestDecodeByteArrayStop decodes a BYTE_ARRAY_STOP series.
func TestDecodeByteArrayStop(t *testing.T) {
	// Three NUL-terminated strings.
	blk := []byte("abc\x00de\x00f\x00")
	enc := &Encoding{ID: EncodingByteArrayStop, ExternalID: 3, StopByte: 0}
	got, err := enc.decodeByteArrays(newTestSource(nil, map[int32][]byte{3: blk}), 3)
	if err != nil {
		t.Fatalf("decodeByteArrays: %v", err)
	}
	want := []string{"abc", "de", "f"}
	for i := range want {
		if string(got[i]) != want[i] {
			t.Errorf("value %d: got %q want %q", i, got[i], want[i])
		}
	}
}

// TestDecodeByteArrayStopRunoff checks a missing stop byte errors.
func TestDecodeByteArrayStopRunoff(t *testing.T) {
	enc := &Encoding{ID: EncodingByteArrayStop, ExternalID: 3, StopByte: 0}
	blk := []byte("noterminator")
	if _, err := enc.decodeByteArrays(newTestSource(nil, map[int32][]byte{3: blk}), 1); err == nil {
		t.Errorf("BYTE_ARRAY_STOP without a stop byte should error")
	}
}

// TestDrainByteArrayStop drains a BYTE_ARRAY_STOP series.
func TestDrainByteArrayStop(t *testing.T) {
	blk := []byte("one\x00two\x00")
	enc := &Encoding{ID: EncodingByteArrayStop, ExternalID: 4, StopByte: 0}
	got, err := enc.drainByteArrays(newTestSource(nil, map[int32][]byte{4: blk}))
	if err != nil {
		t.Fatalf("drainByteArrays: %v", err)
	}
	if len(got) != 2 || string(got[0]) != "one" || string(got[1]) != "two" {
		t.Errorf("drained %v, want [one two]", got)
	}
}

// TestDecodeByteArrayLen decodes a BYTE_ARRAY_LEN series whose length
// and values both come from external blocks.
func TestDecodeByteArrayLen(t *testing.T) {
	values := []string{"hello", "hi", ""}
	var lenBlk, valBlk bytes.Buffer
	for _, v := range values {
		lenBlk.Write(encITF8(int32(len(v))))
		valBlk.WriteString(v)
	}
	enc := &Encoding{
		ID:     EncodingByteArrayLen,
		LenEnc: &Encoding{ID: EncodingExternal, ExternalID: 10},
		ValEnc: &Encoding{ID: EncodingExternal, ExternalID: 11},
	}
	s := newTestSource(nil, map[int32][]byte{10: lenBlk.Bytes(), 11: valBlk.Bytes()})
	got, err := enc.decodeByteArrays(s, len(values))
	if err != nil {
		t.Fatalf("decodeByteArrays: %v", err)
	}
	for i := range values {
		if string(got[i]) != values[i] {
			t.Errorf("value %d: got %q want %q", i, got[i], values[i])
		}
	}
}

// TestDrainByteArrayLen drains a BYTE_ARRAY_LEN series.
func TestDrainByteArrayLen(t *testing.T) {
	values := []string{"abc", "d"}
	var lenBlk, valBlk bytes.Buffer
	for _, v := range values {
		lenBlk.Write(encITF8(int32(len(v))))
		valBlk.WriteString(v)
	}
	enc := &Encoding{
		ID:     EncodingByteArrayLen,
		LenEnc: &Encoding{ID: EncodingExternal, ExternalID: 1},
		ValEnc: &Encoding{ID: EncodingExternal, ExternalID: 2},
	}
	s := newTestSource(nil, map[int32][]byte{1: lenBlk.Bytes(), 2: valBlk.Bytes()})
	got, err := enc.drainByteArrayLen(s)
	if err != nil {
		t.Fatalf("drainByteArrayLen: %v", err)
	}
	if len(got) != 2 || string(got[0]) != "abc" || string(got[1]) != "d" {
		t.Errorf("drained %v, want [abc d]", got)
	}
}

// TestDrainByteArrayLenSharedBlock drains a BYTE_ARRAY_LEN series whose
// length and values sub-encodings name the same external block — the
// interleaved length-then-value layout CRAM permits.
func TestDrainByteArrayLenSharedBlock(t *testing.T) {
	values := []string{"foo", "ba", "x"}
	var blk bytes.Buffer
	for _, v := range values {
		blk.Write(encITF8(int32(len(v))))
		blk.WriteString(v)
	}
	enc := &Encoding{
		ID:     EncodingByteArrayLen,
		LenEnc: &Encoding{ID: EncodingExternal, ExternalID: 9},
		ValEnc: &Encoding{ID: EncodingExternal, ExternalID: 9},
	}
	got, err := enc.drainByteArrayLen(newTestSource(nil, map[int32][]byte{9: blk.Bytes()}))
	if err != nil {
		t.Fatalf("drainByteArrayLen (shared block): %v", err)
	}
	if len(got) != 3 || string(got[0]) != "foo" || string(got[2]) != "x" {
		t.Errorf("drained %v, want [foo ba x]", got)
	}
}

// TestDrainByteArrayLenCoreLength rejects a BYTE_ARRAY_LEN whose length
// sub-encoding lives in the CORE bitstream.
func TestDrainByteArrayLenCoreLength(t *testing.T) {
	enc := &Encoding{
		ID:     EncodingByteArrayLen,
		LenEnc: &Encoding{ID: EncodingHuffman, Symbols: []int32{3}, BitLengths: []int32{0}},
		ValEnc: &Encoding{ID: EncodingExternal, ExternalID: 1},
	}
	if _, err := enc.drainByteArrayLen(newTestSource(nil, nil)); err == nil {
		t.Errorf("draining a BYTE_ARRAY_LEN with a CORE length should error")
	}
}

// TestDecodeBeta decodes a BETA integer series from the CORE bitstream.
func TestDecodeBeta(t *testing.T) {
	var w bitWriter
	want := []int32{0, 5, 15, 8}
	for _, v := range want {
		w.writeBits(uint32(v), 4)
	}
	enc := &Encoding{ID: EncodingBeta, NumBits: 4, Offset: 0}
	got, err := enc.decodeInts(newTestSource(w.bytes(), nil), len(want))
	if err != nil {
		t.Fatalf("decodeInts BETA: %v", err)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("BETA value %d: got %d want %d", i, got[i], want[i])
		}
	}
}

// TestDecodeHuffmanInts decodes a HUFFMAN integer series from the CORE
// bitstream.
func TestDecodeHuffmanInts(t *testing.T) {
	enc := &Encoding{
		ID:         EncodingHuffman,
		Symbols:    []int32{10, 20, 30, 40},
		BitLengths: []int32{1, 2, 3, 3},
	}
	var w bitWriter
	w.writeBits(0b0, 1)   // 10
	w.writeBits(0b10, 2)  // 20
	w.writeBits(0b111, 3) // 40
	got, err := enc.decodeInts(newTestSource(w.bytes(), nil), 3)
	if err != nil {
		t.Fatalf("decodeInts HUFFMAN: %v", err)
	}
	want := []int32{10, 20, 40}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("HUFFMAN value %d: got %d want %d", i, got[i], want[i])
		}
	}
}

// TestDecodeRawBytes reads a fixed-length run through an EXTERNAL
// sub-encoding (the values half of a BYTE_ARRAY_LEN decode).
func TestDecodeRawBytes(t *testing.T) {
	enc := &Encoding{ID: EncodingExternal, ExternalID: 7}
	s := newTestSource(nil, map[int32][]byte{7: []byte("0123456789")})
	got, err := enc.decodeRawBytes(s, 4)
	if err != nil || string(got) != "0123" {
		t.Fatalf("decodeRawBytes = %q, %v", got, err)
	}
}

// TestDecodeNull checks the NULL encoding yields nothing and errors when
// asked for a value.
func TestDecodeNull(t *testing.T) {
	var enc *Encoding
	if iv, err := enc.decodeInts(newTestSource(nil, nil), 0); err != nil || len(iv) != 0 {
		t.Errorf("nil encoding decodeInts(0) = %v, %v", iv, err)
	}
	if _, err := enc.decodeInt(newTestSource(nil, nil)); err == nil {
		t.Errorf("nil encoding decodeInt should error")
	}
	null := &Encoding{ID: EncodingNull}
	if _, err := null.decodeByteArray(newTestSource(nil, nil)); err == nil {
		t.Errorf("NULL encoding decodeByteArray should error")
	}
}

// TestDecodeNegativeCount rejects a negative value count.
func TestDecodeNegativeCount(t *testing.T) {
	enc := &Encoding{ID: EncodingExternal, ExternalID: 1}
	if _, err := enc.decodeInts(newTestSource(nil, nil), -1); err == nil {
		t.Errorf("negative count should error")
	}
	if _, err := enc.decodeByteArrays(newTestSource(nil, nil), -1); err == nil {
		t.Errorf("negative count should error")
	}
}

// TestSeriesSourceHasBlock checks the external-block presence helper.
func TestSeriesSourceHasBlock(t *testing.T) {
	s := newTestSource(nil, map[int32][]byte{5: {1, 2, 3}})
	if !s.hasBlock(5) {
		t.Errorf("hasBlock(5) should be true")
	}
	if s.hasBlock(6) {
		t.Errorf("hasBlock(6) should be false")
	}
}

// TestDecodeByteArrayExternal decodes single bytes through a bare
// EXTERNAL byte-array encoding.
func TestDecodeByteArrayExternal(t *testing.T) {
	enc := &Encoding{ID: EncodingExternal, ExternalID: 1}
	s := newTestSource(nil, map[int32][]byte{1: {'x', 'y', 'z'}})
	got, err := enc.decodeByteArrays(s, 3)
	if err != nil {
		t.Fatalf("decodeByteArrays: %v", err)
	}
	if string(got[0]) != "x" || string(got[2]) != "z" {
		t.Errorf("got %v, want [x y z]", got)
	}
}

// TestDecodeRawBytesHuffman reads a fixed-length byte run through a
// HUFFMAN values sub-encoding (each byte is a Huffman symbol).
func TestDecodeRawBytesHuffman(t *testing.T) {
	enc := &Encoding{
		ID:         EncodingHuffman,
		Symbols:    []int32{'A', 'C', 'G', 'T'},
		BitLengths: []int32{1, 2, 3, 3},
	}
	var w bitWriter
	w.writeBits(0b0, 1)   // A
	w.writeBits(0b10, 2)  // C
	w.writeBits(0b110, 3) // G
	got, err := enc.decodeRawBytes(newTestSource(w.bytes(), nil), 3)
	if err != nil {
		t.Fatalf("decodeRawBytes HUFFMAN: %v", err)
	}
	if string(got) != "ACG" {
		t.Errorf("got %q, want ACG", got)
	}
}

// TestDecodeRawBytesStop reads a fixed-length run through a
// BYTE_ARRAY_STOP values sub-encoding.
func TestDecodeRawBytesStop(t *testing.T) {
	enc := &Encoding{ID: EncodingByteArrayStop, ExternalID: 2, StopByte: 0}
	s := newTestSource(nil, map[int32][]byte{2: []byte("abcdef")})
	got, err := enc.decodeRawBytes(s, 4)
	if err != nil || string(got) != "abcd" {
		t.Fatalf("decodeRawBytes STOP = %q, %v", got, err)
	}
}

// TestDecodeRawBytesNull checks NULL supplies zero bytes but errors when
// asked for more.
func TestDecodeRawBytesNull(t *testing.T) {
	null := &Encoding{ID: EncodingNull}
	if b, err := null.decodeRawBytes(newTestSource(nil, nil), 0); err != nil || b != nil {
		t.Errorf("NULL decodeRawBytes(0) = %v, %v", b, err)
	}
	if _, err := null.decodeRawBytes(newTestSource(nil, nil), 3); err == nil {
		t.Errorf("NULL decodeRawBytes(3) should error")
	}
}

// TestDecodeRawBytesUnsupported rejects a values sub-encoding that
// cannot supply raw bytes.
func TestDecodeRawBytesUnsupported(t *testing.T) {
	enc := &Encoding{ID: EncodingGamma}
	if _, err := enc.decodeRawBytes(newTestSource([]byte{0}, nil), 1); err == nil {
		t.Errorf("GAMMA decodeRawBytes should error")
	}
}

// TestDrainRawBytes drains a byte-valued EXTERNAL series.
func TestDrainRawBytes(t *testing.T) {
	enc := &Encoding{ID: EncodingExternal, ExternalID: 3}
	got, err := enc.drainRawBytes(newTestSource(nil, map[int32][]byte{3: []byte("quality")}))
	if err != nil || string(got) != "quality" {
		t.Fatalf("drainRawBytes = %q, %v", got, err)
	}
	// A nil encoding drains to nothing.
	var nilEnc *Encoding
	if b, err := nilEnc.drainRawBytes(newTestSource(nil, nil)); err != nil || b != nil {
		t.Errorf("nil drainRawBytes = %v, %v", b, err)
	}
	// A non-EXTERNAL encoding cannot be raw-drained.
	beta := &Encoding{ID: EncodingBeta}
	if _, err := beta.drainRawBytes(newTestSource(nil, nil)); err == nil {
		t.Errorf("BETA drainRawBytes should error")
	}
}

// TestDrainEmptyEncodings checks the drains tolerate a nil / NULL
// encoding by returning nothing.
func TestDrainEmptyEncodings(t *testing.T) {
	var nilEnc *Encoding
	if v, err := nilEnc.drainInts(newTestSource(nil, nil)); err != nil || v != nil {
		t.Errorf("nil drainInts = %v, %v", v, err)
	}
	if v, err := nilEnc.drainByteArrays(newTestSource(nil, nil)); err != nil || v != nil {
		t.Errorf("nil drainByteArrays = %v, %v", v, err)
	}
	null := &Encoding{ID: EncodingNull}
	if v, err := null.drainByteArrays(newTestSource(nil, nil)); err != nil || v != nil {
		t.Errorf("NULL drainByteArrays = %v, %v", v, err)
	}
}

// TestDrainByteArraysNonStop rejects draining a non-BYTE_ARRAY_STOP
// byte-array encoding.
func TestDrainByteArraysNonStop(t *testing.T) {
	enc := &Encoding{ID: EncodingByteArrayLen}
	if _, err := enc.drainByteArrays(newTestSource(nil, nil)); err == nil {
		t.Errorf("draining a BYTE_ARRAY_LEN via drainByteArrays should error")
	}
}

// TestDecodeByteArrayLenNegativeLength rejects a negative decoded
// length.
func TestDecodeByteArrayLenNegativeLength(t *testing.T) {
	var lenBlk bytes.Buffer
	lenBlk.Write(encITF8(-1))
	enc := &Encoding{
		ID:     EncodingByteArrayLen,
		LenEnc: &Encoding{ID: EncodingExternal, ExternalID: 1},
		ValEnc: &Encoding{ID: EncodingExternal, ExternalID: 2},
	}
	s := newTestSource(nil, map[int32][]byte{1: lenBlk.Bytes(), 2: {}})
	if _, err := enc.decodeByteArrays(s, 1); err == nil {
		t.Errorf("a negative BYTE_ARRAY_LEN length should error")
	}
}

// TestDecodeIntAllEncodings exercises decodeInt for every integer
// encoding through the Encoding dispatch, with each code's parameters.
func TestDecodeIntAllEncodings(t *testing.T) {
	// GAMMA: value 9 with offset 2 -> 7.
	var gw bitWriter
	encGamma(&gw, 9)
	gammaEnc := &Encoding{ID: EncodingGamma, Offset: 2}
	if v, err := gammaEnc.decodeInt(newTestSource(gw.bytes(), nil)); err != nil || v != 7 {
		t.Errorf("GAMMA decodeInt = %d, %v; want 7", v, err)
	}
	// SUBEXP: value 20, k=2, offset 5 -> 15.
	var sw bitWriter
	encSubexp(&sw, 20, 2)
	subEnc := &Encoding{ID: EncodingSubexp, K: 2, Offset: 5}
	if v, err := subEnc.decodeInt(newTestSource(sw.bytes(), nil)); err != nil || v != 15 {
		t.Errorf("SUBEXP decodeInt = %d, %v; want 15", v, err)
	}
	// GOLOMB: value 17, m=5, offset 1 -> 16.
	var gow bitWriter
	encGolomb(&gow, 17, 5)
	golEnc := &Encoding{ID: EncodingGolomb, M: 5, Offset: 1}
	if v, err := golEnc.decodeInt(newTestSource(gow.bytes(), nil)); err != nil || v != 16 {
		t.Errorf("GOLOMB decodeInt = %d, %v; want 16", v, err)
	}
}

// TestDecodeIntTruncated checks decodeInt surfaces a truncated CORE
// stream for every bitstream encoding rather than panicking.
func TestDecodeIntTruncated(t *testing.T) {
	for _, enc := range []*Encoding{
		{ID: EncodingGamma},
		{ID: EncodingSubexp, K: 4},
		{ID: EncodingGolomb, M: 5},
		{ID: EncodingGolombRice, K: 3},
		{ID: EncodingBeta, NumBits: 8},
		{ID: EncodingHuffman, Symbols: []int32{1, 2}, BitLengths: []int32{1, 1}},
	} {
		if _, err := enc.decodeInt(newTestSource(nil, nil)); err == nil {
			t.Errorf("%s decodeInt on an empty CORE stream should error", enc.ID)
		}
	}
}

// TestDecodeWrongKind checks asking an integer encoding for a byte array
// (and vice versa) errors rather than mis-decoding.
func TestDecodeWrongKind(t *testing.T) {
	beta := &Encoding{ID: EncodingBeta, NumBits: 4}
	if _, err := beta.decodeByteArray(newTestSource([]byte{0}, nil)); err == nil {
		t.Errorf("BETA decodeByteArray should error")
	}
	bal := &Encoding{ID: EncodingByteArrayLen}
	if _, err := bal.decodeInt(newTestSource(nil, nil)); err == nil {
		t.Errorf("BYTE_ARRAY_LEN decodeInt should error")
	}
}
