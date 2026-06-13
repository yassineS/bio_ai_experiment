package cram

import (
	"bytes"
	"testing"
)

// encEncoding builds the on-wire bytes of an encoding: ITF-8 codec id,
// ITF-8 parameter length, then the parameter bytes.
func encEncoding(id EncodingID, params []byte) []byte {
	var b bytes.Buffer
	b.Write(encITF8(int32(id)))
	b.Write(encITF8(int32(len(params))))
	b.Write(params)
	return b.Bytes()
}

// TestParseEncodingNull parses a NULL encoding.
func TestParseEncodingNull(t *testing.T) {
	enc, off, err := parseEncoding(newIntReader(3), encEncoding(EncodingNull, nil), 0)
	if err != nil {
		t.Fatalf("parseEncoding NULL: %v", err)
	}
	if enc.ID != EncodingNull {
		t.Errorf("id = %s, want NULL", enc.ID)
	}
	if off != 2 {
		t.Errorf("offset = %d, want 2", off)
	}
}

// TestParseEncodingExternal parses an EXTERNAL encoding and reads its
// content id.
func TestParseEncodingExternal(t *testing.T) {
	enc, _, err := parseEncoding(newIntReader(3), encEncoding(EncodingExternal, encITF8(37)), 0)
	if err != nil {
		t.Fatalf("parseEncoding EXTERNAL: %v", err)
	}
	if enc.ID != EncodingExternal || enc.ExternalID != 37 {
		t.Errorf("got id=%s extID=%d, want EXTERNAL 37", enc.ID, enc.ExternalID)
	}
}

// TestParseEncodingByteArrayStop parses a BYTE_ARRAY_STOP encoding.
func TestParseEncodingByteArrayStop(t *testing.T) {
	params := append([]byte{0x09}, encITF8(13)...) // stop byte 0x09, content id 13
	enc, _, err := parseEncoding(newIntReader(3), encEncoding(EncodingByteArrayStop, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding BYTE_ARRAY_STOP: %v", err)
	}
	if enc.StopByte != 0x09 || enc.ExternalID != 13 {
		t.Errorf("got stop=%#x id=%d, want 0x09 13", enc.StopByte, enc.ExternalID)
	}
}

// TestParseEncodingByteArrayLen parses a BYTE_ARRAY_LEN encoding with
// two EXTERNAL sub-encodings.
func TestParseEncodingByteArrayLen(t *testing.T) {
	var params bytes.Buffer
	params.Write(encEncoding(EncodingExternal, encITF8(42)))
	params.Write(encEncoding(EncodingExternal, encITF8(37)))
	enc, _, err := parseEncoding(newIntReader(3), encEncoding(EncodingByteArrayLen, params.Bytes()), 0)
	if err != nil {
		t.Fatalf("parseEncoding BYTE_ARRAY_LEN: %v", err)
	}
	if enc.LenEnc == nil || enc.LenEnc.ExternalID != 42 {
		t.Errorf("length sub-encoding wrong: %+v", enc.LenEnc)
	}
	if enc.ValEnc == nil || enc.ValEnc.ExternalID != 37 {
		t.Errorf("values sub-encoding wrong: %+v", enc.ValEnc)
	}
}

// TestParseEncodingBeta parses a BETA encoding and its offset/nbits.
func TestParseEncodingBeta(t *testing.T) {
	params := append(encITF8(5), encITF8(4)...)
	enc, _, err := parseEncoding(newIntReader(3), encEncoding(EncodingBeta, params), 0)
	if err != nil {
		t.Fatalf("parseEncoding BETA: %v", err)
	}
	if enc.Offset != 5 || enc.NumBits != 4 {
		t.Errorf("got offset=%d nbits=%d, want 5 4", enc.Offset, enc.NumBits)
	}
}

// TestParseEncodingBetaBadNumBits rejects a BETA width outside 0..32.
func TestParseEncodingBetaBadNumBits(t *testing.T) {
	params := append(encITF8(0), encITF8(33)...)
	if _, _, err := parseEncoding(newIntReader(3), encEncoding(EncodingBeta, params), 0); err == nil {
		t.Errorf("BETA with 33 bits should be rejected")
	}
}

// TestParseEncodingSubexpGammaGolomb parses the remaining integer codes.
func TestParseEncodingSubexpGammaGolomb(t *testing.T) {
	enc, _, err := parseEncoding(newIntReader(3), encEncoding(EncodingSubexp, append(encITF8(1), encITF8(3)...)), 0)
	if err != nil || enc.Offset != 1 || enc.K != 3 {
		t.Errorf("SUBEXP parse wrong: %+v err=%v", enc, err)
	}
	enc, _, err = parseEncoding(newIntReader(3), encEncoding(EncodingGamma, encITF8(2)), 0)
	if err != nil || enc.Offset != 2 {
		t.Errorf("GAMMA parse wrong: %+v err=%v", enc, err)
	}
	enc, _, err = parseEncoding(newIntReader(3), encEncoding(EncodingGolomb, append(encITF8(0), encITF8(7)...)), 0)
	if err != nil || enc.M != 7 {
		t.Errorf("GOLOMB parse wrong: %+v err=%v", enc, err)
	}
	enc, _, err = parseEncoding(newIntReader(3), encEncoding(EncodingGolombRice, append(encITF8(0), encITF8(4)...)), 0)
	if err != nil || enc.K != 4 {
		t.Errorf("GOLOMB_RICE parse wrong: %+v err=%v", enc, err)
	}
}

// TestParseEncodingGolombBadM rejects a GOLOMB divisor that is not
// positive.
func TestParseEncodingGolombBadM(t *testing.T) {
	if _, _, err := parseEncoding(newIntReader(3), encEncoding(EncodingGolomb, append(encITF8(0), encITF8(0)...)), 0); err == nil {
		t.Errorf("GOLOMB with M=0 should be rejected")
	}
}

// TestParseEncodingHuffman parses a HUFFMAN encoding with a multi-symbol
// alphabet.
func TestParseEncodingHuffman(t *testing.T) {
	var params bytes.Buffer
	params.Write(encITF8(3)) // 3 symbols
	for _, s := range []int32{65, 67, 71} {
		params.Write(encITF8(s))
	}
	params.Write(encITF8(3)) // 3 bit lengths
	for _, l := range []int32{1, 2, 2} {
		params.Write(encITF8(l))
	}
	enc, _, err := parseEncoding(newIntReader(3), encEncoding(EncodingHuffman, params.Bytes()), 0)
	if err != nil {
		t.Fatalf("parseEncoding HUFFMAN: %v", err)
	}
	if len(enc.Symbols) != 3 || len(enc.BitLengths) != 3 {
		t.Fatalf("HUFFMAN sizes wrong: %+v", enc)
	}
	if enc.Symbols[0] != 65 || enc.BitLengths[2] != 2 {
		t.Errorf("HUFFMAN values wrong: %+v", enc)
	}
}

// TestParseEncodingHuffmanMismatch rejects a HUFFMAN whose symbol and
// bit-length counts differ.
func TestParseEncodingHuffmanMismatch(t *testing.T) {
	var params bytes.Buffer
	params.Write(encITF8(2))
	params.Write(encITF8(1))
	params.Write(encITF8(2))
	params.Write(encITF8(1)) // 1 bit length, not 2
	params.Write(encITF8(3))
	if _, _, err := parseEncoding(newIntReader(3), encEncoding(EncodingHuffman, params.Bytes()), 0); err == nil {
		t.Errorf("HUFFMAN with mismatched counts should be rejected")
	}
}

// TestParseEncodingUnknown rejects an unknown codec id.
func TestParseEncodingUnknown(t *testing.T) {
	if _, _, err := parseEncoding(newIntReader(3), encEncoding(EncodingID(99), nil), 0); err == nil {
		t.Errorf("unknown encoding id should be rejected")
	}
}

// TestParseEncodingTruncated checks parseEncoding errors — never panics —
// on input truncated at every byte offset.
func TestParseEncodingTruncated(t *testing.T) {
	var params bytes.Buffer
	params.Write(encEncoding(EncodingExternal, encITF8(42)))
	params.Write(encEncoding(EncodingExternal, encITF8(37)))
	full := encEncoding(EncodingByteArrayLen, params.Bytes())
	for n := 0; n < len(full); n++ {
		func() {
			defer func() {
				if r := recover(); r != nil {
					t.Errorf("parseEncoding panicked at truncation %d: %v", n, r)
				}
			}()
			_, _, _ = parseEncoding(newIntReader(3), full[:n], 0)
		}()
	}
}

// TestEncodingIDString checks the encoding-id names.
func TestEncodingIDString(t *testing.T) {
	cases := map[EncodingID]string{
		EncodingNull: "NULL", EncodingExternal: "EXTERNAL",
		EncodingHuffman: "HUFFMAN", EncodingByteArrayLen: "BYTE_ARRAY_LEN",
		EncodingByteArrayStop: "BYTE_ARRAY_STOP", EncodingBeta: "BETA",
		EncodingSubexp: "SUBEXP", EncodingGamma: "GAMMA",
		EncodingGolomb: "GOLOMB", EncodingGolombRice: "GOLOMB_RICE",
	}
	for id, want := range cases {
		if got := id.String(); got != want {
			t.Errorf("EncodingID(%d).String() = %q, want %q", int32(id), got, want)
		}
	}
	if EncodingID(123).String() == "" {
		t.Errorf("unknown encoding id should have a non-empty string")
	}
}
