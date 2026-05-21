package codec

import (
	"bytes"
	"math/rand"
	"testing"
)

// lzmaTestInputs returns a spread of payloads exercising the empty,
// small, large, repetitive and random shapes a CRAM block can carry.
func lzmaTestInputs() map[string][]byte {
	repetitive := bytes.Repeat([]byte("ACGT"), 4096)
	large := make([]byte, 1<<18)
	for i := range large {
		large[i] = byte(i*31 + 7)
	}
	rng := rand.New(rand.NewSource(1))
	random := make([]byte, 1<<16)
	for i := range random {
		random[i] = byte(rng.Intn(256))
	}
	return map[string][]byte{
		"empty":      {},
		"one byte":   {0x42},
		"small text": []byte("the quick brown fox jumps over the lazy dog"),
		"repetitive": repetitive,
		"large":      large,
		"random":     random,
		"binary":     {0, 0, 0, 1, 2, 3, 0xff, 0xfe, 0, 0},
	}
}

// TestLZMARoundTrip encodes each input with LZMAEncode and asserts that
// LZMADecode recovers it byte-for-byte. LZMAEncode emits the same .xz
// container framing that htslib's CRAM writer produces, so a successful
// round-trip is the oracle for the framing: no method-3 LZMA fixture
// exists in the reference_code corpus to test against directly.
func TestLZMARoundTrip(t *testing.T) {
	for name, in := range lzmaTestInputs() {
		t.Run(name, func(t *testing.T) {
			enc, err := LZMAEncode(in)
			if err != nil {
				t.Fatalf("LZMAEncode: %v", err)
			}
			if len(enc) >= 6 && !bytes.Equal(enc[:6], []byte{0xFD, '7', 'z', 'X', 'Z', 0x00}) {
				t.Fatalf("encoded stream does not start with the .xz magic: %x", enc[:6])
			}
			dec, err := LZMADecode(enc)
			if err != nil {
				t.Fatalf("LZMADecode: %v", err)
			}
			if !bytes.Equal(dec, in) {
				t.Fatalf("round-trip mismatch: got %d bytes, want %d", len(dec), len(in))
			}
		})
	}
}

// TestLZMADecodeMalformed feeds corrupt, truncated and non-.xz inputs to
// the decoder. Every case must return an error and never panic.
func TestLZMADecodeMalformed(t *testing.T) {
	valid, err := LZMAEncode([]byte("a representative CRAM block payload"))
	if err != nil {
		t.Fatalf("LZMAEncode: %v", err)
	}
	cases := map[string][]byte{
		"empty":              {},
		"not xz":             []byte("this is plainly not an xz stream"),
		"truncated header":   valid[:3],
		"truncated body":     valid[:len(valid)/2],
		"bad magic":          append([]byte{0x00, 0x00}, valid[2:]...),
		"flipped data byte":  flipByte(valid, len(valid)/2),
		"all zero":           make([]byte, 64),
		"missing stream end": valid[:len(valid)-2],
	}
	for name, in := range cases {
		t.Run(name, func(t *testing.T) {
			defer func() {
				if r := recover(); r != nil {
					t.Fatalf("LZMADecode panicked on %s input: %v", name, r)
				}
			}()
			if _, err := LZMADecode(in); err == nil {
				t.Fatalf("LZMADecode accepted malformed %s input", name)
			}
		})
	}
}

// flipByte returns a copy of in with the byte at i inverted.
func flipByte(in []byte, i int) []byte {
	out := append([]byte(nil), in...)
	if i >= 0 && i < len(out) {
		out[i] ^= 0xff
	}
	return out
}

func TestLZMAEncodeIsDeterministic(t *testing.T) {
	in := []byte(bytes.Repeat([]byte("CRAM"), 1000))
	a, err := LZMAEncode(in)
	if err != nil {
		t.Fatalf("LZMAEncode: %v", err)
	}
	b, err := LZMAEncode(in)
	if err != nil {
		t.Fatalf("LZMAEncode: %v", err)
	}
	if !bytes.Equal(a, b) {
		t.Fatal("LZMAEncode produced differing output for identical input")
	}
}
