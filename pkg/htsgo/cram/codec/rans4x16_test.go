package codec

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// TestRANS4x16_ComplianceVectors decodes the pre-compressed htscodecs
// r4x16 order-0 vectors and asserts byte-for-byte against the expected
// raw data. This is the compliance oracle: it proves our decoder
// matches the reference C implementation's on-wire format exactly.
func TestRANS4x16_ComplianceVectors(t *testing.T) {
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		want, ok := loadCorpus(t, qfile)
		if !ok {
			continue
		}
		comp, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "r4x16", qfile+".0"))
		if err != nil {
			continue
		}
		ran++
		t.Run(qfile+".o0", func(t *testing.T) {
			got, err := RANS4x16Decode(comp)
			if err != nil {
				t.Fatalf("decode: %v", err)
			}
			if !bytes.Equal(got, want) {
				t.Fatalf("decoded %d bytes, want %d; first mismatch at %d",
					len(got), len(want), firstDiff(got, want))
			}
		})
	}
	if ran == 0 {
		t.Skip("htscodecs submodule not initialised — compliance vectors unavailable")
	}
}

// TestRANS4x16_EncodeMatchesHTScodecs checks our encoder produces
// byte-identical output to the htscodecs reference vectors. rANS
// encoding has no freedom once the frequency-normalisation algorithm is
// fixed, so a byte-exact match is both achievable and the strongest
// possible encoder test.
func TestRANS4x16_EncodeMatchesHTScodecs(t *testing.T) {
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		raw, ok := loadCorpus(t, qfile)
		if !ok {
			continue
		}
		comp, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "r4x16", qfile+".0"))
		if err != nil {
			continue
		}
		ran++
		t.Run(qfile+".o0", func(t *testing.T) {
			got, err := RANS4x16Encode(raw, 0)
			if err != nil {
				t.Fatalf("encode: %v", err)
			}
			if !bytes.Equal(got, comp) {
				t.Errorf("encoded %d bytes, htscodecs vector is %d; first mismatch at %d",
					len(got), len(comp), firstDiff(got, comp))
			}
		})
	}
	if ran == 0 {
		t.Skip("htscodecs submodule not initialised")
	}
}

// TestRANS4x16_RoundTrip exercises encode→decode as the identity across
// a spread of input shapes — no external fixtures needed.
func TestRANS4x16_RoundTrip(t *testing.T) {
	inputs := map[string][]byte{
		"empty":         {},
		"single":        {'A'},
		"two":           {'A', 'C'},
		"three":         {'A', 'C', 'G'},
		"four":          {'A', 'C', 'G', 'T'},
		"seven":         []byte("ACGTNNN"),
		"uniform":       bytes.Repeat([]byte{'N'}, 5000),
		"uniform-zero":  bytes.Repeat([]byte{0}, 5000),
		"two-symbol":    twoSymbol(4000),
		"ascii-text":    []byte("the quick brown fox jumps over the lazy dog, repeatedly. " + repeat("ACGTN", 600)),
		"full-alpha":    fullAlphabet(8192),
		"adjacent-syms": adjacentSymbols(6000),
		"random-small":  randomBytes(t, 37),
		"random-large":  randomBytes(t, 200000),
	}
	for name, in := range inputs {
		t.Run(name, func(t *testing.T) {
			comp, err := RANS4x16Encode(in, 0)
			if err != nil {
				t.Fatalf("encode: %v", err)
			}
			got, err := RANS4x16Decode(comp)
			if err != nil {
				t.Fatalf("decode: %v", err)
			}
			if !bytes.Equal(got, in) {
				t.Fatalf("round-trip mismatch: got %d bytes, want %d (first diff at %d)",
					len(got), len(in), firstDiff(got, in))
			}
		})
	}
}

// TestRANS4x16_DecodeErrors pins the error paths for malformed and
// out-of-scope input.
func TestRANS4x16_DecodeErrors(t *testing.T) {
	cases := []struct {
		name string
		in   []byte
	}{
		{"empty", nil},
		{"order-1 rejected", []byte{0x01, 0x00}},
		{"pack rejected", []byte{0x80, 0x00}},
		{"rle rejected", []byte{0x40, 0x00}},
		{"stripe rejected", []byte{0x08, 0x00}},
		{"x32 rejected", []byte{0x04, 0x00}},
		{"nosz rejected", []byte{0x10, 0x00}},
		{"truncated size varint", []byte{0x00, 0x80}},
		{"cat payload too short", []byte{0x20, 0x05, 'A', 'B'}},
		{"order-0 payload too short", []byte{0x00, 0x04, 1, 2, 3}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if _, err := RANS4x16Decode(c.in); err == nil {
				t.Errorf("expected error for %s", c.name)
			}
		})
	}
}

// TestRANS4x16_EncodeRejectsOrder1 confirms the C2 slice declines the
// order-1 model rather than silently producing an order-0 stream.
func TestRANS4x16_EncodeRejectsOrder1(t *testing.T) {
	if _, err := RANS4x16Encode([]byte("ACGT"), 1); err == nil {
		t.Fatal("expected order-1 encode to be rejected")
	}
}

// adjacentSymbols builds data whose alphabet is a contiguous run of
// byte values, exercising the delta-RLE alphabet encoder/decoder.
func adjacentSymbols(n int) []byte {
	out := make([]byte, n)
	for i := range out {
		out[i] = byte(40 + i%50)
	}
	return out
}
