package codec

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// TestRANS4x16_ComplianceVectors decodes the pre-compressed htscodecs
// r4x16 vectors (order 0 and order 1) and asserts byte-for-byte against
// the expected raw data. This is the compliance oracle: it proves our
// decoder matches the reference C implementation's on-wire format
// exactly.
func TestRANS4x16_ComplianceVectors(t *testing.T) {
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		want, ok := loadCorpus(t, qfile)
		if !ok {
			continue
		}
		for _, order := range []int{0, 1} {
			comp, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "r4x16",
				qfile+"."+itoa(order)))
			if err != nil {
				continue
			}
			ran++
			t.Run(qfile+".o"+itoa(order), func(t *testing.T) {
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
		for _, order := range []int{0, 1} {
			comp, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "r4x16",
				qfile+"."+itoa(order)))
			if err != nil {
				continue
			}
			ran++
			t.Run(qfile+".o"+itoa(order), func(t *testing.T) {
				got, err := RANS4x16Encode(raw, order)
				if err != nil {
					t.Fatalf("encode: %v", err)
				}
				if !bytes.Equal(got, comp) {
					t.Errorf("encoded %d bytes, htscodecs vector is %d; first mismatch at %d",
						len(got), len(comp), firstDiff(got, comp))
				}
			})
		}
	}
	if ran == 0 {
		t.Skip("htscodecs submodule not initialised")
	}
}

// TestRANS4x16_RoundTrip exercises encode→decode as the identity across
// a spread of input shapes and both orders — no external fixtures
// needed.
func TestRANS4x16_RoundTrip(t *testing.T) {
	inputs := map[string][]byte{
		"empty":         {},
		"single":        {'A'},
		"two":           {'A', 'C'},
		"three":         {'A', 'C', 'G'},
		"four":          {'A', 'C', 'G', 'T'},
		"seven":         []byte("ACGTNNN"),
		"eight":         []byte("ACGTNNNN"),
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
		for _, order := range []int{0, 1} {
			t.Run(name+".o"+itoa(order), func(t *testing.T) {
				comp, err := RANS4x16Encode(in, order)
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
}

// TestRANS4x16_DecodeErrors pins the error paths for malformed and
// out-of-scope input.
func TestRANS4x16_DecodeErrors(t *testing.T) {
	cases := []struct {
		name string
		in   []byte
	}{
		{"empty", nil},
		{"pack rejected", []byte{0x80, 0x00}},
		{"rle rejected", []byte{0x40, 0x00}},
		{"stripe rejected", []byte{0x08, 0x00}},
		{"x32 rejected", []byte{0x04, 0x00}},
		{"nosz rejected", []byte{0x10, 0x00}},
		{"truncated size varint", []byte{0x00, 0x80}},
		{"cat payload too short", []byte{0x20, 0x05, 'A', 'B'}},
		{"order-0 payload too short", []byte{0x00, 0x04, 1, 2, 3}},
		{"order-1 payload too short", []byte{0x01, 0x04, 1, 2, 3}},
		// payload byte 0xB0 selects table precision 11 (0xB0>>4),
		// which is neither 10 nor 12 and must be rejected.
		{"order-1 invalid shift 11", append([]byte{0x01, 0x08, 0xB0}, make([]byte, 20)...)},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if _, err := RANS4x16Decode(c.in); err == nil {
				t.Errorf("expected error for %s", c.name)
			}
		})
	}
}

// TestRANS4x16_Order1Shift12 exercises the 12-bit order-1 table
// precision. The auto-tune (ransComputeShift) almost always selects the
// 10-bit precision for synthetic data — the fast_log exponent term
// cancels exactly for near-uniform distributions — so the encoder is
// pinned to shift 12 via compressO1RANS4x16's forceShift seam and the
// resulting stream round-tripped through the public decoder. This is
// the only path that fills the frequency tables to 4096 entries and
// masks the rANS state with 0xFFF.
func TestRANS4x16_Order1Shift12(t *testing.T) {
	inputs := map[string][]byte{
		"two-symbol":    twoSymbol(9000),
		"adjacent-syms": adjacentSymbols(9000),
		"ascii-text":    []byte(repeat("the quick brown fox ", 800)),
		"full-alpha":    fullAlphabet(20000),
	}
	for name, in := range inputs {
		t.Run(name, func(t *testing.T) {
			stream := frameRANS4x16(in, compressO1RANS4x16(in, tfShiftO1), 0x01)
			if stream[0] != 0x01 {
				t.Fatalf("expected order-1 format byte, got 0x%02x (X_CAT fallback?)", stream[0])
			}
			_, cp, _ := varGetU32(stream, 1)
			if got := int(stream[cp] >> 4); got != tfShiftO1 {
				t.Fatalf("expected shift %d in stream, got %d", tfShiftO1, got)
			}
			decoded, err := RANS4x16Decode(stream)
			if err != nil {
				t.Fatalf("decode: %v", err)
			}
			if !bytes.Equal(decoded, in) {
				t.Fatalf("shift-12 round-trip mismatch: got %d bytes, want %d (first diff at %d)",
					len(decoded), len(in), firstDiff(decoded, in))
			}
		})
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
