package codec

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// arithVectorSuffixes are the htscodecs arith_dynamic vector order bytes
// the vendored corpus ships for the q-files: plain order 0/1, PACK
// (0x80), RLE (0x40), PACK|RLE (0xC0) and STRIPE (0x08), each with and
// without the order-1 bit. TestArith_ComplianceVectors skips any suffix
// the corpus does not ship for a given q-file.
var arithVectorSuffixes = []int{0, 1, 8, 9, 64, 65, 128, 129, 192, 193}

// TestArith_ComplianceVectors decodes the pre-compressed htscodecs
// arith_dynamic q-file vectors and asserts byte-for-byte against the
// expected raw data. This is the compliance oracle: it proves our
// adaptive range-coder decoder matches the reference C on-wire format
// exactly.
func TestArith_ComplianceVectors(t *testing.T) {
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		want, ok := loadCorpus(t, qfile)
		if !ok {
			continue
		}
		for _, order := range arithVectorSuffixes {
			comp, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "arith",
				qfile+"."+itoa(order)))
			if err != nil {
				continue
			}
			ran++
			t.Run(qfile+".o"+itoa(order), func(t *testing.T) {
				got, err := ArithDecode(comp)
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

// TestArith_EncodeMatchesHTScodecs checks our encoder produces
// byte-identical output to the htscodecs arith_dynamic q-file vectors.
// The adaptive model is fully deterministic — once the model-update and
// renormalisation rules are fixed there is no freedom — so a byte-exact
// match is the strongest possible encoder test.
func TestArith_EncodeMatchesHTScodecs(t *testing.T) {
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		raw, ok := loadCorpus(t, qfile)
		if !ok {
			continue
		}
		for _, order := range arithVectorSuffixes {
			comp, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "arith",
				qfile+"."+itoa(order)))
			if err != nil {
				continue
			}
			ran++
			t.Run(qfile+".o"+itoa(order), func(t *testing.T) {
				got, err := ArithEncode(raw, order)
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

// TestArith_U32ComplianceVectors decodes the htscodecs arith_dynamic u32
// vectors. The u32 corpus file is fed to the codec verbatim (no q-file
// cut/strip transform). Suffix 4 is X_EXT (a bzip2 payload); decode is
// supported, encode is not (Go has no stdlib bzip2 encoder).
func TestArith_U32ComplianceVectors(t *testing.T) {
	raw, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "u32"))
	if err != nil {
		t.Skip("htscodecs submodule not initialised — u32 corpus unavailable")
	}
	ran := 0
	for _, order := range []int{1, 4, 9, 65} {
		comp, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "arith",
			"u32."+itoa(order)))
		if err != nil {
			continue
		}
		ran++
		t.Run("u32.o"+itoa(order), func(t *testing.T) {
			got, err := ArithDecode(comp)
			if err != nil {
				t.Fatalf("decode: %v", err)
			}
			if !bytes.Equal(got, raw) {
				t.Fatalf("decoded %d bytes, want %d; first mismatch at %d",
					len(got), len(raw), firstDiff(got, raw))
			}
		})
	}
	if ran == 0 {
		t.Skip("htscodecs submodule not initialised — u32 vectors unavailable")
	}
}

// TestArith_U32EncodeMatchesHTScodecs checks the encoder reproduces the
// non-X_EXT u32 vectors byte-exactly. Suffix 4 (X_EXT) is excluded: it
// needs a bzip2 encoder, which Go's standard library does not provide.
func TestArith_U32EncodeMatchesHTScodecs(t *testing.T) {
	raw, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "u32"))
	if err != nil {
		t.Skip("htscodecs submodule not initialised — u32 corpus unavailable")
	}
	ran := 0
	for _, order := range []int{1, 9, 65} {
		comp, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "arith",
			"u32."+itoa(order)))
		if err != nil {
			continue
		}
		ran++
		t.Run("u32.o"+itoa(order), func(t *testing.T) {
			got, err := ArithEncode(raw, order)
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
		t.Skip("htscodecs submodule not initialised — u32 vectors unavailable")
	}
}

// TestArith_XEXTEncodeUnsupported pins the documented encoder gap:
// X_EXT (bzip2) encode returns a clear error rather than producing a
// wrong stream.
func TestArith_XEXTEncodeUnsupported(t *testing.T) {
	if _, err := ArithEncode(bytes.Repeat([]byte("ACGT"), 100), x4x16Ext); err == nil {
		t.Fatal("expected an error for X_EXT encode")
	}
}

// TestArith_RoundTrip exercises encode→decode as the identity across a
// spread of input shapes and both orders — no external fixtures needed.
func TestArith_RoundTrip(t *testing.T) {
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
				comp, err := ArithEncode(in, order)
				if err != nil {
					t.Fatalf("encode: %v", err)
				}
				got, err := ArithDecode(comp)
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

// TestArith_TransformRoundTrip exercises encode→decode as the identity
// for every PACK/RLE/STRIPE transform and order combination.
func TestArith_TransformRoundTrip(t *testing.T) {
	inputs := map[string][]byte{
		"empty":         {},
		"single":        {'A'},
		"tiny":          []byte("ACGT"),
		"twentyone":     bytes.Repeat([]byte("ACGTACGT_"), 3)[:21],
		"low-alpha":     []byte("ACGTACGTACGTNNNNACGTACGTNNNN"),
		"runs":          bytes.Repeat([]byte{'A'}, 4000),
		"two-symbol":    twoSymbol(8000),
		"adjacent-syms": adjacentSymbols(8000),
		"sixteen-sym":   sixteenSymbol(9000),
		"ascii-text":    []byte(repeat("the quick brown fox jumps over the lazy dog. ", 400)),
		"full-alpha":    fullAlphabet(60000),
		"random-large":  randomBytes(t, 120000),
		// The final byte 'Z' appears nowhere else: a valid symbol that is
		// never a context, exercising the order-1 path.
		"last-byte-unique": append(bytes.Repeat([]byte("ABC"), 200), 'Z'),
	}
	for _, transform := range []struct {
		name string
		bits int
	}{
		{"pack", x4x16Pack},
		{"rle", x4x16RLE},
		{"pack+rle", x4x16Pack | x4x16RLE},
		{"stripe", x4x16Stripe},
	} {
		for name, in := range inputs {
			for _, order := range []int{0, 1} {
				t.Run(transform.name+"/"+name+".o"+itoa(order), func(t *testing.T) {
					comp, err := ArithEncode(in, transform.bits|order)
					if err != nil {
						t.Fatalf("encode: %v", err)
					}
					got, err := ArithDecode(comp)
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
}

// TestArith_StripeStreamCount checks that an explicit stripe count N,
// passed in bits 8-15 of the order argument, round-trips for N from 1 to
// 8.
func TestArith_StripeStreamCount(t *testing.T) {
	in := fullAlphabet(40000)
	for n := 1; n <= 8; n++ {
		t.Run("N="+itoa(n), func(t *testing.T) {
			comp, err := ArithEncode(in, x4x16Stripe|(n<<8))
			if err != nil {
				t.Fatalf("encode: %v", err)
			}
			got, err := ArithDecode(comp)
			if err != nil {
				t.Fatalf("decode: %v", err)
			}
			if !bytes.Equal(got, in) {
				t.Fatalf("N=%d round-trip mismatch (first diff at %d)", n, firstDiff(got, in))
			}
		})
	}
}

// TestArith_XCATFallback checks that a transform requested on data the
// range coder cannot shrink still round-trips: the encoder falls back to
// X_CAT for the (transformed) payload and the decoder reverses the
// transform regardless.
func TestArith_XCATFallback(t *testing.T) {
	in := []byte{0x11, 0x22, 0x33, 0x44, 0x55, 0x66}
	for _, bits := range []int{0, x4x16Cat, x4x16Pack, x4x16RLE, x4x16Pack | x4x16RLE} {
		comp, err := ArithEncode(in, bits)
		if err != nil {
			t.Fatalf("encode bits 0x%02x: %v", bits, err)
		}
		got, err := ArithDecode(comp)
		if err != nil {
			t.Fatalf("decode bits 0x%02x: %v", bits, err)
		}
		if !bytes.Equal(got, in) {
			t.Fatalf("bits 0x%02x round-trip mismatch", bits)
		}
	}
}

// TestArith_DecodeErrors pins the error paths for malformed input.
func TestArith_DecodeErrors(t *testing.T) {
	cases := []struct {
		name string
		in   []byte
	}{
		{"empty", nil},
		{"pack meta truncated", []byte{0x80, 0x00}},
		{"stripe truncated", []byte{0x08, 0x00}},
		{"stripe zero streams", []byte{0x08, 0x04, 0x00}},
		{"truncated size varint", []byte{0x00, 0x80}},
		{"cat payload too short", []byte{0x20, 0x05, 'A', 'B'}},
		{"order-0 stream too short", []byte{0x00, 0x04, 0x05, 1, 2}},
		{"order-1 stream too short", []byte{0x01, 0x04, 0x05, 1, 2}},
		{"ext bad bzip2", []byte{0x04, 0x04, 'n', 'o', 't', 'b', 'z'}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if _, err := ArithDecode(c.in); err == nil {
				t.Errorf("expected error for %s", c.name)
			}
		})
	}
}
