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
		t.Fatalf("htscodecs submodule not initialised — compliance vectors unavailable; run `git submodule update --init reference_code/htscodecs`")
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
		t.Fatalf("htscodecs submodule not initialised; run `git submodule update --init reference_code/htscodecs`")
	}
}

// transformVectorSuffixes are the htscodecs r4x16 vector order bytes for
// the PACK/RLE/STRIPE transform layer. The suffix is the order argument
// passed to the encoder: 0x08 STRIPE, 0x40 RLE, 0x80 PACK, optionally
// OR'd with 0x01 for order-1.
var transformVectorSuffixes = []int{8, 9, 64, 65, 128, 129, 192, 193}

// x32VectorSuffixes are the htscodecs r4x16 vector order bytes for the
// 32-way coder (X_32, 0x04). Suffix 4 is X_32 order-0, suffix 5 is X_32
// order-1. The combination suffixes 68/69/132/133/196/197 (X_32 + PACK /
// RLE / PACK|RLE) are included too — TestRANS4x16_X32ComplianceVectors
// skips any that the vendored corpus does not ship.
var x32VectorSuffixes = []int{4, 5, 68, 69, 132, 133, 196, 197}

// TestRANS4x16_X32ComplianceVectors decodes the htscodecs r4x16 32-way
// (X_32) vectors and asserts byte-for-byte against the expected raw
// data. It is the compliance oracle for the C2.3 codec — it proves the
// 32-way decoder matches the reference on-wire format exactly.
func TestRANS4x16_X32ComplianceVectors(t *testing.T) {
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		want, ok := loadCorpus(t, qfile)
		if !ok {
			continue
		}
		for _, order := range x32VectorSuffixes {
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
		t.Fatalf("htscodecs submodule not initialised — X_32 compliance vectors unavailable; run `git submodule update --init reference_code/htscodecs`")
	}
}

// TestRANS4x16_X32EncodeMatchesHTScodecs checks the 32-way encoder
// produces byte-identical output to the htscodecs reference vectors.
// rANS encoding has no freedom once frequency normalisation is fixed, so
// a byte-exact match is the strongest possible encoder test.
func TestRANS4x16_X32EncodeMatchesHTScodecs(t *testing.T) {
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		raw, ok := loadCorpus(t, qfile)
		if !ok {
			continue
		}
		for _, order := range x32VectorSuffixes {
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
		t.Fatalf("htscodecs submodule not initialised; run `git submodule update --init reference_code/htscodecs`")
	}
}

// TestRANS4x16_X32RoundTrip exercises encode→decode as the identity for
// the 32-way coder across a spread of input shapes and both orders — no
// external fixtures needed.
func TestRANS4x16_X32RoundTrip(t *testing.T) {
	inputs := map[string][]byte{
		"empty":         {},
		"single":        {'A'},
		"tiny":          []byte("ACGT"),
		"thirtyone":     bytes.Repeat([]byte("ACGTACGT_"), 4)[:31],
		"thirtytwo":     bytes.Repeat([]byte("ACGTACGT_"), 4)[:32],
		"thirtythree":   bytes.Repeat([]byte("ACGTACGT_"), 4)[:33],
		"low-alpha":     []byte("ACGTACGTACGTNNNNACGTACGTNNNN"),
		"runs":          bytes.Repeat([]byte{'A'}, 4000),
		"two-symbol":    twoSymbol(8000),
		"adjacent-syms": adjacentSymbols(8000),
		"ascii-text":    []byte(repeat("the quick brown fox jumps over the lazy dog. ", 400)),
		"full-alpha":    fullAlphabet(60000),
		"random-large":  randomBytes(t, 120000),
		// The final byte 'Z' appears nowhere else, exercising the
		// order-1 path where T[final] == 0.
		"last-byte-unique": append(bytes.Repeat([]byte("ABC"), 200), 'Z'),
	}
	for name, in := range inputs {
		for _, order := range []int{0, 1} {
			t.Run(name+".o"+itoa(order), func(t *testing.T) {
				comp, err := RANS4x16Encode(in, x4x16X32|order)
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

// TestRANS4x16_X32TransformRoundTrip exercises the 32-way coder combined
// with the PACK/RLE/STRIPE transforms. The combination on-wire forms
// have no q-file compliance vectors in the vendored corpus, so a
// round-trip is the available check.
func TestRANS4x16_X32TransformRoundTrip(t *testing.T) {
	inputs := map[string][]byte{
		"low-alpha":    []byte("ACGTACGTACGTNNNNACGTACGTNNNN"),
		"runs":         bytes.Repeat([]byte{'A'}, 4000),
		"two-symbol":   twoSymbol(8000),
		"ascii-text":   []byte(repeat("the quick brown fox jumps over the lazy dog. ", 400)),
		"full-alpha":   fullAlphabet(40000),
		"random-large": randomBytes(t, 90000),
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
					comp, err := RANS4x16Encode(in, x4x16X32|transform.bits|order)
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
}

// TestRANS4x16_TransformComplianceVectors decodes the htscodecs r4x16
// transform vectors (PACK/RLE/STRIPE, with and without order-1) and
// asserts byte-for-byte against the expected raw data. It complements
// TestRANS4x16_ComplianceVectors, which covers the plain order-0/1
// streams.
func TestRANS4x16_TransformComplianceVectors(t *testing.T) {
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		want, ok := loadCorpus(t, qfile)
		if !ok {
			continue
		}
		for _, order := range transformVectorSuffixes {
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
		t.Fatalf("htscodecs submodule not initialised — compliance vectors unavailable; run `git submodule update --init reference_code/htscodecs`")
	}
}

// TestRANS4x16_TransformEncodeMatchesHTScodecs checks our PACK/RLE/STRIPE
// encoder produces byte-identical output to the htscodecs reference
// vectors. The transform pipeline (bit-packing, run-length splitting,
// stripe transposition and the per-stripe method search) is fully
// deterministic, so a byte-exact match is the strongest possible test.
func TestRANS4x16_TransformEncodeMatchesHTScodecs(t *testing.T) {
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		raw, ok := loadCorpus(t, qfile)
		if !ok {
			continue
		}
		for _, order := range transformVectorSuffixes {
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
		t.Fatalf("htscodecs submodule not initialised; run `git submodule update --init reference_code/htscodecs`")
	}
}

// TestRANS4x16_TransformRoundTrip exercises encode→decode as the
// identity for every transform/order combination across a spread of
// input shapes — no external fixtures needed.
func TestRANS4x16_TransformRoundTrip(t *testing.T) {
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
		// The final byte 'Z' appears nowhere else, so it is a valid
		// symbol that is never a context — exercising the order-1
		// path where T[final] == 0 (see encodeFreq1RANS4x16).
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
					comp, err := RANS4x16Encode(in, transform.bits|order)
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
}

// TestRANS4x16_StripeStreamCount checks that an explicit stripe count N,
// passed in bits 8-15 of the order argument, round-trips for N from 1 to
// 8 — the STRIPE format stores N in the stream and the decoder honours
// it.
func TestRANS4x16_StripeStreamCount(t *testing.T) {
	in := fullAlphabet(40000)
	for n := 1; n <= 8; n++ {
		t.Run("N="+itoa(n), func(t *testing.T) {
			comp, err := RANS4x16Encode(in, x4x16Stripe|(n<<8))
			if err != nil {
				t.Fatalf("encode: %v", err)
			}
			got, err := RANS4x16Decode(comp)
			if err != nil {
				t.Fatalf("decode: %v", err)
			}
			if !bytes.Equal(got, in) {
				t.Fatalf("N=%d round-trip mismatch (first diff at %d)", n, firstDiff(got, in))
			}
		})
	}
}

// sixteenSymbol builds data over a 16-symbol alphabet, the largest
// alphabet hts_pack will still bit-pack (2 symbols per byte).
func sixteenSymbol(n int) []byte {
	out := make([]byte, n)
	for i := range out {
		out[i] = byte(i*7%16) + 32
	}
	return out
}

// TestRANS4x16_PackAlphabetWidths round-trips PACK across each
// symbols-per-byte width hts_pack selects: 8 (alphabet 2), 4 (3-4), 2
// (5-16) and 0 (the constant, alphabet-1 case). It also covers a >16
// alphabet, where hts_pack declines and the encoder clears X_PACK.
func TestRANS4x16_PackAlphabetWidths(t *testing.T) {
	mk := func(alpha, n int) []byte {
		out := make([]byte, n)
		for i := range out {
			out[i] = byte(40 + i%alpha)
		}
		return out
	}
	cases := map[string][]byte{
		"constant-1":    bytes.Repeat([]byte{'Q'}, 3000),
		"alpha-2":       mk(2, 5000),
		"alpha-4":       mk(4, 5000),
		"alpha-9":       mk(9, 5000),
		"alpha-16":      mk(16, 5000),
		"alpha-17":      mk(17, 5000), // too wide to pack
		"odd-tail":      mk(9, 5003),
		"single-symbol": {'Z'},
	}
	for name, in := range cases {
		for _, order := range []int{0, 1} {
			t.Run(name+".o"+itoa(order), func(t *testing.T) {
				comp, err := RANS4x16Encode(in, x4x16Pack|order)
				if err != nil {
					t.Fatalf("encode: %v", err)
				}
				got, err := RANS4x16Decode(comp)
				if err != nil {
					t.Fatalf("decode: %v", err)
				}
				if !bytes.Equal(got, in) {
					t.Fatalf("round-trip mismatch (first diff at %d)", firstDiff(got, in))
				}
			})
		}
	}
}

// TestRANS4x16_TransformXCATFallback checks that a transform requested
// on data rANS cannot shrink still round-trips: the encoder falls back
// to X_CAT for the (transformed) payload and the decoder reverses the
// transform regardless.
func TestRANS4x16_TransformXCATFallback(t *testing.T) {
	// Six random bytes: far too short for rANS to beat storing verbatim.
	in := []byte{0x11, 0x22, 0x33, 0x44, 0x55, 0x66}
	for _, bits := range []int{x4x16Pack, x4x16RLE, x4x16Pack | x4x16RLE} {
		comp, err := RANS4x16Encode(in, bits)
		if err != nil {
			t.Fatalf("encode bits 0x%02x: %v", bits, err)
		}
		got, err := RANS4x16Decode(comp)
		if err != nil {
			t.Fatalf("decode bits 0x%02x: %v", bits, err)
		}
		if !bytes.Equal(got, in) {
			t.Fatalf("bits 0x%02x round-trip mismatch", bits)
		}
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
		{"pack meta truncated", []byte{0x80, 0x00}},
		{"rle meta truncated", []byte{0x40, 0x00}},
		{"stripe truncated", []byte{0x08, 0x00}},
		// X_32 with a non-empty declared size but a payload too short
		// for 32 rANS states.
		{"x32 order-0 payload too short", []byte{0x04, 0x04, 1, 2, 3}},
		{"x32 order-1 payload too short", []byte{0x05, 0x04, 1, 2, 3}},
		{"nosz top-level malformed", []byte{0x10, 0x00}},
		{"truncated size varint", []byte{0x00, 0x80}},
		{"cat payload too short", []byte{0x20, 0x05, 'A', 'B'}},
		{"order-0 payload too short", []byte{0x00, 0x04, 1, 2, 3}},
		{"order-1 payload too short", []byte{0x01, 0x04, 1, 2, 3}},
		// payload byte 0xB0 selects table precision 11 (0xB0>>4),
		// which is neither 10 nor 12 and must be rejected.
		{"order-1 invalid shift 11", append([]byte{0x01, 0x08, 0xB0}, make([]byte, 20)...)},
		// Order-1 frequency-table bomb: a compressed-header order-1
		// stream (payload byte 0xA1) declaring a 2^30-byte uncompressed
		// table. The decoder must reject the oversized table size up
		// front rather than allocate ~1 GiB and run a billion-iteration
		// decode loop. Layout: format 0x01, rawSize=100, then a 16-byte
		// payload: 0xA1, uFreqSz varint 2^30, cFreqSz varint 1, padding.
		{"order-1 freq-table bomb", []byte{
			0x01, 0x64,
			0xA1, 0x84, 0x80, 0x80, 0x80, 0x00, 0x01,
			0, 0, 0, 0, 0, 0, 0, 0, 0,
		}},
		// RLE bomb: an X_RLE|X_CAT stream declaring a 10-byte output
		// whose single run-length varint is 0xFFFFFFFF. The decoder
		// must reject the run before expanding it, not hang. Layout:
		// format, osz=10, uMetaSize=15 (raw, odd), rleLen=1, then the
		// 7-byte raw meta [nsyms=1, sym 'A', run varint 0xFFFFFFFF],
		// then the single literal 'A'.
		{"rle run-length bomb", []byte{
			0x60, 0x0A, 0x0F, 0x01,
			0x01, 'A', 0x8F, 0xFF, 0xFF, 0xFF, 0x7F,
			'A',
		}},
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
