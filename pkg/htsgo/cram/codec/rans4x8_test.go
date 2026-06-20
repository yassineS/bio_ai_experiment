package codec

import (
	"bytes"
	"math/rand"
	"os"
	"path/filepath"
	"strconv"
	"testing"
)

// htscodecsDir is the test corpus from the samtools/htscodecs submodule.
// The compliance tests treat the corpus as a hard requirement: when the
// submodule isn't initialised they t.Fatalf with an init hint (run `git
// submodule update --init reference_code/htscodecs`) rather than skipping,
// so a missing corpus can never hide a parity gap. The round-trip and
// property tests need no external fixtures and always run. See
// docs/CRAM_ROADMAP.md §3.
const htscodecsDir = "../../../../reference_code/htscodecs/tests"

// loadCorpus returns the expected decoded bytes for a q-file. The
// htscodecs rans4x8.test harness feeds the codec `cut -f1 | tr -d
// '\012'` — i.e. for each line, the bytes before the first tab,
// concatenated with no newlines. Most q-files have no tabs (so cut -f1
// is the whole line); q40+dir is tab-delimited and the column-1 cut
// matters.
func loadCorpus(t *testing.T, name string) ([]byte, bool) {
	t.Helper()
	raw, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", name))
	if err != nil {
		return nil, false
	}
	var out []byte
	for _, line := range bytes.Split(raw, []byte{'\n'}) {
		if tab := bytes.IndexByte(line, '\t'); tab >= 0 {
			line = line[:tab]
		}
		out = append(out, line...)
	}
	return out, true
}

// TestRANS4x8_ComplianceVectors decodes the pre-compressed htscodecs
// vectors and asserts byte-for-byte against the expected raw data.
// This is the compliance oracle — it proves our decoder matches the
// reference C implementation's on-wire format exactly.
func TestRANS4x8_ComplianceVectors(t *testing.T) {
	cases := []struct {
		qfile string
		order int
	}{
		{"q4", 0}, {"q4", 1},
		{"q8", 0}, {"q8", 1},
		{"qvar", 0}, {"qvar", 1},
		{"q40+dir", 0}, {"q40+dir", 1},
	}
	ran := 0
	for _, c := range cases {
		want, ok := loadCorpus(t, c.qfile)
		if !ok {
			continue
		}
		compName := filepath.Join(htscodecsDir, "dat", "r4x8",
			c.qfile+"."+itoa(c.order))
		comp, err := os.ReadFile(compName)
		if err != nil {
			continue
		}
		ran++
		t.Run(c.qfile+".o"+itoa(c.order), func(t *testing.T) {
			got, err := RANS4x8Decode(comp)
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
		t.Skipf("htscodecs submodule not initialised — compliance vectors unavailable; run `git submodule update --init reference_code/htscodecs`")
	}
}

// TestRANS4x8_EncodeMatchesHTScodecs checks that our encoder produces
// byte-identical output to the htscodecs reference vectors. rANS
// encoding has no freedom once the frequency-normalisation algorithm
// is fixed, so a byte-exact match is both achievable and the strongest
// possible encoder test.
func TestRANS4x8_EncodeMatchesHTScodecs(t *testing.T) {
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		raw, ok := loadCorpus(t, qfile)
		if !ok {
			continue
		}
		for _, order := range []int{0, 1} {
			comp, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "r4x8",
				qfile+"."+itoa(order)))
			if err != nil {
				continue
			}
			ran++
			t.Run(qfile+".o"+itoa(order), func(t *testing.T) {
				got, err := RANS4x8Encode(raw, order)
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
		t.Skipf("htscodecs submodule not initialised; run `git submodule update --init reference_code/htscodecs`")
	}
}

// TestRANS4x8_RoundTrip exercises encode→decode as the identity across
// a spread of input shapes — no external fixtures needed.
func TestRANS4x8_RoundTrip(t *testing.T) {
	inputs := map[string][]byte{
		"empty":        {},
		"single":       {'A'},
		"two":          {'A', 'C'},
		"three":        {'A', 'C', 'G'},
		"four":         {'A', 'C', 'G', 'T'},
		"uniform":      bytes.Repeat([]byte{'N'}, 5000),
		"two-symbol":   twoSymbol(4000),
		"ascii-text":   []byte("the quick brown fox jumps over the lazy dog, repeatedly. " + repeat("ACGTN", 600)),
		"full-alpha":   fullAlphabet(8192),
		"random-small": randomBytes(t, 37),
		"random-large": randomBytes(t, 200000),
	}
	for name, in := range inputs {
		for _, order := range []int{0, 1} {
			t.Run(name+".o"+itoa(order), func(t *testing.T) {
				comp, err := RANS4x8Encode(in, order)
				if err != nil {
					t.Fatalf("encode: %v", err)
				}
				got, err := RANS4x8Decode(comp)
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

// TestRANS4x8_DecodeErrors pins the error paths for malformed input.
func TestRANS4x8_DecodeErrors(t *testing.T) {
	cases := []struct {
		name string
		in   []byte
	}{
		{"too short", []byte{0, 1, 2}},
		{"bad order", []byte{9, 0, 0, 0, 0, 0, 0, 0, 0}},
		{"compressed-size mismatch", []byte{0, 99, 0, 0, 0, 0, 0, 0, 0}},
		{"truncated freq table", append([]byte{0, 1, 0, 0, 0, 4, 0, 0, 0}, 'A')},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if _, err := RANS4x8Decode(c.in); err == nil {
				t.Errorf("expected error for %s", c.name)
			}
		})
	}
}

// --- helpers -----------------------------------------------------------------

// itoa renders a vector order byte for filenames and subtest names. The
// transform suffixes (64, 128, 192, …) are multi-digit, so it must be a
// real base-10 conversion, not a single-digit shortcut.
func itoa(n int) string {
	return strconv.Itoa(n)
}

func firstDiff(a, b []byte) int {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	for i := 0; i < n; i++ {
		if a[i] != b[i] {
			return i
		}
	}
	if len(a) != len(b) {
		return n
	}
	return -1
}

func twoSymbol(n int) []byte {
	out := make([]byte, n)
	for i := range out {
		if i%7 == 0 {
			out[i] = 'X'
		} else {
			out[i] = 'Y'
		}
	}
	return out
}

func fullAlphabet(n int) []byte {
	out := make([]byte, n)
	for i := range out {
		out[i] = byte(i * 31 % 256)
	}
	return out
}

func repeat(s string, n int) string {
	var b bytes.Buffer
	for i := 0; i < n; i++ {
		b.WriteString(s)
	}
	return b.String()
}

func randomBytes(t *testing.T, n int) []byte {
	t.Helper()
	r := rand.New(rand.NewSource(int64(n) * 2654435761))
	out := make([]byte, n)
	// A skewed distribution (squared uniform) — closer to real
	// quality-score data than flat random, exercising the frequency
	// normaliser harder.
	for i := range out {
		x := r.Float64()
		out[i] = byte(x * x * 256)
	}
	return out
}
