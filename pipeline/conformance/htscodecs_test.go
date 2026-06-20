package conformance

import (
	"bytes"
	"os"
	"path/filepath"
	"strconv"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram/codec"
)

// htscodecsCorpus resolves reference_code/htscodecs/tests or SKIPs.
func htscodecsCorpus(t testing.TB) string {
	t.Helper()
	dir, ok := upstream.HtscodecsTestDir()
	if !ok {
		t.Skipf("htscodecs corpus not initialised (%s missing); run:\n"+
			"  git submodule update --init reference_code/htscodecs\n"+
			"see docs/CONFORMANCE.md", dir)
	}
	return dir
}

// loadRawCorpus reads dat/<name>, stripping any tab-delimited trailing column
// per the htscodecs corpus format (raw bytes are the first column of each
// line, concatenated).
func loadRawCorpus(dir, name string) ([]byte, bool) {
	raw, err := os.ReadFile(filepath.Join(dir, "dat", name))
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

// TestHtscodecs_RANS4x16 round-trips the htscodecs rANS-4x16 compliance vectors
// through OUR in-tree CRAM codec: decode the reference vector and assert it
// equals the raw input, then re-encode the raw input and assert it is
// byte-identical to the reference vector. Byte-identity here means our codec is
// bit-compatible with the reference C implementation's on-wire format — a CRAM
// written by us is readable by upstream and vice versa.
func TestHtscodecs_RANS4x16(t *testing.T) {
	dir := htscodecsCorpus(t)
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		raw, ok := loadRawCorpus(dir, qfile)
		if !ok {
			continue
		}
		for _, order := range []int{0, 1} {
			vec, err := os.ReadFile(filepath.Join(dir, "dat", "r4x16", qfile+"."+strconv.Itoa(order)))
			if err != nil {
				continue
			}
			ran++
			t.Run(qfile+".o"+strconv.Itoa(order), func(t *testing.T) {
				dec, err := codec.RANS4x16Decode(vec)
				if err != nil {
					t.Fatalf("decode: %v", err)
				}
				if !bytes.Equal(dec, raw) {
					t.Errorf("decode != raw input (%d vs %d bytes)", len(dec), len(raw))
				}
				enc, err := codec.RANS4x16Encode(raw, order)
				if err != nil {
					t.Fatalf("encode: %v", err)
				}
				if !bytes.Equal(enc, vec) {
					t.Errorf("re-encode != reference vector (%d vs %d bytes)", len(enc), len(vec))
				}
			})
		}
	}
	if ran == 0 {
		t.Skip("no r4x16 vectors present in corpus")
	}
}

// TestHtscodecs_RANS4x8 does the same for the rANS-4x8 vectors (the older CRAM
// 3.0 codec).
func TestHtscodecs_RANS4x8(t *testing.T) {
	dir := htscodecsCorpus(t)
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		raw, ok := loadRawCorpus(dir, qfile)
		if !ok {
			continue
		}
		for _, order := range []int{0, 1} {
			vec, err := os.ReadFile(filepath.Join(dir, "dat", "r4x8", qfile+"."+strconv.Itoa(order)))
			if err != nil {
				continue
			}
			ran++
			t.Run(qfile+".o"+strconv.Itoa(order), func(t *testing.T) {
				dec, err := codec.RANS4x8Decode(vec)
				if err != nil {
					t.Fatalf("decode: %v", err)
				}
				if !bytes.Equal(dec, raw) {
					t.Errorf("decode != raw input (%d vs %d bytes)", len(dec), len(raw))
				}
				enc, err := codec.RANS4x8Encode(raw, order)
				if err != nil {
					t.Fatalf("encode: %v", err)
				}
				if !bytes.Equal(enc, vec) {
					t.Errorf("re-encode != reference vector (%d vs %d bytes)", len(enc), len(vec))
				}
			})
		}
	}
	if ran == 0 {
		t.Skip("no r4x8 vectors present in corpus")
	}
}

// TestHtscodecs_Arith round-trips the adaptive-arithmetic compliance vectors
// (decode-only byte-identity; the arith encoder is range-coder based and the
// reference vectors are the oracle for decode).
func TestHtscodecs_Arith(t *testing.T) {
	dir := htscodecsCorpus(t)
	ran := 0
	for _, qfile := range []string{"q4", "q8", "qvar", "q40+dir"} {
		raw, ok := loadRawCorpus(dir, qfile)
		if !ok {
			continue
		}
		// arith vectors are named dat/arith/<qfile>.<order>; not every order
		// exists for every input.
		matches, _ := filepath.Glob(filepath.Join(dir, "dat", "arith", qfile+".*"))
		for _, vecPath := range matches {
			vec, err := os.ReadFile(vecPath)
			if err != nil {
				continue
			}
			ran++
			t.Run(filepath.Base(vecPath), func(t *testing.T) {
				dec, err := codec.ArithDecode(vec)
				if err != nil {
					t.Fatalf("decode %s: %v", filepath.Base(vecPath), err)
				}
				if !bytes.Equal(dec, raw) {
					t.Errorf("arith decode != raw input (%d vs %d bytes)", len(dec), len(raw))
				}
			})
		}
	}
	if ran == 0 {
		t.Skip("no arith vectors present in corpus")
	}
}
