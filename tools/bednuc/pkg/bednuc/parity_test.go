package bednuc

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// Parity tests for bednuc.
//
// The upstream `bedtools` test suite at reference_code/bedtools/test/ does
// not ship a `nuc/` subdirectory (it's one of a handful of subcommands
// whose tests live only in the README), so the fixtures here were
// hand-computed against the upstream output format and the upstream
// counting rules in `src/utils/sequenceUtilities/sequenceUtils.cpp`:
//
//   * A/a, C/c, G/c, (T/t/U/u), N/n are counted as A/C/G/T/N; everything
//     else is "other".
//   * %AT = (A+T)/seq_len, %GC = (C+G)/seq_len, printed with `%f` (6
//     decimals).
//   * countPattern walks every position with overlapping matches; -C
//     forces case-insensitive matching.
//   * `-s` reverse-complements `-`-strand intervals before counting.
//
// Each fixture stores its own .fa + .bed + expected output so that we can
// diff byte-for-byte without re-deriving the expected text in the test
// helper. When upstream later ships a nuc/ test directory, these cases
// can be replaced by `t.Run` over `read(*.golden)` pairs.

func readFile(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runParity(t *testing.T, faName, bedName string, opts Options) []byte {
	t.Helper()
	// We need a real on-disk FASTA file so the index can be built. The
	// parity fixtures live under testdata/parity/.
	faPath := filepath.Join("..", "..", "testdata", "parity", faName)
	bed := readFile(t, bedName)
	var out bytes.Buffer
	var warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), faPath, &out, &warn, opts); err != nil {
		t.Fatalf("Run failed: %v\nwarn: %s", err, warn.String())
	}
	return out.Bytes()
}

// case1: three BED3 intervals across a periodic ACGT FASTA. Defaults only.
func TestParity_Nuc_Case1_Default(t *testing.T) {
	got := runParity(t, "case1.fa", "case1.bed", Options{})
	want := readFile(t, "case1.expected.txt")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// case2: BED6 plus and minus strands with -s + -seq. Plus strand should
// keep ACGTA; minus strand should RC to TACGT.
func TestParity_Nuc_Case2_StrandSeq(t *testing.T) {
	got := runParity(t, "case1.fa", "case2.bed", Options{ForceStrand: true, PrintSeq: true})
	want := readFile(t, "case2.expected.txt")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// case3: mixed-case FASTA with an N, exercised with -pattern (case
// sensitive — upstream default).
func TestParity_Nuc_Case3_PatternCase(t *testing.T) {
	got := runParity(t, "case3.fa", "case3.bed", Options{Pattern: "GG", HasPattern: true})
	want := readFile(t, "case3.expected.txt")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// Skipped: -fullHeader — our index always keys on the first whitespace
// token, and a best-effort fallback is built in Run. Exercising it would
// require a FASTA with whitespace in the header AND a BED keyed by that
// full header; upstream's `-fullHeader` semantics flip whether the index
// is built on the full header at indexing time, which we don't expose.
// Documented option-tail gap in PARITY_VALIDATION.md.
func TestParity_Nuc_FullHeader(t *testing.T) {
	t.Skip("-fullHeader: index-time semantics differ; tracked in PARITY_VALIDATION.md")
}
