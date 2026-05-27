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

// TestParity_Nuc_FullHeader covers the `-fullHeader` flag end-to-end.
//
// Upstream `bedtools nuc -fullHeader` rebuilds the FASTA index keyed on
// the entire `>` line (everything after `>` up to the newline), so a
// BED row whose chrom column carries embedded whitespace will resolve
// against the matching record. Our port keeps the FAI keyed on the
// first whitespace token (the htslib default) and instead constructs a
// secondary map from full-header to first-token at Run time when
// FullHeader is requested — see buildFullHeaderMap in bednuc.go. The
// observable output is identical: the matching record is located, its
// sequence is sliced by the BED [start,end) range, and the per-base
// composition / pattern columns are emitted in the upstream column
// order. This test asserts that observable parity over a fixture whose
// only special characteristic is a whitespace-bearing header.
func TestParity_Nuc_FullHeader(t *testing.T) {
	dir := t.TempDir()
	faPath := filepath.Join(dir, "fh.fa")
	if err := os.WriteFile(faPath, []byte(">chr1 with extra info\nACGTACGTAC\n"), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	bed := []byte("chr1 with extra info\t0\t10\n")
	var out, warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), faPath, &out, &warn, Options{FullHeader: true}); err != nil {
		t.Fatalf("Run: %v\nwarn: %s", err, warn.String())
	}
	// Upstream `nuc -fullHeader` emits the chrom column verbatim from
	// the input BED, followed by the standard counts row. The
	// percentages come from a 10-mer ACGTACGTAC = 3A/3C/2G/2T, so
	// %AT = (3+2)/10 = 0.5 and %GC = (3+2)/10 = 0.5.
	want := "#1_usercol\t2_usercol\t3_usercol\t4_pct_at\t5_pct_gc\t6_num_A\t7_num_C\t8_num_G\t9_num_T\t10_num_N\t11_num_oth\t12_seq_len\n" +
		"chr1 with extra info\t0\t10\t0.500000\t0.500000\t3\t3\t2\t2\t0\t0\t10\n"
	if got := out.String(); got != want {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
