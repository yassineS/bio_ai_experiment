package bedcoverage

// Parity tests against the upstream `bedtools coverage` test suite.
//
// Cases are mirrored from reference_code/bedtools/test/coverage/test-coverage.sh.
// Inputs are vendored under tools/bedcoverage/testdata/parity/ (literal copies
// of the upstream fixtures). Expected outputs are the inline heredoc strings
// from the upstream script, copied verbatim.
//
// BAM/SAM input is supported (auto-detected) on both -a and -b; the BAM
// database (-b) cases t10..t13 run against a vendored copy of the upstream
// three_blocks_match.bam fixture. The library-level snapshot tests here cover
// the output shapes; the end-to-end CLI cases that exercise the flag surface
// (legacy -abam input, the -sorted no-op, the float32 covered-fraction column,
// and the exact mutually-exclusive-modes stderr text) live in
// upstream_parity_test.go, which diffs this port's binary against a freshly
// built upstream `bedtools` byte-for-byte.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

func readParityFixture(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runParity(t *testing.T, aFile, bFile string, opts Options) []byte {
	t.Helper()
	a := readParityFixture(t, aFile)
	b := readParityFixture(t, bFile)
	var buf bytes.Buffer
	if _, err := Coverage(bytes.NewReader(a), bytes.NewReader(b), &buf, opts); err != nil {
		t.Fatalf("Coverage: %v", err)
	}
	return buf.Bytes()
}

// coverage.t1 — BAM on -a (no -split), BED on -b. We DO support BAM -a: the
// alignment is read through alnbed as a BED12 record and coverage is computed
// against its whole span. The only divergence from upstream's exact bytes is
// cosmetic: upstream emits BED12 block lists with a trailing comma
// ("10,10,10,") while our recordColumns serialiser omits it ("10,10,10").
// The numbers (count/covered/len/fraction) are byte-identical; the semantic
// BAM -a path is exercised below in TestCoverage_BAMQuery and the SAM unit
// test, so this case stays skipped only for the trailing-comma formatting nit.
// coverage.t1 — BAM query (-a): each alignment is echoed as its full BED12
// (CIGAR blocks, trailing-comma block lists) plus the coverage columns,
// byte-for-byte against bedtools v2.31.1.
func TestParity_Coverage_T1_BAMInput(t *testing.T) {
	got := runParity(t, "cov_a.bam", "cov_b.bed", Options{Mode: ModeDefault})
	want := readParityFixture(t, "t1_bam_a.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t2 — defaults: A, B BED; per-A count + bp + len + frac.
func TestParity_Coverage_T2_Default(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{})
	want := readParityFixture(t, "t2_default.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t3 — -counts.
func TestParity_Coverage_T3_Counts(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{Mode: ModeCounts})
	want := readParityFixture(t, "t3_counts.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t4 — -hist (per-A histogram + "all" footer).
func TestParity_Coverage_T4_Hist(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{Mode: ModeHist})
	want := readParityFixture(t, "t4_hist.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t5 — -d (per-base depth).
func TestParity_Coverage_T5_Depth(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{Mode: ModeDepth})
	want := readParityFixture(t, "t5_d.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t6 — -mean. Upstream prints with float32 precision
// ("1.3200001", "5.5599999"); our port uses float64 so we'd emit "1.32" and
// "5.56". Documented divergence; we still cover the mean op via unit tests.
// coverage.t6 — `-mean` prints the per-A mean depth as a float32 with 7
// decimals (float32 rounding noise included), byte-for-byte against upstream.
func TestParity_Coverage_T6_Mean(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{Mode: ModeMean})
	want := readParityFixture(t, "t6_mean.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t7 — -s (same strand).
func TestParity_Coverage_T7_SameStrand(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{SameStrand: true})
	want := readParityFixture(t, "t7_s.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t8 — -S (opposite strand).
func TestParity_Coverage_T8_OppositeStrand(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{OppositeStrand: true})
	want := readParityFixture(t, "t8_S.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t10 — -split with a BAM database (-b three_blocks_match.bam). The
// BAM's spliced alignment (10M10N10M10N10M) becomes a BED12 record whose three
// blocks (30 bp) are expanded by -split, so c.bed [0,50) is covered by 3
// blocks / 30 bp / fraction 0.6. Byte-for-byte against bedtools v2.31.1.
func TestParity_Coverage_T10_Split(t *testing.T) {
	got := runParity(t, "c.bed", "three_blocks_match.bam", Options{Mode: ModeDefault, Split: true})
	want := readParityFixture(t, "t10_split.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t11 — BAM database (-b) WITHOUT -split: the whole 50 bp read span
// covers c.bed [0,50) entirely (count 1 / 50 bp / fraction 1.0).
func TestParity_Coverage_T11_NoSplit(t *testing.T) {
	got := runParity(t, "c.bed", "three_blocks_match.bam", Options{Mode: ModeDefault})
	want := readParityFixture(t, "t11_nosplit.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t12 — BAM database (-b) with -split -d (per-base depth).
func TestParity_Coverage_T12_SplitDepth(t *testing.T) {
	got := runParity(t, "c.bed", "three_blocks_match.bam", Options{Mode: ModeDepth, Split: true})
	want := readParityFixture(t, "t12_split_d.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t13 — BAM database (-b) with -split -hist.
func TestParity_Coverage_T13_SplitHist(t *testing.T) {
	got := runParity(t, "c.bed", "three_blocks_match.bam", Options{Mode: ModeHist, Split: true})
	want := readParityFixture(t, "t13_split_hist.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage -split on a BED12 database (-b): the 3 blocks (30 bp) of the B
// record cover query [0,60), giving count 3 / covered 30 / fraction 0.5.
// Upstream's own -split tests use BAM; this BED-input case is verified
// byte-for-byte against bedtools v2.31.1.
func TestParity_Coverage_SplitBED12Database(t *testing.T) {
	got := runParity(t, "split_a.bed", "split_b12.bed", Options{Mode: ModeDefault, Split: true})
	want := readParityFixture(t, "split_default.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage -split must reject a BED12 query (-a) record (unsupported) with a
// clear error rather than emit a wrong answer.
func TestCoverage_SplitBED12QueryRejected(t *testing.T) {
	a := "chr1\t0\t50\tq\t0\t+\t0\t0\t0\t2\t10,10,\t0,40,\n"
	b := "chr1\t0\t60\tb\n"
	var buf bytes.Buffer
	if _, err := Coverage(bytes.NewReader([]byte(a)), bytes.NewReader([]byte(b)), &buf, Options{Split: true}); err == nil {
		t.Fatal("expected error for BED12 query under -split")
	}
}
