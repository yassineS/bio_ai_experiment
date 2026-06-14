package bedgenomecov

// Parity tests against the upstream bedtools genomecov test suite.
//
// Cases are mirrored from reference_code/bedtools/test/genomecov/test-genomecov.sh.
// Inputs and expected outputs live under tools/bedgenomecov/testdata/parity/.
// BAM/SAM input is supported (RunBAM / the `-ibam` flag), but the upstream
// BAM/CRAM fixtures (built from SAM by htsutil) aren't vendored, so the
// BAM-input cases stay skipped and BAM/SAM is covered instead by the
// TestRunBAM_* unit tests plus a live `-ibam` cross-check against
// bedtools v2.31.1. The BED-input tests (t11/t12/t13) and a BED12 `-split`
// case are exercised directly here.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

func readGenomecovParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runGenomecovParity(t *testing.T, bedFile, genomeFile string, opts Options) []byte {
	t.Helper()
	bed := readGenomecovParity(t, bedFile)
	gen := readGenomecovParity(t, genomeFile)
	g, err := ReadGenome(bytes.NewReader(gen))
	if err != nil {
		t.Fatalf("ReadGenome %s: %v", genomeFile, err)
	}
	var out bytes.Buffer
	if err := Run(bytes.NewReader(bed), g, &out, opts); err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	return out.Bytes()
}

// genomecov.t1..t10 — these use BAM/CRAM input. BAM/SAM input is now
// supported via RunBAM (the `-ibam` CLI flag); the genome is taken from the
// alignment header. These specific upstream cases need BAM fixtures built by
// htsutil, which we don't vendor — RunBAM is instead covered by the
// SAM-input unit tests (TestRunBAM_*) and a live `-ibam` cross-check against
// bedtools v2.31.1. CRAM input is not yet wired.
func TestParity_Genomecov_T1to10_BAMInputs(t *testing.T) {
	t.Skip("BAM/CRAM upstream fixtures not vendored; BAM/SAM input covered by TestRunBAM_* + live -ibam cross-check")
}

// genomecov.t11 — histogram (default) over y.bed, including chroms in the
// genome that have no coverage.
func TestParity_Genomecov_T11_Histogram(t *testing.T) {
	got := runGenomecovParity(t, "y.bed", "genome.txt", Options{Mode: ModeHistogram, Scale: 1.0})
	want := readGenomecovParity(t, "t11_hist.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t12 — `-bg`: non-zero runs of constant depth.
func TestParity_Genomecov_T12_BedGraph(t *testing.T) {
	got := runGenomecovParity(t, "y.bed", "genome.txt", Options{Mode: ModeBedGraph, Scale: 1.0})
	want := readGenomecovParity(t, "t12_bg.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t13 — `-bga`: include zero-depth runs too.
func TestParity_Genomecov_T13_BedGraphAll(t *testing.T) {
	got := runGenomecovParity(t, "y.bed", "genome.txt", Options{Mode: ModeBedGraphAll, Scale: 1.0})
	want := readGenomecovParity(t, "t13_bga.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov -split on BED12 input: a 3-block record contributes coverage
// over each block ([0,10)/[20,30)/[40,50)) rather than the whole [0,50)
// span. Upstream's -split tests use BAM, but -split is equally valid for
// BED12 input; the expected output is generated from bedtools v2.31.1.
func TestParity_Genomecov_SplitBED12(t *testing.T) {
	got := runGenomecovParity(t, "split_blocks.bed", "split.genome",
		Options{Mode: ModeBedGraphAll, Scale: 1.0, Split: true})
	want := readGenomecovParity(t, "split_bga.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t14 — `-pc` paired-end coverage and t15 `-fs` fragment size are
// BAM-only. Both are implemented (RunBAM via -pc / -fs) and covered by the
// SAM-input unit tests TestRunBAM_PairedCoverage / TestRunBAM_FragmentSize
// plus a live cross-check against bedtools v2.31.1; the upstream BAM fixtures
// (built by htsutil) are not vendored.
func TestParity_Genomecov_T14_PairedEnd(t *testing.T) {
	t.Skip("BAM fixture not vendored; -pc covered by TestRunBAM_PairedCoverage + live cross-check")
}
func TestParity_Genomecov_T15_FragmentSize(t *testing.T) {
	t.Skip("BAM fixture not vendored; -fs covered by TestRunBAM_FragmentSize + live cross-check")
}
func TestParity_Genomecov_T16_EmptyBAM(t *testing.T) {
	t.Skip("upstream BAM fixture not vendored; BAM input supported via -ibam (see TestRunBAM_*)")
}
func TestParity_Genomecov_T17_EmptyCRAM(t *testing.T) {
	t.Skip("CRAM input not yet wired into bedgenomecov")
}
func TestParity_Genomecov_T18_DeepSAM(t *testing.T) {
	t.Skip("upstream SAM fixture not vendored; SAM input supported via RunBAM (see TestRunBAM_SAMInput)")
}
