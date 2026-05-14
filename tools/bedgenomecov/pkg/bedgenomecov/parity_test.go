package bedgenomecov

// Parity tests against the upstream bedtools genomecov test suite.
//
// Cases are mirrored from reference_code/bedtools/test/genomecov/test-genomecov.sh.
// Inputs and expected outputs live under tools/bedgenomecov/testdata/parity/.
// bedgenomecov only consumes BED input today (no BAM/CRAM/SAM parser), so
// most upstream tests — which use `bedtools genomecov -ibam` on a BAM built
// from a SAM fixture by htsutil — are skipped. The three BED-input tests
// (t11/t12/t13) are exercised here.

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

// genomecov.t1..t10 — all use BAM/CRAM input. bedgenomecov is BED-only.
func TestParity_Genomecov_T1to10_BAMInputs(t *testing.T) {
	t.Skip("unimplemented: BAM/SAM/CRAM input. bedgenomecov consumes BED only.")
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

// genomecov.t14..t18 — `-pc` paired-end coverage, `-fs` fragment size, BAM
// empty fixtures, deep SAM. All BAM-only.
func TestParity_Genomecov_T14_PairedEnd(t *testing.T) {
	t.Skip("unimplemented: -pc paired-end coverage (BAM-only feature)")
}
func TestParity_Genomecov_T15_FragmentSize(t *testing.T) {
	t.Skip("unimplemented: -fs fragment size (BAM-only feature)")
}
func TestParity_Genomecov_T16_EmptyBAM(t *testing.T) {
	t.Skip("unimplemented: BAM input")
}
func TestParity_Genomecov_T17_EmptyCRAM(t *testing.T) {
	t.Skip("unimplemented: CRAM input")
}
func TestParity_Genomecov_T18_DeepSAM(t *testing.T) {
	t.Skip("unimplemented: SAM input")
}
