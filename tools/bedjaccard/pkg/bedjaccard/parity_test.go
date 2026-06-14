package bedjaccard

// Parity tests against the upstream bedtools jaccard test suite.
//
// Cases are mirrored from reference_code/bedtools/test/jaccard/test-jaccard.sh.
// Inputs and expected outputs live under tools/bedjaccard/testdata/parity/.
//
// As of the column-ops + discrepancies wave, bedjaccard now pre-merges
// both A and B before computing intersection / union (mirroring upstream's
// `setUseMergedIntervals(true)` in ContextJaccard.cpp), so cases against
// b.bed / c.bed / mixedStrands.bed are byte-for-byte parity with upstream.
// The `-S` single-strand filter (t12/t13), `-split` BED12 block-aware overlap
// (t08), and BAM input (t09, auto-detected) are implemented and asserted.
// Cases still wrapped in t.Skip are unrelated (VCF/GFF input, large fixtures).

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

func readJaccardParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runJaccardParity(t *testing.T, aFile, bFile string, opts Options) []byte {
	t.Helper()
	a := readJaccardParity(t, aFile)
	b := readJaccardParity(t, bFile)
	var out bytes.Buffer
	if _, err := Run(bytes.NewReader(a), bytes.NewReader(b), &out, opts); err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	return out.Bytes()
}

// jaccard.t01 — self-intersect on a.bed (no overlapping records).
func TestParity_Jaccard_T01_SelfIntersect(t *testing.T) {
	got := runJaccardParity(t, "a.bed", "a.bed", Options{})
	want := readJaccardParity(t, "t01_self.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// jaccard.t02 — a.bed vs b.bed. b.bed has overlapping records (b2/b3 at
// 90-101 and 100-110, both covered by a2 100-200); upstream pre-merges B
// before counting. With the merge fix in place this is byte-for-byte parity.
func TestParity_Jaccard_T02_AvsB(t *testing.T) {
	got := runJaccardParity(t, "a.bed", "b.bed", Options{})
	want := readJaccardParity(t, "t02_a_b.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// jaccard.t03 — a.bed vs c.bed. c.bed has a single record; the merge
// changes things only on the A side (after pre-merge A is two disjoint
// records, same as the input).
func TestParity_Jaccard_T03_AvsC(t *testing.T) {
	got := runJaccardParity(t, "a.bed", "c.bed", Options{})
	want := readJaccardParity(t, "t03_a_c.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// jaccard.t05 — same input as t02 but conceptually piped via stdin.
// bedjaccard takes both inputs as io.Reader so stdin vs. file is a no-op
// at the package boundary; reuse t02's fixtures.
func TestParity_Jaccard_T05_StdinA(t *testing.T) {
	got := runJaccardParity(t, "a.bed", "b.bed", Options{})
	want := readJaccardParity(t, "t02_a_b.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// jaccard.t06 — symmetry: jaccard(A, B) byte-equals jaccard(B, A). With
// pre-merge in place both orderings produce the same totals.
func TestParity_Jaccard_T06_Symmetry(t *testing.T) {
	fwd := runJaccardParity(t, "a.bed", "b.bed", Options{})
	rev := runJaccardParity(t, "b.bed", "a.bed", Options{})
	if !bytes.Equal(fwd, rev) {
		t.Fatalf("symmetry violated.\nfwd:\n%s\nrev:\n%s", fwd, rev)
	}
}

// jaccard.t07 — three_blocks_match.bed (BED12 single record) vs e.bed.
// Upstream WITHOUT -split treats the BED12 as a single 0..50 interval; that
// matches what our BED reader does, so this case passes.
func TestParity_Jaccard_T07_ThreeBlocksNoSplit(t *testing.T) {
	got := runJaccardParity(t, "three_blocks_match.bed", "e.bed", Options{})
	want := readJaccardParity(t, "t07_three_blocks.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// jaccard.t08 — `-split` interprets the BED12 blocks: three_blocks_match.bed
// becomes blocks [0,10)/[20,30)/[40,50) (total 30) and only the first block
// overlaps e.bed [5,15), giving intersection 5 / union 35.
func TestParity_Jaccard_T08_ThreeBlocksSplit(t *testing.T) {
	got := runJaccardParity(t, "three_blocks_match.bed", "e.bed", Options{Split: true})
	want := readJaccardParity(t, "t08_three_blocks_split.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// jaccard.t09 — BAM input on both -a and -b. bedjaccard auto-detects BAM and
// converts each alignment to a BED12 interval (its CIGAR blocks). Expected
// value is from the upstream binary on the vendored a.bam / three_blocks_match.bam.
func TestParity_Jaccard_T09_BAMInput(t *testing.T) {
	got := runJaccardParity(t, "a.bam", "three_blocks_match.bam", Options{})
	want := readJaccardParity(t, "t09_bam.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// jaccard.t10 — mixed-strand files, no `-s`. Both A and B are pre-merged
// across all strands before the sweep.
func TestParity_Jaccard_T10_MixedStrandsNoFlag(t *testing.T) {
	got := runJaccardParity(t, "aMixedStrands.bed", "bMixedStrands.bed", Options{})
	want := readJaccardParity(t, "t10_mixed_nostrand.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// jaccard.t11 — mixed-strand files with `-s` (same-strand). Pre-merge
// runs per-strand so cross-strand records don't collapse.
func TestParity_Jaccard_T11_MixedStrandsS(t *testing.T) {
	got := runJaccardParity(t, "aMixedStrands.bed", "bMixedStrands.bed",
		Options{SameStrand: true})
	want := readJaccardParity(t, "t11_mixed_s.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// jaccard.t12 — `-S +` restricts both inputs to the forward strand.
func TestParity_Jaccard_T12_MixedStrandsSPlus(t *testing.T) {
	got := runJaccardParity(t, "aMixedStrands.bed", "bMixedStrands.bed",
		Options{StrandFilter: "+"})
	want := readJaccardParity(t, "t12_mixed_s_plus.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// jaccard.t13 — `-S -` restricts both inputs to the reverse strand.
func TestParity_Jaccard_T13_MixedStrandsSMinus(t *testing.T) {
	got := runJaccardParity(t, "aMixedStrands.bed", "bMixedStrands.bed",
		Options{StrandFilter: "-"})
	want := readJaccardParity(t, "t13_mixed_s_minus.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// jaccard.t14 — a645.bed vs b645.bed: each side has disjoint records, no
// pre-merge gap.
func TestParity_Jaccard_T14_A645vsB645(t *testing.T) {
	got := runJaccardParity(t, "a645.bed", "b645.bed", Options{})
	want := readJaccardParity(t, "t14_a645_b645.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// jaccard.t15 — same as t14 but with -a and -b swapped: jaccard is symmetric.
func TestParity_Jaccard_T15_A645vsB645Symmetry(t *testing.T) {
	// The shell test only checks that the two output files are identical
	// (jaccard is symmetric). The expected value is the same as t14.
	got := runJaccardParity(t, "a645.bed", "b645.bed", Options{})
	rev := runJaccardParity(t, "b645.bed", "a645.bed", Options{})
	if !bytes.Equal(got, rev) {
		t.Fatalf("symmetry violated.\nfwd:\n%s\nrev:\n%s", got, rev)
	}
}

// jaccard.t16 — long.bed vs short.bed (the giant CSHL test fixture); too
// large to vendor under testdata/parity for this round.
func TestParity_Jaccard_T16_LongShort(t *testing.T) {
	t.Skip("not vendored: long.bed/short.bed are large fixtures we don't ship in the parity testdata")
}
