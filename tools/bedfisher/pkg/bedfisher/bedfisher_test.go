package bedfisher

import (
	"bytes"
	"math"
	"strings"
	"testing"
)

// runFisher feeds string inputs through Run and returns the resulting Result
// plus the rendered report.
func runFisher(t *testing.T, a, b, g string, opts Options) (*Result, string) {
	t.Helper()
	var out bytes.Buffer
	res, err := Run(strings.NewReader(a), strings.NewReader(b), strings.NewReader(g), &out, opts)
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	return res, out.String()
}

// TestHandComputed_Symmetric2x2:
//
//	Genome: chr1 length 100.
//	A: 0-10, 20-30, 40-50  -> qCount=3, qUnion=30, qMean=11.
//	B: 5-15, 45-55         -> dCount=2, dUnion=20, dMean=11.
//	Overlaps: (0-10,5-15), (40-50,45-55) -> n11=2.
//	n12=3-2=1, n21=2-2=0.
//	bMean = 11 + 11 = 22. genome/bMean = 100/22 = 4 (truncated).
//	n22_full = max(0+1+2, 4) = 4. n22 = 4-1-0-2 = 1.
func TestHandComputed_Symmetric2x2(t *testing.T) {
	a := "chr1\t0\t10\nchr1\t20\t30\nchr1\t40\t50\n"
	b := "chr1\t5\t15\nchr1\t45\t55\n"
	g := "chr1\t100\n"
	res, out := runFisher(t, a, b, g, Options{})
	if res.OverlapPair != 2 || res.N22Full != 4 || res.N22 != 1 || res.N12 != 1 || res.N21 != 0 {
		t.Errorf("contingency mismatch: %+v", res)
	}
	if !strings.Contains(out, "# Number of overlaps: 2") {
		t.Errorf("missing overlap line: %s", out)
	}
}

// TestHandComputed_NoOverlap: disjoint A and B over a tiny genome.
func TestHandComputed_NoOverlap(t *testing.T) {
	a := "chr1\t0\t10\n"
	b := "chr1\t50\t60\n"
	g := "chr1\t100\n"
	res, _ := runFisher(t, a, b, g, Options{})
	if res.OverlapPair != 0 {
		t.Errorf("expected 0 overlaps, got %d", res.OverlapPair)
	}
	// n11=0, n12=1, n21=1. qMean=11, dMean=11, bMean=22. 100/22=4.
	// n22_full = max(2,4) = 4. n22=4-1-1-0=2.
	if res.N22Full != 4 || res.N22 != 2 {
		t.Errorf("expected n22_full=4 n22=2, got %+v", res)
	}
}

// TestHandComputed_OneSideEmpty: empty A is degenerate; we still emit a report
// without panicking. Output has 0 query intervals; ratio is -nan.
func TestHandComputed_OneSideEmpty(t *testing.T) {
	g := "chr1\t100\n"
	res, _ := runFisher(t, "", "chr1\t0\t10\n", g, Options{})
	if res.QueryCount != 0 || res.DBCount != 1 || res.OverlapPair != 0 {
		t.Errorf("expected qc=0 dc=1 n11=0, got %+v", res)
	}
	if !math.IsNaN(res.Ratio) {
		t.Errorf("expected NaN ratio for empty A, got %v", res.Ratio)
	}
}

// TestStrandSame: -s requires same-strand overlap to count.
func TestStrandSame(t *testing.T) {
	a := "chr1\t0\t10\tx\t0\t+\n"
	b := "chr1\t5\t15\ty\t0\t-\n"
	g := "chr1\t100\n"
	res, _ := runFisher(t, a, b, g, Options{SameStrand: true})
	if res.OverlapPair != 0 {
		t.Errorf("expected 0 overlaps with opposite strands under -s, got %d", res.OverlapPair)
	}
	res, _ = runFisher(t, a, "chr1\t5\t15\ty\t0\t+\n", g, Options{SameStrand: true})
	if res.OverlapPair != 1 {
		t.Errorf("expected 1 overlap with same strands under -s, got %d", res.OverlapPair)
	}
}

// TestFractionFilter exercises -f: overlap too short fails the threshold.
func TestFractionFilter(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr1\t99\t101\n"
	g := "chr1\t500\n"
	// 1-base overlap / 100-base A = 0.01; -f 0.5 should drop it.
	res, _ := runFisher(t, a, b, g, Options{FractionA: 0.5})
	if res.OverlapPair != 0 {
		t.Errorf("expected 0 overlaps under -f 0.5, got %d", res.OverlapPair)
	}
}

// TestMergeA: a_merge has overlapping (10-20, 12-19); without -m we count 4
// query intervals, with -m we count 3.
func TestMergeA(t *testing.T) {
	a := "chr1\t10\t20\nchr1\t12\t19\nchr1\t30\t40\nchr1\t51\t52\n"
	b := "chr1\t15\t25\nchr1\t51\t52\n"
	g := "chr1\t60\n"
	resNoM, _ := runFisher(t, a, b, g, Options{})
	if resNoM.QueryCount != 4 {
		t.Errorf("without -m expected qc=4, got %d", resNoM.QueryCount)
	}
	resM, _ := runFisher(t, a, b, g, Options{MergeInputs: true})
	if resM.QueryCount != 3 {
		t.Errorf("with -m expected qc=3, got %d", resM.QueryCount)
	}
}

// TestUnit_OverlapCount_NonMonotonicEnds is a binary-free regression for the
// overlap-counting bug the parity pipeline found. B is sorted by start, but a
// long early-starting B (5-100) extends past the start of A (50-60). The old
// code binary-searched on ChromEnd over the start-sorted slice — invalid,
// because ChromEnd is not monotonic — and skipped that pair, under-counting.
// The correct count here is 2: A(50,60) overlaps both B(5,100) and B(55,65).
func TestUnit_OverlapCount_NonMonotonicEnds(t *testing.T) {
	a := "chr1\t50\t60\n"
	b := "chr1\t5\t100\nchr1\t10\t12\nchr1\t55\t65\n"
	g := "chr1\t1000\n"
	res, _ := runFisher(t, a, b, g, Options{})
	if res.OverlapPair != 2 {
		t.Fatalf("expected 2 overlaps (long early B must still be counted), got %d", res.OverlapPair)
	}
}

// TestUnit_OverlapCount_DuplicateAcrossA verifies a single B that overlaps
// several (self-overlapping) A records is counted once per A — the same
// intersection-pair accounting upstream's chromsweep uses.
func TestUnit_OverlapCount_DuplicateAcrossA(t *testing.T) {
	a := "chr1\t10\t30\nchr1\t20\t40\nchr1\t25\t50\n"
	b := "chr1\t22\t28\n" // overlaps all three A records.
	g := "chr1\t1000\n"
	res, _ := runFisher(t, a, b, g, Options{})
	if res.OverlapPair != 3 {
		t.Fatalf("expected 3 overlap pairs (one per overlapping A), got %d", res.OverlapPair)
	}
}

// TestKtFisherExact_Known2x2 verifies our Fisher's exact against a textbook
// case: table {1,9,11,3}. R's fisher.test gives two-tailed p = 0.00277 (~).
// Reference: also matches upstream htslib.
func TestKtFisherExact_Known2x2(t *testing.T) {
	left, right, two := ktFisherExact(1, 9, 11, 3)
	if math.Abs(two-0.00277) > 1e-4 {
		t.Errorf("two-tail p = %v, expected ~0.00277", two)
	}
	if !(left > 0 && left < 1) || !(right > 0 && right < 1) {
		t.Errorf("left/right out of (0,1): left=%v right=%v", left, right)
	}
}

// TestKtFisherExact_Degenerate covers min==max (a row or column sums to 0).
func TestKtFisherExact_Degenerate(t *testing.T) {
	l, r, two := ktFisherExact(0, 0, 5, 5)
	if l != 1 || r != 1 || two != 1 {
		t.Errorf("expected all 1s for degenerate table, got %v %v %v", l, r, two)
	}
}

// TestBadGenome triggers the genome-parse error path.
func TestBadGenome(t *testing.T) {
	_, err := Run(strings.NewReader(""), strings.NewReader(""), strings.NewReader("badformat\n"), &bytes.Buffer{}, Options{})
	if err == nil {
		t.Fatal("expected error for genome with one column")
	}
}

// TestBadBED triggers the BED-parse error path.
func TestBadBED(t *testing.T) {
	_, err := Run(strings.NewReader("not_a_record\n"), strings.NewReader(""), strings.NewReader("chr1\t100\n"), &bytes.Buffer{}, Options{})
	if err == nil {
		t.Fatal("expected error for malformed A")
	}
}

// TestMutuallyExclusiveFlags ensures -s and -S can't both be set.
func TestMutuallyExclusiveFlags(t *testing.T) {
	_, err := Run(strings.NewReader(""), strings.NewReader(""), strings.NewReader("chr1\t100\n"), &bytes.Buffer{}, Options{SameStrand: true, OppositeStrand: true})
	if err == nil {
		t.Fatal("expected error for -s + -S")
	}
}

// TestFractionOutOfRange covers the -f validation.
func TestFractionOutOfRange(t *testing.T) {
	g := "chr1\t100\n"
	if _, err := Run(strings.NewReader(""), strings.NewReader(""), strings.NewReader(g), &bytes.Buffer{}, Options{FractionA: 1.5}); err == nil {
		t.Fatal("expected error for -f 1.5")
	}
	if _, err := Run(strings.NewReader(""), strings.NewReader(""), strings.NewReader(g), &bytes.Buffer{}, Options{FractionB: -0.1}); err == nil {
		t.Fatal("expected error for -F -0.1")
	}
}
