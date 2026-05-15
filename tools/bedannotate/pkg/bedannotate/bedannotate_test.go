package bedannotate

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// run is the standard harness: feed A as a string, N B files as strings,
// and return Run's stdout.
func run(t *testing.T, a string, bs []string, opts Options) string {
	t.Helper()
	bRs := make([]io.Reader, len(bs))
	for i, s := range bs {
		bRs[i] = strings.NewReader(s)
	}
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), bRs, &out, opts); err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	return out.String()
}

// TestRun_DefaultFraction is a hand-computed end-to-end check.
//
// A = chr1 0 10. B1 = {[2,5),[8,12)} so covers 3+2 = 5 bases of A (frac
// 0.5). B2 = {[0,3)} so covers 3 bases (frac 0.3).
func TestRun_DefaultFraction(t *testing.T) {
	a := "chr1\t0\t10\n"
	b1 := "chr1\t2\t5\nchr1\t8\t12\n"
	b2 := "chr1\t0\t3\n"
	got := run(t, a, []string{b1, b2}, Options{})
	want := "chr1\t0\t10\t0.500000\t0.300000\n"
	if got != want {
		t.Errorf("default mode mismatch.\nwant: %q\n got: %q", want, got)
	}
}

// TestRun_CountsMode emits per-B integer counts in place of fractions.
// Same fixtures as above: B1 has 2 overlapping records, B2 has 1.
func TestRun_CountsMode(t *testing.T) {
	a := "chr1\t0\t10\n"
	b1 := "chr1\t2\t5\nchr1\t8\t12\n"
	b2 := "chr1\t0\t3\n"
	got := run(t, a, []string{b1, b2}, Options{Mode: ModeCounts})
	want := "chr1\t0\t10\t2\t1\n"
	if got != want {
		t.Errorf("counts mode mismatch.\nwant: %q\n got: %q", want, got)
	}
}

// TestRun_BothMode interleaves count then fraction per B.
func TestRun_BothMode(t *testing.T) {
	a := "chr1\t0\t10\n"
	b1 := "chr1\t2\t5\nchr1\t8\t12\n"
	b2 := "chr1\t0\t3\n"
	got := run(t, a, []string{b1, b2}, Options{Mode: ModeBoth})
	want := "chr1\t0\t10\t2\t0.500000\t1\t0.300000\n"
	if got != want {
		t.Errorf("both mode mismatch.\nwant: %q\n got: %q", want, got)
	}
}

// TestRun_HeaderEmittedWhenNamesSet checks the leading '#' line.
func TestRun_HeaderEmittedWhenNamesSet(t *testing.T) {
	a := "chr1\t0\t10\n"
	b1 := "chr1\t2\t5\n"
	got := run(t, a, []string{b1}, Options{Names: []string{"exons"}})
	lines := strings.SplitN(got, "\n", 2)
	if lines[0] != "#\texons" {
		t.Errorf("header line mismatch: %q", lines[0])
	}
	// Both mode header format.
	got2 := run(t, a, []string{b1}, Options{Names: []string{"exons"}, Mode: ModeBoth})
	if !strings.HasPrefix(got2, "#\texons_cnt\texons_pct\n") {
		t.Errorf("both header mismatch: %q", got2)
	}
}

// TestRun_StrandFilters: B records must match A on +/− according to -s
// and -S. A is '+'; B record on '−' contributes 0 with -s and 100% with
// -S.
func TestRun_StrandFilters(t *testing.T) {
	a := "chr1\t0\t10\tname\t0\t+\n"
	b := "chr1\t0\t10\tx\t0\t-\n"
	// -s: same-strand only → no match.
	got := run(t, a, []string{b}, Options{SameStrand: true, Mode: ModeBoth})
	if !strings.HasSuffix(got, "\t0\t0.000000\n") {
		t.Errorf("expected no overlap with -s, got %q", got)
	}
	// -S: opposite-strand → full overlap.
	got2 := run(t, a, []string{b}, Options{OppositeStrand: true, Mode: ModeBoth})
	if !strings.HasSuffix(got2, "\t1\t1.000000\n") {
		t.Errorf("expected full overlap with -S, got %q", got2)
	}
	// Empty strand on A → both -s and -S exclude.
	a2 := "chr1\t0\t10\n"
	got3 := run(t, a2, []string{b}, Options{SameStrand: true})
	if !strings.HasSuffix(got3, "\t0.000000\n") {
		t.Errorf("expected empty strand to be excluded under -s, got %q", got3)
	}
}

// TestRun_NoMatchAndEmptyB exercises zero-overlap edges.
func TestRun_NoMatchAndEmptyB(t *testing.T) {
	a := "chr1\t0\t10\nchr2\t0\t5\n"
	b := "chrX\t0\t100\n" // no chr1/chr2 records.
	got := run(t, a, []string{b}, Options{Mode: ModeBoth})
	wantLines := []string{
		"chr1\t0\t10\t0\t0.000000",
		"chr2\t0\t5\t0\t0.000000",
	}
	for _, line := range wantLines {
		if !strings.Contains(got, line) {
			t.Errorf("expected %q in:\n%s", line, got)
		}
	}
	// Also exercise an empty B file (no records at all).
	got2 := run(t, a, []string{""}, Options{})
	if !strings.Contains(got2, "chr1\t0\t10\t0.000000") {
		t.Errorf("empty B should yield 0.000000, got:\n%s", got2)
	}
}

// TestDefaultNames extracts basenames from paths with various separators.
func TestDefaultNames(t *testing.T) {
	got := DefaultNames([]string{"b1.bed", "dir/b2.bed", "a/b/c/long.bed"})
	want := []string{"b1.bed", "b2.bed", "long.bed"}
	for i, w := range want {
		if got[i] != w {
			t.Errorf("DefaultNames[%d]=%q, want %q", i, got[i], w)
		}
	}
}

// TestRun_RejectsConflictingStrandFlags ensures -s/-S can't be combined.
func TestRun_RejectsConflictingStrandFlags(t *testing.T) {
	_, err := Run(strings.NewReader(""), nil, &bytes.Buffer{},
		Options{SameStrand: true, OppositeStrand: true})
	if err == nil {
		t.Fatalf("expected error for -s + -S, got nil")
	}
}
