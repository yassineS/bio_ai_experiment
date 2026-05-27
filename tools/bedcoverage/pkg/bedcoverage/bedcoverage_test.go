package bedcoverage

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// runCoverage is a tiny helper that runs the package on a pair of
// in-memory inputs and returns the produced output.
func runCoverage(t *testing.T, a, b string, opts Options) string {
	t.Helper()
	var buf bytes.Buffer
	if _, err := Coverage(strings.NewReader(a), strings.NewReader(b), &buf, opts); err != nil {
		t.Fatalf("Coverage: %v", err)
	}
	return buf.String()
}

func TestCoverage_DefaultBasic(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr1\t10\t50\nchr1\t60\t90\n"
	got := runCoverage(t, a, b, Options{})
	// 2 features overlap, 40+30=70 bp covered, len 100, fraction 0.7
	want := "chr1\t0\t100\t2\t70\t100\t0.7000000\n"
	if got != want {
		t.Errorf("default basic:\nwant: %q\ngot:  %q", want, got)
	}
}

func TestCoverage_DefaultNoOverlap(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr2\t10\t50\n"
	got := runCoverage(t, a, b, Options{})
	want := "chr1\t0\t100\t0\t0\t100\t0.0000000\n"
	if got != want {
		t.Errorf("no overlap: want %q, got %q", want, got)
	}
}

func TestCoverage_FullColumnsPreserved(t *testing.T) {
	// 6-field BED record. Output should preserve all 6 columns then append.
	a := "chr1\t20\t70\t6\t25\t+\n"
	b := "chr1\t30\t40\tb\t1\t-\n"
	got := runCoverage(t, a, b, Options{})
	want := "chr1\t20\t70\t6\t25\t+\t1\t10\t50\t0.2000000\n"
	if got != want {
		t.Errorf("6-col preserve: want %q, got %q", want, got)
	}
}

func TestCoverage_CountsMode(t *testing.T) {
	a := "chr1\t20\t70\t6\t25\t+\n"
	b := "chr1\t30\t40\nchr1\t50\t60\nchr1\t100\t110\n"
	got := runCoverage(t, a, b, Options{Mode: ModeCounts})
	want := "chr1\t20\t70\t6\t25\t+\t2\n"
	if got != want {
		t.Errorf("counts: want %q, got %q", want, got)
	}
}

func TestCoverage_DepthMode(t *testing.T) {
	a := "chr1\t0\t5\n"
	b := "chr1\t1\t3\nchr1\t2\t4\n"
	got := runCoverage(t, a, b, Options{Mode: ModeDepth})
	// 5 bases: pos 1..5 ⇒ depths 0,1,2,1,0
	want := "chr1\t0\t5\t1\t0\n" +
		"chr1\t0\t5\t2\t1\n" +
		"chr1\t0\t5\t3\t2\n" +
		"chr1\t0\t5\t4\t1\n" +
		"chr1\t0\t5\t5\t0\n"
	if got != want {
		t.Errorf("depth:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestCoverage_HistMode(t *testing.T) {
	a := "chr1\t0\t5\n"
	b := "chr1\t1\t3\nchr1\t2\t4\n"
	got := runCoverage(t, a, b, Options{Mode: ModeHist})
	// per-A: depth 0 ⇒ 2bp, depth 1 ⇒ 2bp, depth 2 ⇒ 1bp.
	// all: same numbers (one A only).
	want := "chr1\t0\t5\t0\t2\t5\t0.4000000\n" +
		"chr1\t0\t5\t1\t2\t5\t0.4000000\n" +
		"chr1\t0\t5\t2\t1\t5\t0.2000000\n" +
		"all\t0\t2\t5\t0.4000000\n" +
		"all\t1\t2\t5\t0.4000000\n" +
		"all\t2\t1\t5\t0.2000000\n"
	if got != want {
		t.Errorf("hist:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestCoverage_MeanMedian(t *testing.T) {
	a := "chr1\t0\t5\n"
	b := "chr1\t1\t3\nchr1\t2\t4\n"
	// depths: 0,1,2,1,0 ⇒ mean 0.8, median 1, min 0, max 2, sum 4
	cases := []struct {
		mode Mode
		want string
	}{
		{ModeMean, "chr1\t0\t5\t0.8\n"},
		{ModeMedian, "chr1\t0\t5\t1\n"},
		{ModeMin, "chr1\t0\t5\t0\n"},
		{ModeMax, "chr1\t0\t5\t2\n"},
		{ModeSum, "chr1\t0\t5\t4\n"},
	}
	for _, c := range cases {
		got := runCoverage(t, a, b, Options{Mode: c.mode})
		if got != c.want {
			t.Errorf("mode %v: want %q, got %q", c.mode, c.want, got)
		}
	}
}

func TestCoverage_StrandFilters(t *testing.T) {
	a := "chr1\t0\t100\tx\t0\t+\n"
	b := "chr1\t10\t50\tb1\t0\t+\nchr1\t60\t90\tb2\t0\t-\n"

	// Same strand: only b1.
	got := runCoverage(t, a, b, Options{SameStrand: true})
	want := "chr1\t0\t100\tx\t0\t+\t1\t40\t100\t0.4000000\n"
	if got != want {
		t.Errorf("same strand: want %q, got %q", want, got)
	}

	// Opposite strand: only b2.
	got = runCoverage(t, a, b, Options{OppositeStrand: true})
	want = "chr1\t0\t100\tx\t0\t+\t1\t30\t100\t0.3000000\n"
	if got != want {
		t.Errorf("opposite strand: want %q, got %q", want, got)
	}
}

func TestCoverage_FractionFilters(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr1\t90\t110\n" // 10bp overlap → fracA=0.1, fracB=0.5
	got := runCoverage(t, a, b, Options{FractionA: 0.5})
	want := "chr1\t0\t100\t0\t0\t100\t0.0000000\n"
	if got != want {
		t.Errorf("fracA reject: want %q, got %q", want, got)
	}

	got = runCoverage(t, a, b, Options{FractionA: 0.05})
	want = "chr1\t0\t100\t1\t10\t100\t0.1000000\n"
	if got != want {
		t.Errorf("fracA accept: want %q, got %q", want, got)
	}

	// Reciprocal: requires both.
	got = runCoverage(t, a, b, Options{FractionA: 0.5, FractionB: 0.5, Reciprocal: true})
	want = "chr1\t0\t100\t0\t0\t100\t0.0000000\n"
	if got != want {
		t.Errorf("reciprocal reject: want %q, got %q", want, got)
	}
}

func TestCoverage_EmptyB(t *testing.T) {
	got := runCoverage(t, "chr1\t0\t100\n", "", Options{})
	want := "chr1\t0\t100\t0\t0\t100\t0.0000000\n"
	if got != want {
		t.Errorf("empty B: want %q, got %q", want, got)
	}
}

func TestCoverage_EmptyA(t *testing.T) {
	got := runCoverage(t, "", "chr1\t0\t100\n", Options{})
	if got != "" {
		t.Errorf("empty A: want '', got %q", got)
	}
}

func TestCoverage_ErrorOnBadInput(t *testing.T) {
	_, err := Coverage(strings.NewReader("chr1\tbad\t10\n"), strings.NewReader(""), io.Discard, Options{})
	if err == nil {
		t.Errorf("expected parse error from bad A, got nil")
	}
	_, err = Coverage(strings.NewReader(""), strings.NewReader("chr1\tbad\t10\n"), io.Discard, Options{})
	if err == nil {
		t.Errorf("expected parse error from bad B, got nil")
	}
}

func TestCoverage_OverlappingB_DepthAndCoverage(t *testing.T) {
	// Two B records overlap each other; default mode should count both (count=2)
	// and bp_covered should be the union (not double-counted).
	a := "chr1\t0\t100\n"
	b := "chr1\t10\t60\nchr1\t40\t90\n"
	// Union of [10,60) ∪ [40,90) = [10,90) ⇒ 80 bp covered.
	got := runCoverage(t, a, b, Options{})
	want := "chr1\t0\t100\t2\t80\t100\t0.8000000\n"
	if got != want {
		t.Errorf("overlapping B: want %q, got %q", want, got)
	}
}

func TestCoverage_BClipping(t *testing.T) {
	// B record extends beyond A; depth should clip to A's bounds.
	a := "chr1\t10\t20\n"
	b := "chr1\t0\t100\n"
	got := runCoverage(t, a, b, Options{Mode: ModeMax})
	// All 10 bases of A are covered, max depth = 1.
	if got != "chr1\t10\t20\t1\n" {
		t.Errorf("clipping: got %q", got)
	}
}

func TestCoverage_ZeroLengthA(t *testing.T) {
	// Defensive: a zero-length A interval should not crash; covered=0, frac=0.
	a := "chr1\t10\t10\n"
	b := "chr1\t0\t100\n"
	got := runCoverage(t, a, b, Options{})
	// Even though the half-open overlap predicate gives zero overlap with
	// a zero-length interval, the *interval tree* may still surface the B
	// candidate; whatever it does, the output should be sane.
	if got == "" {
		t.Errorf("zero-length A: produced no output")
	}
}

func TestRecordColumns_BED3(t *testing.T) {
	a := "chr1\t0\t100\n"
	got := runCoverage(t, a, "", Options{})
	want := "chr1\t0\t100\t0\t0\t100\t0.0000000\n"
	if got != want {
		t.Errorf("BED3 cols: want %q, got %q", want, got)
	}
}

func TestRecordColumns_ExtraFields(t *testing.T) {
	// 13-field input → ExtraFields populated. Output should keep all fields.
	// We hand-craft a BED12 with one trailing extra column.
	a := "chr1\t0\t100\tn\t10\t+\t0\t100\t0,0,0\t1\t100,\t0,\textra\n"
	got := runCoverage(t, a, "", Options{Mode: ModeCounts})
	// Note: bed.Reader treats the 13th column as an ExtraField. BED12
	// blockSize / blockStart columns are emitted with the trailing
	// comma that upstream bedtools uses.
	want := "chr1\t0\t100\tn\t10\t+\t0\t100\t0,0,0\t1\t100,\t0,\textra\t0\n"
	if got != want {
		t.Errorf("13-field round-trip:\nwant: %q\ngot:  %q", want, got)
	}
}

func TestStrand_BothFlagsRejectMissingStrand(t *testing.T) {
	// A has no strand. -s should reject any B since A.Strand=="".
	a := "chr1\t0\t100\n"
	b := "chr1\t10\t50\tn\t0\t+\n"
	got := runCoverage(t, a, b, Options{SameStrand: true})
	want := "chr1\t0\t100\t0\t0\t100\t0.0000000\n"
	if got != want {
		t.Errorf("strand missing reject: want %q got %q", want, got)
	}
}
