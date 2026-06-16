package bedintersect

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// runMulti runs IntersectMulti over string A and one-or-more string B files,
// defaulting FilePaths to numbered labels so the multi-DB column logic engages.
func runMulti(t *testing.T, a string, bs []string, opts IntersectOptions) string {
	t.Helper()
	readers := make([]io.Reader, len(bs))
	for i, b := range bs {
		readers[i] = strings.NewReader(b)
	}
	if opts.FilePaths == nil {
		paths := make([]string, len(bs))
		for i := range bs {
			paths[i] = "f" + string(rune('1'+i))
		}
		opts.FilePaths = paths
	}
	var buf bytes.Buffer
	if _, err := IntersectMulti(strings.NewReader(a), readers, &buf, opts); err != nil {
		t.Fatalf("IntersectMulti: %v", err)
	}
	return buf.String()
}

// run is the single-B-file convenience wrapper.
func run(t *testing.T, a, b string, opts IntersectOptions) string {
	t.Helper()
	opts.FilePaths = []string{"b"}
	return runMulti(t, a, []string{b}, opts)
}

func TestCountEachSingleAndMulti(t *testing.T) {
	a := "chr1\t10\t20\ta1\nchr1\t100\t200\ta2\n"
	b1 := "chr1\t100\t110\nchr1\t150\t160\n"
	b2 := "chr1\t10\t15\n"

	// Single B: -C prints A then the count (no DB-id column).
	got := run(t, a, b1, IntersectOptions{CountEach: true, FilePaths: []string{"b"}})
	want := "chr1\t10\t20\ta1\t0\nchr1\t100\t200\ta2\t2\n"
	if got != want {
		t.Fatalf("single -C:\n got %q\nwant %q", got, want)
	}

	// Two B files: -C prints one line per file with the numeric DB-id column.
	got = runMulti(t, a, []string{b1, b2}, IntersectOptions{CountEach: true})
	want = "chr1\t10\t20\ta1\t1\t0\n" +
		"chr1\t10\t20\ta1\t2\t1\n" +
		"chr1\t100\t200\ta2\t1\t2\n" +
		"chr1\t100\t200\ta2\t2\t0\n"
	if got != want {
		t.Fatalf("multi -C:\n got %q\nwant %q", got, want)
	}
}

func TestUnique(t *testing.T) {
	a := "chr1\t10\t20\ta1\nchr1\t100\t200\ta2\n"
	b := "chr1\t100\t110\nchr1\t150\t160\n"
	got := run(t, a, b, IntersectOptions{Unique: true})
	want := "chr1\t100\t200\ta2\n" // a2 overlaps twice but is reported once
	if got != want {
		t.Fatalf("-u:\n got %q\nwant %q", got, want)
	}
}

func TestForceOppositeStrand(t *testing.T) {
	a := "chr1\t10\t20\ta1\t0\t+\nchr1\t100\t200\ta2\t0\t-\n"
	b := "chr1\t10\t20\tb1\t0\t-\nchr1\t100\t200\tb2\t0\t+\n"
	// -S: only opposite-strand hits. a1(+)/b1(-) and a2(-)/b2(+) both qualify.
	got := run(t, a, b, IntersectOptions{ForceOpposite: true, WriteA: true})
	want := "chr1\t10\t20\ta1\t0\t+\nchr1\t100\t200\ta2\t0\t-\n"
	if got != want {
		t.Fatalf("-S:\n got %q\nwant %q", got, want)
	}
	// -s on the same data yields nothing (no same-strand overlaps here).
	got = run(t, a, b, IntersectOptions{StrandSpec: true, WriteA: true})
	if got != "" {
		t.Fatalf("-s should be empty, got %q", got)
	}
}

func TestEitherFraction(t *testing.T) {
	// A is 0..100 (len 100), B is 40..60 (len 20). Overlap 20.
	// fA = 20/100 = 0.2, fB = 20/20 = 1.0.
	a := "chr1\t0\t100\n"
	b := "chr1\t40\t60\n"
	// Without -e: -f 0.5 fails (0.2 < 0.5) AND -F 0.5 passes -> both required -> no hit.
	if got := run(t, a, b, IntersectOptions{FractionA: 0.5, FractionB: 0.5, WriteA: true}); got != "" {
		t.Fatalf("both-required should be empty, got %q", got)
	}
	// With -e: either suffices -> -F passes -> hit.
	got := run(t, a, b, IntersectOptions{FractionA: 0.5, FractionB: 0.5, EitherFraction: true, WriteA: true})
	if got != "chr1\t0\t100\n" {
		t.Fatalf("-e should report A, got %q", got)
	}
}

func TestVCFStructuralVariantEnd(t *testing.T) {
	// A SV deletion whose span comes from SVLEN, and a plain SNV (REF length).
	a := "##fileformat=VCFv4.1\n" +
		"19\t100\tsv1\tG\t<DEL>\t.\t.\tSVLEN=-50;END=149\n"
	// B overlaps only the SV span (100-1 .. 100-1+50 = 99..149) at 120..130.
	b := "19\t120\t130\n"
	got := run(t, a, b, IntersectOptions{WriteA: true})
	if !strings.Contains(got, "sv1") {
		t.Fatalf("SV deletion should overlap B via SVLEN span; got %q", got)
	}
	// An insertion is zero-length, so the same B does not overlap.
	aIns := "##fileformat=VCFv4.1\n19\t100\tins1\tG\t<INS>\t.\t.\tSVLEN=50\n"
	if got := run(t, aIns, b, IntersectOptions{WriteA: true}); got != "" {
		t.Fatalf("insertion is zero-length and should not overlap; got %q", got)
	}
}

func TestSortedValidationError(t *testing.T) {
	// Same chromosome, decreasing start -> out-of-order error under -sorted.
	a := "chr1\t100\t200\nchr1\t50\t60\n"
	b := "chr1\t10\t20\n"
	var buf bytes.Buffer
	_, err := IntersectMulti(strings.NewReader(a), []io.Reader{strings.NewReader(b)}, &buf,
		IntersectOptions{Sorted: true, NameA: "a.bed", FilePaths: []string{"b.bed"}})
	if err == nil || !IsSortError(err) {
		t.Fatalf("expected sort error, got %v", err)
	}
	if !strings.Contains(err.Error(), "out of order record") || !strings.Contains(err.Error(), "chr1\t50\t60") {
		t.Fatalf("sort error message wrong: %q", err.Error())
	}
}

func TestSortedValidationOKWhenSorted(t *testing.T) {
	a := "chr1\t10\t20\nchr1\t30\t40\n"
	b := "chr1\t15\t35\n"
	var buf bytes.Buffer
	if _, err := IntersectMulti(strings.NewReader(a), []io.Reader{strings.NewReader(b)}, &buf,
		IntersectOptions{Sorted: true, WriteA: true, NameA: "a", FilePaths: []string{"b"}}); err != nil {
		t.Fatalf("sorted-but-valid should not error: %v", err)
	}
	if buf.String() != "chr1\t10\t20\nchr1\t30\t40\n" {
		t.Fatalf("sorted output wrong: %q", buf.String())
	}
}

func TestHeaderEcho(t *testing.T) {
	a := "# a comment header\ntrack name=foo\nchr1\t10\t20\n"
	b := "chr1\t15\t25\n"
	// With -header the leading header lines are echoed verbatim before results.
	got := run(t, a, b, IntersectOptions{Header: true, WriteA: true})
	want := "# a comment header\ntrack name=foo\nchr1\t10\t20\n"
	if got != want {
		t.Fatalf("-header:\n got %q\nwant %q", got, want)
	}
	// Without -header the header is dropped.
	got = run(t, a, b, IntersectOptions{WriteA: true})
	if got != "chr1\t10\t20\n" {
		t.Fatalf("no -header should drop header, got %q", got)
	}
}

func TestNameConventionWarningMessage(t *testing.T) {
	// A uses chr1 (no leading zero), B uses chr01 (leading zero): inconsistent.
	a := "chr1\t10\t20\n"
	b := "chr01\t15\t25\n"
	var warn bytes.Buffer
	var out bytes.Buffer
	_, err := IntersectMulti(strings.NewReader(a), []io.Reader{strings.NewReader(b)}, &out,
		IntersectOptions{NameA: "a.bed", FilePaths: []string{"b.bed"}, Warnings: &warn})
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !strings.Contains(warn.String(), "naming convention (leading zero) is inconsistent") {
		t.Fatalf("expected naming-convention warning, got %q", warn.String())
	}
	// -nonamecheck suppresses it.
	warn.Reset()
	out.Reset()
	_, _ = IntersectMulti(strings.NewReader(a), []io.Reader{strings.NewReader(b)}, &out,
		IntersectOptions{NameA: "a.bed", FilePaths: []string{"b.bed"}, NoNameCheck: true, Warnings: &warn})
	if warn.Len() != 0 {
		t.Fatalf("-nonamecheck should suppress the warning, got %q", warn.String())
	}
}

func TestSortOutOrdersHitsAcrossFiles(t *testing.T) {
	// One A record overlapping hits from two B files at scattered positions.
	a := "chr1\t100\t300\ta1\n"
	b1 := "chr1\t250\t260\tz1\n"
	b2 := "chr1\t150\t160\ty1\n"
	// Without -sortout: hits grouped by file (file 1 then file 2).
	got := runMulti(t, a, []string{b1, b2}, IntersectOptions{WriteA: true, WriteB: true})
	if !strings.HasPrefix(got, "chr1\t100\t300\ta1\t1\tchr1\t250\t260\tz1\n") {
		t.Fatalf("default multi order wrong:\n%s", got)
	}
	// With -sortout: hits sorted by position across files (y1@150 before z1@250).
	got = runMulti(t, a, []string{b1, b2}, IntersectOptions{WriteA: true, WriteB: true, SortOut: true})
	want := "chr1\t100\t300\ta1\t2\tchr1\t150\t160\ty1\n" +
		"chr1\t100\t300\ta1\t1\tchr1\t250\t260\tz1\n"
	if got != want {
		t.Fatalf("-sortout order wrong:\n got %q\nwant %q", got, want)
	}
}

func TestFieldCountErrorMessage(t *testing.T) {
	// First line locks 3 fields; a later line with a trailing tab has 4.
	a := "chr1\t10\t20\nchr1\t30\t40\t\n"
	b := "chr1\t15\t25\n"
	var buf bytes.Buffer
	_, err := IntersectMulti(strings.NewReader(a), []io.Reader{strings.NewReader(b)}, &buf,
		IntersectOptions{FilePaths: []string{"b"}})
	if err == nil {
		t.Fatal("expected a field-count error")
	}
	if !IsVerbatimError(err) || !strings.Contains(VerbatimMessage(err), "wrong number of fields") {
		t.Fatalf("field-count error wrong: %v", err)
	}
}
