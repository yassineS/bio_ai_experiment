package bedclosest

import (
	"bytes"
	"strings"
	"testing"
)

func runClosest(t *testing.T, a, b string, opts Options) (string, int, error) {
	t.Helper()
	var buf bytes.Buffer
	n, err := Closest(strings.NewReader(a), strings.NewReader(b), &buf, opts)
	return buf.String(), n, err
}

func TestClosestBasic(t *testing.T) {
	a := "chr1\t10\t20\n"
	b := "chr1\t30\t40\nchr1\t100\t200\n"
	opts := Options{PrintDistance: true}
	got, n, err := runClosest(t, a, b, opts)
	if err != nil {
		t.Fatal(err)
	}
	if n != 1 {
		t.Errorf("n=%d, want 1", n)
	}
	// B[30,40) is closer; gap = 30 - 20 = 10, positive (downstream).
	want := "chr1\t10\t20\tchr1\t30\t40\t10\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestClosestOverlap(t *testing.T) {
	a := "chr1\t10\t30\n"
	b := "chr1\t20\t25\n"
	got, _, err := runClosest(t, a, b, Options{PrintDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t10\t30\tchr1\t20\t25\t0\n"
	if got != want {
		t.Errorf("got %q", got)
	}
}

func TestClosestUpstreamNegative(t *testing.T) {
	a := "chr1\t100\t200\n"
	b := "chr1\t10\t20\n"
	got, _, err := runClosest(t, a, b, Options{PrintDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	// B is upstream of A: gap = 100-20 = 80, sign negative.
	want := "chr1\t100\t200\tchr1\t10\t20\t-80\n"
	if got != want {
		t.Errorf("got %q", got)
	}
}

func TestClosestTiesAll(t *testing.T) {
	a := "chr1\t100\t110\n"
	// Two B's equidistant (both 10bp away).
	b := "chr1\t80\t90\nchr1\t120\t130\n"
	got, n, err := runClosest(t, a, b, Options{PrintDistance: true, TieBreak: TieAll})
	if err != nil {
		t.Fatal(err)
	}
	if n != 2 {
		t.Errorf("n=%d want 2", n)
	}
	want := "chr1\t100\t110\tchr1\t80\t90\t-10\n" +
		"chr1\t100\t110\tchr1\t120\t130\t10\n"
	if got != want {
		t.Errorf("got %q\nwant %q", got, want)
	}
}

func TestClosestTiesFirstLast(t *testing.T) {
	a := "chr1\t100\t110\n"
	b := "chr1\t80\t90\nchr1\t120\t130\n"

	gotF, _, err := runClosest(t, a, b, Options{PrintDistance: true, TieBreak: TieFirst})
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(gotF, "80\t90") || strings.Contains(gotF, "120\t130") {
		t.Errorf("TieFirst output unexpected: %q", gotF)
	}

	gotL, _, err := runClosest(t, a, b, Options{PrintDistance: true, TieBreak: TieLast})
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(gotL, "120\t130") || strings.Contains(gotL, "80\t90") {
		t.Errorf("TieLast output unexpected: %q", gotL)
	}
}

func TestClosestNoDistanceColumn(t *testing.T) {
	a := "chr1\t10\t20\n"
	b := "chr1\t30\t40\n"
	got, _, err := runClosest(t, a, b, Options{PrintDistance: false})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t10\t20\tchr1\t30\t40\n"
	if got != want {
		t.Errorf("got %q", got)
	}
}

func TestClosestRequireOverlap(t *testing.T) {
	a := "chr1\t10\t20\nchr1\t100\t200\n"
	b := "chr1\t150\t160\n"
	got, n, err := runClosest(t, a, b, Options{PrintDistance: true, RequireOverlap: true})
	if err != nil {
		t.Fatal(err)
	}
	// Only the second A overlaps B; the first should be omitted.
	if n != 1 {
		t.Errorf("n=%d want 1", n)
	}
	want := "chr1\t100\t200\tchr1\t150\t160\t0\n"
	if got != want {
		t.Errorf("got %q", got)
	}
}

func TestClosestDifferentChrom(t *testing.T) {
	a := "chr1\t10\t20\n"
	b := "chr2\t10\t20\n"
	got, _, err := runClosest(t, a, b, Options{PrintDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	// No B on chr1: missing sentinel emitted.
	want := "chr1\t10\t20\t.\t-1\t-1\t-1\n"
	if got != want {
		t.Errorf("got %q", got)
	}
}

func TestClosestDifferentChromRequireOverlap(t *testing.T) {
	a := "chr1\t10\t20\n"
	b := "chr2\t10\t20\n"
	got, n, err := runClosest(t, a, b, Options{PrintDistance: true, RequireOverlap: true})
	if err != nil {
		t.Fatal(err)
	}
	if n != 0 || got != "" {
		t.Errorf("got %q n=%d", got, n)
	}
}

func TestClosestStrandA(t *testing.T) {
	// A on '-' strand, B downstream on reference -> with -D a, sign flips.
	a := "chr1\t100\t110\tA\t0\t-\n"
	b := "chr1\t120\t130\tB\t0\t+\n"
	got, _, err := runClosest(t, a, b, Options{PrintDistance: true, DistanceMode: DistanceA})
	if err != nil {
		t.Fatal(err)
	}
	// Ref-signed distance is +10; with -D a and A on '-', flips to -10.
	if !strings.HasSuffix(strings.TrimRight(got, "\n"), "\t-10") {
		t.Errorf("expected -D a to flip sign on minus-strand A; got %q", got)
	}
}

func TestClosestStrandB(t *testing.T) {
	a := "chr1\t100\t110\tA\t0\t+\n"
	b := "chr1\t120\t130\tB\t0\t-\n"
	got, _, err := runClosest(t, a, b, Options{PrintDistance: true, DistanceMode: DistanceB})
	if err != nil {
		t.Fatal(err)
	}
	if !strings.HasSuffix(strings.TrimRight(got, "\n"), "\t-10") {
		t.Errorf("expected -D b to flip sign with B on minus strand; got %q", got)
	}
}

func TestClosestUnsortedError(t *testing.T) {
	a := "chr1\t100\t200\nchr1\t50\t60\n" // not sorted
	b := "chr1\t10\t20\n"
	_, _, err := runClosest(t, a, b, Options{})
	if err == nil {
		t.Error("expected sort error for A")
	}

	bUnsorted := "chr1\t100\t200\nchr1\t50\t60\n"
	_, _, err = runClosest(t, "chr1\t10\t20\n", bUnsorted, Options{})
	if err == nil {
		t.Error("expected sort error for B")
	}
}

func TestClosestParseError(t *testing.T) {
	cases := []string{
		"chr1\t10\n",       // too few fields
		"chr1\tBAD\t20\n",  // bad start
		"chr1\t10\tNOPE\n", // bad end
		"chr1\t50\t10\n",   // end<start
	}
	for _, c := range cases {
		_, _, err := runClosest(t, c, "", Options{})
		if err == nil {
			t.Errorf("expected error for %q", c)
		}
	}
	// Bad B should error too.
	_, _, err := runClosest(t, "chr1\t0\t10\n", "chr1\tBAD\t10\n", Options{})
	if err == nil {
		t.Error("expected B parse error")
	}
}

func TestClosestSweepCorrectness(t *testing.T) {
	// A long B upstream that overlaps A; should still be picked up despite
	// later B's having smaller Start.
	a := "chr1\t500\t510\n"
	b := "chr1\t0\t1000\nchr1\t600\t610\n"
	got, _, err := runClosest(t, a, b, Options{PrintDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	// The long B fully covers A -> dist 0. It's the unique closest.
	want := "chr1\t500\t510\tchr1\t0\t1000\t0\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestClosestSkipsCommentsAndBlankLines(t *testing.T) {
	a := "# comment\n\ntrack name=foo\nchr1\t10\t20\n"
	b := "chr1\t30\t40\n"
	got, _, err := runClosest(t, a, b, Options{PrintDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t10\t20\tchr1\t30\t40\t10\n"
	if got != want {
		t.Errorf("got %q", got)
	}
}

func TestSignedDistanceModes(t *testing.T) {
	aPlus := &Row{Chrom: "chr1", Start: 100, End: 110, Strand: "+"}
	aMinus := &Row{Chrom: "chr1", Start: 100, End: 110, Strand: "-"}
	bDown := &Row{Chrom: "chr1", Start: 120, End: 130, Strand: "+"}
	bUp := &Row{Chrom: "chr1", Start: 50, End: 60, Strand: "-"}

	if d := signedDistance(aPlus, bDown, Options{}); d != 10 {
		t.Errorf("ref downstream = %d, want 10", d)
	}
	if d := signedDistance(aPlus, bUp, Options{}); d != -40 {
		t.Errorf("ref upstream = %d, want -40", d)
	}
	if d := signedDistance(aMinus, bDown, Options{DistanceMode: DistanceA}); d != -10 {
		t.Errorf("a-strand flip = %d, want -10", d)
	}
	if d := signedDistance(aPlus, bUp, Options{DistanceMode: DistanceB}); d != 40 {
		t.Errorf("b-strand flip = %d, want 40", d)
	}
}

func TestCheckSorted(t *testing.T) {
	rows := []*Row{
		{Chrom: "chr1", Start: 0, End: 10},
		{Chrom: "chr1", Start: 5, End: 15},
		{Chrom: "chr2", Start: 0, End: 10}, // new chrom is OK
	}
	if err := CheckSorted(rows, "X"); err != nil {
		t.Errorf("expected no error: %v", err)
	}
	bad := []*Row{
		{Chrom: "chr1", Start: 10, End: 20},
		{Chrom: "chr1", Start: 5, End: 6},
	}
	if err := CheckSorted(bad, "X"); err == nil {
		t.Error("expected sort error")
	}
}

func TestReadAllError(t *testing.T) {
	if _, err := ReadAll(strings.NewReader("chr1\t10\n")); err == nil {
		t.Error("expected error")
	}
}

func TestClosestEmptyA(t *testing.T) {
	got, n, err := runClosest(t, "", "chr1\t10\t20\n", Options{PrintDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	if n != 0 || got != "" {
		t.Errorf("got %q n=%d", got, n)
	}
}

func TestClosestMultipleAOnSameChrom(t *testing.T) {
	a := "chr1\t10\t20\nchr1\t50\t60\nchr1\t100\t110\n"
	b := "chr1\t40\t45\nchr1\t90\t95\n"
	got, n, err := runClosest(t, a, b, Options{PrintDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	if n != 3 {
		t.Errorf("n=%d want 3", n)
	}
	// A[10,20) closest: B[40,45) at distance 20 (only one closer; B[90,95) is 70 away).
	// A[50,60) closest: B[40,45) at -5 (gap 5 upstream) and B[90,95) at 30; -> [40,45).
	// A[100,110) closest: B[90,95) at -5.
	want := "chr1\t10\t20\tchr1\t40\t45\t20\n" +
		"chr1\t50\t60\tchr1\t40\t45\t-5\n" +
		"chr1\t100\t110\tchr1\t90\t95\t-5\n"
	if got != want {
		t.Errorf("got %q\nwant %q", got, want)
	}
}
