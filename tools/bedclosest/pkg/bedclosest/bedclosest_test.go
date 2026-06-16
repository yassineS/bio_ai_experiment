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
	// With -d (unsigned distance): B[30,40) closer; gap = (30-20)+1 = 11.
	got, n, err := runClosest(t, a, b, Options{ReportDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	if n != 1 {
		t.Errorf("n=%d, want 1", n)
	}
	want := "chr1\t10\t20\tchr1\t30\t40\t11\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestClosestUnsignedDistanceAbsolute(t *testing.T) {
	// -d always reports the absolute (unsigned) distance even when B is upstream.
	a := "chr1\t100\t110\n"
	b := "chr1\t10\t20\n"
	got, _, err := runClosest(t, a, b, Options{ReportDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t100\t110\tchr1\t10\t20\t81\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestClosestSignedDistanceRef(t *testing.T) {
	// -D ref signs the distance: upstream (left) B is negative.
	a := "chr1\t100\t110\n"
	b := "chr1\t10\t20\n"
	got, _, err := runClosest(t, a, b, Options{ReportDistance: true, DistanceMode: DistanceSignedRef})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t100\t110\tchr1\t10\t20\t-81\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestClosestOverlap(t *testing.T) {
	a := "chr1\t10\t30\n"
	b := "chr1\t20\t25\n"
	got, _, err := runClosest(t, a, b, Options{ReportDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t10\t30\tchr1\t20\t25\t0\n"
	if got != want {
		t.Errorf("got %q", got)
	}
}

func TestClosestTiesAll(t *testing.T) {
	a := "chr1\t100\t110\n"
	// Two B's equidistant (both 11bp under bedtools' (gap+1) convention).
	b := "chr1\t80\t90\nchr1\t120\t130\n"
	got, n, err := runClosest(t, a, b, Options{ReportDistance: true, DistanceMode: DistanceSignedRef, TieBreak: TieAll})
	if err != nil {
		t.Fatal(err)
	}
	if n != 2 {
		t.Errorf("n=%d want 2", n)
	}
	// Upstream emits the upstream tie first, then the downstream tie.
	want := "chr1\t100\t110\tchr1\t80\t90\t-11\n" +
		"chr1\t100\t110\tchr1\t120\t130\t11\n"
	if got != want {
		t.Errorf("got %q\nwant %q", got, want)
	}
}

func TestClosestTiesFirstLast(t *testing.T) {
	a := "chr1\t100\t110\n"
	b := "chr1\t80\t90\nchr1\t120\t130\n"

	gotF, _, err := runClosest(t, a, b, Options{TieBreak: TieFirst})
	if err != nil {
		t.Fatal(err)
	}
	// -t first: upstream tie wins (the upstream group is consulted before the
	// downstream group when both tie and the mode is not "last").
	if !strings.Contains(gotF, "80\t90") || strings.Contains(gotF, "120\t130") {
		t.Errorf("TieFirst output unexpected: %q", gotF)
	}

	gotL, _, err := runClosest(t, a, b, Options{TieBreak: TieLast})
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
	got, _, err := runClosest(t, a, b, Options{})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t10\t20\tchr1\t30\t40\n"
	if got != want {
		t.Errorf("got %q", got)
	}
}

func TestClosestKClosest(t *testing.T) {
	a := "chr1\t100\t110\n"
	// Three downstream B's at increasing distance.
	b := "chr1\t120\t130\nchr1\t140\t150\nchr1\t160\t170\n"
	got, n, err := runClosest(t, a, b, Options{KClosest: 2})
	if err != nil {
		t.Fatal(err)
	}
	if n != 2 {
		t.Errorf("n=%d want 2", n)
	}
	want := "chr1\t100\t110\tchr1\t120\t130\n" +
		"chr1\t100\t110\tchr1\t140\t150\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestClosestIgnoreOverlaps(t *testing.T) {
	// -io drops the overlapping B; the next-closest non-overlapping B wins.
	a := "chr1\t100\t110\n"
	b := "chr1\t105\t115\nchr1\t200\t210\n"
	got, _, err := runClosest(t, a, b, Options{IgnoreOverlaps: true})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t100\t110\tchr1\t200\t210\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestClosestDifferentNames(t *testing.T) {
	// A "g1" sits between B "g1" (nearer) and B "g2" (farther). With -N the
	// same-named B is skipped, so the farther different-named B is reported.
	a := "chr1\t100\t110\tg1\n"
	b := "chr1\t115\t125\tg1\nchr1\t140\t150\tg2\n"

	got, _, err := runClosest(t, a, b, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(got, "\tg1\tchr1\t115\t125\tg1") {
		t.Errorf("without -N expected nearest same-name B, got %q", got)
	}

	gotN, _, err := runClosest(t, a, b, Options{DifferentNames: true})
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(gotN, "\tg1\tchr1\t140\t150\tg2") {
		t.Errorf("with -N expected the different-name B g2, got %q", gotN)
	}
}

func TestClosestDifferentChromNullShape(t *testing.T) {
	// No B on chr1: a BED3 null placeholder is emitted (no distance column).
	a := "chr1\t10\t20\n"
	b := "chr2\t10\t20\n"
	got, _, err := runClosest(t, a, b, Options{})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t10\t20\t.\t-1\t-1\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}

	// With -d the trailing distance -1 is appended.
	gotD, _, err := runClosest(t, a, b, Options{ReportDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	wantD := "chr1\t10\t20\t.\t-1\t-1\t-1\n"
	if gotD != wantD {
		t.Errorf("got %q want %q", gotD, wantD)
	}
}

func TestClosestNullShapesPerRecordType(t *testing.T) {
	a := "chrZ\t10\t20\n"
	cases := []struct {
		name string
		b    string
		want string
	}{
		{"bed3", "chr1\t100\t200\n", "chrZ\t10\t20\t.\t-1\t-1\n"},
		{"bed4", "chr1\t100\t200\tname\n", "chrZ\t10\t20\t.\t-1\t-1\t.\n"},
		{"bed5", "chr1\t100\t200\tname\t50\n", "chrZ\t10\t20\t.\t-1\t-1\t.\t-1\n"},
		{"bed6", "chr1\t100\t200\tname\t50\t+\n", "chrZ\t10\t20\t.\t-1\t-1\t.\t-1\t.\n"},
		{"bed12", "chr1\t100\t200\tn\t50\t+\t100\t200\t0\t1\t100,\t0,\n",
			"chrZ\t10\t20\t.\t-1\t-1\t.\t-1\t.\t.\t.\t.\t.\t.\t.\n"},
		{"bedgraph", "chr1\t100\t200\t3.5\n", "chrZ\t10\t20\t.\t-1\t-1\t.\n"},
		{"bed4plus", "chr1\t100\t200\tb1\tbx\n", "chrZ\t10\t20\t.\t-1\t-1\t.\t.\n"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got, _, err := runClosest(t, a, c.b, Options{})
			if err != nil {
				t.Fatal(err)
			}
			if got != c.want {
				t.Errorf("got %q want %q", got, c.want)
			}
		})
	}
}

func TestClosestStrandA(t *testing.T) {
	// A on '-' strand, B downstream on reference -> with -D a, sign flips.
	a := "chr1\t100\t110\tA\t0\t-\n"
	b := "chr1\t120\t130\tB\t0\t+\n"
	got, _, err := runClosest(t, a, b, Options{ReportDistance: true, DistanceMode: DistanceA})
	if err != nil {
		t.Fatal(err)
	}
	if !strings.HasSuffix(strings.TrimRight(got, "\n"), "\t-11") {
		t.Errorf("expected -D a to flip sign on minus-strand A; got %q", got)
	}
}

func TestClosestStrandB(t *testing.T) {
	a := "chr1\t100\t110\tA\t0\t+\n"
	b := "chr1\t120\t130\tB\t0\t-\n"
	got, _, err := runClosest(t, a, b, Options{ReportDistance: true, DistanceMode: DistanceB})
	if err != nil {
		t.Fatal(err)
	}
	// -D b flips only on a FORWARD-strand B; B here is '-', so no flip: +11.
	if !strings.HasSuffix(strings.TrimRight(got, "\n"), "\t11") {
		t.Errorf("expected -D b NOT to flip sign with B on minus strand; got %q", got)
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
	_, _, err := runClosest(t, "chr1\t0\t10\n", "chr1\tBAD\t10\n", Options{})
	if err == nil {
		t.Error("expected B parse error")
	}
}

func TestClosestSweepCorrectness(t *testing.T) {
	// A long B upstream that overlaps A; should still be picked up despite a
	// later B having a smaller gap.
	a := "chr1\t500\t510\n"
	b := "chr1\t0\t1000\nchr1\t600\t610\n"
	got, _, err := runClosest(t, a, b, Options{ReportDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t500\t510\tchr1\t0\t1000\t0\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestClosestSkipsCommentsAndBlankLines(t *testing.T) {
	a := "# comment\n\ntrack name=foo\nchr1\t10\t20\n"
	b := "chr1\t30\t40\n"
	got, _, err := runClosest(t, a, b, Options{ReportDistance: true})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t10\t20\tchr1\t30\t40\t11\n"
	if got != want {
		t.Errorf("got %q", got)
	}
}

func TestClosestHeader(t *testing.T) {
	a := "#Header for file a.bed\nchr1\t10\t20\n"
	b := "chr1\t20\t21\n"
	got, _, err := runClosest(t, a, b, Options{PrintHeader: true})
	if err != nil {
		t.Fatal(err)
	}
	want := "#Header for file a.bed\nchr1\t10\t20\tchr1\t20\t21\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestClosestClassify(t *testing.T) {
	aPlus := &Row{Chrom: "chr1", Start: 100, End: 110, StrandV: strandForward}
	aMinus := &Row{Chrom: "chr1", Start: 100, End: 110, StrandV: strandReverse}
	bDown := &Row{Chrom: "chr1", Start: 120, End: 130, StrandV: strandForward}
	bUp := &Row{Chrom: "chr1", Start: 50, End: 60, StrandV: strandReverse}

	if s, d := classify(aPlus, bDown, Options{DistanceMode: DistanceSignedRef}); s != streamDownstream || d != 11 {
		t.Errorf("ref downstream = (%v,%d), want (down,11)", s, d)
	}
	if s, d := classify(aPlus, bUp, Options{DistanceMode: DistanceSignedRef}); s != streamUpstream || d != 41 {
		t.Errorf("ref upstream = (%v,%d), want (up,41)", s, d)
	}
	// -D a flips on reverse-strand A: bDown becomes upstream.
	if s, _ := classify(aMinus, bDown, Options{DistanceMode: DistanceA}); s != streamUpstream {
		t.Errorf("a-strand flip: got %v want upstream", s)
	}
	// -D b flips on a FORWARD-strand B.
	bUpFwd := &Row{Chrom: "chr1", Start: 50, End: 60, StrandV: strandForward}
	if s, _ := classify(aPlus, bUpFwd, Options{DistanceMode: DistanceB}); s != streamDownstream {
		t.Errorf("b-strand (forward) flip: got %v want downstream", s)
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
	if _, _, err := ReadAll(strings.NewReader("chr1\t10\n")); err == nil {
		t.Error("expected error")
	}
}

func TestReadAllHeader(t *testing.T) {
	rows, header, err := ReadAll(strings.NewReader("#h1\n#h2\nchr1\t1\t2\n#mid\n"))
	if err != nil {
		t.Fatal(err)
	}
	if len(rows) != 1 {
		t.Errorf("rows=%d want 1", len(rows))
	}
	if len(header) != 2 || header[0] != "#h1" || header[1] != "#h2" {
		t.Errorf("header=%v want [#h1 #h2]", header)
	}
}

func TestClosestEmptyA(t *testing.T) {
	got, n, err := runClosest(t, "", "chr1\t10\t20\n", Options{ReportDistance: true})
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
	got, n, err := runClosest(t, a, b, Options{ReportDistance: true, DistanceMode: DistanceSignedRef})
	if err != nil {
		t.Fatal(err)
	}
	if n != 3 {
		t.Errorf("n=%d want 3", n)
	}
	want := "chr1\t10\t20\tchr1\t40\t45\t21\n" +
		"chr1\t50\t60\tchr1\t40\t45\t-6\n" +
		"chr1\t100\t110\tchr1\t90\t95\t-6\n"
	if got != want {
		t.Errorf("got %q\nwant %q", got, want)
	}
}

func TestClosestSameStrandFilterSelection(t *testing.T) {
	a := "chr1\t100\t110\ta\t0\t+\n"
	// B[120,130) on '-' is closest by distance but wrong strand; B[200,210) on
	// '+' is the closest eligible candidate under -s.
	b := "chr1\t120\t130\tnear\t0\t-\nchr1\t200\t210\tfar\t0\t+\n"
	got, _, err := runClosest(t, a, b, Options{ReportDistance: true, DistanceMode: DistanceSignedRef, SameStrand: true})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t100\t110\ta\t0\t+\tchr1\t200\t210\tfar\t0\t+\t91\n"
	if got != want {
		t.Errorf("got %q\nwant %q", got, want)
	}
}

func TestSameAndOppositeStrandMutuallyExclusive(t *testing.T) {
	opts := Options{SameStrand: true, OppositeStrand: true}
	if err := opts.Validate(); err == nil {
		t.Fatal("Validate: expected error for -s and -S together, got nil")
	}
	_, _, err := runClosest(t, "chr1\t10\t20\ta\t0\t+\n", "chr1\t30\t40\tb\t0\t+\n", opts)
	if err == nil {
		t.Fatal("Closest: expected error for -s and -S together, got nil")
	}
}

func TestStrandBad(t *testing.T) {
	tests := []struct {
		name     string
		aSV, bSV int
		opts     Options
		wantBad  bool
	}{
		{"no filter", strandForward, strandReverse, Options{}, false},
		{"same equal", strandForward, strandForward, Options{SameStrand: true}, false},
		{"same differ", strandForward, strandReverse, Options{SameStrand: true}, true},
		{"opp differ", strandForward, strandReverse, Options{OppositeStrand: true}, false},
		{"opp equal", strandForward, strandForward, Options{OppositeStrand: true}, true},
		{"same unknown A", strandUnknown, strandReverse, Options{SameStrand: true}, true},
		{"opp unknown B", strandForward, strandUnknown, Options{OppositeStrand: true}, true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			a := &Row{StrandV: tt.aSV}
			b := &Row{StrandV: tt.bSV}
			if got := strandBad(a, b, tt.opts); got != tt.wantBad {
				t.Errorf("strandBad = %v, want %v", got, tt.wantBad)
			}
		})
	}
}

func TestDetectRecordType(t *testing.T) {
	cases := []struct {
		fields []string
		want   RecordType
	}{
		{[]string{"c", "1", "2"}, recBed3},
		{[]string{"c", "1", "2", "n"}, recBed4},
		{[]string{"c", "1", "2", "3.5"}, recBedGraph},
		{[]string{"c", "1", "2", "n", "50"}, recBed5},
		{[]string{"c", "1", "2", "n", "50", "+"}, recBed6},
		{[]string{"c", "1", "2", "n", "50", "+", "1", "2", "0", "1", "1,", "0,"}, recBed12},
		{[]string{"c", "1", "2", "n", "x"}, recBedPlus},
	}
	for _, c := range cases {
		rt, _ := detectRecordType(&Row{Fields: c.fields})
		if rt != c.want {
			t.Errorf("detectRecordType(%v) = %v want %v", c.fields, rt, c.want)
		}
	}
}
