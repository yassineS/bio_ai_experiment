package bedflank

import (
	"bytes"
	"strings"
	"testing"
)

func runFlank(t *testing.T, in string, sizes ChromSizes, opts Options) (string, string, int, error) {
	t.Helper()
	var out, warn bytes.Buffer
	n, err := Flank(strings.NewReader(in), &out, &warn, sizes, opts)
	return out.String(), warn.String(), n, err
}

func TestFlankBoth(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	in := "chr1\t100\t200\n"
	got, _, n, err := runFlank(t, in, sizes, Options{Both: true, BothAdd: 50})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t50\t100\nchr1\t200\t250\n"
	if got != want {
		t.Errorf("want %q got %q", want, got)
	}
	if n != 2 {
		t.Errorf("count %d, want 2", n)
	}
}

func TestFlankAsymmetric(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	in := "chr1\t100\t200\n"
	got, _, n, err := runFlank(t, in, sizes, Options{LeftAdd: 20, RightAdd: 5})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t80\t100\nchr1\t200\t205\n"
	if got != want || n != 2 {
		t.Errorf("got %q n=%d", got, n)
	}
}

func TestFlankClippedAtZero(t *testing.T) {
	sizes := ChromSizes{"chr1": 100}
	in := "chr1\t5\t50\n"
	got, _, _, err := runFlank(t, in, sizes, Options{Both: true, BothAdd: 20})
	if err != nil {
		t.Fatal(err)
	}
	// Left flank clipped at 0; right flank up to 70.
	want := "chr1\t0\t5\nchr1\t50\t70\n"
	if got != want {
		t.Errorf("got %q", got)
	}
}

func TestFlankClippedAtChromEnd(t *testing.T) {
	sizes := ChromSizes{"chr1": 100}
	in := "chr1\t60\t90\n"
	got, _, _, err := runFlank(t, in, sizes, Options{Both: true, BothAdd: 30})
	if err != nil {
		t.Fatal(err)
	}
	// Left: [30, 60). Right: [90, 100).
	want := "chr1\t30\t60\nchr1\t90\t100\n"
	if got != want {
		t.Errorf("got %q", got)
	}
}

func TestFlankSkipsEmpty(t *testing.T) {
	sizes := ChromSizes{"chr1": 100}
	// At position 0, left flank is empty. At position 100, right flank empty.
	in := "chr1\t0\t10\nchr1\t90\t100\n"
	got, _, n, err := runFlank(t, in, sizes, Options{Both: true, BothAdd: 5})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t10\t15\nchr1\t85\t90\n"
	if got != want || n != 2 {
		t.Errorf("got %q n=%d", got, n)
	}
}

func TestFlankStrandSwapsLeftRight(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	// Same input, with strand-specific flank: -l 10 -r 100. On '+' strand the
	// left flank is 10bp, right is 100bp. On '-' strand they're swapped.
	in := "chr1\t100\t200\tg+\t0\t+\nchr1\t100\t200\tg-\t0\t-\n"
	got, _, _, err := runFlank(t, in, sizes, Options{LeftAdd: 10, RightAdd: 100, StrandSpec: true})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t90\t100\tg+\t0\t+\n" +
		"chr1\t200\t300\tg+\t0\t+\n" +
		"chr1\t0\t100\tg-\t0\t-\n" +
		"chr1\t200\t210\tg-\t0\t-\n"
	if got != want {
		t.Errorf("got %q\nwant %q", got, want)
	}
}

func TestFlankPct(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	// Interval length 100, pct=0.1 -> 10bp each side.
	in := "chr1\t100\t200\n"
	got, _, _, err := runFlank(t, in, sizes, Options{Both: true, BothAdd: 0.1, Pct: true})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t90\t100\nchr1\t200\t210\n"
	if got != want {
		t.Errorf("got %q", got)
	}
}

func TestFlankUnknownChrom(t *testing.T) {
	sizes := ChromSizes{"chr1": 100}
	in := "chr1\t10\t20\nchrX\t10\t20\n"
	got, warn, n, err := runFlank(t, in, sizes, Options{Both: true, BothAdd: 5})
	if err != nil {
		t.Fatal(err)
	}
	if n != 2 {
		t.Errorf("expected 2 records emitted (only chr1), got %d", n)
	}
	if !strings.Contains(warn, "chrX") {
		t.Errorf("expected warning about chrX, got %q", warn)
	}
	if strings.Contains(got, "chrX") {
		t.Errorf("chrX should not appear in output: %q", got)
	}
}

func TestFlankErrors(t *testing.T) {
	sizes := ChromSizes{"chr1": 100}
	cases := []string{
		"chr1\t10\n",       // too few fields
		"chr1\tBAD\t20\n",  // bad start
		"chr1\t10\tNOPE\n", // bad end
		"chr1\t50\t10\n",   // end<start
	}
	for _, in := range cases {
		_, err := Flank(strings.NewReader(in), &bytes.Buffer{}, &bytes.Buffer{}, sizes, Options{Both: true, BothAdd: 5})
		if err == nil {
			t.Errorf("expected error for input %q", in)
		}
	}
}

func TestFlankSkipsHeadersAndBlankLines(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	in := "# comment\n\ntrack name=foo\nbrowser hide all\nchr1\t100\t200\n"
	got, _, n, err := runFlank(t, in, sizes, Options{Both: true, BothAdd: 5})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t95\t100\nchr1\t200\t205\n"
	if got != want || n != 2 {
		t.Errorf("got %q n=%d", got, n)
	}
}

func TestReadChromSizes(t *testing.T) {
	in := "chr1\t100\nchr2\t200\n# comment\n\nchr3 300\n"
	sizes, err := ReadChromSizes(strings.NewReader(in))
	if err != nil {
		t.Fatal(err)
	}
	want := ChromSizes{"chr1": 100, "chr2": 200, "chr3": 300}
	for k, v := range want {
		if sizes[k] != v {
			t.Errorf("%s = %d, want %d", k, sizes[k], v)
		}
	}
}

func TestReadChromSizesErrors(t *testing.T) {
	cases := []string{
		"chr1\n",      // missing size
		"chr1\tBAD\n", // non-numeric size
		"chr1\t-5\n",  // negative size
	}
	for _, c := range cases {
		if _, err := ReadChromSizes(strings.NewReader(c)); err == nil {
			t.Errorf("expected error for %q", c)
		}
	}
}

func TestFlankNegativeValuesClamped(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	in := "chr1\t100\t200\n"
	// Negative flank requests are clamped to 0 (no flank emitted).
	got, _, n, err := runFlank(t, in, sizes, Options{LeftAdd: -10, RightAdd: 20})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t200\t220\n"
	if got != want || n != 1 {
		t.Errorf("got %q n=%d", got, n)
	}
}
