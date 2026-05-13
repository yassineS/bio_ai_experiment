package bedslop

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

func runSlop(t *testing.T, input string, sizes ChromSizes, opts Options) (string, string) {
	t.Helper()
	var out, warn bytes.Buffer
	if _, err := Slop(strings.NewReader(input), &out, &warn, sizes, opts); err != nil {
		t.Fatalf("Slop returned error: %v", err)
	}
	return out.String(), warn.String()
}

func TestSlopSymmetricBoth(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\t200\n"
	got, _ := runSlop(t, input, sizes, Options{Both: true, BothAdd: 10})
	want := "chr1\t90\t210\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestSlopAsymmetricLR(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\t200\n"
	got, _ := runSlop(t, input, sizes, Options{LeftAdd: 5, RightAdd: 20})
	want := "chr1\t95\t220\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestSlopClipsToChromBounds(t *testing.T) {
	sizes := ChromSizes{"chr1": 150}
	// Slop by 200 on each side should clip to [0, 150].
	input := "chr1\t100\t130\n"
	got, _ := runSlop(t, input, sizes, Options{Both: true, BothAdd: 200})
	want := "chr1\t0\t150\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestSlopNegativeShrinks(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\t200\n"
	got, _ := runSlop(t, input, sizes, Options{Both: true, BothAdd: -10})
	want := "chr1\t110\t190\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestSlopDropsEmptyInterval(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	// 100..120 shrunk by 50 each side -> 150..70 -> empty, drop.
	input := "chr1\t100\t120\n"
	out, warn := runSlop(t, input, sizes, Options{Both: true, BothAdd: -50})
	if out != "" {
		t.Errorf("expected no output, got %q", out)
	}
	if !strings.Contains(warn, "empty interval") {
		t.Errorf("expected warning about empty interval, got %q", warn)
	}
}

func TestSlopDropsUnknownChrom(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t10\t20\nchrZ\t10\t20\n"
	out, warn := runSlop(t, input, sizes, Options{Both: true, BothAdd: 5})
	if !strings.Contains(out, "chr1\t5\t25") {
		t.Errorf("missing chr1 output, got %q", out)
	}
	if strings.Contains(out, "chrZ") {
		t.Errorf("chrZ should have been dropped, got %q", out)
	}
	if !strings.Contains(warn, "chrZ") {
		t.Errorf("expected warning mentioning chrZ, got %q", warn)
	}
}

func TestSlopStrandSpecSwapsLeftRight(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\t200\tg1\t0\t+\n" +
		"chr1\t100\t200\tg2\t0\t-\n"
	got, _ := runSlop(t, input, sizes, Options{LeftAdd: 10, RightAdd: 100, StrandSpec: true})
	// '+' strand: extend chromStart by 10 -> 90, chromEnd by 100 -> 300
	// '-' strand: swap; "left" on minus strand is at the high end, so extend
	// chromStart by 100 -> 0, chromEnd by 10 -> 210.
	want := "chr1\t90\t300\tg1\t0\t+\n" +
		"chr1\t0\t210\tg2\t0\t-\n"
	if got != want {
		t.Errorf("strand swap mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}

func TestSlopStrandSpecIgnoredWithoutFlag(t *testing.T) {
	// Without StrandSpec, the strand column is ignored.
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\t200\tg1\t0\t-\n"
	got, _ := runSlop(t, input, sizes, Options{LeftAdd: 10, RightAdd: 100})
	want := "chr1\t90\t300\tg1\t0\t-\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestSlopPercentageMode(t *testing.T) {
	sizes := ChromSizes{"chr1": 10000}
	// 100bp interval, -b 0.5 -> add 50 on each side -> [50, 250]
	input := "chr1\t100\t200\n"
	got, _ := runSlop(t, input, sizes, Options{Both: true, BothAdd: 0.5, Pct: true})
	want := "chr1\t50\t250\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestSlopPercentageRounding(t *testing.T) {
	sizes := ChromSizes{"chr1": 10000}
	// 10bp interval, -b 0.25 -> 2.5, rounds to 2 or 3 (math.Round rounds to even? actually away from zero in Go).
	input := "chr1\t100\t110\n"
	got, _ := runSlop(t, input, sizes, Options{Both: true, BothAdd: 0.25, Pct: true})
	// math.Round(2.5) -> 3 in Go (round half away from zero).
	want := "chr1\t97\t113\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestSlopPreservesExtraColumns(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\t200\tname\t900\t+\textra1\textra2\n"
	got, _ := runSlop(t, input, sizes, Options{Both: true, BothAdd: 5})
	want := "chr1\t95\t205\tname\t900\t+\textra1\textra2\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestSlopSkipsHeadersAndBlanks(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "# comment\n" +
		"track name=foo\n" +
		"browser hide all\n" +
		"\n" +
		"chr1\t100\t200\n"
	got, _ := runSlop(t, input, sizes, Options{Both: true, BothAdd: 5})
	want := "chr1\t95\t205\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestSlopErrorsOnTooFewFields(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\n"
	var out, warn bytes.Buffer
	_, err := Slop(strings.NewReader(input), &out, &warn, sizes, Options{Both: true, BothAdd: 5})
	if err == nil {
		t.Errorf("expected error on too-few fields")
	}
}

func TestSlopErrorsOnBadStart(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\tNOPE\t200\n"
	var out, warn bytes.Buffer
	_, err := Slop(strings.NewReader(input), &out, &warn, sizes, Options{Both: true, BothAdd: 5})
	if err == nil {
		t.Errorf("expected error on bad chromStart")
	}
}

func TestSlopErrorsOnBadEnd(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\tNOPE\n"
	var out, warn bytes.Buffer
	_, err := Slop(strings.NewReader(input), &out, &warn, sizes, Options{Both: true, BothAdd: 5})
	if err == nil {
		t.Errorf("expected error on bad chromEnd")
	}
}

func TestSlopNilWarnerOK(t *testing.T) {
	// Passing a nil warn writer must not panic for warning cases.
	sizes := ChromSizes{"chr1": 1000}
	input := "chrZ\t10\t20\n"
	var out bytes.Buffer
	if _, err := Slop(strings.NewReader(input), &out, io.Discard, sizes, Options{Both: true, BothAdd: 5}); err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
}

func TestReadChromSizes(t *testing.T) {
	input := "# comment\n\nchr1\t1000\nchr2\t500\nchr3 200\n"
	sizes, err := ReadChromSizes(strings.NewReader(input))
	if err != nil {
		t.Fatalf("ReadChromSizes: %v", err)
	}
	if sizes["chr1"] != 1000 || sizes["chr2"] != 500 || sizes["chr3"] != 200 {
		t.Errorf("unexpected sizes: %v", sizes)
	}
}

func TestReadChromSizesBadFormat(t *testing.T) {
	if _, err := ReadChromSizes(strings.NewReader("onlyname\n")); err == nil {
		t.Errorf("expected error for single-field line")
	}
	if _, err := ReadChromSizes(strings.NewReader("chr1\tNOT_A_NUMBER\n")); err == nil {
		t.Errorf("expected error for non-numeric size")
	}
	if _, err := ReadChromSizes(strings.NewReader("chr1\t-1\n")); err == nil {
		t.Errorf("expected error for negative size")
	}
}

func TestSlopEmptyInput(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	got, warn := runSlop(t, "", sizes, Options{Both: true, BothAdd: 5})
	if got != "" || warn != "" {
		t.Errorf("expected no output and no warnings, got out=%q warn=%q", got, warn)
	}
}

func TestSlopStrandSpecWithBoth(t *testing.T) {
	// With Both=true the symmetric add is applied; swapping sides is a no-op
	// in that case, but the code path must still work.
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\t200\tg1\t0\t-\n"
	got, _ := runSlop(t, input, sizes, Options{Both: true, BothAdd: 10, StrandSpec: true})
	want := "chr1\t90\t210\tg1\t0\t-\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}
