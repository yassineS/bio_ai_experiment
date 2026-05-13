package bedcomplement

import (
	"bytes"
	"strings"
	"testing"
)

func runComplement(t *testing.T, input string, sizes ChromSizes, order []string) (string, string) {
	t.Helper()
	var out, warn bytes.Buffer
	if _, err := Complement(strings.NewReader(input), &out, &warn, sizes, order); err != nil {
		t.Fatalf("Complement returned error: %v", err)
	}
	return out.String(), warn.String()
}

func TestComplementBasic(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\t200\nchr1\t500\t600\n"
	got, _ := runComplement(t, input, sizes, []string{"chr1"})
	want := "chr1\t0\t100\nchr1\t200\t500\nchr1\t600\t1000\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestComplementStartAtZero(t *testing.T) {
	// First interval starts at 0 -> no leading gap.
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t0\t300\n"
	got, _ := runComplement(t, input, sizes, []string{"chr1"})
	want := "chr1\t300\t1000\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestComplementEndsAtChromSize(t *testing.T) {
	// Last interval ends at chromSize -> no trailing gap.
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\t1000\n"
	got, _ := runComplement(t, input, sizes, []string{"chr1"})
	want := "chr1\t0\t100\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestComplementFullChromCovered(t *testing.T) {
	// Single interval covers the whole chromosome -> no output.
	sizes := ChromSizes{"chr1": 100}
	input := "chr1\t0\t100\n"
	got, _ := runComplement(t, input, sizes, []string{"chr1"})
	if got != "" {
		t.Errorf("expected empty output, got %q", got)
	}
}

func TestComplementChromOnlyInSizes(t *testing.T) {
	// chr2 has no intervals: emit one full-length record.
	sizes := ChromSizes{"chr1": 1000, "chr2": 500}
	input := "chr1\t100\t200\n"
	got, _ := runComplement(t, input, sizes, []string{"chr1", "chr2"})
	want := "chr1\t0\t100\n" +
		"chr1\t200\t1000\n" +
		"chr2\t0\t500\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestComplementSkipsUnknownChrom(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\t200\nchrZ\t10\t20\nchrZ\t30\t40\n"
	got, warn := runComplement(t, input, sizes, []string{"chr1"})
	want := "chr1\t0\t100\nchr1\t200\t1000\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
	if !strings.Contains(warn, "chrZ") {
		t.Errorf("expected warning about chrZ, got %q", warn)
	}
	// Only one warning per chromosome.
	if strings.Count(warn, "chrZ") != 1 {
		t.Errorf("expected exactly one warning about chrZ, got %q", warn)
	}
}

func TestComplementUnsortedDifferentChromError(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000, "chr2": 1000}
	input := "chr1\t100\t200\nchr2\t100\t200\nchr1\t300\t400\n"
	var out, warn bytes.Buffer
	_, err := Complement(strings.NewReader(input), &out, &warn, sizes, nil)
	if err == nil {
		t.Errorf("expected sort error")
	}
}

func TestComplementUnsortedSameChromError(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t300\t400\nchr1\t100\t200\n"
	var out, warn bytes.Buffer
	_, err := Complement(strings.NewReader(input), &out, &warn, sizes, nil)
	if err == nil {
		t.Errorf("expected sort error")
	}
}

func TestComplementOverlappingMerged(t *testing.T) {
	// Overlapping or touching intervals should be merged before complementing.
	sizes := ChromSizes{"chr1": 1000}
	input := "chr1\t100\t200\nchr1\t150\t300\nchr1\t300\t400\n"
	got, _ := runComplement(t, input, sizes, []string{"chr1"})
	want := "chr1\t0\t100\nchr1\t400\t1000\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestComplementSkipsHeadersAndBlanks(t *testing.T) {
	sizes := ChromSizes{"chr1": 1000}
	input := "# comment\ntrack name=foo\nbrowser hide all\n\nchr1\t100\t200\n"
	got, _ := runComplement(t, input, sizes, []string{"chr1"})
	want := "chr1\t0\t100\nchr1\t200\t1000\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestComplementEmptyInput(t *testing.T) {
	sizes := ChromSizes{"chr1": 100, "chr2": 50}
	got, _ := runComplement(t, "", sizes, []string{"chr1", "chr2"})
	want := "chr1\t0\t100\nchr2\t0\t50\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestComplementInvalidStart(t *testing.T) {
	sizes := ChromSizes{"chr1": 100}
	var out, warn bytes.Buffer
	_, err := Complement(strings.NewReader("chr1\tNOPE\t10\n"), &out, &warn, sizes, nil)
	if err == nil {
		t.Errorf("expected error on bad start")
	}
}

func TestComplementInvalidEnd(t *testing.T) {
	sizes := ChromSizes{"chr1": 100}
	var out, warn bytes.Buffer
	_, err := Complement(strings.NewReader("chr1\t10\tNOPE\n"), &out, &warn, sizes, nil)
	if err == nil {
		t.Errorf("expected error on bad end")
	}
}

func TestComplementTooFewFields(t *testing.T) {
	sizes := ChromSizes{"chr1": 100}
	var out, warn bytes.Buffer
	_, err := Complement(strings.NewReader("chr1\t10\n"), &out, &warn, sizes, nil)
	if err == nil {
		t.Errorf("expected error on too few fields")
	}
}

func TestComplementInvalidEndBeforeStart(t *testing.T) {
	sizes := ChromSizes{"chr1": 100}
	var out, warn bytes.Buffer
	_, err := Complement(strings.NewReader("chr1\t50\t10\n"), &out, &warn, sizes, nil)
	if err == nil {
		t.Errorf("expected error on end<start")
	}
}

func TestComplementIntervalClippedToChromSize(t *testing.T) {
	// Interval extends past chromSize: clip and emit the trailing region as
	// fully covered.
	sizes := ChromSizes{"chr1": 100}
	input := "chr1\t50\t1000\n"
	got, _ := runComplement(t, input, sizes, []string{"chr1"})
	want := "chr1\t0\t50\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestReadChromSizes(t *testing.T) {
	in := "# comment\n\nchr1\t100\nchr2\t200\nchr3 50\n"
	sizes, order, err := ReadChromSizes(strings.NewReader(in))
	if err != nil {
		t.Fatalf("ReadChromSizes: %v", err)
	}
	if sizes["chr1"] != 100 || sizes["chr2"] != 200 || sizes["chr3"] != 50 {
		t.Errorf("unexpected sizes: %v", sizes)
	}
	want := []string{"chr1", "chr2", "chr3"}
	if len(order) != len(want) {
		t.Fatalf("order mismatch: %v want %v", order, want)
	}
	for i, c := range want {
		if order[i] != c {
			t.Errorf("order[%d]=%q want %q", i, order[i], c)
		}
	}
}

func TestReadChromSizesBadFormat(t *testing.T) {
	if _, _, err := ReadChromSizes(strings.NewReader("only\n")); err == nil {
		t.Errorf("expected error")
	}
	if _, _, err := ReadChromSizes(strings.NewReader("chr1\tNOPE\n")); err == nil {
		t.Errorf("expected error")
	}
	if _, _, err := ReadChromSizes(strings.NewReader("chr1\t-5\n")); err == nil {
		t.Errorf("expected error for negative size")
	}
}

func TestComplementOrderingLexicographicFallback(t *testing.T) {
	// chromOrder=nil -> chromosomes sorted lexicographically.
	sizes := ChromSizes{"chrZ": 10, "chr1": 10, "chr2": 10}
	got, _ := runComplement(t, "", sizes, nil)
	want := "chr1\t0\t10\nchr2\t0\t10\nchrZ\t0\t10\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestComplementOrderingPrefersUserOrder(t *testing.T) {
	sizes := ChromSizes{"chrZ": 10, "chr1": 10, "chr2": 10}
	got, _ := runComplement(t, "", sizes, []string{"chrZ", "chr2"})
	// chrZ then chr2 from explicit order, then chr1 (the unlisted one) at
	// the end in lex order.
	want := "chrZ\t0\t10\nchr2\t0\t10\nchr1\t0\t10\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestComplementOrderingDedupsAndIgnoresUnknown(t *testing.T) {
	sizes := ChromSizes{"chr1": 10, "chr2": 10}
	got, _ := runComplement(t, "", sizes, []string{"chr1", "chr1", "chrUnknown", "chr2"})
	want := "chr1\t0\t10\nchr2\t0\t10\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}
