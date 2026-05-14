package bedsort

// Parity tests against the upstream bedtools sort test suite.
//
// Cases are mirrored from reference_code/bedtools/test/sort/test-sort.sh.
// Inputs live under tools/bedsort/testdata/parity/<case>.bed and expected
// outputs under <case>.expected.bed. Each case feeds the input through the
// library's Run function with the same options the upstream test used and
// asserts byte-for-byte equality.
//
// Cases for upstream features bedsort does not support yet (notably -header,
// which would preserve the leading `#` line) are wrapped in t.Skip with a
// short rationale.

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func readParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runParity(t *testing.T, inputFile string, opts Options) []byte {
	t.Helper()
	in := readParity(t, inputFile)
	var buf bytes.Buffer
	if err := Run(bytes.NewReader(in), &buf, opts); err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	return buf.Bytes()
}

// faidxOrder loads names.txt-style chromosome order using the same parser the
// CLI uses, so the parity test matches `bedtools sort -faidx names.txt`.
func faidxOrder(t *testing.T, name string) []string {
	t.Helper()
	data := readParity(t, name)
	order, err := ReadFaidx(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("ReadFaidx %s: %v", name, err)
	}
	return order
}

func TestParity_Sort_T01_Default(t *testing.T) {
	got := runParity(t, "a.bed", Options{})
	want := readParity(t, "t01_default.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestParity_Sort_T02_SizeA(t *testing.T) {
	got := runParity(t, "a.bed", Options{Mode: ModeSizeA})
	want := readParity(t, "t02_sizeA.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestParity_Sort_T03_SizeD(t *testing.T) {
	got := runParity(t, "a.bed", Options{Mode: ModeSizeD})
	want := readParity(t, "t03_sizeD.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestParity_Sort_T04_ChrThenSizeA(t *testing.T) {
	got := runParity(t, "a.bed", Options{Mode: ModeChrThenSizeA})
	want := readParity(t, "t04_chrThenSizeA.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestParity_Sort_T05_ChrThenSizeD(t *testing.T) {
	got := runParity(t, "a.bed", Options{Mode: ModeChrThenSizeD})
	want := readParity(t, "t05_chrThenSizeD.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestParity_Sort_T06_ChrThenScoreA(t *testing.T) {
	got := runParity(t, "a.bed", Options{Mode: ModeChrThenScoreA})
	want := readParity(t, "t06_chrThenScoreA.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestParity_Sort_T07_ChrThenScoreD(t *testing.T) {
	got := runParity(t, "a.bed", Options{Mode: ModeChrThenScoreD})
	want := readParity(t, "t07_chrThenScoreD.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestParity_Sort_T08_Faidx(t *testing.T) {
	order := faidxOrder(t, "names.txt")
	got := runParity(t, "a.bed", Options{ChromOrder: order})
	want := readParity(t, "t08_faidx.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestParity_Sort_T09_Header(t *testing.T) {
	t.Skip("upstream -header preserves the leading '#' comment line; bedsort " +
		"strips header/comment lines like the default upstream behaviour. " +
		"Tracked as an unimplemented option.")
}

func TestParity_Sort_T10_ZeroLength(t *testing.T) {
	got := runParity(t, "b.bed", Options{})
	want := readParity(t, "t10_zerolen.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// Sanity check that the parity tests resolve their fixtures relative to the
// package directory; protects against a future refactor that moves the files.
func TestParity_Sort_FixtureLayout(t *testing.T) {
	for _, name := range []string{"a.bed", "b.bed", "names.txt"} {
		if data := readParity(t, name); len(data) == 0 {
			t.Fatalf("fixture %s is empty", name)
		}
	}
	if !strings.Contains(string(readParity(t, "a.bed")), "chr7") {
		t.Fatal("fixture a.bed missing expected content")
	}
}
