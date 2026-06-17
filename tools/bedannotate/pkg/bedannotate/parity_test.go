package bedannotate

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// Fixture-based parity tests for bedannotate. The expected outputs are
// generated directly from the upstream `bedtools annotate` binary (bedtools
// 2.31.1), so these assert byte-for-byte equality against real upstream output:
//
//   * no header is emitted unless -names is given;
//   * the header pads the leading '#' with bedType-1 tabs;
//   * records are reported grouped by chromosome then UCSC bin;
//   * `-counts` / `-both` column semantics match upstream.

func readFile(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func parityRun(t *testing.T, aFile string, bFiles []string, opts Options) []byte {
	t.Helper()
	a := readFile(t, aFile)
	bRs := make([]io.Reader, len(bFiles))
	for i, name := range bFiles {
		bRs[i] = bytes.NewReader(readFile(t, name))
	}
	var out bytes.Buffer
	if _, err := Run(bytes.NewReader(a), bRs, &out, opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	return out.Bytes()
}

// case1.default — three A intervals across two B files, default fractions, no
// header (no -names).
func TestParity_Annotate_Case1_Default(t *testing.T) {
	got := parityRun(t, "case1.a.bed",
		[]string{"case1.b1.bed", "case1.b2.bed"}, Options{})
	want := readFile(t, "case1.expected.default.txt")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\n got:\n%s", want, got)
	}
}

// case1.counts — -counts swaps fractions for record counts (still no header).
func TestParity_Annotate_Case1_Counts(t *testing.T) {
	got := parityRun(t, "case1.a.bed",
		[]string{"case1.b1.bed", "case1.b2.bed"}, Options{Mode: ModeCounts})
	want := readFile(t, "case1.expected.counts.txt")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\n got:\n%s", want, got)
	}
}

// case1.both — -both interleaves count + fraction; -names emits a header padded
// with bedType-1 tabs.
func TestParity_Annotate_Case1_Both(t *testing.T) {
	got := parityRun(t, "case1.a.bed",
		[]string{"case1.b1.bed", "case1.b2.bed"},
		Options{Mode: ModeBoth, Names: []string{"b1", "b2"}})
	want := readFile(t, "case1.expected.both.txt")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\n got:\n%s", want, got)
	}
}
