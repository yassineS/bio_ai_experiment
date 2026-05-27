package bedmultiinter

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// Upstream ships no `multiinter/` test subdir under
// `reference_code/bedtools/test/`. The expected outputs here are the
// upstream's own example from `multiIntersectBed/multiIntersectBedMain.cpp`
// (the `multiintersect_examples()` heredoc); the fixtures are byte-for-byte
// copies of the input files documented there.

func parityFixture(t *testing.T, name string) []byte {
	t.Helper()
	p := filepath.Join("..", "..", "testdata", "parity", name)
	data, err := os.ReadFile(p)
	if err != nil {
		t.Fatalf("read fixture %s: %v", name, err)
	}
	return data
}

func parityReader(t *testing.T, name string) io.Reader {
	t.Helper()
	return bytes.NewReader(parityFixture(t, name))
}

// Parity 1: the canonical default-mode example. Three BED files, eight
// segments emitted.
func TestParity_DefaultExample(t *testing.T) {
	want := parityFixture(t, "example.default.expected")
	var got bytes.Buffer
	if _, err := Run([]io.Reader{
		parityReader(t, "a.bed"),
		parityReader(t, "b.bed"),
		parityReader(t, "c.bed"),
	}, &got, Options{Filenames: []string{"a.bed", "b.bed", "c.bed"}}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity 2: -header + -names. Header row, then the same 8 segments using
// the supplied labels in the 'list' column.
func TestParity_HeaderAndNames(t *testing.T) {
	want := parityFixture(t, "example.header_names.expected")
	var got bytes.Buffer
	if _, err := Run([]io.Reader{
		parityReader(t, "a.bed"),
		parityReader(t, "b.bed"),
		parityReader(t, "c.bed"),
	}, &got, Options{
		Filenames: []string{"a.bed", "b.bed", "c.bed"},
		Names:     []string{"A", "B", "C"},
		Header:    true,
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity 3: -header + -names + -empty + -g. Two extra 0-count gap rows
// frame the eight data rows (chr1:0..6 and chr1:34..5000).
func TestParity_EmptyWithGenome(t *testing.T) {
	want := parityFixture(t, "example.empty.expected")
	var got bytes.Buffer
	sizes := map[string]int{"chr1": 5000}
	if _, err := Run([]io.Reader{
		parityReader(t, "a.bed"),
		parityReader(t, "b.bed"),
		parityReader(t, "c.bed"),
	}, &got, Options{
		Filenames:  []string{"a.bed", "b.bed", "c.bed"},
		Names:      []string{"A", "B", "C"},
		Header:     true,
		Empty:      true,
		ChromSizes: sizes,
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// TestParity_DocumentedFixtureProvenance is an executable footnote. The
// upstream bedtools repo does not ship a test/multiinter/ subdir, so the
// parity fixtures used above are derived from the multiintersect_examples()
// heredoc in src/multiIntersectBed/multiIntersectBedMain.cpp. The test
// logs that provenance so `go test -v` makes the rationale discoverable
// without polluting the skip count.
func TestParity_DocumentedFixtureProvenance(t *testing.T) {
	t.Log("upstream ships no test/multiinter/ subdir; parity fixtures above " +
		"come from multiintersect_examples() in " +
		"src/multiIntersectBed/multiIntersectBedMain.cpp.")
}
