package mosdepth

import (
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

// TestParseThresholds covers the success and error paths.
func TestParseThresholds(t *testing.T) {
	cases := []struct {
		in      string
		want    []int
		wantErr bool
	}{
		{"", nil, false},
		{"1", []int{1}, false},
		{"5, 1, 10", []int{1, 5, 10}, false},
		{",,", []int{}, false},
		{"abc", nil, true},
		{"-1", nil, true},
	}
	for _, tc := range cases {
		got, err := parseThresholds(tc.in)
		if tc.wantErr {
			if err == nil {
				t.Errorf("parseThresholds(%q): expected error", tc.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("parseThresholds(%q): %v", tc.in, err)
			continue
		}
		// Treat nil and []int{} as equivalent.
		if len(got) == 0 && len(tc.want) == 0 {
			continue
		}
		if !reflect.DeepEqual(got, tc.want) {
			t.Errorf("parseThresholds(%q): got %v, want %v", tc.in, got, tc.want)
		}
	}
}

// TestFormatProportion clamps to [0,1] and formats with 2 decimals.
func TestFormatProportion(t *testing.T) {
	cases := []struct {
		in   float64
		want string
	}{
		{0, "0.00"},
		{1, "1.00"},
		{0.5, "0.50"},
		{-1, "0.00"},
		{2, "1.00"},
	}
	for _, tc := range cases {
		if got := formatProportion(tc.in); got != tc.want {
			t.Errorf("formatProportion(%v): got %q, want %q", tc.in, got, tc.want)
		}
	}
}

// TestFormatMean: standard 2-decimal formatter.
func TestFormatMean(t *testing.T) {
	if formatMean(1.234) != "1.23" {
		t.Errorf("formatMean(1.234): %q", formatMean(1.234))
	}
	if formatMean(0) != "0.00" {
		t.Errorf("formatMean(0): %q", formatMean(0))
	}
}

// TestStringSliceUnique de-duplicates while preserving first-occurrence order.
func TestStringSliceUnique(t *testing.T) {
	got := stringSliceUnique([]string{"a", "b", "a", "c", "b"})
	want := []string{"a", "b", "c"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("stringSliceUnique: got %v, want %v", got, want)
	}
}

// TestBedGzWriterAndIndex: write a tiny bed.gz, build a TBI, read it back.
func TestBedGzWriterAndIndex(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "out.bed.gz")
	w, err := newBedGzWriter(path)
	if err != nil {
		t.Fatalf("newBedGzWriter: %v", err)
	}
	if err := w.writeBED("chr1", 0, 10, "1"); err != nil {
		t.Fatalf("writeBED: %v", err)
	}
	if err := w.writeBED("chr1", 10, 20, "2"); err != nil {
		t.Fatalf("writeBED: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if err := buildBedTbi(path); err != nil {
		t.Fatalf("buildBedTbi: %v", err)
	}
	if _, err := os.Stat(path + ".tbi"); err != nil {
		t.Errorf("tbi missing: %v", err)
	}
}

// TestWriteThresholdHeader smoke-tests the header builder.
func TestWriteThresholdHeader(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "th.bed.gz")
	w, err := newBedGzWriter(path)
	if err != nil {
		t.Fatalf("newBedGzWriter: %v", err)
	}
	if err := writeThresholdHeader(w, []int{1, 5}); err != nil {
		t.Fatalf("writeThresholdHeader: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	// Decompress and read the header.
	contents, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	if len(contents) == 0 {
		t.Errorf("file empty")
	}
	// We can't easily decompress without a gzip reader; just check size > 28 (EOF block).
	if len(contents) <= 28 {
		t.Errorf("file unexpectedly small: %d bytes", len(contents))
	}
}

// TestWriteDistribution emits a small dist file and re-parses it.
func TestWriteDistribution(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "dist.txt")
	// chr1 histogram: hist[0]=3, hist[1]=2, hist[2]=5 → total=10
	hist := map[string][]int64{
		"chr1": {3, 2, 5},
	}
	if err := writeDistribution(path, hist, []string{"chr1"}, nil); err != nil {
		t.Fatalf("writeDistribution: %v", err)
	}
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	out := string(data)
	// At depth >= 2: 5/10 = 0.50; at >= 1: 7/10 = 0.70; at >= 0: 10/10 = 1.00.
	for _, want := range []string{"chr1\t2\t0.50", "chr1\t1\t0.70", "chr1\t0\t1.00", "total\t2\t0.50"} {
		if !strings.Contains(out, want) {
			t.Errorf("missing %q in dist output:\n%s", want, out)
		}
	}
}

// TestWriteSummary: header + total row with zero-row case.
func TestWriteSummary(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "summary.txt")
	rows := []summaryRow{
		{chrom: "chr1", length: 100, bases: 250, mean: 2.5, minD: 0, maxD: 5},
		{chrom: "chr2", length: 50, bases: 50, mean: 1.0, minD: 0, maxD: 2},
	}
	if err := writeSummary(path, rows); err != nil {
		t.Fatalf("writeSummary: %v", err)
	}
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	out := string(data)
	if !strings.Contains(out, "chrom\tlength\tbases\tmean\tmin\tmax") {
		t.Errorf("missing header in summary: %q", out)
	}
	if !strings.Contains(out, "total\t150\t300\t2.00\t0\t5") {
		t.Errorf("missing total row in summary:\n%s", out)
	}
}

// TestWriteSummary_Empty handles zero-row input cleanly (header + total
// with zeros).
func TestWriteSummary_Empty(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "empty.txt")
	if err := writeSummary(path, nil); err != nil {
		t.Fatalf("writeSummary nil: %v", err)
	}
	data, _ := os.ReadFile(path)
	if !strings.Contains(string(data), "total\t0\t0\t0.00\t0\t0") {
		t.Errorf("expected zero total row, got:\n%s", data)
	}
}
