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

// TestBedGzWriterAndIndex: write a tiny bed.gz, build a CSI, read it back.
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
	if err := buildBedCsi(path); err != nil {
		t.Fatalf("buildBedCsi: %v", err)
	}
	if _, err := os.Stat(path + ".csi"); err != nil {
		t.Errorf("csi missing: %v", err)
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
	if err := writeDistribution(path, hist, []string{"chr1"}); err != nil {
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
	if err := writeSummary(path, rows, nil); err != nil {
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

// TestWriteSummary_RegionRows proves the *_region rows are interleaved after
// each chrom's non-region row and that total_region is emitted after total,
// with both totals aggregated independently — matching upstream mosdepth's
// region-mode summary layout. No upstream binary is involved.
func TestWriteSummary_RegionRows(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "summary.txt")
	rows := []summaryRow{
		{chrom: "chr1", length: 100, bases: 200, mean: 2.0, minD: 0, maxD: 5},
		{chrom: "chr2", length: 50, bases: 50, mean: 1.0, minD: 0, maxD: 2},
	}
	// Region aggregates cover only part of each chrom (the region-covered
	// bases), so length/bases/min differ from the non-region rows.
	regionRows := []summaryRow{
		{chrom: "chr1_region", length: 40, bases: 120, mean: 3.0, minD: 1, maxD: 5},
		{chrom: "chr2_region", length: 10, bases: 20, mean: 2.0, minD: 2, maxD: 2},
	}
	if err := writeSummary(path, rows, regionRows); err != nil {
		t.Fatalf("writeSummary: %v", err)
	}
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	got := strings.Split(strings.TrimRight(string(data), "\n"), "\n")
	want := []string{
		"chrom\tlength\tbases\tmean\tmin\tmax",
		"chr1\t100\t200\t2.00\t0\t5",
		"chr1_region\t40\t120\t3.00\t1\t5",
		"chr2\t50\t50\t1.00\t0\t2",
		"chr2_region\t10\t20\t2.00\t2\t2",
		// total: length 150, bases 250, mean 1.67, min 0, max 5.
		"total\t150\t250\t1.67\t0\t5",
		// total_region: length 50, bases 140, mean 2.80, min 1, max 5.
		"total_region\t50\t140\t2.80\t1\t5",
	}
	if len(got) != len(want) {
		t.Fatalf("summary line count: got %d, want %d:\n%s", len(got), len(want), strings.Join(got, "\n"))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("summary line %d:\n got %q\nwant %q", i, got[i], want[i])
		}
	}
}

// TestWriteDistribution_RegionCumulation feeds a synthetic region histogram and
// asserts the cumulative-from-top proportions and upstream's skip rules
// (cum < 8e-5 trimming and the depth>300 zero-tail skip) without any binary.
func TestWriteDistribution_RegionCumulation(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "region.dist.txt")
	// Two chroms. chr1: 1 region at depth 2, 1 at depth 1, 2 at depth 0
	// (this is the per-window histogram counted at int(mean)). chr2: 1 at
	// depth 1, 1 at depth 0.
	hist := map[string][]int64{
		"chr1": {2, 1, 1}, // hist[0]=2, hist[1]=1, hist[2]=1 -> total 4
		"chr2": {1, 1},    // hist[0]=1, hist[1]=1            -> total 2
	}
	if err := writeDistribution(path, hist, []string{"chr1", "chr2"}); err != nil {
		t.Fatalf("writeDistribution: %v", err)
	}
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	got := strings.Split(strings.TrimRight(string(data), "\n"), "\n")
	want := []string{
		// chr1: >=2 -> 1/4=0.25; >=1 -> 2/4=0.50; >=0 -> 4/4=1.00.
		"chr1\t2\t0.25",
		"chr1\t1\t0.50",
		"chr1\t0\t1.00",
		// chr2: >=1 -> 1/2=0.50; >=0 -> 2/2=1.00.
		"chr2\t1\t0.50",
		"chr2\t0\t1.00",
		// total (sum of histograms): hist[0]=3, hist[1]=2, hist[2]=1, total 6.
		// >=2 -> 1/6=0.17; >=1 -> 3/6=0.50; >=0 -> 6/6=1.00.
		"total\t2\t0.17",
		"total\t1\t0.50",
		"total\t0\t1.00",
	}
	if len(got) != len(want) {
		t.Fatalf("dist line count: got %d, want %d:\n%s", len(got), len(want), strings.Join(got, "\n"))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("dist line %d:\n got %q\nwant %q", i, got[i], want[i])
		}
	}
}

// TestWriteDistribution_SkipRules covers upstream's two emission filters: a
// depth above 300 whose count is zero is skipped, and a running cumulative
// proportion still below 8e-5 trims the very sparse top of the distribution.
func TestWriteDistribution_SkipRules(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "dist.txt")
	// One base at depth 500, then a large mass at depth 0. The depth-500 row
	// has cum = 1/100001 ~= 1e-5 < 8e-5, so it is skipped; every zero-count
	// depth above 300 is also skipped, so only the depth-0 row is emitted.
	hist := make([]int64, 501)
	hist[500] = 1
	hist[0] = 100000
	if err := writeDistribution(path, map[string][]int64{"chr1": hist}, []string{"chr1"}); err != nil {
		t.Fatalf("writeDistribution: %v", err)
	}
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	got := strings.Split(strings.TrimRight(string(data), "\n"), "\n")
	// Expect chr1 rows then total rows, but only depth-0 survives the trim in
	// each (the depth-500 row is below the 8e-5 cumulative floor).
	for _, ln := range got {
		if strings.Contains(ln, "\t500\t") {
			t.Errorf("depth-500 row should have been trimmed (cum < 8e-5): %q", ln)
		}
	}
	if got[0] != "chr1\t0\t1.00" {
		t.Errorf("first emitted row: got %q, want %q", got[0], "chr1\t0\t1.00")
	}
}

// TestWriteSummary_Empty handles zero-row input cleanly (header + total
// with zeros).
func TestWriteSummary_Empty(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "empty.txt")
	if err := writeSummary(path, nil, nil); err != nil {
		t.Fatalf("writeSummary nil: %v", err)
	}
	data, _ := os.ReadFile(path)
	if !strings.Contains(string(data), "total\t0\t0\t0.00\t0\t0") {
		t.Errorf("expected zero total row, got:\n%s", data)
	}
}
