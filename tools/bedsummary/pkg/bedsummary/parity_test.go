package bedsummary

// Parity tests for `bedtools summary`. The upstream bedtools test corpus
// (reference_code/bedtools/test/) does NOT ship a `summary/` subdirectory,
// so these are spec-driven cases derived from the bedtools v2.31 manual:
//
//   https://bedtools.readthedocs.io/en/latest/content/tools/summary.html
//
// Output columns: chrom, num_ivls, total_ivl_bp, min_ivl_bp, max_ivl_bp,
// mean_ivl_bp, median_ivl_bp, plus a trailing `all` aggregate row
// (suppressible with --skip-all). The column NAMES we emit follow the
// upstream tool's "interval/bp" naming used in newer bedtools releases.
// Numeric values use the project-wide 3-digit precision rule via
// formatNum.

import (
	"bytes"
	"strings"
	"testing"
)

// TestParity_Basic_TwoChroms exercises the canonical 3-record / 2-chrom
// example from the upstream manual.
func TestParity_Basic_TwoChroms(t *testing.T) {
	in := strings.NewReader("chr1\t0\t10\nchr1\t100\t200\nchr2\t0\t50\n")
	var out bytes.Buffer
	if err := Run(in, &out, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chrom\tnum_ivls\ttotal_ivl_bp\tmin_ivl_bp\tmax_ivl_bp\tmean_ivl_bp\tmedian_ivl_bp\n" +
		"chr1\t2\t110\t10\t100\t55\t55\n" +
		"chr2\t1\t50\t50\t50\t50\t50\n" +
		"all\t3\t160\t10\t100\t53.333\t50\n"
	if got := out.String(); got != want {
		t.Errorf("summary mismatch.\ngot:\n%s\nwant:\n%s", got, want)
	}
}

// TestParity_NoHeader confirms --no-header drops the column-header line
// only — the per-chrom rows and the trailing aggregate are unchanged.
func TestParity_NoHeader(t *testing.T) {
	in := strings.NewReader("chr1\t0\t10\n")
	var out bytes.Buffer
	if err := Run(in, &out, Options{NoHeader: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	got := out.String()
	if strings.HasPrefix(got, "chrom\t") {
		t.Errorf("expected no header line, got prefix %q", got[:50])
	}
	if !strings.Contains(got, "chr1\t1\t10\t10\t10\t10\t10\n") {
		t.Errorf("missing chr1 row:\n%s", got)
	}
	if !strings.Contains(got, "all\t1\t10\t10\t10\t10\t10\n") {
		t.Errorf("missing all aggregate row:\n%s", got)
	}
}

// TestParity_SkipAll confirms --skip-all suppresses the trailing
// aggregate row but keeps the column header and per-chrom rows.
func TestParity_SkipAll(t *testing.T) {
	in := strings.NewReader("chr1\t0\t10\nchr2\t0\t30\n")
	var out bytes.Buffer
	if err := Run(in, &out, Options{SkipAll: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	got := out.String()
	if strings.Contains(got, "\nall\t") {
		t.Errorf("expected no all row, got:\n%s", got)
	}
	if !strings.Contains(got, "chr1\t1\t10") || !strings.Contains(got, "chr2\t1\t30") {
		t.Errorf("missing per-chrom rows:\n%s", got)
	}
}

// TestParity_OddCountMedian sanity-checks that an odd-count chromosome
// uses the middle element (not the average of two middles), matching
// upstream's `median()` helper in `summary/SummaryFile.cpp`.
func TestParity_OddCountMedian(t *testing.T) {
	in := strings.NewReader("chr1\t0\t1\nchr1\t0\t5\nchr1\t0\t10\n")
	var out bytes.Buffer
	if err := Run(in, &out, Options{SkipAll: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	// Sorted lengths: 1,5,10 → median 5.
	got := out.String()
	// Sorted lengths: 1,5,10 → min=1 max=10 mean=5.333 median=5.
	if !strings.Contains(got, "chr1\t3\t16\t1\t10\t5.333\t5\n") {
		t.Errorf("expected min=1 max=10 median=5 row, got:\n%s", got)
	}
}
