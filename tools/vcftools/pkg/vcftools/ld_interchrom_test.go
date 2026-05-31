package vcftools

import (
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// TestInterchromGenoR2_BasicOutput: two SNPs on chr1 and two on chr2.
// --interchrom-geno-r2 should emit only cross-chromosome pairs (4 pairs of
// the 6 total). All four samples are 0/0|0/1|0/1|1/1 at every site, so the
// allele count vectors are identical: r²=1 on every pair.
func TestInterchromGenoR2_BasicOutput(t *testing.T) {
	vcfText := minimalVCF([]string{"s1", "s2", "s3", "s4"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"1", 200, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"2", 100, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"2", 200, []string{"0/0", "0/1", "0/1", "1/1"}},
		})
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, InterchromGenoR2: true}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, err := os.ReadFile(prefix + ".interchrom.geno.ld")
	if err != nil {
		t.Fatalf("read .interchrom.geno.ld: %v", err)
	}
	lines := strings.Split(strings.TrimRight(string(data), "\n"), "\n")
	if lines[0] != "CHR1\tPOS1\tCHR2\tPOS2\tN_INDV\tR^2" {
		t.Errorf("bad header: %q", lines[0])
	}
	// 4 pairs.
	if len(lines) != 5 {
		t.Fatalf("expected 4 data rows, got %d:\n%s", len(lines)-1, data)
	}
	for _, ln := range lines[1:] {
		fields := strings.Split(ln, "\t")
		if len(fields) != 6 {
			t.Errorf("row %q: want 6 fields, got %d", ln, len(fields))
			continue
		}
		if fields[0] == fields[2] {
			t.Errorf("row %q: chromosomes match (cross-chrom only expected)", ln)
		}
		if fields[5] != "1" {
			t.Errorf("row %q: r²=%s want 1", ln, fields[5])
		}
	}
}

// TestInterchromHapR2_OnlyCrossChrom: same as above but for haplotype LD.
// All GTs phased.
func TestInterchromHapR2_OnlyCrossChrom(t *testing.T) {
	vcfText := minimalVCF([]string{"s1", "s2", "s3", "s4"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0|0", "0|1", "0|1", "1|1"}},
			{"2", 100, []string{"0|0", "0|1", "0|1", "1|1"}},
		})
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, InterchromHapR2: true}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, _ := os.ReadFile(prefix + ".interchrom.hap.ld")
	lines := strings.Split(strings.TrimRight(string(data), "\n"), "\n")
	if lines[0] != "CHR1\tPOS1\tCHR2\tPOS2\tN_CHR\tR^2\tD\tDprime" {
		t.Errorf("bad header: %q", lines[0])
	}
	// 1 pair (1 chr1 site × 1 chr2 site).
	if len(lines) != 2 {
		t.Fatalf("expected 1 data row, got %d:\n%s", len(lines)-1, data)
	}
}

// TestGenoChiSq_PerfectAssociation hand-computes chi^2 for two perfectly
// correlated sites.
//
// 4 samples, sites: A=[0,0,2,2], B=[0,0,2,2]. The 3×3 table is:
//
//	g_B = 0  1  2  rowTot
//
// g_A=0   2  0  0   2
// g_A=1   0  0  0   0
// g_A=2   0  0  2   2
// colTot   2  0  2
//
// rows=2, cols=2 (g_A=1 and g_B=1 are empty), df = (2-1)*(2-1) = 1.
//
// Expected counts in the 2×2 sub-table:
// E(0,0) = 2*2/4 = 1, E(0,2) = 2*2/4 = 1, E(2,0) = 1, E(2,2) = 1.
// chi^2 = (2-1)^2/1 + (0-1)^2/1 + (0-1)^2/1 + (2-1)^2/1 = 4.
// p-value at chi^2=4, df=1 ≈ 0.0455.
func TestGenoChiSq_PerfectAssociation(t *testing.T) {
	a := makeLDSite("1", 100, "0/0", "0/0", "1/1", "1/1")
	b := makeLDSite("1", 200, "0/0", "0/0", "1/1", "1/1")
	n, chi2, df, p, ok := computeGenoChiSq(a, b)
	if !ok {
		t.Fatalf("expected ok=true")
	}
	if n != 4 {
		t.Errorf("n=%d want 4", n)
	}
	if df != 1 {
		t.Errorf("df=%d want 1", df)
	}
	if math.Abs(chi2-4) > 1e-9 {
		t.Errorf("chi2=%g want 4", chi2)
	}
	if math.Abs(p-0.04550026389635842) > 1e-6 {
		t.Errorf("p=%g want ≈0.04550026", p)
	}
}

// TestGenoChiSq_NoAssociation: hand-computed independent sites give chi^2=0.
// Counts at site A and B set so g_A and g_B are unrelated.
func TestGenoChiSq_NoAssociation(t *testing.T) {
	// 4 samples: A=[0,1,0,1], B=[0,0,1,1]. Marginal: A has 2 of each cat;
	// B has 2 of each cat. The full table is:
	//          g_B=0  g_B=1
	//  g_A=0     1      1
	//  g_A=1     1      1
	// Expected = 1.0 in every cell -> chi^2 = 0.
	a := makeLDSite("1", 100, "0/0", "0/1", "0/0", "0/1")
	b := makeLDSite("1", 200, "0/0", "0/0", "0/1", "0/1")
	n, chi2, df, _, ok := computeGenoChiSq(a, b)
	if !ok {
		t.Fatalf("expected ok=true")
	}
	if n != 4 || df != 1 {
		t.Errorf("n=%d df=%d", n, df)
	}
	if math.Abs(chi2) > 1e-9 {
		t.Errorf("chi2=%g want 0", chi2)
	}
}

// TestGenoChiSq_Monomorphic: monomorphic site has only 1 row -> df=0 -> ok=false.
func TestGenoChiSq_Monomorphic(t *testing.T) {
	a := makeLDSite("1", 100, "0/0", "0/0", "0/0", "0/0")
	b := makeLDSite("1", 200, "0/0", "0/1", "1/1", "0/0")
	_, _, _, _, ok := computeGenoChiSq(a, b)
	if ok {
		t.Errorf("expected ok=false for monomorphic site")
	}
}

// TestGenoChiSq_Integration: end-to-end. Upstream's
// output_genotype_chisq only emits same-chromosome pairs; for this
// 2-chrom fixture that yields 1 pair per chromosome (2 total).
func TestGenoChiSq_Integration(t *testing.T) {
	vcfText := minimalVCF([]string{"s1", "s2", "s3", "s4"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0/0", "0/0", "1/1", "1/1"}},
			{"1", 200, []string{"0/0", "0/0", "1/1", "1/1"}},
			{"2", 100, []string{"0/0", "0/0", "1/1", "1/1"}},
			{"2", 200, []string{"0/0", "0/0", "1/1", "1/1"}},
		})
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, GenoChiSq: true}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, err := os.ReadFile(prefix + ".geno.chisq")
	if err != nil {
		t.Fatalf("read .geno.chisq: %v", err)
	}
	lines := strings.Split(strings.TrimRight(string(data), "\n"), "\n")
	if lines[0] != "CHR\tPOS1\tPOS2\tN_INDV\tCHI^2\tDOF\tPVAL" {
		t.Errorf("bad header: %q", lines[0])
	}
	if len(lines) != 3 {
		t.Fatalf("expected 2 data rows, got %d:\n%s", len(lines)-1, data)
	}
	// Every pair should have chi^2=4, df=1 (since all sites are identical).
	for _, ln := range lines[1:] {
		fields := strings.Split(ln, "\t")
		if len(fields) != 7 {
			t.Errorf("row %q wrong field count", ln)
			continue
		}
		chi2, err := strconv.ParseFloat(fields[4], 64)
		if err != nil {
			t.Errorf("parse chi2 from %q: %v", ln, err)
			continue
		}
		if math.Abs(chi2-4) > 1e-9 {
			t.Errorf("chi2=%v want 4 (row %q)", chi2, ln)
		}
		if fields[5] != "1" {
			t.Errorf("df=%v want 1 (row %q)", fields[5], ln)
		}
	}
}

// TestChiSquareSurvival sanity-checks the chi-square survival function
// against well-known values.
func TestChiSquareSurvival(t *testing.T) {
	cases := []struct {
		x    float64
		df   int
		want float64
	}{
		// Chi^2 = 3.841 at df=1 -> p ≈ 0.05 (1.96^2 normal-95%).
		{3.841, 1, 0.05},
		// Chi^2 = 5.991 at df=2 -> p ≈ 0.05.
		{5.991, 2, 0.05},
		// Chi^2 = 0 -> p = 1.
		{0, 1, 1.0},
		// Large x -> p ≈ 0.
		{50, 1, 0.0},
	}
	for _, tc := range cases {
		got := chiSquareSurvival(tc.x, tc.df)
		// Use a 1% absolute tolerance.
		if math.Abs(got-tc.want) > 0.01 {
			t.Errorf("chiSquareSurvival(%g, %d) = %g, want ≈%g",
				tc.x, tc.df, got, tc.want)
		}
	}
}
