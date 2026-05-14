package vcftools

import (
	"bytes"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// approxEqual compares two floats with a tolerance suitable for the LD math
// (which involves divisions but no transcendental functions).
func approxEqual(a, b, tol float64) bool {
	if math.IsNaN(a) || math.IsNaN(b) {
		return false
	}
	return math.Abs(a-b) <= tol
}

// makeLDSite builds a synthetic site from per-sample GT strings, parsing them
// the same way the production code does.
func makeLDSite(chrom string, pos int, gts ...string) *ldSite {
	v := &vcf.Variant{
		Chrom:  chrom,
		Pos:    pos,
		Ref:    "A",
		Alt:    []string{"G"},
		Format: []string{"GT"},
	}
	for i, gt := range gts {
		v.Samples = append(v.Samples, vcf.Sample{
			Name: "s" + string(rune('1'+i)),
			Data: map[string]string{"GT": gt},
		})
	}
	site, ok := extractLDSite(v)
	if !ok {
		panic("extractLDSite failed")
	}
	return site
}

// TestParseGTForLD covers the GT parser used by both --geno-r2 and --hap-r2.
// Hand-computed expectations: only alleles {0,1} count; phased GTs ("a|b")
// expose haplotypes; unphased GTs ("a/b") expose only the diploid count.
func TestParseGTForLD(t *testing.T) {
	tests := []struct {
		name                         string
		gt                           string
		wantGeno, wantHapA, wantHapB int
	}{
		{"missing dot", ".", -1, -1, -1},
		{"missing diploid", "./.", -1, -1, -1},
		{"missing phased", ".|.", -1, -1, -1},
		{"hom ref unphased", "0/0", 0, -1, -1},
		{"het unphased", "0/1", 1, -1, -1},
		{"hom alt unphased", "1/1", 2, -1, -1},
		{"hom ref phased", "0|0", 0, 0, 0},
		{"het phased 01", "0|1", 1, 0, 1},
		{"het phased 10", "1|0", 1, 1, 0},
		{"hom alt phased", "1|1", 2, 1, 1},
		{"multi-allelic 0/2 skipped", "0/2", -1, -1, -1},
		{"multi-allelic 1|2 skipped", "1|2", -1, -1, -1},
		{"haploid skipped", "0", -1, -1, -1},
		{"garbage skipped", "x/y", -1, -1, -1},
		{"half-missing skipped", "0/.", -1, -1, -1},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			gc, ha, hb := parseGTForLD(tc.gt)
			if gc != tc.wantGeno || ha != tc.wantHapA || hb != tc.wantHapB {
				t.Errorf("parseGTForLD(%q) = (%d,%d,%d), want (%d,%d,%d)",
					tc.gt, gc, ha, hb, tc.wantGeno, tc.wantHapA, tc.wantHapB)
			}
		})
	}
}

// TestComputeGenoR2 verifies r² for a hand-computed example.
//
// Four samples; site A counts g = [0, 1, 1, 2]; site B counts h = [0, 1, 1, 2].
//
//	mean(g) = mean(h) = 1
//	var(g)  = var(h)  = ((1) + 0 + 0 + (1))/4 = 0.5
//	cov(g,h) = ((-1)(-1) + 0 + 0 + (1)(1))/4 = 0.5
//	r² = 0.5² / (0.5 * 0.5) = 1.0
func TestComputeGenoR2_PerfectCorrelation(t *testing.T) {
	a := makeLDSite("1", 100, "0/0", "0/1", "0/1", "1/1")
	b := makeLDSite("1", 200, "0/0", "0/1", "0/1", "1/1")
	n, r2, ok := computeGenoR2(a, b)
	if !ok {
		t.Fatalf("expected ok=true, got false")
	}
	if n != 4 {
		t.Errorf("n=%d want 4", n)
	}
	if !approxEqual(r2, 1.0, 1e-12) {
		t.Errorf("r2=%v want 1.0", r2)
	}
}

// TestComputeGenoR2_NoCorrelation uses orthogonal vectors:
//
//	g = [0,0,2,2], h = [0,2,0,2] -> cov(g,h) = 0 -> r² = 0.
func TestComputeGenoR2_NoCorrelation(t *testing.T) {
	a := makeLDSite("1", 100, "0/0", "0/0", "1/1", "1/1")
	b := makeLDSite("1", 200, "0/0", "1/1", "0/0", "1/1")
	n, r2, ok := computeGenoR2(a, b)
	if !ok {
		t.Fatalf("expected ok=true, got false")
	}
	if n != 4 {
		t.Errorf("n=%d want 4", n)
	}
	if !approxEqual(r2, 0.0, 1e-12) {
		t.Errorf("r2=%v want 0", r2)
	}
}

// TestComputeGenoR2_MonomorphicSkipped: var(g)=0 must return ok=false.
func TestComputeGenoR2_MonomorphicSkipped(t *testing.T) {
	a := makeLDSite("1", 100, "0/0", "0/0", "0/0", "0/0")
	b := makeLDSite("1", 200, "0/0", "0/1", "1/1", "0/1")
	_, _, ok := computeGenoR2(a, b)
	if ok {
		t.Errorf("expected ok=false for monomorphic site")
	}
}

// TestComputeGenoR2_NLessThan2: when fewer than 2 samples are non-missing at
// both sites, skip.
func TestComputeGenoR2_NLessThan2(t *testing.T) {
	a := makeLDSite("1", 100, "0/0", "./.", "./.")
	b := makeLDSite("1", 200, "1/1", "0/0", "0/1")
	n, _, ok := computeGenoR2(a, b)
	if ok {
		t.Errorf("expected ok=false when n<2, got ok=true (n=%d)", n)
	}
	if n != 1 {
		t.Errorf("n=%d want 1", n)
	}
}

// TestComputeGenoR2_MissingPairwise: a sample missing at one site only is
// dropped from the calculation entirely, not counted as zero.
//
// site A = [0/0, 0/1, 1/1, ./.], site B = [0/0, 0/1, 1/1, 0/1]
// Sample 4 is missing at A so it's excluded. Restricted sample set:
//
//	g = [0,1,2], h = [0,1,2] -> mean=1, var=2/3 each, cov=2/3 -> r²=1.
func TestComputeGenoR2_MissingPairwise(t *testing.T) {
	a := makeLDSite("1", 100, "0/0", "0/1", "1/1", "./.")
	b := makeLDSite("1", 200, "0/0", "0/1", "1/1", "0/1")
	n, r2, ok := computeGenoR2(a, b)
	if !ok {
		t.Fatalf("expected ok=true")
	}
	if n != 3 {
		t.Errorf("n=%d want 3", n)
	}
	if !approxEqual(r2, 1.0, 1e-12) {
		t.Errorf("r2=%v want 1.0", r2)
	}
}

// TestComputeHapR2_PerfectLD: two phased sites with all "0|0" and "1|1"
// haplotypes -> pA=pB=0.5, pAB=0.5, D=0.25, Dmax=0.25 (since D>=0,
// min(pA*(1-pB), (1-pA)*pB) = min(0.25, 0.25) = 0.25), Dprime=1,
// r² = 0.25²/(0.5*0.5*0.5*0.5) = 0.0625/0.0625 = 1.
func TestComputeHapR2_PerfectLD(t *testing.T) {
	a := makeLDSite("1", 100, "0|0", "0|0", "1|1", "1|1")
	b := makeLDSite("1", 200, "0|0", "0|0", "1|1", "1|1")
	n, r2, D, Dp, ok := computeHapR2(a, b)
	if !ok {
		t.Fatalf("expected ok=true")
	}
	if n != 8 {
		t.Errorf("nChr=%d want 8", n)
	}
	if !approxEqual(r2, 1.0, 1e-12) {
		t.Errorf("r2=%v want 1.0", r2)
	}
	if !approxEqual(D, 0.25, 1e-12) {
		t.Errorf("D=%v want 0.25", D)
	}
	if !approxEqual(Dp, 1.0, 1e-12) {
		t.Errorf("Dprime=%v want 1.0", Dp)
	}
}

// TestComputeHapR2_HandComputed: a known 4-haplotype example.
// Phased GTs (two individuals): "0|1" and "1|0". Site A haplotypes: 0,1,1,0;
// site B haplotypes: 1,0,0,1. nChr=4, pA=0.5, pB=0.5, pAB=freq("0|1") at
// (A,B). Pair (A=0,B=1): individual1.hapA=0,b=1 -> AB=01 (no). Let's enumerate:
//
//	ind1 hapA: A=0, B=1
//	ind1 hapB: A=1, B=0
//	ind2 hapA: A=1, B=0
//	ind2 hapB: A=0, B=1
//
// pAB = freq(A=0,B=0) = 0/4 = 0. D = 0 - 0.5*0.5 = -0.25.
// |D|=0.25; D<0 so Dmax = min(0.5*0.5, 0.5*0.5) = 0.25; Dprime = -1.
// r² = 0.0625 / (0.5*0.5*0.5*0.5) = 1.
func TestComputeHapR2_PerfectNegativeLD(t *testing.T) {
	a := makeLDSite("1", 100, "0|1", "1|0")
	b := makeLDSite("1", 200, "1|0", "0|1")
	n, r2, D, Dp, ok := computeHapR2(a, b)
	if !ok {
		t.Fatalf("expected ok=true")
	}
	if n != 4 {
		t.Errorf("nChr=%d want 4", n)
	}
	if !approxEqual(r2, 1.0, 1e-12) {
		t.Errorf("r2=%v want 1.0", r2)
	}
	if !approxEqual(D, -0.25, 1e-12) {
		t.Errorf("D=%v want -0.25", D)
	}
	if !approxEqual(Dp, -1.0, 1e-12) {
		t.Errorf("Dprime=%v want -1.0", Dp)
	}
}

// TestComputeHapR2_IgnoresUnphased: an unphased GT should not contribute to
// --hap-r2 nChr.
func TestComputeHapR2_IgnoresUnphased(t *testing.T) {
	a := makeLDSite("1", 100, "0|0", "0/1", "1|1")
	b := makeLDSite("1", 200, "0|0", "0/1", "1|1")
	n, _, _, _, _ := computeHapR2(a, b)
	if n != 4 {
		t.Errorf("nChr=%d want 4 (unphased middle sample skipped)", n)
	}
}

// TestComputeHapR2_MonomorphicSkipped: pA == 0 or 1 (or pB) -> ok=false.
func TestComputeHapR2_MonomorphicSkipped(t *testing.T) {
	a := makeLDSite("1", 100, "0|0", "0|0", "0|0")
	b := makeLDSite("1", 200, "0|1", "1|0", "0|0")
	_, _, _, _, ok := computeHapR2(a, b)
	if ok {
		t.Errorf("expected ok=false when one site is monomorphic")
	}
}

// TestComputeHapR2_NLessThan2: only one phased non-missing individual -> nChr=2 still
// fine, but with zero phased individuals nChr<2 -> ok=false.
func TestComputeHapR2_NLessThan2(t *testing.T) {
	a := makeLDSite("1", 100, "0/1", "0/1") // both unphased
	b := makeLDSite("1", 200, "0|1", "1|0")
	_, _, _, _, ok := computeHapR2(a, b)
	if ok {
		t.Errorf("expected ok=false when no phased haplotypes overlap")
	}
}

// TestWithinLDWindow exercises the window-distance logic.
func TestWithinLDWindow(t *testing.T) {
	a := &ldSite{chrom: "1", pos: 100, chromIdx: 1}
	b := &ldSite{chrom: "1", pos: 1100, chromIdx: 5}
	cross := &ldSite{chrom: "2", pos: 200, chromIdx: 2}

	// Different chromosome -> always false.
	if withinLDWindow(a, cross, &Params{}) {
		t.Error("expected false for different chroms")
	}

	// Unbounded -> true.
	if !withinLDWindow(a, b, &Params{}) {
		t.Error("expected true for unbounded window")
	}

	// SNP window too tight.
	if withinLDWindow(a, b, &Params{LDWindow: 2}) {
		t.Error("expected false when SNP dist (4) > ld-window (2)")
	}

	// bp window too tight.
	if withinLDWindow(a, b, &Params{LDWindowBp: 500}) {
		t.Error("expected false when bp dist (1000) > ld-window-bp (500)")
	}

	// Minimum SNP distance not met.
	if withinLDWindow(a, b, &Params{LDWindowMin: 10}) {
		t.Error("expected false when SNP dist (4) < ld-window-min (10)")
	}

	// Minimum bp distance not met.
	if withinLDWindow(a, b, &Params{LDWindowBpMin: 2000}) {
		t.Error("expected false when bp dist (1000) < ld-window-bp-min (2000)")
	}

	// All constraints satisfied.
	if !withinLDWindow(a, b, &Params{LDWindow: 10, LDWindowBp: 2000, LDWindowMin: 1, LDWindowBpMin: 100}) {
		t.Error("expected true within all bounds")
	}
}

// TestLDPositionAllowed covers --geno-r2-positions / --hap-r2-positions
// gating.
func TestLDPositionAllowed(t *testing.T) {
	pos := positionSet{"1": {100: true}}
	a := &ldSite{chrom: "1", pos: 100}
	b := &ldSite{chrom: "1", pos: 200}
	c := &ldSite{chrom: "1", pos: 300}

	if !ldPositionAllowed(a, b, pos, true) {
		t.Error("expected true when a is in positions set")
	}
	if !ldPositionAllowed(b, a, pos, true) {
		t.Error("expected true when b is in positions set")
	}
	if ldPositionAllowed(b, c, pos, true) {
		t.Error("expected false when neither endpoint is in set")
	}
	if !ldPositionAllowed(b, c, nil, false) {
		t.Error("expected true when restricted is false")
	}
}

// runLDIntegration runs Run() on an in-memory VCF and returns the LD output
// files as strings (whichever exist). Failure to read a file is fatal.
func runLDIntegration(t *testing.T, vcfText string, params *Params) (geno, hap string) {
	t.Helper()
	dir := t.TempDir()
	params.OutPrefix = filepath.Join(dir, "out")
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run() failed: %v", err)
	}
	if data, err := os.ReadFile(params.OutPrefix + ".geno.ld"); err == nil {
		geno = string(data)
	}
	if data, err := os.ReadFile(params.OutPrefix + ".hap.ld"); err == nil {
		hap = string(data)
	}
	return geno, hap
}

// minimalVCF builds a small VCF with N samples and the given (chrom,pos,GTs)
// rows. GTs is one string per sample.
func minimalVCF(samples []string, rows []struct {
	chrom string
	pos   int
	gts   []string
}) string {
	var b bytes.Buffer
	b.WriteString("##fileformat=VCFv4.2\n")
	b.WriteString("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
	for _, s := range samples {
		b.WriteString("\t")
		b.WriteString(s)
	}
	b.WriteString("\n")
	for _, r := range rows {
		b.WriteString(r.chrom)
		b.WriteString("\t")
		// Pos as string via Sprintf-free path
		b.WriteString(itoa(r.pos))
		b.WriteString("\t.\tA\tG\t.\tPASS\t.\tGT")
		for _, gt := range r.gts {
			b.WriteString("\t")
			b.WriteString(gt)
		}
		b.WriteString("\n")
	}
	return b.String()
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	neg := false
	if n < 0 {
		neg = true
		n = -n
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
}

// TestRun_GenoR2_BasicOutput: three sites in perfect LD with each other
// produce three pairs (1-2, 1-3, 2-3) each with r²=1.
func TestRun_GenoR2_BasicOutput(t *testing.T) {
	vcfText := minimalVCF([]string{"s1", "s2", "s3", "s4"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"1", 200, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"1", 300, []string{"0/0", "0/1", "0/1", "1/1"}},
		})
	params := &Params{GenoR2: true, MaxMissing: 1}
	geno, hap := runLDIntegration(t, vcfText, params)
	if hap != "" {
		t.Errorf("expected no .hap.ld output, got %q", hap)
	}
	lines := strings.Split(strings.TrimSpace(geno), "\n")
	if len(lines) != 4 { // header + 3 pairs
		t.Fatalf("expected 4 lines (header + 3 pairs), got %d:\n%s", len(lines), geno)
	}
	if !strings.HasPrefix(lines[0], "CHR\tPOS1\tPOS2\tN_INDV\tR^2") {
		t.Errorf("missing/incorrect header: %q", lines[0])
	}
	for i, line := range lines[1:] {
		fields := strings.Split(line, "\t")
		if len(fields) != 5 {
			t.Errorf("row %d: want 5 fields, got %d (%q)", i, len(fields), line)
			continue
		}
		if fields[3] != "4" {
			t.Errorf("row %d: N_INDV=%s want 4", i, fields[3])
		}
		if fields[4] != "1" {
			t.Errorf("row %d: R^2=%s want 1", i, fields[4])
		}
	}
}

// TestRun_LDWindow restricts pairs by SNP and bp window. With three sites at
// 100, 200, 1500 and ld-window-bp=500, only the (100,200) pair survives.
func TestRun_LDWindow(t *testing.T) {
	vcfText := minimalVCF([]string{"s1", "s2", "s3", "s4"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"1", 200, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"1", 1500, []string{"0/0", "0/1", "0/1", "1/1"}},
		})
	params := &Params{GenoR2: true, LDWindowBp: 500, MaxMissing: 1}
	geno, _ := runLDIntegration(t, vcfText, params)
	lines := strings.Split(strings.TrimSpace(geno), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected header + 1 pair, got %d lines:\n%s", len(lines), geno)
	}
	fields := strings.Split(lines[1], "\t")
	if fields[1] != "100" || fields[2] != "200" {
		t.Errorf("expected pair (100,200), got (%s,%s)", fields[1], fields[2])
	}
}

// TestRun_MinR2_Filter: with --min-r2=0.5 the uncorrelated pair drops out.
func TestRun_MinR2_Filter(t *testing.T) {
	vcfText := minimalVCF([]string{"s1", "s2", "s3", "s4"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0/0", "0/0", "1/1", "1/1"}},
			{"1", 200, []string{"0/0", "1/1", "0/0", "1/1"}}, // orthogonal -> r²=0
			{"1", 300, []string{"0/0", "0/0", "1/1", "1/1"}}, // perfect with site 1
		})
	params := &Params{GenoR2: true, MinR2: 0.5, MaxMissing: 1}
	geno, _ := runLDIntegration(t, vcfText, params)
	lines := strings.Split(strings.TrimSpace(geno), "\n")
	// Header + only the (100,300) pair (r²=1); pairs involving site 200 have r²=0.
	if len(lines) != 2 {
		t.Fatalf("expected header + 1 pair, got %d lines:\n%s", len(lines), geno)
	}
	fields := strings.Split(lines[1], "\t")
	if fields[1] != "100" || fields[2] != "300" {
		t.Errorf("expected pair (100,300), got (%s,%s)", fields[1], fields[2])
	}
}

// TestRun_GenoR2Positions restricts emission to pairs touching a listed
// position. Listing position 200 should emit (100,200) and (200,300), but not
// (100,300).
func TestRun_GenoR2Positions(t *testing.T) {
	posFile := filepath.Join(t.TempDir(), "pos.txt")
	if err := os.WriteFile(posFile, []byte("1\t200\n"), 0o600); err != nil {
		t.Fatalf("write pos file: %v", err)
	}
	vcfText := minimalVCF([]string{"s1", "s2", "s3", "s4"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"1", 200, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"1", 300, []string{"0/0", "0/1", "0/1", "1/1"}},
		})
	params := &Params{GenoR2Positions: posFile, MaxMissing: 1}
	geno, _ := runLDIntegration(t, vcfText, params)
	lines := strings.Split(strings.TrimSpace(geno), "\n")
	if len(lines) != 3 {
		t.Fatalf("expected header + 2 pairs (touching pos 200), got %d:\n%s", len(lines), geno)
	}
	got := map[string]bool{}
	for _, l := range lines[1:] {
		f := strings.Split(l, "\t")
		got[f[1]+"-"+f[2]] = true
	}
	if !got["100-200"] || !got["200-300"] {
		t.Errorf("missing expected pairs, got: %v", got)
	}
	if got["100-300"] {
		t.Errorf("pair (100,300) should be excluded (neither touches pos 200)")
	}
}

// TestRun_HapR2_BasicOutput verifies the .hap.ld file shape and that
// unphased sites are skipped from nChr.
func TestRun_HapR2_BasicOutput(t *testing.T) {
	vcfText := minimalVCF([]string{"s1", "s2", "s3", "s4"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0|0", "0|0", "1|1", "1|1"}},
			{"1", 200, []string{"0|0", "0|0", "1|1", "1|1"}},
		})
	params := &Params{HapR2: true, MaxMissing: 1}
	_, hap := runLDIntegration(t, vcfText, params)
	lines := strings.Split(strings.TrimSpace(hap), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected header + 1 pair, got %d:\n%s", len(lines), hap)
	}
	if !strings.HasPrefix(lines[0], "CHR\tPOS1\tPOS2\tN_CHR\tR^2\tD\tDprime") {
		t.Errorf("missing/incorrect header: %q", lines[0])
	}
	f := strings.Split(lines[1], "\t")
	if len(f) != 7 {
		t.Fatalf("want 7 fields, got %d (%q)", len(f), lines[1])
	}
	if f[3] != "8" {
		t.Errorf("N_CHR=%s want 8", f[3])
	}
	if f[4] != "1" {
		t.Errorf("R^2=%s want 1", f[4])
	}
}

// TestRun_NoLDFlags_NoOutput: when no LD flag is set, no .geno.ld / .hap.ld
// file is produced.
func TestRun_NoLDFlags_NoOutput(t *testing.T) {
	vcfText := minimalVCF([]string{"s1", "s2"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0/0", "0/1"}},
			{"1", 200, []string{"0/0", "0/1"}},
		})
	params := &Params{MaxMissing: 1}
	geno, hap := runLDIntegration(t, vcfText, params)
	if geno != "" || hap != "" {
		t.Errorf("expected no LD output without flags, got geno=%q hap=%q", geno, hap)
	}
}

// TestRun_CrossChromosomePairsSkipped: vcftools does not emit --geno-r2 /
// --hap-r2 across chromosomes (the dedicated --interchrom-geno-r2 /
// --interchrom-hap-r2 options exist for that). Ensure two sites on different
// chromosomes produce no pair.
func TestRun_CrossChromosomePairsSkipped(t *testing.T) {
	vcfText := minimalVCF([]string{"s1", "s2", "s3"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0/0", "0/1", "1/1"}},
			{"2", 100, []string{"0/0", "0/1", "1/1"}},
		})
	params := &Params{GenoR2: true, MaxMissing: 1}
	geno, _ := runLDIntegration(t, vcfText, params)
	lines := strings.Split(strings.TrimSpace(geno), "\n")
	if len(lines) != 1 {
		t.Fatalf("expected only header, got %d lines:\n%s", len(lines), geno)
	}
}

// TestRun_LDWindowMin: --ld-window-min suppresses adjacent SNPs. With min=2,
// only the (100,300) pair survives among three sites.
func TestRun_LDWindowMin(t *testing.T) {
	vcfText := minimalVCF([]string{"s1", "s2", "s3", "s4"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"1", 200, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"1", 300, []string{"0/0", "0/1", "0/1", "1/1"}},
		})
	params := &Params{GenoR2: true, LDWindowMin: 2, MaxMissing: 1}
	geno, _ := runLDIntegration(t, vcfText, params)
	lines := strings.Split(strings.TrimSpace(geno), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected header + 1 pair, got %d:\n%s", len(lines), geno)
	}
	fields := strings.Split(lines[1], "\t")
	if fields[1] != "100" || fields[2] != "300" {
		t.Errorf("expected pair (100,300), got (%s,%s)", fields[1], fields[2])
	}
}

// TestRun_LDWindowBpMin: 200-100=100, 1500-100=1400, 1500-200=1300. With
// LDWindowBpMin=500 only (100,1500) and (200,1500) pass.
func TestRun_LDWindowBpMin(t *testing.T) {
	vcfText := minimalVCF([]string{"s1", "s2", "s3", "s4"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"1", 200, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"1", 1500, []string{"0/0", "0/1", "0/1", "1/1"}},
		})
	params := &Params{GenoR2: true, LDWindowBpMin: 500, MaxMissing: 1}
	geno, _ := runLDIntegration(t, vcfText, params)
	lines := strings.Split(strings.TrimSpace(geno), "\n")
	if len(lines) != 3 {
		t.Fatalf("expected header + 2 pairs, got %d:\n%s", len(lines), geno)
	}
}

// TestExtractLDSite_RejectsEmpty: variant with no samples or no ALT must
// return ok=false (cannot participate in any LD pair).
func TestExtractLDSite_RejectsEmpty(t *testing.T) {
	if _, ok := extractLDSite(&vcf.Variant{Chrom: "1", Pos: 1}); ok {
		t.Error("expected ok=false for empty variant")
	}
	if _, ok := extractLDSite(&vcf.Variant{
		Chrom:   "1",
		Pos:     1,
		Alt:     []string{"G"},
		Samples: []vcf.Sample{},
	}); ok {
		t.Error("expected ok=false for zero samples")
	}
}

// TestExtractLDSite_NoGTField: a sample missing the GT field is marked as
// fully missing (-1 everywhere).
func TestExtractLDSite_NoGTField(t *testing.T) {
	v := &vcf.Variant{
		Chrom:  "1",
		Pos:    1,
		Alt:    []string{"G"},
		Format: []string{"GT"},
		Samples: []vcf.Sample{
			{Name: "s1", Data: map[string]string{}},
			{Name: "s2", Data: map[string]string{"GT": "0/1"}},
		},
	}
	s, ok := extractLDSite(v)
	if !ok {
		t.Fatal("expected ok=true")
	}
	if s.genoCounts[0] != -1 || s.hapA[0] != -1 || s.hapB[0] != -1 {
		t.Errorf("expected sample 0 to be all -1, got %+v / %+v / %+v",
			s.genoCounts[0], s.hapA[0], s.hapB[0])
	}
	if s.genoCounts[1] != 1 {
		t.Errorf("expected sample 1 genoCount=1, got %d", s.genoCounts[1])
	}
}

// TestLDRunner_PruneWindow: with --ld-window=1 only adjacent SNP pairs are
// kept; window should be aggressively pruned so it never grows unbounded.
func TestLDRunner_PruneWindow(t *testing.T) {
	r, err := newLDRunner(&Params{
		OutPrefix: filepath.Join(t.TempDir(), "out"),
		GenoR2:    true,
		LDWindow:  1,
	})
	if err != nil {
		t.Fatalf("newLDRunner: %v", err)
	}
	defer r.close()

	for i := 0; i < 10; i++ {
		v := &vcf.Variant{
			Chrom:  "1",
			Pos:    100 + i,
			Alt:    []string{"G"},
			Format: []string{"GT"},
			Samples: []vcf.Sample{
				{Name: "s1", Data: map[string]string{"GT": "0/0"}},
				{Name: "s2", Data: map[string]string{"GT": "0/1"}},
				{Name: "s3", Data: map[string]string{"GT": "1/1"}},
			},
		}
		r.addVariant(v)
		if len(r.window) > 2 {
			t.Errorf("window grew past 2 (LDWindow=1) at step %d: len=%d", i, len(r.window))
		}
	}
}

// TestRun_GenoAndHap_Together: enabling both flags writes both files.
func TestRun_GenoAndHap_Together(t *testing.T) {
	vcfText := minimalVCF([]string{"s1", "s2", "s3", "s4"},
		[]struct {
			chrom string
			pos   int
			gts   []string
		}{
			{"1", 100, []string{"0|0", "0|0", "1|1", "1|1"}},
			{"1", 200, []string{"0|0", "0|0", "1|1", "1|1"}},
		})
	params := &Params{GenoR2: true, HapR2: true, MaxMissing: 1}
	geno, hap := runLDIntegration(t, vcfText, params)
	if geno == "" {
		t.Errorf("expected .geno.ld output")
	}
	if hap == "" {
		t.Errorf("expected .hap.ld output")
	}
}
