package vcftools

import (
	"bytes"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// TestRelatedness_HandComputed builds a tiny VCF with 4 individuals and 3
// biallelic SNPs and verifies the Yang 2010 unadjusted A_jk values match
// hand-computed expectations.
//
// Setup (rows are samples s1..s4, columns are SNPs):
//
//	SNP1: 0/0  0/1  0/1  1/1   -> x = [0,1,1,2], sumX=4, n=4, p=0.5
//	SNP2: 0/0  0/0  1/1  1/1   -> x = [0,0,2,2], sumX=4, n=4, p=0.5
//	SNP3: 0/0  0/1  1/1  0/1   -> x = [0,1,2,1], sumX=4, n=4, p=0.5
//
// At each SNP p=0.5 so 2p(1-p)=0.5 and the per-pair contribution is
// (x_i - 1)(x_j - 1) / 0.5 = 2(x_i - 1)(x_j - 1).
//
// Pair (s1,s2):
//
//	SNP1: 2*(0-1)*(1-1) = 0
//	SNP2: 2*(0-1)*(0-1) = 2
//	SNP3: 2*(0-1)*(1-1) = 0
//	A_jk = (0 + 2 + 0) / 3 = 2/3
//
// Pair (s1,s4):
//
//	SNP1: 2*(0-1)*(2-1) = -2
//	SNP2: 2*(0-1)*(2-1) = -2
//	SNP3: 2*(0-1)*(1-1) =  0
//	A_jk = -4/3
//
// Diagonal s1: (x1-1)^2 / 0.5 = 2*(x1-1)^2 averaged:
//
//	SNP1: 2*(0-1)^2 = 2
//	SNP2: 2*(0-1)^2 = 2
//	SNP3: 2*(0-1)^2 = 2
//	A_jj = 2
func TestRelatedness_HandComputed(t *testing.T) {
	vcfText := buildMinimalVCFRel(t, []string{"s1", "s2", "s3", "s4"},
		[]relRow{
			{"1", 100, []string{"0/0", "0/1", "0/1", "1/1"}},
			{"1", 200, []string{"0/0", "0/0", "1/1", "1/1"}},
			{"1", 300, []string{"0/0", "0/1", "1/1", "0/1"}},
		})
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, Relatedness: true}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, err := os.ReadFile(prefix + ".relatedness")
	if err != nil {
		t.Fatalf("read .relatedness: %v", err)
	}
	got := parseRelatedness(t, string(data))
	want := map[[2]string]float64{
		// diagonal: average of (x_i - 2p)^2/(2p(1-p)) with p=0.5
		{"s1", "s1"}: 2.0,       // 2 + 2 + 2 / 3 = 2
		{"s2", "s2"}: 2.0 / 3.0, // 0 + 2 + 0 / 3
		{"s3", "s3"}: 4.0 / 3.0, // 0 + 2 + 2 / 3
		{"s4", "s4"}: 4.0 / 3.0, // 2 + 2 + 0 / 3
		// off-diagonal
		{"s1", "s2"}: 2.0 / 3.0,
		{"s1", "s3"}: -4.0 / 3.0, // 0 + -2 + -2 / 3
		{"s1", "s4"}: -4.0 / 3.0,
		{"s2", "s3"}: -2.0 / 3.0,
		{"s2", "s4"}: -2.0 / 3.0,
		{"s3", "s4"}: 2.0 / 3.0,
	}
	if len(got) != len(want) {
		t.Errorf("got %d pairs, want %d (output:\n%s)", len(got), len(want), data)
	}
	for k, v := range want {
		g, ok := got[k]
		if !ok {
			t.Errorf("missing pair %v", k)
			continue
		}
		// Output is C++ defaultfloat (6 significant digits), so compare
		// with a 6-sig-fig tolerance rather than full double precision.
		if math.Abs(g-v) > 1e-5 {
			t.Errorf("pair %v: got %g, want %g", k, g, v)
		}
	}
}

// TestRelatedness_SkipsMonomorphic ensures a SNP that is monomorphic across
// the kept individuals (p=0 or p=1) does NOT contribute to the average. With
// no informative SNPs, every pair has N_sites == 0, so the relatedness is a
// positive 0.0/0.0 NaN, which libstdc++ (and so the upstream binary) prints
// as "nan" — verified against the live `vcftools --relatedness` oracle. We
// assert the output is all "nan" rather than numeric values, matching
// upstream. (This previously asserted "-nan", which encoded an earlier
// hard-coded NaN sign rather than upstream's actual output.)
func TestRelatedness_SkipsMonomorphic(t *testing.T) {
	vcfText := buildMinimalVCFRel(t, []string{"s1", "s2"},
		[]relRow{
			{"1", 100, []string{"0/0", "0/0"}},
			{"1", 200, []string{"1/1", "1/1"}},
		})
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, Relatedness: true}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, _ := os.ReadFile(prefix + ".relatedness")
	for _, ln := range strings.Split(strings.TrimSpace(string(data)), "\n")[1:] {
		fields := strings.Split(ln, "\t")
		if len(fields) != 3 {
			continue
		}
		if fields[2] != "nan" {
			t.Errorf("expected nan for monomorphic-only data, got %q in line %q", fields[2], ln)
		}
	}
}

// TestRelatedness_SkipsMultiAllelic confirms multi-allelic sites are ignored.
func TestRelatedness_SkipsMultiAllelic(t *testing.T) {
	// SNP1 is biallelic -> contributes; SNP2 has 2 ALTs -> skipped.
	vcfText := "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n" +
		"1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\n" +
		"1\t200\t.\tA\tG,C\t.\tPASS\t.\tGT\t0/1\t1/2\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, Relatedness: true}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, _ := os.ReadFile(prefix + ".relatedness")
	// With only SNP1 contributing: p=0.5, 2p(1-p)=0.5.
	// s1: (0-1)*(0-1)/0.5 = 2; s2: (2-1)*(2-1)/0.5 = 2;
	// (s1,s2): (0-1)*(2-1)/0.5 = -2.
	got := parseRelatedness(t, string(data))
	if math.Abs(got[[2]string{"s1", "s2"}]+2) > 1e-9 {
		t.Errorf("got %v want -2", got[[2]string{"s1", "s2"}])
	}
	if math.Abs(got[[2]string{"s1", "s1"}]-2) > 1e-9 {
		t.Errorf("got %v want 2", got[[2]string{"s1", "s1"}])
	}
}

const lrohHeader = "CHROM\tAUTO_START\tAUTO_END\tMIN_START\tMAX_END\tN_VARIANTS_BETWEEN_MAX_BOUNDARIES\tN_MISMATCHES\tINDV"

// TestLROH_HMMDetectsRun checks the forward-backward HMM reports an autozygous
// region for an individual that is homozygous across a tight cluster of sites
// (with other samples het, so the site heterozygosity h makes the homozygous
// emissions favour the autozygous state) and reports none for a het-rich one.
func TestLROH_HMMDetectsRun(t *testing.T) {
	var rows []relRow
	pos := 1000
	for i := 0; i < 60; i++ {
		pos += 200
		rows = append(rows, relRow{"1", pos, []string{"1/1", "0/1", "0/1"}})
	}
	vcfText := buildMinimalVCFRel(t, []string{"s1", "s2", "s3"}, rows)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	if err := Run(strings.NewReader(vcfText), &Params{OutPrefix: prefix, LROH: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, err := os.ReadFile(prefix + ".LROH")
	if err != nil {
		t.Fatalf("read .LROH: %v", err)
	}
	lines := strings.Split(strings.TrimRight(string(data), "\n"), "\n")
	if lines[0] != lrohHeader {
		t.Errorf("bad header: %q", lines[0])
	}
	s1Rows, s2Rows := 0, 0
	for _, ln := range lines[1:] {
		cols := strings.Split(ln, "\t")
		if len(cols) != 8 {
			t.Errorf("row should have 8 columns: %q", ln)
			continue
		}
		switch cols[7] {
		case "s1":
			s1Rows++
		case "s2":
			s2Rows++
		}
	}
	if s1Rows == 0 {
		t.Errorf("expected at least one autozygous region for s1:\n%s", data)
	}
	if s2Rows != 0 {
		t.Errorf("het-rich s2 should have no autozygous region, got %d:\n%s", s2Rows, data)
	}
}

// TestRelatedness2_KingFormula tests the KING-robust kinship for a
// hand-computed pair of "twins" (identical genotypes) and a pair of
// "unrelated" individuals (every site differs).
//
// Setup: 3 samples, 4 biallelic SNPs.
//
//	SNP1: 0/0  0/0  1/1
//	SNP2: 0/1  0/1  0/0
//	SNP3: 1/1  1/1  0/1
//	SNP4: 0/1  0/1  1/1
//
// Upstream's KING-robust estimator is
// phi = (N_AaAa - 2*N_AAaa) / (N_Aa[i] + N_Aa[j]).
//
// Pair (s1,s2) = identical at every site: N_AaAa=2 (SNPs 2,4),
// N_AAaa=0, N_Aa[s1]=2, N_Aa[s2]=2. phi = (2 - 0)/(2 + 2) = 0.5.
//
// Pair (s1,s3): genotypes (0/0,1/1), (0/1,0/0), (1/1,0/1), (0/1,1/1).
//
//	N_AaAa = 0 (no site where both are het)
//	N_AAaa = 1 (SNP1: 0/0 vs 1/1)
//	N_Aa[s1] = 2 (SNPs 2 and 4), N_Aa[s3] = 1 (SNP3)
//	phi = (0 - 2*1)/(2 + 1) = -2/3 = -0.666667
func TestRelatedness2_KingFormula(t *testing.T) {
	rows := []relRow{
		{"1", 100, []string{"0/0", "0/0", "1/1"}},
		{"1", 200, []string{"0/1", "0/1", "0/0"}},
		{"1", 300, []string{"1/1", "1/1", "0/1"}},
		{"1", 400, []string{"0/1", "0/1", "1/1"}},
	}
	vcfText := buildMinimalVCFRel(t, []string{"s1", "s2", "s3"}, rows)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, Relatedness2: true}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, err := os.ReadFile(prefix + ".relatedness2")
	if err != nil {
		t.Fatalf("read .relatedness2: %v", err)
	}
	lines := strings.Split(strings.TrimRight(string(data), "\n"), "\n")
	if lines[0] != "INDV1\tINDV2\tN_AaAa\tN_AAaa\tN1_Aa\tN2_Aa\tRELATEDNESS_PHI" {
		t.Errorf("bad header: %q", lines[0])
	}
	got := make(map[[2]string][]string)
	for _, ln := range lines[1:] {
		fields := strings.Split(ln, "\t")
		got[[2]string{fields[0], fields[1]}] = fields[2:]
	}
	// (s1,s2): phi=0.5
	g, ok := got[[2]string{"s1", "s2"}]
	if !ok {
		t.Fatalf("missing pair (s1,s2): %v", got)
	}
	phi, _ := strconv.ParseFloat(g[4], 64)
	if math.Abs(phi-0.5) > 1e-9 {
		t.Errorf("(s1,s2) phi=%g, want 0.5", phi)
	}
	// (s1,s3): phi=-2/3
	g = got[[2]string{"s1", "s3"}]
	phi, _ = strconv.ParseFloat(g[4], 64)
	if math.Abs(phi-(-2.0/3.0)) > 1e-5 {
		t.Errorf("(s1,s3) phi=%g, want -0.666667", phi)
	}
	// Self phi = 0.5.
	g = got[[2]string{"s1", "s1"}]
	phi, _ = strconv.ParseFloat(g[4], 64)
	if math.Abs(phi-0.5) > 1e-9 {
		t.Errorf("(s1,s1) phi=%g, want 0.5", phi)
	}
}

// TestPhasedBlocks_Basic: emit a single run for s1 (3 phased sites) and
// nothing for s2 (mixed phased/unphased breaks the run too short).
func TestPhasedBlocks_Basic(t *testing.T) {
	rows := []relRow{
		{"1", 100, []string{"0|0", "0/0"}},
		{"1", 200, []string{"0|1", "0|0"}},
		{"1", 300, []string{"1|1", "0|1"}},
	}
	vcfText := buildMinimalVCFRel(t, []string{"s1", "s2"}, rows)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	params := &Params{OutPrefix: prefix, PhasedBlocks: true}
	if err := Run(strings.NewReader(vcfText), params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, err := os.ReadFile(prefix + ".blocks")
	if err != nil {
		t.Fatalf("read .blocks: %v", err)
	}
	lines := strings.Split(strings.TrimRight(string(data), "\n"), "\n")
	if lines[0] != "CHROM\tBLOCK_START\tBLOCK_END\tN_VARIANTS\tINDV" {
		t.Errorf("bad header: %q", lines[0])
	}
	// s1: a single 3-site run (100..300). s2: 2-site run (200..300).
	wantRows := map[string]bool{
		"1\t100\t300\t3\ts1": true,
		"1\t200\t300\t2\ts2": true,
	}
	for _, ln := range lines[1:] {
		if !wantRows[ln] {
			t.Errorf("unexpected row %q", ln)
		}
		delete(wantRows, ln)
	}
	for r := range wantRows {
		t.Errorf("missing row %q", r)
	}
}

// TestIsPhasedDiploid covers the per-call phased-detection helper.
func TestIsPhasedDiploid(t *testing.T) {
	cases := []struct {
		gt   string
		want bool
	}{
		{"0|0", true},
		{"1|0", true},
		{"0/1", false},
		{".|0", false},
		{"0|.", false},
		{".", false},
		{"", false},
		{"0|0|0", true}, // tolerant: stops at first '|'
	}
	for _, tc := range cases {
		if got := isPhasedDiploid(tc.gt); got != tc.want {
			t.Errorf("isPhasedDiploid(%q)=%v want %v", tc.gt, got, tc.want)
		}
	}
}

// TestLROH_NoRunOnHighHet verifies the HMM reports no autozygous region (only
// the header) when every individual is heterozygous everywhere — the
// non-autozygous emissions dominate.
func TestLROH_NoRunOnHighHet(t *testing.T) {
	var rows []relRow
	pos := 1000
	for i := 0; i < 40; i++ {
		pos += 500
		rows = append(rows, relRow{"1", pos, []string{"0/1", "0/1"}})
	}
	vcfText := buildMinimalVCFRel(t, []string{"s1", "s2"}, rows)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	if err := Run(strings.NewReader(vcfText), &Params{OutPrefix: prefix, LROH: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	data, _ := os.ReadFile(prefix + ".LROH")
	lines := strings.Split(strings.TrimRight(string(data), "\n"), "\n")
	if len(lines) != 1 || lines[0] != lrohHeader {
		t.Errorf("expected header-only output, got %d lines:\n%s", len(lines), data)
	}
}

// relRow is a tiny helper struct local to this test file.
type relRow struct {
	chrom string
	pos   int
	gts   []string
}

func buildMinimalVCFRel(t *testing.T, samples []string, rows []relRow) string {
	t.Helper()
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

func parseRelatedness(t *testing.T, data string) map[[2]string]float64 {
	t.Helper()
	out := make(map[[2]string]float64)
	lines := strings.Split(strings.TrimRight(data, "\n"), "\n")
	if len(lines) < 1 || lines[0] != "INDV1\tINDV2\tRELATEDNESS_AJK" {
		t.Fatalf("bad relatedness header: %q", lines[0])
	}
	for _, ln := range lines[1:] {
		fields := strings.Split(ln, "\t")
		if len(fields) != 3 {
			t.Fatalf("bad row %q", ln)
		}
		v, err := strconv.ParseFloat(fields[2], 64)
		if err != nil {
			t.Fatalf("parse %q: %v", fields[2], err)
		}
		out[[2]string{fields[0], fields[1]}] = v
	}
	return out
}
