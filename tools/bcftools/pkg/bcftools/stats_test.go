package bcftools

import (
	"bytes"
	"io"
	"math"
	"os"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// 10-record fixture mixing SNPs, MNPs, indels, and a multi-allelic site.
const statsFixtureVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FILTER=<ID=q10,Description="below 10">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	.	A	G	30	PASS	DP=30;AF=0.05	GT:DP	0/0:10	0/1:10	0/1:10
chr1	150	.	C	T	60	PASS	DP=60;AF=0.15	GT:DP	0/1:20	1/1:20	0/0:20
chr1	200	.	A	C	120	PASS	DP=45;AF=0.95	GT:DP	1/1:15	1/1:15	0/1:15
chr1	250	.	G	T	250	PASS	DP=30;AF=0.5	GT:DP	0/0:10	0/0:10	0/1:10
chr1	300	.	ATG	A	600	PASS	DP=60;AF=0.3	GT:DP	0/1:20	0/0:20	0/0:20
chr1	350	.	A	AT	1200	PASS	DP=30;AF=0.15	GT:DP	0/1:10	0/1:10	0/0:10
chr1	400	.	ACGT	A	80	PASS	DP=45;AF=0.05	GT:DP	0/1:15	0/0:15	0/0:15
chr1	450	.	AC	GT	40	PASS	DP=30;AF=0.5	GT:DP	0/0:10	0/1:10	0/1:10
chr1	500	.	A	AT,ATT	75	q10	DP=60;AF=0.1,0.2	GT:DP	0/1:20	0/2:20	1/2:20
chr2	100	.	C	A	300	PASS	DP=30;AF=0.1	GT:DP	0/0:10	0/1:10	0/0:10
`

func runStats(t *testing.T, input string, opts StatsOptions) (*statsResult, string) {
	t.Helper()
	var out bytes.Buffer
	res, err := Stats(strings.NewReader(input), &out, opts)
	if err != nil {
		t.Fatalf("Stats: %v", err)
	}
	return res, out.String()
}

// helper that hand-counts the fixture sections.
func TestStatsSNHandCount(t *testing.T) {
	res, out := runStats(t, statsFixtureVCF, StatsOptions{InputFile: "test.vcf"})
	// records: 10 total.
	if res.numRecords != 10 {
		t.Errorf("numRecords = %d, want 10", res.numRecords)
	}
	if res.numNoALTs != 0 {
		t.Errorf("numNoALTs = %d, want 0", res.numNoALTs)
	}
	// SNPs: chr1@100 A>G, chr1@150 C>T, chr1@200 A>C, chr1@250 G>T, chr2@100 C>A = 5.
	if res.numSNPs != 5 {
		t.Errorf("numSNPs = %d, want 5", res.numSNPs)
	}
	// MNPs: chr1@450 AC>GT = 1.
	if res.numMNPs != 1 {
		t.Errorf("numMNPs = %d, want 1", res.numMNPs)
	}
	// indels: chr1@300, chr1@350, chr1@400, chr1@500 = 4.
	if res.numIndels != 4 {
		t.Errorf("numIndels = %d, want 4", res.numIndels)
	}
	if res.numMA != 1 {
		t.Errorf("numMA = %d, want 1", res.numMA)
	}
	if res.numMASNP != 0 {
		t.Errorf("numMASNP = %d, want 0", res.numMASNP)
	}
	// Verify the always-present section headers appear in the text. The
	// per-sample sections (PSC/PSI/HWE) are emitted only with -s/-S, mirroring
	// upstream, so they are absent here.
	for _, hdr := range []string{"# SN,", "# AF,", "# QUAL,", "# IDD,", "# ST,", "# DP,"} {
		if !strings.Contains(out, hdr) {
			t.Errorf("missing section header %q in:\n%s", hdr, out)
		}
	}
	// Match the per-section header text exactly as upstream emits it: PSC and
	// PSI carry verbose comma-prefixed comments, while HWE is the bare "# HWE".
	for _, hdr := range []string{"# PSC, Per-sample counts", "# PSI, Per-Sample Indels", "# HWE\n"} {
		if strings.Contains(out, hdr) {
			t.Errorf("unexpected per-sample section %q without -s:\n%s", hdr, out)
		}
	}
	// With -s -, the per-sample sections should appear.
	_, outS := runStats(t, statsFixtureVCF, StatsOptions{InputFile: "test.vcf", EnableSamples: true, Samples: []string{"-"}})
	for _, hdr := range []string{"# PSC, Per-sample counts", "# PSI, Per-Sample Indels", "# HWE\n"} {
		if !strings.Contains(outS, hdr) {
			t.Errorf("missing per-sample section %q with -s -:\n%s", hdr, outS)
		}
	}
	if !strings.Contains(out, "number of samples:\t3") {
		t.Errorf("expected 3 samples line, got:\n%s", out)
	}
}

func TestStatsIDDLengths(t *testing.T) {
	res, _ := runStats(t, statsFixtureVCF, StatsOptions{})
	want := map[int]int{
		-2: 1, // chr1@300 ATG>A  (delta -2)
		1:  1, // chr1@350 A>AT   (delta +1)
		-3: 1, // chr1@400 ACGT>A (delta -3)
		2:  1, // chr1@500 A>ATT  (delta +2 — the 2nd alt of the multi-allele site)
		// Note chr1@500 also has A>AT giving +1, so +1 is bumped twice.
	}
	want[1] = 2
	for length, expected := range want {
		if got := res.indelLen[length]; got != expected {
			t.Errorf("indelLen[%d] = %d, want %d", length, got, expected)
		}
	}
	// Should have exactly 4 distinct length bins for this fixture.
	if len(res.indelLen) != 4 {
		t.Errorf("expected 4 indel-length bins, got %d (%v)", len(res.indelLen), res.indelLen)
	}
}

func TestStatsSTSubstitutions(t *testing.T) {
	res, _ := runStats(t, statsFixtureVCF, StatsOptions{})
	// Hand-counted per-allele substitutions:
	//   A>G (chr1@100), C>T (chr1@150), A>C (chr1@200), G>T (chr1@250), C>A (chr2@100).
	want := map[string]int{
		"A>G": 1,
		"C>T": 1,
		"A>C": 1,
		"G>T": 1,
		"C>A": 1,
	}
	for k, expected := range want {
		if got := res.subst[k]; got != expected {
			t.Errorf("subst[%s] = %d, want %d", k, got, expected)
		}
	}
	// All others should be zero.
	for _, k := range []string{"A>T", "C>G", "G>A", "G>C", "T>A", "T>C", "T>G"} {
		if got := res.subst[k]; got != 0 {
			t.Errorf("subst[%s] = %d, want 0", k, got)
		}
	}
}

func TestStatsAFBinning(t *testing.T) {
	// AF is derived from the allele counts (bcf_calc_ac), not INFO/AF.
	// Each single-sample 0/1 site has AC=1, AN=2; an AC==1 allele is a
	// singleton and lands in bucket 0 (the singleton bin), regardless of the
	// INFO/AF tag value.
	mini := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	A	G	30	PASS	AF=0.05	GT	0/1
chr1	200	.	C	T	30	PASS	AF=0.15	GT	0/1
chr1	300	.	G	A	30	PASS	AF=0.95	GT	0/1
`
	res, _ := runStats(t, mini, StatsOptions{})
	if res.afSNPs[0] != 3 {
		t.Errorf("AF singleton bucket SNPs = %d, want 3", res.afSNPs[0])
	}
	// A non-singleton, non-INFO/AF example: AC=2 over AN=4 → af=0.5.
	mini2 := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	30	PASS	.	GT	0/1	0/1
`
	res2, _ := runStats(t, mini2, StatsOptions{})
	// af=0.5 → bucket int(0.5*99)+1 = 50.
	if res2.afSNPs[50] != 1 {
		t.Errorf("AF bucket 50 SNPs = %d, want 1", res2.afSNPs[50])
	}
}

func TestStatsQUALBinning(t *testing.T) {
	res, _ := runStats(t, statsFixtureVCF, StatsOptions{})
	// QUAL is bucketed as 1 + int(qual*10). The first-ALT ts/tv split is what
	// the QUAL section reports. SNP records and their first-ALT QUAL buckets:
	//   30 (A>G ts), 60 (C>T ts), 120 (A>C tv), 250 (G>T tv), 300 (C>A tv).
	wantTs := map[int]int{1 + 30*10: 1, 1 + 60*10: 1}
	wantTv := map[int]int{1 + 120*10: 1, 1 + 250*10: 1, 1 + 300*10: 1}
	for q, c := range wantTs {
		if got := res.qualTs[q]; got != c {
			t.Errorf("qualTs[%d] = %d, want %d", q, got, c)
		}
	}
	for q, c := range wantTv {
		if got := res.qualTv[q]; got != c {
			t.Errorf("qualTv[%d] = %d, want %d", q, got, c)
		}
	}
}

func TestStatsPSCThreeSampleFixture(t *testing.T) {
	fx := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	.	A	G	30	PASS	DP=30	GT:DP	0/0:10	0/1:10	1/1:10
chr1	200	.	C	T	30	PASS	DP=30	GT:DP	0/1:5	0/1:5	0/0:5
chr1	300	.	G	A	30	PASS	DP=30	GT:DP	1/1:8	0/0:8	0/1:8
chr1	400	.	A	AT	30	PASS	DP=30	GT:DP	0/0:6	0/1:6	1/1:6
chr1	500	.	C	G	30	PASS	DP=30	GT:DP	0/1:4	1/1:4	0/1:4
`
	res, _ := runStats(t, fx, StatsOptions{EnableSamples: true, Samples: []string{"-"}})
	// nRefHom/nNonRefHom/nHets count SNP genotypes only; the chr1@400 A>AT
	// indel genotype is excluded from these (it feeds nIndels / PSI instead).
	// The expectations below match the live upstream `bcftools stats -s -`
	// output for this fixture.
	want := []struct{ refHom, altHom, het int }{
		{2, 1, 2},
		{1, 1, 2},
		{1, 1, 2},
	}
	for i, w := range want {
		if res.pscNRefHom[i] != w.refHom {
			t.Errorf("S%d nRefHom = %d, want %d", i+1, res.pscNRefHom[i], w.refHom)
		}
		if res.pscNNonRefHom[i] != w.altHom {
			t.Errorf("S%d nNonRefHom = %d, want %d", i+1, res.pscNNonRefHom[i], w.altHom)
		}
		if res.pscNHets[i] != w.het {
			t.Errorf("S%d nHets = %d, want %d", i+1, res.pscNHets[i], w.het)
		}
	}
}

func TestStatsSampleRestriction(t *testing.T) {
	res, _ := runStats(t, statsFixtureVCF, StatsOptions{EnableSamples: true, Samples: []string{"S1", "S3"}})
	if len(res.samples) != 2 || res.samples[0] != "S1" || res.samples[1] != "S3" {
		t.Errorf("samples = %v, want [S1 S3]", res.samples)
	}
	// PSC arrays should have length 2 (only restricted samples counted).
	if len(res.pscNRefHom) != 2 {
		t.Errorf("len(pscNRefHom) = %d, want 2", len(res.pscNRefHom))
	}
}

func TestStatsIncludeExpression(t *testing.T) {
	resAll, _ := runStats(t, statsFixtureVCF, StatsOptions{})
	resPass, _ := runStats(t, statsFixtureVCF, StatsOptions{IncludeExpr: `FILTER="PASS"`})
	// The fixture has 1 non-PASS record (chr1@500 q10) — filtering with -i 'FILTER="PASS"' drops it.
	if resPass.numRecords != resAll.numRecords-1 {
		t.Errorf("expected -i FILTER=PASS to drop the q10 record (got %d, all=%d)",
			resPass.numRecords, resAll.numRecords)
	}
}

func TestStatsExcludeExpression(t *testing.T) {
	res, _ := runStats(t, statsFixtureVCF, StatsOptions{ExcludeExpr: `FILTER="q10"`})
	if res.numRecords != 9 {
		t.Errorf("with -e FILTER=q10 expected 9 records, got %d", res.numRecords)
	}
}

func TestStatsApplyFilters(t *testing.T) {
	res, _ := runStats(t, statsFixtureVCF, StatsOptions{ApplyFilters: []string{"PASS"}})
	if res.numRecords != 9 {
		t.Errorf("with -f PASS expected 9 records, got %d", res.numRecords)
	}
}

func TestStatsRegionRestriction(t *testing.T) {
	// -r chr1:100-200 should keep records at chr1@100, chr1@150, chr1@200.
	res, _ := runStats(t, statsFixtureVCF, StatsOptions{Regions: []string{"chr1:100-200"}})
	if res.numRecords != 3 {
		t.Errorf("expected 3 records in chr1:100-200, got %d", res.numRecords)
	}
}

func TestStatsCommandLineHeader(t *testing.T) {
	_, out := runStats(t, statsFixtureVCF, StatsOptions{InputFile: "myfile.vcf"})
	if !strings.Contains(out, "command line was:") {
		t.Errorf("missing command-line header:\n%s", out)
	}
	if !strings.Contains(out, "myfile.vcf") {
		t.Errorf("InputFile not echoed:\n%s", out)
	}
}

func TestStatsParseDepthSpec(t *testing.T) {
	tests := []struct {
		in             string
		min, max, step int
		err            bool
	}{
		{"", 0, 500, 1, false},
		{"0,100,2", 0, 100, 2, false},
		{",,5", 0, 500, 5, false},
		{"1,2", 0, 0, 0, true},
		{"x,2,1", 0, 0, 0, true},
	}
	for _, tt := range tests {
		min, max, step, err := ParseDepthSpec(tt.in)
		if tt.err {
			if err == nil {
				t.Errorf("ParseDepthSpec(%q): expected error", tt.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseDepthSpec(%q): %v", tt.in, err)
			continue
		}
		if min != tt.min || max != tt.max || step != tt.step {
			t.Errorf("ParseDepthSpec(%q) = (%d,%d,%d), want (%d,%d,%d)",
				tt.in, min, max, step, tt.min, tt.max, tt.step)
		}
	}
}

func TestStatsParseAFBins(t *testing.T) {
	bins, err := ParseAFBins("0,0.1,0.5,1.0")
	if err != nil {
		t.Fatalf("ParseAFBins: %v", err)
	}
	want := []float64{0, 0.1, 0.5, 1.0}
	if len(bins) != len(want) {
		t.Fatalf("len = %d, want %d", len(bins), len(want))
	}
	for i := range bins {
		if math.Abs(bins[i]-want[i]) > 1e-9 {
			t.Errorf("bin[%d] = %f, want %f", i, bins[i], want[i])
		}
	}
	if _, err := ParseAFBins("0"); err == nil {
		t.Errorf("expected error for single-bin AF spec")
	}
	if _, err := ParseAFBins("bad,0.5"); err == nil {
		t.Errorf("expected error for non-numeric AF bin")
	}
}

func TestStatsAFBinEdges(t *testing.T) {
	bins := []float64{0.0, 0.1, 0.5, 1.0}
	cases := map[float64]int{
		-0.1: 0,
		0.0:  0,
		0.05: 0,
		0.1:  1,
		0.49: 1,
		0.5:  2,
		0.99: 2,
		1.0:  2,
		1.5:  2,
	}
	for af, want := range cases {
		if got := afBinIndex(bins, af); got != want {
			t.Errorf("afBinIndex(%f) = %d, want %d", af, got, want)
		}
	}
}

func TestStatsClassifyVariant(t *testing.T) {
	tests := []struct {
		ref, alt, want string
	}{
		{"A", "G", "snp"},
		{"AC", "GT", "mnp"},
		{"A", "AT", "indel"},
		{"AT", "A", "indel"},
		{"A", "<DEL>", "other"},
		{"A", ".", "other"},
		{"A", "*", "other"},
		{"N", "A", "other"}, // N is not A/C/G/T
	}
	for _, tt := range tests {
		if got := classifyVariant(tt.ref, tt.alt); got != tt.want {
			t.Errorf("classifyVariant(%q,%q) = %q, want %q", tt.ref, tt.alt, got, tt.want)
		}
	}
}

func TestStatsTransitionType(t *testing.T) {
	cases := map[string]string{
		"A>G": "ts", "G>A": "ts", "C>T": "ts", "T>C": "ts",
		"A>C": "tv", "A>T": "tv", "G>T": "tv", "C>G": "tv",
	}
	for k, want := range cases {
		parts := strings.Split(k, ">")
		if got := transitionType(parts[0], parts[1]); got != want {
			t.Errorf("transitionType(%s) = %s, want %s", k, got, want)
		}
	}
}

func TestStatsDPBin(t *testing.T) {
	opts := StatsOptions{DepthMin: 0, DepthMax: 100, DepthStep: 10}
	cases := map[int]int{
		-5:  0,
		0:   0,
		9:   0,
		10:  10,
		11:  10,
		99:  90,
		100: 100,
		200: 100,
	}
	for in, want := range cases {
		if got := dpBin(opts, in); got != want {
			t.Errorf("dpBin(%d) = %d, want %d", in, got, want)
		}
	}
	// Defaults when step <= 0.
	def := StatsOptions{}
	if got := dpBin(def, 5); got != 5 {
		t.Errorf("default dpBin(5) = %d, want 5", got)
	}
}

func TestStatsFilterSampleSet(t *testing.T) {
	got := filterSampleSet([]string{"a", "b", "c"}, []string{"c", "a"})
	if len(got) != 2 || got[0] != "c" || got[1] != "a" {
		t.Errorf("filterSampleSet returned %v", got)
	}
	got = filterSampleSet([]string{"a", "b"}, nil)
	if len(got) != 2 || got[0] != "a" || got[1] != "b" {
		t.Errorf("nil filter returned %v", got)
	}
}

func TestStatsParseDiploidGT(t *testing.T) {
	cases := []struct {
		gt   string
		a, b int8
		ok   bool
	}{
		{"0/0", 0, 0, true},
		{"0|1", 0, 1, true},
		{"1/2", 1, 2, true},
		{"1", 1, 1, true},
		{".", 0, 0, false},
		{"./.", 0, 0, false},
		{"", 0, 0, false},
		{"x/y", 0, 0, false},
	}
	for _, tt := range cases {
		a, b, _, ok := parseDiploidGT(tt.gt)
		if ok != tt.ok {
			t.Errorf("parseDiploidGT(%q): ok=%v, want %v", tt.gt, ok, tt.ok)
			continue
		}
		if !ok {
			continue
		}
		if a != tt.a || b != tt.b {
			t.Errorf("parseDiploidGT(%q) = (%d,%d), want (%d,%d)", tt.gt, a, b, tt.a, tt.b)
		}
	}
}

func TestStatsHWE(t *testing.T) {
	// Perfect HWE: p=0.5, n=4 → AA=1, Aa=2, aa=1. The het fraction is 2/4 =
	// 0.5; with AC=2/AN=8... here AN=8, AC=4 → af=0.5 → first-ALT AF bucket
	// 50. The het-fraction histogram bin is int(0.5*(nafHWE-1)) = 49.
	mini := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	100	.	A	G	30	PASS	.	GT	0/0	0/1	0/1	1/1
`
	res, _ := runStats(t, mini, StatsOptions{EnableSamples: true, Samples: []string{"-"}})
	bucket := 50 // int(0.5*99)+1
	ihet := int(0.5 * float32(res.nafHWE-1))
	if got := res.afHWE[bucket*res.nafHWE+ihet]; got != 1 {
		t.Errorf("afHWE[bucket=%d,ihet=%d] = %d, want 1", bucket, ihet, got)
	}
}

func TestStatsAFTag(t *testing.T) {
	// With --af-tag we read the named INFO tag instead of deriving AF from the
	// genotypes. The singleton rule does not apply to --af-tag values.
	mini := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##INFO=<ID=MAF,Number=A,Type=Float,Description="MAF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	A	G	30	PASS	AF=0.95;MAF=0.05	GT	0/1
`
	res, _ := runStats(t, mini, StatsOptions{AFTag: "MAF"})
	// MAF=0.05 → bucket int(0.05*99)+1 = 5.
	if res.afSNPs[5] != 1 {
		t.Errorf("AF bucket 5 should have 1 SNP (AFTag=MAF), got %d", res.afSNPs[5])
	}
}

func TestStatsNoALTRecord(t *testing.T) {
	mini := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	A	.	30	PASS	.	GT	0/0
chr1	200	.	C	T	30	PASS	.	GT	0/1
`
	res, _ := runStats(t, mini, StatsOptions{})
	if res.numNoALTs != 1 {
		t.Errorf("numNoALTs = %d, want 1", res.numNoALTs)
	}
	if res.numSNPs != 1 {
		t.Errorf("numSNPs = %d, want 1", res.numSNPs)
	}
}

func TestStatsFirstAlleleOnly(t *testing.T) {
	mini := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	A	G,C	30	PASS	.	GT	0/1
`
	// Without --1st-allele-only: A>G and A>C both counted.
	resAll, _ := runStats(t, mini, StatsOptions{})
	if resAll.subst["A>G"]+resAll.subst["A>C"] != 2 {
		t.Errorf("multi-allelic site should produce 2 substitutions, got A>G=%d A>C=%d",
			resAll.subst["A>G"], resAll.subst["A>C"])
	}
	// With --1st-allele-only: only A>G.
	resFirst, _ := runStats(t, mini, StatsOptions{FirstAlleleOnly: true})
	if resFirst.subst["A>G"] != 1 || resFirst.subst["A>C"] != 0 {
		t.Errorf("1st-allele-only: A>G=%d, A>C=%d, want 1,0",
			resFirst.subst["A>G"], resFirst.subst["A>C"])
	}
}

func TestStatsMultiAllelicSNP(t *testing.T) {
	mini := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	A	G,C	30	PASS	.	GT	0/1
`
	res, _ := runStats(t, mini, StatsOptions{})
	if res.numMA != 1 || res.numMASNP != 1 {
		t.Errorf("MA=%d MASNP=%d, want 1,1", res.numMA, res.numMASNP)
	}
}

func TestStatsStatsFileVCF(t *testing.T) {
	// Round-trip Stats through StatsFile by writing a temp VCF.
	dir := t.TempDir()
	path := dir + "/in.vcf"
	if err := writeTestFile(path, statsFixtureVCF); err != nil {
		t.Fatalf("write fixture: %v", err)
	}
	var buf bytes.Buffer
	res, err := StatsFile(path, &buf, StatsOptions{})
	if err != nil {
		t.Fatalf("StatsFile: %v", err)
	}
	if res.numRecords != 10 {
		t.Errorf("numRecords = %d, want 10", res.numRecords)
	}
	if !strings.Contains(buf.String(), path) {
		t.Errorf("InputFile not auto-filled to %q:\n%s", path, buf.String())
	}
}

func TestStatsFileNotFound(t *testing.T) {
	var buf bytes.Buffer
	if _, err := StatsFile("/nonexistent/file.vcf", &buf, StatsOptions{}); err == nil {
		t.Errorf("expected error for missing file")
	}
}

func TestStatsBadIncludeExpr(t *testing.T) {
	var buf bytes.Buffer
	_, err := Stats(strings.NewReader(statsFixtureVCF), &buf, StatsOptions{IncludeExpr: "INFO/x !! 1"})
	if err == nil {
		t.Errorf("expected error parsing bad include expression")
	}
}

func TestStatsBadRegion(t *testing.T) {
	var buf bytes.Buffer
	_, err := Stats(strings.NewReader(statsFixtureVCF), &buf, StatsOptions{Regions: []string{"chr1:bad"}})
	if err == nil {
		t.Errorf("expected error parsing bad region")
	}
}

func TestStatsDepthSpecRoundtrip(t *testing.T) {
	res, _ := runStats(t, statsFixtureVCF, StatsOptions{DepthMin: 0, DepthMax: 100, DepthStep: 10})
	// dpSites is keyed by the idist bucket index. Site DP=30 with min=0,
	// step=10 maps to bucket 1+(30-0)/10 = 4 (label "30").
	if got := res.dpSites[res.dpIdx(30)]; got < 1 {
		t.Errorf("expected DP bucket for depth 30 to have sites, got %d", got)
	}
	// With step=10 the idist label is the bucket index (i-1+min), matching
	// upstream: depth 30 → bucket index 4 → label "3".
	if lbl := res.dpLabel(res.dpIdx(30)); lbl != "3" {
		t.Errorf("dpLabel for depth 30 (step 10) = %q, want \"3\"", lbl)
	}
}

func TestStatsBCFStream(t *testing.T) {
	// Magic-byte sanity check: a stream that begins with "BCF\x02\x02" but
	// has no real BCF payload should produce a header read error, not panic.
	var buf bytes.Buffer
	bad := []byte{'B', 'C', 'F', 0x02, 0x02}
	_, err := Stats(bytes.NewReader(bad), &buf, StatsOptions{})
	if err == nil {
		t.Errorf("expected error reading malformed BCF stream")
	}
}

func TestStatsBCFRoundtrip(t *testing.T) {
	// Build a real BCF stream from a small VCF fixture using the bcf writer,
	// then re-read it through Stats and verify the SN counters match.
	mini := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	A	G	30	PASS	DP=10	GT	0/1
chr1	200	.	C	T	40	PASS	DP=20	GT	1/1
chr1	300	.	A	AT	50	PASS	DP=30	GT	0/1
`
	// Stats over text VCF first.
	resVCF, _ := runStats(t, mini, StatsOptions{})
	// Now write through the BCF writer.
	var bcfBuf bytes.Buffer
	if err := writeFixtureAsBCF(t, mini, &bcfBuf); err != nil {
		t.Fatalf("write BCF: %v", err)
	}
	var out bytes.Buffer
	resBCF, err := Stats(&bcfBuf, &out, StatsOptions{})
	if err != nil {
		t.Fatalf("Stats(BCF): %v", err)
	}
	if resBCF.numRecords != resVCF.numRecords {
		t.Errorf("BCF numRecords %d != VCF %d", resBCF.numRecords, resVCF.numRecords)
	}
	if resBCF.numSNPs != resVCF.numSNPs {
		t.Errorf("BCF numSNPs %d != VCF %d", resBCF.numSNPs, resVCF.numSNPs)
	}
	if resBCF.numIndels != resVCF.numIndels {
		t.Errorf("BCF numIndels %d != VCF %d", resBCF.numIndels, resVCF.numIndels)
	}
}

// writeTestFile is a tiny helper duplicating the snippet used in the view tests.
func writeTestFile(path, content string) error {
	return os.WriteFile(path, []byte(content), 0o644)
}

// writeFixtureAsBCF reads a VCF string, writes it through the bcf.Writer, and
// returns the BCF bytes in `out`. Used to exercise the BCF stats path with
// a known-good encoder.
func writeFixtureAsBCF(t *testing.T, vcfText string, out *bytes.Buffer) error {
	t.Helper()
	r := vcf.NewReader(strings.NewReader(vcfText))
	hdr, err := r.ReadHeader()
	if err != nil {
		return err
	}
	w, err := bcf.NewWriterFromVCFHeader(out, hdr)
	if err != nil {
		return err
	}
	for {
		v, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if err := w.Write(v); err != nil {
			return err
		}
	}
	return w.Flush()
}
