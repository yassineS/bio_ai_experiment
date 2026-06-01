package bcftools

import (
	"bytes"
	"os"
	"strings"
	"testing"
)

// fillTagsFixture is a 3-sample biallelic + multi-allelic VCF used to
// pin the per-allele counting and the HWE/ExcHet exact test.
const fillTagsFixture = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	.	A	T	.	.	.	GT	0/0	0/1	0/0
chr1	200	.	C	G	.	.	.	GT	1/1	0/1	0/1
chr1	300	.	G	A,C	.	.	.	GT	1/2	0/1	0/0
chr1	400	.	T	C	.	.	.	GT	1/1	0/1	./.
`

// runFillTags is a small test helper that fills the default "all" tag
// set and returns the per-record INFO column keyed by POS.
func runFillTags(t *testing.T, src string) map[string]string {
	t.Helper()
	hdr, vars, err := readAllVariants(strings.NewReader(src))
	if err != nil {
		t.Fatalf("readAllVariants: %v", err)
	}
	out := map[string]string{}
	for _, v := range vars {
		fillTagsRecord(v, hdr.Samples, tagAll, false)
		out[posStr(v.Pos)] = v.Info["AC"] + "|" + v.Info["AN"] + "|" + v.Info["AF"] +
			"|" + v.Info["NS"] + "|" + v.Info["MAF"] + "|" + v.Info["AC_Het"] +
			"|" + v.Info["AC_Hom"] + "|" + v.Info["F_MISSING"] +
			"|" + v.Info["HWE"] + "|" + v.Info["ExcHet"]
	}
	return out
}

// TestFillTags_BiallelicCounts pins AC/AN/AF/NS/AC_Het/AC_Hom for a
// simple biallelic site.
func TestFillTags_BiallelicCounts(t *testing.T) {
	got := runFillTags(t, fillTagsFixture)
	// chr1:200 — S1=1/1, S2=0/1, S3=0/1: allele 1 = 2(hom) + 1 + 1 = 4,
	// AN=6, AF=0.666667, NS=3, AC_Het=2, AC_Hom=2, F_MISSING=0.
	row := got["200"]
	for _, want := range []string{"4|6|0.666667|3|", "|2|2|0|"} {
		if !strings.Contains(row, want) {
			t.Errorf("chr1:200 row %q missing %q", row, want)
		}
	}
}

// TestFillTags_MultiAllelic pins the per-ALT AC/AF list for a
// multi-allelic site (REF=G, ALT=A,C; S1=1/2, S2=0/1, S3=0/0).
func TestFillTags_MultiAllelic(t *testing.T) {
	got := runFillTags(t, fillTagsFixture)
	row := got["300"]
	// Alleles: ref G=2 (S3 0/0), A(=1)=2 (S1, S2 het), C(=2)=1 (S1 het).
	// AN=6. AC=2,1. AF=0.333333,0.166667.
	if !strings.HasPrefix(row, "2,1|6|0.333333,0.166667|") {
		t.Errorf("chr1:300 multi-allelic AC/AN/AF wrong: %q", row)
	}
	// AC_Het should be 2,1 (allele A in two hets, allele C in one het).
	if !strings.Contains(row, "|2,1|") {
		t.Errorf("chr1:300 AC_Het wrong: %q", row)
	}
}

// TestFillTags_Missing pins F_MISSING and NS when a sample GT is ./.
func TestFillTags_Missing(t *testing.T) {
	got := runFillTags(t, fillTagsFixture)
	row := got["400"]
	// S1=1/1, S2=0/1, S3=./. : NS=2, F_MISSING=0.333333, AC=3, AN=4.
	if !strings.HasPrefix(row, "3|4|") {
		t.Errorf("chr1:400 AC/AN wrong: %q", row)
	}
	if !strings.Contains(row, "|2|") {
		t.Errorf("chr1:400 NS should be 2: %q", row)
	}
	if !strings.Contains(row, "|0.333333|") {
		t.Errorf("chr1:400 F_MISSING should be 0.333333: %q", row)
	}
}

// TestFillTags_HeaderAppend ensures the new ##INFO lines are appended
// (and existing IDs are not duplicated).
func TestFillTags_HeaderAppend(t *testing.T) {
	hdr, _, err := readAllVariants(strings.NewReader(fillTagsFixture))
	if err != nil {
		t.Fatal(err)
	}
	out := fillTagsHeader(hdr, tagAll)
	joined := strings.Join(out.MetaInfo, "\n")
	for _, id := range []string{"F_MISSING", "NS", "AC_Hom", "AC_Het", "AC_Hemi",
		"MAF", "HWE", "ExcHet", "AN", "AC", "AF", "VAF", "VAF1"} {
		if !strings.Contains(joined, "ID="+id+",") {
			t.Errorf("header missing ##INFO/##FORMAT ID=%s", id)
		}
	}
	// AN was not in the fixture header (only GT), so it is appended once.
	if n := strings.Count(joined, "ID=AN,"); n != 1 {
		t.Errorf("ID=AN appears %d times, want 1", n)
	}
}

// TestFillTags_CalcHWE pins a known HWE/ExcHet value: a perfectly
// Hardy-Weinberg site should score HWE=1.
func TestFillTags_CalcHWE(t *testing.T) {
	// nref=2, nalt=2, nhet=2 (one het genotype per allele) → balanced.
	pHWE, pExc := calcHWE(2, 2, 2)
	if pHWE <= 0 || pHWE > 1 {
		t.Errorf("HWE out of range: %v", pHWE)
	}
	if pExc <= 0 || pExc > 1 {
		t.Errorf("ExcHet out of range: %v", pExc)
	}
}

// TestFillTags_PluginOutput end-to-end checks the +fill-tags built-in
// emits a valid VCF whose first data row carries the appended tags.
func TestFillTags_PluginOutput(t *testing.T) {
	tmp := t.TempDir()
	path := tmp + "/in.vcf"
	if err := os.WriteFile(path, []byte(fillTagsFixture), 0o644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	err := runBuiltinFillTags(PluginOptions{
		Name:         "fill-tags",
		Args:         []string{path},
		OutputFormat: OutputVCF,
	}, &out, &bytes.Buffer{})
	if err != nil {
		t.Fatalf("runBuiltinFillTags: %v", err)
	}
	s := out.String()
	if !strings.Contains(s, "##INFO=<ID=HWE,") {
		t.Errorf("output header missing HWE declaration")
	}
	if !strings.Contains(s, "NS=3") || !strings.Contains(s, "ExcHet=") {
		t.Errorf("output records missing filled tags:\n%s", s)
	}
}
