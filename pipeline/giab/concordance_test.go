package giab

import (
	"strings"
	"testing"
)

// vcfHeader is a minimal valid VCF header for building synthetic call sets.
const vcfHeader = `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
`

// rec builds a VCF data line. gtPL is the SAMPLE column (e.g. "0/1:0,30,255").
func rec(chrom string, pos int, ref, alt, qual, filter, format, sample string) string {
	return strings.Join([]string{chrom, itoa(pos), ".", ref, alt, qual, filter, ".", format, sample}, "\t")
}

func itoa(i int) string {
	// Tiny local helper to avoid importing strconv in many test files.
	if i == 0 {
		return "0"
	}
	neg := i < 0
	if neg {
		i = -i
	}
	var b [20]byte
	p := len(b)
	for i > 0 {
		p--
		b[p] = byte('0' + i%10)
		i /= 10
	}
	if neg {
		p--
		b[p] = '-'
	}
	return string(b[p:])
}

func vcf(lines ...string) string {
	return vcfHeader + strings.Join(lines, "\n") + "\n"
}

func mustParse(t *testing.T, s string) []VCFRecord {
	t.Helper()
	recs, err := ParseVCF(strings.NewReader(s))
	if err != nil {
		t.Fatalf("ParseVCF: %v", err)
	}
	return recs
}

func TestCompare_Identical(t *testing.T) {
	ours := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", "PASS", "GT:PL", "0/1:255,0,255")))
	up := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", "PASS", "GT:PL", "0/1:255,0,255")))
	c := CompareCallSets(ours, up, nil, DefaultQualULP)
	if c.Common != 1 || c.Identical != 1 || c.Differ != 0 {
		t.Fatalf("expected 1 identical, got %+v", c)
	}
	if c.GenotypeOrFilterFlips != 0 {
		t.Fatalf("identical records should not flip: %+v", c)
	}
}

func TestCompare_QualULPDoesNotFlip(t *testing.T) {
	// QUAL differs by 0.0001, GT and FILTER identical -> ULP-only, no flip.
	ours := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0000", "PASS", "GT:PL", "0/1:255,0,255")))
	up := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0001", "PASS", "GT:PL", "0/1:255,0,255")))
	c := CompareCallSets(ours, up, nil, DefaultQualULP)
	if c.Differ != 1 {
		t.Fatalf("expected 1 differing record, got %+v", c)
	}
	if c.QualULPOnly != 1 {
		t.Fatalf("expected the QUAL diff to be ULP-only, got %+v", c)
	}
	if c.GenotypeOrFilterFlips != 0 {
		t.Fatalf("ULP QUAL diff must not flip a genotype/FILTER, got %+v", c)
	}
	if len(c.Diffs) != 1 || c.Diffs[0].Flips() {
		t.Fatalf("diff should be recorded as non-flipping: %+v", c.Diffs)
	}
	if !c.Diffs[0].QualULP {
		t.Fatalf("diff should be marked QualULP: %+v", c.Diffs[0])
	}
}

func TestCompare_PLULPDoesNotFlip(t *testing.T) {
	// PL last element differs by 1 (a Phred ULP), same GT/FILTER.
	ours := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", "PASS", "GT:PL", "0/1:255,0,254")))
	up := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", "PASS", "GT:PL", "0/1:255,0,255")))
	c := CompareCallSets(ours, up, nil, DefaultQualULP)
	if c.Differ != 1 || c.QualULPOnly != 1 || c.GenotypeOrFilterFlips != 0 {
		t.Fatalf("PL ULP diff should be ULP-only non-flip, got %+v", c)
	}
}

func TestCompare_GenotypeFlip(t *testing.T) {
	// GT differs: 0/1 vs 1/1 -> a genuine flip.
	ours := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", "PASS", "GT:PL", "0/1:255,0,255")))
	up := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", "PASS", "GT:PL", "1/1:255,10,0")))
	c := CompareCallSets(ours, up, nil, DefaultQualULP)
	if c.Differ != 1 {
		t.Fatalf("expected 1 differing record, got %+v", c)
	}
	if c.GenotypeOrFilterFlips != 1 {
		t.Fatalf("expected a genotype flip, got %+v", c)
	}
	if c.QualULPOnly != 0 {
		t.Fatalf("a GT flip is not a ULP-only diff, got %+v", c)
	}
	if !c.Diffs[0].GTFlip || !c.Diffs[0].Flips() {
		t.Fatalf("diff should be marked as GT flip: %+v", c.Diffs[0])
	}
}

func TestCompare_FilterFlip(t *testing.T) {
	// FILTER differs: PASS vs LowQual -> a PASS/FAIL flip.
	ours := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", "PASS", "GT:PL", "0/1:255,0,255")))
	up := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", "LowQual", "GT:PL", "0/1:255,0,255")))
	c := CompareCallSets(ours, up, nil, DefaultQualULP)
	if c.GenotypeOrFilterFlips != 1 || !c.Diffs[0].FilterFlip {
		t.Fatalf("expected a FILTER flip, got %+v", c)
	}
}

func TestCompare_PhasingNotAFlip(t *testing.T) {
	// 0|1 (phased) vs 0/1 (unphased) is the same genotype -> identical.
	ours := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", "PASS", "GT:PL", "0|1:255,0,255")))
	up := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", "PASS", "GT:PL", "0/1:255,0,255")))
	c := CompareCallSets(ours, up, nil, DefaultQualULP)
	if c.Identical != 1 {
		t.Fatalf("phasing-only difference should be identical, got %+v", c)
	}
}

func TestCompare_DotFilterIsPass(t *testing.T) {
	// "." FILTER is treated as PASS by GIAB tooling; PASS vs "." should not flip.
	ours := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", ".", "GT:PL", "0/1:255,0,255")))
	up := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", "PASS", "GT:PL", "0/1:255,0,255")))
	c := CompareCallSets(ours, up, nil, DefaultQualULP)
	if c.Identical != 1 {
		t.Fatalf(". vs PASS should be identical, got %+v", c)
	}
}

func TestCompare_RegionRestriction(t *testing.T) {
	// Two sites: one inside the BED, one outside. Only the inside one counts.
	ours := mustParse(t, vcf(
		rec("chr1", 100, "A", "G", "200.0", "PASS", "GT", "0/1"),
		rec("chr1", 5000, "C", "T", "200.0", "PASS", "GT", "1/1"),
	))
	up := mustParse(t, vcf(
		rec("chr1", 100, "A", "G", "200.0", "PASS", "GT", "0/1"),
		rec("chr1", 5000, "C", "T", "200.0", "PASS", "GT", "0/1"), // differs, but outside
	))
	// BED covers [0,200) on chr1 -> includes POS 100 (0-based 99), excludes 5000.
	region, err := ParseBED(strings.NewReader("chr1\t0\t200\n"))
	if err != nil {
		t.Fatal(err)
	}
	c := CompareCallSets(ours, up, region, DefaultQualULP)
	if c.Common != 1 {
		t.Fatalf("only the in-region site should be common, got %+v", c)
	}
	if c.Differ != 0 || c.GenotypeOrFilterFlips != 0 {
		t.Fatalf("the differing site is outside the BED and must be ignored, got %+v", c)
	}
}

func TestCompare_OnlyOursOnlyUp(t *testing.T) {
	ours := mustParse(t, vcf(rec("chr1", 100, "A", "G", "200.0", "PASS", "GT", "0/1")))
	up := mustParse(t, vcf(rec("chr1", 200, "C", "T", "200.0", "PASS", "GT", "0/1")))
	c := CompareCallSets(ours, up, nil, DefaultQualULP)
	if c.OnlyOurs != 1 || c.OnlyUp != 1 || c.Common != 0 {
		t.Fatalf("expected disjoint call sets, got %+v", c)
	}
}

func TestCompare_QualBeyondULPNotULPButNoFlip(t *testing.T) {
	// QUAL differs by 50 (way beyond ULP) but GT/FILTER unchanged: differ,
	// not ULP-only, not a flip.
	ours := mustParse(t, vcf(rec("chr1", 100, "A", "G", "100.0", "PASS", "GT:PL", "0/1:255,0,255")))
	up := mustParse(t, vcf(rec("chr1", 100, "A", "G", "150.0", "PASS", "GT:PL", "0/1:255,0,255")))
	c := CompareCallSets(ours, up, nil, DefaultQualULP)
	if c.Differ != 1 || c.QualULPOnly != 0 || c.GenotypeOrFilterFlips != 0 {
		t.Fatalf("large QUAL gap w/ same GT: differ but not ULP and not flip, got %+v", c)
	}
}

func TestHeadline(t *testing.T) {
	c := Concordance{Common: 10, Identical: 8, Differ: 2, QualULPOnly: 2, GenotypeOrFilterFlips: 0}
	h := c.Headline()
	if !strings.Contains(h, "0 of which flip") {
		t.Fatalf("headline should report zero flips: %q", h)
	}
}
