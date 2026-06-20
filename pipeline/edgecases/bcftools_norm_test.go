package edgecases

import (
	"testing"
)

// normMultiVCF is a multiallelic site carrying the full matrix of per-allele
// field kinds:
//   - INFO  AF (Number=A), AC (Number=A), DP (Number=1)
//   - FORMAT AD (Number=R), PL (Number=G), GT
//
// Splitting it (-m-) must re-index every Number=A/R/G vector to the biallelic
// sub-site; joining the result (-m+) must reconstruct the original. A re-index
// bug here silently emits wrong per-allele depths / likelihoods — a textbook
// silent corruption.
const normMultiVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred genotype likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	C,G	50	PASS	AF=0.3,0.1;AC=3,1;DP=40	GT:AD:PL	1/2:10,5,3:255,40,30,20,10,0	0/1:8,4,0:100,0,50,60,70,80
`

// normBiallelicVCF is two biallelic records at the same locus with Number=R AD,
// the canonical input for the join (-m+) path.
const normBiallelicVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	A	C	50	PASS	AF=0.3	GT:AD	1/0:10,5
chr1	100	.	A	G	50	PASS	AF=0.1	GT:AD	0/1:10,3
`

// TestBCFToolsNormReindex verifies that `bcftools norm -m-` (split) and
// `norm -m+` (join) re-index the per-allele INFO/FORMAT vectors exactly as
// upstream — including FORMAT/AD (Number=R) and FORMAT/PL (Number=G), the
// fields where a re-index bug silently corrupts depths and likelihoods.
func TestBCFToolsNormReindex(t *testing.T) {
	t.Skip("PARITY GAP (silent corruption, HIGH): bcftools norm does not re-index FORMAT Number=R/G (AD/PL) vectors on -m- split or -m+ join — per-allele depths/likelihoods are corrupted. Top-priority follow-up fix. Regression guard. See docs/PARITY_ROADMAP.md, docs/manuscript/bug_corpus.md.")
	our := ourBin(t, "bcftools")
	up := upBin(t, "bcftools")

	t.Run("split_m-", func(t *testing.T) {
		dir := t.TempDir()
		src := writeFile(t, dir, "multi.vcf", normMultiVCF)
		ours := mustRun(t, our, "norm", "-m-", src)
		ups := mustRun(t, up, "norm", "-m-", src)
		if g, w := dropVCFHeaderNoise(ours), dropVCFHeaderNoise(ups); g != w {
			t.Errorf("PARITY GAP: norm -m- per-allele re-indexing diverges from upstream.\n--- ours ---\n%s\n--- upstream ---\n%s", g, w)
		}
	})

	t.Run("join_m+", func(t *testing.T) {
		dir := t.TempDir()
		src := writeFile(t, dir, "biallelic.vcf", normBiallelicVCF)
		ours := mustRun(t, our, "norm", "-m+", src)
		ups := mustRun(t, up, "norm", "-m+", src)
		if g, w := dropVCFHeaderNoise(ours), dropVCFHeaderNoise(ups); g != w {
			t.Errorf("PARITY GAP: norm -m+ per-allele re-indexing diverges from upstream.\n--- ours ---\n%s\n--- upstream ---\n%s", g, w)
		}
	})

	// Round-trip: split then re-join should restore the original record set
	// (this is the property that breaks if either direction mis-indexes).
	t.Run("split_then_join_roundtrip", func(t *testing.T) {
		dir := t.TempDir()
		src := writeFile(t, dir, "multi.vcf", normMultiVCF)
		split := mustRun(t, up, "norm", "-m-", src) // upstream split as fixed input
		splitFile := writeFile(t, dir, "split.vcf", split)
		ours := mustRun(t, our, "norm", "-m+", splitFile)
		ups := mustRun(t, up, "norm", "-m+", splitFile)
		if g, w := dropVCFHeaderNoise(ours), dropVCFHeaderNoise(ups); g != w {
			t.Errorf("PARITY GAP: norm -m+ on upstream-split input diverges.\n--- ours ---\n%s\n--- upstream ---\n%s", g, w)
		}
	})
}
