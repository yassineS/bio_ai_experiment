package bcftools

// Live upstream-binary parity for two fixes:
//
//  1. `bcftools view -c/-q` (and -C/-Q/-x/-X) recompute and append INFO/AC and
//     INFO/AN even without a sample subset, matching upstream vcfview.c's
//     calc_ac path. The -I/--no-update variant suppresses it.
//
//  2. `bcftools norm -m+` joins biallelic records at the same position into a
//     multiallelic using upstream's exact allele-merge and FORMAT/GT merge
//     rules (vcfnorm.c). In particular "0/1" + "0/1" joins to "2/1" (the second
//     record's allele lands in the first free strand), and indels with
//     differing REF lengths share a common, padded REF.
//
// Both run the genuine upstream C binary (located/built via requireLive) and
// our freshly built port on the same input, then compare stdout byte-for-byte
// after stripping the provenance header lines. They t.Fatalf (never t.Skip)
// when the upstream binary or our port binary is unavailable.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// writeTmpVCF writes content to a temp .vcf file and returns its path.
func writeTmpVCF(t *testing.T, content string) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write temp vcf: %v", err)
	}
	return p
}

// oracleEqualVCF runs both binaries with the given args (the input path is
// appended) and asserts byte-equality of stdout modulo provenance lines.
func oracleEqualVCF(t *testing.T, live, ours, path string, args ...string) {
	t.Helper()
	upArgs := append([]string{}, args...)
	upArgs = append(upArgs, path)
	want := stripProvenanceBytes(runBin(t, live, upArgs...))
	got := stripProvenanceBytes(runBin(t, ours, upArgs...))
	if !bytes.Equal(want, got) {
		t.Fatalf("byte mismatch for args %v\n--- upstream ---\n%s\n--- ours ---\n%s", args, want, got)
	}
}

const acanParityVCF = `##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
1	10	.	A	C	.	.	.	GT	0/1	0/0	1/1
1	20	.	A	C,G	.	.	.	GT	0/1	1/2	2/2
1	30	.	G	T	.	.	.	GT	0/0	0/0	0/1
`

// acanParityWithInfoVCF carries pre-existing INFO/AC and INFO/AN so the test
// covers upstream's preference for the existing tags (bcf_calc_ac reads
// INFO/AC,AN first when no sample subset is applied).
const acanParityWithInfoVCF = `##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##INFO=<ID=AC,Number=A,Type=Integer,Description="orig">
##INFO=<ID=AN,Number=1,Type=Integer,Description="orig">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	10	.	A	C	.	.	AC=99;AN=99	GT	0/1	1/1
1	20	.	A	C,G	.	.	AC=2,3;AN=6	GT	0/1	1/2
`

// TestView_ACANRecompute_UpstreamParity verifies that `view -c/-q/-C/-x` add
// recomputed INFO/AC and INFO/AN to the output (and -I suppresses it), exactly
// matching the upstream C binary byte-for-byte.
func TestView_ACANRecompute_UpstreamParity(t *testing.T) {
	live, ours := requireLive(t)
	path := writeTmpVCF(t, acanParityVCF)
	infoPath := writeTmpVCF(t, acanParityWithInfoVCF)

	cases := [][]string{
		{"view", "--no-version", "-c", "1"},
		{"view", "--no-version", "-C", "4"},
		{"view", "--no-version", "-q", "0.3"},
		{"view", "--no-version", "-Q", "0.6"},
		{"view", "--no-version", "-c", "1", "-I"}, // -I suppresses the update
		{"view", "--no-version", "-x", "-s", "S1,S2"},
		{"view", "--no-version", "-s", "S1,S2"}, // subset triggers update
		{"view", "--no-version", "-s", "S1,S2", "-I"},
	}
	for _, args := range cases {
		t.Run(args[2], func(t *testing.T) {
			oracleEqualVCF(t, live, ours, path, args...)
		})
	}

	// Pre-existing INFO/AC,AN must be preserved (not recomputed) without a
	// subset, matching bcf_calc_ac's INFO-first behaviour.
	t.Run("preexisting-info", func(t *testing.T) {
		oracleEqualVCF(t, live, ours, infoPath, "view", "--no-version", "-c", "1")
	})
}

const normJoinParityVCF = `##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
1	100	.	A	C	.	.	.	GT	0/1
1	100	.	A	G	.	.	.	GT	0/1
1	100	.	AT	A	.	.	.	GT	0/1
1	200	.	C	T	.	.	.	GT	1/1
1	200	.	CG	C	.	.	.	GT	0/1
1	300	.	G	GA	.	.	.	GT	0/1
1	300	.	G	GAA	.	.	.	GT	0/1
`

// normJoinDiffRefVCF exercises joining indels with differing REF lengths into
// a common padded REF (AT>A + ATG>A -> ATG>{AG,A}).
const normJoinDiffRefVCF = `##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
1	100	.	AT	A	.	.	.	GT	0/1
1	100	.	ATG	A	.	.	.	GT	0/1
`

// normJoinPhasedVCF exercises phased multi-ploidy join (three biallelic SNPs
// joined into one record with phased GT preserved).
const normJoinPhasedVCF = `##fileformat=VCFv4.2
##contig=<ID=20,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
20	1	.	A	C	.	.	.	GT	0|0|1|0
20	1	.	A	G	.	.	.	GT	0|1|0|0
20	1	.	A	T	.	.	.	GT	1|0|0|0
`

// TestNorm_JoinMultiallelic_UpstreamParity verifies `norm -m+` join byte-for-byte
// against the upstream C binary across snps/indels/both/any modes, including
// the GT-ordering fix and the common-REF padding for differing-length indels.
func TestNorm_JoinMultiallelic_UpstreamParity(t *testing.T) {
	live, ours := requireLive(t)

	fixtures := map[string]string{
		"mixed":   normJoinParityVCF,
		"diffref": normJoinDiffRefVCF,
		"phased":  normJoinPhasedVCF,
	}
	// -m +both spelled with a space is also valid upstream; the last entry
	// exercises that separate-argument form.
	modeArgs := [][]string{
		{"norm", "--no-version", "-N", "-m+"},
		{"norm", "--no-version", "-N", "-m+snps"},
		{"norm", "--no-version", "-N", "-m+indels"},
		{"norm", "--no-version", "-N", "-m+any"},
		{"norm", "--no-version", "-N", "-m+both"},
		{"norm", "--no-version", "-N", "-m", "+both"},
	}
	for name, content := range fixtures {
		path := writeTmpVCF(t, content)
		for _, args := range modeArgs {
			label := name + "_" + args[len(args)-1]
			t.Run(label, func(t *testing.T) {
				oracleEqualVCF(t, live, ours, path, args...)
			})
		}
	}
}
