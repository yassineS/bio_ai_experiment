package bcftools

import (
	"os"
	"path/filepath"
	"testing"
)

// Live-oracle parity tests for the native setGT plugin, exercising the binomial
// (-t b), random (-t r with -s), and read-depth (-n X) modes that previously
// reported "not supported". Each case is driven through BOTH the real upstream
// bcftools binary (1.23.x, via assertPluginParity's CLI-to-CLI oracle) and our
// port and is asserted byte-identical (modulo provenance). The random cases use
// a fixed -s seed so the deterministic drand48 stream matches upstream exactly.

// TestNativePluginSetGTReadDepth covers -n X (set every allele to the allele
// with the largest FORMAT/AD) across the target selectors.
func TestNativePluginSetGTReadDepth(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	cases := [][]string{
		{"-t", "a", "-n", "X"},
		{"-t", ".", "-n", "X"},
		{"-t", "./.", "-n", "X"},
		{"-t", "./x", "-n", "X"},
		{"-t", "a", "-n", "c:0/X"},
		{"-t", "a", "-n", "c:X/X"},
		{"-t", "a", "-n", "c:X|X"},
	}
	for _, args := range cases {
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, gt, "setGT", args...)
		})
	}
}

// TestNativePluginSetGTBinom covers the two-tailed binomial selector over
// FORMAT/AD, including every comparison operator and several new-gt targets.
func TestNativePluginSetGTBinom(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	skew := writeSetGTSkewFixture(t)
	cases := []struct {
		fixture string
		args    []string
	}{
		{gt, []string{"-t", "b:AD<1e-3", "-n", "0"}},
		{gt, []string{"-t", "b:AD<0.8", "-n", "0"}},
		{gt, []string{"-t", "b:AD>0.1", "-n", "."}},
		{gt, []string{"-t", "b:AD>=0.5", "-n", "c:m/M"}},
		{gt, []string{"-t", "b:AD<=0.9", "-n", "X"}},
		{gt, []string{"-t", "b:AD==1", "-n", "."}},
		// The skewed fixture has hets that actually fail the test, so these
		// thresholds exercise the genotype-changing path, not just no-ops.
		{skew, []string{"-t", "b:AD<0.001", "-n", "0"}},
		{skew, []string{"-t", "b:AD<0.01", "-n", "0"}},
		{skew, []string{"-t", "b:AD<0.05", "-n", "0"}},
		{skew, []string{"-t", "b:AD<0.1", "-n", "c:0/0"}},
		{skew, []string{"-t", "b:AD<0.5", "-n", "."}},
		{skew, []string{"-t", "b:AD<0.9", "-n", "0"}},
	}
	for _, tc := range cases {
		t.Run(filepath.Base(tc.fixture)+" "+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "setGT", tc.args...)
		})
	}
}

// TestNativePluginSetGTRandom covers the random selector. Because setGT seeds
// htslib's deterministic drand48 from -s (default 0) and nothing else, a fixed
// seed yields a byte-reproducible result that must match upstream exactly.
func TestNativePluginSetGTRandom(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	cases := [][]string{
		{"-t", "r:0.5", "-n", "."}, // default seed 0
		{"-t", "r:0.5", "-s", "1", "-n", "."},
		{"-t", "r:0.3", "-s", "7", "-n", "0"},
		{"-t", "r:0.7", "-s", "99", "-n", "c:m/m"},
		{"-t", "r:0.5", "-s", "42", "-n", "X"},
		{"-t", "r:0.25", "-s", "12345", "-n", "0p"},
	}
	for _, args := range cases {
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, gt, "setGT", args...)
		})
	}
}

// writeSetGTSkewFixture writes a small VCF whose heterozygous genotypes have
// strongly skewed FORMAT/AD, so the binomial test changes genotypes rather than
// passing them through unchanged. It returns the absolute path.
func writeSetGTSkewFixture(t *testing.T) string {
	t.Helper()
	const body = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	.	A	T	.	.	.	GT:AD	0/1:1,19	0/1:10,10	0/1:20,1
chr1	150	.	A	T	.	.	.	GT:AD	1|0:3,17	0/0:18,0	./.:0,0
chr1	200	.	A	T,C	.	.	.	GT:AD	1/2:2,18,0	0|1:30,2,0	1/2:0,1,25
`
	dir := t.TempDir()
	path := filepath.Join(dir, "setgt_skew.vcf")
	if err := os.WriteFile(path, []byte(body), 0o644); err != nil {
		t.Fatalf("write skew fixture: %v", err)
	}
	return path
}
