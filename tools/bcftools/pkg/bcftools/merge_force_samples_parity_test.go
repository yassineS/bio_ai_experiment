package bcftools

// Live upstream-binary parity for `bcftools merge --force-samples`.
//
// Upstream merge errors on duplicate sample names unless --force-samples is
// given, in which case the clashing name from input i (0-based) is prefixed
// with "<i+1>:" (repeatedly, if the prefixed form also clashes). This test
// bgzips+indexes the fixtures with the vendored htslib tools (upstream merge
// requires indexed input), runs the genuine upstream binary and our port on
// the same inputs, and asserts byte-equality of stdout modulo provenance.
// It t.Fatalf's (never t.Skip) when the binaries are unavailable.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// writeMergeInput writes content to dir/<name> then bgzips+tabix-indexes it via
// the shared bgzipIndex helper, returning the .gz path (upstream merge needs an
// indexed input).
func writeMergeInput(t *testing.T, htslibDir, dir, name, content string) string {
	t.Helper()
	plain := filepath.Join(dir, name)
	if err := os.WriteFile(plain, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return bgzipIndex(t, htslibDir, plain, "-p", "vcf", "-f")
}

const mergeForceA = `##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B
1	10	.	A	C	.	.	.	GT	0/1	1/1
`

const mergeForceB = `##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	C
1	10	.	A	C	.	.	.	GT	0/0	0/1
`

const mergeForceC = `##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	D
1	10	.	A	C	.	.	.	GT	1/1	0/0
`

// TestMerge_ForceSamples_UpstreamParity checks the two- and three-way clash
// renaming (A, 2:A, 3:A) against the upstream C binary byte-for-byte.
func TestMerge_ForceSamples_UpstreamParity(t *testing.T) {
	live, ours := requireLive(t)
	tools := upstreamBcftoolsConsensus(t)
	htslibDir := filepath.Dir(tools.bgzip)

	dir := t.TempDir()
	a := writeMergeInput(t, htslibDir, dir, "a.vcf", mergeForceA)
	b := writeMergeInput(t, htslibDir, dir, "b.vcf", mergeForceB)
	c := writeMergeInput(t, htslibDir, dir, "c.vcf", mergeForceC)

	cases := [][]string{
		{a, b},
		{a, b, c},
	}
	for i, inputs := range cases {
		args := append([]string{"merge", "--no-version", "--force-samples"}, inputs...)
		t.Run([]string{"two-way", "three-way"}[i], func(t *testing.T) {
			want := stripProvenanceBytes(runBin(t, live, args...))
			got := stripProvenanceBytes(runBin(t, ours, args...))
			if !bytes.Equal(want, got) {
				t.Fatalf("merge --force-samples mismatch\n--- upstream ---\n%s\n--- ours ---\n%s", want, got)
			}
		})
	}
}
