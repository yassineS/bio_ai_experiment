package bcftools

// Live upstream-binary parity for `bcftools concat -a` overlap-merge ordering.
//
// htslib's synced reader (bcf_sr_sort.c) emits records that share a position
// grouped by their REF>ALT variant string, with the groups ordered by
// descending pre-dedup record count and ties broken by first-appearance. Our
// port reproduces that in orderPositionGroup. These tests bgzip+index the
// fixtures (concat -a requires indexed input) and assert byte-equality against
// the genuine upstream binary, including the -d (rm-dup) modes whose collapse
// keeps the same emission order. t.Fatalf (never t.Skip) when a binary is
// missing.

import (
	"bytes"
	"path/filepath"
	"testing"
)

// concatOverlapA / B overlap at POS 20 with the same variants in different
// per-file order, exercising the count- and first-appearance-based ordering.
const concatOverlapA = `##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
1	20	.	A	G	.	.	.	GT	0/1
1	20	.	A	ATTT	.	.	.	GT	0/1
`

const concatOverlapB = `##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
1	20	.	C	T	.	.	.	GT	0/1
1	20	.	A	ATTT	.	.	.	GT	1/1
1	20	.	A	G	.	.	.	GT	1/1
`

// TestConcat_OverlapOrder_UpstreamParity asserts the overlap-merge ordering
// matches the upstream C binary byte-for-byte, in both input orders and under
// the exact/none rm-dup modes (which preserve the same ordering).
func TestConcat_OverlapOrder_UpstreamParity(t *testing.T) {
	live, ours := requireLive(t)
	tools := upstreamBcftoolsConsensus(t)
	htslibDir := filepath.Dir(tools.bgzip)

	dir := t.TempDir()
	a := writeMergeInput(t, htslibDir, dir, "ca.vcf", concatOverlapA)
	b := writeMergeInput(t, htslibDir, dir, "cb.vcf", concatOverlapB)

	orders := [][2]string{{a, b}, {b, a}}
	flagSets := [][]string{
		{"-a"},
		{"-a", "-d", "exact"},
		{"-a", "-d", "none"},
	}
	for _, ab := range orders {
		for _, flags := range flagSets {
			args := append([]string{"concat", "--no-version"}, flags...)
			args = append(args, ab[0], ab[1])
			name := filepath.Base(ab[0]) + "_" + filepath.Base(ab[1])
			for _, f := range flags {
				name += "_" + f
			}
			t.Run(name, func(t *testing.T) {
				want := stripProvenanceBytes(runBin(t, live, args...))
				got := stripProvenanceBytes(runBin(t, ours, args...))
				if !bytes.Equal(want, got) {
					t.Fatalf("concat overlap order mismatch\n--- upstream ---\n%s\n--- ours ---\n%s", want, got)
				}
			})
		}
	}
}
