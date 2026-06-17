package bcftools

// Live upstream-binary parity for `bcftools annotate -x` field removal.
//
// Covers the FILTER header-ordering fix (bare FILTER drops every ##FILTER line
// including PASS; FILTER/NAME drops only that line and rewrites a now-empty
// record to PASS), plus INFO and FORMAT removal (FORMAT keeps GT). Each case
// runs the genuine upstream C binary and our freshly built port on the same
// fixture and asserts byte-equality of stdout modulo provenance lines.
// t.Fatalf (never t.Skip) when a binary is unavailable.

import (
	"testing"
)

const annotateRemoveVCF = `##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele freq">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
1	10	rs1	A	C	5	q10	DP=10;AF=0.5	GT:DP	0/1:7
1	20	.	A	G	60	PASS	DP=20;AF=0.1	GT:DP	1/1:9
`

// TestAnnotate_Remove_UpstreamParity asserts every -x removal mode matches the
// upstream C binary byte-for-byte, with particular attention to FILTER header
// line ordering (the originally reported divergence).
func TestAnnotate_Remove_UpstreamParity(t *testing.T) {
	live, ours := requireLive(t)
	path := writeTmpVCF(t, annotateRemoveVCF)

	specs := []string{
		"FILTER",
		"FILTER/q10",
		"FILTER/s50",
		"INFO",
		"INFO/DP",
		"INFO/AF",
		"ID",
		"FORMAT",
		"FMT",
		"FORMAT/DP",
		"FORMAT/GT",
		"FILTER,INFO/DP",
		"FILTER/q10,INFO,ID",
	}
	for _, spec := range specs {
		t.Run(spec, func(t *testing.T) {
			oracleEqualVCF(t, live, ours, path, "annotate", "--no-version", "-x", spec)
		})
	}
}
