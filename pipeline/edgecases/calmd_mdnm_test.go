package edgecases

import (
	"sort"
	"strings"
	"testing"
)

// calmdRefFASTA has a run of N (ambiguity) bases in the middle so reads aligned
// there exercise the ambiguous-reference path of MD/NM computation.
const calmdRefFASTA = `>chr1
ACGTACGTACGTACGTACGTNNNNNACGTACGTACGTACGTACGTACGTACGT
`

// calmdInputSAM exercises plain matches, a substitution, explicit =/X CIGAR
// ops, an N (reference skip) op, and a read over the N-run — every path where
// MD/NM can be miscomputed.
const calmdInputSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:53
plain	0	chr1	1	60	8M	*	0	0	ACGTACGT	FFFFFFFF
mismatch	0	chr1	1	60	8M	*	0	0	ACCTACGT	FFFFFFFF
eqx	0	chr1	1	60	4=1X3=	*	0	0	ACGTACGT	FFFFFFFF
nskip	0	chr1	1	60	4M3N4M	*	0	0	ACGTACGT	FFFFFFFF
ambig	0	chr1	21	60	5M	*	0	0	ACGTA	FFFFF
`

// mdnmTags extracts a sorted "QNAME MD NM" line per record so the comparison is
// on the recomputed tags, independent of column order or unrelated aux tags.
func mdnmTags(sam string) []string {
	var out []string
	for _, ln := range strings.Split(sam, "\n") {
		if ln == "" || strings.HasPrefix(ln, "@") {
			continue
		}
		cols := strings.Split(ln, "\t")
		if len(cols) < 11 {
			continue
		}
		var md, nm string
		for _, c := range cols[11:] {
			switch {
			case strings.HasPrefix(c, "MD:Z:"):
				md = c
			case strings.HasPrefix(c, "NM:i:"):
				nm = c
			}
		}
		out = append(out, cols[0]+" "+md+" "+nm)
	}
	sort.Strings(out)
	return out
}

// TestCalmdMDNMTags verifies `samtools calmd` recomputes MD:Z and NM:i exactly
// as upstream across =/X CIGAR ops, N (reference-skip) ops, and ambiguity
// codes (N bases in the reference). A wrong MD/NM is silent corruption of
// downstream variant calling that relies on these tags.
func TestCalmdMDNMTags(t *testing.T) {
	our := ourBin(t, "samtools")
	up := upBin(t, "samtools")

	dir := t.TempDir()
	ref := writeFile(t, dir, "ref.fa", calmdRefFASTA)
	src := writeFile(t, dir, "in.sam", calmdInputSAM)

	ours := mustRun(t, our, "calmd", src, ref)
	// Upstream calmd auto-detects SAM but accept -S explicitly for clarity.
	ups := mustRun(t, up, "calmd", "-S", src, ref)

	g := strings.Join(mdnmTags(ours), "\n")
	w := strings.Join(mdnmTags(ups), "\n")
	if g != w {
		t.Errorf("calmd MD/NM diverges from upstream:\n--- ours ---\n%s\n--- upstream ---\n%s", g, w)
	}
}
