package edgecases

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// sortInputSAM mixes:
//   - records that tie on (RNAME, POS) but differ in QNAME (stable-order check),
//   - QNAMEs with numeric suffixes that must sort by strnum_cmp natural order
//     (x2 < x9 < x10 < x100, NOT lexical),
//   - unmapped reads (FLAG 4, RNAME '*') that must be placed consistently.
const sortInputSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:1000
x10	0	chr1	100	60	5M	*	0	0	ACGTA	FFFFF
x2	0	chr1	100	60	5M	*	0	0	ACGTA	FFFFF
x100	0	chr1	100	60	5M	*	0	0	ACGTA	FFFFF
x9	0	chr1	50	60	5M	*	0	0	ACGTA	FFFFF
read1a	0	chr2	10	60	5M	*	0	0	ACGTA	FFFFF
read1	0	chr2	10	60	5M	*	0	0	ACGTA	FFFFF
um2	4	*	0	0	*	*	0	0	ACGTA	FFFFF
um1	4	*	0	0	*	*	0	0	ACGTA	FFFFF
`

// TestSortStabilityStrnumCmp verifies coordinate sort and queryname (-n) sort
// at -@1 are byte-for-byte identical to upstream, including tie-break order and
// unmapped/'*' placement. Multi-thread tie order is allowed to differ between
// implementations, so the test pins -@1 (single thread) where order is
// deterministic.
func TestSortStabilityStrnumCmp(t *testing.T) {
	our := ourBin(t, "samtools")
	up := upBin(t, "samtools")

	cases := []struct {
		name string
		args []string
	}{
		{"coordinate", []string{"sort", "-@1", "-O", "sam"}},
		{"queryname_strnum", []string{"sort", "-n", "-@1", "-O", "sam"}},
	}

	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			dir := t.TempDir()
			src := writeFile(t, dir, "in.sam", sortInputSAM)
			ours := mustRun(t, our, append(append([]string{}, tc.args...), src)...)
			ups := mustRun(t, up, append(append([]string{}, tc.args...), src)...)
			// Compare with @PG dropped / @SQ provenance stripped, but WITHOUT
			// sorting records — the whole point is to verify record order.
			g := upstream.NormalizeSAM(ours, false)
			w := upstream.NormalizeSAM(ups, false)
			if g != w {
				t.Errorf("%s sort at -@1 diverges from upstream (record order / tie-break / unmapped placement):\n--- ours ---\n%s\n--- upstream ---\n%s",
					tc.name, g, w)
			}
		})
	}
}
