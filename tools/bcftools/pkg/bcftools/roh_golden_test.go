package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// rohGoldenDir is where the upstream bcftools `roh` fixtures are
// vendored: roh.1.vcf.gz, roh.1.tab.gz and the roh.1.*.out goldens.
const rohGoldenDir = "../../testdata/roh"

// stripComments drops the leading "#"-prefixed header lines, matching
// the `| grep -v ^#` in bcftools' test.pl test_roh helper.
func stripComments(s string) string {
	var b strings.Builder
	for _, line := range strings.Split(s, "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		b.WriteString(line)
		b.WriteByte('\n')
	}
	return b.String()
}

// TestRoh_UpstreamGoldens replays the four `roh` invocations from
// reference_code/bcftools/test/test.pl and compares our output
// byte-for-byte (after dropping the comment header, exactly as
// test.pl does) against the upstream roh.1.*.out goldens.
//
// test.pl invocations (lines 1094-1100):
//
//	roh.1.1.out: -Or -G30 --AF-dflt 0.4
//	roh.1.1.out: -Or -G30 --AF-file roh.1.tab.gz
//	roh.1.1.out: -Or -G30 --AF-file roh.1.tab.gz --ignore-homref
//	roh.1.2.out: -G30 --AF-dflt 0.4 -r 1:100174876-100318245
//	roh.1.3.out: -G30 --AF-dflt 0.4 -r 1:100174876-100318245 --ignore-homref
//	roh.1.3.out: -G30 --AF-dflt 0.4 -r 1:... --ignore-homref --include-noalt
//	roh.1.4.out: -G30 --AF-dflt 0.4 -r 1:100174876-100318245 --include-noalt
func TestRoh_UpstreamGoldens(t *testing.T) {
	vcf := filepath.Join(rohGoldenDir, "roh.1.vcf.gz")
	tab := filepath.Join(rohGoldenDir, "roh.1.tab.gz")
	region := "1:100174876-100318245"

	cases := []struct {
		name   string
		golden string
		opts   RohOptions
	}{
		{
			name:   "Or_AFdflt",
			golden: "roh.1.1.out",
			opts:   RohOptions{GTsOnly: 30, AFDflt: AFDfltPtr(0.4), OutputTypes: "r"},
		},
		{
			name:   "Or_AFfile",
			golden: "roh.1.1.out",
			opts:   RohOptions{GTsOnly: 30, AFFile: tab, OutputTypes: "r"},
		},
		{
			name:   "Or_AFfile_ignoreHomref",
			golden: "roh.1.1.out",
			opts:   RohOptions{GTsOnly: 30, AFFile: tab, IgnoreHomRef: true, OutputTypes: "r"},
		},
		{
			name:   "region_AFdflt",
			golden: "roh.1.2.out",
			opts:   RohOptions{GTsOnly: 30, AFDflt: AFDfltPtr(0.4), Regions: []string{region}},
		},
		{
			name:   "region_AFdflt_ignoreHomref",
			golden: "roh.1.3.out",
			opts:   RohOptions{GTsOnly: 30, AFDflt: AFDfltPtr(0.4), Regions: []string{region}, IgnoreHomRef: true},
		},
		{
			name:   "region_AFdflt_ignoreHomref_noalt",
			golden: "roh.1.3.out",
			opts:   RohOptions{GTsOnly: 30, AFDflt: AFDfltPtr(0.4), Regions: []string{region}, IgnoreHomRef: true, IncludeNoalt: true},
		},
		{
			name:   "region_AFdflt_noalt",
			golden: "roh.1.4.out",
			opts:   RohOptions{GTsOnly: 30, AFDflt: AFDfltPtr(0.4), Regions: []string{region}, IncludeNoalt: true},
		},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			want, err := os.ReadFile(filepath.Join(rohGoldenDir, c.golden))
			if err != nil {
				t.Fatalf("read golden: %v", err)
			}
			var out bytes.Buffer
			if _, err := RohFile(vcf, &out, c.opts); err != nil {
				t.Fatalf("RohFile: %v", err)
			}
			got := stripComments(out.String())
			if got != string(want) {
				t.Errorf("output mismatch for %s\n--- got ---\n%s\n--- want ---\n%s",
					c.golden, got, want)
			}
		})
	}
}

// TestRoh_DistanceScaledTransitions confirms that the transition
// probability for a gap between markers scales with the physical
// distance: the precomputed matrix-power chain means a 2-bp gap is
// not the same as a 1-bp gap.
func TestRoh_DistanceScaledTransitions(t *testing.T) {
	h := newHMM(2, baseTprob(DefaultHWtoAZ, DefaultAZtoHW), 10000)
	// base^1
	h.setTprobForDiff(0)
	base := append([]float64(nil), h.currTprob...)
	// base^3 (gap of 3 bp -> posDiff 2)
	h.setTprobForDiff(2)
	far := append([]float64(nil), h.currTprob...)
	if matAt(far, 2, stateHW, stateAZ) <= matAt(base, 2, stateHW, stateAZ) {
		t.Errorf("HW->AZ should grow with distance: base=%g far=%g",
			matAt(base, 2, stateHW, stateAZ), matAt(far, 2, stateHW, stateAZ))
	}
}

// TestRoh_RecRateModifier checks that -M/--rec-rate scales the
// off-diagonal transitions by the interval cross-over probability.
func TestRoh_RecRateModifier(t *testing.T) {
	m := newTprobModifier(RohOptions{RecRate: 1e-3}, nil, "1")
	if m == nil {
		t.Fatal("expected a tprob modifier for RecRate>0")
	}
	tp := baseTprob(DefaultHWtoAZ, DefaultAZtoHW)
	before := matAt(tp, 2, stateHW, stateAZ)
	m.apply(1000, 1100, tp) // 100-bp gap, ci = 100*1e-3 = 0.1
	after := matAt(tp, 2, stateHW, stateAZ)
	if after >= before {
		t.Errorf("rec-rate modifier should shrink HW->AZ here: before=%g after=%g", before, after)
	}
}

// TestRoh_ViterbiTrainingConverges runs --viterbi-training over a
// fixture and confirms the Baum-Welch loop terminates and still
// decodes per-site states.
func TestRoh_ViterbiTrainingConverges(t *testing.T) {
	var out bytes.Buffer
	r, err := Roh(strings.NewReader(fixtureRoh), &out, RohOptions{
		GTsOnly:         30,
		ViterbiTraining: 1e-8,
		HWtoAZ:          1e-3,
		AZtoHW:          1e-3,
	})
	if err != nil {
		t.Fatalf("Roh --viterbi-training: %v", err)
	}
	if len(r.Sites) != 8 {
		t.Errorf("expected 8 decoded sites, got %d", len(r.Sites))
	}
}

// TestRoh_EstimateAF estimates AF from the cohort genotypes and
// confirms a site is processed when --estimate-AF is active.
func TestRoh_EstimateAF(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	1000	.	A	T	.	.	.	GT	1/1	0/1	0/0
chr1	2000	.	C	G	.	.	.	GT	1/1	0/1	0/0
`
	var out bytes.Buffer
	r, err := Roh(strings.NewReader(src), &out, RohOptions{
		GTsOnly:    30,
		EstimateAF: "-",
		Samples:    []string{"S1"},
	})
	if err != nil {
		t.Fatalf("Roh --estimate-AF: %v", err)
	}
	// 6 alt + 6 ref alleles per site over 3 samples -> AF 0.5; both
	// sites must be analysed (not skipped for missing AF).
	if len(r.Sites) != 2 {
		t.Errorf("expected 2 sites with --estimate-AF, got %d", len(r.Sites))
	}
}

// TestRoh_BufferSizeParse pins parseBufferSize against upstream's
// argument grammar.
func TestRoh_BufferSizeParse(t *testing.T) {
	cases := []struct {
		spec     string
		nsamples int
		wantMax  int
		wantOlap int
		wantErr  bool
	}{
		{"", 1, 0, 0, false},
		{"1000", 1, 1000, 10, false},
		{"1000,50", 1, 1000, 50, false},
		{"bad", 1, 0, 0, true},
		{"1000,-1", 1, 0, 0, true},
	}
	for _, c := range cases {
		gotMax, gotOlap, err := parseBufferSize(c.spec, c.nsamples)
		if (err != nil) != c.wantErr {
			t.Errorf("parseBufferSize(%q): err=%v wantErr=%v", c.spec, err, c.wantErr)
			continue
		}
		if err != nil {
			continue
		}
		if gotMax != c.wantMax || gotOlap != c.wantOlap {
			t.Errorf("parseBufferSize(%q) = (%d,%d), want (%d,%d)",
				c.spec, gotMax, gotOlap, c.wantMax, c.wantOlap)
		}
	}
}

// TestRoh_BufferSizeMatchesUnbuffered verifies that running with a
// modest --buffer-size still yields the same region as the unbuffered
// pass for the upstream fixture.
func TestRoh_BufferSizeMatchesUnbuffered(t *testing.T) {
	vcf := filepath.Join(rohGoldenDir, "roh.1.vcf.gz")
	var unbuf, buf bytes.Buffer
	if _, err := RohFile(vcf, &unbuf, RohOptions{GTsOnly: 30, AFDflt: AFDfltPtr(0.4), OutputTypes: "r"}); err != nil {
		t.Fatalf("unbuffered: %v", err)
	}
	if _, err := RohFile(vcf, &buf, RohOptions{GTsOnly: 30, AFDflt: AFDfltPtr(0.4), OutputTypes: "r", BufferSize: "500,20"}); err != nil {
		t.Fatalf("buffered: %v", err)
	}
	if stripComments(unbuf.String()) != stripComments(buf.String()) {
		t.Errorf("buffered output differs from unbuffered:\nunbuf:\n%s\nbuf:\n%s",
			stripComments(unbuf.String()), stripComments(buf.String()))
	}
}
