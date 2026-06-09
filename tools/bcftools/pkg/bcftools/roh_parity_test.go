package bcftools

import (
	"bytes"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// rohGoldenDir is where the upstream bcftools `roh` fixtures are
// vendored: roh.1.vcf.gz and roh.1.tab.gz. The expected regions are
// produced live by the upstream `bcftools roh` binary at test time, not
// from committed golden files.
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

// TestRoh_UpstreamParity replays the seven `roh` invocations from
// reference_code/bcftools/test/test.pl (lines 1094-1100), running BOTH the
// live upstream `bcftools roh` binary and the Go RohFile port on the same
// fixture and comparing their output in-process (after dropping the comment
// header, exactly as test.pl's `| grep -v ^#` does). The upstream binary is
// built on demand; a build failure is fatal, never skipped.
//
// test.pl invocations:
//
//	-Or -G30 --AF-dflt 0.4
//	-Or -G30 --AF-file roh.1.tab.gz
//	-Or -G30 --AF-file roh.1.tab.gz --ignore-homref
//	-G30 --AF-dflt 0.4 -r 1:100174876-100318245
//	-G30 --AF-dflt 0.4 -r 1:100174876-100318245 --ignore-homref
//	-G30 --AF-dflt 0.4 -r 1:... --ignore-homref --include-noalt
//	-G30 --AF-dflt 0.4 -r 1:100174876-100318245 --include-noalt
func TestRoh_UpstreamParity(t *testing.T) {
	bin := upstreamBcftools(t)
	vcf := filepath.Join(rohGoldenDir, "roh.1.vcf.gz")
	tab := filepath.Join(rohGoldenDir, "roh.1.tab.gz")
	region := "1:100174876-100318245"

	// The upstream binary needs a .csi index alongside the VCF to honour
	// `-r`, and a .tbi alongside the AF tab file. The Go port reads
	// sequentially and needs neither. Build indexed copies in a temp dir so
	// the committed fixture tree stays index-free.
	idxVCF := indexedVCFCopy(t, bin, vcf)
	idxTab := indexedTabCopy(t, bin, tab)

	cases := []struct {
		name string
		args []string // upstream CLI args
		opts RohOptions
	}{
		{
			name: "Or_AFdflt",
			args: []string{"-Or", "-G30", "--AF-dflt", "0.4"},
			opts: RohOptions{GTsOnly: 30, AFDflt: AFDfltPtr(0.4), OutputTypes: "r"},
		},
		{
			name: "Or_AFfile",
			args: []string{"-Or", "-G30", "--AF-file", idxTab},
			opts: RohOptions{GTsOnly: 30, AFFile: tab, OutputTypes: "r"},
		},
		{
			name: "Or_AFfile_ignoreHomref",
			args: []string{"-Or", "-G30", "--AF-file", idxTab, "--ignore-homref"},
			opts: RohOptions{GTsOnly: 30, AFFile: tab, IgnoreHomRef: true, OutputTypes: "r"},
		},
		{
			name: "region_AFdflt",
			args: []string{"-G30", "--AF-dflt", "0.4", "-r", region},
			opts: RohOptions{GTsOnly: 30, AFDflt: AFDfltPtr(0.4), Regions: []string{region}},
		},
		{
			name: "region_AFdflt_ignoreHomref",
			args: []string{"-G30", "--AF-dflt", "0.4", "-r", region, "--ignore-homref"},
			opts: RohOptions{GTsOnly: 30, AFDflt: AFDfltPtr(0.4), Regions: []string{region}, IgnoreHomRef: true},
		},
		{
			name: "region_AFdflt_ignoreHomref_noalt",
			args: []string{"-G30", "--AF-dflt", "0.4", "-r", region, "--ignore-homref", "--include-noalt"},
			opts: RohOptions{GTsOnly: 30, AFDflt: AFDfltPtr(0.4), Regions: []string{region}, IgnoreHomRef: true, IncludeNoalt: true},
		},
		{
			name: "region_AFdflt_noalt",
			args: []string{"-G30", "--AF-dflt", "0.4", "-r", region, "--include-noalt"},
			opts: RohOptions{GTsOnly: 30, AFDflt: AFDfltPtr(0.4), Regions: []string{region}, IncludeNoalt: true},
		},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			// Live upstream invocation (uses the indexed VCF copy).
			cmd := exec.Command(bin, append([]string{"roh"}, append(c.args, idxVCF)...)...)
			var upOut, upErr bytes.Buffer
			cmd.Stdout = &upOut
			cmd.Stderr = &upErr
			if err := cmd.Run(); err != nil {
				t.Fatalf("upstream bcftools roh %v: %v\n%s", c.args, err, upErr.String())
			}
			want := stripComments(upOut.String())

			// Go port.
			var out bytes.Buffer
			if _, err := RohFile(vcf, &out, c.opts); err != nil {
				t.Fatalf("RohFile: %v", err)
			}
			got := stripComments(out.String())

			if got != want {
				t.Errorf("output mismatch for %v\n--- got (go) ---\n%s\n--- want (upstream) ---\n%s",
					c.args, got, want)
			}
		})
	}
}

// indexedVCFCopy copies a bgzipped VCF into a temp dir and runs
// `bcftools index` so the upstream binary can honour region (`-r`)
// queries. It returns the temp-dir VCF path.
func indexedVCFCopy(t *testing.T, bin, vcf string) string {
	t.Helper()
	dir := t.TempDir()
	dst := filepath.Join(dir, filepath.Base(vcf))
	if err := copyFile(vcf, dst); err != nil {
		t.Fatalf("copy %s: %v", vcf, err)
	}
	cmd := exec.Command(bin, "index", dst)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("bcftools index %s: %v\n%s", dst, err, out)
	}
	return dst
}

// indexedTabCopy copies a bgzipped CHR\tPOS\tREF,ALT\tAF allele-frequency
// table into a temp dir and tabix-indexes it (`-s1 -b2 -e2`, matching
// upstream's bcf_sr_set_targets), returning the temp-dir path. The tabix
// binary is the one built next to the vendored htslib.
func indexedTabCopy(t *testing.T, bcftoolsBin, tab string) string {
	t.Helper()
	dir := t.TempDir()
	dst := filepath.Join(dir, filepath.Base(tab))
	if err := copyFile(tab, dst); err != nil {
		t.Fatalf("copy %s: %v", tab, err)
	}
	referenceCode := filepath.Dir(filepath.Dir(bcftoolsBin)) // .../reference_code
	tabix := filepath.Join(referenceCode, "htslib", "tabix")
	cmd := exec.Command(tabix, "-s1", "-b2", "-e2", dst)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("tabix %s: %v\n%s", dst, err, out)
	}
	return dst
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
