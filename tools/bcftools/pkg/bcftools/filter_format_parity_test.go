package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Live upstream-binary parity for the FORMAT/sample-level filter engine.
//
// Each case runs `bcftools view -i/-e <expr>` and `bcftools filter -i/-e
// <expr>` on the SAME fixture through BOTH the real upstream C binary
// (built on demand by upstreamBcftools) and our Go port, then diffs stdout
// byte-for-byte after stripping provenance (##bcftools_*) lines. No committed
// golden file is involved. The expression set covers FMT/DP and FMT/GQ
// comparisons, every GT class keyword, exact GT strings, mixed INFO+FORMAT
// expressions, the QUAL||FMT guard case, regex (~/!~) and negations (-e).

// runUpstreamCmd runs the upstream bcftools binary with argv and returns
// stdout, failing loudly (never skipping) per the project's live-oracle rule.
func runUpstreamCmd(t *testing.T, bin string, argv ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, argv...)
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream %v: %v\nstderr: %s", argv, err, errBuf.String())
	}
	return out.Bytes()
}

// runOursCmd runs OUR built bcftools binary with argv and returns stdout.
func runOursCmd(t *testing.T, argv ...string) []byte {
	t.Helper()
	if ourBinPath == "" {
		t.Fatalf("local bcftools port binary not built; cannot run CLI oracle")
	}
	cmd := exec.Command(ourBinPath, argv...)
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("ours %v: %v\nstderr: %s", argv, err, errBuf.String())
	}
	return out.Bytes()
}

// formatFilterExprs is the shared expression set exercised against view and
// filter. Each must produce byte-identical output (modulo provenance) on both
// the upstream binary and our port.
var formatFilterExprs = []string{
	// FORMAT numeric comparisons.
	`FMT/DP<10`,
	`FMT/DP>=15`,
	`FORMAT/DP>20`,
	`FMT/GQ>38`,
	`FMT/GQ<=20`,
	// GT class keywords.
	`GT="het"`,
	`GT="hom"`,
	`GT="alt"`,
	`GT="ref"`,
	`GT="mis"`,
	`GT="hap"`,
	`GT="HET"`, // case-insensitive
	`GT="rr"`,
	`GT="ra"`,
	`GT="aa"`,
	`GT="aA"`,
	`GT!="het"`,
	// Exact genotype strings.
	`GT="0/1"`,
	`GT="1|1"`,
	`GT="0|0"`,
	// Regex over genotype text.
	`GT~"0"`,
	`GT!~"1"`,
	// Mixed INFO + FORMAT.
	`INFO/AC>1`,
	`INFO/AC>1 && FMT/DP>10`,
	`INFO/AC>1 || FMT/GQ>50`,
	`GT="het" && FMT/DP>15`,
	`GT="alt" && FMT/GQ>=38`,
	`GT="het" || GT="hom"`,
	`QUAL>30 || FMT/GQ>50`, // the QUAL||FMT guard
	`FMT/GQ>30 && FMT/DP>15`,
	`GT="mis" || FMT/DP<6`,
}

func TestFilterFormatViewParity(t *testing.T) {
	bin := upstreamBcftools(t)
	fixture := parityFixture(t, "gt_plugins.vcf")
	for _, mode := range []string{"-i", "-e"} {
		for _, expr := range formatFilterExprs {
			mode, expr := mode, expr
			t.Run(mode+" "+expr, func(t *testing.T) {
				argv := []string{"view", "--no-version", mode, expr, fixture}
				up := runUpstreamCmd(t, bin, argv...)
				ours := runOursCmd(t, argv...)
				if !bytes.Equal(stripProvenanceBytes(up), stripProvenanceBytes(ours)) {
					t.Fatalf("view %s %q diverges from upstream\n--- upstream ---\n%s\n--- ours ---\n%s",
						mode, expr, up, ours)
				}
			})
		}
	}
}

func TestFilterFormatFilterParity(t *testing.T) {
	bin := upstreamBcftools(t)
	fixture := parityFixture(t, "gt_plugins.vcf")
	for _, mode := range []string{"-i", "-e"} {
		for _, expr := range formatFilterExprs {
			mode, expr := mode, expr
			t.Run(mode+" "+expr, func(t *testing.T) {
				argv := []string{"filter", "--no-version", mode, expr, fixture}
				up := runUpstreamCmd(t, bin, argv...)
				ours := runOursCmd(t, argv...)
				if !bytes.Equal(stripProvenanceBytes(up), stripProvenanceBytes(ours)) {
					t.Fatalf("filter %s %q diverges from upstream\n--- upstream ---\n%s\n--- ours ---\n%s",
						mode, expr, up, ours)
				}
				// Also check the soft-filter (-s) variant: header Description
				// and FILTER annotation must match too.
				argvS := []string{"filter", "--no-version", "-s", "LowX", mode, expr, fixture}
				upS := runUpstreamCmd(t, bin, argvS...)
				oursS := runOursCmd(t, argvS...)
				if !bytes.Equal(stripProvenanceBytes(upS), stripProvenanceBytes(oursS)) {
					t.Fatalf("filter -s %s %q diverges from upstream\n--- upstream ---\n%s\n--- ours ---\n%s",
						mode, expr, upS, oursS)
				}
			})
		}
	}
}

// TestFilterFormatRemoveOverlapsParity checks the re-enabled remove-overlaps
// -i/-e filter pre-selection against upstream.
func TestFilterFormatRemoveOverlapsParity(t *testing.T) {
	bin := upstreamBcftools(t)
	fixture := parityFixture(t, "overlaps.vcf")
	cases := [][]string{
		{"-i", "QUAL>10"},
		{"-e", "QUAL<40"},
		{"-i", `GT="het"`},
		{"-m", "dup", "-i", "QUAL>10"},
		{"-m", "min(QUAL)", "-e", "QUAL<25"},
		{"-i", "DP>10"},
		{"-e", "DP<10"},
		{"-M", "FOO", "-i", "QUAL>30"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "remove-overlaps", args...)
		})
	}
}

// TestFilterFormatPruneParity checks the re-enabled prune -i/-e filter
// pre-selection (with deterministic supported window modes) against upstream.
func TestFilterFormatPruneParity(t *testing.T) {
	bin := upstreamBcftools(t)
	fixture := parityFixture(t, "prune.vcf")
	cases := [][]string{
		{"-n", "1", "-N", "1st", "-w", "1000", "-i", "QUAL>10"},
		{"-n", "1", "-N", "1st", "-w", "1000", "-e", "QUAL<40"},
		{"-m", "count=1", "-w", "1000", "-i", "QUAL>20"},
		{"-m", "count=2", "-w", "5000", "-e", "QUAL<30"},
		{"-n", "1", "-N", "1st", "-w", "1000", "-i", `GT="het"`},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "prune", args...)
		})
	}
}

// --- unit tests for the per-sample mask API and GT classification ----------

func fmtVariant() *vcf.Variant {
	// Four samples mirroring gt_plugins.vcf pos 100.
	mk := func(gt, dp, gq string) vcf.Sample {
		return vcf.Sample{Data: map[string]string{"GT": gt, "DP": dp, "GQ": gq}}
	}
	return &vcf.Variant{
		Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"T"}, Qual: 30,
		Filter: []string{"PASS"},
		Info:   map[string]string{"AC": "2", "AN": "8"},
		Format: []string{"GT", "DP", "GQ"},
		Samples: []vcf.Sample{
			mk("0/0", "10", "30"), // S1: ref
			mk("0/1", "20", "30"), // S2: het/alt
			mk("./.", "5", "10"),  // S3: missing
			mk("1|1", "15", "40"), // S4: hom alt
		},
	}
}

func TestEvalSamplesMask(t *testing.T) {
	v := fmtVariant()
	cases := []struct {
		expr     string
		passSite bool
		mask     []bool
	}{
		{`FMT/GQ>38`, true, []bool{false, false, false, true}},
		{`FMT/DP<10`, true, []bool{false, false, true, false}},
		{`GT="het"`, true, []bool{false, true, false, false}},
		{`GT="alt"`, true, []bool{false, true, false, true}},
		{`GT="ref"`, true, []bool{true, false, false, false}},
		{`GT="mis"`, true, []bool{false, false, true, false}},
		// && combines per-sample masks by OR ("can be true in different samples").
		{`GT="het" && FMT/GQ>38`, true, []bool{false, true, false, true}},
		// site INFO term gates the per-sample mask through.
		{`INFO/AC>1 && FMT/DP<10`, true, []bool{false, false, true, false}},
		// site INFO term that fails => no samples.
		{`INFO/AC>5 && FMT/DP<10`, false, []bool{false, false, false, false}},
	}
	for _, c := range cases {
		f, err := CompileFilter(c.expr)
		if err != nil {
			t.Fatalf("CompileFilter(%q): %v", c.expr, err)
		}
		pass, mask := f.EvalSamples(v)
		if pass != c.passSite {
			t.Errorf("EvalSamples(%q) passSite=%v, want %v", c.expr, pass, c.passSite)
		}
		if mask == nil {
			t.Fatalf("EvalSamples(%q) returned nil mask (expected per-sample)", c.expr)
		}
		for i := range c.mask {
			if mask[i] != c.mask[i] {
				t.Errorf("EvalSamples(%q) mask[%d]=%v, want %v (mask=%v)", c.expr, i, mask[i], c.mask[i], mask)
				break
			}
		}
	}
}

func TestEvalSamplesSiteLevelNilMask(t *testing.T) {
	v := fmtVariant()
	f, err := CompileFilter(`INFO/AC>1`)
	if err != nil {
		t.Fatal(err)
	}
	pass, mask := f.EvalSamples(v)
	if !pass {
		t.Fatal("INFO/AC>1 should pass site")
	}
	if mask != nil {
		t.Fatalf("site-level expression should yield a nil per-sample mask, got %v", mask)
	}
}

func TestGTClassification(t *testing.T) {
	cases := []struct {
		gt  string
		cls string
		out bool
	}{
		{"0/1", "het", true}, {"0/0", "het", false}, {"1|1", "het", false},
		{"0/0", "hom", true}, {"0/1", "hom", false}, {"1/1", "hom", true},
		{"0/1", "alt", true}, {"1/1", "alt", true}, {"0/0", "alt", false},
		{"0/0", "ref", true}, {"0/1", "ref", false}, {"./.", "ref", false},
		{"./.", "mis", true}, {"./1", "mis", true}, {"0/1", "mis", false},
		{"0", "hap", true}, {"0/1", "hap", false},
		{"0/0", "rr", true}, {"0/1", "ra", true}, {"1/1", "aa", true}, {"1/2", "aA", true},
	}
	for _, c := range cases {
		got := gtMatchesClass(c.gt, c.cls)
		if got != c.out {
			t.Errorf("gtMatchesClass(%q,%q)=%v, want %v", c.gt, c.cls, got, c.out)
		}
	}
}
