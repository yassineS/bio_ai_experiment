package bcftools

import (
	"bytes"
	"testing"
)

// Live-oracle parity tests for the batch-2 native plugins (allele-length,
// dosage, impute-info, variant-distance, check-ploidy, vcf2table, and the
// tag2tag --QR-QA-to-QS mode). Each builds the genuine upstream bcftools via
// buildBcftools() and compares its stdout, byte-for-byte after stripping
// provenance, against the native pipeline.

func TestNativePluginDosage(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	multi := parityFixture(t, "multi.vcf")
	cases := []struct {
		fixture string
		args    []string
	}{
		{gt, nil}, // default PL,GL,GT -> PL handler
		{gt, []string{"-t", "GT"}},
		{gt, []string{"-t", "PL"}},
		{gt, []string{"-t", "GL,GT"}}, // GL absent -> falls to GT
		{gt, []string{"-t", "GT,PL"}},
		{multi, []string{"-t", "GT"}},
		{multi, nil}, // no PL/GL header -> GT handler
	}
	for _, tc := range cases {
		tc := tc
		t.Run(joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "dosage", tc.args...)
		})
	}
}

func TestNativePluginImputeInfo(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gp := parityFixture(t, "impute_info.vcf")
	gt := parityFixture(t, "gt_plugins.vcf") // no GP -> all sites unchanged
	t.Run("gp", func(t *testing.T) {
		assertPluginParity(t, bin, gp, "impute-info")
	})
	t.Run("no-gp", func(t *testing.T) {
		assertPluginParity(t, bin, gt, "impute-info")
	})
}

func TestNativePluginAlleleLength(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	for _, fx := range []string{"basic.vcf", "multi.vcf", "gt_plugins.vcf"} {
		fx := fx
		t.Run(fx, func(t *testing.T) {
			fixture := parityFixture(t, fx)
			// allele-length takes no plugin args; invoke as `+allele-length file`.
			argv := pluginCLIArgs("allele-length", fixture, nil)
			up := runUpstreamPlugin(t, bin, argv...)
			ours := runOursPlugin(t, argv...)
			if !bytes.Equal(up, ours) {
				t.Fatalf("+allele-length diverges from upstream (argv=%v)\n--- upstream ---\n%s\n--- ours ---\n%s", argv, snippet(up, 1200), snippet(ours, 1200))
			}
		})
	}
}

func TestNativePluginCheckPloidy(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	xy := parityFixture(t, "fixploidy_xy.vcf")
	for _, fx := range []string{gt, xy} {
		fx := fx
		for _, args := range [][]string{nil, {"-m"}} {
			args := args
			t.Run(joinArgs(args), func(t *testing.T) {
				assertPluginParity(t, bin, fx, "check-ploidy", args...)
			})
		}
	}
}

func TestNativePluginVariantDistance(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	multi := parityFixture(t, "multi.vcf")
	dup := parityFixture(t, "dup.vcf")
	cases := []struct {
		fixture string
		args    []string
	}{
		{gt, nil},
		{gt, []string{"-d", "nearest"}},
		{gt, []string{"-d", "fwd"}},
		{gt, []string{"-d", "rev"}},
		{gt, []string{"-d", "both"}},
		{gt, []string{"-n", "NEAR"}},
		{gt, []string{"-d", "both", "-n", "DD"}},
		{multi, nil},
		{multi, []string{"-d", "both"}},
		{dup, []string{"-d", "both"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			// variant-distance is a run()-style plugin: its options precede the
			// input file rather than following `--`. assertPluginParity builds
			// that exact upstream-accepted form via pluginCLIArgs and drives
			// both binaries through their CLIs.
			assertPluginParity(t, bin, tc.fixture, "variant-distance", tc.args...)
		})
	}
}

func TestNativePluginVcf2Table(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	basic := parityFixture(t, "basic.vcf")
	multi := parityFixture(t, "multi.vcf")
	gt := parityFixture(t, "gt_plugins.vcf")
	cases := []struct {
		fixture string
		args    []string
	}{
		{basic, nil},
		{multi, nil},
		{gt, nil},
		{gt, []string{"-x", "INFO"}},
		{gt, []string{"-x", "GT"}},
		{gt, []string{"-x", "GTTYPES"}},
		{gt, []string{"-x", "HOM_REF,NO_CALL"}},
		{gt, []string{"-x", "VC"}},
		{multi, []string{"-x", "HET"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "vcf2table", tc.args...)
		})
	}
}

func TestNativePluginTag2TagQRQA(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	qrqa := parityFixture(t, "tag2tag_qrqa.vcf")
	for _, args := range [][]string{{"--QR-QA-to-QS"}, {"--QR-QA-to-QS", "-r"}} {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, qrqa, "tag2tag", args...)
		})
	}
}

// shortName returns the base file name without extension for subtest labels.
func shortName(path string) string {
	for i := len(path) - 1; i >= 0; i-- {
		if path[i] == '/' {
			path = path[i+1:]
			break
		}
	}
	for i := len(path) - 1; i >= 0; i-- {
		if path[i] == '.' {
			return path[:i]
		}
	}
	return path
}

// TestNativePluginUnsupportedModes asserts that the modes deliberately left
// unsupported fail with a clean error from Init rather than diverging silently.
func TestNativePluginUnsupportedModes(t *testing.T) {
	gt := parityFixture(t, "gt_plugins.vcf")
	cases := []struct {
		name string
		args []string
	}{
		// frameshifts is now natively supported (see
		// native_plugin_gvcfz_frameshifts_oracle_test.go); only its missing-option
		// error remains a clean failure, exercised elsewhere.
		{"tag2tag", []string{"--LXX-to-XX"}},
		{"tag2tag", []string{"--XX-to-LXX"}},
		{"tag2tag", []string{"--PL-to-LPL"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name+"_"+joinArgs(tc.args), func(t *testing.T) {
			var out, errBuf bytes.Buffer
			err := RunPlugin(PluginOptions{
				Name:         tc.name,
				Args:         tc.args,
				InputFile:    gt,
				OutputFormat: OutputVCF,
			}, &out, &errBuf)
			if err == nil {
				t.Fatalf("expected an unsupported error for +%s %v, got nil", tc.name, tc.args)
			}
		})
	}
}
