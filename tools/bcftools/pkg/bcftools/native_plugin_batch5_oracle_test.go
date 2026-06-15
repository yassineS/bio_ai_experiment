package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// Live-oracle parity tests for the batch-5 sample-group / per-sample statistics
// native plugins (guess-ploidy, smpl-stats, indel-stats, contrast, ad-bias,
// GTisec, GTsubset). Each builds the genuine upstream bcftools via
// buildBcftools() and drives BOTH that binary and OUR port through their CLIs
// with the SAME upstream-accepted argv, diffing stdout (and, for contrast, the
// stderr summary) byte-for-byte after stripping provenance — the same harness
// as batch 2-4 (pluginCLIArgs / runUpstreamPlugin / runOursPlugin /
// assertPluginParity / stripProvenanceBytes), kept strictly CLI-to-CLI.
//
// Dispatch styles (detected from each plugin's .c, and mirrored by the host via
// IsRunStyleNativePlugin):
//   - guess-ploidy run()-style  => `+guess-ploidy <opts> FILE`
//   - smpl-stats   run()-style  => `+smpl-stats <opts> FILE`
//   - indel-stats  run()-style  => `+indel-stats <opts> FILE`
//   - contrast     run()-style  => `+contrast <opts> FILE`
//   - ad-bias      generic      => `+ad-bias FILE -- <opts>`
//   - GTisec       generic      => `+GTisec FILE -- <opts>`
//   - GTsubset     generic      => `+GTsubset FILE -- <opts>`

// TestNativePluginGuessPloidy checks the per-sample sex/ploidy report from the
// PL, GL and GT likelihood sources, plus the verbose table. The GL verbose mode
// is exercised only for the predicted-sex columns (the final "$4-$5" column can
// differ by one ULP from upstream because libm's pow/log are not bit-identical
// to Go's math; PL and GT verbose match exactly and cover the verbose path).
func TestNativePluginGuessPloidy(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fixture := parityFixture(t, "guess_ploidy.vcf")
	cases := [][]string{
		nil,                         // default PL
		{"-t", "PL"},                // explicit PL
		{"-t", "GL"},                // GL source
		{"-t", "GT"},                // GT source
		{"-v"},                      // verbose PL table
		{"-v2"},                     // verbose PL + DBG lines (getopt attached form)
		{"-t", "GT", "-v"},          // verbose GT table
		{"-t", "GT", "-v2"},         // verbose GT + DBG lines
		{"-t", "GT", "-e", "0.01"},  // custom error rate
		{"--AF-tag", "AF"},          // AF from INFO tag
		{"--AF-tag", "AF", "-v"},    // AF tag, verbose
		{"--AF-dflt", "0.3"},        // custom default AF
		{"-i"},                      // include indels (none here)
		{"-t", "GT", "-e", "0.001"}, // default error rate explicit
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "guess-ploidy", args...)
		})
	}
}

// TestNativePluginSmplStats checks the per-sample and per-site statistics table
// (default "all" filter). The CMD report line embeds the verbatim argv, which is
// identical for both binaries since they share the same upstream-accepted argv.
func TestNativePluginSmplStats(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	for _, fx := range []string{"gt_plugins.vcf", "indels.vcf", "multi.vcf", "basic.vcf"} {
		fx := fx
		t.Run(fx, func(t *testing.T) {
			assertPluginParity(t, bin, parityFixture(t, fx), "smpl-stats")
		})
	}
}

// TestNativePluginIndelStats checks the indel summary, length, VAF and DFRAC
// distributions (default mode). The --max-len/--nvaf knobs are covered with the
// only values upstream's (buggy) --nvaf [0,1] validation accepts.
func TestNativePluginIndelStats(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	indels := parityFixture(t, "indels.vcf")
	cases := []struct {
		fixture string
		args    []string
	}{
		{indels, nil},
		{indels, []string{"--max-len", "5"}},
		{indels, []string{"--max-len", "50"}},
		{indels, []string{"--nvaf", "1"}},
		{indels, []string{"-c", "ANN"}}, // alternate (absent) CSQ tag
		{parityFixture(t, "gt_plugins.vcf"), nil},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "indel-stats", tc.args...)
		})
	}
}

// TestNativePluginContrast checks the per-site association annotations and novel
// allele/genotype lists across two sample groups. Both the annotated VCF
// (stdout) and the summary line (stderr) are parity-checked.
func TestNativePluginContrast(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fixture := parityFixture(t, "gt_plugins.vcf")
	cases := [][]string{
		{"-0", "S1,S2", "-1", "S3,S4"},
		{"-a", "PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT", "-0", "S1,S2", "-1", "S3,S4"},
		{"-a", "NASSOC", "-0", "S1", "-1", "S2,S3,S4"},
		{"-a", "NOVELAL,NOVELGT", "-0", "S3,S4", "-1", "S1,S2"},
		{"-a", "NOVELGT", "-0", "S1", "-1", "S2,S3,S4"},
		{"-0", "S1,S3", "-1", "S2,S4"},
		{"-a", "PASSOC", "-0", "S1,S2,S3", "-1", "S4"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "contrast", args...)
			assertPluginStderrParity(t, bin, fixture, "contrast", args...)
		})
	}

	// @file sample lists for the two groups.
	dir := t.TempDir()
	ctrl := filepath.Join(dir, "ctrl.txt")
	cas := filepath.Join(dir, "case.txt")
	if err := os.WriteFile(ctrl, []byte("S1\nS2\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(cas, []byte("S3\nS4\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	t.Run("file_lists", func(t *testing.T) {
		assertPluginParity(t, bin, fixture, "contrast", "-0", ctrl, "-1", cas)
	})
}

// TestNativePluginAdBias checks the per-pair Fisher-test hit table and the
// trailing summary, covering the depth/alt-depth thresholds and the snp/indel
// variant-type restriction. The sample-pair file is created under testdata.
func TestNativePluginAdBias(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fixture := parityFixture(t, "gt_plugins.vcf")
	pairs := parityFixture(t, "adbias_pairs.txt")
	cases := [][]string{
		nil,
		{"-s", pairs},
		{"-s", pairs, "-t", "1"},
		{"-s", pairs, "-t", "1e-3"},
		{"-s", pairs, "-t", "1", "-d", "20"},
		{"-s", pairs, "-t", "1", "-a", "5"},
		{"-s", pairs, "-t", "1", "-v", "snp"},
		{"-s", pairs, "-t", "1", "-v", "indel"},
	}
	for _, args := range cases {
		args := args
		// ad-bias always requires -s; prepend it for the bare case.
		full := args
		if len(full) == 0 {
			full = []string{"-s", pairs}
		}
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "ad-bias", full...)
		})
	}
}

// TestNativePluginGTisec checks the subset genotype-intersection counts in
// banker's sequence order, across the default, verbose (-v), human-readable (-H)
// and missing-count (-m) modes and their combinations.
func TestNativePluginGTisec(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	for _, fx := range []string{"gt_plugins.vcf", "multi.vcf", "basic.vcf"} {
		fx := fx
		fixture := parityFixture(t, fx)
		for _, args := range [][]string{nil, {"-v"}, {"-m"}, {"-H"}, {"-m", "-v"}, {"-m", "-H"}} {
			args := args
			t.Run(fx+"_"+joinArgs(args), func(t *testing.T) {
				assertPluginParity(t, bin, fixture, "GTisec", args...)
			})
		}
	}
}

// TestNativePluginGTsubset checks the exclusive-genotype-sharing record filter
// across several sample subsets, exercising the raw-GT (phase-aware) comparison.
func TestNativePluginGTsubset(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fixture := parityFixture(t, "gt_plugins.vcf")
	for _, sel := range []string{
		"S1", "S2", "S3", "S4",
		"S1,S2", "S1,S3", "S1,S4", "S2,S3", "S2,S4", "S3,S4",
		"S2,S3,S4", "S1,S2,S3,S4",
	} {
		sel := sel
		t.Run(sel, func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "GTsubset", "-s", sel)
		})
	}
}

// TestNativePluginBatch5Unsupported asserts the deliberately unsupported modes
// fail with a clean Init error rather than diverging silently. These go through
// RunPlugin directly so host-level region/flag interception does not mask the
// plugin's own rejection.
func TestNativePluginBatch5Unsupported(t *testing.T) {
	gt := parityFixture(t, "gt_plugins.vcf")
	gp := parityFixture(t, "guess_ploidy.vcf")
	pairs := parityFixture(t, "adbias_pairs.txt")
	cases := []struct {
		name string
		args []string
	}{
		// guess-ploidy: region/genome jumps and filter expressions are unsupported.
		{"guess-ploidy", []string{"-g", "b37"}},
		{"guess-ploidy", []string{"-r", "X:2699521-154931043"}},
		{"guess-ploidy", []string{"-R", "regions.txt"}},
		{"guess-ploidy", []string{"--include", "QUAL>10"}},
		{"guess-ploidy", []string{"--exclude", "QUAL<10"}},
		// smpl-stats: filter, region and -o file modes are unsupported.
		{"smpl-stats", []string{"-i", "GQ>30"}},
		{"smpl-stats", []string{"-e", "GQ<30"}},
		{"smpl-stats", []string{"-r", "chr1"}},
		{"smpl-stats", []string{"-o", "out.txt"}},
		// indel-stats: filter, PED, region and -o file modes are unsupported.
		{"indel-stats", []string{"-i", "GQ>30"}},
		{"indel-stats", []string{"-p", "trios.ped"}},
		{"indel-stats", []string{"-r", "chr1"}},
		{"indel-stats", []string{"--nvaf", "10"}}, // upstream's [0,1] validation rejects it too
		// contrast: filter, rare-allele enrichment and region modes are unsupported.
		{"contrast", []string{"-0", "S1", "-1", "S2", "-f", "0.001"}},
		{"contrast", []string{"-0", "S1", "-1", "S2", "-i", "QUAL>10"}},
		{"contrast", []string{"-0", "S1", "-1", "S2", "-r", "chr1"}},
		{"contrast", nil}, // missing -0/-1
		// ad-bias: clean-vcf and convert-format modes are unsupported.
		{"ad-bias", []string{"-s", pairs, "-c"}},
		{"ad-bias", []string{"-s", pairs, "-f", "%CHROM"}},
		{"ad-bias", nil}, // missing -s
	}
	fixtureFor := func(name string) string {
		if name == "guess-ploidy" {
			return gp
		}
		return gt
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name+"_"+joinArgs(tc.args), func(t *testing.T) {
			var out, errBuf bytes.Buffer
			err := RunPlugin(PluginOptions{
				Name:         tc.name,
				Args:         tc.args,
				InputFile:    fixtureFor(tc.name),
				OutputFormat: OutputVCF,
			}, &out, &errBuf)
			if err == nil {
				t.Fatalf("expected an unsupported error for +%s %v, got nil", tc.name, tc.args)
			}
		})
	}
}
