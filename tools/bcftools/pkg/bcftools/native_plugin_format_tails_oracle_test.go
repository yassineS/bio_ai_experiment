package bcftools

import (
	"bytes"
	"path/filepath"
	"testing"
)

// Live-oracle parity tests for the remaining format/output-mode tails closed in
// this wave: ad-bias --clean-vcf / -f, remove-overlaps --missing / -Ot,
// tag2tag --LXX-to-XX localized expansion, guess-ploidy -g/--genome, and the
// af-dist bin-from-file input. Each drives BOTH the upstream binary and our
// port through their CLIs with the same upstream-accepted argv and diffs stdout
// byte-for-byte (after stripping provenance), exactly the established harness.

// TestNativePluginAdBiasCleanVCF checks the -c/--clean-vcf allele-subsetting
// mode (the VCF subset to only ALT alleles passing -t) and the -f convert-format
// mode (a query-style column appended to each FT report line). ad-bias is a
// generic init/process plugin, so its options follow `--`.
func TestNativePluginAdBiasCleanVCF(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	multi := parityFixture(t, "adbias_multi.vcf") // multiallelic + AD/PL/GT
	gt := parityFixture(t, "gt_plugins.vcf")
	pairs := parityFixture(t, "adbias_pairs.txt")
	cases := []struct {
		fixture string
		args    []string
	}{
		// --clean-vcf across thresholds: drops failing ALT alleles (remapping
		// AC/AD/PL/GT) and whole sites where nothing passes.
		{multi, []string{"-s", pairs, "-c", "-t", "1e-9"}},
		{multi, []string{"-s", pairs, "-c", "-t", "1e-6"}},
		{multi, []string{"-s", pairs, "-c", "-t", "1e-3"}},
		{multi, []string{"-s", pairs, "-c", "-t", "0.5"}},
		{multi, []string{"-s", pairs, "-c", "-t", "0.9"}},
		{multi, []string{"-s", pairs, "-c", "-t", "1"}},
		{gt, []string{"-s", pairs, "-c"}},
		{gt, []string{"-s", pairs, "-c", "-t", "1"}},
		// --clean-vcf long form requires (and ignores) an argument upstream.
		{multi, []string{"-s", pairs, "--clean-vcf", "ignored", "-t", "0.9"}},
		// -f convert-format appended to each FT line (report mode).
		{gt, []string{"-s", pairs, "-t", "1", "-f", "%CHROM:%POS"}},
		{gt, []string{"-s", pairs, "-t", "1", "-f", "%INFO/AC"}},
		{gt, []string{"-s", pairs, "-t", "1", "-f", "[%SAMPLE=%GT ]"}},
		{gt, []string{"-s", pairs, "-t", "1", "-f", `%REF\t%ALT\t%QUAL`}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "ad-bias", tc.args...)
		})
	}
}

// TestNativePluginRemoveOverlapsMissing checks the min(QUAL) --missing modes
// (the default scalar 0 and the DP coverage heuristic) and the -Ot/-Otz text
// list output, both stdout and the Processed/Removed stderr summary. The
// fixture deliberately uses co-located / cleanly-overlapping records so the
// upstream oracle is self-consistent (see UPSTREAM_BUGS.md
// bcftools-remove-overlaps-minqual-stale-mark).
func TestNativePluginRemoveOverlapsMissing(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	dp := parityFixture(t, "overlaps_dp.vcf")
	overlaps := parityFixture(t, "overlaps.vcf")
	cases := [][]string{
		{"-m", "min(QUAL)", "--missing", "0"},
		{"-m", "min(QUAL)", "--missing", "DP"},
		{"-m", "min(QUAL)", "--missing", "DP", "-M", "MQ"},
		{"-m", "min(QUAL)", "--missing", "DP", "--reverse"},
		// -Ot text list (chr,pos) in every removal/marking mode.
		{"-O", "t"},
		{"-m", "dup", "-O", "t"},
		{"-M", "OLAP", "-O", "t"},
		{"-m", "min(QUAL)", "-O", "t"},
		{"-m", "min(QUAL)", "--missing", "DP", "-O", "t"},
	}
	for _, args := range cases {
		args := args
		fixture := dp
		// The pure -Ot cases also run on the original overlaps fixture.
		t.Run("dp_"+joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "remove-overlaps", args...)
			assertPluginStderrParity(t, bin, fixture, "remove-overlaps", args...)
		})
	}
	for _, args := range [][]string{{"-O", "t"}, {"-m", "dup", "-O", "t"}, {"-M", "X", "-O", "t"}} {
		args := args
		t.Run("overlaps_"+joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, overlaps, "remove-overlaps", args...)
			assertPluginStderrParity(t, bin, overlaps, "remove-overlaps", args...)
		})
	}
}

// TestNativePluginTag2TagLocalized checks the localized-allele expansion
// (--LXX-to-XX, --LPL-to-PL, --LAD-to-AD) using FORMAT/LAA to map per-sample
// localized indices back to global PL (Number=G) and AD (Number=R), with -r
// (replace), -d (defaults) and -s (skip-nalt). tag2tag is a generic plugin.
func TestNativePluginTag2TagLocalized(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	laa := parityFixture(t, "tag2tag_laa.vcf")
	cases := [][]string{
		{"--LXX-to-XX"},
		{"--LXX-to-XX", "-r"},
		{"--LAD-to-AD"},
		{"--LAD-to-AD", "-r"},
		{"--LPL-to-PL"},
		{"--LPL-to-PL", "-r"},
		{"--LXX-to-XX", "-d", "AD:0,PL:0"},
		{"--LXX-to-XX", "-d", "AD:.,PL:."},
		{"--LXX-to-XX", "-s", "3"},
		{"--LXX-to-XX", "-s", "3", "-r"},
		{"--LXX-to-XX", "-s", "2"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, laa, "tag2tag", args...)
		})
	}
}

// TestNativePluginGuessPloidyGenome checks the -g/--genome shortcut. Upstream
// reads the bgzipped+indexed input (a region requires the index); our port
// streams the plain file. The b37/b38 presets (no-chr prefix) match the X
// contig; hg19/hg38 (chr prefix) match nothing, leaving every sample "U" for
// both binaries.
func TestNativePluginGuessPloidyGenome(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	gp := parityFixture(t, "guess_ploidy.vcf")
	dir := t.TempDir()
	gpGz := bgzipAndIndex(t, bin, gp, filepath.Join(dir, "guess_ploidy.vcf.gz"))
	cases := []struct {
		name string
		args []string
	}{
		{"b37", []string{"-g", "b37"}},
		{"b38", []string{"-g", "b38"}},
		{"hg19", []string{"-g", "hg19"}},
		{"hg38", []string{"-g", "hg38"}},
		{"b37_verbose", []string{"-g", "b37", "-v"}},
		{"b37_gt", []string{"-g", "b37", "-t", "GT"}},
		{"genome_long", []string{"--genome", "b37"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			upArgv := pluginCLIArgs("guess-ploidy", gpGz, tc.args)
			ourArgv := pluginCLIArgs("guess-ploidy", gp, tc.args)
			assertStdoutParity(t, bin, upArgv, ourArgv)
		})
	}
}

// TestNativePluginAfDistBinFile checks reading bin boundaries from a file (one
// float per line), the af-dist -p/-d file form that bin_init reads via
// hts_readlist. af-dist is a generic plugin.
func TestNativePluginAfDistBinFile(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	af := parityFixture(t, "afdist.vcf")
	probBins := parityFixture(t, "afdist_probbins.txt")
	devBins := parityFixture(t, "afdist_devbins.txt")
	cases := [][]string{
		{"-p", probBins},
		{"-d", devBins},
		{"-p", probBins, "-d", devBins},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, af, "af-dist", args...)
		})
	}
}

// TestNativePluginFormatTailsUnsupported asserts the deliberately unsupported
// localized direction (--XX-to-LXX, a todo upstream too) and the bad
// guess-ploidy genome value still fail with a clean Init error.
func TestNativePluginFormatTailsUnsupported(t *testing.T) {
	laa := parityFixture(t, "tag2tag_laa.vcf")
	gp := parityFixture(t, "guess_ploidy.vcf")
	cases := []struct {
		name string
		args []string
	}{
		{"tag2tag", []string{"--XX-to-LXX"}},
		{"tag2tag", []string{"--PL-to-LPL"}},
		{"tag2tag", []string{"--AD-to-LAD"}},
		{"guess-ploidy", []string{"-g", "notagenome"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name+"_"+joinArgs(tc.args), func(t *testing.T) {
			fixture := laa
			if tc.name == "guess-ploidy" {
				fixture = gp
			}
			var out, errBuf bytes.Buffer
			err := RunPlugin(PluginOptions{
				Name:         tc.name,
				Args:         tc.args,
				InputFile:    fixture,
				OutputFormat: OutputVCF,
			}, &out, &errBuf)
			if err == nil {
				t.Fatalf("expected an unsupported error for +%s %v, got nil", tc.name, tc.args)
			}
		})
	}
}
