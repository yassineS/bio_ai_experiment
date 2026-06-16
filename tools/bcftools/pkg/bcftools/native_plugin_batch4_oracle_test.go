package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"testing"
)

// Live-oracle parity tests for the batch-4 stateful/windowed native plugins
// (af-dist, check-sparsity, remove-overlaps, prune) and the deliberately
// unsupported gvcfz registration. Each builds the genuine upstream bcftools via
// buildBcftools() and drives BOTH that binary and OUR port through their CLIs
// with the SAME upstream-accepted argv, diffing stdout byte-for-byte after
// stripping provenance — exactly the batch-2/3 harness (pluginCLIArgs /
// runUpstreamPlugin / runOursPlugin / assertPluginParity / stripProvenanceBytes),
// kept strictly CLI-to-CLI.
//
// Dispatch styles (detected from each plugin's .c):
//   - af-dist        generic init/process => `+af-dist FILE -- <opts>`
//   - check-sparsity run()-style          => `+check-sparsity <opts> FILE`
//   - remove-overlaps run()-style         => `+remove-overlaps <opts> FILE`
//   - prune          run()-style          => `+prune <opts> FILE`
//   - gvcfz          run()-style (natively supported; errors only on missing -g)
// pluginCLIArgs builds the matching form automatically from IsRunStyleNativePlugin.

// TestNativePluginAfDist checks the AF/GT probability-distribution tables, which
// are printed to stdout while the VCF output is suppressed. af-dist is a generic
// init/process plugin, so its options follow `--`.
func TestNativePluginAfDist(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	af := parityFixture(t, "afdist.vcf")
	gt := parityFixture(t, "gt_plugins.vcf") // no AF tag -> empty distributions
	cases := []struct {
		fixture string
		args    []string
	}{
		{af, nil},                                             // default AF tag, default bins
		{af, []string{"-s"}},                                  // per-sample HWE log probability
		{af, []string{"-l", "0.3,0.5"}},                       // list genotypes in a probability band
		{af, []string{"-t", "EUR_AF"}},                        // alternate AF tag
		{af, []string{"-d", "0,0.25,0.5,0.75,1"}},             // custom deviation bins
		{af, []string{"-p", "0,0.5,1"}},                       // custom probability bins
		{af, []string{"-s", "-l", "0.3,0.5", "-t", "EUR_AF"}}, // combined
		{gt, nil}, // no AF tag present at all
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "af-dist", tc.args...)
		})
	}
}

// TestNativePluginCheckSparsity checks the per-chromosome "missing samples"
// report, printed to stdout. check-sparsity is a run()-style plugin, so its
// options precede the file.
func TestNativePluginCheckSparsity(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	sparsity := parityFixture(t, "sparsity.vcf")
	gt := parityFixture(t, "gt_plugins.vcf")
	cases := []struct {
		fixture string
		args    []string
	}{
		{sparsity, nil}, // default min_sites=1, group by chromosome
		{sparsity, []string{"-n", "1"}},
		{sparsity, []string{"-n", "2"}},
		{sparsity, []string{"-n", "3"}},
		{sparsity, []string{"-n", "5"}},
		{gt, nil},
		{gt, []string{"-n", "2"}},
		{gt, []string{"-n", "3"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "check-sparsity", tc.args...)
		})
	}
}

// TestNativePluginRemoveOverlaps checks the overlap/dup/min(QUAL) removal and
// marking modes; the emitted VCF is parity-checked on stdout and the
// "Processed/Removed" summary is parity-checked on stderr.
func TestNativePluginRemoveOverlaps(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	overlaps := parityFixture(t, "overlaps.vcf")
	dup := parityFixture(t, "dup.vcf")
	basic := parityFixture(t, "basic.vcf")
	multi := parityFixture(t, "multi.vcf")
	cases := []struct {
		fixture string
		args    []string
	}{
		{overlaps, nil},                                 // default: remove overlaps
		{overlaps, []string{"-m", "overlap"}},           // explicit overlap
		{overlaps, []string{"-m", "dup"}},               // remove duplicates
		{overlaps, []string{"-M", "OLAP"}},              // mark overlaps, keep all
		{overlaps, []string{"-m", "dup", "-M", "DUPM"}}, // mark duplicates
		{overlaps, []string{"--reverse"}},               // invert: keep overlaps
		{overlaps, []string{"-m", "dup", "--reverse"}},  // keep only duplicates
		{overlaps, []string{"-m", "min(QUAL)"}},         // resolve overlaps by QUAL
		{overlaps, []string{"-m", "min(QUAL)", "-M", "MQ"}},
		{dup, nil},
		{dup, []string{"-m", "dup"}},
		{basic, nil},
		{multi, nil},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "remove-overlaps", tc.args...)
			assertPluginStderrParity(t, bin, tc.fixture, "remove-overlaps", tc.args...)
		})
	}
}

// TestNativePluginPrune checks the supported window/count pruning modes: -n with
// the deterministic "1st" selection, -m count= cluster removal, and maxAF with
// an explicit --AF-tag. prune is run()-style.
func TestNativePluginPrune(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	fixture := parityFixture(t, "prune.vcf")
	cases := [][]string{
		{"-n", "1", "-N", "1st", "-w", "100bp"},
		{"-n", "2", "-N", "1st", "-w", "200bp"},
		{"-n", "1", "-N", "1st", "-w", "2"},   // site-count window
		{"-n", "2", "-N", "1st", "-w", "3"},   // site-count window
		{"-n", "1", "-N", "1st", "-w", "1Mb"}, // whole-chromosome window
		{"-m", "count=2", "-w", "100bp"},
		{"-m", "count=3", "-w", "1000bp"},
		{"-m", "count=2", "-w", "60bp"},
		{"-m", "count=1", "-w", "200bp"},
		{"-m", "count=1", "-w", "100bp"},
		{"-m", "count=4", "-w", "2000bp"},
		{"--AF-tag", "AF", "-n", "1", "-w", "100bp"},
		{"--AF-tag", "AF", "-n", "2", "-w", "1000bp"},
		{"--AF-tag", "AF", "-n", "1", "-N", "maxAF", "-w", "50000bp"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "prune", args...)
		})
	}
}

// TestNativePluginBatch4Unsupported asserts the deliberately unsupported paths
// fail with a clean error from Init rather than diverging silently. These go
// through RunPlugin directly (rather than the CLI) so that host-level flag
// interception of region options does not mask the plugin's own rejection.
func TestNativePluginBatch4Unsupported(t *testing.T) {
	prune := parityFixture(t, "prune.vcf")
	sparsity := parityFixture(t, "sparsity.vcf")
	overlaps := parityFixture(t, "overlaps.vcf")
	gt := parityFixture(t, "gt_plugins.vcf")
	cases := []struct {
		name string
		args []string
	}{
		// gvcfz: now natively supported (the FORMAT/GT filter engine is wired in;
		// see native_plugin_gvcfz_frameshifts_oracle_test.go). The only remaining
		// hard error is the missing-required -g option.
		{"gvcfz", nil},
		// prune: every mode (LD/annotation/rand/keep-sites/default-maxAF) is now
		// supported and validated against upstream in
		// native_plugin_prune_oracle_test.go; the only remaining hard error is a
		// genuinely-invalid combination (--keep-sites with --nsites-per-win).
		{"prune", []string{"-n", "1", "-k"}}, // upstream: -k cannot combine with -n
		// remove-overlaps: the --missing DP heuristic and the -Ot/-Otz text-list
		// output are now supported (parity-checked in
		// TestNativePluginRemoveOverlapsMissing); -i/-e filters and -r/-R/-t/-T
		// region/target selection are covered elsewhere. The only remaining hard
		// errors are a bad --mark expression and the --missing DP + non-min(QUAL)
		// combination upstream rejects.
		{"remove-overlaps", []string{"-m", "frobnicate"}},
		{"remove-overlaps", []string{"-m", "dup", "--missing", "DP"}},
	}
	fixtureFor := func(name string) string {
		switch name {
		case "prune":
			return prune
		case "check-sparsity":
			return sparsity
		case "remove-overlaps":
			return overlaps
		default:
			return gt
		}
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

// assertPluginStderrParity asserts that OUR port and upstream produce identical
// stderr when driven through their CLIs with the same argv. It is used for
// remove-overlaps, whose "Processed/Removed" summary is emitted on stderr while
// the VCF goes to stdout. The few provenance/log lines that legitimately differ
// are not produced by remove-overlaps, so a plain byte comparison suffices.
func assertPluginStderrParity(t *testing.T, bin, fixture, name string, args ...string) {
	t.Helper()
	argv := pluginCLIArgs(name, fixture, args)
	up := runPluginStderr(t, exec.Command(bin, argv...), true)
	ours := runPluginStderr(t, exec.Command(ourBinPath, argv...), false)
	if !bytes.Equal(up, ours) {
		t.Fatalf("+%s %v stderr diverges from upstream (argv=%v)\n--- upstream ---\n%s\n--- ours ---\n%s",
			name, args, argv, snippet(up, 800), snippet(ours, 800))
	}
}

// runPluginStderr runs a prepared command with the vendored plugins on the
// environment and returns its stderr, discarding stdout. upstream==true marks
// the upstream binary (used only for the failure label).
func runPluginStderr(t *testing.T, cmd *exec.Cmd, upstream bool) []byte {
	t.Helper()
	if ourBinPath == "" && !upstream {
		t.Fatalf("local bcftools port binary not built; cannot run CLI oracle")
	}
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
	var errBuf bytes.Buffer
	cmd.Stdout = nil
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		side := "ours"
		if upstream {
			side = "upstream"
		}
		t.Fatalf("%s %v: %v\nstderr: %s", side, cmd.Args, err, errBuf.String())
	}
	return errBuf.Bytes()
}
