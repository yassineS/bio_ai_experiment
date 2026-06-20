package bcftools

import (
	"bytes"
	"testing"
)

// Live-oracle parity tests for the batch-6 trio / pedigree native plugins
// (trio-stats, trio-switch-rate, mendelian2). Each builds the genuine upstream
// bcftools via buildBcftools() and drives BOTH that binary and OUR port through
// their CLIs with the SAME upstream-accepted argv, diffing stdout byte-for-byte
// after stripping provenance — the same harness as batches 2-5 (pluginCLIArgs /
// runUpstreamPlugin / runOursPlugin / assertPluginParity / stripProvenanceBytes),
// kept strictly CLI-to-CLI, no committed goldens.
//
// Dispatch styles (detected from each plugin's .c, mirrored by the host via
// IsRunStyleNativePlugin):
//   - trio-stats        run()-style => `+trio-stats <opts> FILE`
//   - mendelian2        run()-style => `+mendelian2 <opts> FILE`
//   - trio-switch-rate  generic     => `+trio-switch-rate FILE -- <opts>`
//
// The four batch-6 plugins that cannot be made byte-reproducible against
// upstream — parental-origin, color-chrs, trio-dnm3 (and the unsupported
// region/target modes) — are exercised in TestNativePluginBatch6Unsupported,
// which asserts a clean Init/Run error rather than silently wrong output.

// TestNativePluginTrioStats checks the per-trio transmission / Mendelian-error
// table, including the streamed MERR and TRANSMITTED debug lines (-d), across
// the PED and -P pfm trio sources and a two-trio PED (which also exercises the
// cmp_trios sort, the PED being listed in reverse index order).
func TestNativePluginTrioStats(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "trio.vcf")
	ped := parityFixture(t, "trio.ped")
	multi := parityFixture(t, "trio_multi.vcf")
	multiPed := parityFixture(t, "trio_multi.ped")
	cases := []struct {
		fixture string
		args    []string
	}{
		{trio, []string{"-p", ped}},
		{trio, []string{"-P", "CHILD,FATHER,MOTHER"}},
		{trio, []string{"-p", ped, "-d", "mendel-errors"}},
		{trio, []string{"-p", ped, "-d", "transmitted"}},
		{trio, []string{"-p", ped, "-d", "mendel-errors,transmitted"}},
		{trio, []string{"-P", "CHILD,FATHER,MOTHER", "-d", "transmitted"}},
		{trio, []string{"-p", ped, "-v", "2"}}, // verbosity passthrough, no effect
		{multi, []string{"-p", multiPed}},      // two trios, reverse-listed -> sort
		{multi, []string{"-p", multiPed, "-d", "mendel-errors,transmitted"}},
		// -i/-e site and FORMAT pre-filter parity. A FORMAT expression is folded to
		// a per-trio verdict (include: all three members match; exclude: none of
		// the three match), matching trio-stats.c. The trio fixture carries GQ.
		{trio, []string{"-p", ped, "-i", "QUAL>10"}},
		{trio, []string{"-p", ped, "-e", "QUAL<10"}},
		{trio, []string{"-p", ped, "-i", `GT="het"`}},
		{trio, []string{"-p", ped, "-e", `GT="het"`}},
		{trio, []string{"-p", ped, "-i", "FMT/GQ>30"}},
		{trio, []string{"-p", ped, "-e", "FMT/GQ>30"}},
		{trio, []string{"-P", "CHILD,FATHER,MOTHER", "-i", `GT="het"`}},
		{multi, []string{"-p", multiPed, "-i", `GT="het"`}},
		{multi, []string{"-p", multiPed, "-e", "FMT/GQ>30"}},
		// Curly-brace multi-threshold expansion: each {a,b,c} element becomes its
		// own FLT* section (and multiple groups combine as a cartesian product),
		// with the per-filter MERR / TRANSMITTED debug lines streamed interleaved
		// per record, matching upstream's parse_filters() + run() loop.
		{trio, []string{"-p", ped, "-i", "FMT/GQ>{10,30}"}},                                    // single list
		{trio, []string{"-p", ped, "-e", "FMT/GQ>{10,30}"}},                                    // EXCLUDE, expanded
		{trio, []string{"-p", ped, "-i", "FMT/GQ>{10,30}", "-d", "mendel-errors,transmitted"}}, // braces + debug
		{trio, []string{"-p", ped, "-i", "QUAL>{10,30} && FMT/GQ>{20,40}"}},                    // two groups: cartesian
		{trio, []string{"-p", ped, "-i", "FMT/GQ>{}"}},                                         // empty list -> "all"
		{multi, []string{"-p", multiPed, "-i", "FMT/GQ>{20,40}"}},                              // two trios, expanded
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "trio-stats", tc.args...)
		})
	}
}

// TestNativePluginTrioSwitchRate checks the per-trio phase-switch table and the
// optional per-population (PED 7th column) rollup. The version/command-line
// banner is removed by the provenance stripper.
func TestNativePluginTrioSwitchRate(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "trio.vcf")
	ped := parityFixture(t, "trio.ped")
	popPed := parityFixture(t, "trio_pop.ped")
	multi := parityFixture(t, "trio_multi.vcf")
	multiPed := parityFixture(t, "trio_multi.ped")
	cases := []struct {
		fixture string
		args    []string
	}{
		{trio, []string{"-p", ped}},
		{trio, []string{"-p", popPed}},    // 7th column -> POP lines
		{multi, []string{"-p", multiPed}}, // two trios + two populations
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "trio-switch-rate", tc.args...)
		})
	}
}

// TestNativePluginMendelian2Plugin checks the `+mendelian2` plugin form, which
// reuses the existing Mendelian2 engine. The text count summary (-m c) and the
// VCF/BCF-emitting modes (a/d/e/E/g/m/M/S and combinations) are parity-checked
// against upstream across the -p pfm and -P ped trio sources and a two-trio
// PED. The annotate mode (-m a) verifies the four MERR/MGOOD/MMISS/MNORULE INFO
// annotations now emitted by the engine.
func TestNativePluginMendelian2Plugin(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "trio.vcf")
	ped := parityFixture(t, "trio.ped")
	multi := parityFixture(t, "trio_multi.vcf")
	multiPed := parityFixture(t, "trio_multi.ped")
	pfm := []string{"-p", "CHILD,FATHER,MOTHER"}
	cases := []struct {
		fixture string
		args    []string
	}{
		{trio, append(append([]string{}, pfm...), "-m", "c")},
		{trio, []string{"-P", ped, "-m", "c"}},
		{trio, append(append([]string{}, pfm...), "-m", "a")},
		{trio, append(append([]string{}, pfm...), "-m", "d")},
		{trio, append(append([]string{}, pfm...), "-m", "e")},
		{trio, append(append([]string{}, pfm...), "-m", "E")},
		{trio, append(append([]string{}, pfm...), "-m", "g")},
		{trio, append(append([]string{}, pfm...), "-m", "m")},
		{trio, append(append([]string{}, pfm...), "-m", "M")},
		{trio, append(append([]string{}, pfm...), "-m", "S")},
		{trio, append(append([]string{}, pfm...), "-m", "ad")},
		{trio, append(append([]string{}, pfm...), "-m", "aE")},
		{trio, append([]string{}, pfm...)}, // default mode (c)
		{multi, []string{"-P", multiPed, "-m", "c"}},
		{multi, []string{"-P", multiPed, "-m", "a"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "mendelian2", tc.args...)
		})
	}
}

// TestNativePluginBatch6Unsupported asserts the deliberately unsupported
// plugins / modes fail with a clean error rather than diverging silently. These
// go through RunPlugin directly so host-level region/flag interception does not
// mask the plugin's own rejection.
func TestNativePluginBatch6Unsupported(t *testing.T) {
	trio := parityFixture(t, "trio.vcf")
	cases := []struct {
		name string
		args []string
	}{
		// trio-dnm3: the NAIVE model is parity-checked in TestNativePluginTrioDNM3Naive
		// and the DMM/ALM/DNG float models in TestNativePluginTrioDNM3FloatModels (with
		// the tolerance-aware proximity helper). parental-origin and color-chrs are
		// fully ported and parity-checked in their own tests. The only remaining
		// trio-dnm3 error path is the missing -p/-P option, plus the AD+PL-only DMM
		// site that has no FORMAT/AD/QM (asserted in the float-model test).
		{"trio-dnm3", nil}, // missing -p/-P
		// trio-stats: the -a/--alt-trios accounting and -o file output are now
		// supported (parity-checked in TestNativePluginTrioStatsAltTrios and
		// TestNativePluginTrioStatsOutputFile); the -i/-e filter modes, the
		// curly-brace multi-threshold expansion and -t/-T streaming targets are
		// parity-checked elsewhere. Only the missing-trio-source error remains.
		{"trio-stats", nil}, // missing -p/-P
		// trio-switch-rate: only -p is supported.
		{"trio-switch-rate", nil}, // missing -p
		{"trio-switch-rate", []string{"-x"}},
		// mendelian2: -W remains unsupported; -r/-R/-t/-T region/target selection
		// is now supported (see TestNativePluginRegionTarget).
		{"mendelian2", nil}, // missing -p/-P
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name+"_"+joinArgs(tc.args), func(t *testing.T) {
			var out, errBuf bytes.Buffer
			err := RunPlugin(PluginOptions{
				Name:         tc.name,
				Args:         tc.args,
				InputFile:    trio,
				OutputFormat: OutputVCF,
			}, &out, &errBuf)
			if err == nil {
				t.Fatalf("expected an unsupported/clean error for +%s %v, got nil", tc.name, tc.args)
			}
		})
	}
}
