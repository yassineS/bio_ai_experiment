package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// Live-oracle parity tests for the "stats wave": the remaining flags of the
// four bcftools stats/contrast plugins closed in this change —
//
//   - indel-stats -p/--ped (de-novo indel mode) and -o/--output FILE
//   - trio-stats  -a/--alt-trios (deferred singleton/doubleton) and -o
//   - smpl-stats  -o/--output FILE
//   - contrast    -f/--max-allele-freq (rare-allele enrichment)
//
// Every case drives BOTH the genuine upstream bcftools (built via
// buildBcftools(), BCFTOOLS_PLUGINS pointed at the vendored .so dir) and OUR
// port through their CLIs with the SAME upstream-accepted argv, diffing the
// report byte-for-byte after stripping provenance — the same harness as
// batches 2-6 (pluginCLIArgs / runUpstreamPlugin / runOursPlugin /
// assertPluginParity / assertPluginStderrParity / stripProvenanceBytes), kept
// strictly CLI-to-CLI with no committed goldens.

// assertStatsOutputFileParity drives both binaries with the SAME -o FILE path
// and compares the FILE contents (not stdout) byte-for-byte after stripping
// provenance. Because the CMD report line echoes the verbatim argv, using the
// identical -o path for both makes the files byte-identical when the report is
// correct. It mirrors upstream report_stats()'s "write to FILE instead of
// stdout" behaviour for the three stdout-style stats plugins.
func assertStatsOutputFileParity(t *testing.T, bin, fixture, name string, args ...string) {
	t.Helper()
	dir := t.TempDir()
	outPath := filepath.Join(dir, "report.txt")

	full := append(append([]string{}, args...), "-o", outPath)

	// Upstream run.
	if err := os.Remove(outPath); err != nil && !os.IsNotExist(err) {
		t.Fatal(err)
	}
	_ = runUpstreamPlugin(t, bin, pluginCLIArgs(name, fixture, full)...)
	up, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatalf("upstream did not write -o file: %v", err)
	}

	// Our run, same path.
	if err := os.Remove(outPath); err != nil && !os.IsNotExist(err) {
		t.Fatal(err)
	}
	_ = runOursPlugin(t, pluginCLIArgs(name, fixture, full)...)
	ours, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatalf("our port did not write -o file: %v", err)
	}

	if !bytes.Equal(stripProvenanceBytes(up), stripProvenanceBytes(ours)) {
		t.Fatalf("+%s %v -o FILE diverges from upstream\n--- upstream (%d bytes) ---\n%s\n--- ours (%d bytes) ---\n%s",
			name, args, len(up), snippet(up, 1200), len(ours), snippet(ours, 1200))
	}

	// The -o FILE bytes must also equal the stdout form modulo the CMD line
	// (which echoes the -o argument): a stronger invariant from report_stats().
	std := runOursPlugin(t, pluginCLIArgs(name, fixture, args)...)
	if !bytes.Equal(stripProvenanceBytes(dropCMD(ours)), stripProvenanceBytes(dropCMD(std))) {
		t.Fatalf("+%s %v: -o FILE bytes differ from stdout (modulo CMD)\n--- file ---\n%s\n--- stdout ---\n%s",
			name, args, snippet(ours, 1200), snippet(std, 1200))
	}
}

// dropCMD removes the single "CMD\t..." report line so the -o-file and stdout
// forms (which legitimately differ only in the echoed output path) can be
// compared for equality.
func dropCMD(b []byte) []byte {
	var keep [][]byte
	for _, line := range bytes.Split(b, []byte("\n")) {
		if bytes.HasPrefix(line, []byte("CMD\t")) {
			continue
		}
		keep = append(keep, line)
	}
	return bytes.Join(keep, []byte("\n"))
}

// TestNativePluginStatsOutputFile checks the -o/--output FILE report parity for
// the three stdout-style stats plugins. The file contents must be byte-identical
// to upstream and (modulo the echoed CMD line) to the stdout form.
func TestNativePluginStatsOutputFile(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	indels := parityFixture(t, "indels.vcf")
	trioIndels := parityFixture(t, "trio_indels.vcf")
	trioIndelsPed := parityFixture(t, "trio_indels.ped")
	trio := parityFixture(t, "trio_multi.vcf")
	trioPed := parityFixture(t, "trio_multi.ped")
	cases := []struct {
		name    string
		fixture string
		args    []string
	}{
		{"smpl-stats", gt, nil},
		{"smpl-stats", indels, nil},
		{"smpl-stats", gt, []string{"-i", `GT="het"`}},
		{"indel-stats", indels, nil},
		{"indel-stats", trioIndels, []string{"-p", trioIndelsPed}},
		{"trio-stats", trio, []string{"-p", trioPed}},
		{"trio-stats", trio, []string{"-p", trioPed, "-a", "1"}},
		{"trio-stats", trio, []string{"-p", trioPed, "-d", "mendel-errors,transmitted"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name+"_"+shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertStatsOutputFileParity(t, bin, tc.fixture, tc.name, tc.args...)
		})
	}
}

// TestNativePluginIndelStatsPED checks the PED de-novo indel mode (-p): the
// stats are restricted to de-novo indels in each trio's child, the SN* "number
// of samples" column reports the trio count, and --alt2ref-DNM widens the DNM
// definition. The trio_indels fixture carries two trios with inherited indels, a
// het DNM insertion, a homozygous DNM, a multiallelic indel DNM and CSQ-tagged
// frameshift/inframe sites.
func TestNativePluginIndelStatsPED(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	fixture := parityFixture(t, "trio_indels.vcf")
	ped := parityFixture(t, "trio_indels.ped")
	cases := [][]string{
		{"-p", ped},
		{"-p", ped, "--alt2ref-DNM"},
		{"-p", ped, "--max-len", "3"},
		{"-p", ped, "--max-len", "50"},
		// Per-trio FORMAT/site filter folding in PED mode (include: all three
		// members pass; exclude: none of the three match), matching
		// indel-stats.c. A site INCLUDE is safe; a per-sample EXCLUDE on FORMAT/GQ
		// exercises the inverted mask without hitting upstream's site-EXCLUDE
		// NULL-smpl_pass segfault.
		{"-p", ped, "-i", "QUAL>=45"},
		{"-p", ped, "-i", "FMT/GQ>30"},
		{"-p", ped, "-e", "FMT/GQ<30"},
		{"-p", ped, "-i", `GT="het"`},
		// Curly-brace expansion combined with -p.
		{"-p", ped, "-i", "QUAL>{40,45,50}"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "indel-stats", args...)
		})
	}
}

// TestNativePluginTrioStatsAltTrios checks the -a/--alt-trios deferred
// singleton/doubleton accounting on a two-trio fixture (where an allele can be
// shared across trios): a singleton/doubleton is counted only if its allele
// appears in at most -a trios at the site. The verbose TRANSMITTED lines are
// also emitted from the deferred loop, so -d transmitted is parity-checked too.
func TestNativePluginTrioStatsAltTrios(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	multi := parityFixture(t, "trio_multi.vcf")
	ped := parityFixture(t, "trio_multi.ped")
	cases := [][]string{
		{"-p", ped, "-a", "1"},
		{"-p", ped, "-a", "2"},
		{"-p", ped, "-a", "3"},
		{"-p", ped, "-a", "0"}, // 0 == unlimited (== default)
		{"-p", ped, "-a", "1", "-d", "transmitted"},
		{"-p", ped, "-a", "2", "-d", "mendel-errors,transmitted"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, multi, "trio-stats", args...)
		})
	}
}

// TestNativePluginContrastEnrichment checks the rare-allele enrichment mode
// (-f): the per-site VCF + annotations are unchanged, and the extra region-wide
// "max_AC/PASSOC/FASSOC/NASSOC:" stderr summary line (Fisher's exact probability
// and control/case non-REF fractions over the pooled minor-allele counts) must
// match upstream. Both an integer (allele-count) threshold and float
// (allele-frequency) thresholds are covered.
func TestNativePluginContrastEnrichment(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	fixture := parityFixture(t, "gt_plugins.vcf")
	cases := [][]string{
		{"-f", "0.5", "-0", "S1,S2", "-1", "S3,S4"},
		{"-f", "0.25", "-0", "S1,S2", "-1", "S3,S4"},
		{"-f", "0.75", "-0", "S1,S2", "-1", "S3,S4"},
		{"-f", "0.001", "-0", "S1,S2", "-1", "S3,S4"}, // floor to 1
		{"-f", "1", "-0", "S1,S2", "-1", "S3,S4"},     // integer count
		{"-f", "2", "-0", "S1,S2", "-1", "S3,S4"},
		{"-f", "3", "-0", "S1,S2", "-1", "S3,S4"},
		{"-f", "0.5", "-0", "S1", "-1", "S2,S3,S4"},
		{"-f", "2", "-a", "NASSOC", "-0", "S1,S2", "-1", "S3,S4"},
		{"-f", "2", "-a", "PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT", "-0", "S1,S2", "-1", "S3,S4"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, fixture, "contrast", args...)
			assertPluginStderrParity(t, bin, fixture, "contrast", args...)
		})
	}
}
