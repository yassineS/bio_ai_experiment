package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// Live-oracle parity tests for the shared -r/-R/-t/-T region/target selection in
// the native plugins. Every case drives BOTH the genuine upstream bcftools
// binary and OUR port through their CLIs with an upstream-accepted argv, then
// diffs the provenance-stripped stdout byte-for-byte (assertStdoutParity).
//
// Two upstream realities shape the harness:
//
//   - -r/--regions and -R/--regions-file use the index, so upstream REQUIRES a
//     bgzip-compressed, tabix-indexed input. OUR port streams a fully decoded
//     VCF and applies an equivalent overlap filter, so it reads the plain VCF.
//     Region cases therefore feed upstream the indexed copy (bgzipAndIndex) and
//     ours the plain fixture.
//   - -t/--targets and -T/--targets-file are a streaming positional filter that
//     upstream applies to plain VCF, so both binaries read the same plain file.
//
// The semantic difference the cases pin down (verified against upstream
// 1.23.1): -r is span-OVERLAP based while -t is record-START based. In
// overlaps.vcf the deletion at chr1:150 (REF=CGT) spans 150..152; `-r
// chr1:151-151` INCLUDES it (overlap) but `-t chr1:151-151` EXCLUDES it (its
// start, 150, is not 151).

// rtRegionFile writes a region-list / BED file with the given lines into dir and
// returns its path.
func rtRegionFile(t *testing.T, dir, name string, lines ...string) string {
	t.Helper()
	path := filepath.Join(dir, name)
	var buf []byte
	for _, l := range lines {
		buf = append(buf, l...)
		buf = append(buf, '\n')
	}
	if err := os.WriteFile(path, buf, 0o644); err != nil {
		t.Fatalf("write region file %s: %v", path, err)
	}
	return path
}

// TestNativePluginRegionSelection exercises the -r/-R region (overlap) selection
// across every native plugin that now honours it, comparing OUR plain-input run
// against upstream's indexed-input run.
func TestNativePluginRegionSelection(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	overlaps := parityFixture(t, "overlaps.vcf")
	gt := parityFixture(t, "gt_plugins.vcf")

	dir := t.TempDir()
	overlapsGz := bgzipAndIndex(t, bin, overlaps, filepath.Join(dir, "overlaps.vcf.gz"))
	gtGz := bgzipAndIndex(t, bin, gt, filepath.Join(dir, "gt_plugins.vcf.gz"))

	// The synced reader (the plugin region/target path) parses -R/-T files as
	// TSV: a .bed path is 0-based half-open, any other file is 1-based inclusive
	// (chr<TAB>beg<TAB>end). The chr:beg-end colon syntax is only valid for
	// inline -r/-t strings, not files.
	bedAll := rtRegionFile(t, dir, "all.bed", "chr1\t0\t10000", "chr2\t0\t10000")
	regionList := rtRegionFile(t, dir, "regions.txt", "chr1\t150\t405", "chr2\t100\t200")
	bedChr1 := rtRegionFile(t, dir, "chr1.bed", "chr1\t99\t405")

	// Each plugin is run-style, so pluginCLIArgs builds `+name <opts> <file>`.
	// The region options are appended to the plugin's option list; the file is
	// the indexed copy for upstream and the plain copy for ours.
	type plug struct {
		name    string
		fixture string // plain fixture for ours
		gz      string // indexed fixture for upstream
		base    []string
	}
	// check-sparsity is covered separately (TestNativePluginCheckSparsityRegion):
	// its -R file format is verbatim region-list strings, not the synced-reader
	// TSV the plugins below use.
	plugs := []plug{
		{"remove-overlaps", overlaps, overlapsGz, []string{"-M", "OLAP"}},
		{"prune", overlaps, overlapsGz, []string{"-n", "1", "-N", "1st", "-w", "100bp"}},
		{"smpl-stats", gt, gtGz, nil},
		{"indel-stats", gt, gtGz, nil},
		{"contrast", gt, gtGz, []string{"-0", "S1,S2", "-1", "S3,S4"}},
	}
	regionCases := []struct {
		name string
		args []string // region option(s), appended after base; file substituted in
	}{
		{"r_single", []string{"-r", "chr1:151-151"}}, // overlap: indel 150..152 included
		{"r_range", []string{"-r", "chr1:200-405"}},  // a mid-file window
		{"r_commalist", []string{"-r", "chr1:150,chr2"}},
		{"r_openended", []string{"-r", "chr1:300-"}},   // chr:beg-
		{"r_wholefile", []string{"-r", "chr1,chr2"}},   // spans the whole file
		{"r_nothing", []string{"-r", "chr9:1-100"}},    // matches nothing
		{"R_bed_all", []string{"-R", "__BED_ALL__"}},   // BED file, whole genome
		{"R_bed_chr1", []string{"-R", "__BED_CHR1__"}}, // BED file, subset
		{"R_regionlist", []string{"-R", "__REGLIST__"}},
	}
	subst := func(args []string) []string {
		out := make([]string, len(args))
		for i, a := range args {
			switch a {
			case "__BED_ALL__":
				out[i] = bedAll
			case "__BED_CHR1__":
				out[i] = bedChr1
			case "__REGLIST__":
				out[i] = regionList
			default:
				out[i] = a
			}
		}
		return out
	}
	for _, p := range plugs {
		for _, rc := range regionCases {
			p, rc := p, rc
			t.Run(p.name+"_"+rc.name, func(t *testing.T) {
				opts := append(append([]string{}, p.base...), subst(rc.args)...)
				// Feed BOTH binaries the same indexed .vcf.gz: upstream's -r
				// requires the index, and OUR port reads bgzipped VCF
				// transparently, so the report plugins (smpl-stats, indel-stats)
				// echo an identical input path in their CMD line.
				argv := pluginCLIArgs(p.name, p.gz, opts)
				assertStdoutParity(t, bin, argv, argv)
			})
		}
	}
}

// TestNativePluginCheckSparsityRegion checks check-sparsity's -r/-R selection,
// which (uniquely) groups and LABELS its per-region report by the verbatim
// region token. Its -R file is a verbatim region-list (colon syntax) upstream,
// NOT the synced-reader TSV. A colon line is matched byte-for-byte against
// upstream; a TSV/BED line is the fix-on-port case (see
// TestNativePluginCheckSparsityRegionBEDFixOnPort below). Upstream needs the
// indexed input; ours reads the same .vcf.gz.
func TestNativePluginCheckSparsityRegion(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	dir := t.TempDir()
	gtGz := bgzipAndIndex(t, bin, gt, filepath.Join(dir, "gt_plugins.vcf.gz"))
	// Verbatim region-list file (colon syntax — what check-sparsity's
	// hts_readlist + tbx_itr_querys accepts).
	regionList := rtRegionFile(t, dir, "cs_regions.txt", "chr1:200-405", "chr2")

	cases := []struct {
		name string
		args []string
	}{
		{"r_range", []string{"-n", "1", "-r", "chr1:200-405"}},
		{"r_bare_chrom", []string{"-n", "1", "-r", "chr2"}},
		{"r_commalist", []string{"-n", "1", "-r", "chr1:100-100,chr2"}},
		{"r_nothing", []string{"-n", "1", "-r", "chr9"}},
		{"r_n2", []string{"-n", "2", "-r", "chr1"}},
		{"R_regionlist", []string{"-n", "1", "-R", regionList}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			argv := pluginCLIArgs("check-sparsity", gtGz, tc.args)
			assertStdoutParity(t, bin, argv, argv)
		})
	}
}

// TestNativePluginCheckSparsityRegionBEDFixOnPort pins the fix-on-port for the
// documented upstream bug (docs/UPSTREAM_BUGS.md#bcftools-check-sparsity-regions-file):
// upstream reads -R with hts_readlist + tbx_itr_querys, so a BED/TSV line fails
// to parse and upstream emits NOTHING. Our port parses BED/TSV the synced-reader
// way. We validate against the REAL upstream binary by feeding upstream the
// equivalent colon region-list it DOES parse (a BED "chr 0 10000" line is
// 0-based half-open => 1-based "chr:1-10000"); ours given the BED file must match
// upstream given that colon-equivalent, label and all.
func TestNativePluginCheckSparsityRegionBEDFixOnPort(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	dir := t.TempDir()
	gtGz := bgzipAndIndex(t, bin, gt, filepath.Join(dir, "gt_plugins.vcf.gz"))
	bedFile := rtRegionFile(t, dir, "cs.bed", "chr1\t0\t10000", "chr2\t0\t10000")
	// The colon region-list upstream parses correctly; equivalent 1-based windows.
	colonEquiv := rtRegionFile(t, dir, "cs_bed_equiv.txt", "chr1:1-10000", "chr2:1-10000")

	// 1) Document the bug: upstream with the BED -R file produces no output.
	upBED := runUpstreamPluginCmd(t, bin, pluginCLIArgs("check-sparsity", gtGz, []string{"-n", "1", "-R", bedFile}))
	if len(bytes.TrimSpace(stripProvenanceBytes(upBED))) != 0 {
		t.Fatalf("expected upstream to silently emit nothing for a BED -R file, got:\n%s", upBED)
	}
	// 2) Fix-on-port: ours given the BED file equals upstream given the
	//    colon-equivalent region-list it can actually parse.
	oursArgv := pluginCLIArgs("check-sparsity", gtGz, []string{"-n", "1", "-R", bedFile})
	upArgv := pluginCLIArgs("check-sparsity", gtGz, []string{"-n", "1", "-R", colonEquiv})
	assertStdoutParity(t, bin, upArgv, oursArgv)
	// Sanity: the fixed output is actually non-empty (otherwise the equality is
	// trivially the empty string and would not prove the fix).
	if len(bytes.TrimSpace(stripProvenanceBytes(runOursPluginCmd(t, oursArgv)))) == 0 {
		t.Fatalf("fix-on-port produced empty output; expected a per-region sparse report")
	}
}

// TestNativePluginTargetSelection exercises the -t/-T target (start-position)
// selection, including the leading-^ negation, across every native plugin that
// now honours it. Targets work on plain VCF for both binaries.
func TestNativePluginTargetSelection(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	overlaps := parityFixture(t, "overlaps.vcf")
	gt := parityFixture(t, "gt_plugins.vcf")

	dir := t.TempDir()
	// Synced-reader TSV target files: a 3-column region-list (1-based
	// inclusive), a .bed file (0-based half-open), and a two-column chr<TAB>pos
	// file (single 1-based position).
	regionList := rtRegionFile(t, dir, "targets.txt", "chr1\t150\t405", "chr2\t100\t200")
	bedFile := rtRegionFile(t, dir, "targets.bed", "chr1\t149\t405")
	twoCol := rtRegionFile(t, dir, "twocol.txt", "chr1\t150", "chr1\t300")

	type plug struct {
		name    string
		fixture string
		base    []string
	}
	plugs := []plug{
		{"remove-overlaps", overlaps, []string{"-M", "OLAP"}},
		{"prune", overlaps, []string{"-n", "1", "-N", "1st", "-w", "100bp"}},
		{"smpl-stats", gt, nil},
		{"indel-stats", gt, nil},
		{"contrast", gt, []string{"-0", "S1,S2", "-1", "S3,S4"}},
	}
	targetCases := []struct {
		name string
		args []string
	}{
		{"t_single", []string{"-t", "chr1:151-151"}}, // start-based: indel start 150 not in 151 => excluded
		{"t_startboundary", []string{"-t", "chr1:150-152"}},
		{"t_commalist", []string{"-t", "chr1:150,chr2"}},
		{"t_wholefile", []string{"-t", "chr1,chr2"}},
		{"t_nothing", []string{"-t", "chr9:1-100"}},
		{"t_negate", []string{"-t", "^chr1"}}, // exclude all of chr1
		{"t_negate_range", []string{"-t", "^chr1:100-305"}},
		{"T_regionlist", []string{"-T", "__REGLIST__"}},
		{"T_bed", []string{"-T", "__BED__"}},
		{"T_twocol", []string{"-T", "__TWOCOL__"}},
		{"T_negate_file", []string{"-T", "^__REGLIST__"}},
	}
	subst := func(args []string) []string {
		out := make([]string, len(args))
		for i, a := range args {
			switch a {
			case "__REGLIST__":
				out[i] = regionList
			case "^__REGLIST__":
				out[i] = "^" + regionList
			case "__BED__":
				out[i] = bedFile
			case "__TWOCOL__":
				out[i] = twoCol
			default:
				out[i] = a
			}
		}
		return out
	}
	for _, p := range plugs {
		for _, tc := range targetCases {
			p, tc := p, tc
			t.Run(p.name+"_"+tc.name, func(t *testing.T) {
				opts := append(append([]string{}, p.base...), subst(tc.args)...)
				argv := pluginCLIArgs(p.name, p.fixture, opts)
				assertStdoutParity(t, bin, argv, argv)
			})
		}
	}
}

// TestNativePluginGuessPloidyRegion checks that guess-ploidy honours -r/-R
// region selection. guess-ploidy's -t is --tag (an INFO field), NOT targets, so
// only -r/-R are region/target cases here. Upstream needs the indexed input.
func TestNativePluginGuessPloidyRegion(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	gp := parityFixture(t, "guess_ploidy.vcf")
	dir := t.TempDir()
	gpGz := bgzipAndIndex(t, bin, gp, filepath.Join(dir, "guess_ploidy.vcf.gz"))
	regionList := rtRegionFile(t, dir, "gp_regions.txt", "X\t2699000\t2800000")

	cases := []struct {
		name string
		args []string
	}{
		{"r_par1", []string{"-r", "X:2699521-2800000"}},
		{"r_nothing", []string{"-r", "Y"}},
		{"r_commalist", []string{"-r", "X:2700000-2700200,X:2700300-2700400"}},
		{"R_regionlist", []string{"-R", regionList}},
		// guess-ploidy's -t is --tag, not targets: pass it through to the plugin.
		{"t_is_tag", []string{"-t", "GT"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			// The -t GT case uses the plain file for both (no region jump).
			if tc.name == "t_is_tag" {
				argv := pluginCLIArgs("guess-ploidy", gp, tc.args)
				assertStdoutParity(t, bin, argv, argv)
				return
			}
			upArgv := pluginCLIArgs("guess-ploidy", gpGz, tc.args)
			ourArgv := pluginCLIArgs("guess-ploidy", gp, tc.args)
			assertStdoutParity(t, bin, upArgv, ourArgv)
		})
	}
}

// TestNativePluginMendelian2RegionTarget checks that the +mendelian2 plugin
// honours -r/-R/-t/-T. The shared engine streams plain VCF, so target cases use
// the plain fixture for both; region cases give upstream the indexed copy.
func TestNativePluginMendelian2RegionTarget(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "color_chrs_trio.vcf")
	dir := t.TempDir()
	trioGz := bgzipAndIndex(t, bin, trio, filepath.Join(dir, "trio.vcf.gz"))
	pfm := []string{"-p", "CHILD,FATHER,MOTHER"}

	t.Run("t_chr1", func(t *testing.T) {
		args := append(append([]string{}, pfm...), "-t", "chr1")
		argv := pluginCLIArgs("mendelian2", trio, args)
		assertStdoutParity(t, bin, argv, argv)
	})
	t.Run("t_negate", func(t *testing.T) {
		args := append(append([]string{}, pfm...), "-t", "^chr1:1-6000000")
		argv := pluginCLIArgs("mendelian2", trio, args)
		assertStdoutParity(t, bin, argv, argv)
	})
	t.Run("r_chr1", func(t *testing.T) {
		args := append(append([]string{}, pfm...), "-r", "chr1:1-6000000")
		up := pluginCLIArgs("mendelian2", trioGz, args)
		ours := pluginCLIArgs("mendelian2", trio, args)
		assertStdoutParity(t, bin, up, ours)
	})
}

// TestNativePluginTrioStatsTargets checks that +trio-stats honours -t/-T
// streaming targets (previously rejected). trio-stats needs a PED/PFM; the trio
// fixture's three samples form a single trio via -P/--pfm.
func TestNativePluginTrioStatsTargets(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "color_chrs_trio.vcf")
	pfm := []string{"-P", "CHILD,FATHER,MOTHER"}

	cases := []struct {
		name string
		args []string
	}{
		{"t_chr1", []string{"-t", "chr1"}},
		{"t_negate", []string{"-t", "^chr1:1-6000000"}},
		{"t_range", []string{"-t", "chr1:1-6000000"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			args := append(append([]string{}, pfm...), tc.args...)
			argv := pluginCLIArgs("trio-stats", trio, args)
			assertStdoutParity(t, bin, argv, argv)
		})
	}
}

// TestNativePluginIsecGTRegionTarget checks that +isecGT applies region/target
// selection to BOTH input streams (matching upstream's synced reader). Upstream
// requires both inputs bgzipped+indexed; ours streams the plain files.
func TestNativePluginIsecGTRegionTarget(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	aPlain := parityFixture(t, "isecgt_a.vcf")
	bPlain := parityFixture(t, "isecgt_b.vcf")
	dir := t.TempDir()
	aGz := bgzipAndIndex(t, bin, aPlain, filepath.Join(dir, "a.vcf.gz"))
	bGz := bgzipAndIndex(t, bin, bPlain, filepath.Join(dir, "b.vcf.gz"))

	// -r region: both files indexed for upstream, plain for ours.
	t.Run("r_chr1", func(t *testing.T) {
		assertStdoutParity(t, bin,
			[]string{"+isecGT", "-r", "chr1", aGz, bGz},
			[]string{"+isecGT", "-r", "chr1", aPlain, bPlain})
	})
	// -t targets work on plain VCF for both.
	t.Run("t_chr1", func(t *testing.T) {
		assertStdoutParity(t, bin,
			[]string{"+isecGT", "-t", "chr1", aGz, bGz},
			[]string{"+isecGT", "-t", "chr1", aPlain, bPlain})
	})
}

// TestNativePluginSplitScatterRegionTarget checks that +split and +scatter, the
// multi-output plugins, apply region/target selection to the records before
// fanning them into per-file outputs. Both write into their own temp dirs which
// are then diffed file-by-file (assertDirParity). Upstream needs the indexed
// input for -r; ours reads the plain file.
func TestNativePluginSplitScatterRegionTarget(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	dir := t.TempDir()
	gtGz := bgzipAndIndex(t, bin, gt, filepath.Join(dir, "gt_plugins.vcf.gz"))

	// split: -t target (plain VCF both).
	t.Run("split_t_chr1", func(t *testing.T) {
		assertDirParity(t, bin,
			func(d string) []byte {
				return runUpstreamPluginCmd(t, bin, []string{"+split", "-o", d, "-t", "chr1", gt})
			},
			func(d string) []byte {
				return runOursPluginCmd(t, []string{"+split", "-o", d, "-t", "chr1", gt})
			})
	})
	// split: -r region (indexed upstream, plain ours).
	t.Run("split_r_chr1", func(t *testing.T) {
		assertDirParity(t, bin,
			func(d string) []byte {
				return runUpstreamPluginCmd(t, bin, []string{"+split", "-o", d, "-r", "chr1", gtGz})
			},
			func(d string) []byte {
				return runOursPluginCmd(t, []string{"+split", "-o", d, "-r", "chr1", gt})
			})
	})
	// scatter by region label, with a -t pre-filter applied first.
	t.Run("scatter_t_chr1", func(t *testing.T) {
		assertDirParity(t, bin,
			func(d string) []byte {
				return runUpstreamPluginCmd(t, bin, []string{"+scatter", "-o", d, "-s", "chr1,chr2", "-t", "chr1", gt})
			},
			func(d string) []byte {
				return runOursPluginCmd(t, []string{"+scatter", "-o", d, "-s", "chr1,chr2", "-t", "chr1", gt})
			})
	})
	t.Run("scatter_r_chr1", func(t *testing.T) {
		assertDirParity(t, bin,
			func(d string) []byte {
				return runUpstreamPluginCmd(t, bin, []string{"+scatter", "-o", d, "-s", "chr1,chr2", "-r", "chr1", gtGz})
			},
			func(d string) []byte {
				return runOursPluginCmd(t, []string{"+scatter", "-o", d, "-s", "chr1,chr2", "-r", "chr1", gt})
			})
	})
}
