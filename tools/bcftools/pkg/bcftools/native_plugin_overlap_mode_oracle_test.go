package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// readFileT reads a file, failing the test on error.
func readFileT(t *testing.T, path string) []byte {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return b
}

// upstreamPluginSucceeds reports whether the upstream bcftools binary exits zero
// for argv. It is used to document a fix-on-port where upstream errors outright.
func upstreamPluginSucceeds(t *testing.T, bin string, argv []string) bool {
	t.Helper()
	cmd := exec.Command(bin, argv...)
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
	return cmd.Run() == nil
}

// Live-oracle parity tests for the shared --regions-overlap / --targets-overlap
// MODE option (pos|0, record|1, variant|2) now honoured by every region/target-
// aware native plugin, and for the gtisec ploidy>2 fix-on-port. Every case drives
// the genuine upstream bcftools binary and OUR port through their CLIs and diffs
// the provenance-stripped stdout byte-for-byte.
//
// The discriminating record in overlaps.vcf is the deletion at chr1:150
// (REF=CGT, ALT=C) spanning 150..152. Under the three modes htslib derives the
// record interval [beg,end] as:
//
//   - pos     (0): beg=end=150              (POS only)
//   - record  (1): beg=150, end=152         (POS..POS+rlen-1)
//   - variant (2): beg=151, end=152         (POS+off..POS+rlen-1; off=len of the
//                  shared "C" prefix between REF and ALT, =1)
//
// so the windows below separate the modes:
//
//   -r chr1:150-150 : pos keeps it, record keeps it, variant DROPS it (151..152
//                     does not reach 150).
//   -r chr1:151-151 : pos DROPS it (150!=151), record keeps it, variant keeps it.
//
// Upstream applies --regions-overlap via the index (so the input must be
// bgzipped+indexed); --targets-overlap is a streaming positional filter on plain
// VCF. The default modes are regions=record(1) and targets=pos(0), matching
// htslib's bcf_sr_init.

// TestNativePluginRegionsOverlapMode pins --regions-overlap pos|record|variant
// against upstream on the `contrast` plugin (the per-record plugin that upstream
// genuinely wires to bcf_sr_set_opt(BCF_SR_REGIONS_OVERLAP)). The discriminating
// deletion at chr1:150 (CGT>C, span 150..152, variant span 151..152) separates
// the three modes against the windows below.
func TestNativePluginRegionsOverlapMode(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	fx := parityFixture(t, "overlap_modes_contrast.vcf")
	dir := t.TempDir()
	gz := bgzipAndIndex(t, bin, fx, filepath.Join(dir, "ov.vcf.gz"))
	base := []string{"-0", "S1,S2", "-1", "S3,S4"}

	cases := []struct {
		name string
		args []string
	}{
		// chr1:150-150 distinguishes variant from {pos,record}.
		{"w150_pos", []string{"-r", "chr1:150-150", "--regions-overlap", "pos"}},
		{"w150_record", []string{"-r", "chr1:150-150", "--regions-overlap", "record"}},
		{"w150_variant", []string{"-r", "chr1:150-150", "--regions-overlap", "variant"}},
		// chr1:151-151 distinguishes pos from {record,variant}.
		{"w151_pos", []string{"-r", "chr1:151-151", "--regions-overlap", "0"}},
		{"w151_record", []string{"-r", "chr1:151-151", "--regions-overlap", "1"}},
		{"w151_variant", []string{"-r", "chr1:151-151", "--regions-overlap", "2"}},
		// Insertion at chr2:400 (G>GAAA): record span 400..400, variant span empty.
		{"ins_w400_record", []string{"-r", "chr2:400-400", "--regions-overlap", "record"}},
		{"ins_w400_variant", []string{"-r", "chr2:400-400", "--regions-overlap", "variant"}},
		// A wider window with the variant mode (exercises plumbing, no separation).
		{"w_wide_variant", []string{"-r", "chr1:140-160", "--regions-overlap", "variant"}},
	}
	for _, c := range cases {
		c := c
		t.Run(c.name, func(t *testing.T) {
			argv := pluginCLIArgs("contrast", gz, append(append([]string{}, base...), c.args...))
			assertStdoutParity(t, bin, argv, argv)
		})
	}
}

// TestNativePluginTargetsOverlapMode pins --targets-overlap pos|record|variant on
// contrast. Targets are a streaming positional filter on plain VCF for both
// binaries; the default (pos) is the start-position semantics, while
// record/variant extend the match to the record / true-variant span.
func TestNativePluginTargetsOverlapMode(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	fx := parityFixture(t, "overlap_modes_contrast.vcf")
	base := []string{"-0", "S1,S2", "-1", "S3,S4"}

	cases := []struct {
		name string
		args []string
	}{
		// -t chr1:151-151: pos drops the deletion (start 150), record/variant keep
		// it (its span reaches 151).
		{"w151_pos", []string{"-t", "chr1:151-151", "--targets-overlap", "pos"}},
		{"w151_record", []string{"-t", "chr1:151-151", "--targets-overlap", "record"}},
		{"w151_variant", []string{"-t", "chr1:151-151", "--targets-overlap", "variant"}},
		// -t chr1:150-150: pos/record keep, variant drops (variant span 151..152).
		{"w150_pos", []string{"-t", "chr1:150-150", "--targets-overlap", "0"}},
		{"w150_record", []string{"-t", "chr1:150-150", "--targets-overlap", "1"}},
		{"w150_variant", []string{"-t", "chr1:150-150", "--targets-overlap", "2"}},
		// Negated target with a non-default mode.
		{"negate_record", []string{"-t", "^chr1:150-150", "--targets-overlap", "record"}},
	}
	for _, c := range cases {
		c := c
		t.Run(c.name, func(t *testing.T) {
			argv := pluginCLIArgs("contrast", fx, append(append([]string{}, base...), c.args...))
			assertStdoutParity(t, bin, argv, argv)
		})
	}
}

// TestNativePluginSplitScatterOverlapMode pins --regions-overlap on +split and
// +scatter (multi-output plugins that upstream wires to BCF_SR_REGIONS_OVERLAP).
// Output dirs are diffed file-by-file.
func TestNativePluginSplitScatterOverlapMode(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	fx := parityFixture(t, "overlap_modes_contrast.vcf")
	dir := t.TempDir()
	gz := bgzipAndIndex(t, bin, fx, filepath.Join(dir, "ov.vcf.gz"))

	t.Run("split_r_variant", func(t *testing.T) {
		assertDirParity(t, bin,
			func(d string) []byte {
				return runUpstreamPluginCmd(t, bin, []string{"+split", "-o", d, "-r", "chr1:150-150", "--regions-overlap", "variant", gz})
			},
			func(d string) []byte {
				return runOursPluginCmd(t, []string{"+split", "-o", d, "-r", "chr1:150-150", "--regions-overlap", "variant", fx})
			})
	})
	t.Run("scatter_r_pos", func(t *testing.T) {
		assertDirParity(t, bin,
			func(d string) []byte {
				return runUpstreamPluginCmd(t, bin, []string{"+scatter", "-o", d, "-s", "chr1,chr2", "-r", "chr1:150-150", "--regions-overlap", "pos", gz})
			},
			func(d string) []byte {
				return runOursPluginCmd(t, []string{"+scatter", "-o", d, "-s", "chr1,chr2", "-r", "chr1:150-150", "--regions-overlap", "pos", fx})
			})
	})
}

// TestNativePluginMendelian2OverlapMode is a FIX-ON-PORT: upstream mendelian2
// advertises --regions-overlap/--targets-overlap (and registers the long options
// in getopt) but has NO matching case in its switch, so passing either prints
// usage and exits non-zero. Our port honours them via the shared filter. We
// document the upstream bug, then validate the fix: ours with the DEFAULT-
// equivalent mode (regions=record, targets=pos) equals upstream WITHOUT the
// option (which applies those very defaults). See docs/UPSTREAM_BUGS.md.
func TestNativePluginMendelian2OverlapMode(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "color_chrs_trio.vcf")
	dir := t.TempDir()
	trioGz := bgzipAndIndex(t, bin, trio, filepath.Join(dir, "trio.vcf.gz"))
	pfm := []string{"-p", "CHILD,FATHER,MOTHER"}

	// 1) Document the upstream getopt bug.
	t.Run("upstream_rejects_regions_overlap", func(t *testing.T) {
		argv := pluginCLIArgs("mendelian2", trioGz, append(append([]string{}, pfm...), "-r", "chr1", "--regions-overlap", "record"))
		if upstreamPluginSucceeds(t, bin, argv) {
			t.Fatal("expected upstream mendelian2 to reject --regions-overlap (getopt bug)")
		}
	})

	// 2) Fix-on-port: ours with the default-equivalent mode matches upstream
	//    without the option (upstream defaults are regions=record, targets=pos).
	t.Run("regions_record_equals_default", func(t *testing.T) {
		upArgv := pluginCLIArgs("mendelian2", trioGz, append(append([]string{}, pfm...), "-r", "chr1:1-6000000"))
		ourArgv := pluginCLIArgs("mendelian2", trio, append(append([]string{}, pfm...), "-r", "chr1:1-6000000", "--regions-overlap", "record"))
		assertStdoutParity(t, bin, upArgv, ourArgv)
	})
	t.Run("targets_pos_equals_default", func(t *testing.T) {
		upArgv := pluginCLIArgs("mendelian2", trio, append(append([]string{}, pfm...), "-t", "chr1"))
		ourArgv := pluginCLIArgs("mendelian2", trio, append(append([]string{}, pfm...), "-t", "chr1", "--targets-overlap", "pos"))
		assertStdoutParity(t, bin, upArgv, ourArgv)
	})
}

// TestNativePluginTrioDNM3Targets pins trio-dnm3's newly wired -t/-T streaming
// targets (previously rejected) against upstream. trio-dnm3 needs a trio (-P)
// and FORMAT/GT; the NAIVE model is used (--use-NAIVE) so the output is
// byte-exact. Targets work on plain VCF for both binaries.
func TestNativePluginTrioDNM3Targets(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "color_chrs_trio.vcf")
	dir := t.TempDir()
	regionList := rtRegionFile(t, dir, "dnm_targets.txt", "chr1\t1\t6000000")
	base := []string{"-p", "CHILD,FATHER,MOTHER", "--use-NAIVE"}

	cases := []struct {
		name string
		args []string
	}{
		{"t_chr1", []string{"-t", "chr1"}},
		{"t_range", []string{"-t", "chr1:1-6000000"}},
		{"t_negate", []string{"-t", "^chr1:1-6000000"}},
		{"T_regionlist", []string{"-T", regionList}},
		{"t_overlap_record", []string{"-t", "chr1:1-6000000", "--targets-overlap", "record"}},
	}
	for _, c := range cases {
		c := c
		t.Run(c.name, func(t *testing.T) {
			args := append(append([]string{}, base...), c.args...)
			argv := pluginCLIArgs("trio-dnm3", trio, args)
			assertStdoutParity(t, bin, argv, argv)
		})
	}
}

// TestNativePluginTrioDNM3Regions pins trio-dnm3's -r/-R region selection now
// routed through the shared filter (including --regions-overlap). Upstream needs
// the indexed input; ours reads the plain file.
func TestNativePluginTrioDNM3Regions(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "color_chrs_trio.vcf")
	dir := t.TempDir()
	trioGz := bgzipAndIndex(t, bin, trio, filepath.Join(dir, "trio.vcf.gz"))
	base := []string{"-p", "CHILD,FATHER,MOTHER", "--use-NAIVE"}

	cases := []struct {
		name string
		args []string
	}{
		{"r_chr1", []string{"-r", "chr1:1-6000000"}},
		{"r_variant", []string{"-r", "chr1:1-6000000", "--regions-overlap", "variant"}},
	}
	for _, c := range cases {
		c := c
		t.Run(c.name, func(t *testing.T) {
			args := append(append([]string{}, base...), c.args...)
			up := pluginCLIArgs("trio-dnm3", trioGz, args)
			ours := pluginCLIArgs("trio-dnm3", trio, args)
			assertStdoutParity(t, bin, up, ours)
		})
	}
}

// TestNativePluginTrioDNM3OutputFile pins trio-dnm3's -o/--output FILE (previously
// rejected): the bytes written to the file must equal what would go to stdout.
// We run ours twice — once to stdout, once with -o FILE — and require the file
// contents to match the stdout run AND the upstream stdout run (provenance
// stripped).
func TestNativePluginTrioDNM3OutputFile(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "color_chrs_trio.vcf")
	dir := t.TempDir()
	base := []string{"-p", "CHILD,FATHER,MOTHER", "--use-NAIVE"}

	// Upstream stdout (the reference).
	upArgv := pluginCLIArgs("trio-dnm3", trio, base)
	up := stripProvenanceBytes(runUpstreamPluginCmd(t, bin, upArgv))

	// Ours to stdout.
	oursStdout := stripProvenanceBytes(runOursPluginCmd(t, pluginCLIArgs("trio-dnm3", trio, base)))
	if !bytes.Equal(up, oursStdout) {
		t.Fatalf("ours stdout diverges from upstream\n--- upstream ---\n%s\n--- ours ---\n%s",
			snippet(up, 2000), snippet(oursStdout, 2000))
	}

	// Ours with -o FILE: the file must contain identical bytes.
	outFile := filepath.Join(dir, "dnm_out.vcf")
	oFileArgv := pluginCLIArgs("trio-dnm3", trio, append(append([]string{}, base...), "-o", outFile))
	if got := runOursPluginCmd(t, oFileArgv); len(bytes.TrimSpace(got)) != 0 {
		t.Fatalf("expected -o to write nothing to stdout, got:\n%s", got)
	}
	fileBytes := stripProvenanceBytes(readFileT(t, outFile))
	if !bytes.Equal(up, fileBytes) {
		t.Fatalf("-o FILE contents differ from upstream stdout\n--- upstream ---\n%s\n--- file ---\n%s",
			snippet(up, 2000), snippet(fileBytes, 2000))
	}
}

// TestNativePluginGTisecPloidyGt2 pins the fix-on-port for gtisec ploidy>2.
// Upstream errors ("gtisec does not support ploidy higher than 2"); our port
// generalises the unordered-multiset intersection to arbitrary ploidy. Because
// GTisec's output is keyed purely by the sample-sharing PARTITION (counts per
// sample subset, sample names — never the genotype values), a partition-
// preserving diploid remap of the polyploid fixture is a valid upstream oracle:
// upstream on the diploid-equivalent must byte-match ours on the polyploid file.
// See docs/UPSTREAM_BUGS.md#bcftools-gtisec-ploidy.
func TestNativePluginGTisecPloidyGt2(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	poly := parityFixture(t, "gtisec_polyploid.vcf")
	diploidEquiv := parityFixture(t, "gtisec_polyploid_diploid_equiv.vcf")

	// 1) Document the bug: upstream errors on the polyploid input (non-zero exit).
	t.Run("upstream_rejects_ploidy_gt2", func(t *testing.T) {
		if upstreamPluginSucceeds(t, bin, pluginCLIArgs("GTisec", poly, nil)) {
			t.Fatal("expected upstream GTisec to fail on ploidy>2 input")
		}
	})

	// 2) Fix-on-port: ours on the polyploid file equals upstream on the
	//    partition-equivalent diploid file, across the output modes.
	for _, args := range [][]string{nil, {"-v"}, {"-m"}, {"-H"}, {"-m", "-v"}, {"-m", "-H"}} {
		args := args
		t.Run("fix_"+joinArgs(args), func(t *testing.T) {
			upArgv := pluginCLIArgs("GTisec", diploidEquiv, args)
			ourArgv := pluginCLIArgs("GTisec", poly, args)
			assertStdoutParity(t, bin, upArgv, ourArgv)
		})
	}
}
