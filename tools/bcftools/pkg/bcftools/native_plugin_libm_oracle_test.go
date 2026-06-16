package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// Live-oracle parity tests for the three previously-"unsupported" libm-adjacent
// trio/CNV plugins now fully (or, for trio-dnm3, partially) ported:
//
//   - parental-origin: byte parity on the %e DBG probabilities and %f quality.
//     Its incomplete-beta binomial tail goes through the in-tree kfunc port
//     (native_kfunc.go / callc.go::kfBetai), which uses htslib's deterministic
//     AS245 kf_lgamma rather than libm, so the dup/del scores are bit-identical
//     to upstream on linux/amd64.
//
//   - color-chrs: byte parity on the <prefix>.dat file. Its HMM Viterbi decode
//     is all IEEE-754 +,-,*,/ (the matrix-power transition chain and the
//     constant-product emissions), so the segmentation reproduces the C HMM
//     exactly. Because it writes a file rather than stdout, this test drives the
//     two binaries into separate prefixes and diffs the produced .dat files.
//
//   - trio-dnm3 (NAIVE model only): byte parity on the FORMAT/DNM and FORMAT/VA
//     annotations. The NAIVE verdict is a pure integer Mendelian-consistency
//     table lookup (no floating point), so it is bit-reproducible. The DMM/ALM/
//     DNG float models remain unsupported (asserted in
//     TestNativePluginBatch6Unsupported).
//
// All comparisons are CLI-to-CLI against the live upstream binary with
// provenance stripped; no committed goldens.

// TestNativePluginParentalOrigin checks the parental-origin summary and DBG
// lines across the dup/del CNV types, the -d debug listing, -g greedy, the -b
// min-binom threshold, and the -i/-e per-trio filters, all byte-for-byte.
func TestNativePluginParentalOrigin(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "trio.vcf")
	del := parityFixture(t, "parental_origin_del.vcf")
	cases := []struct {
		fixture string
		args    []string
	}{
		{trio, []string{"-p", "CHILD,FATHER,MOTHER", "-t", "dup"}},
		{trio, []string{"-p", "CHILD,FATHER,MOTHER", "-t", "dup", "-d"}},
		{trio, []string{"-p", "CHILD,FATHER,MOTHER", "-t", "dup", "-g", "-d"}},
		{trio, []string{"-p", "CHILD,FATHER,MOTHER", "-t", "dup", "-b", "0.3", "-d"}},
		{trio, []string{"-p", "CHILD,FATHER,MOTHER", "-t", "dup", "-b", "0", "-d"}},
		{trio, []string{"-p", "CHILD,FATHER,MOTHER", "-t", "del"}},
		{trio, []string{"-p", "CHILD,FATHER,MOTHER", "-t", "del", "-d"}},
		{del, []string{"-p", "CHILD,FATHER,MOTHER", "-t", "del", "-d"}},
		{del, []string{"-p", "CHILD,FATHER,MOTHER", "-t", "del", "-g", "-d"}},
		// Per-sample filter (does not trigger the upstream site-only-exclude NULL
		// deref; see docs/UPSTREAM_BUGS.md).
		{trio, []string{"-p", "CHILD,FATHER,MOTHER", "-t", "dup", "-d", "-i", "QUAL>10"}},
		{trio, []string{"-p", "CHILD,FATHER,MOTHER", "-t", "dup", "-d", "-i", "FMT/GQ>30"}},
		{trio, []string{"-p", "CHILD,FATHER,MOTHER", "-t", "dup", "-d", "-e", "FMT/GQ>50"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "parental-origin", tc.args...)
		})
	}
}

// TestNativePluginParentalOriginNoNullDeref asserts our fix-on-port for the
// upstream segfault: a site-only -e expression that matches no site must NOT
// crash and must run to completion (upstream 1.23.1 NULL-derefs here, see
// docs/UPSTREAM_BUGS.md). We drive only our port; upstream produces no
// comparable output.
func TestNativePluginParentalOriginNoNullDeref(t *testing.T) {
	if ourBinPath == "" {
		t.Fatalf("local bcftools port binary not built")
	}
	trio := parityFixture(t, "trio.vcf")
	cmd := exec.Command(ourBinPath, "+parental-origin", "-p", "CHILD,FATHER,MOTHER", "-t", "dup", "-e", "QUAL<10", trio)
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("ours +parental-origin -e QUAL<10 should not crash, got: %v\nstderr: %s", err, errBuf.String())
	}
	if !bytes.Contains(out.Bytes(), []byte("\tpaternal\t")) && !bytes.Contains(out.Bytes(), []byte("\tuncertain\t")) && !bytes.Contains(out.Bytes(), []byte("\tmaternal\t")) {
		t.Fatalf("expected a summary line in the output, got:\n%s", out.String())
	}
}

// TestNativePluginColorChrs checks the <prefix>.dat segmentation for both the
// trio (-t) and unrelated (-u) modes, including the large-gap matrix-power
// transition chain and the unphased/missing-GT skip, byte-for-byte against the
// .dat upstream writes.
func TestNativePluginColorChrs(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "color_chrs_trio.vcf")
	unrel := parityFixture(t, "color_chrs_unrel.vcf")
	cases := []struct {
		name    string
		fixture string
		args    []string // plugin options, WITHOUT the -p prefix (added per-run)
	}{
		{"trio", trio, []string{"-t", "MOTHER,FATHER,CHILD"}},
		{"trio_grch38", trio, []string{"-t", "MOTHER,FATHER,CHILD"}},
		{"unrelated", unrel, []string{"-u", "S1,S2"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			dir := t.TempDir()
			upPrefix := filepath.Join(dir, "up")
			ourPrefix := filepath.Join(dir, "our")

			// color-chrs is a generic init/process plugin: options follow `--`.
			upArgs := append([]string{"+color-chrs", tc.fixture, "--", "-p", upPrefix}, tc.args...)
			ourArgs := append([]string{"+color-chrs", tc.fixture, "--", "-p", ourPrefix}, tc.args...)

			runColorChrs(t, bin, upArgs, pluginDirAbs(t))
			if ourBinPath == "" {
				t.Fatalf("local bcftools port binary not built")
			}
			runColorChrs(t, ourBinPath, ourArgs, pluginDirAbs(t))

			upDat, err := os.ReadFile(upPrefix + ".dat")
			if err != nil {
				t.Fatalf("read upstream .dat: %v", err)
			}
			ourDat, err := os.ReadFile(ourPrefix + ".dat")
			if err != nil {
				t.Fatalf("read ours .dat: %v", err)
			}
			if !bytes.Equal(upDat, ourDat) {
				t.Fatalf("+color-chrs %v .dat diverges from upstream\n--- upstream (%d bytes) ---\n%s\n--- ours (%d bytes) ---\n%s",
					tc.args, len(upDat), upDat, len(ourDat), ourDat)
			}
		})
	}
}

// runColorChrs runs a bcftools binary (upstream or ours) with the color-chrs
// argv and the vendored plugin dir on the environment. color-chrs writes its
// output to a file, so stdout is ignored.
func runColorChrs(t *testing.T, bin string, argv []string, pluginDir string) {
	t.Helper()
	cmd := exec.Command(bin, argv...)
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDir)
	var errBuf bytes.Buffer
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("%s %v: %v\nstderr: %s", bin, argv, err, errBuf.String())
	}
}

// TestNativePluginTrioDNM3Naive checks the trio-dnm3 NAIVE model's FORMAT/DNM
// and FORMAT/VA annotations across the PFM and PED trio sources, the
// strictly-novel flag, custom tags, chrX (male/female proband, GRCh38 list),
// multi-allelic sites, and the per-trio -i/-e filters, byte-for-byte.
func TestNativePluginTrioDNM3Naive(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "trio.vcf")
	ped := parityFixture(t, "trio.ped")
	multi := parityFixture(t, "trio_multi.vcf")
	multiPed := parityFixture(t, "trio_multi.ped")
	chrX := parityFixture(t, "trio_dnm_x.vcf")
	ma := parityFixture(t, "trio_dnm_ma.vcf")
	cases := []struct {
		fixture string
		args    []string
	}{
		{trio, []string{"--use-NAIVE", "-p", "CHILD,FATHER,MOTHER"}},
		{trio, []string{"--dnm-tag", "DNM:flag", "-p", "CHILD,FATHER,MOTHER"}},
		{trio, []string{"--use-NAIVE", "-P", ped}},
		{trio, []string{"--use-NAIVE", "-n", "-p", "CHILD,FATHER,MOTHER"}},
		{trio, []string{"--use-NAIVE", "--va", "DA", "-p", "CHILD,FATHER,MOTHER"}},
		{trio, []string{"--dnm-tag", "MYDNM:flag", "-p", "CHILD,FATHER,MOTHER"}},
		{multi, []string{"--use-NAIVE", "-P", multiPed}},
		{chrX, []string{"--use-NAIVE", "-p", "1X:CHILD,FATHER,MOTHER"}},
		{chrX, []string{"--use-NAIVE", "-p", "2X:CHILD,FATHER,MOTHER"}},
		{chrX, []string{"--use-NAIVE", "-p", "CHILD,FATHER,MOTHER"}},
		{chrX, []string{"--use-NAIVE", "-X", "GRCh38", "-p", "1X:CHILD,FATHER,MOTHER"}},
		{ma, []string{"--use-NAIVE", "-p", "CHILD,FATHER,MOTHER"}},
		{trio, []string{"--use-NAIVE", "-p", "CHILD,FATHER,MOTHER", "-i", "QUAL>10"}},
		{trio, []string{"--use-NAIVE", "-p", "CHILD,FATHER,MOTHER", "-e", "FMT/GQ<30"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "trio-dnm3", tc.args...)
		})
	}
}
