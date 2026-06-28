package bcftools

import (
	"bytes"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
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
//     two binaries into separate prefixes and diffs the produced .dat files. The
//     SW lines' nSwitches/rate columns are the one exception: upstream inflates
//     them via an out-of-bounds hap_switch[state][-1] read on each chromosome's
//     first segment (undefined behaviour, build-dependent garbage), so we assert
//     ours' correct deterministic counts against a per-fixture golden there
//     rather than byte-matching the buggy upstream binary. Everything else on
//     the SW line (SW/sample/chrom/nHets) and all SG/header lines stay byte-for-
//     byte against upstream. See docs/UPSTREAM_BUGS.md#bcftools-color-chrs-oob-switch.
//
//   - trio-dnm3 NAIVE model: byte parity on the FORMAT/DNM and FORMAT/VA
//     annotations. The NAIVE verdict is a pure integer Mendelian-consistency
//     table lookup (no floating point), so it is bit-reproducible.
//
//   - trio-dnm3 DMM/ALM/DNG float models: proximity parity (string fields exact,
//     DNM/VA/VAF numbers within the documented libm tolerance) in
//     TestNativePluginTrioDNM3FloatModels. The de-novo log score is a long
//     log/exp/pow/lgamma reduction; on linux/amd64 these inputs land byte-for-byte
//     because the incomplete-beta / lgamma kernels go through the bit-stable
//     in-tree kfunc port, but the proximity bar is the contract on any platform.
//
// All comparisons are CLI-to-CLI against the live upstream binary with
// provenance stripped; no committed goldens.

// TestNativePluginParentalOrigin checks the parental-origin summary and DBG
// lines across the dup/del CNV types, the -d debug listing, -g greedy, the -b
// min-binom threshold, and the -i/-e per-trio filters, all byte-for-byte.
func TestNativePluginParentalOrigin(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
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
		t.Skipf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "color_chrs_trio.vcf")
	unrel := parityFixture(t, "color_chrs_unrel.vcf")
	// swGolden maps "sample\tchrom" -> the correct, deterministic
	// "nSwitches\trate" pair our port emits on the SW line. These columns are NOT
	// byte-checked against upstream: upstream inflates them on each chromosome's
	// first segment via an out-of-bounds hap_switch[state][-1] read (undefined
	// behaviour; the garbage value is build/layout-dependent, here 3 = SW_MOTHER|
	// SW_FATHER, so +1 mother and +1 father per chromosome). Ours counts only the
	// real phase switches in the decoded segment path. See
	// docs/UPSTREAM_BUGS.md#bcftools-color-chrs-oob-switch. An empty golden means
	// "no SW lines carry UB-tainted columns" (unrelated mode never reads
	// hap_switch), so every SW line stays fully byte-checked against upstream.
	cases := []struct {
		name     string
		fixture  string
		args     []string          // plugin options, WITHOUT the -p prefix (added per-run)
		swGolden map[string]string // "sample\tchrom" -> "nSwitches\trate"
	}{
		// chr1: father :1->:2->:1 = 2 real switches; mother stays :1 = 0.
		// chr2: a single segment = 0 switches for both. rate = nSwitches/(nHets-1).
		{"trio", trio, []string{"-t", "MOTHER,FATHER,CHILD"}, map[string]string{
			"MOTHER\tchr1": "0\t0.000000",
			"FATHER\tchr1": "2\t0.666667",
			"MOTHER\tchr2": "0\t0.000000",
			"FATHER\tchr2": "0\t0.000000",
		}},
		{"trio_grch38", trio, []string{"-t", "MOTHER,FATHER,CHILD"}, map[string]string{
			"MOTHER\tchr1": "0\t0.000000",
			"FATHER\tchr1": "2\t0.666667",
			"MOTHER\tchr2": "0\t0.000000",
			"FATHER\tchr2": "0\t0.000000",
		}},
		// Unrelated mode never touches hap_switch, so its SW lines (all-zero here)
		// are already byte-identical to upstream: no golden override needed.
		{"unrelated", unrel, []string{"-u", "S1,S2"}, nil},
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
			compareColorChrsDat(t, tc.args, upDat, ourDat, tc.swGolden)
		})
	}
}

// compareColorChrsDat asserts ours' .dat output matches upstream's, with the SW
// line's UB-tainted columns handled specially. For every line that is not an SW
// data line (SG segments, "# ..." headers, and any blank trailing line) we
// require byte-equality ours==upstream — those are fully deterministic and the
// oracle must hold. For SW data lines ("SW\t<sample>\t<chrom>\t<nHets>\t
// <nSwitches>\t<rate>") we require ours==upstream on the SW/sample/chrom/nHets
// fields (also deterministic) but check the nSwitches and switch-rate fields
// against the supplied per-fixture golden instead, because upstream inflates
// them via an out-of-bounds hap_switch read on each chromosome's first segment
// (see docs/UPSTREAM_BUGS.md#bcftools-color-chrs-oob-switch). If swGolden has no
// entry for a given "sample\tchrom" key, that SW line is byte-checked against
// upstream in full (used by unrelated mode, which never reads hap_switch).
func compareColorChrsDat(t *testing.T, args []string, upDat, ourDat []byte, swGolden map[string]string) {
	t.Helper()
	upLines := strings.Split(string(upDat), "\n")
	ourLines := strings.Split(string(ourDat), "\n")
	if len(upLines) != len(ourLines) {
		t.Fatalf("+color-chrs %v .dat line count diverges: upstream %d vs ours %d\n--- upstream ---\n%s\n--- ours ---\n%s",
			args, len(upLines), len(ourLines), upDat, ourDat)
	}
	for i := range upLines {
		up, our := upLines[i], ourLines[i]
		if !strings.HasPrefix(up, "SW\t") {
			// SG segments, "# ..." headers, trailing blank line: pure oracle.
			if up != our {
				t.Fatalf("+color-chrs %v .dat line %d diverges from upstream\n--- upstream ---\n%q\n--- ours ---\n%q",
					args, i+1, up, our)
			}
			continue
		}
		// SW line layout: SW \t sample \t chrom \t nHets \t nSwitches \t rate
		upF := strings.Split(up, "\t")
		ourF := strings.Split(our, "\t")
		if len(upF) != 6 || len(ourF) != 6 {
			t.Fatalf("+color-chrs %v malformed SW line %d: upstream %q ours %q", args, i+1, up, our)
		}
		// Deterministic prefix (SW/sample/chrom/nHets) must byte-match upstream.
		for f := 0; f < 4; f++ {
			if upF[f] != ourF[f] {
				t.Fatalf("+color-chrs %v SW line %d field %d (deterministic) diverges: upstream %q ours %q",
					args, i+1, f, upF[f], ourF[f])
			}
		}
		key := upF[1] + "\t" + upF[2] // "sample\tchrom"
		want, ok := swGolden[key]
		if !ok {
			// No golden override: this SW line carries no UB-tainted columns, so
			// the nSwitches/rate columns are byte-checked against upstream too.
			if up != our {
				t.Fatalf("+color-chrs %v SW line %d diverges from upstream\n--- upstream ---\n%q\n--- ours ---\n%q",
					args, i+1, up, our)
			}
			continue
		}
		// UB-tainted columns: assert ours equals the recorded correct golden
		// rather than the buggy upstream value.
		got := ourF[4] + "\t" + ourF[5]
		if got != want {
			t.Fatalf("+color-chrs %v SW line %d (%s) nSwitches/rate = %q, want golden %q "+
				"(upstream is %q/%q — UB-inflated, see docs/UPSTREAM_BUGS.md#bcftools-color-chrs-oob-switch)",
				args, i+1, key, got, want, upF[4], upF[5])
		}
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
		t.Skipf("build upstream bcftools: %v", err)
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

// TestNativePluginTrioDNM3FloatModels checks the DMM/ALM/DNG float models'
// FORMAT/DNM (phred/log/prob), FORMAT/VA and FORMAT/VAF against upstream using
// the tolerance-aware proximity comparison (numeric_parity_test.go): string
// fields must match exactly, while the de-novo score is allowed the last-ULP
// libm slack documented there. The de-novo log score is a long log/exp/pow/
// lgamma reduction; on linux/amd64 these inputs actually land byte-for-byte (the
// incomplete-beta / lgamma kernels go through the bit-stable in-tree kfunc port),
// but the proximity bar is the correct contract for the float models on any
// platform. The test reports the maximum observed per-model deviation.
func TestNativePluginTrioDNM3FloatModels(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	trio := parityFixture(t, "trio.vcf")         // AD+PL only
	fl := parityFixture(t, "trio_dnm_float.vcf") // AD+PL+QS+QM+SP, incl. a 5-allele site
	multi := parityFixture(t, "trio_multi.vcf")  // multi-trio AD+PL
	cases := []struct {
		model   string
		fixture string
		args    []string
	}{
		// DNG (PL only) works on every fixture.
		{"DNG", trio, []string{"--use-DNG", "-p", "CHILD,FATHER,MOTHER"}},
		{"DNG", trio, []string{"--use-DNG", "-n", "-p", "CHILD,FATHER,MOTHER"}},
		{"DNG", trio, []string{"--use-DNG", "--dnm-tag", "DNM:phred", "-p", "CHILD,FATHER,MOTHER"}},
		{"DNG", trio, []string{"--use-DNG", "--dnm-tag", "DNM:prob", "-p", "CHILD,FATHER,MOTHER"}},
		{"DNG", trio, []string{"--use-DNG", "--dng-priors", "-p", "CHILD,FATHER,MOTHER"}},
		{"DNG", trio, []string{"--use-DNG", "--va", "DA", "--vaf", "MYVAF", "-p", "CHILD,FATHER,MOTHER"}},
		{"DNG", multi, []string{"--use-DNG", "-p", "K1,F1,M1"}},
		{"DNG", fl, []string{"--use-DNG", "--strand-bias", "0.05", "-p", "CHILD,FATHER,MOTHER"}},
		// ALM over fake-QS-from-AD (--with-pAD) and over PL (--with-pPL) on trio.vcf,
		// and over real FORMAT/QS on the rich fixture.
		{"ALM", trio, []string{"--use-ALM", "--with-pAD", "-p", "CHILD,FATHER,MOTHER"}},
		{"ALM", trio, []string{"--use-ALM", "--with-pPL", "-p", "CHILD,FATHER,MOTHER"}},
		{"ALM", trio, []string{"--use-ALM", "--with-pAD", "-n", "-p", "CHILD,FATHER,MOTHER"}},
		{"ALM", trio, []string{"--use-ALM", "--with-pAD", "--ad", "0.5", "-p", "CHILD,FATHER,MOTHER"}},
		{"ALM", fl, []string{"--use-ALM", "-p", "CHILD,FATHER,MOTHER"}},
		{"ALM", fl, []string{"--use-ALM", "--strand-bias", "0.05", "-p", "CHILD,FATHER,MOTHER"}},
		// DMM needs FORMAT/AD+QM (the rich fixture) or --max-QM negative on AD+PL.
		{"DMM", trio, []string{"--use-DMM", "--max-QM", "-1", "-p", "CHILD,FATHER,MOTHER"}},
		{"DMM", fl, []string{"--use-DMM", "-p", "CHILD,FATHER,MOTHER"}},
		{"DMM", fl, []string{"--use-DMM", "--with-cAD", "-p", "CHILD,FATHER,MOTHER"}},
		{"DMM", fl, []string{"--use-DMM", "--min-vaf", "0.3", "-p", "CHILD,FATHER,MOTHER"}},
		{"DMM", fl, []string{"--use-DMM", "--strand-bias", "0.05", "-p", "CHILD,FATHER,MOTHER"}},
		{"DMM", fl, []string{"--use-DMM", "--pn", "0.02,1:snv", "--pns", "0.05,2:snv", "-p", "CHILD,FATHER,MOTHER"}},
		{"DMM", fl, []string{"--use-DMM", "--dnm-tag", "DNM:phred", "-p", "CHILD,FATHER,MOTHER"}},
	}
	maxDev := map[string]float64{"DNG": 0, "ALM": 0, "DMM": 0}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.model+"_"+shortName(tc.fixture)+"_"+joinArgs(tc.args), func(t *testing.T) {
			argv := pluginCLIArgs("trio-dnm3", tc.fixture, tc.args)
			up := string(stripProvenanceBytes(runUpstreamPlugin(t, bin, argv...)))
			ours := string(stripProvenanceBytes(runOursPlugin(t, argv...)))
			if diffs := compareProximityDefault(up, ours); diffs != nil {
				var b strings.Builder
				for _, d := range diffs {
					b.WriteString("  " + d.String() + "\n")
				}
				t.Fatalf("+trio-dnm3 %v diverges beyond tolerance:\n%s", tc.args, b.String())
			}
			if d := maxFieldDeviation(up, ours); d > maxDev[tc.model] {
				maxDev[tc.model] = d
			}
		})
	}
	t.Logf("max observed per-model DNM/VA/VAF field deviation (within tolerance): DNG=%g ALM=%g DMM=%g",
		maxDev["DNG"], maxDev["ALM"], maxDev["DMM"])
}

// maxFieldDeviation returns the largest absolute difference between any pair of
// numeric fields in two proximity-equal texts, for reporting how close the
// libm-bound scores actually land.
func maxFieldDeviation(want, got string) float64 {
	wl := strings.Split(want, "\n")
	gl := strings.Split(got, "\n")
	maxd := 0.0
	n := len(wl)
	if len(gl) < n {
		n = len(gl)
	}
	for i := 0; i < n; i++ {
		wf := strings.FieldsFunc(wl[i], func(r rune) bool { return r == '\t' || r == ' ' || r == ':' || r == ',' })
		gf := strings.FieldsFunc(gl[i], func(r rune) bool { return r == '\t' || r == ' ' || r == ':' || r == ',' })
		m := len(wf)
		if len(gf) < m {
			m = len(gf)
		}
		for j := 0; j < m; j++ {
			a, aok := parseNumericField(wf[j])
			b, bok := parseNumericField(gf[j])
			if aok && bok && !math.IsInf(a, 0) && !math.IsInf(b, 0) && !math.IsNaN(a) && !math.IsNaN(b) {
				if d := math.Abs(a - b); d > maxd {
					maxd = d
				}
			}
		}
	}
	return maxd
}

// TestNativePluginTrioDNM3FloatErrors asserts the float models still report the
// clean fatal-data errors upstream does: DMM on an AD/QM-less input, and a
// non-positive --phi.
func TestNativePluginTrioDNM3FloatErrors(t *testing.T) {
	gtOnly := parityFixture(t, "trio_dnm_x.vcf") // GT only, no PL/AD
	cases := [][]string{
		{"--use-DMM", "-p", "CHILD,FATHER,MOTHER"}, // no FORMAT/AD/QM/PL
		{"--use-DNG", "-p", "CHILD,FATHER,MOTHER"}, // no FORMAT/PL
		{"--use-DMM", "--phi", "0", "-p", "CHILD,FATHER,MOTHER"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			var out, errBuf bytes.Buffer
			err := RunPlugin(PluginOptions{
				Name:         "trio-dnm3",
				Args:         args,
				InputFile:    gtOnly,
				OutputFormat: OutputVCF,
			}, &out, &errBuf)
			if err == nil {
				t.Fatalf("expected a clean error for +trio-dnm3 %v, got nil", args)
			}
		})
	}
}
