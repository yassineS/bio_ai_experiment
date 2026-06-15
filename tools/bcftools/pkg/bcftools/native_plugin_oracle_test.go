package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// Live-oracle parity tests for the native (in-process, pure-Go) plugins.
//
// Each test builds the genuine upstream bcftools binary from the vendored
// submodule via buildBcftools() (t.Fatalf — not t.Skip — when the build is
// impossible) and runs `bcftools +<name> <fixture> -- <args>` on BOTH the
// upstream binary (BCFTOOLS_PLUGINS pointed at the vendored .so directory) and
// the native Go pipeline (RunPlugin, which dispatches to the native registry).
// The two outputs are compared byte-for-byte after stripping the provenance
// lines (##bcftools_*, version banners) that legitimately differ between the
// two — reusing stripProvenanceBytes from the existing live-oracle suite.

// pluginDirAbs returns the absolute path of the vendored compiled-plugin
// directory that the upstream binary dlopen's via BCFTOOLS_PLUGINS.
func pluginDirAbs(t *testing.T) string {
	t.Helper()
	dir, err := filepath.Abs(filepath.Join("..", "..", "..", "..", "reference_code", "bcftools", "plugins"))
	if err != nil {
		t.Fatalf("plugin dir abs: %v", err)
	}
	return dir
}

// pluginCLIArgs builds the argv for a `+name` plugin invocation in the form
// that the REAL upstream bcftools accepts, so the very same argv can be handed
// to both the upstream binary and our port. There are two upstream forms:
//
//   - run()-style plugins (those whose upstream .so exports a `run` symbol,
//     e.g. variant-distance) take their own options BEFORE the input file with
//     no `--` separator: `+name <opts...> <file>`.
//   - the generic init/process plugins take host options before the file and
//     the plugin's own options after a `--`: `+name <file> -- <opts...>`.
//
// Mirroring upstream here is what makes the oracle a genuine CLI-to-CLI test:
// if our port rejects a form upstream accepts, the subprocess fails and the
// test fails (which is exactly how the variant-distance bug would have been
// caught before the fix).
func pluginCLIArgs(name, fixture string, args []string) []string {
	if IsRunStyleNativePlugin(name) {
		argv := append([]string{"+" + name}, args...)
		return append(argv, fixture)
	}
	argv := []string{"+" + name, fixture}
	if len(args) > 0 {
		argv = append(argv, "--")
		argv = append(argv, args...)
	}
	return argv
}

// runUpstreamPlugin runs the upstream bcftools binary with the vendored
// plugins directory on PATH and returns its stdout.
func runUpstreamPlugin(t *testing.T, bin string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
	var out bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = nil
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream %v: %v", args, err)
	}
	return out.Bytes()
}

// runOursPlugin runs OUR built bcftools binary (ourBinPath, produced by
// TestMain) through its actual CLI with the given argv, returning stdout. It
// drives the port the same way a user would, so CLI argument-routing bugs are
// visible to the oracle. A non-zero exit is a hard failure — e.g. our port
// rejecting an upstream-accepted command form.
func runOursPlugin(t *testing.T, args ...string) []byte {
	t.Helper()
	if ourBinPath == "" {
		t.Fatalf("local bcftools port binary not built; cannot run CLI oracle")
	}
	cmd := exec.Command(ourBinPath, args...)
	cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDirAbs(t))
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("ours %v: %v\nstderr: %s", args, err, errBuf.String())
	}
	return out.Bytes()
}

// assertPluginParity asserts that OUR port and upstream produce byte-identical
// stdout (modulo provenance) when BOTH are driven through their CLIs with the
// SAME upstream-accepted command form for `+name`. This is a true CLI-to-CLI
// oracle: it fails if our port rejects (or mis-routes) a form upstream accepts.
func assertPluginParity(t *testing.T, bin, fixture, name string, args ...string) {
	t.Helper()
	argv := pluginCLIArgs(name, fixture, args)
	up := runUpstreamPlugin(t, bin, argv...)
	ours := runOursPlugin(t, argv...)
	if !bytes.Equal(stripProvenanceBytes(up), stripProvenanceBytes(ours)) {
		t.Fatalf("+%s %v diverges from upstream (argv=%v)\n--- upstream (%d bytes) ---\n%s\n--- ours (%d bytes) ---\n%s",
			name, args, argv, len(up), snippet(up, 1200), len(ours), snippet(ours, 1200))
	}
}

func parityFixture(t *testing.T, name string) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("abs: %v", err)
	}
	if _, err := os.Stat(abs); err != nil {
		t.Fatalf("fixture %s: %v", name, err)
	}
	return abs
}

func TestNativePluginFillTags(t *testing.T) {
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
		{basic, []string{"-t", "all"}},
		{basic, []string{"-t", "AN,AC"}},
		{basic, []string{"-t", "AN,AC,AF,MAF,NS"}},
		{basic, []string{"-t", "AC_Hom,AC_Het,AC_Hemi"}},
		{basic, []string{"-t", "HWE,ExcHet"}},
		{basic, []string{"-t", "END,TYPE"}},
		{basic, []string{"-t", "VAF,VAF1"}},
		{basic, []string{"-t", "F_MISSING"}},
		{basic, []string{"-d", "-t", "all"}},
		{multi, []string{"-t", "all"}},
		{gt, []string{"-t", "all"}},
		{gt, []string{"-t", "VAF,VAF1"}},
	}
	for _, tc := range cases {
		t.Run(filepath.Base(tc.fixture)+" "+joinArgs(tc.args), func(t *testing.T) {
			assertPluginParity(t, bin, tc.fixture, "fill-tags", tc.args...)
		})
	}
}

func TestNativePluginFillANAC(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	for _, fx := range []string{"basic.vcf", "multi.vcf", "gt_plugins.vcf"} {
		fx := fx
		t.Run(fx, func(t *testing.T) {
			assertPluginParity(t, bin, parityFixture(t, fx), "fill-AN-AC")
		})
	}
}

func TestNativePluginMissing2Ref(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	for _, args := range [][]string{nil, {"-p"}, {"-m"}, {"-p", "-m"}} {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, gt, "missing2ref", args...)
		})
	}
}

func TestNativePluginSetGT(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	cases := [][]string{
		{"-t", ".", "-n", "0"},
		{"-t", ".", "-n", "0p"},
		{"-t", "./.", "-n", "0"},
		{"-t", "./x", "-n", "."},
		{"-t", "a", "-n", "M"},
		{"-t", "a", "-n", "Mp"},
		{"-t", "a", "-n", "m"},
		{"-t", "a", "-n", "p"},
		{"-t", "a", "-n", "u"},
		{"-t", "a", "-n", "i"},
		{"-t", "a", "-n", "c:0/0"},
		{"-t", ".", "-n", "c:0/0"},
		{"-t", "./x", "-n", "c:m/M"},
		// Filter-query mode (-t q with -i/-e): per-sample genotype setting.
		{"-t", "q", "-n", "0", "-i", `GT="het"`},
		{"-t", "q", "-n", ".", "-e", "FMT/DP>10"},
		{"-t", "q", "-n", "0", "-i", "FMT/GQ>30"},
		{"-t", "q", "-n", "M", "-i", `GT="alt"`},
		{"-t", "q", "-n", "0", "-e", `GT="mis"`},
		{"-t", "q", "-n", ".", "-i", "INFO/AC>1"},
		{"-t", "q", "-n", "0", "-i", `GT="het" && FMT/DP>15`},
		{"-t", "q", "-n", "0", "-e", `GT="het" && FMT/DP>15`},
		{"-t", "q", "-n", "c:0/0", "-i", `GT="alt"`},
		{"-t", "q", "-n", "0", "-e", "INFO/AC>1"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, gt, "setGT", args...)
		})
	}
}

func TestNativePluginTag2Tag(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gl := parityFixture(t, "gt_plugins.vcf")
	gtfx := parityFixture(t, "tag2tag_gt.vcf")
	glCases := [][]string{
		{"--PL-to-GL"},
		{"--PL-to-GP"},
		{"--PL-to-GL", "-r"},
	}
	for _, args := range glCases {
		args := args
		t.Run("gl "+joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, gl, "tag2tag", args...)
		})
	}
	gtCases := [][]string{
		{"--PL-to-GT"},
		{"--GP-to-GT"},
		{"--PL-to-GT", "-t", "0.3"},
		{"--PL-to-GT", "-r"},
		{"--PL-to-GP"},
		{"--PL-to-GL"},
	}
	for _, args := range gtCases {
		args := args
		t.Run("gt "+joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, gtfx, "tag2tag", args...)
		})
	}

	// GL-source conversions exercised against a GL-bearing fixture.
	glfx := parityFixture(t, "tag2tag_gl.vcf")
	for _, args := range [][]string{{"--GL-to-PL"}, {"--GL-to-GP"}, {"--GL-to-PL", "-r"}} {
		args := args
		t.Run("glsrc "+joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, glfx, "tag2tag", args...)
		})
	}
}

func TestNativePluginFixPloidy(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	gt := parityFixture(t, "gt_plugins.vcf")
	xy := parityFixture(t, "fixploidy_xy.vcf")
	for _, args := range [][]string{nil, {"-f", "1"}, {"-f", "2"}, {"-f", "3"}, {"-d", "1"}} {
		args := args
		t.Run("gt "+joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, gt, "fixploidy", args...)
		})
	}
	// X/Y/MT default-table and explicit -p/-s region+sex paths.
	t.Run("xy default", func(t *testing.T) {
		assertPluginParity(t, bin, xy, "fixploidy")
	})

	sexFile := filepath.Join(t.TempDir(), "sex.txt")
	if err := os.WriteFile(sexFile, []byte("male1 M\nmale2 M\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	t.Run("xy sex", func(t *testing.T) {
		assertPluginParity(t, bin, xy, "fixploidy", "-s", sexFile)
	})

	ploidyFile := filepath.Join(t.TempDir(), "ploidy.txt")
	pl := "X 1 60000 M 1\nX 2699521 154931043 M 1\nY 1 59373566 M 1\nY 1 59373566 F 0\nMT 1 16569 M 1\nMT 1 16569 F 1\n"
	if err := os.WriteFile(ploidyFile, []byte(pl), 0o644); err != nil {
		t.Fatal(err)
	}
	t.Run("xy ploidy+sex", func(t *testing.T) {
		assertPluginParity(t, bin, xy, "fixploidy", "-p", ploidyFile, "-s", sexFile)
	})
}

func TestNativePluginCounts(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Fatalf("build upstream bcftools: %v", err)
	}
	for _, fx := range []string{"basic.vcf", "multi.vcf", "gt_plugins.vcf"} {
		fx := fx
		t.Run(fx, func(t *testing.T) {
			fixture := parityFixture(t, fx)
			// counts is a generic plugin taking no plugin options: `+counts FILE`.
			argv := pluginCLIArgs("counts", fixture, nil)
			up := runUpstreamPlugin(t, bin, argv...)
			ours := runOursPlugin(t, argv...)
			if !bytes.Equal(up, ours) {
				t.Fatalf("+counts diverges from upstream (argv=%v)\n--- upstream ---\n%s\n--- ours ---\n%s", argv, up, ours)
			}
		})
	}
}

// TestNativePluginParallelDeterminism asserts that running a per-record native
// plugin with multiple worker threads produces byte-identical output to the
// single-threaded run (the ordering guarantee of processRecords).
func TestNativePluginParallelDeterminism(t *testing.T) {
	fixture, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", "basic.vcf"))
	if err != nil {
		t.Fatal(err)
	}
	var single bytes.Buffer
	if err := RunPlugin(PluginOptions{Name: "fill-tags", Args: []string{"-t", "all"}, InputFile: fixture, OutputFormat: OutputVCF, Threads: 1}, &single, &bytes.Buffer{}); err != nil {
		t.Fatalf("single-threaded: %v", err)
	}
	for _, nthreads := range []int{2, 4, 8} {
		var multi bytes.Buffer
		if err := RunPlugin(PluginOptions{Name: "fill-tags", Args: []string{"-t", "all"}, InputFile: fixture, OutputFormat: OutputVCF, Threads: nthreads}, &multi, &bytes.Buffer{}); err != nil {
			t.Fatalf("threads=%d: %v", nthreads, err)
		}
		if !bytes.Equal(single.Bytes(), multi.Bytes()) {
			t.Fatalf("threads=%d output differs from single-threaded run", nthreads)
		}
	}
}

// joinArgs renders args for a subtest name, using "default" for the empty set.
func joinArgs(args []string) string {
	if len(args) == 0 {
		return "default"
	}
	out := ""
	for i, a := range args {
		if i > 0 {
			out += "_"
		}
		out += a
	}
	return out
}
