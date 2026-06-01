package mosdepth

// Live-binary oracle tests.
//
// Each test in this file invokes the genuine upstream mosdepth binary
// vendored under reference_code/mosdepth/mosdepth (a statically-linked
// Linux build of mosdepth 0.3.14) AND the local Go port (built once in
// TestMain into a temp dir) on the same fixture, then asserts
// byte-equality of the emitted output files. Intentional structural
// deltas (the upstream emits .csi, we emit .tbi; D4 unimplemented) are
// either compared via the underlying decoded BED or t.Skip'd with a
// pointer to docs/PARITY_ROADMAP.md.
//
// Divergences are surfaced via t.Errorf so they remain visible in CI.
// In-slice fixes for small divergences live in their own commits; large
// ones leave the t.Errorf hot and add a roadmap entry.
//
// The whole file t.Skip's when either binary is unavailable, so CI
// without the vendored upstream binary still passes.

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// liveOurBin is set by TestMain to the path of the locally-built mosdepth
// port binary. Empty when the build failed.
var liveOurBin string

// TestMain builds the local mosdepth port binary into a per-suite temp
// dir exactly once, so each TestLive* case can shell out to it cheaply.
// If the build fails, liveOurBin stays "" and every live-oracle test
// will t.Skip.
func TestMain(m *testing.M) {
	tmp, err := os.MkdirTemp("", "mosdepth-live-oracle-*")
	if err != nil {
		fmt.Fprintln(os.Stderr, "live-oracle: failed to make tempdir:", err)
		os.Exit(m.Run())
	}
	bin := filepath.Join(tmp, "mosdepth")
	cmd := exec.Command("go", "build", "-o", bin, "../../cmd/mosdepth")
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		fmt.Fprintln(os.Stderr, "live-oracle: go build failed:", err)
	} else {
		liveOurBin = bin
	}
	code := m.Run()
	_ = os.RemoveAll(tmp)
	os.Exit(code)
}

// liveUpstreamBin returns the absolute path to the vendored upstream
// mosdepth binary, or "" if it cannot be located / is not executable.
func liveUpstreamBin(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "..", "..",
		"reference_code", "mosdepth", "mosdepth"))
	if err != nil {
		return ""
	}
	fi, err := os.Stat(abs)
	if err != nil || fi.IsDir() || fi.Mode()&0111 == 0 {
		return ""
	}
	return abs
}

// requireLive skips the test when either binary is unavailable, and
// returns (upstream-path, our-port-path) otherwise.
func requireLive(t *testing.T) (live, ours string) {
	t.Helper()
	live = liveUpstreamBin(t)
	if live == "" {
		t.Skip("upstream mosdepth binary not found at reference_code/mosdepth/mosdepth; skipping live oracle")
	}
	if liveOurBin == "" {
		t.Skip("local mosdepth port binary not built; skipping live oracle")
	}
	return live, liveOurBin
}

// runMosdepth runs the binary with the given args and per-test prefix.
// stdout/stderr are captured and discarded unless the run fails (then
// stderr is surfaced into the test log). The args slice is prepended
// with the prefix and BAM by callers; this helper assumes args is the
// full argv.
func runMosdepth(t *testing.T, bin string, env []string, args ...string) {
	t.Helper()
	cmd := exec.Command(bin, args...)
	if env != nil {
		// Always inherit the parent env so PATH etc. still work.
		cmd.Env = append(os.Environ(), env...)
	}
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("mosdepth (%s) %v failed: %v\nstderr:\n%s",
			filepath.Base(bin), args, err, stderr.String())
	}
}

// readMaybeGz reads a file, transparently decompressing it if the path
// ends in .gz. Returns the raw (decoded) bytes so callers can compare
// the underlying BED records — BGZF framing differs across writers (we
// emit .tbi-tabix-compatible BGZF, upstream emits htslib BGZF) but the
// decoded bytes must match.
func readMaybeGz(t *testing.T, path string) []byte {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	if !strings.HasSuffix(path, ".gz") {
		return data
	}
	gr, err := gzip.NewReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("gzip.NewReader %s: %v", path, err)
	}
	gr.Multistream(true)
	out, err := io.ReadAll(gr)
	if err != nil {
		t.Fatalf("gunzip %s: %v", path, err)
	}
	return out
}

// diffBytes returns a short, human-readable diff of two byte buffers.
// When inputs are textual it shows the first differing line. For binary
// inputs it shows the first differing offset.
func diffBytes(a, b []byte) string {
	if bytes.Equal(a, b) {
		return ""
	}
	// Find first differing offset.
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	off := n
	for i := 0; i < n; i++ {
		if a[i] != b[i] {
			off = i
			break
		}
	}
	// Try line-based diff for textual content.
	if isMostlyText(a) && isMostlyText(b) {
		al := strings.Split(string(a), "\n")
		bl := strings.Split(string(b), "\n")
		var buf strings.Builder
		max := len(al)
		if len(bl) > max {
			max = len(bl)
		}
		shown := 0
		for i := 0; i < max && shown < 8; i++ {
			var av, bv string
			if i < len(al) {
				av = al[i]
			}
			if i < len(bl) {
				bv = bl[i]
			}
			if av != bv {
				fmt.Fprintf(&buf, "line %d:\n  want: %q\n  got:  %q\n", i+1, av, bv)
				shown++
			}
		}
		return fmt.Sprintf("len(want)=%d len(got)=%d first byte differs at offset %d\n%s",
			len(a), len(b), off, buf.String())
	}
	return fmt.Sprintf("len(want)=%d len(got)=%d first byte differs at offset %d", len(a), len(b), off)
}

// isMostlyText returns true if the slice is mostly printable ASCII /
// UTF-8, so diffBytes can switch to line-based output.
func isMostlyText(b []byte) bool {
	if len(b) == 0 {
		return true
	}
	n := len(b)
	if n > 4096 {
		n = 4096
	}
	nonText := 0
	for i := 0; i < n; i++ {
		c := b[i]
		if c == '\n' || c == '\t' || c == '\r' {
			continue
		}
		if c < 0x20 || c == 0x7f {
			nonText++
		}
	}
	return nonText*10 < n
}

// liveFixtureDir is the directory of BAM / BED fixtures shared with
// parity_test.go. It is resolved the same way (..,..,testdata,parity).
func liveFixtureDir(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity"))
	if err != nil {
		t.Fatalf("Abs: %v", err)
	}
	return abs
}

// ensureBAI generates a `.bai` next to bamPath if missing. The upstream
// mosdepth binary requires an index; our port doesn't.
func ensureBAI(t *testing.T, bamPath string) {
	t.Helper()
	if _, err := os.Stat(bamPath + ".bai"); err == nil {
		return
	}
	// Best-effort: try the vendored samtools binary.
	st, err := filepath.Abs(filepath.Join("..", "..", "..", "..",
		"reference_code", "samtools", "samtools"))
	if err != nil {
		t.Skipf("cannot resolve samtools path: %v", err)
	}
	if fi, err := os.Stat(st); err != nil || fi.IsDir() {
		t.Skipf("samtools binary not found at %s; cannot create BAI for %s", st, bamPath)
	}
	cmd := exec.Command(st, "index", bamPath)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Skipf("samtools index %s failed: %v\n%s", bamPath, err, string(out))
	}
}

// oracleCompareFile compares one named output between live and ours,
// applying gzip decode when appropriate. Reports the divergence via
// t.Errorf (not t.Fatalf) so subsequent file comparisons in the same
// subtest still run.
func oracleCompareFile(t *testing.T, name, livePath, oursPath string) {
	t.Helper()
	_, liveErr := os.Stat(livePath)
	_, oursErr := os.Stat(oursPath)
	if liveErr != nil && oursErr != nil {
		// Neither side produced the file -- nothing to do.
		return
	}
	if liveErr != nil {
		t.Errorf("DIVERGENCE: %s: live did not produce %s but ours produced %s", name, livePath, oursPath)
		return
	}
	if oursErr != nil {
		t.Errorf("DIVERGENCE: %s: ours did not produce %s but live produced %s", name, oursPath, livePath)
		return
	}
	want := readMaybeGz(t, livePath)
	got := readMaybeGz(t, oursPath)
	if !bytes.Equal(want, got) {
		t.Errorf("DIVERGENCE in %s:\n%s", name, diffBytes(want, got))
	}
}

// liveRunPair drives both binaries with identical positional/flag args
// (plus an env-var slice that is applied identically to both runs) and
// returns the two output prefixes.
func liveRunPair(t *testing.T, bam string, env []string, flags []string) (livePrefix, oursPrefix string) {
	t.Helper()
	live, ours := requireLive(t)

	bamPath := filepath.Join(liveFixtureDir(t), bam)
	ensureBAI(t, bamPath)

	tmp := t.TempDir()
	livePrefix = filepath.Join(tmp, "live")
	oursPrefix = filepath.Join(tmp, "ours")

	liveArgs := append([]string{}, flags...)
	liveArgs = append(liveArgs, livePrefix, bamPath)
	runMosdepth(t, live, env, liveArgs...)

	oursArgs := append([]string{}, flags...)
	oursArgs = append(oursArgs, oursPrefix, bamPath)
	runMosdepth(t, ours, env, oursArgs...)
	return livePrefix, oursPrefix
}

// oracleCompareCommon compares the three always-emitted files (plus
// optional regions / quantized / thresholds when present on either side)
// between live and ours. Callers pass in the list of expected file
// suffixes to check; absent files are still compared so divergences in
// "which files were emitted" also surface.
func oracleCompareCommon(t *testing.T, livePrefix, oursPrefix string, suffixes ...string) {
	t.Helper()
	for _, sfx := range suffixes {
		oracleCompareFile(t, sfx, livePrefix+sfx, oursPrefix+sfx)
	}
}

// commonSuffixes returns the list of output suffixes we always compare.
// Index files (.csi / .tbi) are excluded — that delta is intentional
// and tracked separately in docs/PARITY_ROADMAP.md#mosdepth.
func commonSuffixes(extra ...string) []string {
	base := []string{
		".mosdepth.global.dist.txt",
		".mosdepth.summary.txt",
	}
	return append(base, extra...)
}

// ---- Per-mode subtests ----

// TestLive_Bare covers the bare invocation with no extra flags: per-base
// emission across every chromosome in the BAM header.
func TestLive_Bare(t *testing.T) {
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil, []string{})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".per-base.bed.gz")...)
}

// TestLive_NoPerBase covers `-n / --no-per-base`: per-base file omitted.
func TestLive_NoPerBase(t *testing.T) {
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil, []string{"-n"})
	oracleCompareCommon(t, livePrefix, oursPrefix, commonSuffixes()...)
}

// TestLive_Chrom covers `-c MT`: restrict to one chromosome.
func TestLive_Chrom(t *testing.T) {
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil, []string{"-c", "MT"})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".per-base.bed.gz")...)
}

// TestLive_ByWindow covers `-b 5000`: fixed-window regions.
func TestLive_ByWindow(t *testing.T) {
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil,
		[]string{"-c", "MT", "-n", "-b", "5000"})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".regions.bed.gz")...)
}

// TestLive_ByBED covers `--by <BED>`: per-region coverage on a BED file.
func TestLive_ByBED(t *testing.T) {
	bed := filepath.Join(liveFixtureDir(t), "track.bed")
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil,
		[]string{"-c", "MT", "-n", "--by", bed})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".regions.bed.gz")...)
}

// TestLive_Quantize_DefaultLabels covers `-q 0:1:5:150` with no env
// override. Upstream emits "cutoff[i]:cutoff[i+1]" range labels by
// default (e.g. "0:1", "1:5"); we emit the NO_COVERAGE / LOW_COVERAGE
// mnemonics. The divergence is surfaced via DIVERGENCE: in t.Errorf and
// tracked in docs/PARITY_ROADMAP.md#mosdepth.
func TestLive_Quantize_DefaultLabels(t *testing.T) {
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil,
		[]string{"-c", "MT", "--fast-mode", "-q", "0:1:5:150"})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".quantized.bed.gz")...)
}

// TestLive_Quantize_EnvOverride covers `-q 0:1:5:150` with the
// MOSDEPTH_Q{i} environment variables set to override default labels.
// Both binaries should honour the env vars and emit identical labels.
func TestLive_Quantize_EnvOverride(t *testing.T) {
	env := []string{
		"MOSDEPTH_Q0=NO_COVERAGE",
		"MOSDEPTH_Q1=LOW_COVERAGE",
		"MOSDEPTH_Q2=CALLABLE",
		"MOSDEPTH_Q3=HIGH_COVERAGE",
	}
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", env,
		[]string{"-c", "MT", "--fast-mode", "-q", "0:1:5:150"})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".quantized.bed.gz")...)
}

// TestLive_Thresholds covers `-T 1,5,10` with `--by <BED>`.
func TestLive_Thresholds(t *testing.T) {
	bed := filepath.Join(liveFixtureDir(t), "track.bed")
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil,
		[]string{"-c", "MT", "--fast-mode", "--by", bed, "-T", "0,1,2"})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".regions.bed.gz", ".thresholds.bed.gz")...)
}

// TestLive_MapQ covers `-Q 20`: minimum MAPQ filter.
func TestLive_MapQ(t *testing.T) {
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil,
		[]string{"-c", "MT", "--fast-mode", "-Q", "20"})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".per-base.bed.gz")...)
}

// TestLive_ExcludeFlag covers `-F 1796`: exclude-flag mask (the
// default; we set it explicitly to ensure the CLI parses correctly).
func TestLive_ExcludeFlag(t *testing.T) {
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil,
		[]string{"-c", "MT", "--fast-mode", "-F", "1796"})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".per-base.bed.gz")...)
}

// TestLive_IncludeFlag covers `-i 2`: require-flag mask (proper pair).
func TestLive_IncludeFlag(t *testing.T) {
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil,
		[]string{"-c", "MT", "--fast-mode", "-i", "2"})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".per-base.bed.gz")...)
}

// TestLive_FastMode covers `-x / --fast-mode`: skip CIGAR walking.
func TestLive_FastMode(t *testing.T) {
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil,
		[]string{"-c", "MT", "--fast-mode"})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".per-base.bed.gz")...)
}

// TestLive_FastModeDefaultOff covers running without --fast-mode, which
// engages overlap-pair detection.
func TestLive_FastModeDefaultOff(t *testing.T) {
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil,
		[]string{"-c", "MT"})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".per-base.bed.gz")...)
}

// TestLive_UseMedian covers `-m / --use-median`. Upstream emits median
// depth instead of mean for `--by` regions. Our port does not implement
// `--use-median` and the `-m` short alias is bound to --fragment-mode,
// so this test is skipped with a TODO.
func TestLive_UseMedian(t *testing.T) {
	t.Skip("not implemented: --use-median; TODO in docs/PARITY_ROADMAP.md#mosdepth")
}

// TestLive_MeanMAPQ covers `--mean-mapq`. Not implemented in our port.
func TestLive_MeanMAPQ(t *testing.T) {
	t.Skip("not implemented: --mean-mapq; TODO in docs/PARITY_ROADMAP.md#mosdepth")
}

// TestLive_D4 covers `-d / --d4`. Our port rejects --d4 with a clear
// error; upstream emits a `.per-base.d4` file. This is a documented
// non-goal for v1.
func TestLive_D4(t *testing.T) {
	t.Skip("not implemented: --d4 (ErrD4NotImplemented); see docs/PARITY_ROADMAP.md#mosdepth")
}

// TestLive_FragmentMode covers `-a / --fragment-mode`: paired-end
// fragment depth. Note our port binds `-m` to --fragment-mode whereas
// upstream binds `-a` — the long flag --fragment-mode is the portable
// form and what we pass here.
func TestLive_FragmentMode(t *testing.T) {
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil,
		[]string{"-c", "MT", "--fragment-mode"})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".per-base.bed.gz")...)
}

// TestLive_ReadGroups covers `-R / --read-groups`. The ovl.bam test
// fixture has a single RG, so naming it should be a no-op; naming a
// nonexistent RG should leave depth at zero across the chromosome.
func TestLive_ReadGroups(t *testing.T) {
	// Inspect the BAM header via the live binary to find an RG.
	// Failing that, just try "x" which won't match any RG.
	live, _ := requireLive(t)
	bamPath := filepath.Join(liveFixtureDir(t), "ovl.bam")
	ensureBAI(t, bamPath)

	// Run a quick "list RG" by reading the SAM header.
	// We don't have samtools view wired in here; use the upstream
	// mosdepth binary's behaviour: name a guaranteed-nonexistent RG.
	_ = live
	livePrefix, oursPrefix := liveRunPair(t, "ovl.bam", nil,
		[]string{"-c", "MT", "--fast-mode", "--read-groups", "no-such-rg-xyz"})
	oracleCompareCommon(t, livePrefix, oursPrefix,
		commonSuffixes(".per-base.bed.gz")...)
}
