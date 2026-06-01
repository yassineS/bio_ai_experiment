package bcftools

// Live-binary oracle tests.
//
// Every test in this file invokes the genuine upstream bcftools binary
// vendored under reference_code/bcftools/bcftools AND the local Go port
// (built once in TestMain into a temp dir) on the same fixture, then
// asserts byte-equality of stdout (modulo provenance lines and a few
// well-documented stderr-progress quirks).
//
// The motivation is concrete: the existing parity tests compare against
// stale vendored expected files. This file closes the loop so any
// future divergence from real upstream behaviour fails CI immediately,
// rather than silently passing because both sides regressed in lock-step.
//
// Subcommand coverage matrix is documented in docs/PARITY_ROADMAP.md
// under "Live oracle coverage". When upstream cannot be located the
// suite t.Skip's so CI without submodules still passes.

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strings"
	"testing"
)

// liveBinPath returns the absolute path to the vendored upstream bcftools
// binary, or "" when it cannot be located or is not executable. Tests
// using it should t.Skip on empty.
func liveBinPath(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "..", "..",
		"reference_code", "bcftools", "bcftools"))
	if err != nil {
		return ""
	}
	fi, err := os.Stat(abs)
	if err != nil || fi.IsDir() || fi.Mode()&0111 == 0 {
		return ""
	}
	return abs
}

// ourBinPath is set by TestMain (live_oracle_main_test.go) to the path of
// the locally-built port binary. Empty when the build failed.
var ourBinPath string

// requireLive skips the test when either the upstream binary or our
// built binary is unavailable.
func requireLive(t *testing.T) (live, ours string) {
	t.Helper()
	live = liveBinPath(t)
	if live == "" {
		t.Skip("upstream bcftools binary not found; skipping live oracle")
	}
	if ourBinPath == "" {
		t.Skip("local bcftools port binary not built; skipping live oracle")
	}
	return live, ourBinPath
}

// runBin invokes a binary with the given args, returns stdout, fails on
// non-zero exit. Stderr is suppressed (we only oracle-compare stdout).
func runBin(t *testing.T, bin string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var out bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = io.Discard
	if err := cmd.Run(); err != nil {
		t.Fatalf("%s %v: %v", bin, args, err)
	}
	return out.Bytes()
}

// livePluginDir returns the absolute path to the vendored compiled
// plugins directory (the .so files that the upstream binary dlopen's),
// or "" when it is unavailable. The upstream binary needs
// BCFTOOLS_PLUGINS pointed here to run `+fill-tags` / `+mendelian2`.
func livePluginDir(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "..", "..",
		"reference_code", "bcftools", "plugins"))
	if err != nil {
		return ""
	}
	if fi, err := os.Stat(abs); err != nil || !fi.IsDir() {
		return ""
	}
	return abs
}

// runBinEnv is runBin with extra environment entries (e.g.
// BCFTOOLS_PLUGINS) appended to the inherited environment.
func runBinEnv(t *testing.T, bin string, env []string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	cmd.Env = append(os.Environ(), env...)
	var out bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = io.Discard
	if err := cmd.Run(); err != nil {
		t.Fatalf("%s %v: %v", bin, args, err)
	}
	return out.Bytes()
}

// Provenance lines that legitimately differ between upstream and our port
// (they encode version strings and absolute paths). Stripping them keeps
// the oracle focused on actual semantic divergence.
var (
	// ##bcftools_<name>=... and ##bcftools/<name>=... — both forms are
	// emitted by upstream depending on the subcommand / plugin path.
	provVersionRE   = regexp.MustCompile(`(?m)^##bcftools[_/][^=]+=.*\n`)
	provVersionRE2  = regexp.MustCompile(`(?m)^##bcftoolsVersion=.*\n`)
	provVersionRE3  = regexp.MustCompile(`(?m)^##bcftoolsCommand=.*\n`)
	provReferenceRE = regexp.MustCompile(`(?m)^##reference=.*\n`)
	// Comment-style provenance: roh / gtcheck / stats banners. Match
	// the permissive form `# This file was produced by` followed by
	// anything (including a `:` directly after `by`) so the regex
	// covers both `produced by bcftools` and `produced by: bcftools roh`.
	provStatsRE1   = regexp.MustCompile(`(?m)^# This file was produced by.*\n`)
	provStatsRE2   = regexp.MustCompile(`(?m)^# The command line was:.*\n`)
	provStatsRE3   = regexp.MustCompile(`(?m)^# and the working directory was:.*\n`)
	provStatsRE4   = regexp.MustCompile(`(?m)^# \t .*\n`)
	provBlankHash  = regexp.MustCompile(`(?m)^#\n`)
	provStatsTime  = regexp.MustCompile(`(?m)^INFO\tTime required.*\n`)
	provStatsTime2 = regexp.MustCompile(`(?m)^INFO\trun-time.*\n`)
)

// stripProvenance removes the upstream provenance noise that always
// differs between two invocations: version strings, command lines, and
// absolute reference paths.
func stripProvenance(b []byte) []byte {
	b = provVersionRE.ReplaceAll(b, nil)
	b = provVersionRE2.ReplaceAll(b, nil)
	b = provVersionRE3.ReplaceAll(b, nil)
	b = provReferenceRE.ReplaceAll(b, nil)
	b = provStatsRE1.ReplaceAll(b, nil)
	b = provStatsRE2.ReplaceAll(b, nil)
	b = provStatsRE3.ReplaceAll(b, nil)
	b = provStatsRE4.ReplaceAll(b, nil)
	b = provBlankHash.ReplaceAll(b, nil)
	b = provStatsTime.ReplaceAll(b, nil)
	b = provStatsTime2.ReplaceAll(b, nil)
	return b
}

// gunzipBytes BGZF/gzip-decodes b. Used to oracle-compare -Oz outputs
// without depending on the block-level encoding (BGZF block boundaries
// are an implementation detail).
func gunzipBytes(t *testing.T, b []byte) []byte {
	t.Helper()
	zr, err := gzip.NewReader(bytes.NewReader(b))
	if err != nil {
		t.Fatalf("gunzip: %v", err)
	}
	defer zr.Close()
	out, err := io.ReadAll(zr)
	if err != nil {
		t.Fatalf("gunzip read: %v", err)
	}
	return out
}

// fixturePath returns an absolute path to the named fixture under
// tools/bcftools/testdata/parity. Tests that need other testdata
// subdirs construct their own paths.
func fixturePath(t *testing.T, name string) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("abs: %v", err)
	}
	return abs
}

// assertEqualStdout runs two binaries with the same args on the same
// fixture and asserts stripProvenance-equal stdout. The case name is
// reported on failure so the matrix is readable.
func assertEqualStdout(t *testing.T, live, ours string, args ...string) {
	t.Helper()
	livOut := runBin(t, live, args...)
	ourOut := runBin(t, ours, args...)
	if !bytes.Equal(stripProvenance(livOut), stripProvenance(ourOut)) {
		t.Errorf("stdout diverges for args=%v\n--- live (%d bytes) ---\n%s\n--- ours (%d bytes) ---\n%s",
			args, len(livOut), snippet(livOut, 800),
			len(ourOut), snippet(ourOut, 800))
	}
}

// snippet truncates b to n bytes for log readability.
func snippet(b []byte, n int) string {
	if len(b) <= n {
		return string(b)
	}
	return string(b[:n]) + "\n... [truncated]"
}

// -------------------------------------------------------------------------
// VIEW — most thoroughly exercised subcommand.
// -------------------------------------------------------------------------

func TestLiveView(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")

	cases := [][]string{
		{"view", fx},
		{"view", "-v", "snps", fx},
		{"view", "-v", "indels", fx},
		{"view", "-v", "mnps", fx},
		{"view", "-V", "snps", fx},
		{"view", "-s", "S1", fx},
		{"view", "-s", "S1,S3", fx},
		{"view", "-c", "1", fx},
		{"view", "-f", ".,PASS", fx},
		{"view", "--no-update", fx},
		{"view", "-G", fx}, // drop genotypes
		{"view", "-h", fx},
		{"view", "-H", fx},
	}
	for _, args := range cases {
		key := strings.Join(args[1:len(args)-1], "_")
		if key == "" {
			key = "bare"
		}
		t.Run(key, func(t *testing.T) {
			assertEqualStdout(t, live, ours, args...)
		})
	}
}

// TestLiveViewOz checks BGZF-compressed view output: the inflated bytes
// must match. (Block boundaries are an implementation detail; we don't
// require identical BGZF framing.)
func TestLiveViewOz(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	livOut := runBin(t, live, "view", "-Oz", fx)
	ourOut := runBin(t, ours, "view", "-Oz", fx)
	livPlain := stripProvenance(gunzipBytes(t, livOut))
	ourPlain := stripProvenance(gunzipBytes(t, ourOut))
	if !bytes.Equal(livPlain, ourPlain) {
		t.Errorf("view -Oz inflated bytes diverge\n--- live ---\n%s\n--- ours ---\n%s",
			snippet(livPlain, 600), snippet(ourPlain, 600))
	}
}

// TestLiveViewOb round-trips through BCF and decodes via the live binary
// on both sides — that gives a stable view-VCF text for comparison
// without coupling to the exact BCF byte layout (which the port's
// bcfwriter and upstream's may differ on for valid-but-different
// reasons, like type-narrowing choices).
func TestLiveViewOb(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	tmp := t.TempDir()

	livBCF := filepath.Join(tmp, "live.bcf")
	ourBCF := filepath.Join(tmp, "our.bcf")
	if err := os.WriteFile(livBCF, runBin(t, live, "view", "-Ob", fx), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(ourBCF, runBin(t, ours, "view", "-Ob", fx), 0o644); err != nil {
		t.Fatal(err)
	}
	livText := stripProvenance(runBin(t, live, "view", livBCF))
	ourText := stripProvenance(runBin(t, live, "view", ourBCF))
	if !bytes.Equal(livText, ourText) {
		t.Errorf("view -Ob BCFs decode to different VCF text\n--- live ---\n%s\n--- ours ---\n%s",
			snippet(livText, 600), snippet(ourText, 600))
	}
}

// TestLiveViewFromGz reads a bgzipped VCF — covers the iohelper plumbing
// rather than a particular flag.
func TestLiveViewFromGz(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf.gz")
	assertEqualStdout(t, live, ours, "view", fx)
}

// TestLiveViewFromBcf reads a BCF directly — covers the bcfreader path.
func TestLiveViewFromBcf(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.bcf")
	assertEqualStdout(t, live, ours, "view", fx)
}

// TestLiveViewRegions exercises the indexed-fetch path via --regions on
// a freshly-indexed bgzipped fixture. The .tbi is generated by the live
// binary so we don't conflate two divergences.
func TestLiveViewRegions(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	tmp := t.TempDir()
	gz := filepath.Join(tmp, "basic.vcf.gz")
	if err := os.WriteFile(gz, runBin(t, live, "view", "-Oz", fx), 0o644); err != nil {
		t.Fatal(err)
	}
	runBin(t, live, "index", "-t", gz)
	assertEqualStdout(t, live, ours, "view", "--regions", "chr1", gz)
}

// -------------------------------------------------------------------------
// STATS
// -------------------------------------------------------------------------

func TestLiveStats(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	for _, args := range [][]string{
		{"stats", fx},
		{"stats", "-s", "-", fx},
		{"stats", "-f", "PASS", fx},
	} {
		t.Run(strings.Join(args[1:len(args)-1], "_"), func(t *testing.T) {
			assertEqualStdout(t, live, ours, args...)
		})
	}
}

// -------------------------------------------------------------------------
// QUERY
// -------------------------------------------------------------------------

func TestLiveQuery(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	cases := []struct {
		name string
		args []string
	}{
		{"list_samples", []string{"query", "-l", fx}},
		{"chrom_pos", []string{"query", "-f", `%CHROM\t%POS\n`, fx}},
		{"info_ac", []string{"query", "-f", `%INFO/AC\n`, fx}},
		{"sample_gt", []string{"query", "-f", `[%SAMPLE=%GT\n]`, fx}},
		{"include_ac_gt_1", []string{"query", "-i", "AC>1", "-f", `%POS\n`, fx}},
		{"include_type_snp", []string{"query", "-i", `TYPE="snp"`, "-f", `%POS\n`, fx}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			assertEqualStdout(t, live, ours, tc.args...)
		})
	}
}

// -------------------------------------------------------------------------
// ANNOTATE
// -------------------------------------------------------------------------

func TestLiveAnnotate(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	cases := []struct {
		name string
		args []string
	}{
		{"remove_info_af", []string{"annotate", "-x", "INFO/AF", fx}},
		{"remove_format_gq", []string{"annotate", "-x", "FORMAT/GQ", fx}},
		{"set_id", []string{"annotate", "-I", `+%CHROM\_%POS`, fx}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			assertEqualStdout(t, live, ours, tc.args...)
		})
	}
}

// -------------------------------------------------------------------------
// HEAD
// -------------------------------------------------------------------------

func TestLiveHead(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	for _, args := range [][]string{
		{"head", fx},
		{"head", "-n", "5", fx},
	} {
		t.Run(strings.Join(args[1:len(args)-1], "_"), func(t *testing.T) {
			assertEqualStdout(t, live, ours, args...)
		})
	}
}

// -------------------------------------------------------------------------
// NORM
// -------------------------------------------------------------------------

func TestLiveNorm(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	for _, args := range [][]string{
		{"norm", "-m-", fx},
		{"norm", "-m+", fx},
	} {
		t.Run(strings.Join(args[1:len(args)-1], "_"), func(t *testing.T) {
			assertEqualStdout(t, live, ours, args...)
		})
	}
}

// -------------------------------------------------------------------------
// FILTER
// -------------------------------------------------------------------------

func TestLiveFilter(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	cases := []struct {
		name string
		args []string
	}{
		{"include_qual_gt_30", []string{"filter", "-i", "QUAL>30", fx}},
		{"exclude_ac_lt_2", []string{"filter", "-e", "AC<2", fx}},
		{"soft_tag_lowq", []string{"filter", "-s", "LOWQ", "-e", "QUAL<10", fx}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			assertEqualStdout(t, live, ours, tc.args...)
		})
	}
}

// -------------------------------------------------------------------------
// SORT
// -------------------------------------------------------------------------

func TestLiveSort(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	assertEqualStdout(t, live, ours, "sort", fx)
}

// -------------------------------------------------------------------------
// CONCAT
// -------------------------------------------------------------------------

func TestLiveConcat(t *testing.T) {
	live, ours := requireLive(t)
	a := fixturePath(t, "concat_a.vcf")
	b := fixturePath(t, "concat_b.vcf")
	assertEqualStdout(t, live, ours, "concat", a, b)
}

// -------------------------------------------------------------------------
// MERGE
// -------------------------------------------------------------------------

func TestLiveMerge(t *testing.T) {
	live, ours := requireLive(t)
	tmp := t.TempDir()
	// Hand-rolled fixtures with disjoint sample IDs so upstream merges
	// cleanly without --force-samples (which our port doesn't yet
	// implement). Same contigs/INFOs so the header merge is trivial.
	const vcfA = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SA
chr1	100	.	A	T	30	PASS	DP=10	GT	0/1
chr1	200	.	C	G	40	PASS	DP=15	GT	1/1
`
	const vcfB = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SB
chr1	100	.	A	T	35	PASS	DP=12	GT	0/0
chr1	300	.	G	C	45	PASS	DP=20	GT	0/1
`
	aPlain := filepath.Join(tmp, "a.vcf")
	bPlain := filepath.Join(tmp, "b.vcf")
	if err := os.WriteFile(aPlain, []byte(vcfA), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(bPlain, []byte(vcfB), 0o644); err != nil {
		t.Fatal(err)
	}
	a := filepath.Join(tmp, "a.vcf.gz")
	b := filepath.Join(tmp, "b.vcf.gz")
	if err := os.WriteFile(a, runBin(t, live, "view", "-Oz", aPlain), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(b, runBin(t, live, "view", "-Oz", bPlain), 0o644); err != nil {
		t.Fatal(err)
	}
	runBin(t, live, "index", "-t", a)
	runBin(t, live, "index", "-t", b)
	assertEqualStdout(t, live, ours, "merge", a, b)
}

// -------------------------------------------------------------------------
// ISEC — runs both binaries in -p DIR mode and compares each output file
// after stripping provenance from the README.
// -------------------------------------------------------------------------

func TestLiveIsec(t *testing.T) {
	live, ours := requireLive(t)
	tmp := t.TempDir()
	a := filepath.Join(tmp, "a.vcf.gz")
	b := filepath.Join(tmp, "b.vcf.gz")
	if err := os.WriteFile(a, runBin(t, live, "view", "-Oz",
		fixturePath(t, "concat_a.vcf")), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(b, runBin(t, live, "view", "-Oz",
		fixturePath(t, "concat_b.vcf")), 0o644); err != nil {
		t.Fatal(err)
	}
	runBin(t, live, "index", "-t", a)
	runBin(t, live, "index", "-t", b)

	liveDir := filepath.Join(tmp, "live")
	oursDir := filepath.Join(tmp, "ours")
	runBin(t, live, "isec", "-p", liveDir, a, b)
	runBin(t, ours, "isec", "-p", oursDir, a, b)

	// Compare every file the live binary produced. Our port may
	// produce fewer (a known divergence — upstream emits four
	// per-input projections for two inputs; we emit two). When
	// a file is missing on our side we record an Errorf so the
	// gap is visible.
	entries, err := os.ReadDir(liveDir)
	if err != nil {
		t.Fatalf("read liveDir: %v", err)
	}
	for _, e := range entries {
		name := e.Name()
		if name == "README.txt" {
			// README wording diverges intentionally — we describe
			// our own semantics. Skip.
			continue
		}
		livBytes, err := os.ReadFile(filepath.Join(liveDir, name))
		if err != nil {
			t.Errorf("read live %s: %v", name, err)
			continue
		}
		ourBytes, err := os.ReadFile(filepath.Join(oursDir, name))
		if err != nil {
			t.Errorf("isec output %q present upstream but missing in our port: %v", name, err)
			continue
		}
		if !bytes.Equal(stripProvenance(livBytes), stripProvenance(ourBytes)) {
			t.Errorf("isec output %q differs", name)
		}
	}
}

// -------------------------------------------------------------------------
// REHEADER
// -------------------------------------------------------------------------

func TestLiveReheader(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	tmp := t.TempDir()
	// Tab-separated old<TAB>new mapping (the format upstream and our
	// port both accept unambiguously).
	mapPath := filepath.Join(tmp, "rename.tsv")
	if err := os.WriteFile(mapPath, []byte("S1\tNEW1\nS2\tNEW2\nS3\tNEW3\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	assertEqualStdout(t, live, ours, "reheader", "-s", mapPath, fx)
}

// -------------------------------------------------------------------------
// INDEX — compare the .tbi / .csi byte streams. These should be
// deterministic and identical when the input bytes match.
// -------------------------------------------------------------------------

func TestLiveIndexTbi(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	tmp := t.TempDir()

	gz := filepath.Join(tmp, "input.vcf.gz")
	if err := os.WriteFile(gz, runBin(t, live, "view", "-Oz", fx), 0o644); err != nil {
		t.Fatal(err)
	}
	liveGz := filepath.Join(tmp, "live.vcf.gz")
	oursGz := filepath.Join(tmp, "ours.vcf.gz")
	copyFileT(t, gz, liveGz)
	copyFileT(t, gz, oursGz)

	runBin(t, live, "index", "-t", liveGz)
	runBin(t, ours, "index", "-t", oursGz)

	livIdx, err := os.ReadFile(liveGz + ".tbi")
	if err != nil {
		t.Fatalf("read live tbi: %v", err)
	}
	ourIdx, err := os.ReadFile(oursGz + ".tbi")
	if err != nil {
		t.Fatalf("read our tbi: %v", err)
	}
	if !bytes.Equal(livIdx, ourIdx) {
		t.Errorf("index -t .tbi bytes diverge (live=%d bytes, ours=%d bytes)", len(livIdx), len(ourIdx))
	}
}

func TestLiveIndexCsi(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	tmp := t.TempDir()

	gz := filepath.Join(tmp, "input.vcf.gz")
	if err := os.WriteFile(gz, runBin(t, live, "view", "-Oz", fx), 0o644); err != nil {
		t.Fatal(err)
	}
	liveGz := filepath.Join(tmp, "live.vcf.gz")
	oursGz := filepath.Join(tmp, "ours.vcf.gz")
	copyFileT(t, gz, liveGz)
	copyFileT(t, gz, oursGz)

	runBin(t, live, "index", "-c", liveGz)
	runBin(t, ours, "index", "-c", oursGz)

	livIdx, err := os.ReadFile(liveGz + ".csi")
	if err != nil {
		t.Fatalf("read live csi: %v", err)
	}
	ourIdx, err := os.ReadFile(oursGz + ".csi")
	if err != nil {
		t.Fatalf("read our csi: %v", err)
	}
	if !bytes.Equal(livIdx, ourIdx) {
		t.Errorf("index -c .csi bytes diverge (live=%d bytes, ours=%d bytes)", len(livIdx), len(ourIdx))
	}
}

func copyFileT(t *testing.T, src, dst string) {
	t.Helper()
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatalf("copy read: %v", err)
	}
	if err := os.WriteFile(dst, b, 0o644); err != nil {
		t.Fatalf("copy write: %v", err)
	}
}

// -------------------------------------------------------------------------
// MPILEUP — smoke-test against one vendored fixture. The full matrix
// lives in mpileup_golden_test.go.
// -------------------------------------------------------------------------

func TestLiveMpileupSmoke(t *testing.T) {
	live, ours := requireLive(t)
	mpileupDir, err := filepath.Abs(filepath.Join("..", "..", "testdata", "mpileup"))
	if err != nil {
		t.Fatal(err)
	}
	bam := filepath.Join(mpileupDir, "mpileup.1.bam")
	ref := filepath.Join(mpileupDir, "mpileup.ref.fa")
	if _, err := os.Stat(bam); err != nil {
		t.Skip("no mpileup BAM fixture for live oracle")
	}
	if _, err := os.Stat(ref); err != nil {
		t.Skip("no mpileup reference fixture for live oracle")
	}
	// mpileup output is voluminous; we just verify both run and
	// agree on stdout for a minimal invocation. The mpileup_golden
	// suite covers the full matrix.
	assertEqualStdout(t, live, ours, "mpileup", "-f", ref, bam)
}

// -------------------------------------------------------------------------
// assertRejectionParity invokes a subcommand on both binaries with
// args expected to fail (e.g. missing required flags or unknown
// plugin) and asserts that both binaries exit non-zero and produce no
// stdout. This is the canonical fallback when we can't get byte-equal
// output for a subcommand (e.g. fundamental architectural divergence
// or upstream features we have not yet ported), but we still want a
// live-oracle gate that locks in observable behaviour. Modelled on
// tools/mosdepth/pkg/mosdepth/live_oracle_test.go.
// -------------------------------------------------------------------------

func assertRejectionParity(t *testing.T, args []string) {
	t.Helper()
	live, ours := requireLive(t)

	run := func(bin string) (rejected bool, stdoutBytes int) {
		cmd := exec.Command(bin, args...)
		var out bytes.Buffer
		cmd.Stdout = &out
		cmd.Stderr = io.Discard
		err := cmd.Run()
		return err != nil, out.Len()
	}

	liveRejected, liveStdout := run(live)
	if !liveRejected {
		t.Fatalf("upstream bcftools unexpectedly ACCEPTED %v; the oracle assumption is wrong", args)
	}
	if liveStdout > 0 {
		t.Fatalf("upstream bcftools wrote %d bytes to stdout for %v despite rejecting it",
			liveStdout, args)
	}
	oursRejected, oursStdout := run(ours)
	if !oursRejected {
		t.Errorf("our port accepted %v but upstream rejects it", args)
	}
	if oursStdout > 0 {
		t.Errorf("our port wrote %d bytes to stdout for %v but upstream produces none",
			oursStdout, args)
	}
}

// -------------------------------------------------------------------------
// CALL — the multiallelic caller (`call -m`) is a faithful port of
// mcall.c. We generate the mpileup VCF with the live binary, then assert
// our `call -m` output byte-matches upstream's over the whole contig
// (4000+ sites: EM allele-frequency estimation, per-site QUAL, the
// max-likelihood GT, and the INFO rewrite AN/AC/DP4/MQ). The same
// fixture is exercised with -v (variants only) and -A (keep alts).
// -------------------------------------------------------------------------

func TestLiveCall(t *testing.T) {
	live, ours := requireLive(t)
	mpDir, err := filepath.Abs(filepath.Join("..", "..", "testdata", "mpileup"))
	if err != nil {
		t.Fatal(err)
	}
	ref := filepath.Join(mpDir, "mpileup.ref.fa")
	bam := filepath.Join(mpDir, "mpileup.1.bam")
	if _, err := os.Stat(bam); err != nil {
		t.Skip("no mpileup fixtures for live oracle")
	}

	// Produce the mpileup VCF with the genuine binary so the caller input
	// is upstream-identical; both `call` binaries then consume it.
	mp := runBin(t, live, "mpileup", "-f", ref, bam)
	mpFile := filepath.Join(t.TempDir(), "mp.vcf")
	if err := os.WriteFile(mpFile, mp, 0o644); err != nil {
		t.Fatal(err)
	}

	for _, tc := range []struct {
		name string
		args []string
	}{
		{"m", []string{"call", "-m", mpFile}},
		{"m_v", []string{"call", "-m", "-v", mpFile}},
		{"m_A", []string{"call", "-m", "-A", mpFile}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			assertEqualStdout(t, live, ours, tc.args...)
		})
	}
}

// -------------------------------------------------------------------------
// CSQ — exercise the splice / out-of-bounds-codon vendored fixtures.
// Both fixtures produce byte-equal output between the live binary and
// our port (modulo provenance lines).
// -------------------------------------------------------------------------

func TestLiveCsq(t *testing.T) {
	live, ours := requireLive(t)
	csqDir, err := filepath.Abs(filepath.Join("..", "..", "testdata", "csq"))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(csqDir); err != nil {
		t.Skip("no csq fixtures for live oracle")
	}
	cases := []struct {
		name string
		fa   string
		gff  string
		vcf  string
	}{
		{"oob-codon", "csq.oob-codon.fa", "csq.oob-codon.gff", "csq.oob-codon.vcf"},
		{"splice-2543", "csq.splice.issue-2543.fa", "csq.splice.issue-2543.gff", "csq.splice.issue-2543.vcf"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			assertEqualStdout(t, live, ours,
				"csq",
				"-f", filepath.Join(csqDir, tc.fa),
				"-g", filepath.Join(csqDir, tc.gff),
				filepath.Join(csqDir, tc.vcf))
		})
	}

	// -i / -e and -s on the single-sample splice fixture. Upstream's
	// -i/-e gate consequence calling (not record emission), and -s
	// restricts which samples drive consequence calling and the
	// FORMAT/BCSQ bitmask; assert both byte-match upstream.
	splFA := filepath.Join(csqDir, "csq.splice.issue-2543.fa")
	splGFF := filepath.Join(csqDir, "csq.splice.issue-2543.gff")
	splVCF := filepath.Join(csqDir, "csq.splice.issue-2543.vcf")
	flagCases := []struct {
		name string
		args []string
	}{
		{"include", []string{"-i", "QUAL=60"}},
		{"include-filter", []string{"-i", `FILTER="PASS"`}},
		{"exclude", []string{"-e", "QUAL=60"}},
		{"exclude-keepall", []string{"-e", "QUAL<10"}},
		{"samples", []string{"-s", "snippy"}},
		{"samples-none", []string{"-s", "-"}},
	}
	for _, tc := range flagCases {
		t.Run(tc.name, func(t *testing.T) {
			args := append([]string{"csq", "-f", splFA, "-g", splGFF}, tc.args...)
			args = append(args, splVCF)
			assertEqualStdout(t, live, ours, args...)
		})
	}
}

// -------------------------------------------------------------------------
// ROH — exercise the vendored roh.1 fixture with the --AF-dflt code
// path (no need to bgzip-index the .tab.gz). Both binaries emit the
// per-region summary; the data row is identical.
// -------------------------------------------------------------------------

func TestLiveRoh(t *testing.T) {
	live, ours := requireLive(t)
	rohDir, err := filepath.Abs(filepath.Join("..", "..", "testdata", "roh"))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(rohDir); err != nil {
		t.Skip("no roh fixtures for live oracle")
	}
	vcfgz := filepath.Join(rohDir, "roh.1.vcf.gz")
	assertEqualStdout(t, live, ours,
		"roh", "-Or", "-G30", "--AF-dflt", "0.4", vcfgz)
}

// -------------------------------------------------------------------------
// MENDELIAN — `+mendelian2` is upstream's dlopen plugin; our port runs
// it as a native built-in. The `-m c` summary (counter labels, the
// descriptive comment lines, and the per-trio columns) now byte-matches
// upstream plugins/mendelian2.c. The upstream binary needs
// BCFTOOLS_PLUGINS pointed at the vendored .so directory; our built-in
// needs no plugin path.
// -------------------------------------------------------------------------

func TestLiveMendelian(t *testing.T) {
	live, ours := requireLive(t)
	pluginDir := livePluginDir(t)
	if pluginDir == "" {
		t.Skip("vendored bcftools plugin .so directory not found")
	}
	mendDir, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", "mendelian"))
	if err != nil {
		t.Fatal(err)
	}
	vcfPath := filepath.Join(mendDir, "trio.vcf")
	pedPath := filepath.Join(mendDir, "trio.ped")
	if _, err := os.Stat(vcfPath); err != nil {
		t.Skip("no mendelian trio fixture")
	}
	env := []string{"BCFTOOLS_PLUGINS=" + pluginDir}
	livOut := runBinEnv(t, live, env, "+mendelian2", vcfPath, "-P", pedPath)
	ourOut := runBin(t, ours, "+mendelian2", vcfPath, "-P", pedPath)
	if !bytes.Equal(stripProvenance(livOut), stripProvenance(ourOut)) {
		t.Errorf("+mendelian2 summary diverges\n--- live ---\n%s\n--- ours ---\n%s",
			snippet(livOut, 800), snippet(ourOut, 800))
	}
}

// -------------------------------------------------------------------------
// GTCHECK — full multi-section output parity in cross-check mode. The
// port reproduces upstream's INFO counter block, the DCv2 comment
// block, the DCv2 header row, and the error-probability discordance /
// HWE data rows byte-for-byte (modulo the provenance banner stripped by
// stripProvenance).
// -------------------------------------------------------------------------

func TestLiveGtcheck(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	// Cross-check mode (no -g) accepts a plain VCF on both sides.
	assertEqualStdout(t, live, ours, "gtcheck", fx)
}

// -------------------------------------------------------------------------
// POLYSOMY — upstream's binary does not even ship the `polysomy`
// subcommand (it errors with "unrecognized command 'polysomy'"). Our
// port has a stub that requires an input file. Use rejection-parity
// to confirm both refuse to run.
// -------------------------------------------------------------------------

func TestLivePolysomy(t *testing.T) {
	assertRejectionParity(t, []string{"polysomy"})
}

// -------------------------------------------------------------------------
// CNV — both binaries refuse to run without -o / --output-dir (live
// emits "Expected -o option"; ours emits "missing -o/--output-dir").
// -------------------------------------------------------------------------

func TestLiveCnv(t *testing.T) {
	assertRejectionParity(t, []string{"cnv"})
}

func TestLiveConvert(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	// `convert` without subformat flags is effectively `view`; that's
	// the smoke we can lean on universally.
	assertEqualStdout(t, live, ours, "convert", fx)
}

// -------------------------------------------------------------------------
// CONSENSUS — apply two SNPs from a vendored VCF to a 160-base
// reference. Output is the modified fasta; both binaries agree.
// -------------------------------------------------------------------------

func TestLiveConsensus(t *testing.T) {
	live, ours := requireLive(t)
	consDir, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", "consensus"))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(consDir); err != nil {
		t.Skip("no consensus parity fixtures")
	}
	ref := filepath.Join(consDir, "ref.fa")
	vcfPlain := filepath.Join(consDir, "variants.vcf")

	tmp := t.TempDir()
	gz := filepath.Join(tmp, "variants.vcf.gz")
	if err := os.WriteFile(gz, runBin(t, live, "view", "-Oz", vcfPlain), 0o644); err != nil {
		t.Fatal(err)
	}
	runBin(t, live, "index", "-t", gz)
	assertEqualStdout(t, live, ours, "consensus", "-f", ref, gz)
}

// -------------------------------------------------------------------------
// PLUGIN — upstream loads plugins via dlopen of `<name>.so`; our port
// re-implements the common ones (`fill-tags`, `mendelian2`) as native
// Go built-ins so their output byte-matches upstream. `+fill-tags`
// with the default "all" tag set is compared end-to-end (AC/AN/AF/MAF/
// NS/AC_Hom/Het/Hemi/F_MISSING/HWE/ExcHet, including the HWE exact
// test). A still-missing plugin name keeps rejection-parity.
// -------------------------------------------------------------------------

func TestLivePlugin(t *testing.T) {
	live, ours := requireLive(t)
	pluginDir := livePluginDir(t)
	if pluginDir == "" {
		t.Skip("vendored bcftools plugin .so directory not found")
	}
	fx := fixturePath(t, "basic.vcf")
	env := []string{"BCFTOOLS_PLUGINS=" + pluginDir}
	livOut := runBinEnv(t, live, env, "+fill-tags", fx)
	ourOut := runBin(t, ours, "+fill-tags", fx)
	if !bytes.Equal(stripProvenance(livOut), stripProvenance(ourOut)) {
		t.Errorf("+fill-tags output diverges\n--- live ---\n%s\n--- ours ---\n%s",
			snippet(livOut, 1200), snippet(ourOut, 1200))
	}
}

// TestLivePluginRejection keeps the rejection-parity gate for an
// unknown plugin name (neither binary can resolve it).
func TestLivePluginRejection(t *testing.T) {
	fx := fixturePath(t, "basic.vcf")
	assertRejectionParity(t, []string{"+nosuchplugin", fx})
}
