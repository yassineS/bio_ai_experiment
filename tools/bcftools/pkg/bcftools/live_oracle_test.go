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

// Provenance lines that legitimately differ between upstream and our port
// (they encode version strings and absolute paths). Stripping them keeps
// the oracle focused on actual semantic divergence.
var (
	provVersionRE   = regexp.MustCompile(`(?m)^##bcftools_[^=]+=.*\n`)
	provVersionRE2  = regexp.MustCompile(`(?m)^##bcftoolsVersion=.*\n`)
	provVersionRE3  = regexp.MustCompile(`(?m)^##bcftoolsCommand=.*\n`)
	provReferenceRE = regexp.MustCompile(`(?m)^##reference=.*\n`)
	provStatsRE1    = regexp.MustCompile(`(?m)^# This file was produced by .*\n`)
	provStatsRE2    = regexp.MustCompile(`(?m)^# The command line was:.*\n`)
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
// Skip-only stubs for subcommands without committed live-oracle
// fixtures. Each carries a TODO pointing at what's needed to enable it.
// -------------------------------------------------------------------------

func TestLiveCall(t *testing.T) {
	requireLive(t)
	t.Skip("TODO: needs an mpileup VCF input fixture committed to testdata/parity/")
}

func TestLiveCsq(t *testing.T) {
	requireLive(t)
	csqDir, err := filepath.Abs(filepath.Join("..", "..", "testdata", "csq"))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(csqDir); err != nil {
		t.Skip("no csq fixtures for live oracle")
	}
	t.Skip("TODO: csq requires (gff, fa, vcf) triple; wire up to a single representative case under testdata/csq/")
}

func TestLiveRoh(t *testing.T) {
	requireLive(t)
	rohDir, err := filepath.Abs(filepath.Join("..", "..", "testdata", "roh"))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(rohDir); err != nil {
		t.Skip("no roh fixtures for live oracle")
	}
	t.Skip("TODO: roh needs --AF-tag/--AF-file plumbing; smoke-test under testdata/roh/ to be added")
}

func TestLiveMendelian(t *testing.T) {
	requireLive(t)
	t.Skip("TODO: needs a trio VCF + PED fixture under testdata/parity/")
}

func TestLiveGtcheck(t *testing.T) {
	requireLive(t)
	t.Skip("TODO: needs two indexed VCFs with overlapping samples for cross-check")
}

func TestLivePolysomy(t *testing.T) {
	requireLive(t)
	t.Skip("TODO: needs a BAF-distribution input fixture")
}

func TestLiveCnv(t *testing.T) {
	requireLive(t)
	t.Skip("TODO: needs a CNV-ready input (tumour/normal pair)")
}

func TestLiveConvert(t *testing.T) {
	live, ours := requireLive(t)
	fx := fixturePath(t, "basic.vcf")
	// `convert` without subformat flags is effectively `view`; that's
	// the smoke we can lean on universally.
	assertEqualStdout(t, live, ours, "convert", fx)
}

func TestLiveConsensus(t *testing.T) {
	requireLive(t)
	t.Skip("TODO: needs a small fasta + matching vcf.gz under testdata/parity/consensus/")
}

func TestLivePlugin(t *testing.T) {
	requireLive(t)
	t.Skip("TODO: our plugin surface diverges (+fill-tags etc.); to be wired once plugin discovery is unified")
}
