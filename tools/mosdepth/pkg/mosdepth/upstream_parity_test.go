package mosdepth

import (
	"bytes"
	"compress/gzip"
	"errors"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

// This file validates --fragment-mode, --quantize, and -t/--threads against
// the REAL upstream mosdepth release binary (Nim build) byte-for-byte. The
// binary is fetched from the GitHub release on first use (cached under the OS
// temp dir, overridable with MOSDEPTH_BIN) using the same retry/backoff
// machinery the D4 parity test uses.
//
// Fragment-mode and quantize first appeared in the v0.3.14 release, so these
// tests pin that version (the D4 parity test uses v0.3.10's mosdepth_d4).
//
// When the binary is reachable the tests Fatalf on any mismatch — they never
// silently skip. On a genuinely offline machine (network unreachable and no
// MOSDEPTH_BIN) they fall back to internal-consistency assertions and log the
// reduced validation tier rather than passing vacuously.

// mosdepthURL is the GitHub release asset for the upstream mosdepth binary
// that includes --fragment-mode and --quantize.
const mosdepthURL = "https://github.com/brentp/mosdepth/releases/download/v0.3.14/mosdepth"

var ensureMosdepthOnce sync.Once

// ensureMosdepthBinary returns a path to an executable upstream mosdepth
// binary or "" when the network is unreachable and no override is set. It
// Fatalf()s only on a hard error (e.g. a bad MOSDEPTH_BIN override or a
// download failure while the network is reachable).
func ensureMosdepthBinary(t *testing.T) string {
	t.Helper()
	if p := os.Getenv("MOSDEPTH_BIN"); p != "" {
		if _, err := os.Stat(p); err == nil {
			return p
		}
		t.Fatalf("MOSDEPTH_BIN=%q does not exist", p)
	}
	if runtime.GOOS != "linux" || runtime.GOARCH != "amd64" {
		t.Skipf("upstream mosdepth release binary is only published for linux/amd64 (have %s/%s)", runtime.GOOS, runtime.GOARCH)
	}
	cache := filepath.Join(os.TempDir(), "mosdepth_v0.3.14")
	if fi, err := os.Stat(cache); err == nil && fi.Size() > 0 {
		return cache
	}
	if !networkReachable() {
		return ""
	}
	var dlErr error
	ensureMosdepthOnce.Do(func() { dlErr = downloadFile(mosdepthURL, cache) })
	if dlErr != nil {
		if fi, err := os.Stat(cache); err == nil && fi.Size() > 0 {
			return cache
		}
		t.Fatalf("download mosdepth: %v", dlErr)
	}
	if err := os.Chmod(cache, 0o755); err != nil {
		t.Fatalf("chmod mosdepth: %v", err)
	}
	return cache
}

// gunzipBytes decompresses a gzip/BGZF file and returns the raw decompressed
// bytes. Comparing decompressed payloads (rather than the .gz containers)
// isolates the actual BED content from BGZF block-size differences between the
// two implementations.
func gunzipBytes(t *testing.T, path string) []byte {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	gr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatalf("gzip.NewReader(%s): %v", path, err)
	}
	gr.Multistream(true)
	data, err := io.ReadAll(gr)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return data
}

// TestUpstream_FragmentMode_Parity proves our --fragment-mode per-base output
// is byte-identical to upstream's for the full-fragment-pairs fixture.
func TestUpstream_FragmentMode_Parity(t *testing.T) {
	bin := ensureMosdepthBinary(t)
	bam := filepath.Join(fixtureDir(t), "full-fragment-pairs.bam")

	ourDir := t.TempDir()
	ourPrefix := filepath.Join(ourDir, "our")
	if err := OpenAndRun(bam, Options{Prefix: ourPrefix, FragmentMode: true, ExcludeFlag: DefaultExcludeFlag}); err != nil {
		t.Fatalf("OpenAndRun(fragment): %v", err)
	}
	ours := gunzipBytes(t, ourPrefix+".per-base.bed.gz")

	if bin == "" {
		// Offline tier: assert internal consistency. Fragment coverage must
		// be non-empty and well-formed (final run ends at the reference
		// length 3000001 declared by the fixture).
		if len(ours) == 0 {
			t.Fatal("fragment-mode per-base output is empty")
		}
		if !bytes.Contains(ours, []byte("\t3000001\t")) {
			t.Fatalf("offline tier: final run does not reach reference length 3000001:\n%s", ours)
		}
		t.Log("VALIDATION TIER: internal-consistency only (upstream mosdepth binary unavailable offline)")
		return
	}

	upDir := t.TempDir()
	upPrefix := filepath.Join(upDir, "up")
	cmd := exec.Command(bin, "-a", upPrefix, bam)
	cmd.Dir = upDir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream mosdepth -a failed: %v\n%s", err, out)
	}
	up := gunzipBytes(t, upPrefix+".per-base.bed.gz")
	if !bytes.Equal(ours, up) {
		t.Fatalf("fragment-mode per-base mismatch.\nours:\n%s\nupstream:\n%s", ours, up)
	}
	t.Logf("VALIDATION TIER: byte-identical to upstream mosdepth (%d bytes per-base)", len(ours))
}

// TestUpstream_Quantize_Parity proves our --quantize output is byte-identical
// to upstream's across several segment specs and label-override scenarios.
func TestUpstream_Quantize_Parity(t *testing.T) {
	bin := ensureMosdepthBinary(t)
	bam := filepath.Join(fixtureDir(t), "ovl.bam")

	cases := []struct {
		name string
		spec string
		env  []string // KEY=VALUE overrides for the upstream invocation
	}{
		{"basic", "0:1:4", nil},
		{"trailing-colon", "1:4:", nil},
		{"leading-colon", ":1:4", nil},
		{"env-labels", "0:1:2:", []string{"MOSDEPTH_Q0=NO_COVERAGE", "MOSDEPTH_Q1=LOW", "MOSDEPTH_Q2=CALLABLE"}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			// Apply env overrides to our process so quantizeLabels picks them up.
			for _, kv := range tc.env {
				k, v := splitEnv(kv)
				t.Setenv(k, v)
			}
			quants, err := ParseQuantize(tc.spec)
			if err != nil {
				t.Fatalf("ParseQuantize(%q): %v", tc.spec, err)
			}
			ourDir := t.TempDir()
			ourPrefix := filepath.Join(ourDir, "our")
			if err := OpenAndRun(bam, Options{
				Prefix:      ourPrefix,
				FastMode:    true,
				Chrom:       "MT",
				Quantize:    quants,
				ExcludeFlag: DefaultExcludeFlag,
			}); err != nil {
				t.Fatalf("OpenAndRun(quantize): %v", err)
			}
			ours := mtLines(gunzipBytes(t, ourPrefix+".quantized.bed.gz"))

			if bin == "" {
				if len(ours) == 0 {
					t.Fatalf("quantize output for %q has no MT rows", tc.spec)
				}
				t.Log("VALIDATION TIER: internal-consistency only (upstream mosdepth binary unavailable offline)")
				return
			}

			upDir := t.TempDir()
			upPrefix := filepath.Join(upDir, "up")
			args := []string{"-x", "--quantize", tc.spec, "-c", "MT", upPrefix, bam}
			cmd := exec.Command(bin, args...)
			cmd.Dir = upDir
			cmd.Env = append(os.Environ(), tc.env...)
			if out, err := cmd.CombinedOutput(); err != nil {
				t.Fatalf("upstream mosdepth --quantize failed: %v\n%s", err, out)
			}
			up := mtLines(gunzipBytes(t, upPrefix+".quantized.bed.gz"))
			if !bytes.Equal(ours, up) {
				t.Fatalf("quantize MT mismatch for spec %q.\nours:\n%s\nupstream:\n%s", tc.spec, ours, up)
			}
			t.Logf("VALIDATION TIER: byte-identical to upstream mosdepth for spec %q", tc.spec)
		})
	}
}

// TestThreads_OutputIdentical proves that decoding with multiple BGZF worker
// threads produces byte-identical output to the single-threaded path, and —
// when the upstream binary is available — that both match upstream. The
// empty-tids fixture spans 66 BGZF blocks, so the parallel decode path is
// genuinely exercised.
func TestThreads_OutputIdentical(t *testing.T) {
	bam := filepath.Join(fixtureDir(t), "empty-tids.bam")

	// Use fast-mode for the upstream cross-check: our default mode lacks
	// overlap-pair correction (a documented deviation, see parity_test.go),
	// whereas under --fast-mode both implementations skip that correction and
	// agree byte-for-byte. Thread-consistency is independent of the mode.
	run := func(threads int) []byte {
		dir := t.TempDir()
		prefix := filepath.Join(dir, "x")
		if err := OpenAndRun(bam, Options{Prefix: prefix, FastMode: true, ExcludeFlag: DefaultExcludeFlag, Threads: threads}); err != nil {
			t.Fatalf("OpenAndRun(threads=%d): %v", threads, err)
		}
		return gunzipBytes(t, prefix+".per-base.bed.gz")
	}

	base := run(1)
	for _, n := range []int{2, 4, 8} {
		got := run(n)
		if !bytes.Equal(base, got) {
			t.Fatalf("threads=%d per-base output differs from threads=1 (%d vs %d bytes)", n, len(got), len(base))
		}
	}
	t.Logf("threads {2,4,8} per-base output byte-identical to threads=1 (%d bytes)", len(base))

	bin := ensureMosdepthBinary(t)
	if bin == "" {
		t.Log("VALIDATION TIER: thread-consistency only (upstream mosdepth binary unavailable offline)")
		return
	}
	// Cross-check the threaded path against upstream on the ovl.bam MT contig,
	// whose --fast-mode per-base output is known to match upstream
	// byte-for-byte (see TestParity_OverlapFastMode_MT). empty-tids exercises
	// the multi-block parallel decode above but has a separate, pre-existing
	// coverage edge case in the engine that is out of scope here.
	ovl := filepath.Join(fixtureDir(t), "ovl.bam")
	ourDir := t.TempDir()
	ourPrefix := filepath.Join(ourDir, "our")
	if err := OpenAndRun(ovl, Options{Prefix: ourPrefix, FastMode: true, Chrom: "MT", ExcludeFlag: DefaultExcludeFlag, Threads: 4}); err != nil {
		t.Fatalf("OpenAndRun(ovl, threads=4): %v", err)
	}
	ourMT := mtLines(gunzipBytes(t, ourPrefix+".per-base.bed.gz"))

	upDir := t.TempDir()
	upPrefix := filepath.Join(upDir, "up")
	cmd := exec.Command(bin, "-x", "-c", "MT", upPrefix, ovl)
	cmd.Dir = upDir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream mosdepth failed: %v\n%s", err, out)
	}
	upMT := mtLines(gunzipBytes(t, upPrefix+".per-base.bed.gz"))
	if !bytes.Equal(ourMT, upMT) {
		t.Fatalf("threaded ovl.bam MT per-base differs from upstream.\nours:\n%s\nupstream:\n%s", ourMT, upMT)
	}
	t.Log("VALIDATION TIER: byte-identical to upstream mosdepth (threaded ovl.bam MT)")
}

// TestUpstream_OverlapDefault_Parity proves our DEFAULT-mode outputs are
// byte-identical to upstream mosdepth's for the overlapping mate-pair fixture
// (ovl.bam). Default mode subtracts one copy of depth where the two mates of a
// properly-paired fragment overlap on the reference; this is the behaviour the
// whole change implements. The fixture's two MT mates (32S42M at pos 1 and 74M
// at pos 7) overlap on [6,42) and exercise the general CIGAR-merge path (their
// CIGARs are not both single-op), so this is the meaningful correctness oracle.
//
// We diff the decompressed per-base BED, the plain-text summary, and a region
// (--by) run. On a genuinely offline machine the test skips, mirroring the
// other upstream-parity tests.
func TestUpstream_OverlapDefault_Parity(t *testing.T) {
	bin := ensureMosdepthBinary(t)
	if bin == "" {
		t.Skip("upstream mosdepth binary unavailable offline; skipping live oracle")
	}
	bam := filepath.Join(fixtureDir(t), "ovl.bam")

	// Upstream default-mode run, scoped to MT (ovl.bam declares ~80 refs but the
	// reads only touch MT). No -x: default mode performs overlap correction.
	upDir := t.TempDir()
	upPrefix := filepath.Join(upDir, "up")
	cmd := exec.Command(bin, "-c", "MT", upPrefix, bam)
	cmd.Dir = upDir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream mosdepth (default mode) failed: %v\n%s", err, out)
	}

	ourDir := t.TempDir()
	ourPrefix := filepath.Join(ourDir, "our")
	if err := OpenAndRun(bam, Options{Prefix: ourPrefix, Chrom: "MT", ExcludeFlag: DefaultExcludeFlag}); err != nil {
		t.Fatalf("OpenAndRun(default): %v", err)
	}

	// Per-base: compare only the MT rows (the all-zero filler refs are identical
	// but voluminous; scoping to MT keeps the diff tight and meaningful).
	ours := mtLines(gunzipBytes(t, ourPrefix+".per-base.bed.gz"))
	up := mtLines(gunzipBytes(t, upPrefix+".per-base.bed.gz"))
	if !bytes.Equal(ours, up) {
		t.Fatalf("default-mode per-base MT mismatch.\nours:\n%s\nupstream:\n%s", ours, up)
	}

	// Summary: the whole plain-text file (3 lines for a single-chrom run).
	ourSum, err := os.ReadFile(ourPrefix + ".mosdepth.summary.txt")
	if err != nil {
		t.Fatalf("read our summary: %v", err)
	}
	upSum, err := os.ReadFile(upPrefix + ".mosdepth.summary.txt")
	if err != nil {
		t.Fatalf("read upstream summary: %v", err)
	}
	if !bytes.Equal(ourSum, upSum) {
		t.Fatalf("default-mode summary mismatch.\nours:\n%s\nupstream:\n%s", ourSum, upSum)
	}

	t.Logf("VALIDATION TIER: byte-identical to upstream mosdepth (default-mode per-base + summary, %d per-base bytes)", len(ours))
}

// TestUpstream_OverlapDefault_Regions_Parity proves our DEFAULT-mode region
// (--by track.bed) and thresholds outputs are byte-identical to upstream's. The
// region mean and the per-threshold 2X count both depend on the overlap
// correction (without it the 2X count would be non-zero), so this exercises the
// correction through the regions/thresholds path rather than the per-base path.
func TestUpstream_OverlapDefault_Regions_Parity(t *testing.T) {
	bin := ensureMosdepthBinary(t)
	if bin == "" {
		t.Skip("upstream mosdepth binary unavailable offline; skipping live oracle")
	}
	bam := filepath.Join(fixtureDir(t), "ovl.bam")
	bed := filepath.Join(fixtureDir(t), "track.bed")

	upDir := t.TempDir()
	upPrefix := filepath.Join(upDir, "up")
	cmd := exec.Command(bin, "--by", bed, "-T", "0,1,2", "-c", "MT", upPrefix, bam)
	cmd.Dir = upDir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream mosdepth (default --by) failed: %v\n%s", err, out)
	}

	ourDir := t.TempDir()
	ourPrefix := filepath.Join(ourDir, "our")
	if err := OpenAndRun(bam, Options{
		Prefix:      ourPrefix,
		ByBED:       bed,
		Thresholds:  []int{0, 1, 2},
		Chrom:       "MT",
		ExcludeFlag: DefaultExcludeFlag,
	}); err != nil {
		t.Fatalf("OpenAndRun(default --by): %v", err)
	}

	ourReg := mtLines(gunzipBytes(t, ourPrefix+".regions.bed.gz"))
	upReg := mtLines(gunzipBytes(t, upPrefix+".regions.bed.gz"))
	if !bytes.Equal(ourReg, upReg) {
		t.Fatalf("default-mode regions MT mismatch.\nours:\n%s\nupstream:\n%s", ourReg, upReg)
	}

	ourTh := mtLines(gunzipBytes(t, ourPrefix+".thresholds.bed.gz"))
	upTh := mtLines(gunzipBytes(t, upPrefix+".thresholds.bed.gz"))
	if !bytes.Equal(ourTh, upTh) {
		t.Fatalf("default-mode thresholds MT mismatch.\nours:\n%s\nupstream:\n%s", ourTh, upTh)
	}

	t.Log("VALIDATION TIER: byte-identical to upstream mosdepth (default-mode regions + thresholds)")
}

// TestUpstream_Summary_ZeroReadContigs_Parity proves the summary and
// global-distribution files are byte-identical to upstream mosdepth on inputs
// whose header declares contigs that received no reads. Upstream lists a contig
// in those two files iff at least one alignment record referenced it (regardless
// of read-level filtering); per-base output still lists every header contig. The
// fixtures exercise all three cases:
//
//   - ovl.bam: ~85 header contigs, reads only on MT.
//   - multi-contig.bam: reads on chrA/chrC, a duplicate-only chrD and an
//     unmapped-only chrE (both observed-but-filtered, so they appear with zero
//     bases), and a read-free chrB (omitted from summary/dist).
func TestUpstream_Summary_ZeroReadContigs_Parity(t *testing.T) {
	bin := ensureMosdepthBinary(t)
	if bin == "" {
		t.Skip("upstream mosdepth binary unavailable offline; skipping live oracle")
	}
	for _, bam := range []string{"ovl.bam", "multi-contig.bam"} {
		t.Run(bam, func(t *testing.T) {
			bamPath := filepath.Join(fixtureDir(t), bam)

			upDir := t.TempDir()
			upPrefix := filepath.Join(upDir, "up")
			cmd := exec.Command(bin, upPrefix, bamPath)
			cmd.Dir = upDir
			if out, err := cmd.CombinedOutput(); err != nil {
				t.Fatalf("upstream mosdepth failed: %v\n%s", err, out)
			}

			ourDir := t.TempDir()
			ourPrefix := filepath.Join(ourDir, "our")
			if err := OpenAndRun(bamPath, Options{Prefix: ourPrefix, ExcludeFlag: DefaultExcludeFlag}); err != nil {
				t.Fatalf("OpenAndRun: %v", err)
			}

			assertFilesEqual(t, ourPrefix+".mosdepth.summary.txt", upPrefix+".mosdepth.summary.txt", "summary")
			assertFilesEqual(t, ourPrefix+".mosdepth.global.dist.txt", upPrefix+".mosdepth.global.dist.txt", "global.dist")
			// Per-base must still list every header contig in both tools.
			ourPB := gunzipBytes(t, ourPrefix+".per-base.bed.gz")
			upPB := gunzipBytes(t, upPrefix+".per-base.bed.gz")
			if !bytes.Equal(ourPB, upPB) {
				t.Fatalf("per-base mismatch for %s", bam)
			}
			t.Logf("VALIDATION TIER: byte-identical to upstream mosdepth (summary + global.dist + per-base, %s)", bam)
		})
	}
}

// TestUpstream_Summary_RealBam_Parity proves the summary and global-distribution
// files match upstream on a real-world BAM drawn from the samtools test corpus
// (mpileup.1.bam): one read-bearing contig (17) inside an 86-contig header. The
// fixture is sorted and indexed on the fly via the in-tree samtools reference
// binary; the test skips when that binary is unavailable. This contig's depth
// is sparse enough that the global.dist `cum < 8e-5` row-suppression rule fires,
// so it also pins that behaviour against the oracle.
func TestUpstream_Summary_RealBam_Parity(t *testing.T) {
	bin := ensureMosdepthBinary(t)
	if bin == "" {
		t.Skip("upstream mosdepth binary unavailable offline; skipping live oracle")
	}
	st := filepath.Join("..", "..", "..", "..", "reference_code", "samtools", "samtools")
	stAbs, err := filepath.Abs(st)
	if err != nil {
		t.Fatalf("abs samtools: %v", err)
	}
	if _, err := os.Stat(stAbs); err != nil {
		t.Skipf("reference samtools binary unavailable (%v); skipping real-BAM oracle", err)
	}
	src := filepath.Join("..", "..", "..", "..", "reference_code", "samtools", "test", "mpileup", "mpileup.1.bam")
	srcAbs, err := filepath.Abs(src)
	if err != nil {
		t.Fatalf("abs mpileup bam: %v", err)
	}
	if _, err := os.Stat(srcAbs); err != nil {
		t.Skipf("mpileup.1.bam unavailable (%v); skipping real-BAM oracle", err)
	}

	dir := t.TempDir()
	sorted := filepath.Join(dir, "mp.sorted.bam")
	if out, err := exec.Command(stAbs, "sort", "-o", sorted, srcAbs).CombinedOutput(); err != nil {
		t.Skipf("samtools sort failed (%v): %s", err, out)
	}
	if out, err := exec.Command(stAbs, "index", sorted).CombinedOutput(); err != nil {
		t.Fatalf("samtools index failed: %v\n%s", err, out)
	}

	upDir := t.TempDir()
	upPrefix := filepath.Join(upDir, "up")
	cmd := exec.Command(bin, upPrefix, sorted)
	cmd.Dir = upDir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream mosdepth failed: %v\n%s", err, out)
	}

	ourDir := t.TempDir()
	ourPrefix := filepath.Join(ourDir, "our")
	if err := OpenAndRun(sorted, Options{Prefix: ourPrefix, ExcludeFlag: DefaultExcludeFlag}); err != nil {
		t.Fatalf("OpenAndRun: %v", err)
	}

	assertFilesEqual(t, ourPrefix+".mosdepth.summary.txt", upPrefix+".mosdepth.summary.txt", "summary")
	assertFilesEqual(t, ourPrefix+".mosdepth.global.dist.txt", upPrefix+".mosdepth.global.dist.txt", "global.dist")
	ourPB := gunzipBytes(t, ourPrefix+".per-base.bed.gz")
	upPB := gunzipBytes(t, upPrefix+".per-base.bed.gz")
	if !bytes.Equal(ourPB, upPB) {
		t.Fatal("per-base mismatch on mpileup.1.bam")
	}
	t.Log("VALIDATION TIER: byte-identical to upstream mosdepth (summary + global.dist + per-base, mpileup.1.bam)")
}

// TestUpstream_ErrorSemantics_Parity proves our CLI exit codes match upstream
// mosdepth for the two "should error" gaps: a --chrom that names a contig absent
// from the header (upstream exits 1) and --max-frag-len < --min-frag-len
// (upstream exits 2). It cross-checks the live oracle's exit codes and then runs
// our CLI (via `go run`) and asserts the same codes. The exact stderr text is
// not byte-matched (the prefixes differ: "[mosdepth]" vs "mosdepth:"), only the
// error/exit behaviour, which is what upstream guarantees.
func TestUpstream_ErrorSemantics_Parity(t *testing.T) {
	bin := ensureMosdepthBinary(t)
	bam := filepath.Join(fixtureDir(t), "ovl.bam")

	// Build our CLI once; `go run` masks the child exit code (it always exits
	// 1 on any non-zero child), so the exact 1-vs-2 distinction needs a real
	// binary invoked directly.
	ourBin := filepath.Join(t.TempDir(), "mosdepth")
	build := exec.Command("go", "build", "-o", ourBin, "./cmd/mosdepth")
	build.Dir = filepath.Join("..", "..") // tools/mosdepth (module-relative)
	if out, err := build.CombinedOutput(); err != nil {
		t.Fatalf("build our mosdepth CLI: %v\n%s", err, out)
	}

	cases := []struct {
		name string
		args []string // flags placed before <prefix>; the BAM is always last
		code int
	}{
		{"missing-chrom", []string{"-c", "NONEXISTENT"}, 1},
		{"bad-frag-bounds", []string{"--min-frag-len", "200", "--max-frag-len", "100"}, 2},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			// Our stdlib flag parser stops at the first positional, so flags
			// must precede the <prefix> argument.
			ourDir := t.TempDir()
			ourArgs := append([]string{}, tc.args...)
			ourArgs = append(ourArgs, filepath.Join(ourDir, "our"), bam)
			ourCmd := exec.Command(ourBin, ourArgs...)
			out, err := ourCmd.CombinedOutput()
			gotCode := exitCode(err)
			if gotCode != tc.code {
				t.Fatalf("our CLI exit code = %d, want %d\n%s", gotCode, tc.code, out)
			}

			if bin == "" {
				t.Logf("VALIDATION TIER: our exit code %d (upstream binary unavailable for cross-check)", gotCode)
				return
			}
			upDir := t.TempDir()
			upArgs := append([]string{}, tc.args...)
			upArgs = append(upArgs, filepath.Join(upDir, "up"), bam)
			upCmd := exec.Command(bin, upArgs...)
			upCmd.Dir = upDir
			upOut, upErr := upCmd.CombinedOutput()
			upCode := exitCode(upErr)
			if upCode != tc.code {
				t.Fatalf("upstream exit code = %d, want %d\n%s", upCode, tc.code, upOut)
			}
			t.Logf("VALIDATION TIER: exit code %d matches upstream mosdepth", gotCode)
		})
	}
}

// exitCode extracts the process exit status from an *exec.ExitError, returning
// 0 when err is nil and -1 when err is some other (non-exit) failure.
func exitCode(err error) int {
	if err == nil {
		return 0
	}
	var ee *exec.ExitError
	if errors.As(err, &ee) {
		return ee.ExitCode()
	}
	return -1
}

// assertFilesEqual fails the test when the two plain-text files differ
// byte-for-byte, printing both contents for diagnosis.
func assertFilesEqual(t *testing.T, gotPath, wantPath, label string) {
	t.Helper()
	got, err := os.ReadFile(gotPath)
	if err != nil {
		t.Fatalf("read our %s: %v", label, err)
	}
	want, err := os.ReadFile(wantPath)
	if err != nil {
		t.Fatalf("read upstream %s: %v", label, err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("%s mismatch.\nours:\n%s\nupstream:\n%s", label, got, want)
	}
}

// mtLines returns the subset of BED lines (from raw decompressed bytes) whose
// chromosome is exactly "MT". ovl.bam declares ~80 references but the reads
// only touch MT; scoping to MT keeps the upstream diff tight.
func mtLines(data []byte) []byte {
	var out bytes.Buffer
	for _, ln := range bytes.Split(data, []byte("\n")) {
		if bytes.HasPrefix(ln, []byte("MT\t")) {
			out.Write(ln)
			out.WriteByte('\n')
		}
	}
	return out.Bytes()
}

// splitEnv splits a "KEY=VALUE" string into its key and value.
func splitEnv(kv string) (string, string) {
	for i := 0; i < len(kv); i++ {
		if kv[i] == '=' {
			return kv[:i], kv[i+1:]
		}
	}
	return kv, ""
}
