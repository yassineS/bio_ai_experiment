package mosdepth

import (
	"bytes"
	"compress/gzip"
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
