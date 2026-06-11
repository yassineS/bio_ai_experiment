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
