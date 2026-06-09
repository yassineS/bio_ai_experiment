package mosdepth

import (
	"bytes"
	"io"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
	"time"
)

// TestD4_UpstreamBinaryParity proves that our D4 writer produces a
// `.per-base.d4` that is byte-identical to the one the REAL upstream mosdepth
// binary emits for the same BAM.
//
// It runs the official `mosdepth_d4` release binary (a build of mosdepth with
// D4 support linked in) on a small fixture BAM, runs our implementation on the
// same input, and asserts the two `.per-base.d4` files are byte-for-byte equal.
//
// The upstream binary is downloaded from the GitHub release on first use and
// cached under the OS temp dir; the download is gated behind a reachability
// check so a genuinely offline machine fails with a clear message rather than
// hanging. This test never silently skips when it can reach the network.
func TestD4_UpstreamBinaryParity(t *testing.T) {
	if runtime.GOOS != "linux" || runtime.GOARCH != "amd64" {
		t.Skipf("upstream mosdepth_d4 release binary is only published for linux/amd64 (have %s/%s)", runtime.GOOS, runtime.GOARCH)
	}

	bin := ensureMosdepthD4Binary(t)

	fixtures := fixtureDir(t)
	bam := filepath.Join(fixtures, "d4-small.bam")
	if _, err := os.Stat(bam); err != nil {
		t.Fatalf("fixture BAM missing: %v", err)
	}

	// Run the upstream binary. It writes <prefix>.per-base.d4 alongside the
	// summary/dist files. -x (fast-mode) matches our writer's pipeline.
	upDir := t.TempDir()
	upPrefix := filepath.Join(upDir, "up")
	cmd := exec.Command(bin, "--d4", "-x", upPrefix, bam)
	cmd.Dir = upDir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream mosdepth_d4 failed: %v\n%s", err, out)
	}
	upD4 := upPrefix + ".per-base.d4"
	upBytes, err := os.ReadFile(upD4)
	if err != nil {
		t.Fatalf("read upstream d4: %v", err)
	}

	// Run our implementation on the same BAM.
	ourDir := t.TempDir()
	ourPrefix := filepath.Join(ourDir, "our")
	if err := OpenAndRun(bam, Options{
		Prefix:   ourPrefix,
		FastMode: true,
		D4Output: true,
	}); err != nil {
		t.Fatalf("OpenAndRun: %v", err)
	}
	ourD4 := ourPrefix + ".per-base.d4"
	ourBytes, err := os.ReadFile(ourD4)
	if err != nil {
		t.Fatalf("read our d4: %v", err)
	}

	if len(ourBytes) != len(upBytes) {
		t.Fatalf("D4 size mismatch: ours=%d upstream=%d bytes", len(ourBytes), len(upBytes))
	}
	if !bytes.Equal(ourBytes, upBytes) {
		// Report the first differing byte to aid debugging.
		for i := range ourBytes {
			if ourBytes[i] != upBytes[i] {
				t.Fatalf("D4 byte mismatch at offset %d: ours=0x%02x upstream=0x%02x (size %d, identical)",
					i, ourBytes[i], upBytes[i], len(ourBytes))
			}
		}
	}
	t.Logf("D4 byte-identical to upstream: %d bytes", len(ourBytes))

	// Cross-check: our reader decodes the upstream file to the same per-base
	// depths it decodes from ours, confirming the encoding is semantically
	// correct and not merely coincidentally equal byte-wise.
	ur, err := openD4Reader(upD4)
	if err != nil {
		t.Fatalf("openD4Reader(upstream): %v", err)
	}
	or, err := openD4Reader(ourD4)
	if err != nil {
		t.Fatalf("openD4Reader(ours): %v", err)
	}
	for _, c := range []string{"chr1", "chr2"} {
		ud, err := ur.chromDepths(c)
		if err != nil {
			t.Fatalf("upstream chromDepths(%s): %v", c, err)
		}
		od, err := or.chromDepths(c)
		if err != nil {
			t.Fatalf("our chromDepths(%s): %v", c, err)
		}
		if !equalInt32(ud, od) {
			t.Fatalf("decoded depths differ for %s", c)
		}
	}
}

// mosdepthD4URL is the GitHub release asset for the D4-enabled mosdepth build.
const mosdepthD4URL = "https://github.com/brentp/mosdepth/releases/download/v0.3.10/mosdepth_d4"

var ensureMosdepthD4Once sync.Once

// ensureMosdepthD4Binary returns a path to an executable upstream mosdepth_d4
// binary, downloading and caching it on first use. It t.Fatalf()s (never
// silently skips) if the binary cannot be obtained while the network is
// reachable.
func ensureMosdepthD4Binary(t *testing.T) string {
	t.Helper()

	// Allow an override so CI / local runs can point at a pre-fetched binary.
	if p := os.Getenv("MOSDEPTH_D4_BIN"); p != "" {
		if _, err := os.Stat(p); err == nil {
			return p
		}
		t.Fatalf("MOSDEPTH_D4_BIN=%q does not exist", p)
	}

	cache := filepath.Join(os.TempDir(), "mosdepth_d4_v0.3.10")
	if fi, err := os.Stat(cache); err == nil && fi.Size() > 0 {
		return cache
	}

	// Gate the download behind a reachability check: if the network is truly
	// unavailable, fail loudly with an actionable message instead of hanging.
	if !networkReachable() {
		t.Fatalf("cannot download upstream mosdepth_d4 binary: network unreachable. " +
			"Set MOSDEPTH_D4_BIN to a local copy of the binary to run this parity test.")
	}

	var dlErr error
	ensureMosdepthD4Once.Do(func() {
		dlErr = downloadFile(mosdepthD4URL, cache)
	})
	if dlErr != nil {
		// Re-check in case a concurrent run populated the cache.
		if fi, err := os.Stat(cache); err == nil && fi.Size() > 0 {
			return cache
		}
		t.Fatalf("download mosdepth_d4: %v", dlErr)
	}
	if err := os.Chmod(cache, 0o755); err != nil {
		t.Fatalf("chmod mosdepth_d4: %v", err)
	}
	return cache
}

// networkReachable reports whether github.com is reachable over HTTPS within a
// short timeout.
func networkReachable() bool {
	client := &http.Client{Timeout: 8 * time.Second}
	resp, err := client.Head("https://github.com")
	if err != nil {
		return false
	}
	_ = resp.Body.Close()
	return true
}

// downloadFile fetches url into dst with a bounded retry/backoff, following
// redirects (the default http.Client behaviour).
func downloadFile(url, dst string) error {
	var lastErr error
	backoff := []time.Duration{0, 2 * time.Second, 4 * time.Second, 8 * time.Second, 16 * time.Second}
	for _, d := range backoff {
		if d > 0 {
			time.Sleep(d)
		}
		client := &http.Client{Timeout: 5 * time.Minute}
		resp, err := client.Get(url)
		if err != nil {
			lastErr = err
			continue
		}
		if resp.StatusCode != http.StatusOK {
			_ = resp.Body.Close()
			lastErr = &httpStatusError{code: resp.StatusCode}
			continue
		}
		tmp := dst + ".part"
		f, err := os.Create(tmp)
		if err != nil {
			_ = resp.Body.Close()
			return err
		}
		_, err = io.Copy(f, resp.Body)
		_ = resp.Body.Close()
		_ = f.Close()
		if err != nil {
			lastErr = err
			continue
		}
		if err := os.Rename(tmp, dst); err != nil {
			return err
		}
		return nil
	}
	return lastErr
}

type httpStatusError struct{ code int }

func (e *httpStatusError) Error() string {
	return "unexpected HTTP status " + http.StatusText(e.code)
}
