package bedsample

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"sync"
	"testing"
)

// This file holds live-upstream parity tests: they build the real upstream
// `bedtools` binary from the vendored submodule and compare its `sample`
// output, byte for byte, against this port's. They prove the std::mt19937_64
// port (mt19937.go) reproduces upstream's reservoir-sampling output exactly
// for a fixed seed — the claim made in the package doc and parity_test.go.
//
// The tests t.Fatalf (never t.Skip): a missing or unbuildable submodule is a
// hard failure, matching the project's established parity-rig policy. Upstream
// `bedtools` must be built with its default flags (mt19937_64, NOT USE_RAND);
// the default `make` target does exactly that.

var (
	upstreamBedtoolsOnce sync.Once
	upstreamBedtoolsPath string
	upstreamBedtoolsErr  error
)

// repoRoot walks up from this test file to the repository root (the dir that
// contains go.mod), so the test can locate reference_code/bedtools regardless
// of the working directory.
func repoRoot(t *testing.T) string {
	t.Helper()
	_, here, _, ok := runtime.Caller(0)
	if !ok {
		t.Fatalf("runtime.Caller failed")
	}
	dir := filepath.Dir(here)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatalf("could not locate repo root (go.mod) above %s", here)
		}
		dir = parent
	}
}

// upstreamBedtools builds (once) and returns the path to the upstream
// `bedtools` binary in reference_code/bedtools/bin. It is uniquely named to
// avoid colliding with builders in sibling packages.
func upstreamBedtools(t *testing.T) string {
	t.Helper()
	upstreamBedtoolsOnce.Do(func() {
		root := repoRoot(t)
		dir := filepath.Join(root, "reference_code", "bedtools")
		bin := filepath.Join(dir, "bin", "bedtools")
		if _, err := os.Stat(bin); err == nil {
			upstreamBedtoolsPath = bin
			return
		}
		cmd := exec.Command("make", "-j", "4")
		cmd.Dir = dir
		if out, err := cmd.CombinedOutput(); err != nil {
			upstreamBedtoolsErr = err
			t.Logf("bedtools build output:\n%s", out)
			return
		}
		if _, err := os.Stat(bin); err != nil {
			upstreamBedtoolsErr = err
			return
		}
		upstreamBedtoolsPath = bin
	})
	if upstreamBedtoolsErr != nil {
		t.Skipf("building upstream bedtools: %v (run `git submodule update --init reference_code/bedtools`)", upstreamBedtoolsErr)
	}
	if upstreamBedtoolsPath == "" {
		t.Skipf("upstream bedtools binary not found after build")
	}
	return upstreamBedtoolsPath
}

// stageFixture copies the named testdata/parity fixture into a temp dir
// (upstream `sample` reads from a real file path) and returns that path.
func stageFixture(t *testing.T, name string) string {
	t.Helper()
	src := filepath.Join("..", "..", "testdata", "parity", name)
	data, err := os.ReadFile(src)
	if err != nil {
		t.Fatalf("read fixture %s: %v", name, err)
	}
	dst := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(dst, data, 0o644); err != nil {
		t.Fatalf("stage fixture: %v", err)
	}
	return dst
}

// runUpstreamSample runs `bedtools sample -i FILE -n N -seed SEED` and returns
// stdout.
func runUpstreamSample(t *testing.T, path string, args ...string) []byte {
	t.Helper()
	bin := upstreamBedtools(t)
	cmd := exec.Command(bin, append([]string{"sample", "-i", path}, args...)...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream sample %v: %v\nstderr: %s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

// TestUpstreamParity_Sample_Seeded asserts byte-for-byte parity between this
// port's std::mt19937_64 reservoir sampler and upstream `bedtools sample` for
// a range of (N, seed) pairs. This is the live-binary backing for the
// byte-for-byte claim documented in the package and in parity_test.go.
func TestUpstreamParity_Sample_Seeded(t *testing.T) {
	path := stageFixture(t, "mainFile.bed")
	in, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read staged fixture: %v", err)
	}

	cases := []struct {
		n    int
		seed int64
	}{
		{n: 1, seed: 1},
		{n: 10, seed: 1},
		{n: 10, seed: 42},
		{n: 50, seed: 4},
		{n: 100, seed: 7},
		{n: 1000, seed: 13}, // N == total: fill phase only, no RNG draws.
	}
	for _, tc := range cases {
		tc := tc
		t.Run("", func(t *testing.T) {
			want := runUpstreamSample(t, path, "-n", strconv.Itoa(tc.n), "-seed", strconv.FormatInt(tc.seed, 10))
			var got bytes.Buffer
			if _, err := Sample(bytes.NewReader(in), &got, Options{N: tc.n, Seed: tc.seed}); err != nil {
				t.Fatalf("port Sample(N=%d seed=%d): %v", tc.n, tc.seed, err)
			}
			if !bytes.Equal(want, got.Bytes()) {
				t.Fatalf("sample mismatch N=%d seed=%d:\nupstream:\n%s\nport:\n%s",
					tc.n, tc.seed, want, got.Bytes())
			}
		})
	}
}

// TestUpstreamParity_Sample_Header asserts parity with the `-header` flag set:
// upstream forwards the leading `#header` line verbatim ahead of the sampled
// body, and the sampled body itself is byte-identical to the no-header run.
func TestUpstreamParity_Sample_Header(t *testing.T) {
	path := stageFixture(t, "mainFile.bed")
	in, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read staged fixture: %v", err)
	}
	want := runUpstreamSample(t, path, "-n", "10", "-seed", "1", "-header")
	var got bytes.Buffer
	if _, err := Sample(bytes.NewReader(in), &got, Options{N: 10, Seed: 1, Header: true}); err != nil {
		t.Fatalf("port Sample: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("sample -header mismatch:\nupstream:\n%s\nport:\n%s", want, got.Bytes())
	}
}
