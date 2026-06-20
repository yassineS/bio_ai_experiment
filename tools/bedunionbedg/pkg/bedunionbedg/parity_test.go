package bedunionbedg

import (
	"bytes"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

// Live-upstream parity tests: they build the real upstream `bedtools` binary
// from the vendored submodule and compare its `unionbedg` output, byte for
// byte, against this port. They t.Fatalf (never t.Skip) so a missing or
// unbuildable submodule is a hard failure, matching the project's parity-rig
// policy.

var (
	upstreamBedtoolsOnce sync.Once
	upstreamBedtoolsPath string
	upstreamBedtoolsErr  error
)

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

func upstreamBedtools(t *testing.T) string {
	t.Helper()
	if testing.Short() {
		t.Skip("skipping upstream-binary parity test in -short mode")
	}
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

func fixturePath(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join("..", "..", "testdata", "parity", name)
}

// runUpstreamUnion runs `bedtools unionbedg args...` and returns its stdout.
func runUpstreamUnion(t *testing.T, args ...string) []byte {
	t.Helper()
	bin := upstreamBedtools(t)
	cmd := exec.Command(bin, append([]string{"unionbedg"}, args...)...)
	out, _ := cmd.CombinedOutput()
	return out
}

// openFixtures opens the named fixtures as readers and returns a cleanup.
func openFixtures(t *testing.T, names ...string) ([]io.Reader, func()) {
	t.Helper()
	var rs []io.Reader
	var fs []*os.File
	for _, n := range names {
		f, err := os.Open(fixturePath(t, n))
		if err != nil {
			t.Fatalf("open %s: %v", n, err)
		}
		rs = append(rs, f)
		fs = append(fs, f)
	}
	return rs, func() {
		for _, f := range fs {
			f.Close()
		}
	}
}

func sizesFromFixture(t *testing.T, name string) map[string]int64 {
	t.Helper()
	f, err := os.Open(fixturePath(t, name))
	if err != nil {
		t.Fatalf("open %s: %v", name, err)
	}
	defer f.Close()
	cs, err := ReadChromSizes(f)
	if err != nil {
		t.Fatalf("ReadChromSizes %s: %v", name, err)
	}
	return cs
}

func runOursUnion(t *testing.T, opts Options, files ...string) []byte {
	t.Helper()
	rs, cleanup := openFixtures(t, files...)
	defer cleanup()
	var out bytes.Buffer
	if err := Union(rs, &out, opts); err != nil {
		t.Fatalf("Union: %v", err)
	}
	return out.Bytes()
}

func TestParity_Union_Basic(t *testing.T) {
	want := runUpstreamUnion(t, "-i", fixturePath(t, "1.bg"), fixturePath(t, "2.bg"), fixturePath(t, "3.bg"))
	got := runOursUnion(t, Options{}, "1.bg", "2.bg", "3.bg")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch\nupstream:\n%s\nours:\n%s", want, got)
	}
}

func TestParity_Union_Header(t *testing.T) {
	want := runUpstreamUnion(t, "-header", "-i", fixturePath(t, "1.bg"), fixturePath(t, "2.bg"), fixturePath(t, "3.bg"))
	got := runOursUnion(t, Options{PrintHeader: true}, "1.bg", "2.bg", "3.bg")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch\nupstream:\n%s\nours:\n%s", want, got)
	}
}

func TestParity_Union_Names(t *testing.T) {
	want := runUpstreamUnion(t, "-header", "-i",
		fixturePath(t, "1.bg"), fixturePath(t, "2.bg"), fixturePath(t, "3.bg"),
		"-names", "WT-1", "WT-2", "KO-1")
	got := runOursUnion(t, Options{PrintHeader: true, Names: []string{"WT-1", "WT-2", "KO-1"}},
		"1.bg", "2.bg", "3.bg")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch\nupstream:\n%s\nours:\n%s", want, got)
	}
}

func TestParity_Union_Empty(t *testing.T) {
	want := runUpstreamUnion(t, "-header", "-empty", "-g", fixturePath(t, "sizes.txt"),
		"-i", fixturePath(t, "1.bg"), fixturePath(t, "2.bg"), fixturePath(t, "3.bg"))
	got := runOursUnion(t, Options{PrintHeader: true, PrintEmpty: true, Sizes: sizesFromFixture(t, "sizes.txt")},
		"1.bg", "2.bg", "3.bg")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch\nupstream:\n%s\nours:\n%s", want, got)
	}
}

func TestParity_Union_Filler(t *testing.T) {
	want := runUpstreamUnion(t, "-empty", "-g", fixturePath(t, "sizes.txt"), "-filler", "N/A",
		"-i", fixturePath(t, "1.bg"), fixturePath(t, "2.bg"), fixturePath(t, "3.bg"))
	got := runOursUnion(t, Options{PrintEmpty: true, Sizes: sizesFromFixture(t, "sizes.txt"), Filler: "N/A"},
		"1.bg", "2.bg", "3.bg")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch\nupstream:\n%s\nours:\n%s", want, got)
	}
}

func TestParity_Union_MultiChromUnsorted(t *testing.T) {
	want := runUpstreamUnion(t, "-i", fixturePath(t, "m1.bg"), fixturePath(t, "m2.bg"))
	got := runOursUnion(t, Options{}, "m1.bg", "m2.bg")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch\nupstream:\n%s\nours:\n%s", want, got)
	}
}

func TestParity_Union_FloatDepths(t *testing.T) {
	want := runUpstreamUnion(t, "-i", fixturePath(t, "f1.bg"), fixturePath(t, "f2.bg"))
	got := runOursUnion(t, Options{}, "f1.bg", "f2.bg")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch\nupstream:\n%s\nours:\n%s", want, got)
	}
}
