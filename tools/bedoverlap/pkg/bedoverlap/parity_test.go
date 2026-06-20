package bedoverlap

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

// Live-upstream parity tests: they build the real upstream `bedtools` binary
// from the vendored submodule and compare its `overlap` output, byte for byte,
// against this port. They t.Fatalf (never t.Skip) so a missing or unbuildable
// submodule is a hard failure, matching the project's parity-rig policy.

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

// runUpstream runs `bedtools overlap -i <fixture> -cols <spec>` and returns
// combined stdout/stderr.
func runUpstream(t *testing.T, fixture, spec string) []byte {
	t.Helper()
	bin := upstreamBedtools(t)
	cmd := exec.Command(bin, "overlap", "-i", fixturePath(t, fixture), "-cols", spec)
	out, _ := cmd.CombinedOutput()
	return out
}

// runOurs runs this port's Overlap over the same fixture and column spec.
func runOurs(t *testing.T, fixture, spec string) []byte {
	t.Helper()
	cols, err := ParseCols(spec)
	if err != nil {
		t.Fatalf("ParseCols(%q): %v", spec, err)
	}
	data, err := os.ReadFile(fixturePath(t, fixture))
	if err != nil {
		t.Fatalf("read fixture %s: %v", fixture, err)
	}
	var out bytes.Buffer
	if err := Overlap(bytes.NewReader(data), &out, cols); err != nil {
		t.Fatalf("Overlap(%s): %v", fixture, err)
	}
	return out.Bytes()
}

func assertParity(t *testing.T, fixture, spec string) {
	t.Helper()
	want := runUpstream(t, fixture, spec)
	got := runOurs(t, fixture, spec)
	if !bytes.Equal(got, want) {
		t.Fatalf("parity mismatch for %s -cols %s\nupstream:\n%s\nours:\n%s", fixture, spec, want, got)
	}
}

func TestParity_Overlap_Window(t *testing.T)   { assertParity(t, "window.txt", "2,3,6,7") }
func TestParity_Overlap_Touching(t *testing.T) { assertParity(t, "touching.txt", "2,3,6,7") }
func TestParity_Overlap_Nested(t *testing.T)   { assertParity(t, "nested.txt", "2,3,6,7") }
func TestParity_Overlap_Disjoint(t *testing.T) { assertParity(t, "disjoint.txt", "2,3,6,7") }
