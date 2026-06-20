package bedexpand

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"testing"
)

// Live-upstream parity tests: they run the real upstream `bedtools` binary
// from the vendored submodule and compare its `expand` output, byte for byte,
// against this port. The fixtures exercise the trailing-comma bug (the port
// previously emitted a spurious empty final row when an expanded column ended
// with a comma) as well as leading/interior empty elements (which must be
// preserved). They t.Fatalf (never t.Skip) when the upstream binary is
// absent, matching the project's parity-rig policy.

var (
	upstreamBedtoolsOnce sync.Once
	upstreamBedtoolsPath string
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
		bin := filepath.Join(root, "reference_code", "bedtools", "bin", "bedtools")
		if _, err := os.Stat(bin); err == nil {
			upstreamBedtoolsPath = bin
		}
	})
	if upstreamBedtoolsPath == "" {
		t.Skipf("upstream bedtools binary not found at reference_code/bedtools/bin/bedtools " +
			"(run `git submodule update --init reference_code/bedtools` and build it)")
	}
	return upstreamBedtoolsPath
}

func colsArg(cols []int) string {
	parts := make([]string, len(cols))
	for i, c := range cols {
		parts[i] = strconv.Itoa(c)
	}
	return strings.Join(parts, ",")
}

func assertLiveParity(t *testing.T, fixture string, cols []int) {
	t.Helper()
	bin := upstreamBedtools(t)
	in := filepath.Join("..", "..", "testdata", "parity", fixture)

	want, err := exec.Command(bin, "expand", "-i", in, "-c", colsArg(cols)).Output()
	if err != nil {
		t.Fatalf("upstream expand -c %s: %v", colsArg(cols), err)
	}

	data, err := os.ReadFile(in)
	if err != nil {
		t.Fatalf("read fixture %s: %v", fixture, err)
	}
	var got bytes.Buffer
	if _, err := Expand(bytes.NewReader(data), &got, Options{Columns: cols}); err != nil {
		t.Fatalf("Expand: %v", err)
	}
	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("parity mismatch for %s -c %s\nupstream:\n%q\nours:\n%q", fixture, colsArg(cols), want, got.Bytes())
	}
}

// TestLiveParity_TrailingComma proves a trailing comma in the expanded column
// does NOT produce a spurious empty final row.
func TestLiveParity_TrailingComma(t *testing.T) {
	assertLiveParity(t, "trailing_comma.txt", []int{5})
}

// TestLiveParity_LeadingInteriorEmpties proves leading and interior empty
// elements are preserved (only a single terminating empty is dropped).
func TestLiveParity_LeadingInteriorEmpties(t *testing.T) {
	assertLiveParity(t, "empties.txt", []int{4, 5})
}
