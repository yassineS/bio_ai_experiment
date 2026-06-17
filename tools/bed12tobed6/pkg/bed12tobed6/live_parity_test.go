package bed12tobed6

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

// Live-upstream parity tests: they run the real upstream `bedtools` binary
// from the vendored submodule and compare its `bed12tobed6` output, byte for
// byte, against this port. The fixtures carry non-zero scores so they
// exercise the score-propagation bug (the port previously emitted "0" in the
// score column instead of carrying the parent record's score onto each
// emitted BED6 block). They t.Fatalf (never t.Skip) when the upstream binary
// is absent, matching the project's parity-rig policy.

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
		t.Fatalf("upstream bedtools binary not found at reference_code/bedtools/bin/bedtools " +
			"(run `git submodule update --init reference_code/bedtools` and build it)")
	}
	return upstreamBedtoolsPath
}

func fixturePath(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join("..", "..", "testdata", "parity", name)
}

func assertLiveParity(t *testing.T, fixture string, opts Options, upstreamArgs ...string) {
	t.Helper()
	bin := upstreamBedtools(t)
	in := fixturePath(t, fixture)

	args := append([]string{"bed12tobed6", "-i", in}, upstreamArgs...)
	want, err := exec.Command(bin, args...).Output()
	if err != nil {
		t.Fatalf("upstream bed12tobed6 %v: %v", upstreamArgs, err)
	}

	data, err := os.ReadFile(in)
	if err != nil {
		t.Fatalf("read fixture %s: %v", fixture, err)
	}
	var got bytes.Buffer
	if _, err := Convert(bytes.NewReader(data), &got, opts); err != nil {
		t.Fatalf("Convert: %v", err)
	}
	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("parity mismatch for %s (%v)\nupstream:\n%s\nours:\n%s", fixture, upstreamArgs, want, got.Bytes())
	}
}

// TestLiveParity_ScorePropagation proves the parent record's score (column 5)
// is carried onto each emitted BED6 block, matching upstream.
func TestLiveParity_ScorePropagation(t *testing.T) {
	assertLiveParity(t, "scored.bed", Options{})
}

// TestLiveParity_ScoredNumbered proves -n still overrides the score with the
// (strand-aware) 1-based block number, byte-for-byte against upstream.
func TestLiveParity_ScoredNumbered(t *testing.T) {
	assertLiveParity(t, "scored.bed", Options{NumberBlocks: true}, "-n")
}
