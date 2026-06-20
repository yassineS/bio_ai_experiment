package bedsplit

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
// from the vendored submodule and compare its `split -a size` output, byte
// for byte, against this port — both the stdout manifest AND every emitted
// shard file. The fixture uses records of distinct sizes so the
// size-descending sort order is fully determined (upstream's std::sort is
// unstable for equal-length records). They cover the bin-packing heuristic
// across several -n values, including -n > record-count. They t.Fatalf
// (never t.Skip) when the upstream binary is absent.

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

func fixtureAbs(t *testing.T, name string) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("abs path: %v", err)
	}
	return abs
}

// readShards reads <prefix>.NNNNN.bed shard files from dir (in numeric order)
// keyed by their base filename.
func readShards(t *testing.T, dir string) map[string][]byte {
	t.Helper()
	entries, err := os.ReadDir(dir)
	if err != nil {
		t.Fatalf("readdir %s: %v", dir, err)
	}
	out := make(map[string][]byte)
	for _, e := range entries {
		if e.IsDir() || !strings.HasSuffix(e.Name(), ".bed") {
			continue
		}
		data, err := os.ReadFile(filepath.Join(dir, e.Name()))
		if err != nil {
			t.Fatalf("read shard %s: %v", e.Name(), err)
		}
		out[e.Name()] = data
	}
	return out
}

func assertSizeParity(t *testing.T, fixture string, n int) {
	t.Helper()
	bin := upstreamBedtools(t)
	in := fixtureAbs(t, fixture)

	// Upstream writes shards into its working directory using the prefix.
	upDir := t.TempDir()
	cmd := exec.Command(bin, "split", "-i", in, "-n", strconv.Itoa(n), "-p", "s", "-a", "size")
	cmd.Dir = upDir
	wantManifest, err := cmd.Output()
	if err != nil {
		t.Fatalf("upstream split -n %d: %v", n, err)
	}

	// Our port: run with a prefix inside our own temp dir, then strip the
	// directory prefix from the manifest so the filenames line up with
	// upstream's "s.NNNNN.bed".
	ourDir := t.TempDir()
	prefix := filepath.Join(ourDir, "s")
	data, err := os.ReadFile(in)
	if err != nil {
		t.Fatalf("read fixture: %v", err)
	}
	var manifest bytes.Buffer
	if _, err := Split(bytes.NewReader(data), &manifest, Options{N: n, Prefix: prefix, Algorithm: AlgSize}); err != nil {
		t.Fatalf("Split: %v", err)
	}
	gotManifest := strings.ReplaceAll(manifest.String(), ourDir+string(os.PathSeparator), "")

	if gotManifest != string(wantManifest) {
		t.Fatalf("manifest mismatch (-n %d)\nupstream:\n%s\nours:\n%s", n, wantManifest, gotManifest)
	}

	wantShards := readShards(t, upDir)
	gotShards := readShards(t, ourDir)
	if len(wantShards) != len(gotShards) {
		t.Fatalf("shard count mismatch (-n %d): upstream=%d ours=%d", n, len(wantShards), len(gotShards))
	}
	for name, want := range wantShards {
		got, ok := gotShards[name]
		if !ok {
			t.Fatalf("missing shard %s (-n %d)", name, n)
		}
		if !bytes.Equal(got, want) {
			t.Fatalf("shard %s mismatch (-n %d)\nupstream:\n%s\nours:\n%s", name, n, want, got)
		}
	}
}

func TestLiveParity_Size_N1(t *testing.T)  { assertSizeParity(t, "distinct_sizes.bed", 1) }
func TestLiveParity_Size_N2(t *testing.T)  { assertSizeParity(t, "distinct_sizes.bed", 2) }
func TestLiveParity_Size_N3(t *testing.T)  { assertSizeParity(t, "distinct_sizes.bed", 3) }
func TestLiveParity_Size_N5(t *testing.T)  { assertSizeParity(t, "distinct_sizes.bed", 5) }
func TestLiveParity_Size_N10(t *testing.T) { assertSizeParity(t, "distinct_sizes.bed", 10) }
