package bedfisher

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

// Live-upstream parity tests run the real vendored `bedtools fisher` binary and
// compare its full report, byte for byte, against this port. The overlap_*
// fixtures contain heavily self-overlapping A and B intervals (long records
// whose end coordinates are non-monotonic over a start-sorted slice). They
// exercise the overlap-counting bug the parity pipeline found: the port used to
// binary-search on ChromEnd over a start-sorted B slice — invalid, because a
// long B can start before A yet extend past A.Start — which dropped those pairs
// and under-counted overlaps (skewing n11/n12/n21/n22 and the p-values). These
// tests t.Fatalf (never t.Skip) when the upstream binary is absent.

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

func fixturePath(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join("..", "..", "testdata", "parity", name)
}

// runFisherCLIArgs builds the upstream argument list for `bedtools fisher`.
func runFisherCLIArgs(aFile, bFile, gFile string, opts Options) []string {
	args := []string{"fisher", "-a", aFile, "-b", bFile, "-g", gFile}
	if opts.MergeInputs {
		args = append(args, "-m")
	}
	if opts.SameStrand {
		args = append(args, "-s")
	}
	if opts.OppositeStrand {
		args = append(args, "-S")
	}
	if opts.Reciprocal {
		args = append(args, "-r")
	}
	if opts.FractionA > 0 {
		args = append(args, "-f", strconv.FormatFloat(opts.FractionA, 'g', -1, 64))
	}
	if opts.FractionB > 0 {
		args = append(args, "-F", strconv.FormatFloat(opts.FractionB, 'g', -1, 64))
	}
	return args
}

// assertLiveFisherParity runs upstream `bedtools fisher` and the port over the
// same fixtures and asserts byte-for-byte equality.
func assertLiveFisherParity(t *testing.T, aFile, bFile, gFile string, opts Options) {
	t.Helper()
	bin := upstreamBedtools(t)
	aPath := fixturePath(t, aFile)
	bPath := fixturePath(t, bFile)
	gPath := fixturePath(t, gFile)

	want, err := exec.Command(bin, runFisherCLIArgs(aPath, bPath, gPath, opts)...).Output()
	if err != nil {
		t.Fatalf("upstream fisher: %v", err)
	}

	aData, err := os.ReadFile(aPath)
	if err != nil {
		t.Fatalf("read %s: %v", aFile, err)
	}
	bData, err := os.ReadFile(bPath)
	if err != nil {
		t.Fatalf("read %s: %v", bFile, err)
	}
	gData, err := os.ReadFile(gPath)
	if err != nil {
		t.Fatalf("read %s: %v", gFile, err)
	}

	var got bytes.Buffer
	if _, err := Run(bytes.NewReader(aData), bytes.NewReader(bData), bytes.NewReader(gData), &got, opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("fisher parity mismatch (%s vs %s)\nupstream:\n%s\nours:\n%s", aFile, bFile, want, got.Bytes())
	}
}

// TestLiveParity_Fisher_HeavyOverlap proves the overlap counter matches upstream
// byte-for-byte on a dataset with thousands of self-overlapping intervals (the
// regression that previously under-counted overlaps).
func TestLiveParity_Fisher_HeavyOverlap(t *testing.T) {
	assertLiveFisherParity(t, "overlap_a.bed", "overlap_b.bed", "overlap.genome", Options{})
}

// TestLiveParity_Fisher_HeavyOverlapMerged proves the same dataset matches with
// -m (pre-merged A intervals).
func TestLiveParity_Fisher_HeavyOverlapMerged(t *testing.T) {
	assertLiveFisherParity(t, "overlap_a.bed", "overlap_b.bed", "overlap.genome", Options{MergeInputs: true})
}
