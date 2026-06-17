package bedannotate

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

// Live-upstream parity tests run the real vendored `bedtools annotate` binary
// over multi-file -files inputs and compare its output, byte for byte, against
// this port. They lock in the two regressions the parity pipeline found:
//   - no spurious "# <file>" header (upstream emits a header only with -names);
//   - upstream's record order (grouped per chromosome, then by UCSC bin).
// They t.Fatalf (never t.Skip) when the upstream binary is absent.

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

// assertLiveAnnotateParity runs upstream `bedtools annotate` and the port over
// the same fixtures with the given extra upstream args, mapping those args to
// the equivalent Options.
func assertLiveAnnotateParity(t *testing.T, aFile string, bFiles []string, opts Options, upstreamArgs ...string) {
	t.Helper()
	bin := upstreamBedtools(t)
	aPath := fixturePath(t, aFile)

	args := []string{"annotate", "-i", aPath, "-files"}
	bPaths := make([]string, len(bFiles))
	for i, f := range bFiles {
		bPaths[i] = fixturePath(t, f)
	}
	args = append(args, bPaths...)
	args = append(args, upstreamArgs...)

	want, err := exec.Command(bin, args...).Output()
	if err != nil {
		t.Fatalf("upstream annotate %v: %v", upstreamArgs, err)
	}

	aData, err := os.ReadFile(aPath)
	if err != nil {
		t.Fatalf("read %s: %v", aFile, err)
	}
	bRs := make([]io.Reader, len(bPaths))
	for i, p := range bPaths {
		d, err := os.ReadFile(p)
		if err != nil {
			t.Fatalf("read %s: %v", bFiles[i], err)
		}
		bRs[i] = bytes.NewReader(d)
	}

	var got bytes.Buffer
	if _, err := Run(bytes.NewReader(aData), bRs, &got, opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("annotate parity mismatch (%s %v)\nupstream:\n%s\nours:\n%s",
			aFile, upstreamArgs, want, got.Bytes())
	}
}

// TestLiveParity_Annotate_NoHeaderDefault proves no header is emitted (and the
// record order matches) on a multi-file default run.
func TestLiveParity_Annotate_NoHeaderDefault(t *testing.T) {
	assertLiveAnnotateParity(t, "order_a.bed", []string{"order_b1.bed", "order_b2.bed"}, Options{})
}

// TestLiveParity_Annotate_Counts proves -counts parity.
func TestLiveParity_Annotate_Counts(t *testing.T) {
	assertLiveAnnotateParity(t, "order_a.bed", []string{"order_b1.bed", "order_b2.bed"},
		Options{Mode: ModeCounts}, "-counts")
}

// TestLiveParity_Annotate_BothNamed proves -both with -names (header) parity.
func TestLiveParity_Annotate_BothNamed(t *testing.T) {
	assertLiveAnnotateParity(t, "order_a.bed", []string{"order_b1.bed", "order_b2.bed"},
		Options{Mode: ModeBoth, Names: []string{"b1", "b2"}}, "-both", "-names", "b1", "b2")
}

// TestLiveParity_Annotate_LargeMultiChrom stresses the chromosome/bin ordering
// over a 500-interval, 4-chromosome (chr1, chr2, chr10, chrX) dataset.
func TestLiveParity_Annotate_LargeMultiChrom(t *testing.T) {
	assertLiveAnnotateParity(t, "multi_a.bed", []string{"multi_b1.bed", "multi_b2.bed"}, Options{})
}

// TestLiveParity_Annotate_SameStrand proves -s parity over stranded fixtures.
func TestLiveParity_Annotate_SameStrand(t *testing.T) {
	assertLiveAnnotateParity(t, "multi_a.bed", []string{"multi_b1.bed", "multi_b2.bed"},
		Options{SameStrand: true}, "-s")
}

// TestLiveParity_Annotate_OppositeStrand proves -S parity over stranded
// fixtures.
func TestLiveParity_Annotate_OppositeStrand(t *testing.T) {
	assertLiveAnnotateParity(t, "multi_a.bed", []string{"multi_b1.bed", "multi_b2.bed"},
		Options{OppositeStrand: true}, "-S")
}
