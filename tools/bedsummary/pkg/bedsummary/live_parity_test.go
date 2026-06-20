package bedsummary

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

// Live-upstream parity tests run the real vendored `bedtools summary` binary
// and compare its output, byte for byte, against this port. They exercise the
// full upstream column set (chrom_length / the genome-fraction columns) and the
// required -g genome file — the features the port previously lacked. They
// t.Fatalf (never t.Skip) when the upstream binary is absent.

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

func assertLiveSummaryParity(t *testing.T, bedFile, genomeFile string) {
	t.Helper()
	bin := upstreamBedtools(t)
	bedPath := fixturePath(t, bedFile)
	genPath := fixturePath(t, genomeFile)

	want, err := exec.Command(bin, "summary", "-i", bedPath, "-g", genPath).Output()
	if err != nil {
		t.Fatalf("upstream summary: %v", err)
	}

	bedData, err := os.ReadFile(bedPath)
	if err != nil {
		t.Fatalf("read %s: %v", bedFile, err)
	}
	genData, err := os.ReadFile(genPath)
	if err != nil {
		t.Fatalf("read %s: %v", genomeFile, err)
	}
	g, err := ParseGenome(bytes.NewReader(genData))
	if err != nil {
		t.Fatalf("ParseGenome: %v", err)
	}

	var got bytes.Buffer
	if err := Run(bytes.NewReader(bedData), g, &got, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("summary parity mismatch (%s, %s)\nupstream:\n%s\nours:\n%s",
			bedFile, genomeFile, want, got.Bytes())
	}
}

// TestLiveParity_Summary_MultiChrom proves byte-for-byte parity on a small
// multi-chromosome input whose genome file includes a chromosome with no
// intervals (exercising the -1 default row) and the per-row trailing tab.
func TestLiveParity_Summary_MultiChrom(t *testing.T) {
	assertLiveSummaryParity(t, "multi.bed", "multi.genome")
}

// TestLiveParity_Summary_Large proves byte-for-byte parity on a 4000-interval,
// 5-chromosome dataset (fraction columns, fixed-9 precision, genome ordering).
func TestLiveParity_Summary_Large(t *testing.T) {
	assertLiveSummaryParity(t, "big.bed", "big.genome")
}
