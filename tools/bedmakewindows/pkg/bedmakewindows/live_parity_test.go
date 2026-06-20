package bedmakewindows

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

// Live-upstream parity tests: they run the real upstream `bedtools` binary
// from the vendored submodule and compare its `makewindows` output, byte for
// byte, against this port. They exercise the default (no -i, ID_NONE / BED3)
// path plus every -i value (src, winnum, srcwinnum) over both a BED file and
// a genome file. The genome cases cover the regression where `-i src` /
// `-i srcwinnum` must annotate windows with the chromosome name (upstream
// builds each genome interval with name == chrom). They t.Fatalf (never
// t.Skip) when the upstream binary is absent.

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

func fixture(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join("..", "..", "testdata", "parity", name)
}

func oursWindows(t *testing.T, ivs []Interval, opts Options) []byte {
	t.Helper()
	var out, warn bytes.Buffer
	if _, err := MakeWindows(ivs, &out, &warn, opts); err != nil {
		t.Fatalf("MakeWindows: %v", err)
	}
	return out.Bytes()
}

// assertBEDParity compares against upstream `makewindows -b <bed> -n <count>`
// with the given -i naming (empty iArg means the default, no -i).
func assertBEDParity(t *testing.T, count int, iArg string, naming Naming) {
	t.Helper()
	bin := upstreamBedtools(t)
	bed := fixture(t, "named.bed")
	args := []string{"makewindows", "-b", bed, "-n", strconv.Itoa(count)}
	if iArg != "" {
		args = append(args, "-i", iArg)
	}
	want, err := exec.Command(bin, args...).Output()
	if err != nil {
		t.Fatalf("upstream makewindows %v: %v", args, err)
	}
	f, err := os.Open(bed)
	if err != nil {
		t.Fatalf("open bed: %v", err)
	}
	defer f.Close()
	ivs, err := FromBED(f)
	if err != nil {
		t.Fatalf("FromBED: %v", err)
	}
	got := oursWindows(t, ivs, Options{Count: count, Naming: naming})
	if !bytes.Equal(got, want) {
		t.Fatalf("BED parity mismatch (-i %q)\nupstream:\n%s\nours:\n%s", iArg, want, got)
	}
}

func assertGenomeParity(t *testing.T, count int, iArg string, naming Naming) {
	t.Helper()
	bin := upstreamBedtools(t)
	gen := fixture(t, "sizes.genome")
	args := []string{"makewindows", "-g", gen, "-n", strconv.Itoa(count)}
	if iArg != "" {
		args = append(args, "-i", iArg)
	}
	want, err := exec.Command(bin, args...).Output()
	if err != nil {
		t.Fatalf("upstream makewindows %v: %v", args, err)
	}
	f, err := os.Open(gen)
	if err != nil {
		t.Fatalf("open genome: %v", err)
	}
	defer f.Close()
	ivs, err := FromGenome(f)
	if err != nil {
		t.Fatalf("FromGenome: %v", err)
	}
	got := oursWindows(t, ivs, Options{Count: count, Naming: naming})
	if !bytes.Equal(got, want) {
		t.Fatalf("genome parity mismatch (-i %q)\nupstream:\n%s\nours:\n%s", iArg, want, got)
	}
}

func TestLiveParity_BED_Default(t *testing.T)   { assertBEDParity(t, 3, "", NoName) }
func TestLiveParity_BED_Src(t *testing.T)       { assertBEDParity(t, 3, "src", NameSrc) }
func TestLiveParity_BED_WinNum(t *testing.T)    { assertBEDParity(t, 3, "winnum", NameWinNum) }
func TestLiveParity_BED_SrcWinNum(t *testing.T) { assertBEDParity(t, 3, "srcwinnum", NameSrcWinNum) }
func TestLiveParity_Genome_Src(t *testing.T)    { assertGenomeParity(t, 2, "src", NameSrc) }
func TestLiveParity_Genome_SrcWinNum(t *testing.T) {
	assertGenomeParity(t, 2, "srcwinnum", NameSrcWinNum)
}
