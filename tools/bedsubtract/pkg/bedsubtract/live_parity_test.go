package bedsubtract

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

// Live-upstream parity tests run the real vendored `bedtools subtract` binary
// and compare its output, byte for byte, against this port. They focus on the
// newly added reciprocal (-r) flag: the recip_* fixtures contain B intervals
// whose overlap covers enough of A but not enough of B (and vice versa), so -r
// changes the result versus plain -f. They t.Fatalf (never t.Skip) when the
// upstream binary is absent.

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

func subtractCLIArgs(aFile, bFile string, opts Options) []string {
	args := []string{"subtract", "-a", aFile, "-b", bFile}
	if opts.RemoveEntire {
		args = append(args, "-A")
	}
	if opts.RemoveSum {
		args = append(args, "-N")
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
	if opts.MinFraction > 0 {
		args = append(args, "-f", strconv.FormatFloat(opts.MinFraction, 'g', -1, 64))
	}
	return args
}

func assertLiveSubtractParity(t *testing.T, aFile, bFile string, opts Options) {
	t.Helper()
	bin := upstreamBedtools(t)
	aPath := fixturePath(t, aFile)
	bPath := fixturePath(t, bFile)

	want, err := exec.Command(bin, subtractCLIArgs(aPath, bPath, opts)...).Output()
	if err != nil {
		t.Fatalf("upstream subtract: %v", err)
	}

	aData, err := os.ReadFile(aPath)
	if err != nil {
		t.Fatalf("read %s: %v", aFile, err)
	}
	bData, err := os.ReadFile(bPath)
	if err != nil {
		t.Fatalf("read %s: %v", bFile, err)
	}

	var got bytes.Buffer
	if _, err := Subtract(bytes.NewReader(aData), bytes.NewReader(bData), &got, opts); err != nil {
		t.Fatalf("Subtract: %v", err)
	}
	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("subtract parity mismatch (%s vs %s, %+v)\nupstream:\n%s\nours:\n%s",
			aFile, bFile, opts, want, got.Bytes())
	}
}

// TestLiveParity_Subtract_Reciprocal proves -r matches upstream byte-for-byte:
// a B record that covers enough of A but too little of B must NOT trigger
// subtraction under -r (whereas plain -f would subtract it).
func TestLiveParity_Subtract_Reciprocal(t *testing.T) {
	assertLiveSubtractParity(t, "recip_a.bed", "recip_b.bed", Options{MinFraction: 0.1, Reciprocal: true})
}

// TestLiveParity_Subtract_ReciprocalHigh proves -r with a higher threshold.
func TestLiveParity_Subtract_ReciprocalHigh(t *testing.T) {
	assertLiveSubtractParity(t, "recip_a.bed", "recip_b.bed", Options{MinFraction: 0.5, Reciprocal: true})
}

// TestLiveParity_Subtract_FractionOnly is the same fixtures without -r, proving
// the reciprocal flag genuinely changes the outcome.
func TestLiveParity_Subtract_FractionOnly(t *testing.T) {
	assertLiveSubtractParity(t, "recip_a.bed", "recip_b.bed", Options{MinFraction: 0.1})
}

// TestLiveParity_Subtract_ReciprocalNoFraction proves -r with no -f behaves
// like plain subtract (fraction defaults to ~0 on both sides).
func TestLiveParity_Subtract_ReciprocalNoFraction(t *testing.T) {
	assertLiveSubtractParity(t, "recip_a.bed", "recip_b.bed", Options{Reciprocal: true})
}
