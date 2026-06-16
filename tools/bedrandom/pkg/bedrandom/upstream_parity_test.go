package bedrandom

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

// This file holds live-upstream parity tests: they build the real upstream
// `bedtools` binary from the vendored submodule and compare its `random`
// output, byte for byte, against this port's. They prove the std::mt19937_64
// port (mt19937.go), the rand_range rejection bound, and the genome-projection
// + redraw draw order reproduce upstream `bedtools random -seed N` exactly.
//
// The tests t.Fatalf (never t.Skip): a missing or unbuildable submodule is a
// hard failure, matching the project's established parity-rig policy. Upstream
// `bedtools` must be built with its default flags (mt19937_64, NOT USE_RAND);
// the default `make` target does exactly that.

var (
	upstreamBedtoolsOnce sync.Once
	upstreamBedtoolsPath string
	upstreamBedtoolsErr  error
)

// repoRoot walks up from this test file to the repository root (the dir that
// contains go.mod), so the test can locate reference_code/bedtools regardless
// of the working directory.
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

// upstreamBedtools builds (once) and returns the path to the upstream
// `bedtools` binary in reference_code/bedtools/bin.
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
		t.Fatalf("building upstream bedtools: %v (run `git submodule update --init reference_code/bedtools`)", upstreamBedtoolsErr)
	}
	if upstreamBedtoolsPath == "" {
		t.Fatalf("upstream bedtools binary not found after build")
	}
	return upstreamBedtoolsPath
}

// genomePath returns the absolute path to a vendored testdata/parity genome.
func genomePath(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join("..", "..", "testdata", "parity", name)
}

// runUpstreamRandom runs `bedtools random -g GENOME -l L -n N -seed SEED` and
// returns stdout.
func runUpstreamRandom(t *testing.T, genome string, l, n int, seed int64) []byte {
	t.Helper()
	bin := upstreamBedtools(t)
	cmd := exec.Command(bin, "random",
		"-g", genome,
		"-l", strconv.Itoa(l),
		"-n", strconv.Itoa(n),
		"-seed", strconv.FormatInt(seed, 10),
	)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream random -g %s -l %d -n %d -seed %d: %v\nstderr: %s",
			genome, l, n, seed, err, stderr.String())
	}
	return stdout.Bytes()
}

// TestUpstreamParity_Random_Seeded asserts byte-for-byte parity between this
// port and upstream `bedtools random` across several (genome, L, N, seed)
// tuples. The small genome with a large -l exercises the off-the-end redraw
// path heavily; the hg19 genome exercises projection across many chromosomes.
func TestUpstreamParity_Random_Seeded(t *testing.T) {
	cases := []struct {
		genome string
		l, n   int
		seed   int64
	}{
		{"human.hg19.genome", 100, 50, 42},
		{"human.hg19.genome", 100, 25, 1},
		{"human.hg19.genome", 5000, 30, 123},
		{"human.hg19.genome", 1, 40, 7},
		{"small.genome", 100, 60, 42},
		{"small.genome", 400, 50, 7}, // large L vs 500bp chr3 → many redraws
		{"small.genome", 100, 30, 99},
		{"tiny.genome", 100, 40, 5},
	}
	for _, tc := range cases {
		tc := tc
		t.Run("", func(t *testing.T) {
			gpath := genomePath(t, tc.genome)
			want := runUpstreamRandom(t, gpath, tc.l, tc.n, tc.seed)

			gf, err := os.Open(gpath)
			if err != nil {
				t.Fatalf("open genome: %v", err)
			}
			genome, err := ParseGenome(gf)
			gf.Close()
			if err != nil {
				t.Fatalf("ParseGenome: %v", err)
			}

			var got bytes.Buffer
			if _, err := Generate(genome, &got, Options{
				Length: tc.l, Num: tc.n, Seed: int(tc.seed), HaveSeed: true,
			}); err != nil {
				t.Fatalf("port Generate: %v", err)
			}
			if !bytes.Equal(want, got.Bytes()) {
				t.Fatalf("random mismatch genome=%s l=%d n=%d seed=%d:\nupstream:\n%s\nport:\n%s",
					tc.genome, tc.l, tc.n, tc.seed, want, got.Bytes())
			}
		})
	}
}
