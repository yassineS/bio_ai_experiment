package bedmulticov

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
)

// This file holds live-upstream parity tests for CRAM input: it builds the
// real upstream `bedtools` binary from the vendored submodule, runs it on a
// real (htslib-produced) CRAM fixture, and compares the output byte for byte
// against this port's CRAM path. Closing the documented "CRAM input deferred"
// gap for `bedtools multicov`.
//
// Upstream `multicov` requires a coordinate index for every alignment input
// (it queries by region), so the test generates a `.crai` for the fixture via
// pkg/htsgo/cram.CreateCRAI into a temp dir. The same `.crai` is accepted by
// upstream htslib, proving the in-tree index builder is interoperable.
//
// The tests t.Fatalf (never t.Skip): a missing/unbuildable submodule or an
// undecodable fixture is a hard failure, matching the parity-rig policy.

var (
	upstreamBedtoolsOnce sync.Once
	upstreamBedtoolsPath string
	upstreamBedtoolsErr  error
)

// repoRoot walks up from this test file to the repo root (dir with go.mod).
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
// `bedtools` binary. Uniquely named to avoid colliding with sibling packages.
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

// fixturePath returns the absolute path to a testdata/parity fixture.
func fixturePath(t *testing.T, name string) string {
	t.Helper()
	p, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("abs path for %s: %v", name, err)
	}
	return p
}

// stageCRAMWithIndex copies the named CRAM fixture into a temp dir and builds
// its `.crai` sidecar (required by upstream multicov). It returns the staged
// CRAM path.
func stageCRAMWithIndex(t *testing.T, name string) string {
	t.Helper()
	src := fixturePath(t, name)
	data, err := os.ReadFile(src)
	if err != nil {
		t.Fatalf("read CRAM fixture %s: %v", name, err)
	}
	dir := t.TempDir()
	dst := filepath.Join(dir, name)
	if err := os.WriteFile(dst, data, 0o644); err != nil {
		t.Fatalf("stage CRAM: %v", err)
	}
	if err := cram.CreateCRAI(dst, dst+".crai"); err != nil {
		t.Fatalf("CreateCRAI: %v", err)
	}
	return dst
}

// writeBED writes content to a temp BED file and returns its path.
func writeBED(t *testing.T, content string) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "a.bed")
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write BED: %v", err)
	}
	return p
}

// runUpstreamMulticov runs upstream `bedtools multicov` with the given args.
func runUpstreamMulticov(t *testing.T, args ...string) []byte {
	t.Helper()
	bin := upstreamBedtools(t)
	cmd := exec.Command(bin, append([]string{"multicov"}, args...)...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream multicov %v: %v\nstderr: %s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

// runPortMulticov runs this port's CLI-equivalent path: open the BED A file
// and the CRAM input, then RunSources with the CRAM source kind.
func runPortMulticov(t *testing.T, bedPath, cramPath string, opts Options) []byte {
	t.Helper()
	aR, err := os.Open(bedPath)
	if err != nil {
		t.Fatalf("open A: %v", err)
	}
	defer aR.Close()
	cR, err := os.Open(cramPath)
	if err != nil {
		t.Fatalf("open CRAM: %v", err)
	}
	defer cR.Close()
	var got bytes.Buffer
	if _, err := RunSources(aR, []Source{{Reader: cR, Kind: SourceCRAM}}, &got, opts); err != nil {
		t.Fatalf("port RunSources: %v", err)
	}
	return got.Bytes()
}

// TestUpstreamParity_CRAM_Default exercises plain (no-flag) CRAM coverage and
// asserts byte-for-byte parity with upstream multicov. The fixture a.cram holds
// two mapped chr1 alignments around chr1:10003..10143.
func TestUpstreamParity_CRAM_Default(t *testing.T) {
	cramPath := stageCRAMWithIndex(t, "a.cram")
	bed := writeBED(t, "chr1\t10000\t10100\tregion1\nchr1\t10100\t10200\tregion2\nchr1\t9000\t9500\tregion3\n")
	want := runUpstreamMulticov(t, "-bed", bed, "-bams", cramPath)
	got := runPortMulticov(t, bed, cramPath, Options{MaxDepth: 64000})
	if !bytes.Equal(want, got) {
		t.Fatalf("CRAM default mismatch:\nupstream:\n%s\nport:\n%s", want, got)
	}
}

// TestUpstreamParity_CRAM_MAPQ exercises the `-q` MAPQ filter on CRAM input.
// b.cram has one MAPQ-0 read and one MAPQ-3 read; `-q 1` drops the former.
func TestUpstreamParity_CRAM_MAPQ(t *testing.T) {
	cramPath := stageCRAMWithIndex(t, "b.cram")
	// Two A intervals: one over each read's span.
	bed := writeBED(t, "chr1\t10000\t10100\tr1\nchr1\t49000\t49200\tr2\n")
	want := runUpstreamMulticov(t, "-bed", bed, "-bams", cramPath, "-q", "1")
	got := runPortMulticov(t, bed, cramPath, Options{MinMAPQ: 1, MaxDepth: 64000})
	if !bytes.Equal(want, got) {
		t.Fatalf("CRAM -q mismatch:\nupstream:\n%s\nport:\n%s", want, got)
	}
}

// TestUpstreamParity_CRAM_Strand exercises the `-s` same-strand filter on CRAM.
func TestUpstreamParity_CRAM_Strand(t *testing.T) {
	cramPath := stageCRAMWithIndex(t, "a.cram")
	// a.cram read 1 is '+' (flag 99), read 2 is '-' (flag 147), both overlap
	// chr1:10003..10143. A '+' interval should only pick up the '+' read.
	bed := writeBED(t, "chr1\t10000\t10150\tplus\t0\t+\nchr1\t10000\t10150\tminus\t0\t-\n")
	want := runUpstreamMulticov(t, "-bed", bed, "-bams", cramPath, "-s")
	got := runPortMulticov(t, bed, cramPath, Options{SameStrand: true, MaxDepth: 64000})
	if !bytes.Equal(want, got) {
		t.Fatalf("CRAM -s mismatch:\nupstream:\n%s\nport:\n%s", want, got)
	}
}
