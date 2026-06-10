package bedmultiinter

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

// This file holds live-upstream parity tests: they build the real upstream
// `bedtools` binary from the vendored submodule and compare its output, byte
// for byte, against this port's. They exercise the VCF/GFF autodetection
// added to `multiinter`, which upstream gets from its shared BedFile parser.
//
// The tests t.Fatalf (never t.Skip) so a missing or unbuildable submodule is
// a hard failure, matching the project's established parity-rig policy.

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
// `bedtools` binary in reference_code/bedtools/bin. It is uniquely named to
// avoid colliding with builders in sibling packages.
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

// fixturePath returns the absolute path to a testdata/parity fixture.
func fixturePath(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join("..", "..", "testdata", "parity", name)
}

// runUpstreamMultiinter runs `bedtools multiinter -i <files...>` and returns
// its stdout.
func runUpstreamMultiinter(t *testing.T, files ...string) []byte {
	t.Helper()
	bin := upstreamBedtools(t)
	args := append([]string{"multiinter", "-i"}, files...)
	cmd := exec.Command(bin, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream multiinter %v: %v\nstderr: %s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

// runPortMultiinter runs this port's Run over the same files using default
// options (no -names, no -header), matching the upstream invocation above.
func runPortMultiinter(t *testing.T, files ...string) []byte {
	t.Helper()
	bRs := make([]io.Reader, len(files))
	for i, f := range files {
		fh, err := os.Open(f)
		if err != nil {
			t.Fatalf("open fixture %s: %v", f, err)
		}
		defer fh.Close()
		bRs[i] = fh
	}
	var got bytes.Buffer
	if _, err := Run(bRs, &got, Options{Filenames: files}); err != nil {
		t.Fatalf("port Run: %v", err)
	}
	return got.Bytes()
}

// TestUpstreamParity_VCF_GFF_Mix runs a 3-way VCF + GFF + BED intersection and
// asserts byte-for-byte parity with upstream `bedtools multiinter`. This is the
// closing case for the documented "VCF/GFF input not implemented" gap.
func TestUpstreamParity_VCF_GFF_Mix(t *testing.T) {
	vcf := fixturePath(t, "multi.vcf")
	gff := fixturePath(t, "multi.gff")
	bed := fixturePath(t, "multi.bed")
	want := runUpstreamMultiinter(t, vcf, gff, bed)
	got := runPortMultiinter(t, vcf, gff, bed)
	if !bytes.Equal(want, got) {
		t.Fatalf("multiinter VCF/GFF/BED mismatch:\nupstream:\n%s\nport:\n%s", want, got)
	}
}

// TestUpstreamParity_VCF_GFF_SingleRecord uses the upstream issue311 fixtures
// (one VCF row and one GFF row that both cover chr1:31-32) to lock the
// coordinate-conversion semantics (VCF POS-1 + len(REF); GFF start-1, end).
func TestUpstreamParity_VCF_GFF_SingleRecord(t *testing.T) {
	vcf := fixturePath(t, "issue311.vcf")
	gff := fixturePath(t, "issue311.gff")
	want := runUpstreamMultiinter(t, vcf, gff)
	got := runPortMultiinter(t, vcf, gff)
	if !bytes.Equal(want, got) {
		t.Fatalf("issue311 VCF/GFF mismatch:\nupstream:\n%s\nport:\n%s", want, got)
	}
}
