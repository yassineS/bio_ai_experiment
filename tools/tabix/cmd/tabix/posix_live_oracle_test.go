package main

// Live-binary oracle test for POSIX getopt-style short-flag bundling in
// tabix, now that the CLI is routed through cliflag.Parse.
//
// The cross-binary assertion runs the genuine upstream `tabix` and our port
// on the same bgzipped VCF, building the index with a value-concatenated
// `-pvcf` (== `-p vcf`) cluster and querying a region; the emitted records
// must match line-for-line. The intra-binary assertion proves the bundled /
// value-concatenated spellings behave identically to the canonical
// spelled-out forms within our own binary.
//
// Per the project's testing rules the helpers t.Fatalf rather than t.Skip
// when the upstream binary cannot be built.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"testing"
)

var (
	tabixOurOnce sync.Once
	tabixOurPath string
	tabixOurErr  error
)

func buildOurTabixBinary(t *testing.T) string {
	t.Helper()
	tabixOurOnce.Do(func() {
		dir, err := os.MkdirTemp("", "our-tabix-")
		if err != nil {
			tabixOurErr = err
			return
		}
		bin := filepath.Join(dir, "tabix")
		cmd := exec.Command("go", "build", "-o", bin, ".")
		if out, err := cmd.CombinedOutput(); err != nil {
			tabixOurErr = tabixBuildErr{err: err, out: out}
			return
		}
		tabixOurPath = bin
	})
	if tabixOurErr != nil {
		t.Fatalf("build our tabix: %v", tabixOurErr)
	}
	return tabixOurPath
}

type tabixBuildErr struct {
	err error
	out []byte
}

func (e tabixBuildErr) Error() string { return e.err.Error() + ": " + string(e.out) }

func tabixRepoRoot(t *testing.T) string {
	t.Helper()
	dir, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatal("could not locate repo root (go.mod)")
		}
		dir = parent
	}
}

func upstreamTabixBinary(t *testing.T) string {
	t.Helper()
	bin := filepath.Join(tabixRepoRoot(t), "reference_code", "htslib", "tabix")
	if _, err := os.Stat(bin); err != nil {
		t.Skipf("upstream tabix not found at %s (build reference_code/htslib first): %v", bin, err)
	}
	return bin
}

func upstreamBgzip(t *testing.T) string {
	t.Helper()
	bin := filepath.Join(tabixRepoRoot(t), "reference_code", "htslib", "bgzip")
	if _, err := os.Stat(bin); err != nil {
		t.Skipf("upstream bgzip not found at %s: %v", bin, err)
	}
	return bin
}

func runTabix(t *testing.T, bin, dir string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	cmd.Dir = dir
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	if err := cmd.Run(); err != nil {
		t.Fatalf("%s %v failed: %v\nstderr: %s", bin, args, err, errb.String())
	}
	return out.Bytes()
}

const tabixFixtureVCF = "##fileformat=VCFv4.2\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
	"chr1\t100\t.\tA\tT\t.\t.\t.\n" +
	"chr1\t200\t.\tC\tG\t.\t.\t.\n" +
	"chr1\t300\t.\tG\tA\t.\t.\t.\n" +
	"chr2\t150\t.\tG\tA\t.\t.\t.\n"

// makeBgzippedVCF writes the fixture VCF and bgzips it (via upstream bgzip)
// into dir, returning the .vcf.gz path. It does NOT build the .tbi index;
// each test builds it with the spelling under test.
func makeBgzippedVCF(t *testing.T, dir, bgzipBin string) string {
	t.Helper()
	vcf := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(vcf, []byte(tabixFixtureVCF), 0o644); err != nil {
		t.Fatal(err)
	}
	gz := vcf + ".gz"
	out, err := exec.Command(bgzipBin, "-c", vcf).Output()
	if err != nil {
		t.Fatalf("bgzip fixture: %v", err)
	}
	if err := os.WriteFile(gz, out, 0o644); err != nil {
		t.Fatal(err)
	}
	return gz
}

// TestLiveTabixPosixBundledBuildAndQuery is the upstream-parity gate: build
// the index with the value-concatenated `-pvcf` cluster and query a region;
// our port's records must match upstream's line-for-line.
func TestLiveTabixPosixBundledBuildAndQuery(t *testing.T) {
	ours := buildOurTabixBinary(t)
	up := upstreamTabixBinary(t)
	bg := upstreamBgzip(t)

	// Separate fixtures so the two .tbi files don't collide.
	upDir := t.TempDir()
	ourDir := t.TempDir()
	upGz := makeBgzippedVCF(t, upDir, bg)
	ourGz := makeBgzippedVCF(t, ourDir, bg)

	// Build with the bundled `-pvcf` (== `-p vcf`) form, exactly as upstream
	// getopt accepts it — the form that must now parse in our port too.
	runTabix(t, up, upDir, "-pvcf", upGz)
	runTabix(t, ours, ourDir, "-pvcf", ourGz)

	upRecords := runTabix(t, up, upDir, upGz, "chr1:100-250")
	ourRecords := runTabix(t, ours, ourDir, ourGz, "chr1:100-250")
	if !bytes.Equal(upRecords, ourRecords) {
		t.Fatalf("region query mismatch:\nupstream:\n%s\nours:\n%s", upRecords, ourRecords)
	}
}

// TestTabixPosixBundlingEquivalentToCanonical proves, within our binary, that
// the value-concatenated build flag (`-pvcf`) and a bundled query header flag
// produce identical results to the canonical spelled-out forms.
func TestTabixPosixBundlingEquivalentToCanonical(t *testing.T) {
	ours := buildOurTabixBinary(t)
	bg := upstreamBgzip(t)

	// Build-flag equivalence: -pvcf == -p vcf (verified by querying the
	// resulting index).
	bundledDir := t.TempDir()
	canonDir := t.TempDir()
	bundledGz := makeBgzippedVCF(t, bundledDir, bg)
	canonGz := makeBgzippedVCF(t, canonDir, bg)
	runTabix(t, ours, bundledDir, "-pvcf", bundledGz)
	runTabix(t, ours, canonDir, "-p", "vcf", canonGz)

	bundled := runTabix(t, ours, bundledDir, bundledGz, "chr1")
	canonical := runTabix(t, ours, canonDir, canonGz, "chr1")
	if !bytes.Equal(bundled, canonical) {
		t.Fatalf("build -pvcf query differs from -p vcf:\n-pvcf:\n%s\n-p vcf:\n%s", bundled, canonical)
	}

	// Query-flag equivalence: -h (print-header) on a query stream.
	withHdrBundled := runTabix(t, ours, bundledDir, bundledGz, "-h", "chr1:100-250")
	withHdrCanon := runTabix(t, ours, canonDir, canonGz, "--print-header", "chr1:100-250")
	if !bytes.Equal(withHdrBundled, withHdrCanon) {
		t.Fatalf("-h query differs from --print-header:\n-h:\n%s\n--print-header:\n%s", withHdrBundled, withHdrCanon)
	}
}
