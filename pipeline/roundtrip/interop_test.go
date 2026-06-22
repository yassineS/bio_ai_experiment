package roundtrip

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pipeline/fixtures"
	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// --- pure unit tests (no binaries; always run, including under -short) ------

// TestSecondContigVCF picks the second distinct contig from a VCF body so a
// region query crosses a contig boundary on multi-contig fixtures.
func TestSecondContigVCF(t *testing.T) {
	dir := t.TempDir()
	vcf := filepath.Join(dir, "v.vcf")
	body := "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t10\t.\tA\tC\t.\t.\t.\n" +
		"chr1\t20\t.\tA\tC\t.\t.\t.\n" +
		"chr2\t30\t.\tA\tC\t.\t.\t.\n" +
		"chr3\t40\t.\tA\tC\t.\t.\t.\n"
	if err := os.WriteFile(vcf, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	got, err := secondContigVCF(vcf)
	if err != nil {
		t.Fatalf("secondContigVCF: %v", err)
	}
	if got != "chr2" {
		t.Fatalf("secondContigVCF = %q, want chr2", got)
	}
}

// TestSecondContigVCFSingle falls back to the first contig when only one exists.
func TestSecondContigVCFSingle(t *testing.T) {
	dir := t.TempDir()
	vcf := filepath.Join(dir, "v.vcf")
	body := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chrA\t1\t.\tA\tC\t.\t.\t.\n" +
		"chrA\t2\t.\tA\tC\t.\t.\t.\n"
	if err := os.WriteFile(vcf, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	got, err := secondContigVCF(vcf)
	if err != nil {
		t.Fatalf("secondContigVCF: %v", err)
	}
	if got != "chrA" {
		t.Fatalf("secondContigVCF = %q, want chrA fallback", got)
	}
}

// TestCopyFileRoundTrip checks the scratch-copy helper used by index interop.
func TestCopyFileRoundTrip(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(dir, "src")
	dst := filepath.Join(dir, "dst")
	want := []byte("the quick brown fox\x00\x01\x02")
	if err := os.WriteFile(src, want, 0o644); err != nil {
		t.Fatal(err)
	}
	if err := copyFile(dst, src); err != nil {
		t.Fatalf("copyFile: %v", err)
	}
	got, err := os.ReadFile(dst)
	if err != nil {
		t.Fatal(err)
	}
	if string(got) != string(want) {
		t.Fatalf("copyFile content mismatch")
	}
}

// --- live interop suite (needs upstream binaries; skipped under -short) -----

// requireInterop gates the live interop checks: it skips cleanly under -short
// and when the upstream samtools binary (the common dependency) is unavailable,
// keeping `go test -short ./...` and the macOS CI job green without building any
// upstream binary.
func requireInterop(t *testing.T) {
	t.Helper()
	if testing.Short() {
		t.Skip("interop checks shell out to upstream binaries; skipped under -short")
	}
	if _, err := upstream.Binary("samtools"); err != nil {
		t.Skipf("upstream samtools unavailable: %v", err)
	}
}

// interopManifest generates (or reuses) the multi-contig smoke fixture set.
// Generation itself uses upstream binaries, so callers must gate with
// requireInterop first. The smoke tier has ≥2 contigs, so cross-contig bins,
// RNEXT, and coordinate-sort ordering are exercised.
func interopManifest(t *testing.T) (*fixtures.Manifest, string) {
	t.Helper()
	root, err := upstream.RepoRoot()
	if err != nil {
		t.Fatal(err)
	}
	man, err := fixtures.Generate(fixtures.Options{Scale: fixtures.Smoke})
	if err != nil {
		t.Skipf("fixtures unavailable: %v", err)
	}
	if man.Params.NumContigs < 2 {
		t.Fatalf("fixture is not multi-contig (NumContigs=%d)", man.Params.NumContigs)
	}
	cache := filepath.Join(root, "pipeline", ".fixtures", "bin")
	return man, cache
}

// newEnv builds the interop env over a fresh scratch dir (cleaned up by the
// test) so the live checks have somewhere to write.
func newEnv(t *testing.T, man *fixtures.Manifest, cache string) *env {
	t.Helper()
	tmp := t.TempDir()
	return &env{man: man, cache: cache, tmp: tmp}
}

// expectInterop runs one interop check and requires it to PASS (a SKIP here
// would mean the upstream binary that requireInterop did not gate is missing;
// since samtools is present, the htslib/bcftools binaries built from the same
// submodules are expected too — but we still tolerate a SKIP rather than fail).
func expectInterop(t *testing.T, r Result) {
	t.Helper()
	switch r.Status {
	case Pass:
		t.Logf("%s (%s): PASS — %s", r.Name, r.Format, r.Detail)
	case Skip:
		t.Skipf("%s (%s): SKIP — %s", r.Name, r.Format, r.Detail)
	default:
		t.Fatalf("%s (%s): %s — %s", r.Name, r.Format, r.Status, r.Detail)
	}
}

func TestBGZFInterop(t *testing.T) {
	requireInterop(t)
	m, c := interopManifest(t)
	expectInterop(t, newEnv(t, m, c).bgzfInterop())
}
func TestBAMInterop(t *testing.T) {
	requireInterop(t)
	m, c := interopManifest(t)
	expectInterop(t, newEnv(t, m, c).bamInterop())
}
func TestCRAMInterop(t *testing.T) {
	requireInterop(t)
	m, c := interopManifest(t)
	expectInterop(t, newEnv(t, m, c).cramInterop())
}
func TestVCFGzInterop(t *testing.T) {
	requireInterop(t)
	m, c := interopManifest(t)
	expectInterop(t, newEnv(t, m, c).vcfGzInterop())
}
func TestBCFInterop(t *testing.T) {
	requireInterop(t)
	m, c := interopManifest(t)
	expectInterop(t, newEnv(t, m, c).bcfInterop())
}
func TestFASTQInterop(t *testing.T) {
	requireInterop(t)
	m, c := interopManifest(t)
	expectInterop(t, newEnv(t, m, c).fastqInterop())
}
func TestBAIInterop(t *testing.T) {
	requireInterop(t)
	m, c := interopManifest(t)
	expectInterop(t, newEnv(t, m, c).baiInterop())
}
func TestCSIInterop(t *testing.T) {
	requireInterop(t)
	m, c := interopManifest(t)
	expectInterop(t, newEnv(t, m, c).csiInterop())
}
func TestTBIInterop(t *testing.T) {
	requireInterop(t)
	m, c := interopManifest(t)
	expectInterop(t, newEnv(t, m, c).tbiInterop())
}

// TestInteropSuite runs every interop check through the same env and asserts
// none FAILs — the end-to-end gate that pipeline/cmd/full-validation relies on.
func TestInteropSuite(t *testing.T) {
	requireInterop(t)
	man, cache := interopManifest(t)
	e := newEnv(t, man, cache)
	results := e.interopChecks()
	if len(results) != 9 {
		t.Fatalf("expected 9 interop checks, got %d", len(results))
	}
	var passed int
	for _, r := range results {
		switch r.Status {
		case Fail:
			t.Errorf("%s (%s): FAIL — %s", r.Name, r.Format, r.Detail)
		case Pass:
			passed++
		case Skip:
			t.Logf("%s (%s): SKIP — %s", r.Name, r.Format, r.Detail)
		}
	}
	if passed == 0 {
		t.Skip("no interop check ran (all upstream binaries unavailable)")
	}
}
