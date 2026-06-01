package main

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

// refDir is the htslib reference checkout that ships the genuine `tabix` and
// `bgzip` binaries used for live-oracle comparisons.
const refDir = "../../../../reference_code/htslib"

// sampleVCF is the body shared by most fixtures: two contigs, six records.
const sampleVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr2>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	.	.	.
chr1	150	.	C	G	.	.	.
chr1	200	.	G	A	.	.	.
chr1	250	.	T	C	.	.	.
chr2	100	.	A	T	.	.	.
chr2	300	.	G	C	.	.	.
`

// refBin returns the path to a reference htslib binary, or "" if it is not
// present/executable so the caller can skip the live-oracle assertion.
func refBin(name string) string {
	p, err := filepath.Abs(filepath.Join(refDir, name))
	if err != nil {
		return ""
	}
	if fi, err := os.Stat(p); err != nil || fi.IsDir() {
		return ""
	}
	return p
}

// writeBGZF compresses content with the reference bgzip and builds a .tbi
// index with the reference tabix, returning the path to the .gz file. It
// skips the test when the reference binaries are unavailable.
func writeBGZF(t *testing.T, dir, name, content string) string {
	t.Helper()
	bgzipBin := refBin("bgzip")
	tabixBin := refBin("tabix")
	if bgzipBin == "" || tabixBin == "" {
		t.Skip("reference htslib bgzip/tabix binary not available")
	}
	plain := filepath.Join(dir, name)
	if err := os.WriteFile(plain, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", plain, err)
	}
	gz := plain + ".gz"
	if out, err := exec.Command(bgzipBin, "-f", plain).CombinedOutput(); err != nil {
		t.Fatalf("bgzip: %v: %s", err, out)
	}
	if out, err := exec.Command(tabixBin, "-p", "vcf", "-f", gz).CombinedOutput(); err != nil {
		t.Fatalf("tabix index: %v: %s", err, out)
	}
	return gz
}

// runArgs invokes the tool's run() with args and returns captured stdout.
func runArgs(t *testing.T, args ...string) (string, int) {
	t.Helper()
	var stdout, stderr bytes.Buffer
	code := run(args, strings.NewReader(""), &stdout, &stderr)
	if code != 0 {
		t.Logf("stderr: %s", stderr.String())
	}
	return stdout.String(), code
}

// referenceTabix runs the genuine tabix binary and returns its stdout.
func referenceTabix(t *testing.T, args ...string) string {
	t.Helper()
	bin := refBin("tabix")
	if bin == "" {
		t.Skip("reference tabix binary not available")
	}
	cmd := exec.Command(bin, args...)
	var out bytes.Buffer
	cmd.Stdout = &out
	if err := cmd.Run(); err != nil {
		t.Fatalf("reference tabix %v: %v", args, err)
	}
	return out.String()
}

func TestTargetsMatchesReference(t *testing.T) {
	dir := t.TempDir()
	gz := writeBGZF(t, dir, "sample.vcf", sampleVCF)

	// Tab-format (1-based) and BED-format target files.
	tab := filepath.Join(dir, "targets.txt")
	if err := os.WriteFile(tab, []byte("chr1\t145\t205\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	bed := filepath.Join(dir, "targets.bed")
	if err := os.WriteFile(bed, []byte("chr1\t99\t150\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	cases := [][]string{
		{"-T", tab, gz},
		{"-T", tab, gz, "chr2"}, // region chr2 + chr1 targets => empty
		{"-T", bed, gz},
		{"-R", tab, gz, "chr2"}, // -R differs: regions union
	}
	for _, args := range cases {
		t.Run(strings.Join(args, "_"), func(t *testing.T) {
			want := referenceTabix(t, args...)
			got, code := runArgs(t, args...)
			if code != 0 {
				t.Fatalf("run exited %d", code)
			}
			if got != want {
				t.Errorf("output mismatch\nargs: %v\n got: %q\nwant: %q", args, got, want)
			}
		})
	}
}

// TestTargetsStrictPostFilter checks the key semantic distinction from -R:
// a positional region combined with a non-overlapping targets file yields no
// records, because --targets is a post-filter rather than a region union.
func TestTargetsStrictPostFilter(t *testing.T) {
	dir := t.TempDir()
	gz := writeBGZF(t, dir, "sample.vcf", sampleVCF)
	tab := filepath.Join(dir, "targets.txt")
	if err := os.WriteFile(tab, []byte("chr1\t145\t205\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	// -T targets(chr1) with region chr2: post-filter removes everything.
	got, code := runArgs(t, "-T", tab, gz, "chr2")
	if code != 0 {
		t.Fatalf("run exited %d", code)
	}
	if got != "" {
		t.Errorf("expected empty output, got %q", got)
	}

	// -R targets(chr1) with region chr2: union, so chr2 records are kept.
	got, code = runArgs(t, "-R", tab, gz, "chr2")
	if code != 0 {
		t.Fatalf("run exited %d", code)
	}
	if !strings.Contains(got, "chr2\t100") || !strings.Contains(got, "chr1\t150") {
		t.Errorf("-R union expected chr1+chr2 records, got %q", got)
	}
}

// Ensure tabix.Config remains importable (compile-time guard for the
// post-filter parser that depends on its exported fields).
var _ = tabix.Config{}
