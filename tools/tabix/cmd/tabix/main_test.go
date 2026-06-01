package main

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
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

// TestReheaderBodyPreserved validates the reheader output: the new header
// replaces the old, the body is byte-identical to the original, and the
// result is valid BGZF that the reference tabix can re-index and query.
func TestReheaderBodyPreserved(t *testing.T) {
	dir := t.TempDir()
	gz := writeBGZF(t, dir, "sample.vcf", sampleVCF)

	newHdr := filepath.Join(dir, "newhdr.txt")
	const hdr = `##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr2>
##source=reheadered
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
`
	if err := os.WriteFile(newHdr, []byte(hdr), 0o644); err != nil {
		t.Fatal(err)
	}

	var stdout, stderr bytes.Buffer
	code := run([]string{"-r", newHdr, gz}, strings.NewReader(""), &stdout, &stderr)
	if code != 0 {
		t.Fatalf("reheader exited %d: %s", code, stderr.String())
	}

	got := decompress(t, stdout.Bytes())
	wantHdr, wantBody := splitVCF(hdr + bodyOf(sampleVCF))
	gotHdr, gotBody := splitVCF(got)
	if gotHdr != wantHdr {
		t.Errorf("header mismatch\n got: %q\nwant: %q", gotHdr, wantHdr)
	}
	if gotBody != wantBody {
		t.Errorf("body mismatch\n got: %q\nwant: %q", gotBody, wantBody)
	}

	// The output must be a valid, re-indexable BGZF file.
	out := filepath.Join(dir, "reheadered.vcf.gz")
	if err := os.WriteFile(out, stdout.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
	if bin := refBin("tabix"); bin != "" {
		if o, err := exec.Command(bin, "-p", "vcf", "-f", out).CombinedOutput(); err != nil {
			t.Fatalf("reference re-index of reheadered file failed: %v: %s", err, o)
		}
	}
}

// TestReheaderTrailingNewline ensures a header file lacking a final newline is
// normalised so the first body record stays on its own line.
func TestReheaderTrailingNewline(t *testing.T) {
	dir := t.TempDir()
	gz := writeBGZF(t, dir, "sample.vcf", sampleVCF)
	newHdr := filepath.Join(dir, "nonl.txt")
	if err := os.WriteFile(newHdr, []byte("##a\n##b\n#CHROM\tPOS"), 0o644); err != nil {
		t.Fatal(err)
	}
	var stdout, stderr bytes.Buffer
	if code := run([]string{"-r", newHdr, gz}, strings.NewReader(""), &stdout, &stderr); code != 0 {
		t.Fatalf("reheader exited %d: %s", code, stderr.String())
	}
	got := decompress(t, stdout.Bytes())
	lines := strings.Split(strings.TrimRight(got, "\n"), "\n")
	if lines[2] != "#CHROM\tPOS" {
		t.Errorf("expected header line preserved, got %q", lines[2])
	}
	if lines[3] != "chr1\t100\t.\tA\tT\t.\t.\t." {
		t.Errorf("expected first record on its own line, got %q", lines[3])
	}
}

// TestReheaderMultiBlock reheaders a file large enough to span many BGZF
// blocks and confirms the body is preserved byte-for-byte.
func TestReheaderMultiBlock(t *testing.T) {
	dir := t.TempDir()
	var b strings.Builder
	b.WriteString("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for i := 1; i <= 20000; i++ {
		b.WriteString("chr1\t")
		b.WriteString(itoa(i))
		b.WriteString("\t.\tA\tT\t.\t.\t.\n")
	}
	gz := writeBGZF(t, dir, "big.vcf", b.String())

	newHdr := filepath.Join(dir, "bighdr.txt")
	const hdr = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##new=yes\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	if err := os.WriteFile(newHdr, []byte(hdr), 0o644); err != nil {
		t.Fatal(err)
	}
	var stdout, stderr bytes.Buffer
	if code := run([]string{"-r", newHdr, gz}, strings.NewReader(""), &stdout, &stderr); code != 0 {
		t.Fatalf("reheader exited %d: %s", code, stderr.String())
	}
	got := decompress(t, stdout.Bytes())
	_, gotBody := splitVCF(got)
	_, wantBody := splitVCF(b.String())
	if gotBody != wantBody {
		t.Errorf("multi-block body mismatch (len got=%d want=%d)", len(gotBody), len(wantBody))
	}
}

// --- helpers ---

func decompress(t *testing.T, gz []byte) string {
	t.Helper()
	r, err := bgzip.NewReader(bytes.NewReader(gz))
	if err != nil {
		t.Fatalf("bgzf reader: %v", err)
	}
	defer r.Close()
	var out bytes.Buffer
	if _, err := out.ReadFrom(r); err != nil {
		t.Fatalf("bgzf read: %v", err)
	}
	return out.String()
}

func splitVCF(s string) (header, body string) {
	var h, b strings.Builder
	for _, line := range strings.SplitAfter(s, "\n") {
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, "#") {
			h.WriteString(line)
		} else {
			b.WriteString(line)
		}
	}
	return h.String(), b.String()
}

func bodyOf(s string) string {
	_, b := splitVCF(s)
	return b
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	return string(buf[i:])
}

// Ensure tabix.Config remains importable (compile-time guard for the
// post-filter parser that depends on its exported fields).
var _ = tabix.Config{}
