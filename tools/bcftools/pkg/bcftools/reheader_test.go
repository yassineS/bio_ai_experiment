package bcftools

import (
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const reheaderVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	rs1	A	T	30	PASS	DP=10	GT	0/1	0/0
chr1	200	rs2	C	G	30	PASS	DP=20	GT	1/1	0/1
`

func TestReheaderSamplesPositional(t *testing.T) {
	dir := t.TempDir()
	srename := filepath.Join(dir, "names.txt")
	if err := os.WriteFile(srename, []byte("NEW1\nNEW2\n"), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	n, err := Reheader(strings.NewReader(reheaderVCF), &out, ReheaderOptions{SamplesFile: srename})
	if err != nil {
		t.Fatal(err)
	}
	if n != 2 {
		t.Errorf("record count = %d, want 2", n)
	}
	if !strings.Contains(out.String(), "\tNEW1\tNEW2\n") {
		t.Errorf("samples not renamed:\n%s", out.String())
	}
	if strings.Contains(out.String(), "\tS1\tS2\n") {
		t.Error("old sample names still present")
	}
}

func TestReheaderSamplesMapping(t *testing.T) {
	dir := t.TempDir()
	srename := filepath.Join(dir, "map.txt")
	if err := os.WriteFile(srename, []byte("S2\tRENAMED2\n# comment\nS1\tRENAMED1\n"), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := Reheader(strings.NewReader(reheaderVCF), &out, ReheaderOptions{SamplesFile: srename}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out.String(), "RENAMED1") || !strings.Contains(out.String(), "RENAMED2") {
		t.Errorf("expected renamed names:\n%s", out.String())
	}
	// Order should still be (S1->RENAMED1, S2->RENAMED2) — i.e. the original
	// header sample order, not the file order.
	idx1 := strings.Index(out.String(), "RENAMED1")
	idx2 := strings.Index(out.String(), "RENAMED2")
	if idx1 > idx2 {
		t.Errorf("expected RENAMED1 before RENAMED2 (header order):\n%s", out.String())
	}
}

func TestReheaderHeaderFile(t *testing.T) {
	dir := t.TempDir()
	hpath := filepath.Join(dir, "newhdr.vcf")
	const newHdr = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=99999>
##source=ReheaderTest
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
`
	if err := os.WriteFile(hpath, []byte(newHdr), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := Reheader(strings.NewReader(reheaderVCF), &out, ReheaderOptions{HeaderFile: hpath}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out.String(), "##source=ReheaderTest") {
		t.Errorf("new header not used:\n%s", out.String())
	}
	if !strings.Contains(out.String(), "##contig=<ID=chr1,length=99999>") {
		t.Errorf("new contig length missing:\n%s", out.String())
	}
	// Records should still be present with the original samples (matching
	// the new header's sample list).
	if !strings.Contains(out.String(), "rs1") || !strings.Contains(out.String(), "rs2") {
		t.Errorf("records dropped:\n%s", out.String())
	}
}

func TestReheaderFAIContigs(t *testing.T) {
	dir := t.TempDir()
	fai := filepath.Join(dir, "ref.fa.fai")
	// Two-column FAI: NAME LENGTH OFFSET LINEBASES LINEWIDTH; we only care
	// about the first two for reheader.
	if err := os.WriteFile(fai, []byte("chr1\t250000000\t6\t60\t61\nchrX\t155270560\t12\t60\t61\n"), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := Reheader(strings.NewReader(reheaderVCF), &out, ReheaderOptions{FaiFile: fai}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out.String(), "##contig=<ID=chr1,length=250000000>") {
		t.Errorf("FAI chr1 length not applied:\n%s", out.String())
	}
	if !strings.Contains(out.String(), "##contig=<ID=chrX,length=155270560>") {
		t.Errorf("FAI chrX missing:\n%s", out.String())
	}
}

func TestReheaderFile(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(in, []byte(reheaderVCF), 0644); err != nil {
		t.Fatal(err)
	}
	srename := filepath.Join(dir, "names.txt")
	if err := os.WriteFile(srename, []byte("X1\nX2\n"), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	n, err := ReheaderFile(in, &out, ReheaderOptions{SamplesFile: srename})
	if err != nil {
		t.Fatal(err)
	}
	if n != 2 {
		t.Errorf("ReheaderFile got %d records, want 2", n)
	}
	if !strings.Contains(out.String(), "X1") || !strings.Contains(out.String(), "X2") {
		t.Errorf("ReheaderFile names not applied:\n%s", out.String())
	}
}

// TestReheaderMirrorsInputCompression verifies that ReheaderFile re-emits a
// compressed (.gz) input as compressed output and a plain input as plain
// output, mirroring upstream `bcftools reheader` (which has no -O flag), and
// that an explicit OutputFormat overrides that auto-detection.
func TestReheaderMirrorsInputCompression(t *testing.T) {
	dir := t.TempDir()
	srename := filepath.Join(dir, "names.txt")
	if err := os.WriteFile(srename, []byte("X1\nX2\n"), 0644); err != nil {
		t.Fatal(err)
	}
	isGzip := func(b []byte) bool { return len(b) >= 2 && b[0] == 0x1f && b[1] == 0x8b }

	// Plain input -> plain output.
	plain := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(plain, []byte(reheaderVCF), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := ReheaderFile(plain, &out, ReheaderOptions{SamplesFile: srename}); err != nil {
		t.Fatal(err)
	}
	if isGzip(out.Bytes()) {
		t.Error("plain input produced compressed output; want plain")
	}
	if !strings.Contains(out.String(), "X1") {
		t.Errorf("rename not applied:\n%s", out.String())
	}

	// Gzip input -> compressed output (magic 0x1f 0x8b).
	gzPath := filepath.Join(dir, "in.vcf.gz")
	var gzBuf bytes.Buffer
	gw := gzip.NewWriter(&gzBuf)
	if _, err := gw.Write([]byte(reheaderVCF)); err != nil {
		t.Fatal(err)
	}
	if err := gw.Close(); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(gzPath, gzBuf.Bytes(), 0644); err != nil {
		t.Fatal(err)
	}
	out.Reset()
	if _, err := ReheaderFile(gzPath, &out, ReheaderOptions{SamplesFile: srename}); err != nil {
		t.Fatal(err)
	}
	if !isGzip(out.Bytes()) {
		t.Error("gzip input produced plain output; want compressed")
	}

	// Explicit OutputFormat wins over the input-compression mirror.
	out.Reset()
	if _, err := ReheaderFile(gzPath, &out, ReheaderOptions{
		SamplesFile: srename, OutputFormat: OutputVCF, OutputFormatExplicit: true,
	}); err != nil {
		t.Fatal(err)
	}
	if isGzip(out.Bytes()) {
		t.Error("explicit -O v on gzip input produced compressed output; want plain")
	}
}

func TestReheaderMissingFile(t *testing.T) {
	if _, err := ReheaderFile("/no/such/file", &bytes.Buffer{}, ReheaderOptions{}); err == nil {
		t.Error("expected error for missing input")
	}
}

func TestReheaderBadFAI(t *testing.T) {
	dir := t.TempDir()
	fai := filepath.Join(dir, "bad.fai")
	if err := os.WriteFile(fai, []byte("chr1\tnotanumber\n"), 0644); err != nil {
		t.Fatal(err)
	}
	if _, err := Reheader(strings.NewReader(reheaderVCF), &bytes.Buffer{}, ReheaderOptions{FaiFile: fai}); err == nil {
		t.Error("expected error parsing bad FAI")
	}
}
