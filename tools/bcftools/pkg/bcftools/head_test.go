package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const headTestVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	rs1	A	T	30	PASS	DP=10	GT	0/1	0/0
chr1	200	rs2	C	G	30	PASS	DP=20	GT	1/1	0/1
`

func TestHeadFullHeader(t *testing.T) {
	var out bytes.Buffer
	if err := Head(strings.NewReader(headTestVCF), &out, HeadOptions{}); err != nil {
		t.Fatal(err)
	}
	got := out.String()
	for _, want := range []string{
		"##fileformat=VCFv4.2",
		"##contig=<ID=chr1",
		"##contig=<ID=chr2",
		"##INFO=<ID=DP",
		"##FORMAT=<ID=GT",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
	} {
		if !strings.Contains(got, want) {
			t.Errorf("missing %q in head output:\n%s", want, got)
		}
	}
	// No data rows.
	if strings.Contains(got, "rs1") || strings.Contains(got, "rs2") {
		t.Errorf("unexpected data row in head output:\n%s", got)
	}
}

func TestHeadNumLinesCap(t *testing.T) {
	var out bytes.Buffer
	if err := Head(strings.NewReader(headTestVCF), &out, HeadOptions{NumLines: 2}); err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(strings.TrimRight(out.String(), "\n"), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected 2 lines, got %d:\n%s", len(lines), out.String())
	}
	if lines[0] != "##fileformat=VCFv4.2" {
		t.Errorf("first line = %q", lines[0])
	}
}

func TestHeadSamplesOnly(t *testing.T) {
	var out bytes.Buffer
	if err := Head(strings.NewReader(headTestVCF), &out, HeadOptions{SamplesOnly: true}); err != nil {
		t.Fatal(err)
	}
	got := strings.TrimRight(out.String(), "\n")
	if got != "S1\nS2" {
		t.Errorf("samples-only got %q want \"S1\\nS2\"", got)
	}
}

func TestHeadFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(path, []byte(headTestVCF), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if err := HeadFile(path, &out, HeadOptions{NumLines: 1}); err != nil {
		t.Fatal(err)
	}
	if !strings.HasPrefix(out.String(), "##fileformat=VCFv4.2") {
		t.Errorf("HeadFile content unexpected:\n%s", out.String())
	}
}

func TestHeadEmptyHeaderError(t *testing.T) {
	var out bytes.Buffer
	if err := Head(strings.NewReader("not a vcf"), &out, HeadOptions{}); err == nil {
		t.Error("expected error reading malformed input")
	}
}

func TestHeadNoSamplesNoFormat(t *testing.T) {
	const noSamples = `##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1	.	A	T	.	PASS	DP=1
`
	var out bytes.Buffer
	if err := Head(strings.NewReader(noSamples), &out, HeadOptions{}); err != nil {
		t.Fatal(err)
	}
	if strings.Contains(out.String(), "FORMAT") {
		t.Errorf("FORMAT column should not appear when no samples:\n%s", out.String())
	}
}
