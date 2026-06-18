package htsfile

import (
	"bytes"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const viewVCFText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	rs1	A	T	30	PASS	DP=10	GT	0/1
chr1	200	rs2	C	G	.	PASS	DP=20	GT	1/1
`

// TestViewVCF round-trips a plain VCF through View and checks the header and
// records survive unchanged (the htslib -c canonical text form).
func TestViewVCF(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(p, []byte(viewVCFText), 0644); err != nil {
		t.Fatal(err)
	}
	f, err := Identify(p)
	if err != nil {
		t.Fatal(err)
	}
	if f.Payload != PayloadVCF {
		t.Fatalf("identified as %v, want VCF", f.Payload)
	}
	var out bytes.Buffer
	if err := View(p, f, &out); err != nil {
		t.Fatal(err)
	}
	got := out.String()
	for _, want := range []string{"##fileformat=VCFv4.2", "#CHROM\tPOS", "chr1\t100\trs1\tA\tT", "chr1\t200\trs2\tC\tG"} {
		if !strings.Contains(got, want) {
			t.Errorf("view output missing %q:\n%s", want, got)
		}
	}
}

// TestViewUnsupported confirms non-VCF formats report ErrViewUnsupported rather
// than emitting non-canonical bytes.
func TestViewUnsupported(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "reads.fa")
	if err := os.WriteFile(p, []byte(">seq1\nACGT\n"), 0644); err != nil {
		t.Fatal(err)
	}
	f, err := Identify(p)
	if err != nil {
		t.Fatal(err)
	}
	if err := View(p, f, &bytes.Buffer{}); !errors.Is(err, ErrViewUnsupported) {
		t.Errorf("View(FASTA) error = %v, want ErrViewUnsupported", err)
	}
}
