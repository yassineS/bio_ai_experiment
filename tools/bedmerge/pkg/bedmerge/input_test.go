package bedmerge

import (
	"bytes"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"testing"
)

// repoRootDir walks up from this test file to the module root (the dir holding
// go.mod), used to locate upstream BAM fixtures.
func repoRootDir(t *testing.T) string {
	t.Helper()
	_, file, _, _ := runtime.Caller(0)
	dir := filepath.Dir(file)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatalf("could not locate repo root from %s", file)
		}
		dir = parent
	}
}

// TestVCFStructuralVariantLength covers the SV-length path in vcfEnd/vcfSVLen:
// a <DEL> with SVLEN takes the abs-max magnitude; an <DUP> with only END uses
// END-POS+1; an <INS...> is zero length; a plain ALT uses len(REF).
func TestVCFStructuralVariantLength(t *testing.T) {
	vcf := strings.Join([]string{
		"##fileformat=VCFv4.1",
		"chr1\t100\tdel\tG\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=-50,-200;END=300",
		"chr1\t500\tdup\tC\t<DUP>\t.\t.\tSVTYPE=DUP;END=600",
		"chr1\t800\tins\tA\t<INS:ME>\t.\t.\tSVTYPE=INS;SVLEN=120",
		"chr1\t1000\tsnv\tACGT\tA\t.\t.\t.",
	}, "\n") + "\n"

	var buf bytes.Buffer
	if _, err := Merge(strings.NewReader(vcf), &buf, MergeOptions{}); err != nil {
		t.Fatalf("Merge VCF SV failed: %v", err)
	}
	// del: start 99, len abs-max(50,200)=200 -> 99..299
	// dup: start 499, len END-POS+1 = 600-500+1 = 101 -> 499..600
	// ins: start 799, len 0 -> 799..799
	// snv: start 999, len(REF=ACGT)=4 -> 999..1003
	want := "chr1\t99\t299\n" +
		"chr1\t499\t600\n" +
		"chr1\t799\t799\n" +
		"chr1\t999\t1003\n"
	if buf.String() != want {
		t.Fatalf("VCF SV output mismatch.\nwant:\n%q\ngot:\n%q", want, buf.String())
	}
}

// TestGFFInputAutoDetect covers the GFF branch of parseTextRecord: 1-based
// inclusive coordinates become 0-based half-open, and the feature merges.
func TestGFFInputAutoDetect(t *testing.T) {
	gff := "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=g1\n" +
		"chr1\tsrc\texon\t150\t250\t.\t+\t.\tID=e1\n"
	var buf bytes.Buffer
	if _, err := Merge(strings.NewReader(gff), &buf, MergeOptions{}); err != nil {
		t.Fatalf("Merge GFF failed: %v", err)
	}
	// 100..200 -> 99..200, 150..250 -> 149..250; merged 99..250.
	if got := buf.String(); got != "chr1\t99\t250\n" {
		t.Fatalf("GFF output mismatch: got %q", got)
	}
}

// TestUnexpectedFormatError covers detectFormat's failure path.
func TestUnexpectedFormatError(t *testing.T) {
	if _, err := Merge(strings.NewReader("not\ta\tbed\n"), &bytes.Buffer{}, MergeOptions{}); err == nil {
		t.Fatalf("expected unexpected-format error, got nil")
	}
}

// TestScientificNotationCoords covers parseChrPos's strtod fallback.
func TestScientificNotationCoords(t *testing.T) {
	var buf bytes.Buffer
	if _, err := Merge(strings.NewReader("chr1\t8e02\t830\n"), &buf, MergeOptions{}); err != nil {
		t.Fatalf("Merge scientific coords failed: %v", err)
	}
	if got := buf.String(); got != "chr1\t800\t830\n" {
		t.Fatalf("scientific coord mismatch: got %q", got)
	}
}

// TestWriteHeaderEchoesComments covers WriteHeader: comment/track/browser lines
// are echoed and the first data line stops the copy.
func TestWriteHeaderEchoesComments(t *testing.T) {
	in := "#chr\tstart\tstop\ntrack name=x\nbrowser y\nchr1\t10\t20\n"
	var buf bytes.Buffer
	if err := WriteHeader(strings.NewReader(in), &buf); err != nil {
		t.Fatalf("WriteHeader failed: %v", err)
	}
	want := "#chr\tstart\tstop\ntrack name=x\nbrowser y\n"
	if buf.String() != want {
		t.Fatalf("WriteHeader mismatch.\nwant:\n%q\ngot:\n%q", want, buf.String())
	}
}

// TestValidateBAMColumns exercises the BAM column-constraint checks directly.
func TestValidateBAMColumns(t *testing.T) {
	if err := validateBAMColumns(fmtBED, &ColumnOps{Columns: []int{2}}); err != nil {
		t.Errorf("non-BAM input should not validate columns: %v", err)
	}
	co := &ColumnOps{Columns: []int{12}}
	err := validateBAMColumns(fmtBAM, co)
	be, ok := err.(*BAMColumnError)
	if !ok || be.Column != 12 || be.Flags {
		t.Errorf("expected out-of-range BAMColumnError(12), got %v", err)
	}
	if be.Error() == "" {
		t.Error("BAMColumnError.Error() should be non-empty")
	}
	flagsErr := validateBAMColumns(fmtBAM, &ColumnOps{Columns: []int{2}})
	if fe, ok := flagsErr.(*BAMColumnError); !ok || !fe.Flags {
		t.Errorf("expected Flags BAMColumnError, got %v", flagsErr)
	} else if fe.Error() == "" {
		t.Error("flags BAMColumnError.Error() should be non-empty")
	}
}

// TestBAMInputColumnOps covers the BAM read path (readBAM, bamToRecord,
// bamCigarStr, bamQualStr) using the upstream fullFields.bam fixture. It asserts
// our output equals the recorded expected file bamCol1Collapse.txt.
func TestBAMInputColumnOps(t *testing.T) {
	root := repoRootDir(t)
	dir := filepath.Join(root, "reference_code", "bedtools", "test", "merge")
	bamPath := filepath.Join(dir, "fullFields.bam")
	if _, err := os.Stat(bamPath); err != nil {
		t.Fatalf("upstream BAM fixture unavailable: %v\n"+
			"run: git submodule update --init reference_code/bedtools", err)
	}
	f, err := os.Open(bamPath)
	if err != nil {
		t.Fatalf("open BAM: %v", err)
	}
	defer f.Close()

	co, err := ParseColumnOps("1", "collapse")
	if err != nil {
		t.Fatalf("ParseColumnOps: %v", err)
	}
	var buf bytes.Buffer
	if _, err := Merge(f, &buf, MergeOptions{ColumnOps: co}); err != nil {
		t.Fatalf("Merge BAM failed: %v", err)
	}
	want, err := os.ReadFile(filepath.Join(dir, "bamCol1Collapse.txt"))
	if err != nil {
		t.Fatalf("read expected: %v", err)
	}
	if !bytes.Equal(buf.Bytes(), want) {
		t.Fatalf("BAM col1 collapse mismatch.\nwant:\n%s\ngot:\n%s", want, buf.String())
	}
}

// TestBAMCigarAndQualOps covers the CIGAR (col 6) and QUAL (col 11) BAM fields.
func TestBAMCigarAndQualOps(t *testing.T) {
	root := repoRootDir(t)
	dir := filepath.Join(root, "reference_code", "bedtools", "test", "merge")
	cases := []struct{ col, expected string }{
		{"6", "bamCol6Collapse.txt"},
		{"11", "bamCol11Collapse.txt"},
	}
	for _, c := range cases {
		f, err := os.Open(filepath.Join(dir, "fullFields.bam"))
		if err != nil {
			t.Fatalf("open BAM: %v", err)
		}
		co, _ := ParseColumnOps(c.col, "collapse")
		var buf bytes.Buffer
		if _, err := Merge(f, &buf, MergeOptions{ColumnOps: co}); err != nil {
			f.Close()
			t.Fatalf("Merge BAM col %s: %v", c.col, err)
		}
		f.Close()
		want, err := os.ReadFile(filepath.Join(dir, c.expected))
		if err != nil {
			t.Fatalf("read expected %s: %v", c.expected, err)
		}
		if !bytes.Equal(buf.Bytes(), want) {
			t.Fatalf("BAM col %s mismatch.\nwant:\n%s\ngot:\n%s", c.col, want, buf.String())
		}
	}
}
