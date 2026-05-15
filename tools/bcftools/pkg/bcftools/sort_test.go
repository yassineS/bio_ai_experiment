package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const sortHeaderText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##contig=<ID=chr10,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
`

func makeSortVCF(records ...string) string {
	return sortHeaderText + strings.Join(records, "\n") + "\n"
}

// Hand-computed: contig declaration order is chr1, chr2, chr10 — so chr10
// records must come AFTER chr2 records, even though chr10 < chr2 lexically.
func TestSortBasic(t *testing.T) {
	in := makeSortVCF(
		"chr10\t100\trsX\tA\tT\t.\tPASS\tDP=1\tGT\t0/1",
		"chr2\t50\trsB\tA\tT\t.\tPASS\tDP=1\tGT\t0/1",
		"chr1\t300\trsC\tA\tT\t.\tPASS\tDP=1\tGT\t0/1",
		"chr1\t100\trsA\tA\tT\t.\tPASS\tDP=1\tGT\t0/1",
	)
	var out bytes.Buffer
	n, err := Sort(strings.NewReader(in), &out, SortOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if n != 4 {
		t.Fatalf("expected 4 records, got %d", n)
	}
	wantOrder := []string{"rsA", "rsC", "rsB", "rsX"}
	prev := -1
	for _, w := range wantOrder {
		idx := strings.Index(out.String(), w)
		if idx < 0 {
			t.Fatalf("missing %s:\n%s", w, out.String())
		}
		if idx <= prev {
			t.Errorf("%s before its predecessor (idx %d <= prev %d):\n%s", w, idx, prev, out.String())
		}
		prev = idx
	}
}

// Records that share (CHROM, POS) order by REF then ALT.
func TestSortStableTies(t *testing.T) {
	in := makeSortVCF(
		"chr1\t100\trsB\tA\tG\t.\tPASS\tDP=1\tGT\t0/1",
		"chr1\t100\trsA\tA\tT\t.\tPASS\tDP=1\tGT\t0/1",
	)
	var out bytes.Buffer
	if _, err := Sort(strings.NewReader(in), &out, SortOptions{}); err != nil {
		t.Fatal(err)
	}
	if strings.Index(out.String(), "rsB") > strings.Index(out.String(), "rsA") {
		t.Errorf("ALT-tie not ordered:\n%s", out.String())
	}
}

func TestSortFilePathInputOutput(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(in, []byte(makeSortVCF(
		"chr2\t50\trsX\tA\tT\t.\tPASS\tDP=1\tGT\t0/1",
		"chr1\t100\trsY\tA\tT\t.\tPASS\tDP=1\tGT\t0/1",
	)), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	n, err := SortFile(in, &out, SortOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if n != 2 {
		t.Errorf("expected 2 records, got %d", n)
	}
	if strings.Index(out.String(), "rsY") > strings.Index(out.String(), "rsX") {
		t.Errorf("chr1/rsY should come before chr2/rsX:\n%s", out.String())
	}
}

func TestSortEmptyOK(t *testing.T) {
	var out bytes.Buffer
	n, err := Sort(strings.NewReader(sortHeaderText), &out, SortOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if n != 0 {
		t.Errorf("empty input got %d records", n)
	}
	// Header should still be emitted.
	if !strings.Contains(out.String(), "#CHROM") {
		t.Error("missing header in output")
	}
}

func TestSortFileMissing(t *testing.T) {
	if _, err := SortFile("/no/such/file", &bytes.Buffer{}, SortOptions{}); err == nil {
		t.Error("expected error for missing file")
	}
}

// Parity-style test: feed a small VCF that mimics the layout upstream
// `bcftools sort` produces, verify our output's record order matches the
// hand-computed expected order.
func TestSortParityHandComputed(t *testing.T) {
	in := makeSortVCF(
		"chr1\t300\trs3\tA\tT\t.\tPASS\tDP=3\tGT\t0/1",
		"chr1\t100\trs1\tA\tT\t.\tPASS\tDP=1\tGT\t0/1",
		"chr2\t250\trs4\tA\tT\t.\tPASS\tDP=4\tGT\t0/1",
		"chr1\t200\trs2\tA\tT\t.\tPASS\tDP=2\tGT\t0/1",
		"chr2\t100\trs0\tA\tT\t.\tPASS\tDP=0\tGT\t0/1",
	)
	var out bytes.Buffer
	if _, err := Sort(strings.NewReader(in), &out, SortOptions{}); err != nil {
		t.Fatal(err)
	}
	// Expected: rs1, rs2, rs3, rs0 (chr2/100), rs4 (chr2/250).
	wantOrder := []string{"rs1", "rs2", "rs3", "rs0", "rs4"}
	prev := -1
	for _, w := range wantOrder {
		idx := strings.Index(out.String(), w)
		if idx <= prev {
			t.Errorf("%s out of order at idx %d (prev %d):\n%s", w, idx, prev, out.String())
		}
		prev = idx
	}
}
