package bcftools

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const concatHeader = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
`

func concatVCF(records ...string) string {
	return concatHeader + strings.Join(records, "\n") + "\n"
}

func runConcat(t *testing.T, inputs []string, opts ConcatOptions) (string, int) {
	t.Helper()
	readers := make([]NamedReader, len(inputs))
	for i, s := range inputs {
		readers[i] = NamedReader{Name: "in" + string(rune('0'+i)), Reader: strings.NewReader(s)}
	}
	var out bytes.Buffer
	n, err := Concat(readers, &out, opts)
	if err != nil {
		t.Fatalf("Concat: %v", err)
	}
	return out.String(), n
}

func TestConcatNonOverlapping(t *testing.T) {
	a := concatVCF(
		"chr1\t100\trs1\tA\tT\t30\tPASS\tDP=10\tGT\t0/1",
	)
	b := concatVCF(
		"chr1\t200\trs2\tC\tG\t30\tPASS\tDP=20\tGT\t1/1",
	)
	c := concatVCF(
		"chr2\t50\trs3\tG\tA\t30\tPASS\tDP=30\tGT\t0/0",
	)
	out, n := runConcat(t, []string{a, b, c}, ConcatOptions{})
	if n != 3 {
		t.Fatalf("expected 3 records, got %d", n)
	}
	for _, want := range []string{"rs1", "rs2", "rs3"} {
		if !strings.Contains(out, want) {
			t.Errorf("output missing %q:\n%s", want, out)
		}
	}
	// rs1 should appear before rs3 in the concatenation order.
	if strings.Index(out, "rs1") > strings.Index(out, "rs3") {
		t.Errorf("records out of order:\n%s", out)
	}
}

func TestConcatAllowOverlaps(t *testing.T) {
	a := concatVCF(
		"chr1\t100\trs1\tA\tT\t30\tPASS\tDP=10\tGT\t0/1",
		"chr1\t400\trs4\tC\tG\t30\tPASS\tDP=10\tGT\t0/1",
	)
	b := concatVCF(
		"chr1\t200\trs2\tA\tT\t30\tPASS\tDP=10\tGT\t0/1",
		"chr1\t300\trs3\tA\tT\t30\tPASS\tDP=10\tGT\t0/1",
	)
	out, _ := runConcat(t, []string{a, b}, ConcatOptions{AllowOverlaps: true})
	// Records should now be globally sorted: rs1, rs2, rs3, rs4.
	wantOrder := []string{"rs1", "rs2", "rs3", "rs4"}
	prev := -1
	for _, w := range wantOrder {
		idx := strings.Index(out, w)
		if idx < 0 {
			t.Errorf("missing %s in:\n%s", w, out)
			continue
		}
		if idx <= prev {
			t.Errorf("sort-merge order wrong, %s at %d before previous %d:\n%s", w, idx, prev, out)
		}
		prev = idx
	}
}

func TestConcatRemoveDuplicates(t *testing.T) {
	a := concatVCF(
		"chr1\t100\trs1\tA\tT\t30\tPASS\tDP=10\tGT\t0/1",
	)
	b := concatVCF(
		"chr1\t100\trs1\tA\tT\t30\tPASS\tDP=10\tGT\t0/1",
		"chr1\t200\trs2\tA\tT\t30\tPASS\tDP=10\tGT\t0/1",
	)
	out, n := runConcat(t, []string{a, b}, ConcatOptions{AllowOverlaps: true, RemoveDuplicates: true})
	if n != 2 {
		t.Fatalf("expected 2 records after dedup, got %d:\n%s", n, out)
	}
	if strings.Count(out, "\trs1\t") != 1 {
		t.Errorf("rs1 should appear once after dedup:\n%s", out)
	}
}

func TestConcatHeaderMergeContigUnion(t *testing.T) {
	a := `##fileformat=VCFv4.2
##contig=<ID=chrA,length=10>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chrA	1	.	A	T	.	PASS	DP=10
`
	b := `##fileformat=VCFv4.2
##contig=<ID=chrB,length=10>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chrB	1	.	C	G	.	PASS	DP=20
`
	out, _ := runConcat(t, []string{a, b}, ConcatOptions{})
	if !strings.Contains(out, "##contig=<ID=chrA") || !strings.Contains(out, "##contig=<ID=chrB") {
		t.Errorf("expected both contigs in merged header:\n%s", out)
	}
	// Only one INFO=DP line should appear despite both files declaring it.
	if strings.Count(out, "##INFO=<ID=DP") != 1 {
		t.Errorf("duplicate INFO/DP not deduplicated:\n%s", out)
	}
}

func TestConcatConflictingDefinitionsError(t *testing.T) {
	a := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1	.	A	T	.	PASS	DP=10
`
	// Note: Number=2 instead of Number=1 — this should be flagged.
	b := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10>
##INFO=<ID=DP,Number=2,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	2	.	A	T	.	PASS	DP=10,20
`
	readers := []NamedReader{
		{Name: "a", Reader: strings.NewReader(a)},
		{Name: "b", Reader: strings.NewReader(b)},
	}
	var out bytes.Buffer
	_, err := Concat(readers, &out, ConcatOptions{})
	if err == nil {
		t.Fatal("expected error for conflicting INFO definitions")
	}
	if !strings.Contains(err.Error(), "INFO") || !strings.Contains(err.Error(), "DP") {
		t.Errorf("error should mention INFO and DP: %v", err)
	}
}

func TestConcatSampleMismatchError(t *testing.T) {
	a := concatHeader + "chr1\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
	b := strings.Replace(concatHeader, "S1", "S2", 1) + "chr1\t2\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
	readers := []NamedReader{
		{Name: "a", Reader: strings.NewReader(a)},
		{Name: "b", Reader: strings.NewReader(b)},
	}
	var out bytes.Buffer
	_, err := Concat(readers, &out, ConcatOptions{})
	if err == nil {
		t.Fatal("expected error for sample mismatch")
	}
}

func TestConcatFileList(t *testing.T) {
	dir := t.TempDir()
	a := filepath.Join(dir, "a.vcf")
	b := filepath.Join(dir, "b.vcf")
	list := filepath.Join(dir, "files.txt")
	if err := os.WriteFile(a, []byte(concatVCF("chr1\t1\trsA\tA\tT\t.\tPASS\tDP=1\tGT\t0/1")), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(b, []byte(concatVCF("chr2\t1\trsB\tA\tT\t.\tPASS\tDP=2\tGT\t0/0")), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(list, []byte("# header\n"+a+"\n"+b+"\n"), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	n, err := ConcatFiles(nil, &out, ConcatOptions{FileList: list})
	if err != nil {
		t.Fatal(err)
	}
	if n != 2 {
		t.Fatalf("expected 2 records, got %d", n)
	}
	for _, w := range []string{"rsA", "rsB"} {
		if !strings.Contains(out.String(), w) {
			t.Errorf("missing %s:\n%s", w, out.String())
		}
	}
}

func TestConcatBCFRoundTrip(t *testing.T) {
	// Build two VCF files, concat to BCF, and read back through ViewFile to
	// verify the round-trip.
	dir := t.TempDir()
	a := filepath.Join(dir, "a.vcf")
	b := filepath.Join(dir, "b.vcf")
	if err := os.WriteFile(a, []byte(concatVCF("chr1\t1\trsA\tA\tT\t.\tPASS\tDP=1\tGT\t0/1")), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(b, []byte(concatVCF("chr2\t1\trsB\tA\tT\t.\tPASS\tDP=2\tGT\t0/0")), 0644); err != nil {
		t.Fatal(err)
	}
	outPath := filepath.Join(dir, "merged.bcf")
	f, err := os.Create(outPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := ConcatFiles([]string{a, b}, f, ConcatOptions{OutputFormat: OutputBCF}); err != nil {
		t.Fatalf("ConcatFiles -O b: %v", err)
	}
	f.Close()

	var back bytes.Buffer
	if _, err := ViewFile(outPath, &back, ViewOptions{}, io.Discard); err != nil {
		t.Fatalf("ViewFile(bcf): %v", err)
	}
	for _, w := range []string{"rsA", "rsB", "##contig=<ID=chr1", "##contig=<ID=chr2"} {
		if !strings.Contains(back.String(), w) {
			t.Errorf("BCF round-trip missing %q:\n%s", w, back.String())
		}
	}
}

func TestConcatNoInputs(t *testing.T) {
	if _, err := Concat(nil, &bytes.Buffer{}, ConcatOptions{}); err == nil {
		t.Error("expected error for no inputs")
	}
	if _, err := ConcatFiles(nil, &bytes.Buffer{}, ConcatOptions{}); err == nil {
		t.Error("expected error for no inputs via ConcatFiles")
	}
}

func TestConcatFileListBadPath(t *testing.T) {
	if _, err := ReadFileList("/no/such/file"); err == nil {
		t.Error("expected error for missing file list")
	}
}

func TestStructuredIDParses(t *testing.T) {
	k, id := structuredID(`##INFO=<ID=DP,Number=1,Type=Integer,Description="Read, depth">`)
	if k != "INFO" || id != "DP" {
		t.Errorf("INFO/DP: got %q/%q", k, id)
	}
	k, id = structuredID(`##contig=<ID=chr1,length=1000>`)
	if k != "contig" || id != "chr1" {
		t.Errorf("contig/chr1: got %q/%q", k, id)
	}
	k, _ = structuredID(`##source=v1.0`)
	if k != "" {
		t.Errorf("unstructured line: got kind %q", k)
	}
}

func TestConcatOrderPreserved(t *testing.T) {
	// Without AllowOverlaps the input order is preserved verbatim.
	a := concatVCF("chr2\t1\trsA\tA\tT\t.\tPASS\tDP=1\tGT\t0/1")
	b := concatVCF("chr1\t1\trsB\tA\tT\t.\tPASS\tDP=2\tGT\t0/0")
	out, _ := runConcat(t, []string{a, b}, ConcatOptions{})
	if strings.Index(out, "rsA") > strings.Index(out, "rsB") {
		t.Errorf("input order not preserved:\n%s", out)
	}
}

func TestConcatFilesOpenError(t *testing.T) {
	if _, err := ConcatFiles([]string{"/no/such/file"}, &bytes.Buffer{}, ConcatOptions{}); err == nil {
		t.Error("expected error opening missing file")
	}
}

func TestConcatFilesFileListInvalid(t *testing.T) {
	if _, err := ConcatFiles(nil, &bytes.Buffer{}, ConcatOptions{FileList: "/no/such/list"}); err == nil {
		t.Error("expected error for missing file list")
	}
}

// TestConcatBCFInputStream verifies that BCF inputs are accepted and
// converted back to VCF on output.
func TestConcatBCFInputStream(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1	rsX	A	T	.	PASS	DP=1	GT	0/1
`
	dir := t.TempDir()
	bcfPath := filepath.Join(dir, "x.bcf")
	f, err := os.Create(bcfPath)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := View(strings.NewReader(vcfText), f, ViewOptions{OutputFormat: OutputBCF}); err != nil {
		t.Fatal(err)
	}
	f.Close()
	var out bytes.Buffer
	n, err := ConcatFiles([]string{bcfPath}, &out, ConcatOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if n != 1 {
		t.Fatalf("expected 1 record, got %d", n)
	}
	if !strings.Contains(out.String(), "rsX") {
		t.Errorf("missing rsX in BCF concat output:\n%s", out.String())
	}
}

// TestConcatSortFallback exercises the unknown-contig hash ordering.
func TestConcatSortFallback(t *testing.T) {
	if sortFallback("abc") == sortFallback("xyz") {
		t.Error("sortFallback collision for different names — extremely unlikely")
	}
	if sortFallback("foo") != sortFallback("foo") {
		t.Error("sortFallback should be deterministic")
	}
}

// TestConcatAllowOverlapsUnknownContig exercises the AllowOverlaps path when
// some records use contigs not declared in any header — they should still
// merge in a stable order.
func TestConcatAllowOverlapsUnknownContig(t *testing.T) {
	const noContig = `##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
unkA	5	a	A	T	.	PASS	DP=1
`
	const noContig2 = `##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
unkB	5	b	A	T	.	PASS	DP=1
`
	out, _ := runConcat(t, []string{noContig, noContig2}, ConcatOptions{AllowOverlaps: true})
	if !strings.Contains(out, "unkA") || !strings.Contains(out, "unkB") {
		t.Errorf("expected both unkA and unkB:\n%s", out)
	}
}

func TestSplitTopLevelQuotedComma(t *testing.T) {
	got := splitTopLevel(`ID=DP,Number=1,Description="hello, world",Type=Integer`)
	if len(got) != 4 {
		t.Errorf("got %d parts: %v", len(got), got)
	}
	if got[2] != `Description="hello, world"` {
		t.Errorf("comma-in-quotes split: got %q", got[2])
	}
}

func TestMergeHeadersEmpty(t *testing.T) {
	if _, err := MergeHeaders(nil); err == nil {
		t.Error("expected error for empty input")
	}
}
