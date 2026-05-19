package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func mergeVCFOneSample(sample string, records ...string) string {
	hdr := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	` + sample + "\n"
	return hdr + strings.Join(records, "\n") + "\n"
}

// readMerged is a small helper that runs Merge over in-memory inputs and
// returns the body text plus the parsed merged variants for assertions.
func readMerged(t *testing.T, inputs []string, opts MergeOptions) (string, []*vcf.Variant) {
	t.Helper()
	dir := t.TempDir()
	paths := make([]string, len(inputs))
	for i, s := range inputs {
		p := filepath.Join(dir, "in"+string(rune('0'+i))+".vcf")
		if err := os.WriteFile(p, []byte(s), 0644); err != nil {
			t.Fatal(err)
		}
		paths[i] = p
	}
	var out bytes.Buffer
	if _, err := MergeFiles(paths, &out, opts); err != nil {
		t.Fatalf("MergeFiles: %v", err)
	}
	r := vcf.NewReader(strings.NewReader(out.String()))
	if _, err := r.ReadHeader(); err != nil {
		t.Fatal(err)
	}
	vs, err := r.ReadAll()
	if err != nil {
		t.Fatal(err)
	}
	return out.String(), vs
}

// Hand-computed: two single-sample inputs covering disjoint sites; merged
// header should hold both samples and each record should have GT data for
// both (with `./.` for the missing-from-this-input sample).
func TestMergeBasicTwoSamples(t *testing.T) {
	a := mergeVCFOneSample("S1", "chr1\t100\trsA\tA\tT\t.\tPASS\tDP=10\tGT:DP\t0/1:10")
	b := mergeVCFOneSample("S2", "chr1\t200\trsB\tA\tG\t.\tPASS\tDP=20\tGT:DP\t1/1:20")
	out, recs := readMerged(t, []string{a, b}, MergeOptions{})
	if len(recs) != 2 {
		t.Fatalf("expected 2 merged records, got %d:\n%s", len(recs), out)
	}
	if !strings.Contains(out, "\tS1\tS2\n") {
		t.Errorf("merged header should hold S1 and S2:\n%s", out)
	}
	// First record: rsA, S1 has 0/1, S2 should be ./.
	if recs[0].ID != "rsA" {
		t.Errorf("first record ID = %q, want rsA", recs[0].ID)
	}
	if recs[0].Samples[0].Data["GT"] != "0/1" {
		t.Errorf("S1 GT = %q want 0/1", recs[0].Samples[0].Data["GT"])
	}
	if recs[0].Samples[1].Data["GT"] != "./." {
		t.Errorf("S2 GT = %q want ./.", recs[0].Samples[1].Data["GT"])
	}
}

// Hand-computed: same site (CHROM, POS, REF) in both inputs with different
// ALTs ⇒ multiallelic merge. S1 contributed ALT[0]=T (becomes index 1);
// S2 contributed ALT[0]=G (becomes index 2). GT remapping is verified.
func TestMergeAllelicUnion(t *testing.T) {
	a := mergeVCFOneSample("S1", "chr1\t100\trs1\tA\tT\t.\tPASS\tDP=10\tGT\t0/1")
	b := mergeVCFOneSample("S2", "chr1\t100\trs1\tA\tG\t.\tPASS\tDP=20\tGT\t1/1")
	out, recs := readMerged(t, []string{a, b}, MergeOptions{})
	if len(recs) != 1 {
		t.Fatalf("expected 1 merged record (alt-union), got %d:\n%s", len(recs), out)
	}
	r := recs[0]
	wantAlt := []string{"T", "G"}
	if !sameStringSliceTest(r.Alt, wantAlt) {
		t.Errorf("merged ALT = %v, want %v", r.Alt, wantAlt)
	}
	if r.Samples[0].Data["GT"] != "0/1" {
		t.Errorf("S1 GT = %q (want 0/1, T is allele 1)", r.Samples[0].Data["GT"])
	}
	// S2 used to be 1/1 (ALT=G) but G is now index 2 in the merged record.
	if r.Samples[1].Data["GT"] != "2/2" {
		t.Errorf("S2 GT = %q (want 2/2 after remap)", r.Samples[1].Data["GT"])
	}
}

// Disjoint-sample requirement.
func TestMergeRejectsOverlappingSamples(t *testing.T) {
	a := mergeVCFOneSample("DUP", "chr1\t100\trsA\tA\tT\t.\tPASS\tDP=10\tGT\t0/1")
	b := mergeVCFOneSample("DUP", "chr1\t200\trsB\tA\tT\t.\tPASS\tDP=10\tGT\t0/1")
	dir := t.TempDir()
	p1 := filepath.Join(dir, "a.vcf")
	p2 := filepath.Join(dir, "b.vcf")
	if err := os.WriteFile(p1, []byte(a), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(p2, []byte(b), 0644); err != nil {
		t.Fatal(err)
	}
	if _, err := MergeFiles([]string{p1, p2}, &bytes.Buffer{}, MergeOptions{}); err == nil {
		t.Error("expected error for overlapping sample sets")
	}
}

// FILTER union: PASS only kept when every input was PASS, otherwise the
// non-PASS values are unioned.
func TestMergeFilters(t *testing.T) {
	a := mergeVCFOneSample("S1", "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1")
	b := mergeVCFOneSample("S2", "chr1\t100\trs1\tA\tT\t.\tLowQ\t.\tGT\t0/1")
	_, recs := readMerged(t, []string{a, b}, MergeOptions{})
	if len(recs) != 1 {
		t.Fatalf("expected 1 record")
	}
	if !sameStringSliceTest(recs[0].Filter, []string{"LowQ"}) {
		t.Errorf("FILTER union = %v want [LowQ]", recs[0].Filter)
	}
}

// Test the merge-mode = none path: identical (CHROM, POS, REF, ALT) collapse
// but different ALTs at the same POS DO NOT collapse.
func TestMergeModeNone(t *testing.T) {
	a := mergeVCFOneSample("S1", "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/1")
	b := mergeVCFOneSample("S2", "chr1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1")
	_, recs := readMerged(t, []string{a, b}, MergeOptions{MergeMode: MergeNone})
	if len(recs) != 2 {
		t.Errorf("MergeNone should keep two records (different ALTs); got %d", len(recs))
	}
}

// File-list flag.
func TestMergeFileList(t *testing.T) {
	dir := t.TempDir()
	p1 := filepath.Join(dir, "a.vcf")
	p2 := filepath.Join(dir, "b.vcf")
	list := filepath.Join(dir, "files.txt")
	if err := os.WriteFile(p1, []byte(mergeVCFOneSample("S1", "chr1\t1\trsA\tA\tT\t.\tPASS\t.\tGT\t0/1")), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(p2, []byte(mergeVCFOneSample("S2", "chr1\t2\trsB\tA\tT\t.\tPASS\t.\tGT\t0/1")), 0644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(list, []byte(p1+"\n"+p2+"\n"), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := MergeFiles(nil, &out, MergeOptions{FileList: list}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out.String(), "rsA") || !strings.Contains(out.String(), "rsB") {
		t.Errorf("missing records:\n%s", out.String())
	}
}

func TestMergeNeedsTwoFiles(t *testing.T) {
	if _, err := MergeFiles([]string{"only-one.vcf"}, &bytes.Buffer{}, MergeOptions{}); err == nil {
		t.Error("expected error for single input")
	}
}

func TestParseMergeMode(t *testing.T) {
	cases := []struct {
		in   string
		want MergeMode
		ok   bool
	}{
		{"", MergeBoth, true},
		{"both", MergeBoth, true},
		{"none", MergeNone, true},
		{"snps", MergeSNPs, true},
		{"indels", MergeIndels, true},
		{"all", MergeAll, true},
		{"id", MergeID, true},
		{"bogus", MergeBoth, false},
	}
	for _, c := range cases {
		got, err := ParseMergeMode(c.in)
		if c.ok {
			if err != nil {
				t.Errorf("ParseMergeMode(%q) err: %v", c.in, err)
				continue
			}
			if got != c.want {
				t.Errorf("ParseMergeMode(%q) = %v want %v", c.in, got, c.want)
			}
		} else {
			if err == nil {
				t.Errorf("ParseMergeMode(%q) should have errored", c.in)
			}
		}
	}
}

func TestRemapGTByMap(t *testing.T) {
	m := map[int]int{0: 0, 1: 2}
	if got := remapGTByMap("0/1", m); got != "0/2" {
		t.Errorf("remapGTByMap(0/1) = %q want 0/2", got)
	}
	if got := remapGTByMap("1|0", m); got != "2|0" {
		t.Errorf("remapGTByMap(1|0) = %q want 2|0", got)
	}
	if got := remapGTByMap("./.", m); got != "./." {
		t.Errorf("remapGTByMap(./.) = %q want ./.", got)
	}
}

func TestIsSNPRecord(t *testing.T) {
	v := &vcf.Variant{Ref: "A", Alt: []string{"T"}}
	if !isSNPRecord(v) {
		t.Error("A->T should be a SNP")
	}
	v2 := &vcf.Variant{Ref: "AT", Alt: []string{"A"}}
	if isSNPRecord(v2) {
		t.Error("AT->A should not be a SNP")
	}
}

func sameStringSliceTest(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
