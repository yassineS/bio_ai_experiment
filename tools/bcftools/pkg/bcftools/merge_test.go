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

// TestMergeInfoDPSum verifies the default INFO rule (DP:sum): a site shared by
// two inputs combines their INFO/DP by summing, while a non-ruled tag (here
// none) is unaffected.
func TestMergeInfoDPSum(t *testing.T) {
	a := mergeVCFOneSample("S1", "chr1\t100\trs1\tA\tT\t.\tPASS\tDP=10\tGT:DP\t0/1:10")
	b := mergeVCFOneSample("S2", "chr1\t100\trs1\tA\tT\t.\tPASS\tDP=25\tGT:DP\t1/1:25")
	_, recs := readMerged(t, []string{a, b}, MergeOptions{})
	if len(recs) != 1 {
		t.Fatalf("expected 1 merged record, got %d", len(recs))
	}
	if got := recs[0].Info["DP"]; got != "35" {
		t.Errorf("INFO/DP = %q, want 35 (10+25 summed)", got)
	}
	// -i - disables rules, so DP falls back to first-input value.
	_, recs2 := readMerged(t, []string{a, b}, MergeOptions{InfoRules: "-"})
	if got := recs2[0].Info["DP"]; got != "10" {
		t.Errorf("INFO/DP with rules disabled = %q, want 10 (first input)", got)
	}
}

// TestMergeDuplicatePairing verifies that intra-position duplicate records are
// paired per-input occurrence (kept as distinct merged rows) rather than
// collapsed, and that a SNP/indel pair at one position keeps file order.
func TestMergeDuplicatePairing(t *testing.T) {
	// Each input holds two records at chr1:100: an indel (file order first) then
	// a SNP. A self-merge must yield two rows, indel before SNP.
	dup := mergeVCFOneSample("S1",
		"chr1\t100\trsIndel\tA\tAG\t.\tPASS\tDP=10\tGT:DP\t0/1:10",
		"chr1\t100\trsSnp\tA\tT\t.\tPASS\tDP=20\tGT:DP\t1/1:20")
	b := mergeVCFOneSample("S2",
		"chr1\t100\trsIndel\tA\tAG\t.\tPASS\tDP=5\tGT:DP\t0/1:5",
		"chr1\t100\trsSnp\tA\tT\t.\tPASS\tDP=7\tGT:DP\t1/1:7")
	_, recs := readMerged(t, []string{dup, b}, MergeOptions{})
	if len(recs) != 2 {
		t.Fatalf("expected 2 merged records (duplicates paired), got %d", len(recs))
	}
	if recs[0].ID != "rsIndel" || recs[1].ID != "rsSnp" {
		t.Errorf("record order = %q,%q, want rsIndel,rsSnp (file order)", recs[0].ID, recs[1].ID)
	}
	if recs[0].Info["DP"] != "15" || recs[1].Info["DP"] != "27" {
		t.Errorf("DP sums = %q,%q, want 15,27", recs[0].Info["DP"], recs[1].Info["DP"])
	}
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

// TestUnitMergeHeadersForceSamples checks the duplicate-sample resolution in
// mergeMergeHeaders without invoking any binary: clashing names from input i
// are prefixed with "<i+1>:", repeatedly if needed, mirroring upstream.
func TestUnitMergeHeadersForceSamples(t *testing.T) {
	mk := func(samples ...string) *vcf.Header {
		return &vcf.Header{
			MetaInfo: []string{"##fileformat=VCFv4.2"},
			Samples:  samples,
		}
	}
	// Without --force-samples a clash is an error.
	if _, _, err := mergeMergeHeaders([]*vcf.Header{mk("A", "B"), mk("A", "C")}, false); err == nil {
		t.Fatalf("expected duplicate-sample error without --force-samples")
	}
	// With --force-samples: A + A -> A, 2:A.
	hdr, renames, err := mergeMergeHeaders([]*vcf.Header{mk("A", "B"), mk("A", "C")}, true)
	if err != nil {
		t.Fatalf("force-samples merge: %v", err)
	}
	want := []string{"A", "B", "2:A", "C"}
	if strings.Join(hdr.Samples, ",") != strings.Join(want, ",") {
		t.Errorf("samples=%v want %v", hdr.Samples, want)
	}
	if renames[1]["A"] != "2:A" {
		t.Errorf("renames[1][A]=%q want 2:A", renames[1]["A"])
	}
	// Three-way clash: A + A + A -> A, 2:A, 3:A.
	hdr3, _, err := mergeMergeHeaders([]*vcf.Header{mk("A"), mk("A"), mk("A")}, true)
	if err != nil {
		t.Fatalf("three-way force-samples merge: %v", err)
	}
	want3 := []string{"A", "2:A", "3:A"}
	if strings.Join(hdr3.Samples, ",") != strings.Join(want3, ",") {
		t.Errorf("three-way samples=%v want %v", hdr3.Samples, want3)
	}
}
