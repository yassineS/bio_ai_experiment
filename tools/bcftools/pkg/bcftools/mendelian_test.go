package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// trioFixture: three samples (CHILD, FATHER, MOTHER), two diploid
// autosomal sites. Site 1 (chr1:10) is consistent: child=0/1 with
// parents 0/0 and 0/1. Site 2 (chr1:20) is inconsistent: child=1/1
// with parents 0/0 and 0/0.
func trioFixture() string {
	return `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000>
##contig=<ID=X,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CHILD	FATHER	MOTHER
chr1	10	.	A	T	.	PASS	DP=1	GT	0/1	0/0	0/1
chr1	20	.	G	C	.	PASS	DP=2	GT	1/1	0/0	0/0
chr1	30	.	C	A	.	PASS	DP=3	GT	./.	0/1	0/1
X	100	.	A	T	.	PASS	DP=4	GT	0/0	1/1	0/0
`
}

func TestParseMendelianMode(t *testing.T) {
	cases := []struct {
		in   string
		want MendelianMode
	}{
		{"", MendelianAnnotate},
		{"a", MendelianAnnotate},
		{"c", MendelianCount},
		{"x", MendelianXChrom},
		{"d", MendelianDelete},
		{"+", MendelianPlusPG},
	}
	for _, c := range cases {
		m, err := ParseMendelianMode(c.in)
		if err != nil {
			t.Errorf("ParseMendelianMode(%q): unexpected error %v", c.in, err)
			continue
		}
		if m != c.want {
			t.Errorf("ParseMendelianMode(%q) = %v want %v", c.in, m, c.want)
		}
	}
	if _, err := ParseMendelianMode("z"); err == nil {
		t.Error("expected error on unknown mode")
	}
}

func TestParseTrioFlag(t *testing.T) {
	tr, err := ParseTrioFlag("kid,dad,mom")
	if err != nil {
		t.Fatalf("ParseTrioFlag: %v", err)
	}
	if tr.Child != "kid" || tr.Father != "dad" || tr.Mother != "mom" {
		t.Errorf("ParseTrioFlag fields: %+v", tr)
	}
	if _, err := ParseTrioFlag("kid,dad"); err == nil {
		t.Error("expected error on bad trio")
	}
	if _, err := ParseTrioFlag("kid,,mom"); err == nil {
		t.Error("expected error on empty member")
	}
}

func TestMendelianConsistent(t *testing.T) {
	// Child = 0/1, Father = 0/0, Mother = 0/1 → consistent.
	if !mendelianConsistent([]int{0, 1}, []int{0, 0}, []int{0, 1}, false) {
		t.Error("expected 0/1 from 0/0 x 0/1 to be consistent")
	}
	// Child = 1/1, Father = 0/0, Mother = 0/0 → inconsistent.
	if mendelianConsistent([]int{1, 1}, []int{0, 0}, []int{0, 0}, false) {
		t.Error("expected 1/1 from 0/0 x 0/0 to be INconsistent")
	}
	// Child = 0/1, Father = 1/1, Mother = 1/1 → inconsistent (child must be 1/1).
	if mendelianConsistent([]int{0, 1}, []int{1, 1}, []int{1, 1}, false) {
		t.Error("expected 0/1 from 1/1 x 1/1 to be INconsistent")
	}
	// X-chrom male child: child 0/0 from father 1/1 and mother 0/0 →
	// consistent under haploidFather (only mother matters).
	if !mendelianConsistent([]int{0, 0}, []int{1, 1}, []int{0, 0}, true) {
		t.Error("expected X-chrom 0/0 with mother 0/0 to be consistent under haploidFather")
	}
	// X-chrom male child: child 0/1 from mother 0/0 → INconsistent.
	if mendelianConsistent([]int{0, 1}, []int{1, 1}, []int{0, 0}, true) {
		t.Error("expected X-chrom 0/1 with mother 0/0 to be INconsistent under haploidFather")
	}
}

func TestParseTrioGT(t *testing.T) {
	cases := []struct {
		in   string
		want []int
		ok   bool
	}{
		{"0/1", []int{0, 1}, true},
		{"1|0", []int{1, 0}, true},
		{"./.", nil, false},
		{"1", nil, false},
		{"0/.", nil, false},
		{"", nil, false},
		{"a/b", nil, false},
	}
	for _, c := range cases {
		got, ok := parseTrioGT(c.in)
		if ok != c.ok {
			t.Errorf("parseTrioGT(%q): ok=%v want %v", c.in, ok, c.ok)
			continue
		}
		if !ok {
			continue
		}
		if len(got) != 2 || got[0] != c.want[0] || got[1] != c.want[1] {
			t.Errorf("parseTrioGT(%q) = %v want %v", c.in, got, c.want)
		}
	}
}

func TestMendelian_AnnotateDefault(t *testing.T) {
	var out bytes.Buffer
	sum, err := Mendelian(strings.NewReader(trioFixture()), &out, MendelianOptions{
		Trios:        []Trio{{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"}},
		Mode:         MendelianAnnotate,
		OutputFormat: OutputVCF,
	})
	if err != nil {
		t.Fatalf("Mendelian: %v", err)
	}
	if sum.TotalRecords != 4 {
		t.Errorf("want 4 records, got %d", sum.TotalRecords)
	}
	// chr1:30 (./.) skipped as missing; X:100 is autosomal-mode → counts.
	// Two error sites: chr1:20 (1/1 from 0/0 x 0/0) and X:100 (0/0 from
	// 1/1 x 0/0 in autosomal mode is inconsistent).
	if sum.RecordsWithError != 2 {
		t.Errorf("want 2 records with error in autosomal mode, got %d", sum.RecordsWithError)
	}
	if len(sum.Trios) != 1 || sum.Trios[0].NError != 2 || sum.Trios[0].NTested != 3 {
		t.Errorf("trio stats: %+v", sum.Trios)
	}
	body := out.String()
	for _, want := range []string{
		`##INFO=<ID=MERR,`,
		"chr1\t10",
		";MERR=0",
		"chr1\t20",
		";MERR=1",
	} {
		if !strings.Contains(body, want) {
			t.Errorf("output missing %q:\n%s", want, body)
		}
	}
}

func TestMendelian_XChromMode(t *testing.T) {
	var out bytes.Buffer
	sum, err := Mendelian(strings.NewReader(trioFixture()), &out, MendelianOptions{
		Trios:        []Trio{{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"}},
		Mode:         MendelianXChrom,
		OutputFormat: OutputVCF,
	})
	if err != nil {
		t.Fatalf("Mendelian: %v", err)
	}
	// X:100: child 0/0, mother 0/0, father 1/1. Under haploidFather
	// only the mother needs to source both alleles, which she can
	// (0/0). So this is *consistent* in X mode and *inconsistent* in
	// autosomal mode. Hence only chr1:20 remains an error.
	if sum.RecordsWithError != 1 {
		t.Errorf("want 1 record with error in X mode, got %d:\n%s", sum.RecordsWithError, out.String())
	}
}

func TestMendelian_CountMode(t *testing.T) {
	var out bytes.Buffer
	sum, err := Mendelian(strings.NewReader(trioFixture()), &out, MendelianOptions{
		Trios: []Trio{{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"}},
		Mode:  MendelianCount,
	})
	if err != nil {
		t.Fatalf("Mendelian: %v", err)
	}
	if sum.TotalRecords != 4 {
		t.Errorf("want 4 records, got %d", sum.TotalRecords)
	}
	body := out.String()
	for _, want := range []string{
		"# trio\tchild\tfather\tmother",
		"1\tCHILD\tFATHER\tMOTHER",
		"# totals\trecords=4\twith_error=2",
	} {
		if !strings.Contains(body, want) {
			t.Errorf("summary missing %q:\n%s", want, body)
		}
	}
}

func TestMendelian_DeleteMode(t *testing.T) {
	var out bytes.Buffer
	_, err := Mendelian(strings.NewReader(trioFixture()), &out, MendelianOptions{
		Trios:        []Trio{{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"}},
		Mode:         MendelianDelete,
		OutputFormat: OutputVCF,
	})
	if err != nil {
		t.Fatalf("Mendelian: %v", err)
	}
	body := out.String()
	// Bad records: chr1:20 (autosomal error), X:100 (autosomal error
	// under default mode). chr1:30 is missing but not an error -> kept.
	if strings.Contains(body, "chr1\t20") {
		t.Errorf("chr1:20 should be deleted in delete mode:\n%s", body)
	}
	if strings.Contains(body, "X\t100") {
		t.Errorf("X:100 should be deleted in delete mode (autosomal default):\n%s", body)
	}
	for _, want := range []string{"chr1\t10", "chr1\t30"} {
		if !strings.Contains(body, want) {
			t.Errorf("output missing kept record %q:\n%s", want, body)
		}
	}
}

// TestMendelian_GoodMode verifies the `-m g` (good-only) mode emits the
// same record set as `-m d` (delete bad). Locks in the upstream-parity
// interpretation that "good-only" === "keep consistent" === "drop bad".
func TestMendelian_GoodMode(t *testing.T) {
	var outGood, outDel bytes.Buffer
	opts := MendelianOptions{
		Trios:        []Trio{{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"}},
		OutputFormat: OutputVCF,
	}
	opts.Mode = MendelianGood
	if _, err := Mendelian(strings.NewReader(trioFixture()), &outGood, opts); err != nil {
		t.Fatalf("Mendelian good: %v", err)
	}
	opts.Mode = MendelianDelete
	if _, err := Mendelian(strings.NewReader(trioFixture()), &outDel, opts); err != nil {
		t.Fatalf("Mendelian delete: %v", err)
	}
	if outGood.String() != outDel.String() {
		t.Fatalf("good-mode output differs from delete-mode:\n--- good\n%s\n--- delete\n%s", outGood.String(), outDel.String())
	}
}

func TestMendelian_LegacyCountFlag(t *testing.T) {
	// -c flag is the legacy alias for -m c.
	var out bytes.Buffer
	_, err := Mendelian(strings.NewReader(trioFixture()), &out, MendelianOptions{
		Trios: []Trio{{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"}},
		Count: true,
	})
	if err != nil {
		t.Fatalf("Mendelian -c: %v", err)
	}
	if !strings.Contains(out.String(), "# totals\trecords=4") {
		t.Errorf("expected TSV summary from -c flag:\n%s", out.String())
	}
}

func TestMendelian_LegacyDeleteFlag(t *testing.T) {
	var out bytes.Buffer
	_, err := Mendelian(strings.NewReader(trioFixture()), &out, MendelianOptions{
		Trios:        []Trio{{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"}},
		Delete:       true,
		OutputFormat: OutputVCF,
	})
	if err != nil {
		t.Fatalf("Mendelian -d: %v", err)
	}
	if strings.Contains(out.String(), "chr1\t20") {
		t.Errorf("chr1:20 should be deleted by legacy -d flag:\n%s", out.String())
	}
}

func TestMendelian_TrioFile(t *testing.T) {
	dir := t.TempDir()
	inPath := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(inPath, []byte(trioFixture()), 0644); err != nil {
		t.Fatal(err)
	}
	tf := filepath.Join(dir, "trios.txt")
	if err := os.WriteFile(tf, []byte("# header\nCHILD\tFATHER\tMOTHER\n"), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	sum, err := MendelianFile(inPath, &out, MendelianOptions{
		TrioFile: tf,
		Mode:     MendelianCount,
	})
	if err != nil {
		t.Fatalf("MendelianFile: %v", err)
	}
	if len(sum.Trios) != 1 || sum.Trios[0].NError != 2 {
		t.Errorf("trios: %+v", sum.Trios)
	}
}

func TestMendelian_MissingTrio(t *testing.T) {
	var out bytes.Buffer
	_, err := Mendelian(strings.NewReader(trioFixture()), &out, MendelianOptions{})
	if err == nil {
		t.Error("expected error when no trio is supplied")
	}
}

func TestMendelian_UnknownSample(t *testing.T) {
	var out bytes.Buffer
	_, err := Mendelian(strings.NewReader(trioFixture()), &out, MendelianOptions{
		Trios: []Trio{{Child: "GHOST", Father: "FATHER", Mother: "MOTHER"}},
	})
	if err == nil {
		t.Error("expected error when child sample is missing")
	}
}

func TestMendelian_MERRPreservesExistingValue(t *testing.T) {
	// Pre-existing MERR=99 should be overwritten with the computed
	// value rather than appended.
	in := strings.ReplaceAll(trioFixture(),
		"DP=2\tGT", "DP=2;MERR=99\tGT")
	var out bytes.Buffer
	_, err := Mendelian(strings.NewReader(in), &out, MendelianOptions{
		Trios:        []Trio{{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"}},
		Mode:         MendelianAnnotate,
		OutputFormat: OutputVCF,
	})
	if err != nil {
		t.Fatalf("Mendelian: %v", err)
	}
	// chr1:20 is the error site; MERR must be 1, not 99.
	for _, line := range strings.Split(out.String(), "\n") {
		if strings.HasPrefix(line, "chr1\t20\t") {
			if !strings.Contains(line, "MERR=1") {
				t.Errorf("expected MERR=1 on chr1:20, got %q", line)
			}
			if strings.Contains(line, "MERR=99") {
				t.Errorf("MERR=99 should have been replaced: %q", line)
			}
			return
		}
	}
	t.Errorf("chr1:20 record not found in output:\n%s", out.String())
}
