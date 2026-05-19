package vcftools

import (
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func TestParseBiallelicPL(t *testing.T) {
	cases := []struct {
		in      string
		want    [3]float64
		present bool
	}{
		{"0,10,100", [3]float64{0, 10, 100}, true},
		{"5.5,20,80.5", [3]float64{5.5, 20, 80.5}, true},
		{"", [3]float64{}, false},
		{".", [3]float64{}, false},
		{"0,10", [3]float64{}, false},
		{"0,.,10", [3]float64{}, false},
		{"a,b,c", [3]float64{}, false},
	}
	for _, c := range cases {
		got, ok := parseBiallelicPL(c.in)
		if ok != c.present {
			t.Errorf("parseBiallelicPL(%q): present = %v, want %v", c.in, ok, c.present)
			continue
		}
		if ok && got != c.want {
			t.Errorf("parseBiallelicPL(%q) = %v, want %v", c.in, got, c.want)
		}
	}
}

func TestFormatHelpers(t *testing.T) {
	if got := formatGL(0); got != "0" {
		t.Errorf("formatGL(0) = %q", got)
	}
	if got := formatGL(-1.5); got != "-1.500000" {
		t.Errorf("formatGL(-1.5) = %q", got)
	}
	if got := formatPL(0); got != "0" {
		t.Errorf("formatPL(0) = %q", got)
	}
	if got := formatPL(10); got != "10" {
		t.Errorf("formatPL(10) = %q", got)
	}
	if got := formatPL(2.5); got != "2.5" {
		t.Errorf("formatPL(2.5) = %q", got)
	}
}

func TestIsSimpleSNP(t *testing.T) {
	cases := []struct {
		ref, alt string
		want     bool
	}{
		{"A", "G", true},
		{"a", "g", true},
		{"N", "A", true},
		{"AT", "A", false},
		{"A", "AT", false},
		{"A", "<DEL>", false},
		{"A", "*", false},
	}
	for _, c := range cases {
		if got := isSimpleSNP(c.ref, c.alt); got != c.want {
			t.Errorf("isSimpleSNP(%q,%q) = %v, want %v", c.ref, c.alt, got, c.want)
		}
	}
}

func TestFormatContains(t *testing.T) {
	if !formatContains([]string{"GT", "PL"}, "PL") {
		t.Errorf("expected PL present")
	}
	if formatContains([]string{"GT"}, "PL") {
		t.Errorf("expected PL absent")
	}
}

// writeBeagle emits a single variant through a beagleWriter and returns the
// written file contents.
func writeBeagle(t *testing.T, mode beagleMode, samples []string, vs ...*vcf.Variant) string {
	t.Helper()
	tmp := t.TempDir()
	prefix := filepath.Join(tmp, "b")
	w, err := newBEAGLEWriter(prefix, mode)
	if err != nil {
		t.Fatalf("newBEAGLEWriter: %v", err)
	}
	for _, v := range vs {
		if err := w.write(v, samples); err != nil {
			t.Fatalf("write: %v", err)
		}
	}
	if err := w.close(); err != nil {
		t.Fatalf("close: %v", err)
	}
	suffix := ".BEAGLE.GL"
	if mode == beaglePL {
		suffix = ".BEAGLE.PL"
	}
	got, err := os.ReadFile(prefix + suffix)
	if err != nil {
		t.Fatalf("read: %v", err)
	}
	return string(got)
}

func variantPL(chrom string, pos int, ref string, alt []string, samples []vcf.Sample) *vcf.Variant {
	return &vcf.Variant{
		Chrom:   chrom,
		Pos:     pos,
		ID:      ".",
		Ref:     ref,
		Alt:     alt,
		Qual:    30,
		Filter:  []string{"PASS"},
		Info:    map[string]string{},
		Format:  []string{"GT", "PL"},
		Samples: samples,
	}
}

func TestBeagleGLConversion(t *testing.T) {
	v := variantPL("chr1", 100, "A", []string{"G"}, []vcf.Sample{
		{Name: "s1", Data: map[string]string{"GT": "0/0", "PL": "0,10,100"}},
		{Name: "s2", Data: map[string]string{"GT": "0/1", "PL": "10,0,10"}},
		{Name: "s3", Data: map[string]string{"GT": "./.", "PL": "."}}, // missing
	})
	out := writeBeagle(t, beagleGL, []string{"s1", "s2", "s3"}, v)
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 2 {
		t.Fatalf("got %d lines, want 2:\n%s", len(lines), out)
	}
	headerFields := strings.Fields(lines[0])
	wantHeader := []string{"marker", "allele1", "allele2", "s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3", "s3"}
	if len(headerFields) != len(wantHeader) {
		t.Fatalf("header fields = %d, want %d: %v", len(headerFields), len(wantHeader), headerFields)
	}
	for i := range headerFields {
		if headerFields[i] != wantHeader[i] {
			t.Errorf("header[%d] = %q, want %q", i, headerFields[i], wantHeader[i])
		}
	}
	row := strings.Fields(lines[1])
	if row[0] != "chr1:100" || row[1] != "A" || row[2] != "G" {
		t.Errorf("row prefix = %v", row[:3])
	}
	// s1 PL=0,10,100  -> GL=0, -1, -10
	if row[3] != "0" || row[4] != "-1.000000" || row[5] != "-10.000000" {
		t.Errorf("s1 GL triplet = %v %v %v", row[3], row[4], row[5])
	}
	// s2 PL=10,0,10 -> GL=-1, 0, -1
	if row[6] != "-1.000000" || row[7] != "0" || row[8] != "-1.000000" {
		t.Errorf("s2 GL triplet = %v %v %v", row[6], row[7], row[8])
	}
	// s3 missing -> log10(1/3) ≈ -0.477121
	for i := 9; i < 12; i++ {
		gl, err := strconvFloat(row[i])
		if err != nil {
			t.Fatalf("parse %q: %v", row[i], err)
		}
		if math.Abs(gl-log10OneThird) > 1e-6 {
			t.Errorf("s3 GL[%d] = %v, want ≈ %v", i-9, gl, log10OneThird)
		}
	}
}

func TestBeaglePLRaw(t *testing.T) {
	v := variantPL("chr1", 100, "A", []string{"G"}, []vcf.Sample{
		{Name: "s1", Data: map[string]string{"GT": "0/0", "PL": "0,10,100"}},
		{Name: "s2", Data: map[string]string{"GT": "0/1", "PL": "."}},
	})
	out := writeBeagle(t, beaglePL, []string{"s1", "s2"}, v)
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 2 {
		t.Fatalf("got %d lines: %s", len(lines), out)
	}
	row := strings.Fields(lines[1])
	if row[3] != "0" || row[4] != "10" || row[5] != "100" {
		t.Errorf("s1 PL = %v %v %v", row[3], row[4], row[5])
	}
	if row[6] != "0" || row[7] != "0" || row[8] != "0" {
		t.Errorf("s2 PL (missing) = %v %v %v", row[6], row[7], row[8])
	}
}

func TestBeagleSkipsIndelsMultiAllelicAndNoPL(t *testing.T) {
	tmp := t.TempDir()
	prefix := filepath.Join(tmp, "b")
	w, err := newBEAGLEWriter(prefix, beaglePL)
	if err != nil {
		t.Fatalf("new: %v", err)
	}
	samples := []string{"s1"}
	// indel
	v1 := variantPL("chr1", 100, "A", []string{"AT"}, []vcf.Sample{
		{Name: "s1", Data: map[string]string{"GT": "0/1", "PL": "0,10,100"}},
	})
	// multi-allelic
	v2 := variantPL("chr1", 200, "A", []string{"G", "C"}, []vcf.Sample{
		{Name: "s1", Data: map[string]string{"GT": "0/1", "PL": "0,10,100"}},
	})
	// no PL field
	v3 := &vcf.Variant{
		Chrom: "chr1", Pos: 300, ID: ".", Ref: "A", Alt: []string{"G"},
		Qual: 30, Filter: []string{"PASS"}, Info: map[string]string{},
		Format:  []string{"GT"},
		Samples: []vcf.Sample{{Name: "s1", Data: map[string]string{"GT": "0/1"}}},
	}
	// valid biallelic SNP with PL — should emit
	v4 := variantPL("chr1", 400, "A", []string{"G"}, []vcf.Sample{
		{Name: "s1", Data: map[string]string{"GT": "0/1", "PL": "0,10,100"}},
	})
	for _, v := range []*vcf.Variant{v1, v2, v3, v4} {
		if err := w.write(v, samples); err != nil {
			t.Fatalf("write: %v", err)
		}
	}
	if err := w.close(); err != nil {
		t.Fatalf("close: %v", err)
	}
	body, err := os.ReadFile(prefix + ".BEAGLE.PL")
	if err != nil {
		t.Fatalf("read: %v", err)
	}
	lines := strings.Split(strings.TrimRight(string(body), "\n"), "\n")
	// One header + only v4 should make it through.
	if len(lines) != 2 {
		t.Fatalf("expected 2 lines (header + 1 site), got %d:\n%s", len(lines), body)
	}
	if !strings.HasPrefix(lines[1], "chr1:400\t") && !strings.HasPrefix(lines[1], "chr1:400 ") {
		t.Errorf("expected only chr1:400 row, got: %s", lines[1])
	}
}

func TestBeagleNilCloseSafe(t *testing.T) {
	var b *beagleWriter
	if err := b.close(); err != nil {
		t.Errorf("nil close: %v", err)
	}
}

func TestRunBEAGLEGL(t *testing.T) {
	const vcfIn = `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred GL">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	100	.	A	G	30	PASS	.	GT:PL	0/0:0,10,100	0/1:10,0,10
chr1	200	.	A	AT	30	PASS	.	GT:PL	0/1:0,5,50	0/1:0,5,50
chr1	300	.	A	G,C	30	PASS	.	GT:PL	0/1:0,5,50,5,10,50	0/1:0,5,50,5,10,50
`
	tmp := t.TempDir()
	prefix := filepath.Join(tmp, "x")
	err := Run(strings.NewReader(vcfIn), &Params{
		OutPrefix: prefix,
		BEAGLEGL:  true,
		BEAGLEPL:  true,
	})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	for _, suf := range []string{".BEAGLE.GL", ".BEAGLE.PL"} {
		body, err := os.ReadFile(prefix + suf)
		if err != nil {
			t.Fatalf("read %s: %v", suf, err)
		}
		lines := strings.Split(strings.TrimRight(string(body), "\n"), "\n")
		if len(lines) != 2 {
			t.Errorf("%s: expected header + 1 SNP (chr1:100), got %d lines:\n%s", suf, len(lines), body)
		}
		if !strings.Contains(lines[1], "chr1:100") {
			t.Errorf("%s: expected chr1:100 row, got %s", suf, lines[1])
		}
	}
}

func strconvFloat(s string) (float64, error) {
	return strconv.ParseFloat(s, 64)
}
