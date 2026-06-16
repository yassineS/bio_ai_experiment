package bcftools

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Binary-free table-driven tests for the PURE helpers added in this wave. They
// run with NO upstream binary and NO populated submodules (no exec.Command, no
// BCFTOOLS_PLUGINS), mirroring the established binary-free layer.

// TestUnitGuessPloidyGenomeRegion checks the -g/--genome preset table against
// the exact non-PAR chrX spans upstream hardcodes (guess-ploidy.c case 'g').
func TestUnitGuessPloidyGenomeRegion(t *testing.T) {
	cases := []struct {
		in     string
		want   string
		wantOK bool
	}{
		{"b37", "X:2699521-154931043", true},
		{"B37", "X:2699521-154931043", true}, // case-insensitive (strcasecmp)
		{"b38", "X:2781480-155701381", true},
		{"hg19", "chrX:2699521-154931043", true},
		{"hg38", "chrX:2781480-155701381", true},
		{"HG38", "chrX:2781480-155701381", true},
		{"grch37", "", false},
		{"", "", false},
	}
	for _, tc := range cases {
		got, ok := guessPloidyGenomeRegion(tc.in)
		if ok != tc.wantOK || got != tc.want {
			t.Errorf("guessPloidyGenomeRegion(%q) = (%q,%v), want (%q,%v)", tc.in, got, ok, tc.want, tc.wantOK)
		}
	}
}

// TestUnitGuessPloidyRewriteArgs checks -g PRESET expands to -r REGION in place,
// leaving other args untouched and rejecting an unknown preset.
func TestUnitGuessPloidyRewriteArgs(t *testing.T) {
	p := &guessPloidyPlugin{}
	got, err := p.RewriteArgs([]string{"-t", "GT", "-g", "b37", "-v"})
	if err != nil {
		t.Fatalf("RewriteArgs: %v", err)
	}
	want := []string{"-t", "GT", "-r", "X:2699521-154931043", "-v"}
	if len(got) != len(want) {
		t.Fatalf("RewriteArgs = %v, want %v", got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("RewriteArgs = %v, want %v", got, want)
		}
	}
	if _, err := p.RewriteArgs([]string{"-g", "nope"}); err == nil {
		t.Fatalf("RewriteArgs expected error for unknown preset")
	}
	// Attached forms -gPRESET and --genome=PRESET.
	for _, in := range [][]string{{"-gb38"}, {"--genome=b38"}} {
		got, err := p.RewriteArgs(in)
		if err != nil || len(got) != 2 || got[0] != "-r" || got[1] != "X:2781480-155701381" {
			t.Fatalf("RewriteArgs(%v) = %v, %v", in, got, err)
		}
	}
}

// TestUnitReadBinFile checks the bin-boundary file parser (one verbatim line per
// boundary, blank lines skipped, gzip-transparent via iohelper).
func TestUnitReadBinFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "bins.txt")
	if err := os.WriteFile(path, []byte("0\n0.25\n\n0.5\n1\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	got, err := readBinFile(path)
	if err != nil {
		t.Fatalf("readBinFile: %v", err)
	}
	want := []string{"0", "0.25", "0.5", "1"}
	if len(got) != len(want) {
		t.Fatalf("readBinFile = %v, want %v", got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("readBinFile = %v, want %v", got, want)
		}
	}
}

// TestUnitBinInitFile checks binInit reads from a file (no comma) and applies
// the same min/max boundary fix-up as the inline form.
func TestUnitBinInitFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "b.txt")
	if err := os.WriteFile(path, []byte("0.25\n0.5\n0.75\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	b, err := binInit(path, 0, 1)
	if err != nil {
		t.Fatalf("binInit(file): %v", err)
	}
	// The [0,1] fix-up prepends 0 and appends 1 around 0.25,0.5,0.75.
	want := []float32{0, 0.25, 0.5, 0.75, 1}
	if b.size() != len(want) {
		t.Fatalf("binInit size = %d, want %d (%v)", b.size(), len(want), b.vals)
	}
	for i, w := range want {
		if b.value(i) != w {
			t.Fatalf("bin[%d] = %v, want %v", i, b.value(i), w)
		}
	}
	// Equivalence with the inline comma form.
	inline, err := binInit("0.25,0.5,0.75", 0, 1)
	if err != nil {
		t.Fatalf("binInit(inline): %v", err)
	}
	if inline.size() != b.size() {
		t.Fatalf("file vs inline bin size mismatch: %d vs %d", b.size(), inline.size())
	}
}

// TestUnitTextListLine checks the remove-overlaps -Ot line formatter.
func TestUnitTextListLine(t *testing.T) {
	if got := textListLine("chr1", 200); got != "chr1\t200\n" {
		t.Errorf("textListLine = %q", got)
	}
}

// TestUnitMinQualValues checks the --missing value assignment: scalar default
// for missing QUAL, and the DP coverage scaling maxQual*DP/maxQualDP.
func TestUnitMinQualValues(t *testing.T) {
	mk := func(q float64, dp string) *vcf.Variant {
		return &vcf.Variant{Qual: q, Info: map[string]string{"DP": dp}, InfoOrder: []string{"DP"}}
	}
	run := []*vcf.Variant{mk(50, "40"), mk(-1, "80"), mk(-1, "5"), mk(30, "60")}

	// Scalar mode: present QUAL kept, missing -> the scalar (0 here).
	got := minQualValues(run, markMissingScalar, 0)
	want := []float64{50, 0, 0, 30}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("scalar minQualValues = %v, want %v", got, want)
		}
	}

	// DP mode: maxQual=50 at DP=40; missing -> 50*DP/40.
	gotDP := minQualValues(run, markMissingMaxDP, 0)
	wantDP := []float64{50, float64(float32(50) * 80 / 40), float64(float32(50) * 5 / 40), 30}
	for i := range wantDP {
		if gotDP[i] != wantDP[i] {
			t.Fatalf("DP minQualValues[%d] = %v, want %v", i, gotDP[i], wantDP[i])
		}
	}

	// DP mode with no present QUAL anywhere: maxQualDP==0, missing stays scalar.
	noQual := []*vcf.Variant{mk(-1, "10"), mk(-1, "20")}
	gotNo := minQualValues(noQual, markMissingMaxDP, 7)
	if gotNo[0] != 7 || gotNo[1] != 7 {
		t.Fatalf("DP minQualValues (no present QUAL) = %v, want [7 7]", gotNo)
	}
}

// TestUnitLAARemap checks the localized-allele remap helpers (LAD->AD Number=R,
// LPL->PL Number=G) against the process_LXX layout, including default cells and
// missing LAA.
func TestUnitLAARemap(t *testing.T) {
	// LAD->AD: nals=4, LAA=[2], LAD=[30,18] -> dst[0]=30, dst[2]=18, rest "."
	if got := ladToAD([]string{"30", "18"}, []int{2}, 4, "."); got != "30,.,18,." {
		t.Errorf("ladToAD = %q, want 30,.,18,.", got)
	}
	// LAA=[1,3], LAD=[0,20,15] -> 0,20,.,15
	if got := ladToAD([]string{"0", "20", "15"}, []int{1, 3}, 4, "."); got != "0,20,.,15" {
		t.Errorf("ladToAD = %q, want 0,20,.,15", got)
	}
	// Missing LAA -> only REF depth placed, rest default.
	if got := ladToAD([]string{"40"}, nil, 2, "."); got != "40,." {
		t.Errorf("ladToAD(missing LAA) = %q, want 40,.", got)
	}
	// Custom default 0.
	if got := ladToAD([]string{"30", "18"}, []int{2}, 4, "0"); got != "30,0,18,0" {
		t.Errorf("ladToAD(dflt 0) = %q, want 30,0,18,0", got)
	}

	// LPL->PL: nals=4, LAA=[2], LPL=[0,40,80] -> tmp_laa=[0,2]; G=10 cells.
	if got := lplToPL([]string{"0", "40", "80"}, []int{2}, 4, "."); got != "0,.,.,40,.,80,.,.,.,." {
		t.Errorf("lplToPL = %q", got)
	}
	// LAA=[1,3], LPL=[50,0,60,90,30,255] -> 50,0,60,.,.,.,90,30,.,255
	if got := lplToPL([]string{"50", "0", "60", "90", "30", "255"}, []int{1, 3}, 4, "."); got != "50,0,60,.,.,.,90,30,.,255" {
		t.Errorf("lplToPL = %q", got)
	}
	// Missing LAA -> only the REF/REF cell set.
	if got := lplToPL([]string{"0"}, nil, 2, "."); got != "0,.,." {
		t.Errorf("lplToPL(missing LAA) = %q, want 0,.,.", got)
	}
}

// TestUnitParseLAAList checks LAA parsing (missing field -> empty; stop at a
// missing element).
func TestUnitParseLAAList(t *testing.T) {
	cases := []struct {
		in   string
		want []int
	}{
		{"2", []int{2}},
		{"1,3", []int{1, 3}},
		{".", nil},
		{"", nil},
		{"1,.", []int{1}},
	}
	for _, tc := range cases {
		got := parseLAAList(tc.in)
		if len(got) != len(tc.want) {
			t.Fatalf("parseLAAList(%q) = %v, want %v", tc.in, got, tc.want)
		}
		for i := range tc.want {
			if got[i] != tc.want[i] {
				t.Fatalf("parseLAAList(%q) = %v, want %v", tc.in, got, tc.want)
			}
		}
	}
}

// TestUnitSubsetNumberedList checks the bcf_remove_allele_set value-subsetting
// for Number=A/R/G when ALT allele 2 (of REF,A1,A2,A3) is removed.
func TestUnitSubsetNumberedList(t *testing.T) {
	nROri := 4
	rm := []bool{false, false, true, false} // drop ALT index 2
	mapIdx := []int{0, 1, -1, 2}
	// Number=A: one per ALT (A1,A2,A3) -> drop A2.
	if got, ch := subsetNumberedList("11,22,33", "A", rm, mapIdx, nROri); !ch || got != "11,33" {
		t.Errorf("A subset = %q,%v", got, ch)
	}
	// Number=R: REF + ALTs -> drop A2 entry.
	if got, ch := subsetNumberedList("0,11,22,33", "R", rm, mapIdx, nROri); !ch || got != "0,11,33" {
		t.Errorf("R subset = %q,%v", got, ch)
	}
	// Number=G diploid: 10 cells (b*(b+1)/2+a). Drop any cell touching allele 2.
	// Genotype order: 00,01,11,02,12,22,03,13,23,33 -> keep cells not touching 2.
	gIn := "0,1,2,3,4,5,6,7,8,9"
	gWant := "0,1,2,6,7,9" // 00,01,11,03,13,33
	if got, ch := subsetNumberedList(gIn, "G", rm, mapIdx, nROri); !ch || got != gWant {
		t.Errorf("G subset = %q,%v, want %q", got, ch, gWant)
	}
	// A missing "." value is left unchanged.
	if got, ch := subsetNumberedList(".", "A", rm, mapIdx, nROri); ch || got != "." {
		t.Errorf("A missing = %q,%v", got, ch)
	}
	// A count mismatch leaves it unchanged.
	if got, ch := subsetNumberedList("1,2", "A", rm, mapIdx, nROri); ch || got != "1,2" {
		t.Errorf("A mismatch = %q,%v", got, ch)
	}
}

// TestUnitRemoveAlleleSet checks the end-to-end allele subset on a record: ALT,
// INFO Number=A, FORMAT/AD (R), FORMAT/PL (G) and GT reindexing.
func TestUnitRemoveAlleleSet(t *testing.T) {
	hdr := &vcf.Header{
		Samples: []string{"S1"},
		MetaInfo: []string{
			`##INFO=<ID=AC,Number=A,Type=Integer,Description="x">`,
			`##INFO=<ID=AN,Number=1,Type=Integer,Description="x">`,
			`##FORMAT=<ID=GT,Number=1,Type=String,Description="x">`,
			`##FORMAT=<ID=AD,Number=R,Type=Integer,Description="x">`,
			`##FORMAT=<ID=PL,Number=G,Type=Integer,Description="x">`,
		},
	}
	v := &vcf.Variant{
		Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"T", "G", "C"},
		Info:      map[string]string{"AC": "1,2,1", "AN": "8"},
		InfoOrder: []string{"AC", "AN"},
		Format:    []string{"GT", "AD", "PL"},
		Samples: []vcf.Sample{{Data: map[string]string{
			"GT": "1/3",
			"AD": "0,20,18,15",                    // R: REF,T,G,C
			"PL": "10,11,12,13,14,15,16,17,18,19", // G: 10 cells for 4 alleles
		}}},
	}
	// Remove ALT G (allele index 2) and ALT C (index 3); keep T (index 1).
	rm := []bool{false, false, true, true}
	if err := removeAlleleSet(hdr, v, rm); err != nil {
		t.Fatalf("removeAlleleSet: %v", err)
	}
	if len(v.Alt) != 1 || v.Alt[0] != "T" {
		t.Fatalf("ALT = %v, want [T]", v.Alt)
	}
	if v.Info["AC"] != "1" { // only T's AC survives
		t.Errorf("AC = %q, want 1", v.Info["AC"])
	}
	if v.Info["AN"] != "8" { // Number=1 untouched
		t.Errorf("AN = %q, want 8", v.Info["AN"])
	}
	// AD: REF + T only.
	if v.Samples[0].Data["AD"] != "0,20" {
		t.Errorf("AD = %q, want 0,20", v.Samples[0].Data["AD"])
	}
	// PL: keep cells touching only alleles {0,1}: indices 0(00),1(01),2(11).
	if v.Samples[0].Data["PL"] != "10,11,12" {
		t.Errorf("PL = %q, want 10,11,12", v.Samples[0].Data["PL"])
	}
	// GT 1/3 -> allele 1 stays index 1, allele 3 removed -> "." : "1/."
	if v.Samples[0].Data["GT"] != "1/." {
		t.Errorf("GT = %q, want 1/.", v.Samples[0].Data["GT"])
	}
}
