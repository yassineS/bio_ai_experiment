package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

func newSyntheticVariant(alts, fmts []string) *vcf.Variant {
	return &vcf.Variant{
		Chrom:  "chr1",
		Pos:    1,
		Ref:    "A",
		Alt:    alts,
		Format: fmts,
	}
}

// trio2Fixture: three samples in two families with two trios, plus
// an extra (CHILD2/MISSING1/MOM_X). MOTHER_X is in the VCF as
// MOTHER_X. The variants exercise every per-record code path:
//
//   - chr1:10 — both trios consistent.
//   - chr1:20 — TRIO1 inconsistent (1/1 from 0/0 x 0/0), TRIO2 good.
//   - chr1:30 — TRIO1 missing (./.), TRIO2 good.
//   - chr1:40 — no ALT (skipped as ref-only).
//   - chr1:50 — no FORMAT/GT (skipped as no_gt).
//   - chr1:60 — both trios inconsistent.
func trio2Fixture() string {
	return `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CHILD	FATHER	MOTHER	CHILD2	FATHER2	MOTHER2
chr1	10	.	A	T	.	PASS	DP=1	GT	0/1	0/0	0/1	0/0	0/0	0/0
chr1	20	.	G	C	.	PASS	DP=2	GT	1/1	0/0	0/0	0/1	0/1	0/0
chr1	30	.	C	A	.	PASS	DP=3	GT	./.	0/1	0/1	0/0	0/0	0/0
chr1	40	.	A	.	.	PASS	DP=4	GT	0/0	0/0	0/0	0/0	0/0	0/0
chr1	50	.	A	T	.	PASS	DP=5	DP	1	2	3	4	5	6
chr1	60	.	G	A	.	PASS	DP=6	GT	1/1	0/0	0/0	1/1	0/0	0/0
`
}

func TestParseMendelian2Mode(t *testing.T) {
	cases := []struct {
		in   string
		want Mendelian2Mode
	}{
		{"", Mendelian2Count},
		{"c", Mendelian2Count},
		{"a", Mendelian2Annotate},
		{"d", Mendelian2DeleteGT},
		{"e", Mendelian2ListErr},
		{"E", Mendelian2DropErr},
		{"g", Mendelian2ListGood},
		{"m", Mendelian2ListMiss},
		{"M", Mendelian2DropMiss},
		{"S", Mendelian2DropSkip},
		{"s", Mendelian2DropSkip},
		{"x", Mendelian2ListErr},                        // legacy alias for 'e'
		{"u", Mendelian2ListMiss},                       // legacy alias for 'm'
		{"+", Mendelian2ListGood},                       // legacy alias for 'g'
		{"ad", Mendelian2Annotate | Mendelian2DeleteGT}, // combined
		{"aE", Mendelian2Annotate | Mendelian2DropErr},  // combined
		{"aeg", Mendelian2Annotate | Mendelian2ListErr | Mendelian2ListGood},
	}
	for _, c := range cases {
		got, err := ParseMendelian2Mode(c.in)
		if err != nil {
			t.Errorf("ParseMendelian2Mode(%q): unexpected error %v", c.in, err)
			continue
		}
		if got != c.want {
			t.Errorf("ParseMendelian2Mode(%q) = %d (0x%x) want %d (0x%x)", c.in, got, got, c.want, c.want)
		}
	}
	if _, err := ParseMendelian2Mode("Z"); err == nil {
		t.Error("expected error on unknown mode letter")
	}
}

func TestParseMendelian2PFM(t *testing.T) {
	cases := []struct {
		in      string
		want    Mendelian2PFM
		wantErr bool
	}{
		{"kid,dad,mom", Mendelian2PFM{Child: "kid", Father: "dad", Mother: "mom", Sex: 0}, false},
		{"1X:kid,dad,mom", Mendelian2PFM{Child: "kid", Father: "dad", Mother: "mom", Sex: 1}, false},
		{"2X:kid,dad,mom", Mendelian2PFM{Child: "kid", Father: "dad", Mother: "mom", Sex: 2}, false},
		{"1x:kid,dad,mom", Mendelian2PFM{Child: "kid", Father: "dad", Mother: "mom", Sex: 1}, false}, // case-insensitive
		{"", Mendelian2PFM{}, true},
		{"kid,dad", Mendelian2PFM{}, true},
		{"kid,,mom", Mendelian2PFM{}, true},
		{"1X:kid,,mom", Mendelian2PFM{}, true},
	}
	for _, c := range cases {
		got, err := ParseMendelian2PFM(c.in)
		if c.wantErr {
			if err == nil {
				t.Errorf("ParseMendelian2PFM(%q) expected error, got %+v", c.in, got)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseMendelian2PFM(%q) unexpected error: %v", c.in, err)
			continue
		}
		if got != c.want {
			t.Errorf("ParseMendelian2PFM(%q) = %+v want %+v", c.in, got, c.want)
		}
	}
}

func TestMendelian2_CountModeDefault(t *testing.T) {
	var out bytes.Buffer
	sum, err := Mendelian2(strings.NewReader(trio2Fixture()), &out, Mendelian2Options{
		Trios: []Mendelian2Trio{
			{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"},
			{Child: "CHILD2", Father: "FATHER2", Mother: "MOTHER2"},
		},
		// default Mode -> Count
	})
	if err != nil {
		t.Fatalf("Mendelian2: %v", err)
	}
	if sum.TotalRecords != 6 {
		t.Errorf("want 6 records, got %d", sum.TotalRecords)
	}
	if sum.SitesRefOnly != 1 {
		t.Errorf("want 1 ref-only, got %d", sum.SitesRefOnly)
	}
	if sum.SitesNoGT != 1 {
		t.Errorf("want 1 no-GT, got %d", sum.SitesNoGT)
	}
	// chr1:20 (TRIO1), chr1:60 (both) -> 2 sites with at least one error.
	if sum.SitesMERR != 2 {
		t.Errorf("want 2 MERR sites, got %d", sum.SitesMERR)
	}
	// chr1:30 (TRIO1 ./.) -> 1 site with at least one missing.
	if sum.SitesMissing != 1 {
		t.Errorf("want 1 missing site, got %d", sum.SitesMissing)
	}
	body := out.String()
	for _, want := range []string{
		"# Summary stats",
		"sites_ref_only\t1",
		"sites_no_GT\t1",
		"sites_merr\t2",
		"sites_missing\t1",
		"# TRIO\t[2]id\t[3]child\t[4]father\t[5]mother",
		"TRIO\t1\tCHILD\tFATHER\tMOTHER",
		"TRIO\t2\tCHILD2\tFATHER2\tMOTHER2",
	} {
		if !strings.Contains(body, want) {
			t.Errorf("summary missing %q:\n%s", want, body)
		}
	}
}

func TestMendelian2_AnnotateMode(t *testing.T) {
	var out bytes.Buffer
	_, err := Mendelian2(strings.NewReader(trio2Fixture()), &out, Mendelian2Options{
		Trios: []Mendelian2Trio{
			{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"},
			{Child: "CHILD2", Father: "FATHER2", Mother: "MOTHER2"},
		},
		Mode:         Mendelian2Annotate,
		OutputFormat: OutputVCF,
	})
	if err != nil {
		t.Fatalf("Mendelian2: %v", err)
	}
	body := out.String()
	for _, want := range []string{
		`##INFO=<ID=MERR,`,
		"chr1\t10\t.\tA\tT",
		";MERR=0",
		"chr1\t20",
		";MERR=1", // TRIO1 has an error
		"chr1\t60",
		";MERR=2", // both trios have errors
	} {
		if !strings.Contains(body, want) {
			t.Errorf("annotate output missing %q:\n%s", want, body)
		}
	}
}

func TestMendelian2_DropErrMode(t *testing.T) {
	var out bytes.Buffer
	_, err := Mendelian2(strings.NewReader(trio2Fixture()), &out, Mendelian2Options{
		Trios: []Mendelian2Trio{
			{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"},
			{Child: "CHILD2", Father: "FATHER2", Mother: "MOTHER2"},
		},
		Mode:         Mendelian2DropErr,
		OutputFormat: OutputVCF,
	})
	if err != nil {
		t.Fatalf("Mendelian2: %v", err)
	}
	body := out.String()
	// chr1:20 and chr1:60 must be dropped.
	if strings.Contains(body, "chr1\t20") {
		t.Errorf("chr1:20 should be dropped under -m E:\n%s", body)
	}
	if strings.Contains(body, "chr1\t60") {
		t.Errorf("chr1:60 should be dropped under -m E:\n%s", body)
	}
	if !strings.Contains(body, "chr1\t10") {
		t.Errorf("chr1:10 should be kept under -m E:\n%s", body)
	}
}

func TestMendelian2_ListErrMode(t *testing.T) {
	var out bytes.Buffer
	_, err := Mendelian2(strings.NewReader(trio2Fixture()), &out, Mendelian2Options{
		Trios: []Mendelian2Trio{
			{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"},
			{Child: "CHILD2", Father: "FATHER2", Mother: "MOTHER2"},
		},
		Mode:         Mendelian2ListErr,
		OutputFormat: OutputVCF,
	})
	if err != nil {
		t.Fatalf("Mendelian2: %v", err)
	}
	body := out.String()
	// Only chr1:20 and chr1:60 should appear.
	for _, want := range []string{"chr1\t20", "chr1\t60"} {
		if !strings.Contains(body, want) {
			t.Errorf("expected %q under -m e:\n%s", want, body)
		}
	}
	if strings.Contains(body, "chr1\t10") {
		t.Errorf("chr1:10 (no err) should not appear under -m e:\n%s", body)
	}
}

func TestMendelian2_DeleteGTMode(t *testing.T) {
	var out bytes.Buffer
	_, err := Mendelian2(strings.NewReader(trio2Fixture()), &out, Mendelian2Options{
		Trios: []Mendelian2Trio{
			{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"},
		},
		Mode:         Mendelian2DeleteGT,
		OutputFormat: OutputVCF,
	})
	if err != nil {
		t.Fatalf("Mendelian2: %v", err)
	}
	body := out.String()
	// chr1:20 record is present but TRIO1's GTs are "./.".
	var line20 string
	for _, l := range strings.Split(body, "\n") {
		if strings.HasPrefix(l, "chr1\t20\t") {
			line20 = l
			break
		}
	}
	if line20 == "" {
		t.Fatalf("chr1:20 record missing under -m d:\n%s", body)
	}
	// Expect the trio1 columns (CHILD, FATHER, MOTHER) to be ./.
	parts := strings.Split(line20, "\t")
	if len(parts) < 12 {
		t.Fatalf("chr1:20 line too short: %q", line20)
	}
	for _, p := range parts[9:12] {
		if p != "./." {
			t.Errorf("expected CHILD/FATHER/MOTHER GT == ./. at chr1:20, got %q (line: %s)", p, line20)
		}
	}
}

func TestMendelian2_PFMSingleTrio(t *testing.T) {
	var out bytes.Buffer
	pfm := Mendelian2PFM{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"}
	sum, err := Mendelian2(strings.NewReader(trio2Fixture()), &out, Mendelian2Options{
		PFM:  &pfm,
		Mode: Mendelian2Count,
	})
	if err != nil {
		t.Fatalf("Mendelian2: %v", err)
	}
	if len(sum.Trios) != 1 {
		t.Errorf("want 1 trio (PFM), got %d", len(sum.Trios))
	}
	if sum.Trios[0].NMErr != 2 {
		t.Errorf("trio1 errors: want 2, got %d", sum.Trios[0].NMErr)
	}
}

func TestMendelian2_PEDFile(t *testing.T) {
	dir := t.TempDir()
	inPath := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(inPath, []byte(trio2Fixture()), 0644); err != nil {
		t.Fatal(err)
	}
	pedPath := filepath.Join(dir, "trios.ped")
	ped := "# fam ind dad mom sex pheno\n" +
		"FAM1 CHILD FATHER MOTHER 1 -9\n" +
		"FAM2 CHILD2 FATHER2 MOTHER2 2 -9\n" +
		"FAM3 GHOST GHOSTDAD GHOSTMOM 1 -9\n" + // skipped: not in VCF
		"FAM4 CHILD GHOSTDAD MOTHER 1 -9\n" // skipped: missing parent
	if err := os.WriteFile(pedPath, []byte(ped), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	sum, err := Mendelian2File(inPath, &out, Mendelian2Options{
		PEDFile: pedPath,
		Mode:    Mendelian2Count,
	})
	if err != nil {
		t.Fatalf("Mendelian2File: %v", err)
	}
	if len(sum.Trios) != 2 {
		t.Errorf("want 2 trios from PED, got %d", len(sum.Trios))
	}
	// trios are sorted by child name -> CHILD then CHILD2.
	if sum.Trios[0].Child != "CHILD" || sum.Trios[1].Child != "CHILD2" {
		t.Errorf("trios in wrong order: %+v", sum.Trios)
	}
	if sum.Trios[0].Sex != 1 || sum.Trios[1].Sex != 2 {
		t.Errorf("trios sex mismatch: %+v / %+v", sum.Trios[0], sum.Trios[1])
	}
}

func TestMendelian2_MissingTrioSource(t *testing.T) {
	var out bytes.Buffer
	_, err := Mendelian2(strings.NewReader(trio2Fixture()), &out, Mendelian2Options{})
	if err == nil {
		t.Error("expected error when no trio source is given")
	}
}

func TestMendelian2_UnknownChildName(t *testing.T) {
	var out bytes.Buffer
	_, err := Mendelian2(strings.NewReader(trio2Fixture()), &out, Mendelian2Options{
		Trios: []Mendelian2Trio{{Child: "GHOST", Father: "FATHER", Mother: "MOTHER"}},
	})
	if err == nil {
		t.Error("expected error when child sample is not in input")
	}
}

func TestMendelian2_ClassifyRecord(t *testing.T) {
	// Bit of white-box: locked-in mapping of record shapes to skip
	// reasons. If we ever change this we want the test to break.
	cases := []struct {
		name string
		alts []string
		fmts []string
		want skipReason
	}{
		{"normal", []string{"T"}, []string{"GT"}, skipNone},
		{"no-alt", nil, []string{"GT"}, skipRefOnly},
		{"dot-alt", []string{"."}, []string{"GT"}, skipRefOnly},
		{"no-gt", []string{"T"}, []string{"DP"}, skipNoGT},
		{"many-als", make([]string, 65), []string{"GT"}, skipManyAls},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			v := &fakeVariant{alts: tc.alts, fmts: tc.fmts}
			// minimal stand-in: classifyRecord only reads Alt and
			// Format. We build a vcf.Variant directly.
			realV := newSyntheticVariant(tc.alts, tc.fmts)
			got := classifyRecord(realV)
			if got != tc.want {
				t.Errorf("classifyRecord(%s)=%v want %v (alts=%v fmts=%v)", tc.name, got, tc.want, v.alts, v.fmts)
			}
		})
	}
}

// fakeVariant is just a struct-printable witness for table tests; we
// pass a real *vcf.Variant to classifyRecord (built by
// newSyntheticVariant). Kept separate so the trace inside t.Errorf
// stays readable.
type fakeVariant struct {
	alts []string
	fmts []string
}
