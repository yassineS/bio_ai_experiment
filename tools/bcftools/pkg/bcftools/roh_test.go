package bcftools

import (
	"bytes"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// rohFixture: a single-sample VCF with a clear autozygous stretch in
// the middle. Layout:
//
//	chr1:10  HET    (HW signal)
//	chr1:20  HET    (HW signal)
//	chr1:1000 HOM-ALT  (AZ signal 1)
//	chr1:2000 HOM-ALT  (AZ signal 2)
//	chr1:3000 HOM-ALT  (AZ signal 3)
//	chr1:4000 HOM-ALT  (AZ signal 4)
//	chr1:5000 HOM-ALT  (AZ signal 5)
//	chr1:6000 HOM-REF  (AZ signal 6)
//	chr1:7000 HOM-REF  (AZ signal 7)
//	chr1:8000 HET    (back to HW)
//	chr1:9000 HET    (HW)
//
// AF tags vary so the emission probabilities have signal.
func rohFixture() string {
	return `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=100000>
##INFO=<ID=AF,Number=1,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMP
chr1	10	.	A	T	.	PASS	AF=0.5	GT	0/1
chr1	20	.	G	C	.	PASS	AF=0.5	GT	0/1
chr1	1000	.	C	A	.	PASS	AF=0.3	GT	1/1
chr1	2000	.	T	G	.	PASS	AF=0.3	GT	1/1
chr1	3000	.	A	C	.	PASS	AF=0.3	GT	1/1
chr1	4000	.	G	T	.	PASS	AF=0.3	GT	1/1
chr1	5000	.	C	G	.	PASS	AF=0.3	GT	1/1
chr1	6000	.	T	A	.	PASS	AF=0.3	GT	0/0
chr1	7000	.	A	T	.	PASS	AF=0.3	GT	0/0
chr1	8000	.	G	C	.	PASS	AF=0.5	GT	0/1
chr1	9000	.	C	A	.	PASS	AF=0.5	GT	0/1
`
}

func TestParseRohOutputMode(t *testing.T) {
	cases := []struct {
		in   string
		want RohOutputMode
		err  bool
	}{
		{"", RohOutputRegions, false},
		{"r", RohOutputRegions, false},
		{"R", RohOutputRegions, false},
		{"s", RohOutputSites, false},
		{"sr", RohOutputBoth, false},
		{"rs", RohOutputBoth, false},
		{"bogus", 0, true},
	}
	for _, c := range cases {
		got, err := ParseRohOutputMode(c.in)
		if c.err {
			if err == nil {
				t.Errorf("ParseRohOutputMode(%q) expected error", c.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseRohOutputMode(%q): %v", c.in, err)
			continue
		}
		if got != c.want {
			t.Errorf("ParseRohOutputMode(%q) = %v want %v", c.in, got, c.want)
		}
	}
}

func TestClassifyHardGT(t *testing.T) {
	cases := []struct {
		in     string
		want   int
		wantOK bool
	}{
		{"0/0", 0, true},
		{"0|0", 0, true},
		{"0/1", 1, true},
		{"1|0", 1, true},
		{"1/1", 2, true},
		{"2/2", 2, true}, // any homozygous ALT (>0) counts as HOM-ALT
		{"./.", 0, false},
	}
	for _, c := range cases {
		got, ok := classifyHardGT(c.in)
		if ok != c.wantOK {
			t.Errorf("classifyHardGT(%q) ok = %v want %v", c.in, ok, c.wantOK)
		}
		if ok && got != c.want {
			t.Errorf("classifyHardGT(%q) = %d want %d", c.in, got, c.want)
		}
	}
}

func TestReadAFTag(t *testing.T) {
	mk := func(af string) *vcf.Variant {
		v := &vcf.Variant{Info: map[string]string{}}
		if af != "" {
			v.Info["AF"] = af
		}
		return v
	}
	if got := readAFTag(mk("0.25"), "AF", 0.4); math.Abs(got-0.25) > 1e-9 {
		t.Errorf("readAFTag base = %f want 0.25", got)
	}
	// Comma-list: first value wins.
	if got := readAFTag(mk("0.1,0.9"), "AF", 0.4); math.Abs(got-0.1) > 1e-9 {
		t.Errorf("readAFTag list = %f want 0.1", got)
	}
	// Missing tag -> default.
	v3 := mk("")
	v3.Info = nil
	if got := readAFTag(v3, "AF", 0.7); math.Abs(got-0.7) > 1e-9 {
		t.Errorf("readAFTag missing tag = %f want 0.7", got)
	}
	// 0 and 1 are clamped.
	if got := readAFTag(mk("0"), "AF", 0.4); got <= 0 {
		t.Errorf("readAFTag 0 should be clamped > 0, got %f", got)
	}
	if got := readAFTag(mk("1"), "AF", 0.4); got >= 1 {
		t.Errorf("readAFTag 1 should be clamped < 1, got %f", got)
	}
}

func TestRoh_DetectsAZStretch(t *testing.T) {
	out := &bytes.Buffer{}
	res, err := Roh(strings.NewReader(rohFixture()), out, RohOptions{
		Output: RohOutputRegions,
	})
	if err != nil {
		t.Fatalf("Roh: %v", err)
	}
	// We expect at least one AZ run, and it should cover the HOM
	// stretch (chr1:1000..7000). Loose end check: start before
	// 2000, end after 5000.
	if len(res.Segments) == 0 {
		t.Fatalf("expected at least one AZ segment, got none.\noutput:\n%s", out.String())
	}
	found := false
	for _, s := range res.Segments {
		if s.Chrom != "chr1" {
			continue
		}
		if s.StartPos <= 2000 && s.EndPos >= 5000 && s.NMarkers >= 4 {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("no segment covered the HOM stretch.\nsegments=%+v\noutput:\n%s",
			res.Segments, out.String())
	}

	// The output must contain the RG header banner and at least one
	// "RG\tSAMP\tchr1" data row.
	body := out.String()
	for _, want := range []string{
		"# RG\t",
		"RG\tSAMP\tchr1\t",
	} {
		if !strings.Contains(body, want) {
			t.Errorf("output missing %q\nbody=\n%s", want, body)
		}
	}
}

func TestRoh_SitesModeEmitsST(t *testing.T) {
	out := &bytes.Buffer{}
	_, err := Roh(strings.NewReader(rohFixture()), out, RohOptions{
		Output: RohOutputSites,
	})
	if err != nil {
		t.Fatalf("Roh: %v", err)
	}
	body := out.String()
	if !strings.Contains(body, "# ST\t") {
		t.Errorf("missing ST header banner")
	}
	// Every per-site row begins with "ST\tSAMP\tchr1\t".
	stLines := 0
	for _, ln := range strings.Split(body, "\n") {
		if strings.HasPrefix(ln, "ST\tSAMP\tchr1\t") {
			stLines++
		}
	}
	if stLines != 11 {
		t.Errorf("expected 11 ST rows (one per called site), got %d.\nbody=\n%s", stLines, body)
	}
}

func TestRoh_BothModeEmitsBoth(t *testing.T) {
	out := &bytes.Buffer{}
	_, err := Roh(strings.NewReader(rohFixture()), out, RohOptions{
		Output: RohOutputBoth,
	})
	if err != nil {
		t.Fatalf("Roh: %v", err)
	}
	body := out.String()
	hasST := strings.Contains(body, "ST\tSAMP\tchr1\t")
	hasRG := strings.Contains(body, "RG\tSAMP\tchr1\t")
	if !hasST || !hasRG {
		t.Errorf("both mode: ST=%v RG=%v\nbody=\n%s", hasST, hasRG, body)
	}
}

func TestRoh_UnknownSampleErrors(t *testing.T) {
	_, err := Roh(strings.NewReader(rohFixture()), &bytes.Buffer{}, RohOptions{
		Samples: []string{"NOPE"},
	})
	if err == nil {
		t.Fatal("expected error on unknown sample")
	}
}

func TestRoh_FileRoundTrip(t *testing.T) {
	tmp := t.TempDir()
	in := filepath.Join(tmp, "in.vcf")
	if err := os.WriteFile(in, []byte(rohFixture()), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	out := &bytes.Buffer{}
	res, err := RohFile(in, out, RohOptions{})
	if err != nil {
		t.Fatalf("RohFile: %v", err)
	}
	if res.NSites != 11 {
		t.Errorf("expected 11 sites total, got %d", res.NSites)
	}
}

func TestRoh_NoSegmentsWhenAllHet(t *testing.T) {
	src := `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=100>
##INFO=<ID=AF,Number=1,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
chr1	10	.	A	T	.	PASS	AF=0.5	GT	0/1
chr1	20	.	G	C	.	PASS	AF=0.5	GT	0/1
chr1	30	.	C	A	.	PASS	AF=0.5	GT	0/1
`
	res, err := Roh(strings.NewReader(src), &bytes.Buffer{}, RohOptions{
		Output: RohOutputRegions,
	})
	if err != nil {
		t.Fatalf("Roh: %v", err)
	}
	if len(res.Segments) != 0 {
		t.Errorf("expected no AZ segments for all-het input, got %+v", res.Segments)
	}
}

func TestRoh_EmissionLogAZHomCorrelatesWithAF(t *testing.T) {
	// Invariant: for the AZ state, P(observe HOM-ALT) should grow
	// monotonically with the AF (homozygotes are more probable when
	// the allele is more common).
	prev := math.Inf(-1)
	for _, af := range []float64{0.05, 0.1, 0.3, 0.5, 0.7, 0.9} {
		v := emissionLog(2, af, 1, 1e-3)
		if v <= prev {
			t.Errorf("emissionLog(HOM-ALT|AZ) not monotonic at af=%f: %f <= %f", af, v, prev)
		}
		prev = v
	}
}
