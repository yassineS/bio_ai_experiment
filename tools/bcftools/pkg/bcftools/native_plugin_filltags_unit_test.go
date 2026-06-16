package bcftools

import (
	"math"
	"os"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Binary-free unit tests for the fill-tags pure helpers: the tag formula
// calculators (calcHWE, alleleTotal, formatExprValue, roundHalfAway), the
// tag-list classifier (parseTags), the -S samples-file parser
// (parseFillTagsSamples), and the expression-function evaluators (compileFillExpr
// / evaluate, reduce). These run with NO upstream binary and with the submodules
// unpopulated (no exec.Command, no BCFTOOLS_PLUGINS).

func exprHeader(samples []string, fmtTags ...string) *vcf.Header {
	h := &vcf.Header{Samples: samples}
	for _, t := range fmtTags {
		h.MetaInfo = append(h.MetaInfo, `##FORMAT=<ID=`+t+`,Number=1,Type=Integer,Description="x">`)
	}
	return h
}

func smplVariant(fmtTags []string, perSample []map[string]string) *vcf.Variant {
	v := &vcf.Variant{Format: fmtTags}
	for _, d := range perSample {
		v.Samples = append(v.Samples, vcf.Sample{Data: d})
	}
	return v
}

func TestUnitFillTagsReduce(t *testing.T) {
	nan := math.NaN()
	cases := []struct {
		name string
		op   aggOp
		in   []float64
		want float64
		ok   bool
	}{
		{"sum", aggSum, []float64{10, 20, 5, 15}, 50, true},
		{"sum-with-missing", aggSum, []float64{10, nan, 5}, 15, true},
		{"avg", aggAvg, []float64{10, 20}, 15, true},
		{"max", aggMax, []float64{3, 9, 1}, 9, true},
		{"min", aggMin, []float64{3, 9, 1}, 1, true},
		{"median-odd", aggMedian, []float64{3, 1, 2}, 2, true},
		{"median-even", aggMedian, []float64{1, 2, 3, 4}, 2.5, true},
		{"median-one", aggMedian, []float64{7}, 7, true},
		{"stdev", aggStdev, []float64{2, 4, 4, 4, 5, 5, 7, 9}, 2, true},
		{"stdev-one", aggStdev, []float64{5}, 0, true},
		{"all-missing", aggSum, []float64{nan, nan}, 0, false},
		{"empty", aggSum, nil, 0, false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got, ok := reduce(c.op, c.in)
			if ok != c.ok {
				t.Fatalf("ok=%v want %v", ok, c.ok)
			}
			if ok && math.Abs(got-c.want) > 1e-9 {
				t.Fatalf("got %v want %v", got, c.want)
			}
		})
	}
}

func TestUnitFillTagsExprSiteSum(t *testing.T) {
	hdr := exprHeader([]string{"A", "B", "C"}, "DP")
	e, err := compileFillExpr("sum(FORMAT/DP)", hdr)
	if err != nil {
		t.Fatal(err)
	}
	v := smplVariant([]string{"GT", "DP"}, []map[string]string{
		{"DP": "10"}, {"DP": "20"}, {"DP": "5"},
	})
	r := e.evaluate(v, nil)
	if r.perSample || len(r.site) != 1 || r.site[0] != 35 {
		t.Fatalf("site sum = %#v", r)
	}
}

func TestUnitFillTagsExprSmplSum(t *testing.T) {
	hdr := exprHeader([]string{"A", "B"}, "AD")
	e, err := compileFillExpr("smpl_sum(FORMAT/AD)", hdr)
	if err != nil {
		t.Fatal(err)
	}
	v := smplVariant([]string{"GT", "AD"}, []map[string]string{
		{"AD": "9,1"}, {"AD": "3,5,2"},
	})
	r := e.evaluate(v, nil)
	if !r.perSample {
		t.Fatalf("expected per-sample result: %#v", r)
	}
	if r.values[0][0] != 10 || r.values[1][0] != 10 {
		t.Fatalf("smpl_sum per-sample = %v", r.values)
	}
}

func TestUnitFillTagsExprPopMask(t *testing.T) {
	hdr := exprHeader([]string{"A", "B", "C"}, "DP")
	e, err := compileFillExpr("sum(FORMAT/DP)", hdr)
	if err != nil {
		t.Fatal(err)
	}
	v := smplVariant([]string{"GT", "DP"}, []map[string]string{
		{"DP": "10"}, {"DP": "20"}, {"DP": "5"},
	})
	// Restrict to samples A and C only.
	mask := []bool{true, false, true}
	r := e.evaluate(v, mask)
	if r.site[0] != 15 {
		t.Fatalf("masked sum = %v want 15", r.site[0])
	}
}

func TestUnitFillTagsExprArith(t *testing.T) {
	hdr := exprHeader([]string{"A", "B"}, "DP")
	e, err := compileFillExpr("sum(FORMAT/DP)/2", hdr)
	if err != nil {
		t.Fatal(err)
	}
	v := smplVariant([]string{"GT", "DP"}, []map[string]string{{"DP": "10"}, {"DP": "30"}})
	r := e.evaluate(v, nil)
	if r.site[0] != 20 {
		t.Fatalf("sum/2 = %v want 20", r.site[0])
	}
}

func TestUnitFillTagsExprNMissing(t *testing.T) {
	hdr := exprHeader([]string{"A", "B", "C", "D"})
	e, err := compileFillExpr("F_MISSING", hdr)
	if err != nil {
		t.Fatal(err)
	}
	v := smplVariant([]string{"GT"}, []map[string]string{
		{"GT": "0/0"}, {"GT": "./."}, {"GT": "0/1"}, {"GT": "./1"},
	})
	r := e.evaluate(v, nil)
	// Two of four samples have a missing allele.
	if math.Abs(r.site[0]-0.5) > 1e-9 {
		t.Fatalf("F_MISSING = %v want 0.5", r.site[0])
	}
}

func TestUnitFillTagsFormatExprValue(t *testing.T) {
	cases := []struct {
		x     float64
		isInt bool
		want  string
	}{
		{12.5, false, "12.5"},
		{12.5, true, "13"}, // round half away from zero
		{-12.5, true, "-13"},
		{50, false, "50"},
		{math.NaN(), false, "."},
		{math.NaN(), true, "."},
	}
	for _, c := range cases {
		got := formatExprValue(c.x, c.isInt)
		if got != c.want {
			t.Fatalf("formatExprValue(%v,%v)=%q want %q", c.x, c.isInt, got, c.want)
		}
	}
}

func TestUnitFillTagsCalcHWE(t *testing.T) {
	// A balanced biallelic site (many hets) should give HWE close to 1.
	ph, pe := calcHWE(10, 10, 6)
	if ph <= 0 || ph > 1 || pe <= 0 || pe > 1 {
		t.Fatalf("HWE out of range: ph=%v pe=%v", ph, pe)
	}
	// All-het extreme should yield low excess-het p-value (heterozygote excess).
	_, peExc := calcHWE(8, 8, 8)
	if peExc <= 0 || peExc > 1 {
		t.Fatalf("ExcHet out of range: %v", peExc)
	}
}

func TestUnitFillTagsParseTags(t *testing.T) {
	hdr := exprHeader(nil)
	out := &vcf.Header{}
	p := &fillTagsPlugin{}
	p.pops = []population{{name: "", suffix: ""}}
	flag, err := p.parseTags("AN,AC,AF", hdr, out)
	if err != nil {
		t.Fatal(err)
	}
	if flag != (setAN | setAC | setAF) {
		t.Fatalf("flag=%b want AN|AC|AF", flag)
	}

	// "all" excludes END/TYPE.
	out2 := &vcf.Header{}
	p2 := &fillTagsPlugin{pops: []population{{}}}
	allFlag, err := p2.parseTags("all", hdr, out2)
	if err != nil {
		t.Fatal(err)
	}
	if allFlag&setEND != 0 || allFlag&setType != 0 {
		t.Fatalf("'all' should not set END/TYPE: %b", allFlag)
	}
	if allFlag&setVAF == 0 || allFlag&setHWE == 0 {
		t.Fatalf("'all' should set VAF and HWE: %b", allFlag)
	}

	// Unknown tag errors.
	p3 := &fillTagsPlugin{pops: []population{{}}}
	if _, err := p3.parseTags("BOGUS", hdr, &vcf.Header{}); err == nil {
		t.Fatalf("expected error for unknown tag")
	}
}

func TestUnitFillTagsParseSamples(t *testing.T) {
	dir := t.TempDir()
	fname := filepath.Join(dir, "groups.txt")
	content := "SAMPLE1 GRP1\nSAMPLE2 GRP1,GRP2\nSAMPLE3 GRP2\nMISSINGSAMP GRP3\n"
	if err := os.WriteFile(fname, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	hdr := &vcf.Header{Samples: []string{"SAMPLE1", "SAMPLE2", "SAMPLE3"}}
	pops, err := parseFillTagsSamples(fname, hdr, nil)
	if err != nil {
		t.Fatal(err)
	}
	// GRP1 first (SAMPLE1, SAMPLE2), GRP2 next (SAMPLE2, SAMPLE3); GRP3 skipped
	// because its only member is absent from the VCF.
	if len(pops) != 2 {
		t.Fatalf("got %d pops want 2: %+v", len(pops), pops)
	}
	if pops[0].name != "GRP1" || pops[0].suffix != "_GRP1" {
		t.Fatalf("pop0 = %+v", pops[0])
	}
	wantG1 := []bool{true, true, false}
	if !reflect.DeepEqual(pops[0].mask, wantG1) {
		t.Fatalf("GRP1 mask = %v want %v", pops[0].mask, wantG1)
	}
	wantG2 := []bool{false, true, true}
	if !reflect.DeepEqual(pops[1].mask, wantG2) {
		t.Fatalf("GRP2 mask = %v want %v", pops[1].mask, wantG2)
	}
}

// TestUnitFillTagsParseSamplesShortNames pins the fix-on-port: upstream's
// parse_samples() rejects sample names of length <= 2 (off-by-one in its
// reverse scan, see docs/UPSTREAM_BUGS.md#bcftools-fill-tags-samples-short-names);
// our field-split parser accepts them.
func TestUnitFillTagsParseSamplesShortNames(t *testing.T) {
	dir := t.TempDir()
	fname := filepath.Join(dir, "short.txt")
	if err := os.WriteFile(fname, []byte("S1 GRP1\nS2 GRP1,GRP2\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	hdr := &vcf.Header{Samples: []string{"S1", "S2"}}
	pops, err := parseFillTagsSamples(fname, hdr, nil)
	if err != nil {
		t.Fatalf("short sample names should parse: %v", err)
	}
	if len(pops) != 2 {
		t.Fatalf("got %d pops want 2 (GRP1, GRP2)", len(pops))
	}
}

func TestUnitFillTagsAddExprHeader(t *testing.T) {
	hdr := exprHeader([]string{"A", "B"}, "DP")
	out := &vcf.Header{Samples: hdr.Samples}
	p := &fillTagsPlugin{pops: []population{{name: "", suffix: ""}}}
	if err := p.addExpr("DP2:1=int(sum(FORMAT/DP))", "int(sum(FORMAT/DP))", hdr, out); err != nil {
		t.Fatal(err)
	}
	if len(p.bindings) != 1 {
		t.Fatalf("expected one binding, got %d", len(p.bindings))
	}
	b := p.bindings[0]
	if b.dstTag != "DP2" || !b.isInt || !b.fixedLen || b.count != 1 || b.isFormat {
		t.Fatalf("binding = %+v", b)
	}
	wantHdr := `##INFO=<ID=DP2,Number=1,Type=Integer,Description="Added by +fill-tags expression DP2:1=int(sum(FORMAT/DP))">`
	found := false
	for _, m := range out.MetaInfo {
		if m == wantHdr {
			found = true
		}
	}
	if !found {
		t.Fatalf("header line not added; got %v", out.MetaInfo)
	}
}

func TestUnitFillTagsRoundHalfAway(t *testing.T) {
	cases := []struct {
		in   float64
		want float64
	}{
		{12.5, 13}, {12.4, 12}, {-12.5, -13}, {-12.4, -12}, {0, 0}, {2.5, 3},
	}
	for _, c := range cases {
		if got := roundHalfAway(c.in); got != c.want {
			t.Fatalf("roundHalfAway(%v)=%v want %v", c.in, got, c.want)
		}
	}
}
