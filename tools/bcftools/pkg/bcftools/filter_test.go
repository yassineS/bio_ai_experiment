package bcftools

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func mkVariant() *vcf.Variant {
	return &vcf.Variant{
		Chrom:  "chr1",
		Pos:    100,
		Ref:    "A",
		Alt:    []string{"T"},
		Qual:   30,
		Filter: []string{"PASS"},
		Info: map[string]string{
			"DP": "42",
			"AF": "0.25",
			"H2": "",
		},
	}
}

func TestFilterSimpleCompare(t *testing.T) {
	cases := []struct {
		expr string
		want bool
	}{
		{`INFO/DP > 30`, true},
		{`INFO/DP > 50`, false},
		{`INFO/DP == 42`, true},
		{`INFO/DP != 42`, false},
		{`INFO/DP <= 42`, true},
		{`INFO/DP >= 100`, false},
		{`INFO/DP < 100`, true},
		{`FILTER == "PASS"`, true},
		{`FILTER == "q10"`, false},
		{`FILTER = "PASS"`, true}, // single-equal spelling
		{`INFO/AF < 0.5`, true},
		{`INFO/AF > 0.5`, false},
	}
	v := mkVariant()
	for _, tc := range cases {
		t.Run(tc.expr, func(t *testing.T) {
			f, err := CompileFilter(tc.expr)
			if err != nil {
				t.Fatal(err)
			}
			if got := f.Eval(v); got != tc.want {
				t.Fatalf("expr %q: got %v, want %v", tc.expr, got, tc.want)
			}
		})
	}
}

func TestFilterAndOrNot(t *testing.T) {
	v := mkVariant()
	cases := []struct {
		expr string
		want bool
	}{
		{`INFO/DP > 30 && FILTER="PASS"`, true},
		{`INFO/DP > 100 || FILTER="PASS"`, true},
		{`!(INFO/DP > 30)`, false},
		{`INFO/DP > 30 && !(FILTER="PASS")`, false},
		{`(INFO/DP > 100 || INFO/AF < 0.5) && FILTER="PASS"`, true},
	}
	for _, tc := range cases {
		t.Run(tc.expr, func(t *testing.T) {
			f, err := CompileFilter(tc.expr)
			if err != nil {
				t.Fatal(err)
			}
			if got := f.Eval(v); got != tc.want {
				t.Fatalf("expr %q: got %v, want %v", tc.expr, got, tc.want)
			}
		})
	}
}

func TestFilterFlagInfo(t *testing.T) {
	v := mkVariant()
	// Flag tag H2 is "present" — used as a bare reference.
	f, err := CompileFilter(`INFO/H2`)
	if err != nil {
		t.Fatal(err)
	}
	if !f.Eval(v) {
		t.Fatal("flag H2 should be truthy")
	}
	delete(v.Info, "H2")
	if f.Eval(v) {
		t.Fatal("absent flag should be falsy")
	}
}

func TestFilterNilSafe(t *testing.T) {
	var f *Filter
	if !f.Eval(mkVariant()) {
		t.Fatal("nil filter should accept all variants")
	}
}

func TestFilterBareIdentString(t *testing.T) {
	v := mkVariant()
	// bare identifier (no quotes) is treated as a string for comparison.
	f, err := CompileFilter(`FILTER == PASS`)
	if err != nil {
		t.Fatal(err)
	}
	if !f.Eval(v) {
		t.Fatal("bare ident comparison failed")
	}
}

func TestFilterBoolean(t *testing.T) {
	v := mkVariant()
	f, _ := CompileFilter(`true`)
	if !f.Eval(v) {
		t.Fatal("true should be truthy")
	}
	f, _ = CompileFilter(`false`)
	if f.Eval(v) {
		t.Fatal("false should be falsy")
	}
}

func TestFilterErrors(t *testing.T) {
	bad := []string{
		``,
		`(`,
		`INFO/DP >`,
		`"unterminated`,
		`123abc`,
		`(INFO/DP > 1`,
		`INFO/DP > 1 trailing`,
	}
	for _, b := range bad {
		if _, err := CompileFilter(b); err == nil {
			t.Errorf("expected error for %q", b)
		}
	}
}

func TestFilterNegativeNumber(t *testing.T) {
	v := mkVariant()
	v.Info["X"] = "-5"
	f, err := CompileFilter(`INFO/X < 0`)
	if err != nil {
		t.Fatal(err)
	}
	if !f.Eval(v) {
		t.Fatal("INFO/X (-5) should be < 0")
	}
}

func TestFilterCommaListUsesFirst(t *testing.T) {
	v := mkVariant()
	v.Info["AC"] = "3,5,7"
	f, err := CompileFilter(`INFO/AC > 2`)
	if err != nil {
		t.Fatal(err)
	}
	if !f.Eval(v) {
		t.Fatal("first element of AC should drive comparison")
	}
}

func TestFilterMissingInfoIsFalsy(t *testing.T) {
	v := mkVariant()
	f, err := CompileFilter(`INFO/NOPE > 0`)
	if err != nil {
		t.Fatal(err)
	}
	if f.Eval(v) {
		t.Fatal("missing INFO key should make numeric compare false")
	}
}

func TestFilterStringCompareLexical(t *testing.T) {
	v := mkVariant()
	v.Info["S"] = "banana"
	f, err := CompileFilter(`INFO/S < "cherry"`)
	if err != nil {
		t.Fatal(err)
	}
	if !f.Eval(v) {
		t.Fatal("banana < cherry should be true")
	}
}

func TestFilterScientificNotation(t *testing.T) {
	v := mkVariant()
	v.Info["X"] = "1e-3"
	f, err := CompileFilter(`INFO/X < 1e-2`)
	if err != nil {
		t.Fatal(err)
	}
	if !f.Eval(v) {
		t.Fatal("1e-3 should be less than 1e-2")
	}
}

func TestSplitCommaList(t *testing.T) {
	got := SplitCommaList("a, b ,c")
	want := []string{"a", "b", "c"}
	if len(got) != len(want) {
		t.Fatalf("len: %d vs %d", len(got), len(want))
	}
	for i := range got {
		if got[i] != want[i] {
			t.Errorf("got[%d]=%q want %q", i, got[i], want[i])
		}
	}
	if SplitCommaList("") != nil {
		t.Fatal("empty input should give nil")
	}
}

// mkSampleVariant builds a 3-sample record mirroring basic.vcf pos 400
// (1/2, 0/1, 0/0) for exercising the per-sample GT/FORMAT filter paths.
func mkSampleVariant() *vcf.Variant {
	return &vcf.Variant{
		Chrom:  "chr1",
		Pos:    400,
		ID:     "rs4",
		Ref:    "C",
		Alt:    []string{"T", "G"},
		Qual:   60,
		Filter: []string{"PASS"},
		Info:   map[string]string{"DP": "100"},
		Format: []string{"GT", "DP", "GQ"},
		Samples: []vcf.Sample{
			{Name: "S1", Data: map[string]string{"GT": "1/2", "DP": "40", "GQ": "60"}},
			{Name: "S2", Data: map[string]string{"GT": "0/1", "DP": "30", "GQ": "60"}},
			{Name: "S3", Data: map[string]string{"GT": "0/0", "DP": "30", "GQ": "60"}},
		},
	}
}

func TestFilterColumnOperands(t *testing.T) {
	v := mkSampleVariant()
	cases := []struct {
		expr string
		want bool
	}{
		{`CHROM="chr1"`, true},
		{`CHROM="chr2"`, false},
		{`POS=400`, true},
		{`POS>200`, true},
		{`POS<400`, false},
		{`ID="rs4"`, true},
		{`ID="rs1"`, false},
		{`REF="C"`, true},
		{`REF="A"`, false},
		{`ALT="T"`, true},
		{`ALT="G"`, true},
		{`ALT="A"`, false},
		{`ALT[0]="T"`, true},
		{`ALT[1]="G"`, true},
		{`ALT[1]="T"`, false},
	}
	for _, tc := range cases {
		t.Run(tc.expr, func(t *testing.T) {
			f, err := CompileFilter(tc.expr)
			if err != nil {
				t.Fatal(err)
			}
			if got := f.Eval(v); got != tc.want {
				t.Fatalf("expr %q: got %v want %v", tc.expr, got, tc.want)
			}
		})
	}
}

func TestFilterGTClass(t *testing.T) {
	v := mkSampleVariant() // genotypes: 1/2 (het-alt), 0/1 (het-ref), 0/0 (hom-ref)
	cases := []struct {
		expr string
		want bool
	}{
		{`GT="het"`, true},  // 1/2 and 0/1 are het
		{`GT="hom"`, true},  // 0/0 is hom
		{`GT="ref"`, true},  // 0/0 is ref
		{`GT="alt"`, true},  // 1/2 (and 0/1 has an alt) -> alt
		{`GT="RR"`, true},   // 0/0 -> rr
		{`GT="AA"`, false},  // no hom-alt sample
		{`GT="ra"`, true},   // 0/1 -> ra
		{`GT="Aa"`, true},   // 1/2 -> het-alt
		{`GT="aA"`, true},   // synonym of Aa
		{`GT="mis"`, false}, // no missing genotype
		{`GT="0/1"`, true},  // explicit form present
		{`GT="1/2"`, true},
		{`GT="2/2"`, false},
	}
	for _, tc := range cases {
		t.Run(tc.expr, func(t *testing.T) {
			f, err := CompileFilter(tc.expr)
			if err != nil {
				t.Fatal(err)
			}
			if got := f.Eval(v); got != tc.want {
				t.Fatalf("expr %q: got %v want %v", tc.expr, got, tc.want)
			}
		})
	}
}

func TestFilterFormatPerSample(t *testing.T) {
	v := mkSampleVariant() // DP: 40,30,30  GQ: 60,60,60
	cases := []struct {
		expr string
		want bool
	}{
		{`FORMAT/DP>35`, true}, // S1=40
		{`FMT/DP>50`, false},   // none
		{`FORMAT/GQ=60`, true}, // all
		{`FMT/GQ>60`, false},   // none
		{`FORMAT/DP<35`, true}, // S2,S3=30
	}
	for _, tc := range cases {
		t.Run(tc.expr, func(t *testing.T) {
			f, err := CompileFilter(tc.expr)
			if err != nil {
				t.Fatal(err)
			}
			if got := f.Eval(v); got != tc.want {
				t.Fatalf("expr %q: got %v want %v", tc.expr, got, tc.want)
			}
		})
	}
}

func TestFilterAggregations(t *testing.T) {
	v := mkSampleVariant() // DP: 40,30,30 (sum=100, max=40, min=30, avg=100/3)
	cases := []struct {
		expr string
		want bool
	}{
		{`N_PASS(GT="het")=2`, true}, // 1/2 and 0/1
		{`N_PASS(GT="hom")=1`, true}, // 0/0
		{`COUNT(GT="het")>1`, true},
		{`F_PASS(GT="het")>0.5`, true}, // 2/3
		{`SUM(FMT/DP)=100`, true},
		{`MAX(FMT/DP)=40`, true},
		{`MIN(FMT/DP)=30`, true},
		{`AVG(FMT/GQ)=60`, true},
		{`N_PASS(FMT/DP>35)=1`, true}, // only S1
	}
	for _, tc := range cases {
		t.Run(tc.expr, func(t *testing.T) {
			f, err := CompileFilter(tc.expr)
			if err != nil {
				t.Fatal(err)
			}
			if got := f.Eval(v); got != tc.want {
				t.Fatalf("expr %q: got %v want %v", tc.expr, got, tc.want)
			}
		})
	}
}

func TestFilterStrlenAbs(t *testing.T) {
	v := &vcf.Variant{
		Chrom: "chr1", Pos: 300, Ref: "GAT", Alt: []string{"G"}, Qual: 50,
		Info: map[string]string{},
	}
	cases := []struct {
		expr string
		want bool
	}{
		{`STRLEN(REF)=3`, true},
		{`STRLEN(REF)=1`, false},
		{`STRLEN(ALT)=1`, true},
		{`ABS(POS)=300`, true},
	}
	for _, tc := range cases {
		t.Run(tc.expr, func(t *testing.T) {
			f, err := CompileFilter(tc.expr)
			if err != nil {
				t.Fatal(err)
			}
			if got := f.Eval(v); got != tc.want {
				t.Fatalf("expr %q: got %v want %v", tc.expr, got, tc.want)
			}
		})
	}
}
