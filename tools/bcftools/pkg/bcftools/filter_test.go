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
