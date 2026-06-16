package bedmap

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// runMap is a tiny helper.
func runMap(t *testing.T, a, b string, opts Options) string {
	t.Helper()
	var buf bytes.Buffer
	if _, err := Map(strings.NewReader(a), strings.NewReader(b), &buf, opts); err != nil {
		t.Fatalf("Map: %v", err)
	}
	return buf.String()
}

func TestMap_DefaultSum(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr1\t10\t20\tn\t5\t+\nchr1\t30\t40\tn\t7\t+\n"
	// Default: -c 5 -o sum. 5+7=12.
	got := runMap(t, a, b, Options{})
	want := "chr1\t0\t100\t12\n"
	if got != want {
		t.Errorf("default sum: want %q got %q", want, got)
	}
}

func TestMap_NoOverlapNull(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr2\t10\t20\tn\t5\t+\n"
	got := runMap(t, a, b, Options{})
	if got != "chr1\t0\t100\t.\n" {
		t.Errorf("null placeholder: got %q", got)
	}
}

func TestMap_NoOverlapCustomNull(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr2\t10\t20\tn\t5\t+\n"
	got := runMap(t, a, b, Options{Null: "NA"})
	if got != "chr1\t0\t100\tNA\n" {
		t.Errorf("custom null: got %q", got)
	}
}

func TestMap_CountZeroOnNoOverlap(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr2\t10\t20\tn\t5\t+\n"
	// count and count_distinct should give 0, not null.
	got := runMap(t, a, b, Options{Ops: []string{"count"}})
	if got != "chr1\t0\t100\t0\n" {
		t.Errorf("count zero: got %q", got)
	}
	got = runMap(t, a, b, Options{Ops: []string{"count_distinct"}})
	if got != "chr1\t0\t100\t0\n" {
		t.Errorf("count_distinct zero: got %q", got)
	}
}

func TestMap_CollapseAndDistinct(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr1\t10\t20\ta\t1\t+\nchr1\t30\t40\ta\t2\t+\nchr1\t50\t60\tb\t3\t+\n"
	got := runMap(t, a, b, Options{Columns: []int{4}, Ops: []string{"collapse"}})
	if got != "chr1\t0\t100\ta,a,b\n" {
		t.Errorf("collapse: got %q", got)
	}
	got = runMap(t, a, b, Options{Columns: []int{4}, Ops: []string{"distinct"}})
	if got != "chr1\t0\t100\ta,b\n" {
		t.Errorf("distinct: got %q", got)
	}
	// Custom delim only swaps for collapse/distinct.
	got = runMap(t, a, b, Options{Columns: []int{4}, Ops: []string{"collapse"}, Delim: "|"})
	if got != "chr1\t0\t100\ta|a|b\n" {
		t.Errorf("collapse with delim: got %q", got)
	}
}

func TestMap_MultipleOpsOneColumn(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr1\t10\t20\tn\t1\t+\nchr1\t30\t40\tn\t5\t+\n"
	got := runMap(t, a, b, Options{Columns: []int{5}, Ops: []string{"min", "max"}})
	if got != "chr1\t0\t100\t1\t5\n" {
		t.Errorf("min,max: got %q", got)
	}
}

func TestMap_MultipleColsOneOp(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr1\t10\t20\ta\t2\t+\nchr1\t30\t40\tb\t3\t+\n"
	got := runMap(t, a, b, Options{Columns: []int{4, 5}, Ops: []string{"collapse"}})
	if got != "chr1\t0\t100\ta,b\t2,3\n" {
		t.Errorf("multi-col one op: got %q", got)
	}
}

func TestMap_StrandFilters(t *testing.T) {
	a := "chr1\t0\t100\tx\t0\t+\n"
	b := "chr1\t10\t20\ta\t1\t+\nchr1\t30\t40\tb\t2\t-\n"
	got := runMap(t, a, b, Options{SameStrand: true})
	if got != "chr1\t0\t100\tx\t0\t+\t1\n" {
		t.Errorf("same strand: got %q", got)
	}
	got = runMap(t, a, b, Options{OppositeStrand: true})
	if got != "chr1\t0\t100\tx\t0\t+\t2\n" {
		t.Errorf("opposite strand: got %q", got)
	}
}

func TestMap_FractionFilters(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr1\t90\t110\tn\t5\t+\n" // 10bp overlap → fracA=0.1, fracB=0.5
	got := runMap(t, a, b, Options{FractionA: 0.5})
	if got != "chr1\t0\t100\t.\n" {
		t.Errorf("fracA reject: got %q", got)
	}
	got = runMap(t, a, b, Options{FractionA: 0.05})
	if got != "chr1\t0\t100\t5\n" {
		t.Errorf("fracA accept: got %q", got)
	}
}

func TestMap_PreservesAColumns(t *testing.T) {
	a := "chr1\t0\t100\tnameA\t99\t+\n"
	b := "chr1\t10\t20\tn\t5\t+\n"
	got := runMap(t, a, b, Options{})
	if got != "chr1\t0\t100\tnameA\t99\t+\t5\n" {
		t.Errorf("preserve A: got %q", got)
	}
}

func TestMap_ValidateErrors(t *testing.T) {
	// Column 0 / negative columns are now validated by Map against the
	// database field count (matching upstream's column-range error), not by
	// Validate; verify Map rejects them.
	if _, err := Map(strings.NewReader("chr1\t0\t100\n"), strings.NewReader("chr1\t10\t20\tn\t5\n"), io.Discard, Options{Columns: []int{0}}); err == nil {
		t.Errorf("expected error for col 0")
	}
	opts := Options{Columns: []int{4, 5}, Ops: []string{"sum", "min", "max"}}
	if err := opts.Validate(); err == nil {
		t.Errorf("expected mismatched ops/cols error")
	}
	opts = Options{SameStrand: true, OppositeStrand: true}
	if err := opts.Validate(); err == nil {
		t.Errorf("expected -s/-S mutex error")
	}
}

func TestMap_BadBedRecords(t *testing.T) {
	_, err := Map(strings.NewReader("chr1\tbad\t10\n"), strings.NewReader(""), io.Discard, Options{})
	if err == nil {
		t.Errorf("expected error on bad A")
	}
	_, err = Map(strings.NewReader(""), strings.NewReader("chr1\tbad\t10\n"), io.Discard, Options{})
	if err == nil {
		t.Errorf("expected error on bad B")
	}
	// Column-out-of-range in B.
	_, err = Map(strings.NewReader("chr1\t0\t100\n"), strings.NewReader("chr1\t10\t20\n"), io.Discard, Options{Columns: []int{8}})
	if err == nil {
		t.Errorf("expected error for missing column 8 in B")
	}
}

func TestMap_HeadersAndCommentsSkipped(t *testing.T) {
	a := "#header\nchr1\t0\t100\n"
	b := "track name=foo\nchr1\t10\t20\tn\t5\t+\n"
	got := runMap(t, a, b, Options{})
	if got != "chr1\t0\t100\t5\n" {
		t.Errorf("headers: got %q", got)
	}
}

func TestMap_BEDPlusColumnExtraction(t *testing.T) {
	// BED + extra columns. Extract column 7 with collapse.
	a := "chr1\t0\t100\n"
	b := "chr1\t10\t20\ta\t1\t+\thello\n" +
		"chr1\t30\t40\tb\t2\t+\tworld\n"
	got := runMap(t, a, b, Options{Columns: []int{7}, Ops: []string{"collapse"}})
	if got != "chr1\t0\t100\thello,world\n" {
		t.Errorf("BED+ col 7: got %q", got)
	}
}

func TestMap_Reciprocal(t *testing.T) {
	a := "chr1\t0\t100\n"
	b := "chr1\t90\t110\tn\t5\t+\n" // fracA=0.1, fracB=0.5
	// With -r and -f 0.5 -F 0.5, both must hold ⇒ rejected.
	got := runMap(t, a, b, Options{FractionA: 0.5, FractionB: 0.5, Reciprocal: true})
	if got != "chr1\t0\t100\t.\n" {
		t.Errorf("reciprocal reject: got %q", got)
	}
}

func TestMap_EmptyBColumnRangeError(t *testing.T) {
	// An empty B has zero fields, so the default -c 5 is out of range. Upstream
	// `bedtools map` errors here ("... only has fields 1 - 0."); we match it
	// rather than emitting a null row.
	a := "chr1\t0\t100\n"
	var buf bytes.Buffer
	_, err := Map(strings.NewReader(a), strings.NewReader(""), &buf, Options{})
	if err == nil {
		t.Fatalf("empty B: expected column-range error, got output %q", buf.String())
	}
	if !strings.Contains(err.Error(), "only has fields 1 - 0") {
		t.Errorf("empty B: unexpected error %q", err.Error())
	}
}

func TestMap_EmptyBWithValidColumn(t *testing.T) {
	// When B is empty but A itself supplies enough columns to satisfy -c via a
	// B that does have fields, a no-overlap row still yields the null value.
	a := "chr1\t0\t100\n"
	b := "chr2\t10\t20\tn\t5\n"
	got := runMap(t, a, b, Options{})
	if got != "chr1\t0\t100\t.\n" {
		t.Errorf("no-overlap null: got %q", got)
	}
}
