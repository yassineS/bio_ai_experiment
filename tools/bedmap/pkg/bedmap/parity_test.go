package bedmap

// Parity tests against the upstream `bedtools map` test suite.
//
// Cases mirrored from reference_code/bedtools/test/map/test-map.sh.
// Inputs vendored under tools/bedmap/testdata/parity/.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

func readParityFixture(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runParity(t *testing.T, aFile, bFile string, opts Options) []byte {
	t.Helper()
	a := readParityFixture(t, aFile)
	b := readParityFixture(t, bFile)
	var buf bytes.Buffer
	if _, err := Map(bytes.NewReader(a), bytes.NewReader(b), &buf, opts); err != nil {
		t.Fatalf("Map: %v", err)
	}
	return buf.Bytes()
}

// map.t01 — defaults (-c 5 -o sum).
func TestParity_Map_T01_Default(t *testing.T) {
	got := runParity(t, "ivls.bed", "values.bed", Options{})
	want := readParityFixture(t, "t01_default.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// map.t02 — explicit -o sum.
func TestParity_Map_T02_Sum(t *testing.T) {
	got := runParity(t, "ivls.bed", "values.bed", Options{Ops: []string{"sum"}})
	want := readParityFixture(t, "t02_sum.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// map.t03 — count.
func TestParity_Map_T03_Count(t *testing.T) {
	got := runParity(t, "ivls.bed", "values.bed", Options{Ops: []string{"count"}})
	want := readParityFixture(t, "t03_count.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// map.t04 — mean. Upstream formats integer-valued means as plain integers,
// matching bedmerge's formatNum behaviour. Our values are integer means here.
func TestParity_Map_T04_Mean(t *testing.T) {
	got := runParity(t, "ivls.bed", "values.bed", Options{Ops: []string{"mean"}})
	want := readParityFixture(t, "t04_mean.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// map.t05 — max.
func TestParity_Map_T05_Max(t *testing.T) {
	got := runParity(t, "ivls.bed", "values.bed", Options{Ops: []string{"max"}})
	want := readParityFixture(t, "t05_max.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// map.t06 — min.
func TestParity_Map_T06_Min(t *testing.T) {
	got := runParity(t, "ivls.bed", "values.bed", Options{Ops: []string{"min"}})
	want := readParityFixture(t, "t06_min.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// map.t07 — mode, using values2.bed (which has duplicates).
func TestParity_Map_T07_Mode(t *testing.T) {
	got := runParity(t, "ivls.bed", "values2.bed", Options{Ops: []string{"mode"}})
	want := readParityFixture(t, "t07_mode.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// map.t08 — antimode.
func TestParity_Map_T08_Antimode(t *testing.T) {
	got := runParity(t, "ivls.bed", "values2.bed", Options{Ops: []string{"antimode"}})
	want := readParityFixture(t, "t08_antimode.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// map.t09 — column extraction from BED+: -c 7 -o collapse on values4.bed
// (7th column is a signed integer including negatives).
func TestParity_Map_T09_Collapse_BEDPlus(t *testing.T) {
	got := runParity(t, "ivls.bed", "values4.bed", Options{Columns: []int{7}, Ops: []string{"collapse"}})
	want := readParityFixture(t, "t09_collapse.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// map.t10 — min on the signed col 7 of values4.bed.
func TestParity_Map_T10_MinNegative(t *testing.T) {
	got := runParity(t, "ivls.bed", "values4.bed", Options{Columns: []int{7}, Ops: []string{"min"}})
	want := readParityFixture(t, "t10_min_neg.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// map.t11 — absmin on the signed col 7 of values4.bed.
func TestParity_Map_T11_AbsMin(t *testing.T) {
	got := runParity(t, "ivls.bed", "values4.bed", Options{Columns: []int{7}, Ops: []string{"absmin"}})
	want := readParityFixture(t, "t11_absmin.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// map.t13 — absmax on the signed col 7 of values4.bed.
func TestParity_Map_T13_AbsMax(t *testing.T) {
	got := runParity(t, "ivls.bed", "values4.bed", Options{Columns: []int{7}, Ops: []string{"absmax"}})
	want := readParityFixture(t, "t13_absmax.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// map.t14..t29 — GFF / VCF / BAM input not supported in BED-only port.
func TestParity_Map_T14_GFF(t *testing.T) {
	t.Skip("GFF input not yet supported in bedmap")
}
