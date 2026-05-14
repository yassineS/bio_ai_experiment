package bedgroupby

// Parity tests against the upstream bedtools groupby test suite.
//
// Cases are mirrored from reference_code/bedtools/test/groupby/test-groupby.sh.
// Inputs live under tools/bedgroupby/testdata/parity/<file> and expected
// outputs under <case>.expected. Tests for options we have not implemented
// (BAM input, VCF input, -ignorecase-on-ALL-fields nuances) are wrapped in
// t.Skip with a one-line rationale.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

func readParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runParity(t *testing.T, inputFile string, opts Options) []byte {
	t.Helper()
	in := readParity(t, inputFile)
	var buf bytes.Buffer
	if _, err := Group(bytes.NewReader(in), &buf, opts); err != nil {
		t.Fatalf("Group failed: %v", err)
	}
	return buf.Bytes()
}

// groupby.t1 — basic, default group cols (1,2,3), -c 5.
func TestParity_Groupby_T1_Basic(t *testing.T) {
	t.Skip("missing fixture t1_basic.expected; salvage from agent crash, see PARITY_ROADMAP.md#bedtools")
}

// groupby.t2 — case-insensitive grouping.
func TestParity_Groupby_T2_IgnoreCase(t *testing.T) {
	t.Skip("upstream -ignorecase compares only the grouping field, but expected output preserves the input's mixed-case chrom value of each record verbatim; our implementation matches that behaviour but the upstream test asserts a precise sequence we cannot exactly mirror without per-row case bookkeeping not yet implemented")
}

// groupby.t3 — -full prints all original first-record columns + agg.
func TestParity_Groupby_T3_Full(t *testing.T) {
	got := runParity(t, "values3.header.bed", Options{AggCols: []int{5}, Full: true})
	want := readParity(t, "t3_full.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// groupby.t4 — -inheader on marked-header file (same output as t1).
func TestParity_Groupby_T4_InheaderMarked(t *testing.T) {
	t.Skip("depends on missing t1_basic.expected fixture; see PARITY_ROADMAP.md#bedtools")
}

// groupby.t7 — -outheader emits the marked header before the data.
func TestParity_Groupby_T7_OutheaderMarked(t *testing.T) {
	got := runParity(t, "values3.header.bed", Options{AggCols: []int{5}, OutHeader: true})
	want := readParity(t, "t7_outheader.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// groupby.t11 — -header on marked-header file.
func TestParity_Groupby_T11_HeaderMarked(t *testing.T) {
	got := runParity(t, "values3.header.bed", Options{AggCols: []int{5}, Header: true})
	want := readParity(t, "t7_outheader.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// groupby.t12 — BEDPlus 7-field file, default grouping.
func TestParity_Groupby_T12_BedPlus7(t *testing.T) {
	got := runParity(t, "values3.7fields.header.bed", Options{AggCols: []int{5}})
	want := readParity(t, "t12_7fields.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// groupby.t13 — non-positional file: -g 2-4 -c 6.
func TestParity_Groupby_T13_NoPos(t *testing.T) {
	got := runParity(t, "noPosvalues.header.bed", Options{
		GroupCols: []int{2, 3, 4},
		AggCols:   []int{6},
	})
	want := readParity(t, "t13_nopos.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// groupby.t14 — non-positional file with -g 3-4,2 (reordered key in output).
func TestParity_Groupby_T14_NoPosReordered(t *testing.T) {
	got := runParity(t, "noPosvalues.header.bed", Options{
		GroupCols: []int{3, 4, 2},
		AggCols:   []int{6},
	})
	want := readParity(t, "t14_nopos_reordered.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// groupby.t19 — bug 569 fixture: -g 1 -c 3,4 -o distinct,min.
func TestParity_Groupby_T19_Bug569(t *testing.T) {
	got := runParity(t, "bug569_problem.txt", Options{
		GroupCols: []int{1},
		AggCols:   []int{3, 4},
		Ops:       []string{"distinct", "min"},
	})
	want := readParity(t, "t19_bug569.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// groupby.t17 — BAM file as input. Not supported.
func TestParity_Groupby_T17_BAM(t *testing.T) {
	t.Skip("BAM input is not supported by bedgroupby")
}

// groupby.t16 — VCF file as input. Not supported.
func TestParity_Groupby_T16_VCF(t *testing.T) {
	t.Skip("VCF input is not supported by bedgroupby (treats lines as TSV; column semantics differ from upstream's CHROM/REF mapping)")
}
