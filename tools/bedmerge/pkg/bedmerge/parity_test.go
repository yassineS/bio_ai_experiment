package bedmerge

// Parity tests against the upstream bedtools merge test suite.
//
// Cases are mirrored from reference_code/bedtools/test/merge/test-merge.sh.
// Inputs live under tools/bedmerge/testdata/parity/<case>.bed and expected
// outputs under <case>.expected.bed. Tests that exercise features bedmerge
// does not implement (VCF input) are wrapped in t.Skip with a one-line
// rationale rather than being deleted; the `-S` strand filter (t16/t17) and
// the custom `-delim` list delimiter (t10) are implemented and asserted.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

func readMergeParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runMergeParity(t *testing.T, inputFile string, opts MergeOptions) []byte {
	t.Helper()
	in := readMergeParity(t, inputFile)
	var buf bytes.Buffer
	if _, err := Merge(bytes.NewReader(in), &buf, opts); err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	return buf.Bytes()
}

func mustParseOps(t *testing.T, cols, ops string) *ColumnOps {
	t.Helper()
	co, err := ParseColumnOps(cols, ops)
	if err != nil {
		t.Fatalf("ParseColumnOps(%q, %q): %v", cols, ops, err)
	}
	return co
}

// merge.t1 — basic merge, BED3 output.
func TestParity_Merge_T1_Basic(t *testing.T) {
	got := runMergeParity(t, "a.bed", MergeOptions{})
	want := readMergeParity(t, "t1_basic.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t3 — count of merged intervals via -c 1 -o count.
func TestParity_Merge_T3_Count(t *testing.T) {
	got := runMergeParity(t, "a.bed", MergeOptions{ColumnOps: mustParseOps(t, "1", "count")})
	want := readMergeParity(t, "t3_count.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t5 — collapse names via -c 4 -o collapse.
func TestParity_Merge_T5_CollapseNames(t *testing.T) {
	got := runMergeParity(t, "a.names.bed", MergeOptions{ColumnOps: mustParseOps(t, "4", "collapse")})
	want := readMergeParity(t, "t5_collapse_names.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t6 — collapse names + sum scores.
func TestParity_Merge_T6_CollapseSum(t *testing.T) {
	got := runMergeParity(t, "a.full.bed", MergeOptions{ColumnOps: mustParseOps(t, "4,5", "collapse,sum")})
	want := readMergeParity(t, "t6_collapse_sum.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t7 — count, sum (single op applied per column).
func TestParity_Merge_T7_CountSum(t *testing.T) {
	got := runMergeParity(t, "a.full.bed", MergeOptions{ColumnOps: mustParseOps(t, "5", "count,sum")})
	want := readMergeParity(t, "t7_count_sum.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t8 — collapse + sum + count over three columns.
func TestParity_Merge_T8_ThreeOps(t *testing.T) {
	got := runMergeParity(t, "a.full.bed", MergeOptions{ColumnOps: mustParseOps(t, "4,5,4", "collapse,sum,count")})
	want := readMergeParity(t, "t8_three_ops.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t9a — stranded merge (-s) with collapse, sum, count.
func TestParity_Merge_T9a_Stranded(t *testing.T) {
	got := runMergeParity(t, "a.full.bed", MergeOptions{
		StrandSpec: true,
		ColumnOps:  mustParseOps(t, "4,5,6", "collapse,sum,count"),
	})
	want := readMergeParity(t, "t9a_stranded.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t9b — stranded merge (-s) keeping the strand column via collapse on col 6.
func TestParity_Merge_T9b_StrandedStrand(t *testing.T) {
	got := runMergeParity(t, "a.full.bed", MergeOptions{
		StrandSpec: true,
		ColumnOps:  mustParseOps(t, "4,5,6", "collapse,sum,collapse"),
	})
	want := readMergeParity(t, "t9b_stranded_strands.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t10 — custom delimiter via -delim "|" with -c 4 -o collapse.
func TestParity_Merge_T10_CustomDelim(t *testing.T) {
	co := mustParseOps(t, "4", "collapse")
	co.Delim = "|"
	got := runMergeParity(t, "a.names.bed", MergeOptions{ColumnOps: co})
	want := readMergeParity(t, "t10_custom_delim.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t13 — VCF input. bedmerge only accepts BED.
func TestParity_Merge_T13_VCFInput(t *testing.T) {
	t.Skip("unimplemented: VCF input; bedmerge is BED-only")
}

// merge.t14 — GFF input: features are auto-detected (1-based start/end in
// columns 4/5) and merged into BED3, byte-for-byte against bedtools v2.31.1.
func TestParity_Merge_T14_GFFInput(t *testing.T) {
	got := runMergeParity(t, "features.gff", MergeOptions{})
	want := readMergeParity(t, "t14_gff.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t15 — stranded merge with mixed '.' strands. Upstream's
// FileRecordMergeMgr drops UNKNOWN/`.` records under `-s` (see
// reference_code/bedtools/src/utils/FileRecordTools/FileRecordMergeMgr.cpp
// lines 47-58 + 96-129) and merges `+` and `-` independently, then
// emits the two streams in (chrom, start, end) order. bedmerge now
// matches this behaviour.
func TestParity_Merge_T15_MixedStrandsFanOut(t *testing.T) {
	got := runMergeParity(t, "mixedStrands.bed", MergeOptions{StrandSpec: true})
	want := readMergeParity(t, "t15_mixed_strands_s.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t16 / t17 — `-S +` / `-S -` keep only records on the given strand
// (dropping the other strand and `.`/unknown) and merge the survivors
// positionally, emitting BED3. bedmerge implements this via
// MergeOptions.StrandFilter.
func TestParity_Merge_T16_StrandFilterPlus(t *testing.T) {
	got := runMergeParity(t, "mixedStrands.bed", MergeOptions{StrandFilter: "+"})
	want := readMergeParity(t, "t16_strand_filter_plus.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
func TestParity_Merge_T17_StrandFilterMinus(t *testing.T) {
	got := runMergeParity(t, "mixedStrands.bed", MergeOptions{StrandFilter: "-"})
	want := readMergeParity(t, "t17_strand_filter_minus.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t20 — chromosome change handling (BED3 output, 4-col input ignored).
func TestParity_Merge_T20_ChromChange(t *testing.T) {
	got := runMergeParity(t, "b.bed", MergeOptions{})
	want := readMergeParity(t, "t20_chrom_change.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// merge.t21 — BED3 output from 6-col input under default options.
func TestParity_Merge_T21_BED3FromFull(t *testing.T) {
	got := runMergeParity(t, "a.full.bed", MergeOptions{})
	want := readMergeParity(t, "t21_bed3.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
