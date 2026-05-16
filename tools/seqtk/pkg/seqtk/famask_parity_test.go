package seqtk

// Byte-for-byte parity tests for `seqtk famask` and `seqtk mergefa`
// against the upstream C reference implementation
// (reference_code/seqtk v1.5-r133). Fixtures live under
// tools/seqtk/testdata/parity/ and were produced by running the
// pinned upstream binary with the same input pairs and flags as the
// Go invocations below.

import (
	"bytes"
	"io"
	"testing"
)

// TestParity_Seqtk_Famask_Simple is the smallest readable fixture:
// two records, one with an X+x mix, one with N overwrites.
func TestParity_Seqtk_Famask_Simple(t *testing.T) {
	src := readParityFile(t, "famask_simple_src.fa")
	mask := readParityFile(t, "famask_simple_mask.fa")
	var out bytes.Buffer
	if err := famaskImpl(bytes.NewReader(src), bytes.NewReader(mask), &out, io.Discard); err != nil {
		t.Fatalf("Famask: %v", err)
	}
	want := readParityFile(t, "famask_simple.expected.fa")
	mustEqualBytes(t, "famask simple", out.Bytes(), want)
}

// TestParity_Seqtk_Famask_WrapAndMismatch exercises the 60-base wrap
// and the unequal-length warn path. The expected output came from
// `seqtk famask famask_src.fa famask_mask.fa 2>/dev/null`.
func TestParity_Seqtk_Famask_WrapAndMismatch(t *testing.T) {
	src := readParityFile(t, "famask_src.fa")
	mask := readParityFile(t, "famask_mask.fa")
	var out bytes.Buffer
	// We discard the warn stream here (the parity assertion is on
	// stdout-equivalent only — upstream's stderr text isn't part of
	// the byte-parity bar).
	if err := famaskImpl(bytes.NewReader(src), bytes.NewReader(mask), &out, io.Discard); err != nil {
		t.Fatalf("Famask: %v", err)
	}
	want := readParityFile(t, "famask.expected.fa")
	mustEqualBytes(t, "famask wrap+mismatch", out.Bytes(), want)
}

// TestParity_Seqtk_Mergefa_DefaultFasta covers the default OR-merge
// across A/C/G/T pairs and an N row.
func TestParity_Seqtk_Mergefa_DefaultFasta(t *testing.T) {
	a := readParityFile(t, "mergefa_a.fa")
	b := readParityFile(t, "mergefa_b.fa")
	var out bytes.Buffer
	if err := mergefaImpl(bytes.NewReader(a), bytes.NewReader(b), &out, io.Discard, MergefaOptions{}); err != nil {
		t.Fatalf("Mergefa default: %v", err)
	}
	want := readParityFile(t, "mergefa_default.expected.fa")
	mustEqualBytes(t, "mergefa default", out.Bytes(), want)
}

// TestParity_Seqtk_Mergefa_IntersectFasta covers `-i`.
func TestParity_Seqtk_Mergefa_IntersectFasta(t *testing.T) {
	a := readParityFile(t, "mergefa_a.fa")
	b := readParityFile(t, "mergefa_b.fa")
	var out bytes.Buffer
	if err := mergefaImpl(bytes.NewReader(a), bytes.NewReader(b), &out, io.Discard, MergefaOptions{Intersect: true}); err != nil {
		t.Fatalf("Mergefa -i: %v", err)
	}
	want := readParityFile(t, "mergefa_i.expected.fa")
	mustEqualBytes(t, "mergefa -i", out.Bytes(), want)
}

// TestParity_Seqtk_Mergefa_MaskFasta covers `-m` (lowercases on N).
func TestParity_Seqtk_Mergefa_MaskFasta(t *testing.T) {
	a := readParityFile(t, "mergefa_a.fa")
	b := readParityFile(t, "mergefa_b.fa")
	var out bytes.Buffer
	if err := mergefaImpl(bytes.NewReader(a), bytes.NewReader(b), &out, io.Discard, MergefaOptions{Mask: true}); err != nil {
		t.Fatalf("Mergefa -m: %v", err)
	}
	want := readParityFile(t, "mergefa_m.expected.fa")
	mustEqualBytes(t, "mergefa -m", out.Bytes(), want)
}

// TestParity_Seqtk_Mergefa_HaploidFasta covers `-h` (lowercases hets).
func TestParity_Seqtk_Mergefa_HaploidFasta(t *testing.T) {
	a := readParityFile(t, "mergefa_a.fa")
	b := readParityFile(t, "mergefa_b.fa")
	var out bytes.Buffer
	if err := mergefaImpl(bytes.NewReader(a), bytes.NewReader(b), &out, io.Discard, MergefaOptions{Haploid: true}); err != nil {
		t.Fatalf("Mergefa -h: %v", err)
	}
	want := readParityFile(t, "mergefa_h.expected.fa")
	mustEqualBytes(t, "mergefa -h", out.Bytes(), want)
}

// TestParity_Seqtk_Mergefa_QualLoweringFastq covers `-q INT` with
// FASTQ input: low-quality bases are lowercased before the merge.
func TestParity_Seqtk_Mergefa_QualLoweringFastq(t *testing.T) {
	a := readParityFile(t, "mergefa_a.fq")
	b := readParityFile(t, "mergefa_b.fq")
	var out bytes.Buffer
	if err := mergefaImpl(bytes.NewReader(a), bytes.NewReader(b), &out, io.Discard, MergefaOptions{Quality: 20}); err != nil {
		t.Fatalf("Mergefa -q 20: %v", err)
	}
	want := readParityFile(t, "mergefa_q20.expected.fa")
	mustEqualBytes(t, "mergefa -q 20", out.Bytes(), want)
}

// TestParity_Seqtk_Mergefa_LineWrap pins the 60-column output wrap
// against a 125-base homozygous input.
func TestParity_Seqtk_Mergefa_LineWrap(t *testing.T) {
	a := readParityFile(t, "mergefa_long_a.fa")
	b := readParityFile(t, "mergefa_long_b.fa")
	var out bytes.Buffer
	if err := mergefaImpl(bytes.NewReader(a), bytes.NewReader(b), &out, io.Discard, MergefaOptions{}); err != nil {
		t.Fatalf("Mergefa long: %v", err)
	}
	want := readParityFile(t, "mergefa_long.expected.fa")
	mustEqualBytes(t, "mergefa long wrap", out.Bytes(), want)
}
