package seqtk

// Byte-for-byte parity tests for `seqtk fqchk` and `seqtk hety` against
// the upstream C reference implementation (reference_code/seqtk
// v1.5-r133). Fixtures live under tools/seqtk/testdata/parity/ and were
// produced by running the pinned upstream binary with the same input
// files and flags as the Go invocations below.

import (
	"bytes"
	"testing"
)

// TestParity_Seqtk_Fqchk_Default verifies that Fqchk emits the same
// bytes as `seqtk fqchk` (default qthres = 20) for a small mixed-length
// FASTQ.
func TestParity_Seqtk_Fqchk_Default(t *testing.T) {
	in := readParityFile(t, "fqchk_mixed.fq")
	var out bytes.Buffer
	if err := Fqchk(bytes.NewReader(in), &out, FqchkOptions{QThres: DefaultFqchkQThres}); err != nil {
		t.Fatalf("Fqchk: %v", err)
	}
	want := readParityFile(t, "fqchk_mixed_default.expected.txt")
	mustEqualBytes(t, "fqchk default (mixed.fq)", out.Bytes(), want)
}

// TestParity_Seqtk_Fqchk_Q0 verifies the qthres = 0 path where the
// trailing per-row columns become one %Qk per distinct observed quality.
func TestParity_Seqtk_Fqchk_Q0(t *testing.T) {
	in := readParityFile(t, "fqchk_mixed.fq")
	var out bytes.Buffer
	if err := Fqchk(bytes.NewReader(in), &out, FqchkOptions{QThres: 0}); err != nil {
		t.Fatalf("Fqchk -q0: %v", err)
	}
	want := readParityFile(t, "fqchk_mixed_q0.expected.txt")
	mustEqualBytes(t, "fqchk -q0 (mixed.fq)", out.Bytes(), want)
}

// TestParity_Seqtk_Fqchk_Q30 verifies the qthres = 30 split. Together
// with the default and -q0 cases this covers all three branches of
// upstream's fqc_aux output schema.
func TestParity_Seqtk_Fqchk_Q30(t *testing.T) {
	in := readParityFile(t, "fqchk_mixed.fq")
	var out bytes.Buffer
	if err := Fqchk(bytes.NewReader(in), &out, FqchkOptions{QThres: 30}); err != nil {
		t.Fatalf("Fqchk -q30: %v", err)
	}
	want := readParityFile(t, "fqchk_mixed_q30.expected.txt")
	mustEqualBytes(t, "fqchk -q30 (mixed.fq)", out.Bytes(), want)
}

// TestParity_Seqtk_Fqchk_SmallFqDefault re-uses the shared small.fq
// fixture (also used by the comp parity tests) as a second input —
// the records there have a different length / quality distribution.
func TestParity_Seqtk_Fqchk_SmallFqDefault(t *testing.T) {
	in := readParityFile(t, "small.fq")
	var out bytes.Buffer
	if err := Fqchk(bytes.NewReader(in), &out, FqchkOptions{QThres: DefaultFqchkQThres}); err != nil {
		t.Fatalf("Fqchk small.fq: %v", err)
	}
	want := readParityFile(t, "fqchk_small_default.expected.txt")
	mustEqualBytes(t, "fqchk default (small.fq)", out.Bytes(), want)
}

// TestParity_Seqtk_Fqchk_SmallFqQ0 covers -q0 on the shared small.fq.
// This adds a fixture with 3 distinct quality values (vs mixed.fq's 4),
// exercising the dynamic %Qk-column count.
func TestParity_Seqtk_Fqchk_SmallFqQ0(t *testing.T) {
	in := readParityFile(t, "small.fq")
	var out bytes.Buffer
	if err := Fqchk(bytes.NewReader(in), &out, FqchkOptions{QThres: 0}); err != nil {
		t.Fatalf("Fqchk -q0 small.fq: %v", err)
	}
	want := readParityFile(t, "fqchk_small_q0.expected.txt")
	mustEqualBytes(t, "fqchk -q0 (small.fq)", out.Bytes(), want)
}

// TestParity_Seqtk_Hety_Default exercises the default window (50000 bp)
// which causes the two short records in hety_basic.fa to each emit a
// single window via the i == l flush.
func TestParity_Seqtk_Hety_Default(t *testing.T) {
	in := readParityFile(t, "hety_basic.fa")
	var out bytes.Buffer
	if err := Hety(bytes.NewReader(in), &out, HetyOptions{WinSize: DefaultHetyWinSize, NStart: DefaultHetyNStart}); err != nil {
		t.Fatalf("Hety default: %v", err)
	}
	want := readParityFile(t, "hety_basic_default.expected.txt")
	mustEqualBytes(t, "hety default (hety_basic.fa)", out.Bytes(), want)
}

// TestParity_Seqtk_Hety_W30 covers the "window smaller than the
// sequence" path: multiple windows per record, and the partial-final
// window at i == l.
func TestParity_Seqtk_Hety_W30(t *testing.T) {
	in := readParityFile(t, "hety_basic.fa")
	var out bytes.Buffer
	if err := Hety(bytes.NewReader(in), &out, HetyOptions{WinSize: 30, NStart: DefaultHetyNStart}); err != nil {
		t.Fatalf("Hety -w30: %v", err)
	}
	want := readParityFile(t, "hety_basic_w30.expected.txt")
	mustEqualBytes(t, "hety -w30 (hety_basic.fa)", out.Bytes(), want)
}

// TestParity_Seqtk_Hety_W30_T3 changes the step (-t) so the windows
// overlap differently. Exercises the cnt[buf[y]]-- eviction path with
// a non-default step.
func TestParity_Seqtk_Hety_W30_T3(t *testing.T) {
	in := readParityFile(t, "hety_basic.fa")
	var out bytes.Buffer
	if err := Hety(bytes.NewReader(in), &out, HetyOptions{WinSize: 30, NStart: 3}); err != nil {
		t.Fatalf("Hety -w30 -t3: %v", err)
	}
	want := readParityFile(t, "hety_basic_w30_t3.expected.txt")
	mustEqualBytes(t, "hety -w30 -t3 (hety_basic.fa)", out.Bytes(), want)
}

// TestParity_Seqtk_Hety_W30_M covers the -m (lower-as-N) branch on a
// FASTA with no lowercase bases, confirming that -m is a no-op there
// and the output equals the no-m run.
func TestParity_Seqtk_Hety_W30_M(t *testing.T) {
	in := readParityFile(t, "hety_basic.fa")
	var out bytes.Buffer
	if err := Hety(bytes.NewReader(in), &out, HetyOptions{WinSize: 30, NStart: DefaultHetyNStart, IsLowerMask: true}); err != nil {
		t.Fatalf("Hety -w30 -m: %v", err)
	}
	want := readParityFile(t, "hety_basic_w30_m.expected.txt")
	mustEqualBytes(t, "hety -w30 -m (hety_basic.fa)", out.Bytes(), want)
}

// TestParity_Seqtk_Hety_LowercaseMask exercises -m on a FASTA with
// actual lowercase bases, pinning the islower() -> N substitution
// against upstream.
func TestParity_Seqtk_Hety_LowercaseMask(t *testing.T) {
	in := readParityFile(t, "hety_lowercase.fa")
	var out bytes.Buffer
	if err := Hety(bytes.NewReader(in), &out, HetyOptions{WinSize: 6, NStart: 1, IsLowerMask: true}); err != nil {
		t.Fatalf("Hety -w6 -t1 -m: %v", err)
	}
	want := readParityFile(t, "hety_lowercase_w6_t1_m.expected.txt")
	mustEqualBytes(t, "hety -w6 -t1 -m (hety_lowercase.fa)", out.Bytes(), want)
}
