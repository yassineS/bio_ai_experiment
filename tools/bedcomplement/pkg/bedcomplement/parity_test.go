package bedcomplement

// Parity tests against the upstream bedtools complement test suite.
//
// Cases are mirrored from reference_code/bedtools/test/complement/test-complement.sh.
// Each case has a small `.in.bed` and corresponding `.genome` plus an
// `.expected.bed`. Tests that exercise upstream features bedcomplement does
// not implement (notably -L, which restricts the output to chromosomes seen
// in the input) are wrapped in t.Skip with a one-line rationale.

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"
)

func readComplementParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runComplementParity(t *testing.T, inputFile, genomeFile string) []byte {
	t.Helper()
	in := readComplementParity(t, inputFile)
	g := readComplementParity(t, genomeFile)
	sizes, order, err := ReadChromSizes(bytes.NewReader(g))
	if err != nil {
		t.Fatalf("ReadChromSizes %s: %v", genomeFile, err)
	}
	var buf bytes.Buffer
	if _, err := Complement(bytes.NewReader(in), &buf, io.Discard, sizes, order); err != nil {
		t.Fatalf("Complement failed: %v", err)
	}
	return buf.Bytes()
}

// complement.t1 — basic baseline complement on a 20-bp chrom.
func TestParity_Complement_T1_Baseline(t *testing.T) {
	got := runComplementParity(t, "t1.in.bed", "t1.genome")
	want := readComplementParity(t, "t1.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// complement.t2 — both ends are covered, middle is the gap.
func TestParity_Complement_T2_EndsCovered(t *testing.T) {
	got := runComplementParity(t, "t2.in.bed", "t1.genome")
	want := readComplementParity(t, "t2.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// complement.t3 — middle is covered; emit leading and trailing gaps.
func TestParity_Complement_T3_MiddleCovered(t *testing.T) {
	got := runComplementParity(t, "t3.in.bed", "t1.genome")
	want := readComplementParity(t, "t3.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// complement.t4 — chromosome entirely covered: output is empty.
func TestParity_Complement_T4_EntirelyCovered(t *testing.T) {
	got := runComplementParity(t, "t4.in.bed", "t1.genome")
	want := readComplementParity(t, "t4.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%q\ngot:\n%q", want, got)
	}
}

// complement.t5 — chromosome with no intervals on it appears in full.
func TestParity_Complement_T5_NothingCovered(t *testing.T) {
	got := runComplementParity(t, "t5.in.bed", "t5.genome")
	want := readComplementParity(t, "t5.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// complement.t6 — issue #356: 4-column input, BED3 output, large chrom size.
func TestParity_Complement_T6_Issue356(t *testing.T) {
	got := runComplementParity(t, "t6.in.bed", "t6.genome")
	want := readComplementParity(t, "t6.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// complement.t7 — multiple chromosomes, both partially covered.
func TestParity_Complement_T7_MultipleChroms(t *testing.T) {
	got := runComplementParity(t, "t7.in.bed", "t5.genome")
	want := readComplementParity(t, "t7.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// complement.t8 — multiple chromosomes; chr1 fully covered → no chr1 rows.
func TestParity_Complement_T8_OneFullyCovered(t *testing.T) {
	got := runComplementParity(t, "t8.in.bed", "t5.genome")
	want := readComplementParity(t, "t8.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// complement.t9 — input interval (chr1 90-110) exceeds the chromosome
// length (chr1 = 100). Upstream emits both the warning
// `***** WARNING: chr1:90-110 exceeds the length of chromosome (chr1)`
// to stderr AND the clipped complement gap `chr1\t0\t90` on stdout.
// Our port now mirrors both halves.
func TestParity_Complement_T9_RecordExceedsChrom(t *testing.T) {
	in := []byte("chr1\t90\t110\n")
	g := []byte("chr1\t100\n")
	sizes, order, err := ReadChromSizes(bytes.NewReader(g))
	if err != nil {
		t.Fatalf("ReadChromSizes: %v", err)
	}
	var out, warn bytes.Buffer
	if _, err := Complement(bytes.NewReader(in), &out, &warn, sizes, order); err != nil {
		t.Fatalf("Complement: %v", err)
	}
	wantOut := "chr1\t0\t90\n"
	if got := out.String(); got != wantOut {
		t.Errorf("stdout mismatch.\nwant: %q\ngot:  %q", wantOut, got)
	}
	wantWarn := "***** WARNING: chr1:90-110 exceeds the length of chromosome (chr1)\n"
	if got := warn.String(); got != wantWarn {
		t.Errorf("stderr mismatch.\nwant: %q\ngot:  %q", wantWarn, got)
	}
}

// complement.t9b / t10 (script duplicates 'complement.t9' label) — issue #503,
// the -L flag limits output to chromosomes that had records in the input.
func TestParity_Complement_T9b_DashLLimit(t *testing.T) {
	in := readComplementParity(t, "issue_503.bed")
	g := readComplementParity(t, "issue_503.genome")
	sizes, order, err := ReadChromSizes(bytes.NewReader(g))
	if err != nil {
		t.Fatalf("ReadChromSizes: %v", err)
	}
	var buf bytes.Buffer
	if _, err := ComplementWithOptions(bytes.NewReader(in), &buf, io.Discard, sizes, order, Options{LimitToInput: true}); err != nil {
		t.Fatalf("Complement failed: %v", err)
	}
	want := readComplementParity(t, "t9b_L.expected.bed")
	if !bytes.Equal(buf.Bytes(), want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// complement.t10 — same input as t9b but WITHOUT -L: emit gaps on chr1 plus
// the full chr2/chr3 (since they have no intervals).
func TestParity_Complement_T10_FullGenomeNoL(t *testing.T) {
	got := runComplementParity(t, "issue_503.bed", "issue_503.genome")
	want := readComplementParity(t, "t10_full.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
