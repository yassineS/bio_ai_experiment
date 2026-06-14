package bedcomplement

// Parity tests against the upstream bedtools complement test suite.
//
// Cases are mirrored from reference_code/bedtools/test/complement/test-complement.sh.
// Each case has a small `.in.bed` and corresponding `.genome` plus an
// `.expected.bed`. The `-L` flag (restrict output to chromosomes seen in the
// input, t9b) and the out-of-range-record warning + clamp case (t9) are both
// implemented and asserted byte-for-byte.

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
	return runComplementParityLimit(t, inputFile, genomeFile, false)
}

func runComplementParityLimit(t *testing.T, inputFile, genomeFile string, limitToInput bool) []byte {
	t.Helper()
	in := readComplementParity(t, inputFile)
	g := readComplementParity(t, genomeFile)
	sizes, order, err := ReadChromSizes(bytes.NewReader(g))
	if err != nil {
		t.Fatalf("ReadChromSizes %s: %v", genomeFile, err)
	}
	var buf bytes.Buffer
	if _, err := Complement(bytes.NewReader(in), &buf, io.Discard, sizes, order, limitToInput); err != nil {
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

// complement.t9 — an input interval that exceeds the chromosome length: a
// warning is emitted to stderr and the interval is clamped, then the complement
// is computed. Combined stderr+stdout matches bedtools v2.31.1 byte-for-byte.
func TestParity_Complement_T9_RecordExceedsChrom(t *testing.T) {
	sizes, order, err := ReadChromSizes(bytes.NewReader([]byte("chr1\t100\n")))
	if err != nil {
		t.Fatalf("ReadChromSizes: %v", err)
	}
	// Upstream's `&> obs` captures stderr then stdout; our warn writes during
	// ingestion and out flushes afterward, so a shared buffer preserves that
	// order.
	var combined bytes.Buffer
	if _, err := Complement(bytes.NewReader([]byte("chr1\t90\t110\n")), &combined, &combined, sizes, order, false); err != nil {
		t.Fatalf("Complement: %v", err)
	}
	want := "***** WARNING: chr1:90-110 exceeds the length of chromosome (chr1)\nchr1\t0\t90\n"
	if combined.String() != want {
		t.Fatalf("mismatch.\nwant:\n%q\ngot:\n%q", want, combined.String())
	}
}

// complement.t9b / t10 (script duplicates 'complement.t9' label) — issue #503,
// the -L flag limits output to chromosomes that had records in the input.
// With -L only chr1 (the chromosome with intervals) is emitted.
func TestParity_Complement_T9b_DashLLimit(t *testing.T) {
	got := runComplementParityLimit(t, "issue_503.bed", "issue_503.genome", true)
	want := readComplementParity(t, "t9b_limit.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
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
