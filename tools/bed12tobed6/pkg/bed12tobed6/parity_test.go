package bed12tobed6

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// runParity invokes Convert with the given options on the named input fixture
// and compares the output byte-for-byte against the expected file.
func runParity(t *testing.T, inputFile, expectedFile string, opts Options) {
	t.Helper()
	in, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", inputFile))
	if err != nil {
		t.Fatalf("read input: %v", err)
	}
	want, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", expectedFile))
	if err != nil {
		t.Fatalf("read expected: %v", err)
	}
	var got bytes.Buffer
	if _, err := Convert(bytes.NewReader(in), &got, opts); err != nil {
		t.Fatalf("Convert: %v", err)
	}
	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("parity mismatch:\nwant:\n%s\ngot:\n%s", want, got.Bytes())
	}
}

// TestParity_Bed12ToBed6_T1_OneBlock mirrors upstream test-bed12tobed6.sh t1.
func TestParity_Bed12ToBed6_T1_OneBlock(t *testing.T) {
	runParity(t, "one_blocks.bed", "t1.expected.bed", Options{})
}

// TestParity_Bed12ToBed6_T2_TwoBlocks mirrors upstream t2.
func TestParity_Bed12ToBed6_T2_TwoBlocks(t *testing.T) {
	runParity(t, "two_blocks.bed", "t2.expected.bed", Options{})
}

// TestParity_Bed12ToBed6_T3_ThreeBlocks mirrors upstream t3.
func TestParity_Bed12ToBed6_T3_ThreeBlocks(t *testing.T) {
	runParity(t, "three_blocks.bed", "t3.expected.bed", Options{})
}

// TestParity_Bed12ToBed6_T4_ThreeBlocksNumbered mirrors upstream t4 (-n).
func TestParity_Bed12ToBed6_T4_ThreeBlocksNumbered(t *testing.T) {
	runParity(t, "three_blocks.bed", "t4.expected.bed", Options{NumberBlocks: true})
}

// TestParity_Bed12ToBed6_T5_ReverseStrandNumbered mirrors upstream t5 — same
// record but on the '-' strand, with -n, which causes upstream to reverse the
// block numbering.
func TestParity_Bed12ToBed6_T5_ReverseStrandNumbered(t *testing.T) {
	runParity(t, "t5.input.bed", "t5.expected.bed", Options{NumberBlocks: true})
}
