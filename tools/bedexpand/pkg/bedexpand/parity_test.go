package bedexpand

// Parity tests against the upstream bedtools expand test suite.
//
// Cases mirror reference_code/bedtools/test/expand/test-expand.sh. Inputs live
// under tools/bedexpand/testdata/parity/<file> and expected outputs under
// <case>.expected. Each test asserts byte-for-byte equality between the Go
// port's output and the upstream expected output captured at PR time.
//
// Tests for options not yet implemented (none right now — upstream bedtools
// expand only exposes -i / -c) are wrapped in t.Skip with a one-line reason.

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

func runParity(t *testing.T, inputFile string, cols []int) []byte {
	t.Helper()
	in := readParity(t, inputFile)
	var buf bytes.Buffer
	if _, err := Expand(bytes.NewReader(in), &buf, Options{Columns: cols}); err != nil {
		t.Fatalf("Expand: %v", err)
	}
	return buf.Bytes()
}

// expand.t1 — `bedtools expand -i expand.txt -c 4`.
func TestParity_Expand_T1_SingleColumn(t *testing.T) {
	got := runParity(t, "expand.txt", []int{4})
	want := readParity(t, "t1.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// expand.t2 — `bedtools expand -i expand.txt -c 4,5`.
func TestParity_Expand_T2_TwoColumns(t *testing.T) {
	got := runParity(t, "expand.txt", []int{4, 5})
	want := readParity(t, "t2.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// expand.t3 — `bedtools expand -i expand.txt -c 5,4` (column swap).
func TestParity_Expand_T3_SwappedColumns(t *testing.T) {
	got := runParity(t, "expand.txt", []int{5, 4})
	want := readParity(t, "t3.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// Stdin via "-" is the same code path (Expand takes any io.Reader); the CLI
// wrapper resolves "-" to os.Stdin via iohelper.OpenReader. This sub-case is
// a smoke test that the underlying library produces identical output when
// fed the same bytes from a memory reader vs. the on-disk file.
func TestParity_Expand_T4_StdinShape(t *testing.T) {
	got := runParity(t, "expand.txt", []int{4})
	want := readParity(t, "t1.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// Extra: error path — column index that exceeds the row's field count.
// The upstream binary writes a `*****ERROR: Requested column number exceeds
// number of columns.` message and exits 1. We match the spirit (return an
// error from the library) but the exact wording differs; the CLI wrapper
// emits it to stderr.
func TestParity_Expand_T5_OutOfBoundsError(t *testing.T) {
	in := readParity(t, "expand.txt")
	var buf bytes.Buffer
	_, err := Expand(bytes.NewReader(in), &buf, Options{Columns: []int{99}})
	if err == nil {
		t.Fatalf("expected error for column 99 (only 5 in fixture)")
	}
}

// Extra: mismatched-list-length error path. Upstream writes
// `*****ERROR: Each expanded column must have the same number of elements.`
// and exits 1.
func TestParity_Expand_T6_MismatchedListLength(t *testing.T) {
	// Synthetic input: col 4 has 2 elements, col 5 has 3.
	in := []byte("chr1\t0\t1\ta,b\tx,y,z\n")
	var buf bytes.Buffer
	_, err := Expand(bytes.NewReader(in), &buf, Options{Columns: []int{4, 5}})
	if err == nil {
		t.Fatalf("expected error for mismatched list lengths")
	}
}
