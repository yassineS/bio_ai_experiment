package bedigv

// Parity tests for `bedtools igv`.
//
// Upstream ships no `igv/` test subdir under `reference_code/bedtools/test/`,
// so expected outputs are derived directly from the upstream source
// (`reference_code/bedtools/src/bedToIgv/bedToIgv.cpp`) by mechanically
// playing forward the ProcessBed() routine on the fixture below. Each test
// asserts byte-for-byte equality between this Go port's stdout and the
// fixture's hand-derived expected file.
//
// Fixture: testdata/parity/igv.bed (three BED6 records). Cases:
//
//   - t1: defaults  -- `bedtools igv -i igv.bed`.
//   - t2: `-path /tmp/snaps -sess my.xml -name`.
//   - t3: `-slop 50 -img svg -sort base -clps`.
//
// If the upstream project ever adds an `igv/` test directory, we can drop
// these in favour of the upstream golden files.

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
	if _, err := Run(bytes.NewReader(in), &buf, opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	return buf.Bytes()
}

// t1: defaults — `bedtools igv -i igv.bed`.
func TestParity_Igv_T1_Defaults(t *testing.T) {
	got := runParity(t, "igv.bed", Options{})
	want := readParity(t, "t1.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// t2: `bedtools igv -i igv.bed -path /tmp/snaps -sess my.xml -name`.
func TestParity_Igv_T2_PathSessionName(t *testing.T) {
	got := runParity(t, "igv.bed", Options{
		Path:     "/tmp/snaps",
		Session:  "my.xml",
		UseNames: true,
	})
	want := readParity(t, "t2.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// t3: `bedtools igv -i igv.bed -slop 50 -img svg -sort base -clps`.
func TestParity_Igv_T3_SlopSortCollapseImg(t *testing.T) {
	got := runParity(t, "igv.bed", Options{
		Slop:      50,
		Sort:      SortBase,
		Collapse:  true,
		ImageType: ImageSVG,
	})
	want := readParity(t, "t3.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
