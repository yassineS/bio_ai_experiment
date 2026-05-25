package bedclosest

// Parity tests against the upstream bedtools closest test suite.
//
// Cases are mirrored from reference_code/bedtools/test/closest/test-closest.sh.
// Inputs and expected outputs live under tools/bedclosest/testdata/parity/.
// Tests for upstream features bedclosest does not implement (notably the
// strand filters -s/-S, -N "force different names", -k k-nearest, multi-DB
// input with -names/-filenames) are wrapped in t.Skip with a one-line
// rationale.
//
// Note: bedclosest matches the upstream `(b.start - a.end) + 1` distance
// formula, so two touching records report distance 1 (not 0).

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"
)

func readClosestParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runClosestParity(t *testing.T, aFile, bFile string, opts Options) []byte {
	t.Helper()
	a := readClosestParity(t, aFile)
	b := readClosestParity(t, bFile)
	var out bytes.Buffer
	if _, err := Closest(bytes.NewReader(a), bytes.NewReader(b), &out, opts); err != nil {
		t.Fatalf("Closest failed: %v", err)
	}
	return out.Bytes()
}

func runClosestMultiParity(t *testing.T, aFile string, bFiles, labels []string, opts Options) []byte {
	t.Helper()
	a := readClosestParity(t, aFile)
	readers := make([]io.Reader, len(bFiles))
	for i, f := range bFiles {
		readers[i] = bytes.NewReader(readClosestParity(t, f))
	}
	var out bytes.Buffer
	if _, err := ClosestMulti(bytes.NewReader(a), readers, &out, opts, labels); err != nil {
		t.Fatalf("ClosestMulti failed: %v", err)
	}
	return out.Bytes()
}

// closest.t1 — 1bp apart, off-by-one check; A upstream of B.
func TestParity_Closest_T1_OneBpApart(t *testing.T) {
	got := runClosestParity(t, "a.bed", "b.bed", Options{PrintDistance: true, DistanceMode: DistanceAbsolute})
	want := readClosestParity(t, "t1.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// closest.t2 — reciprocal of t1.
func TestParity_Closest_T2_OneBpApartReverse(t *testing.T) {
	got := runClosestParity(t, "b.bed", "a.bed", Options{PrintDistance: true, DistanceMode: DistanceAbsolute})
	want := readClosestParity(t, "t2.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// closest.t3 — 0bp gap (touching by 1 base of overlap): distance is 0 only
// when intervals truly overlap on a 0-based half-open scale.
func TestParity_Closest_T3_OneBaseOverlap(t *testing.T) {
	got := runClosestParity(t, "a.bed", "b-one-bp-closer.bed", Options{PrintDistance: true, DistanceMode: DistanceAbsolute})
	want := readClosestParity(t, "t3.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// closest.t4 — reciprocal of t3.
func TestParity_Closest_T4_OneBaseOverlapReverse(t *testing.T) {
	got := runClosestParity(t, "b-one-bp-closer.bed", "a.bed", Options{PrintDistance: true, DistanceMode: DistanceAbsolute})
	want := readClosestParity(t, "t4.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// closest.t5 — BED4 input (with name column); distance column appended.
func TestParity_Closest_T5_Named(t *testing.T) {
	got := runClosestParity(t, "a.names.bed", "b.names.bed", Options{PrintDistance: true, DistanceMode: DistanceAbsolute})
	want := readClosestParity(t, "t5.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// closest.t6 — `-N` forces the closest B to have a *different* name than A.
func TestParity_Closest_T6_DifferentNames(t *testing.T) {
	got := runClosestParity(t, "a.names.bed", "b.names.bed", Options{PrintDistance: true, DistanceMode: DistanceAbsolute, RequireDifferentNames: true})
	want := readClosestParity(t, "t6_N.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// closest.t7 — `-s` same-strand filter; A's strand is + but B's strand is -,
// so no B is eligible and the BED6 null-placeholder is emitted.
func TestParity_Closest_T7_SameStrand(t *testing.T) {
	got := runClosestParity(t, "strand-test-a.bed", "strand-test-b.bed", Options{SameStrand: true})
	want := readClosestParity(t, "t7_s.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// closest.t8 — `-S` opposite-strand filter; A's strand is + and the lone B has
// strand -, so it is reported.
func TestParity_Closest_T8_OppositeStrand(t *testing.T) {
	got := runClosestParity(t, "strand-test-a.bed", "strand-test-b.bed", Options{OppositeStrand: true})
	want := readClosestParity(t, "t8_S.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// closest.t9 — report ALL overlapping features when ties (default tie mode).
func TestParity_Closest_T9_TiesAll(t *testing.T) {
	got := runClosestParity(t, "close-a.bed", "close-b.bed", Options{TieBreak: TieAll})
	want := readClosestParity(t, "t9_all.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// closest.t10 — `-t first`: report only the first tied B.
func TestParity_Closest_T10_TiesFirst(t *testing.T) {
	got := runClosestParity(t, "close-a.bed", "close-b.bed", Options{TieBreak: TieFirst})
	want := readClosestParity(t, "t10_first.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// closest.t11 — `-t last`: report only the last tied B.
func TestParity_Closest_T11_TiesLast(t *testing.T) {
	got := runClosestParity(t, "close-a.bed", "close-b.bed", Options{TieBreak: TieLast})
	want := readClosestParity(t, "t11_last.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// closest.t13..t15 — multiple databases (`-b A B C`, `-names ...`,
// `-filenames`). These exercise the ClosestMulti entry point.
func TestParity_Closest_T13_MultipleDatabases(t *testing.T) {
	got := runClosestMultiParity(t, "mq1.bed",
		[]string{"mdb1.bed", "mdb2.bed", "mdb3.bed"},
		[]string{"1", "2", "3"},
		Options{})
	want := readClosestParity(t, "t13_multi.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestParity_Closest_T14_DBNames(t *testing.T) {
	got := runClosestMultiParity(t, "mq1.bed",
		[]string{"mdb1.bed", "mdb2.bed", "mdb3.bed"},
		[]string{"a", "b", "c"},
		Options{})
	want := readClosestParity(t, "t14_names.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestParity_Closest_T15_DBFilenames(t *testing.T) {
	got := runClosestMultiParity(t, "mq1.bed",
		[]string{"mdb1.bed", "mdb2.bed", "mdb3.bed"},
		[]string{"mdb1.bed", "mdb2.bed", "mdb3.bed"},
		Options{})
	want := readClosestParity(t, "t15_filenames.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
