package bedintersect

// Parity tests against the upstream bedtools intersect test suite.
//
// Cases are mirrored from reference_code/bedtools/test/intersect/test-intersect.sh.
// Inputs and expected outputs live under tools/bedintersect/testdata/parity/.
// bedintersect implements the most common subset of `bedtools intersect`
// options; tests that exercise -wo, -wao, -wa+-wb combined output, BAM/VCF
// input, or `-split` are wrapped in t.Skip with a one-line rationale.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

func readIntersectParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runIntersectParity(t *testing.T, aFile, bFile string, opts IntersectOptions) []byte {
	t.Helper()
	a := readIntersectParity(t, aFile)
	b := readIntersectParity(t, bFile)
	var out bytes.Buffer
	if _, err := Intersect(bytes.NewReader(a), bytes.NewReader(b), &out, opts); err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}
	return out.Bytes()
}

// intersect.t01 — basic self intersection: each A overlaps with itself.
func TestParity_Intersect_T01_SelfIntersect(t *testing.T) {
	got := runIntersectParity(t, "a.bed", "a.bed", IntersectOptions{})
	want := readIntersectParity(t, "t01_self.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// intersect.t02 — `-v` self intersect: expected empty.
func TestParity_Intersect_T02_VSelf(t *testing.T) {
	got := runIntersectParity(t, "a.bed", "a.bed", IntersectOptions{NoOverlap: true})
	want := readIntersectParity(t, "t02_v.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%q\ngot:\n%q", want, got)
	}
}

// intersect.t03 — `-c` appends per-A overlap counts.
func TestParity_Intersect_T03_Count(t *testing.T) {
	got := runIntersectParity(t, "a.bed", "b.bed", IntersectOptions{Count: true})
	want := readIntersectParity(t, "t03_c.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// intersect.t04 — `-c -s`: strand-aware count.
func TestParity_Intersect_T04_CountStrand(t *testing.T) {
	got := runIntersectParity(t, "a.bed", "b.bed", IntersectOptions{Count: true, StrandSpec: true})
	want := readIntersectParity(t, "t04_cs.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// intersect.t05 — `-c -s -f 0.1`: strand + fraction-of-A filter.
func TestParity_Intersect_T05_CountStrandFractionA(t *testing.T) {
	got := runIntersectParity(t, "a.bed", "b.bed", IntersectOptions{Count: true, StrandSpec: true, FractionA: 0.1})
	want := readIntersectParity(t, "t05_csf01.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// intersect.t06 — default a vs b intersection: emit clipped A for each hit.
func TestParity_Intersect_T06_Default(t *testing.T) {
	got := runIntersectParity(t, "a.bed", "b.bed", IntersectOptions{})
	want := readIntersectParity(t, "t06_default.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// intersect.t07 — `-wa`: emit the full A record per hit.
func TestParity_Intersect_T07_WriteA(t *testing.T) {
	got := runIntersectParity(t, "a.bed", "b.bed", IntersectOptions{WriteA: true})
	want := readIntersectParity(t, "t07_wa.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// intersect.t08 — `-wa -wb`: emit A and B side-by-side. Not implemented:
// bedintersect's WriteA and WriteB toggles are mutually exclusive, so we
// cannot emit both fields together.
func TestParity_Intersect_T08_WAWB(t *testing.T) {
	t.Skip("unimplemented: -wa with -wb combined output (side-by-side A and B columns)")
}

// intersect.t09 — `-wo`: emit A, B and the overlap length per hit.
func TestParity_Intersect_T09_WO(t *testing.T) {
	t.Skip("unimplemented: -wo (write overlap length alongside A and B)")
}

// intersect.t10 — `-wao`: like -wo but emit A even when there are no hits.
func TestParity_Intersect_T10_WAO(t *testing.T) {
	t.Skip("unimplemented: -wao (write all + overlap length)")
}

// intersect.t11 / t12 — `-wo -s` / `-wao -s`: same writer features as t09/t10.
func TestParity_Intersect_T11_WOStrand(t *testing.T) {
	t.Skip("unimplemented: -wo combined writer (and -s strand filter on it)")
}
func TestParity_Intersect_T12_WAOStrand(t *testing.T) {
	t.Skip("unimplemented: -wao combined writer (and -s strand filter on it)")
}

// intersect.t13..t16 — `-a -`, `-a stdin`, `-b -`, `-b stdin`: these are
// CLI-level stdin redirections; the parity layer drives the library API
// directly which always takes io.Readers, so the behavior is covered by t06.
func TestParity_Intersect_T13_StdinAVariants(t *testing.T) {
	got := runIntersectParity(t, "a.bed", "b.bed", IntersectOptions{})
	want := readIntersectParity(t, "t06_default.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
