package bedsubtract

// Parity tests against the upstream bedtools subtract test suite.
//
// Cases are mirrored from reference_code/bedtools/test/subtract/test-subtract.sh.
// Inputs and expected outputs live under tools/bedsubtract/testdata/parity/.
// Tests for upstream features bedsubtract does not implement (notably `-N`,
// which sums coverage across all B overlaps and drops A as a whole) are
// wrapped in t.Skip with a one-line rationale.

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"
)

func readSubtractParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runSubtractParity(t *testing.T, aFile string, bFiles []string, opts Options) []byte {
	t.Helper()
	a := readSubtractParity(t, aFile)
	// Concatenate B files (mirrors upstream's "-b b.bed b2.bed" multi-DB form).
	var bAll bytes.Buffer
	for _, f := range bFiles {
		bAll.Write(readSubtractParity(t, f))
	}
	var out bytes.Buffer
	if _, err := Subtract(bytes.NewReader(a), bytes.NewReader(bAll.Bytes()), &out, opts); err != nil {
		t.Fatalf("Subtract failed: %v", err)
	}
	// Silence unused import warning when no skipped helpers reference io.
	_ = io.Discard
	return out.Bytes()
}

// subtract.t1 — baseline subtraction, no flags.
func TestParity_Subtract_T1_Baseline(t *testing.T) {
	got := runSubtractParity(t, "a.bed", []string{"b.bed"}, Options{})
	want := readSubtractParity(t, "t1_baseline.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// subtract.t2 — -f 0.5: overlap fraction below threshold ⇒ no subtraction.
func TestParity_Subtract_T2_F05(t *testing.T) {
	got := runSubtractParity(t, "a.bed", []string{"b.bed"}, Options{MinFraction: 0.5})
	want := readSubtractParity(t, "t2_f05.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// subtract.t3 — -f 0.1: threshold met, subtract.
func TestParity_Subtract_T3_F01(t *testing.T) {
	got := runSubtractParity(t, "a.bed", []string{"b.bed"}, Options{MinFraction: 0.1})
	want := readSubtractParity(t, "t3_f01.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// subtract.t4 — -s (same strand): A and B differ on strand ⇒ no subtraction.
func TestParity_Subtract_T4_SameStrand(t *testing.T) {
	got := runSubtractParity(t, "a.bed", []string{"b.bed"}, Options{SameStrand: true})
	want := readSubtractParity(t, "t4_s.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// subtract.t5 — -S (opposite strand): subtract since A=+, B=-.
func TestParity_Subtract_T5_OppositeStrand(t *testing.T) {
	got := runSubtractParity(t, "a.bed", []string{"b.bed"}, Options{OppositeStrand: true})
	want := readSubtractParity(t, "t5_S.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// subtract.t6 — -A removes any A interval that overlaps B at all.
func TestParity_Subtract_T6_RemoveEntire(t *testing.T) {
	got := runSubtractParity(t, "a.bed", []string{"b.bed"}, Options{RemoveEntire: true})
	want := readSubtractParity(t, "t6_A.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// subtract.t7 — -A -f 0.5: threshold blocks the per-B match, A1 survives intact.
func TestParity_Subtract_T7_RemoveEntireF05(t *testing.T) {
	got := runSubtractParity(t, "a.bed", []string{"b.bed"}, Options{RemoveEntire: true, MinFraction: 0.5})
	want := readSubtractParity(t, "t7_A_f05.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// subtract.t8 — -A -f 0.1: threshold met, A1 dropped entirely.
func TestParity_Subtract_T8_RemoveEntireF01(t *testing.T) {
	got := runSubtractParity(t, "a.bed", []string{"b.bed"}, Options{RemoveEntire: true, MinFraction: 0.1})
	want := readSubtractParity(t, "t8_A_f01.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// subtract.t9 — -N sums coverage across all B overlaps and drops A iff that
// union exceeds the fraction. With -f 0.4 the union exactly equals the
// threshold (4/10), the `>` test fails, so A is emitted intact.
func TestParity_Subtract_T9_DashN(t *testing.T) {
	got := runSubtractParity(t, "c.bed", []string{"d.bed"}, Options{UnionDrop: true, MinFraction: 0.4})
	want := readSubtractParity(t, "t9_N_f04.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// subtract.t10 — -f 0.39 just barely under the union coverage (0.4 > 0.39)
// ⇒ A is dropped entirely.
func TestParity_Subtract_T10_DashNStricter(t *testing.T) {
	got := runSubtractParity(t, "c.bed", []string{"d.bed"}, Options{UnionDrop: true, MinFraction: 0.39})
	want := readSubtractParity(t, "t10_N_f039.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// subtract.t19 — multiple -b databases (b.bed AND b2.bed) concatenated.
func TestParity_Subtract_T19_TwoDatabases(t *testing.T) {
	got := runSubtractParity(t, "a.bed", []string{"b.bed", "b2.bed"}, Options{})
	want := readSubtractParity(t, "t19_two_db.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// subtract.t20 — multiple -b databases with -f 0.8 too tight ⇒ no subtraction.
func TestParity_Subtract_T20_TwoDatabasesF08(t *testing.T) {
	got := runSubtractParity(t, "a.bed", []string{"b.bed", "b2.bed"}, Options{MinFraction: 0.8})
	want := readSubtractParity(t, "t20_two_db_f08.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// subtract.t21 — -f 0.6 -N with two DBs. a1 has union coverage 7/10 (0.7 > 0.6
// ⇒ drop); a2 has union coverage 10/20 = 0.5 (0.5 > 0.6 false ⇒ keep).
func TestParity_Subtract_T21_TwoDBsNFraction(t *testing.T) {
	got := runSubtractParity(t, "a.bed", []string{"b.bed", "b2.bed"},
		Options{UnionDrop: true, MinFraction: 0.6})
	want := readSubtractParity(t, "t21_two_db_N_f06.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
