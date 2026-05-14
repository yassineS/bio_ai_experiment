package bedflank

// Parity tests against the upstream bedtools flank test suite.
//
// Cases are mirrored from reference_code/bedtools/test/flank/test-flank.sh.
// Inputs and expected outputs live under tools/bedflank/testdata/parity/.
// All upstream flank test cases (t1..t11) are implementable here because
// bedflank supports all of the options exercised by that test script.

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"
)

func readFlankParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runFlankParity(t *testing.T, inputFile, genomeFile string, opts Options) []byte {
	t.Helper()
	in := readFlankParity(t, inputFile)
	g := readFlankParity(t, genomeFile)
	sizes, err := ReadChromSizes(bytes.NewReader(g))
	if err != nil {
		t.Fatalf("ReadChromSizes %s: %v", genomeFile, err)
	}
	var out bytes.Buffer
	if _, err := Flank(bytes.NewReader(in), &out, io.Discard, sizes, opts); err != nil {
		t.Fatalf("Flank failed: %v", err)
	}
	return out.Bytes()
}

// flank.t1 — -b 5: both flanks symmetric.
func TestParity_Flank_T1_B5(t *testing.T) {
	got := runFlankParity(t, "a.bed", "tiny.genome", Options{Both: true, BothAdd: 5})
	want := readFlankParity(t, "t1_b5.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// flank.t2 — -l 5 -r 5: equivalent to -b 5.
func TestParity_Flank_T2_L5R5(t *testing.T) {
	got := runFlankParity(t, "a.bed", "tiny.genome", Options{LeftAdd: 5, RightAdd: 5})
	want := readFlankParity(t, "t1_b5.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// flank.t3 — -l 5 -r 0: left flank only.
func TestParity_Flank_T3_L5R0(t *testing.T) {
	got := runFlankParity(t, "a.bed", "tiny.genome", Options{LeftAdd: 5, RightAdd: 0})
	want := readFlankParity(t, "t3_l5_r0.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// flank.t4 — -l 0 -r 5: right flank only.
func TestParity_Flank_T4_L0R5(t *testing.T) {
	got := runFlankParity(t, "a.bed", "tiny.genome", Options{LeftAdd: 0, RightAdd: 5})
	want := readFlankParity(t, "t4_l0_r5.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// flank.t5 — -l 5 -r 0 -s: strand-aware. '+' yields left flank, '-' yields right.
func TestParity_Flank_T5_L5R0Strand(t *testing.T) {
	got := runFlankParity(t, "a.bed", "tiny.genome", Options{LeftAdd: 5, RightAdd: 0, StrandSpec: true})
	want := readFlankParity(t, "t5_l5_r0_s.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// flank.t6 — -l 0 -r 5 -s.
func TestParity_Flank_T6_L0R5Strand(t *testing.T) {
	got := runFlankParity(t, "a.bed", "tiny.genome", Options{LeftAdd: 0, RightAdd: 5, StrandSpec: true})
	want := readFlankParity(t, "t6_l0_r5_s.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// flank.t7 — -b 5 -s: symmetric flanks, strand-aware (no observable change on
// symmetric flanks but the code path is exercised).
func TestParity_Flank_T7_B5Strand(t *testing.T) {
	got := runFlankParity(t, "a.bed", "tiny.genome", Options{Both: true, BothAdd: 5, StrandSpec: true})
	want := readFlankParity(t, "t1_b5.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// flank.t8 — -b 200: left flank clipped to chrom start.
func TestParity_Flank_T8_B200(t *testing.T) {
	got := runFlankParity(t, "a.bed", "tiny.genome", Options{Both: true, BothAdd: 200})
	want := readFlankParity(t, "t8_b200.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// flank.t9 — -r 1000 only: right flank clipped to chrom end.
func TestParity_Flank_T9_R1000(t *testing.T) {
	got := runFlankParity(t, "a.bed", "tiny.genome", Options{LeftAdd: 0, RightAdd: 1000})
	want := readFlankParity(t, "t9_r1000.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// flank.t10 — -b 2000: both flanks clipped.
func TestParity_Flank_T10_B2000(t *testing.T) {
	got := runFlankParity(t, "a.bed", "tiny.genome", Options{Both: true, BothAdd: 2000})
	want := readFlankParity(t, "t10_b2000.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// flank.t11 — -b 2000 -s.
func TestParity_Flank_T11_B2000Strand(t *testing.T) {
	got := runFlankParity(t, "a.bed", "tiny.genome", Options{Both: true, BothAdd: 2000, StrandSpec: true})
	want := readFlankParity(t, "t10_b2000.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
