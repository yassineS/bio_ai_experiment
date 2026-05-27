package bedslop

// Parity tests against the upstream bedtools slop test suite.
//
// Cases are mirrored from reference_code/bedtools/test/slop/test-slop.sh.
// Inputs and expected outputs live under tools/bedslop/testdata/parity/.
// Tests that exercise upstream features bedslop does not implement (notably
// reading the giant external human.hg19.genome used by slop.t13/t14) are
// wrapped in t.Skip with a one-line rationale.

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"
)

func readSlopParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runSlopParity(t *testing.T, inputFile, genomeFile string, opts Options) []byte {
	t.Helper()
	in := readSlopParity(t, inputFile)
	g := readSlopParity(t, genomeFile)
	sizes, err := ReadChromSizes(bytes.NewReader(g))
	if err != nil {
		t.Fatalf("ReadChromSizes %s: %v", genomeFile, err)
	}
	var out bytes.Buffer
	if _, err := Slop(bytes.NewReader(in), &out, io.Discard, sizes, opts); err != nil {
		t.Fatalf("Slop failed: %v", err)
	}
	return out.Bytes()
}

// slop.t1 — -b 5: symmetric flanks.
func TestParity_Slop_T1_B5(t *testing.T) {
	got := runSlopParity(t, "a.bed", "tiny.genome", Options{Both: true, BothAdd: 5})
	want := readSlopParity(t, "t1_b5.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t2 — -l 5 -r 5: equivalent to -b 5.
func TestParity_Slop_T2_L5R5(t *testing.T) {
	got := runSlopParity(t, "a.bed", "tiny.genome", Options{LeftAdd: 5, RightAdd: 5})
	want := readSlopParity(t, "t1_b5.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t3 — -l 5 -r 0.
func TestParity_Slop_T3_L5R0(t *testing.T) {
	got := runSlopParity(t, "a.bed", "tiny.genome", Options{LeftAdd: 5, RightAdd: 0})
	want := readSlopParity(t, "t3_l5_r0.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t4 — -l 0 -r 5.
func TestParity_Slop_T4_L0R5(t *testing.T) {
	got := runSlopParity(t, "a.bed", "tiny.genome", Options{LeftAdd: 0, RightAdd: 5})
	want := readSlopParity(t, "t4_l0_r5.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t5 — -l 5 -r 0 -s: strand-aware swap on '-' rows.
func TestParity_Slop_T5_L5R0Strand(t *testing.T) {
	got := runSlopParity(t, "a.bed", "tiny.genome", Options{LeftAdd: 5, RightAdd: 0, StrandSpec: true})
	want := readSlopParity(t, "t5_l5_r0_s.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t6 — -l 0 -r 5 -s.
func TestParity_Slop_T6_L0R5Strand(t *testing.T) {
	got := runSlopParity(t, "a.bed", "tiny.genome", Options{LeftAdd: 0, RightAdd: 5, StrandSpec: true})
	want := readSlopParity(t, "t6_l0_r5_s.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t8 — slop past chrom start; clip to 0.
func TestParity_Slop_T8_PastStart(t *testing.T) {
	got := runSlopParity(t, "a.bed", "tiny.genome", Options{Both: true, BothAdd: 200})
	want := readSlopParity(t, "t8_b200.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t9 — slop past chrom end; clip to chrom size.
func TestParity_Slop_T9_PastEnd(t *testing.T) {
	got := runSlopParity(t, "a.bed", "tiny.genome", Options{LeftAdd: 0, RightAdd: 1000})
	want := readSlopParity(t, "t9_r1000.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t10 — slop past both ends; clip to [0, chromSize].
func TestParity_Slop_T10_PastBoth(t *testing.T) {
	got := runSlopParity(t, "a.bed", "tiny.genome", Options{Both: true, BothAdd: 2000})
	want := readSlopParity(t, "t10_b2000.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t13 / t14 — float-point precision regression tests using
// human.hg19.genome (chr1 length 249,250,621). Now vendored under
// tools/bedslop/testdata/parity/human.hg19.genome.
func TestParity_Slop_T13_FloatPrecision(t *testing.T) {
	in := []byte("chr1\t16778271\t16778571\n")
	g := readSlopParity(t, "human.hg19.genome")
	sizes, err := ReadChromSizes(bytes.NewReader(g))
	if err != nil {
		t.Fatalf("ReadChromSizes: %v", err)
	}
	var out bytes.Buffer
	if _, err := Slop(bytes.NewReader(in), &out, io.Discard, sizes, Options{LeftAdd: 200, RightAdd: 200}); err != nil {
		t.Fatalf("Slop failed: %v", err)
	}
	want := []byte("chr1\t16778071\t16778771\n")
	if !bytes.Equal(out.Bytes(), want) {
		t.Fatalf("mismatch.\nwant: %q\ngot:  %q", want, out.Bytes())
	}
}
func TestParity_Slop_T14_FloatPrecisionB(t *testing.T) {
	in := []byte("chr1\t16778272\t16778572\n")
	g := readSlopParity(t, "human.hg19.genome")
	sizes, err := ReadChromSizes(bytes.NewReader(g))
	if err != nil {
		t.Fatalf("ReadChromSizes: %v", err)
	}
	var out bytes.Buffer
	if _, err := Slop(bytes.NewReader(in), &out, io.Discard, sizes, Options{LeftAdd: 200, RightAdd: 200}); err != nil {
		t.Fatalf("Slop failed: %v", err)
	}
	want := []byte("chr1\t16778072\t16778772\n")
	if !bytes.Equal(out.Bytes(), want) {
		t.Fatalf("mismatch.\nwant: %q\ngot:  %q", want, out.Bytes())
	}
}

// slop.t16 — negative -l (no strand): straight subtraction on the left edge.
func TestParity_Slop_T16_NegLeft(t *testing.T) {
	got := runSlopParity(t, "t16_negl_in.bed", "tiny.genome", Options{LeftAdd: -60, RightAdd: 60})
	want := readSlopParity(t, "t16_negl.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t17 — negative -l with -s on a '-' strand record: swap l/r first.
func TestParity_Slop_T17_NegLeftStrand(t *testing.T) {
	got := runSlopParity(t, "t17_in.bed", "tiny.genome", Options{LeftAdd: -60, RightAdd: 60, StrandSpec: true})
	want := readSlopParity(t, "t17_negl_s.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t20 — both flanks negative, large enough to invert: coordinates swap.
func TestParity_Slop_T20_Crossover(t *testing.T) {
	got := runSlopParity(t, "t20_in.bed", "tiny.genome", Options{LeftAdd: -60, RightAdd: -60, StrandSpec: true})
	want := readSlopParity(t, "t20_crossover.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t21 — negative slop pushes interval past chrom end: collapse to 1bp at
// the right boundary.
func TestParity_Slop_T21_EdgeRight(t *testing.T) {
	got := runSlopParity(t, "t21_in.bed", "tiny.genome", Options{LeftAdd: 60, RightAdd: -60, StrandSpec: true})
	want := readSlopParity(t, "t21_edge.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// slop.t22 — negative slop pushes interval below chrom start: collapse to 1bp
// at the left boundary.
func TestParity_Slop_T22_EdgeLeft(t *testing.T) {
	got := runSlopParity(t, "t22_in.bed", "tiny.genome", Options{LeftAdd: -60, RightAdd: 60, StrandSpec: true})
	want := readSlopParity(t, "t22_edge.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
