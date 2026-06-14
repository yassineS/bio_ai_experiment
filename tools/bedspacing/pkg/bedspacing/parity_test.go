package bedspacing

// Parity tests against the upstream bedtools spacing test suite.
//
// Cases mirror reference_code/bedtools/test/spacing/test-spacing.sh. The
// upstream suite ships only one inline test (spacing.t01); the rest of the
// cases here exercise the documented edge cases that the upstream
// implementation in src/spacingFile/spacingFile.cpp covers (per-chrom
// reset, exact abut, overlap, first record).

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

func runParity(t *testing.T, inputFile string) []byte {
	t.Helper()
	in := readParity(t, inputFile)
	var buf bytes.Buffer
	if _, err := Spacing(bytes.NewReader(in), &buf); err != nil {
		t.Fatalf("Spacing: %v", err)
	}
	return buf.Bytes()
}

// spacing.t01 — the canonical upstream test: covers ".", "-1", "0", N>0 and
// the per-chromosome reset.
func TestParity_Spacing_T01_Basic(t *testing.T) {
	got := runParity(t, "a.bed")
	want := readParity(t, "t01.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// Synthetic t02: two-chromosome file where each chrom's first record is ".".
func TestParity_Spacing_T02_PerChromReset(t *testing.T) {
	in := []byte("chr1\t0\t10\nchr1\t20\t30\nchr2\t100\t200\nchr2\t250\t300\n")
	want := []byte("chr1\t0\t10\t.\nchr1\t20\t30\t10\nchr2\t100\t200\t.\nchr2\t250\t300\t50\n")
	var buf bytes.Buffer
	if _, err := Spacing(bytes.NewReader(in), &buf); err != nil {
		t.Fatalf("Spacing: %v", err)
	}
	if !bytes.Equal(buf.Bytes(), want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// Synthetic t03: exact abut emits "0".
func TestParity_Spacing_T03_ExactAbut(t *testing.T) {
	in := []byte("chr1\t0\t10\nchr1\t10\t20\n")
	want := []byte("chr1\t0\t10\t.\nchr1\t10\t20\t0\n")
	var buf bytes.Buffer
	if _, err := Spacing(bytes.NewReader(in), &buf); err != nil {
		t.Fatalf("Spacing: %v", err)
	}
	if !bytes.Equal(buf.Bytes(), want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// Synthetic t04: overlap emits "-1".
func TestParity_Spacing_T04_Overlap(t *testing.T) {
	in := []byte("chr1\t0\t30\nchr1\t10\t40\n")
	want := []byte("chr1\t0\t30\t.\nchr1\t10\t40\t-1\n")
	var buf bytes.Buffer
	if _, err := Spacing(bytes.NewReader(in), &buf); err != nil {
		t.Fatalf("Spacing: %v", err)
	}
	if !bytes.Equal(buf.Bytes(), want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// Synthetic t05: BED6 input — extra columns preserved verbatim.
func TestParity_Spacing_T05_BED6Preserved(t *testing.T) {
	in := []byte("chr1\t0\t10\tfoo\t100\t+\nchr1\t20\t30\tbar\t200\t-\n")
	want := []byte("chr1\t0\t10\tfoo\t100\t+\t.\nchr1\t20\t30\tbar\t200\t-\t10\n")
	var buf bytes.Buffer
	if _, err := Spacing(bytes.NewReader(in), &buf); err != nil {
		t.Fatalf("Spacing: %v", err)
	}
	if !bytes.Equal(buf.Bytes(), want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// Synthetic t06: single record file.
func TestParity_Spacing_T06_SingleRecord(t *testing.T) {
	in := []byte("chr1\t0\t10\n")
	want := []byte("chr1\t0\t10\t.\n")
	var buf bytes.Buffer
	if _, err := Spacing(bytes.NewReader(in), &buf); err != nil {
		t.Fatalf("Spacing: %v", err)
	}
	if !bytes.Equal(buf.Bytes(), want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// spacing on BAM input — parity with `bedtools spacing -i spacing.bam -bed`.
//
// Upstream `bedtools spacing` accepts BAM (its usage line reads
// "-i <bed/gff/vcf/bam>"); with -bed it prints each alignment as a BED12 line
// followed by the spacing column. The vendored spacing.bam (built from
// spacing.sam) and spacing_bam.expected were produced by the bundled
// reference_code/bedtools/bin/bedtools v2.31.1; this test asserts our output
// matches it byte-for-byte. The fixture exercises ".", "0", "-1", a positive
// gap, a spliced (N-CIGAR, multi-block) alignment, and the per-chromosome
// reset.
func TestParity_Spacing_T07_BAMInput(t *testing.T) {
	in := readParity(t, "spacing.bam")
	want := readParity(t, "spacing_bam.expected")
	var buf bytes.Buffer
	if _, err := Spacing(bytes.NewReader(in), &buf); err != nil {
		t.Fatalf("Spacing on BAM: %v", err)
	}
	if !bytes.Equal(buf.Bytes(), want) {
		t.Fatalf("BAM spacing mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// SAM-text input is auto-detected (leading '@' header) and routed through the
// same alignment-to-BED12 conversion as BAM. The expected output is identical
// to the BAM fixture because spacing.bam was built from spacing.sam.
func TestSpacing_SAMTextInput(t *testing.T) {
	in := readParity(t, "spacing.sam")
	want := readParity(t, "spacing_bam.expected")
	var buf bytes.Buffer
	if _, err := Spacing(bytes.NewReader(in), &buf); err != nil {
		t.Fatalf("Spacing on SAM: %v", err)
	}
	if !bytes.Equal(buf.Bytes(), want) {
		t.Fatalf("SAM spacing mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}
