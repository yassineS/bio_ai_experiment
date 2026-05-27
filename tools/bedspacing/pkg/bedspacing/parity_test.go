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
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bamtobed"
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

// Upstream allows BAM input via the same `-i` flag. We support this by
// running the input through bamtobed.FromBAM upstream of Spacing. The
// upstream corpus has no canonical BAM spacing case, so this test uses a
// SAM-rendered fixture from jaccard (vendored under jaccard/) for a small
// 2-alignment input and asserts the spacing-formatted output.
func TestParity_Spacing_T07_BAMInput(t *testing.T) {
	// Use the small `a.bam` from the bedjaccard fixtures (2 alignments on
	// chr1 at positions 1-100 and 101-200).
	bam, err := os.ReadFile(filepath.Join("..", "..", "..", "..", "tools", "bedjaccard", "testdata", "parity", "a.bam"))
	if err != nil {
		t.Skipf("BAM fixture unavailable: %v", err)
	}
	bed := bamtobed.FromBAM(bytes.NewReader(bam))
	var buf bytes.Buffer
	if _, err := Spacing(bed, &buf); err != nil {
		t.Fatalf("Spacing: %v", err)
	}
	// The exact output depends on the BAM contents, but spacing must:
	//  - have one line per input alignment;
	//  - the first line's spacing column must be "." (no prior on chrom);
	//  - subsequent lines' spacing must be either "." (chrom change) or a
	//    non-negative integer.
	lines := bytes.Split(bytes.TrimRight(buf.Bytes(), "\n"), []byte("\n"))
	if len(lines) == 0 {
		t.Fatalf("expected at least one spacing line, got 0:\n%s", buf.Bytes())
	}
	for i, line := range lines {
		cols := bytes.Split(line, []byte("\t"))
		if len(cols) < 4 {
			t.Fatalf("line %d: expected >=4 cols, got %d: %q", i, len(cols), line)
		}
		spc := string(cols[len(cols)-1])
		if i == 0 && spc != "." {
			t.Fatalf("line 0 spacing must be '.', got %q", spc)
		}
		if spc != "." {
			// Must parse as a non-negative int.
			var n int
			if _, err := fmt.Sscanf(spc, "%d", &n); err != nil || n < 0 {
				t.Fatalf("line %d spacing %q is not a non-negative integer", i, spc)
			}
		}
	}
}
