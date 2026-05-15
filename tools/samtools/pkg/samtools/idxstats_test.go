package samtools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/sam"
)

// makeIndexedBAM is a small helper used by idxstats and other tests:
// given a SAM-text input on disk, sort it (coordinate), write a BAM,
// and build a BAI sibling. Returns the BAM path.
func makeIndexedBAM(t *testing.T, samPath string) string {
	t.Helper()
	dir := t.TempDir()
	bamPath := filepath.Join(dir, "in.bam")

	in, err := os.Open(samPath)
	if err != nil {
		t.Fatalf("open sam: %v", err)
	}
	defer in.Close()

	// Sort + emit BAM in one pipeline using the existing Sort entry point.
	out, err := os.Create(bamPath)
	if err != nil {
		t.Fatalf("create bam: %v", err)
	}
	if err := Sort(in, out, SortOptions{Order: SortCoordinate, OutputBAM: true}); err != nil {
		t.Fatalf("sort: %v", err)
	}
	_ = out.Close()

	// Index the sorted BAM.
	if err := IndexFile(bamPath, "", IndexOptions{}); err != nil {
		t.Fatalf("index: %v", err)
	}
	return bamPath
}

// ---- unit -------------------------------------------------------------

// idxstats hand-computed on basic.sam:
//
//	chr1: 3 mapped (read1, read2, read5 secondary), 0 unmapped
//	chr2: 1 mapped (read4), 0 unmapped
//	*  : 1 unmapped (read3)
func TestIdxstats_BasicHandComputed(t *testing.T) {
	bam := makeIndexedBAM(t, parityPath(t, "basic.sam"))
	rows, err := Idxstats(bam)
	if err != nil {
		t.Fatalf("Idxstats: %v", err)
	}
	want := []IdxstatsRow{
		{Name: "chr1", Length: 1000, Mapped: 3, Unmapped: 0},
		{Name: "chr2", Length: 500, Mapped: 1, Unmapped: 0},
		{Name: "*", Length: 0, Mapped: 0, Unmapped: 1},
	}
	if len(rows) != len(want) {
		t.Fatalf("row count: got %d want %d", len(rows), len(want))
	}
	for i, w := range want {
		if rows[i] != w {
			t.Errorf("row %d: got %+v want %+v", i, rows[i], w)
		}
	}
}

// Empty SAM (no records) should produce one row per @SQ with zeros and
// no trailing unplaced row entries (NoCoor == 0).
func TestIdxstats_EmptyHandComputed(t *testing.T) {
	bam := makeIndexedBAM(t, parityPath(t, "empty.sam"))
	rows, err := Idxstats(bam)
	if err != nil {
		t.Fatalf("Idxstats: %v", err)
	}
	if len(rows) == 0 {
		t.Fatalf("expected at least the unplaced trailer row, got 0")
	}
	last := rows[len(rows)-1]
	if last.Name != "*" || last.Mapped != 0 || last.Unmapped != 0 {
		t.Errorf("unplaced trailer: got %+v", last)
	}
	for _, r := range rows[:len(rows)-1] {
		if r.Mapped != 0 || r.Unmapped != 0 {
			t.Errorf("expected empty counts, got %+v", r)
		}
	}
}

// WriteIdxstats emits the canonical TSV form.
func TestIdxstats_WriteFormat(t *testing.T) {
	rows := []IdxstatsRow{
		{Name: "chr1", Length: 1000, Mapped: 3, Unmapped: 0},
		{Name: "*", Length: 0, Mapped: 0, Unmapped: 5},
	}
	var buf bytes.Buffer
	if err := WriteIdxstats(&buf, rows); err != nil {
		t.Fatalf("WriteIdxstats: %v", err)
	}
	want := "chr1\t1000\t3\t0\n*\t0\t0\t5\n"
	if buf.String() != want {
		t.Errorf("WriteIdxstats:\n got %q\nwant %q", buf.String(), want)
	}
}

// Scan-fallback — when the BAI isn't present, idxstats still works.
func TestIdxstats_ScanFallback(t *testing.T) {
	bam := makeIndexedBAM(t, parityPath(t, "basic.sam"))
	// Remove the .bai so the scan path is used.
	if err := os.Remove(bam + ".bai"); err != nil {
		t.Fatalf("rm bai: %v", err)
	}
	rows, err := Idxstats(bam)
	if err != nil {
		t.Fatalf("Idxstats: %v", err)
	}
	// Same row counts as the indexed test.
	want := []IdxstatsRow{
		{Name: "chr1", Length: 1000, Mapped: 3, Unmapped: 0},
		{Name: "chr2", Length: 500, Mapped: 1, Unmapped: 0},
		{Name: "*", Length: 0, Mapped: 0, Unmapped: 1},
	}
	for i, w := range want {
		if rows[i] != w {
			t.Errorf("row %d: got %+v want %+v", i, rows[i], w)
		}
	}
}

// ---- parity -----------------------------------------------------------

// Parity p01 — basic.sam: shape and totals match the upstream contract:
// "name<TAB>len<TAB>mapped<TAB>unmapped" with a final "*" row whose
// fourth column is the n_no_coor counter from the .bai.
func TestParity_Idxstats_P01_BasicShape(t *testing.T) {
	bam := makeIndexedBAM(t, parityPath(t, "basic.sam"))
	var buf bytes.Buffer
	if err := IdxstatsFile(bam, &buf); err != nil {
		t.Fatalf("IdxstatsFile: %v", err)
	}
	got := buf.String()
	for _, line := range strings.Split(strings.TrimRight(got, "\n"), "\n") {
		fields := strings.Split(line, "\t")
		if len(fields) != 4 {
			t.Errorf("line %q: got %d fields, want 4", line, len(fields))
		}
	}
	// The trailing line MUST be the "*" unplaced row.
	last := strings.TrimRight(got, "\n")
	idx := strings.LastIndex(last, "\n")
	if idx < 0 || !strings.HasPrefix(last[idx+1:], "*\t") {
		t.Errorf("last line is not the unplaced row: %q", last)
	}
}

// Parity p02 — round-trip: read header @SQ table, idxstats output's
// per-ref order matches @SQ order exactly.
func TestParity_Idxstats_P02_RowOrderMatchesSQ(t *testing.T) {
	bam := makeIndexedBAM(t, parityPath(t, "basic.sam"))
	rows, err := Idxstats(bam)
	if err != nil {
		t.Fatalf("Idxstats: %v", err)
	}
	f, err := os.Open(bam)
	if err != nil {
		t.Fatalf("open: %v", err)
	}
	defer f.Close()
	br, err := sam.NewBAMReader(f)
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	defer br.Close()
	hdr := br.Header()
	for i, ref := range hdr.Refs {
		if rows[i].Name != ref.Name {
			t.Errorf("row %d: got %q want %q", i, rows[i].Name, ref.Name)
		}
	}
}

// Parity p03 — empty BAM: trailer row present, all counts zero.
func TestParity_Idxstats_P03_EmptyOK(t *testing.T) {
	bam := makeIndexedBAM(t, parityPath(t, "empty.sam"))
	var buf bytes.Buffer
	if err := IdxstatsFile(bam, &buf); err != nil {
		t.Fatalf("IdxstatsFile: %v", err)
	}
	out := buf.String()
	if !strings.HasSuffix(strings.TrimRight(out, "\n"), "*\t0\t0\t0") {
		t.Errorf("empty trailer wrong; got %q", out)
	}
}
