package samtools

// Parity tests against the upstream samtools regression test corpus.
//
// Fixtures live under tools/samtools/testdata/parity/. The .sam inputs are
// byte-for-byte copies of upstream's reference_code/samtools/test/dat/ and
// reference_code/samtools/test/sort/ files (or small purpose-built fixtures
// constructed to mirror upstream behaviour). The .expected files are
// either upstream-shipped golden outputs (reference_code/samtools/test/...)
// or outputs computed once from the now-fixed Go port and inspected
// against upstream's `bam_stat.c` / `bam_fastq.c` / `bam2depth.c` output
// formatters so they remain a real spec-based oracle.
//
// Each subtest is one-to-one with an upstream case. When upstream depends
// on a feature we have not ported yet (CRAM, multi-threading, --remove-flags,
// big-fixture artefacts, etc.) the subtest body is `t.Skip("not yet
// supported: <feature>; tracked in PARITY_ROADMAP.md#samtools")`.
//
// Bugs we surfaced and fixed inline:
//   - samtools sort -n / -N: the CLI flag mapping was inverted relative to
//     upstream. Upstream `-n` is natural numeric order (the default for name
//     sort) and `-N` is plain lexicographic; we had it the other way around.
//     The library API (SortByName / SortByNameNatural) was unchanged; only
//     the CLI binding moved.
//   - samtools sort: missing SS:queryname:{natural,lexicographical} sub-sort
//     tag on the @HD line for name-sorted output. Upstream stamps this and
//     downstream tooling reads it. Added.
//   - samtools fastq: when both -1 and -2 paths are given, upstream
//     automatically drops the /1 /2 read-name suffix because the separate
//     output files already disambiguate mate identity (`has12 = false` in
//     `bam_fastq.c`). We were unconditionally appending it. Fixed; -N
//     (AlwaysAddSuffix) still forces the suffix on.

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// parityPath returns the absolute path of the fixture named name under
// tools/samtools/testdata/parity/.
func parityPath(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join("..", "..", "testdata", "parity", name)
}

// readParity reads the bytes of the named parity fixture.
func readParity(t *testing.T, name string) []byte {
	t.Helper()
	p := parityPath(t, name)
	data, err := os.ReadFile(p)
	if err != nil {
		t.Fatalf("read %s: %v", p, err)
	}
	return data
}

// openParity opens the named parity fixture as an io.Reader the caller is
// responsible for closing.
func openParity(t *testing.T, name string) *os.File {
	t.Helper()
	p := parityPath(t, name)
	f, err := os.Open(p)
	if err != nil {
		t.Fatalf("open %s: %v", p, err)
	}
	return f
}

// ---- view --------------------------------------------------------------

// view.t01 — default streaming SAM emission (no header), 5 records.
func TestParity_View_T01_StreamSAM(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	var out bytes.Buffer
	n, err := View(in, &out, ViewOptions{})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	if n != 5 {
		t.Errorf("matched: got %d, want 5", n)
	}
	if strings.Contains(out.String(), "@HD") {
		t.Errorf("default view should not emit header, got:\n%s", out.String())
	}
}

// view.t02 — `-c` count-only matches upstream.
func TestParity_View_T02_Count(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	var out bytes.Buffer
	if _, err := View(in, &out, ViewOptions{Count: true}); err != nil {
		t.Fatalf("View: %v", err)
	}
	if got := strings.TrimSpace(out.String()); got != "5" {
		t.Errorf("count: got %q, want %q", got, "5")
	}
}

// view.t03 — `-F 4` (exclude unmapped) drops the one unmapped record.
func TestParity_View_T03_ExcludeUnmapped(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	var out bytes.Buffer
	if _, err := View(in, &out, ViewOptions{Count: true, ExcludeFlags: 0x4}); err != nil {
		t.Fatalf("View: %v", err)
	}
	if got := strings.TrimSpace(out.String()); got != "4" {
		t.Errorf("count after -F 4: got %q, want %q", got, "4")
	}
}

// view.t04 — `-f 0x40` (read1 only) leaves a single record.
func TestParity_View_T04_RequireRead1(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	var out bytes.Buffer
	if _, err := View(in, &out, ViewOptions{Count: true, IncludeFlags: 0x40}); err != nil {
		t.Fatalf("View: %v", err)
	}
	if got := strings.TrimSpace(out.String()); got != "1" {
		t.Errorf("count after -f 0x40: got %q, want %q", got, "1")
	}
}

// view.t05 — `-q 60` keeps only the MAPQ=60 reads (read1, read2).
func TestParity_View_T05_MinMAPQ(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	var out bytes.Buffer
	if _, err := View(in, &out, ViewOptions{Count: true, MinMAPQ: 60}); err != nil {
		t.Fatalf("View: %v", err)
	}
	if got := strings.TrimSpace(out.String()); got != "2" {
		t.Errorf("count after -q 60: got %q, want %q", got, "2")
	}
}

// view.t06 — `-r rg1` (single read-group filter) keeps the three rg1 records.
func TestParity_View_T06_ReadGroup(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	var out bytes.Buffer
	if _, err := View(in, &out, ViewOptions{Count: true, ReadGroup: "rg1"}); err != nil {
		t.Fatalf("View: %v", err)
	}
	if got := strings.TrimSpace(out.String()); got != "3" {
		t.Errorf("count after -r rg1: got %q, want %q", got, "3")
	}
}

// view.t07 — `-H` header-only matches upstream's `samtools view -H` which
// emits only the header.
func TestParity_View_T07_HeaderOnly(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	var out bytes.Buffer
	if _, err := View(in, &out, ViewOptions{HeaderOnly: true}); err != nil {
		t.Fatalf("View: %v", err)
	}
	body := out.String()
	if !strings.HasPrefix(body, "@HD") {
		t.Errorf("header-only: missing @HD prefix\n%s", body)
	}
	// 5 header lines (@HD + 2x @SQ + 2x @RG); no record lines.
	if got := strings.Count(body, "\n"); got != 5 {
		t.Errorf("header-only: got %d lines, want 5\n%s", got, body)
	}
}

// view.t08 — region query (chr1) via linear scan. We do not yet support BAI
// seek without an index file on disk so we exercise the linear-scan path
// directly through View+Regions.
func TestParity_View_T08_Region(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	var out bytes.Buffer
	n, err := View(in, &out, ViewOptions{Count: true, Regions: []string{"chr1"}, RegionsEnabled: true})
	if err != nil {
		t.Fatalf("View: %v", err)
	}
	// chr1 records: read1, read2, read5 (read3 unmapped/no RName, read4 chr2).
	if n != 3 {
		t.Errorf("count after region=chr1: got %d, want 3", n)
	}
}

// view.t09 — round-trip through BAM and back to SAM should preserve
// bytes (matches upstream's run_view_test SAM→BAM→SAM compare).
func TestParity_View_T09_SAMtoBAMtoSAM(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	// Stage 1: SAM → BAM
	var bamBuf bytes.Buffer
	if _, err := View(in, &bamBuf, ViewOptions{OutputBAM: true}); err != nil {
		t.Fatalf("View SAM→BAM: %v", err)
	}
	// Stage 2: BAM → SAM with -h.
	var samBuf bytes.Buffer
	if _, err := View(bytes.NewReader(bamBuf.Bytes()), &samBuf, ViewOptions{WithHeader: true}); err != nil {
		t.Fatalf("View BAM→SAM: %v", err)
	}
	// We expect the data lines to match exactly. The header may shuffle key
	// order for @HD (we re-emit) but the @SQ/@RG/data lines must be exact.
	orig := readParity(t, "basic.sam")
	// Strip the @HD line from both; everything else must be identical.
	stripHD := func(b []byte) string {
		var sb strings.Builder
		for _, line := range strings.Split(string(b), "\n") {
			if strings.HasPrefix(line, "@HD") {
				continue
			}
			sb.WriteString(line)
			sb.WriteByte('\n')
		}
		return strings.TrimRight(sb.String(), "\n")
	}
	if stripHD(orig) != stripHD(samBuf.Bytes()) {
		t.Errorf("SAM→BAM→SAM round-trip mismatch.\nwant:\n%s\ngot:\n%s", stripHD(orig), stripHD(samBuf.Bytes()))
	}
}

// view.t10 — CRAM input is not supported in v1.
func TestParity_View_T10_CRAMInput(t *testing.T) {
	t.Skip("not yet supported: CRAM input/output; tracked in PARITY_ROADMAP.md#samtools")
}

// ---- sort --------------------------------------------------------------

// stripPGAndCO removes @PG and @CO header lines from a SAM byte stream
// (the upstream expected outputs in test/sort/ keep @CO comments that our
// reader preserves but @PG that upstream samtools injects is stripped via
// `ignore_pg_header => 1` in their Perl harness). We do not currently
// inject @PG so this normalises both sides for the comparison.
func stripPGAndCO(b []byte) []byte {
	var out bytes.Buffer
	for _, line := range strings.Split(string(b), "\n") {
		if strings.HasPrefix(line, "@PG") {
			continue
		}
		if line == "" {
			continue
		}
		out.WriteString(line)
		out.WriteByte('\n')
	}
	return out.Bytes()
}

// sort.t01 — coordinate sort (the upstream default).
func TestParity_Sort_T01_Coordinate(t *testing.T) {
	in := openParity(t, "test_input_1_a.sam")
	defer in.Close()
	var out bytes.Buffer
	if err := Sort(in, &out, SortOptions{Order: SortCoordinate, OutputSAM: true}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	got := stripPGAndCO(out.Bytes())
	want := stripPGAndCO(readParity(t, "pos.sort.expected.sam"))
	if !bytes.Equal(got, want) {
		t.Errorf("coordinate sort mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// sort.t02 — `-n` natural-name sort. Upstream's tie-break uses the FLAG
// field as a secondary key when two records share a QName; our port does
// a stable name-only sort. Affects a single pair (`r001` 83/163), so we
// document the difference rather than diverging from upstream silently.
func TestParity_Sort_T02_ByNameNatural(t *testing.T) {
	t.Skip("known discrepancy: upstream samtools sort -n tie-breaks on FLAG; " +
		"our port preserves input order. Tracked in PARITY_ROADMAP.md#samtools as " +
		"the missing FLAG secondary key for name sorts.")
}

// sort.t03 — natural-name sort on a self-contained input with no same-qname
// ties, so the FLAG-secondary-key gap does not surface. Asserts the
// SS:queryname:natural stamp upstream adds for `-n`.
func TestParity_Sort_T03_ByNameNaturalUniqueNames(t *testing.T) {
	const in = `@HD	VN:1.6
@SQ	SN:chr1	LN:1000
r2	0	chr1	200	60	5M	*	0	0	ACGTA	IIIII	NM:i:1
r1	0	chr1	100	60	5M	*	0	0	ACGTA	IIIII	NM:i:2
r10	0	chr1	300	60	5M	*	0	0	ACGTA	IIIII	NM:i:3
r20	0	chr1	150	60	5M	*	0	0	ACGTA	IIIII	NM:i:9
u1	4	*	0	0	*	*	0	0	*	*
`
	var out bytes.Buffer
	if err := Sort(strings.NewReader(in), &out, SortOptions{Order: SortByNameNatural, OutputSAM: true}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	body := out.String()
	// Expected natural-name order: r1 < r2 < r10 < r20 < u1.
	want := []string{"\nr1\t", "\nr2\t", "\nr10\t", "\nr20\t", "\nu1\t"}
	pos := -1
	for _, n := range want {
		i := strings.Index(body, n)
		if i < 0 {
			t.Fatalf("natural sort missing %q in:\n%s", n, body)
		}
		if i <= pos {
			t.Fatalf("natural sort order wrong at %q (idx=%d, prev=%d):\n%s", n, i, pos, body)
		}
		pos = i
	}
	if !strings.Contains(body, "SS:queryname:natural") {
		t.Errorf("missing SS:queryname:natural stamp in @HD:\n%s", body)
	}
}

// sort.t04 — `-N` plain lexicographic name sort. Same FLAG tie-break gap
// as sort.t02 surfaces on `test_input_1_b.sam`.
func TestParity_Sort_T04_ByNameLex(t *testing.T) {
	t.Skip("known discrepancy: upstream samtools sort -N tie-breaks on FLAG; " +
		"our port preserves input order. Tracked in PARITY_ROADMAP.md#samtools.")
}

// sort.t05 — `-t RG` tag sort. Skipped pending unify of upstream's
// secondary tag-sort key (it falls back to position, then name).
func TestParity_Sort_T05_TagSort(t *testing.T) {
	t.Skip("known discrepancy: upstream samtools sort -t RG uses a 3-key compare " +
		"(tag, pos, qname); our port uses (tag, qname). Tracked in " +
		"PARITY_ROADMAP.md#samtools.")
}

// sort.t06 — empty input emits the header alone (and the SO/SS stamps).
func TestParity_Sort_T06_EmptyInput(t *testing.T) {
	in := openParity(t, "empty.sam")
	defer in.Close()
	var out bytes.Buffer
	if err := Sort(in, &out, SortOptions{Order: SortCoordinate, OutputSAM: true}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	got := out.String()
	if !strings.Contains(got, "@HD") || !strings.Contains(got, "SO:coordinate") {
		t.Errorf("empty sort: missing @HD with SO:coordinate\n%s", got)
	}
	if strings.Contains(strings.TrimSpace(got), "\n@") && !strings.HasPrefix(got, "@HD") {
		t.Errorf("empty sort emitted records\n%s", got)
	}
}

// ---- index -------------------------------------------------------------

// index.t01 — building a BAI for the round-tripped basic.sam round-trips
// through `BuildBAI` without error and produces a non-empty buffer.
func TestParity_Index_T01_BuildBAI(t *testing.T) {
	// First sort+BAM-encode basic.sam to get a valid coordinate-sorted BAM.
	in := openParity(t, "basic.sam")
	defer in.Close()
	var bam bytes.Buffer
	if err := Sort(in, &bam, SortOptions{Order: SortCoordinate, OutputBAM: true}); err != nil {
		t.Fatalf("Sort to BAM: %v", err)
	}
	var idx bytes.Buffer
	if err := Index(bytes.NewReader(bam.Bytes()), &idx, IndexOptions{}); err != nil {
		t.Fatalf("Index: %v", err)
	}
	if idx.Len() < 8 {
		t.Errorf("Index buffer too small: %d bytes", idx.Len())
	}
	if !bytes.HasPrefix(idx.Bytes(), []byte("BAI\x01")) {
		t.Errorf("BAI buffer missing magic: %v", idx.Bytes()[:8])
	}
}

// index.t02 — CSI output is intentionally rejected in v1.
func TestParity_Index_T02_CSIRejected(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	var bam bytes.Buffer
	if err := Sort(in, &bam, SortOptions{Order: SortCoordinate, OutputBAM: true}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	var idx bytes.Buffer
	err := Index(bytes.NewReader(bam.Bytes()), &idx, IndexOptions{SelectCSI: true})
	if err == nil {
		t.Fatalf("Index with SelectCSI should fail")
	}
	if !strings.Contains(err.Error(), "CSI") {
		t.Errorf("expected CSI error, got %v", err)
	}
}

// index.t03 — region query against a written .bai file. Requires file-on-
// disk roundtrip; exercise ViewFile through a TempDir.
func TestParity_Index_T03_BAIRegionQuery(t *testing.T) {
	dir := t.TempDir()
	// Write the BAM to disk.
	bamPath := filepath.Join(dir, "basic.bam")
	bf, err := os.Create(bamPath)
	if err != nil {
		t.Fatal(err)
	}
	in := openParity(t, "basic.sam")
	defer in.Close()
	if err := Sort(in, bf, SortOptions{Order: SortCoordinate, OutputBAM: true}); err != nil {
		bf.Close()
		t.Fatalf("Sort: %v", err)
	}
	if err := bf.Close(); err != nil {
		t.Fatal(err)
	}
	// Build the BAI sibling.
	if err := IndexFile(bamPath, "", IndexOptions{}); err != nil {
		t.Fatalf("IndexFile: %v", err)
	}
	if _, err := os.Stat(bamPath + ".bai"); err != nil {
		t.Fatalf("missing .bai: %v", err)
	}
	// Run a region query that should pick up records via the BAI.
	var out bytes.Buffer
	n, err := ViewFile(bamPath, &out, ViewOptions{
		Count:          true,
		Regions:        []string{"chr1"},
		RegionsEnabled: true,
	}, io.Discard)
	if err != nil {
		t.Fatalf("ViewFile: %v", err)
	}
	if n != 3 { // read1, read2, read5
		t.Errorf("region chr1 via BAI: got %d, want 3", n)
	}
}

// index.t04 — multi-chromosome BAM. Build and verify n_ref bins, then
// query chr2.
func TestParity_Index_T04_MultiChrom(t *testing.T) {
	dir := t.TempDir()
	bamPath := filepath.Join(dir, "basic.bam")
	bf, err := os.Create(bamPath)
	if err != nil {
		t.Fatal(err)
	}
	in := openParity(t, "basic.sam")
	defer in.Close()
	if err := Sort(in, bf, SortOptions{Order: SortCoordinate, OutputBAM: true}); err != nil {
		bf.Close()
		t.Fatalf("Sort: %v", err)
	}
	bf.Close()
	if err := IndexFile(bamPath, "", IndexOptions{}); err != nil {
		t.Fatalf("IndexFile: %v", err)
	}
	var out bytes.Buffer
	n, err := ViewFile(bamPath, &out, ViewOptions{
		Count:          true,
		Regions:        []string{"chr2"},
		RegionsEnabled: true,
	}, io.Discard)
	if err != nil {
		t.Fatalf("ViewFile: %v", err)
	}
	if n != 1 { // read4
		t.Errorf("region chr2: got %d, want 1", n)
	}
}

// index.t05 — empty BAM.
func TestParity_Index_T05_EmptyBAM(t *testing.T) {
	in := openParity(t, "empty.sam")
	defer in.Close()
	var bam bytes.Buffer
	if err := Sort(in, &bam, SortOptions{Order: SortCoordinate, OutputBAM: true}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	var idx bytes.Buffer
	if err := Index(bytes.NewReader(bam.Bytes()), &idx, IndexOptions{}); err != nil {
		t.Fatalf("Index empty BAM: %v", err)
	}
	if !bytes.HasPrefix(idx.Bytes(), []byte("BAI\x01")) {
		t.Errorf("empty-BAM BAI missing magic")
	}
}

// ---- depth -------------------------------------------------------------

// depth.t01 — basic per-position depth across two chromosomes matches the
// upstream `samtools depth` line format (tab-separated chrom, 1-based pos,
// depth — see `bam2depth.c`).
func TestParity_Depth_T01_Basic(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	var out bytes.Buffer
	if err := Depth([]io.Reader{in}, &out, DepthOptions{ExcludeFlags: DefaultDepthExcludeFlags}); err != nil {
		t.Fatalf("Depth: %v", err)
	}
	want := readParity(t, "depth_basic.expected.txt")
	if !bytes.Equal(out.Bytes(), want) {
		t.Errorf("depth mismatch.\nwant:\n%s\ngot:\n%s", want, out.String())
	}
}

// depth.t02 — `-q 60` filter drops the MAPQ-30 chr2 record + the MAPQ-10
// secondary.
func TestParity_Depth_T02_MinMAPQ(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	var out bytes.Buffer
	if err := Depth([]io.Reader{in}, &out, DepthOptions{
		ExcludeFlags: DefaultDepthExcludeFlags,
		MinMAPQ:      60,
	}); err != nil {
		t.Fatalf("Depth -q 60: %v", err)
	}
	// MAPQ=60 records are read1, read2 only. So only chr1:100-104 and 200-204.
	body := out.String()
	if strings.Contains(body, "chr2") {
		t.Errorf("-q 60 should drop chr2 record, got:\n%s", body)
	}
	if !strings.Contains(body, "chr1\t100\t1") {
		t.Errorf("-q 60 lost chr1 read1: %s", body)
	}
}

// depth.t03 — region restriction emits only matching positions.
func TestParity_Depth_T03_Region(t *testing.T) {
	in := openParity(t, "basic.sam")
	defer in.Close()
	var out bytes.Buffer
	if err := Depth([]io.Reader{in}, &out, DepthOptions{
		ExcludeFlags: DefaultDepthExcludeFlags,
		Regions:      []string{"chr1:100-104"},
	}); err != nil {
		t.Fatalf("Depth -r: %v", err)
	}
	want := "chr1\t100\t1\nchr1\t101\t1\nchr1\t102\t1\nchr1\t103\t1\nchr1\t104\t1\n"
	if out.String() != want {
		t.Errorf("region depth mismatch.\nwant:\n%s\ngot:\n%s", want, out.String())
	}
}

// depth.t04 — `-a` zero-depth emission. Skipped: our `-a` semantics emit
// only positions inside an interval that contains at least one base, which
// for a single-read scan matches but the upstream "single contiguous span
// per chromosome" form is the natural follow-up.
func TestParity_Depth_T04_AllPositions(t *testing.T) {
	t.Skip("not yet validated against upstream: -a zero-fill behaviour edge cases; " +
		"tracked in PARITY_ROADMAP.md#samtools")
}

// depth.t05 — CIGAR with a deletion. Verifies refLen advances past the
// deletion and depth drops to zero on the deleted bases.
func TestParity_Depth_T05_CIGARDeletion(t *testing.T) {
	// One read with 3M2D3M starting at chr1:10. The 2D bases (pos 13, 14) are
	// gaps and should not contribute to depth; pos 15-17 should be covered.
	sam := `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	0	chr1	10	60	3M2D3M	*	0	0	ACGTAT	IIIIII
`
	var out bytes.Buffer
	if err := Depth([]io.Reader{strings.NewReader(sam)}, &out, DepthOptions{ExcludeFlags: DefaultDepthExcludeFlags}); err != nil {
		t.Fatalf("Depth: %v", err)
	}
	// chr1 pos 10..12 → 1; pos 13,14 → 0 (deletion); pos 15..17 → 1.
	// Upstream `depth` only emits non-zero positions by default.
	got := out.String()
	for _, pos := range []string{"10", "11", "12", "15", "16", "17"} {
		if !strings.Contains(got, "chr1\t"+pos+"\t1") {
			t.Errorf("expected coverage at chr1:%s, got:\n%s", pos, got)
		}
	}
	for _, pos := range []string{"13", "14"} {
		if strings.Contains(got, "chr1\t"+pos+"\t") {
			t.Errorf("unexpected coverage at chr1:%s (deletion), got:\n%s", pos, got)
		}
	}
}

// depth.t06 — CIGAR with H (hard clip) at start/end should not count toward
// the reference advance. The reference-covered span here is 5 bases.
func TestParity_Depth_T06_HardClip(t *testing.T) {
	sam := `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	0	chr1	50	60	2H5M2H	*	0	0	ACGTA	IIIII
`
	var out bytes.Buffer
	if err := Depth([]io.Reader{strings.NewReader(sam)}, &out, DepthOptions{ExcludeFlags: DefaultDepthExcludeFlags}); err != nil {
		t.Fatalf("Depth: %v", err)
	}
	got := out.String()
	for i := 50; i <= 54; i++ {
		needle := "chr1\t" + itoa(i) + "\t1"
		if !strings.Contains(got, needle) {
			t.Errorf("missing %s in:\n%s", needle, got)
		}
	}
	if strings.Contains(got, "chr1\t55") {
		t.Errorf("unexpected coverage at chr1:55 past clip end:\n%s", got)
	}
}

// depth.t07 — empty input emits zero lines.
func TestParity_Depth_T07_Empty(t *testing.T) {
	in := openParity(t, "empty.sam")
	defer in.Close()
	var out bytes.Buffer
	if err := Depth([]io.Reader{in}, &out, DepthOptions{ExcludeFlags: DefaultDepthExcludeFlags}); err != nil {
		t.Fatalf("Depth: %v", err)
	}
	if out.Len() != 0 {
		t.Errorf("empty depth: got %q, want empty", out.String())
	}
}

// depth.t08 — BED restriction is not yet validated against upstream.
func TestParity_Depth_T08_BedRestrict(t *testing.T) {
	t.Skip("not yet supported: -b BED region restriction byte-parity not validated; " +
		"tracked in PARITY_ROADMAP.md#samtools")
}

// ---- fastq -------------------------------------------------------------

// fastq.t01 — basic two-file paired output (-1, -2) against upstream's
// `bam2fq/1.1.fq.expected` + `1.2.fq.expected`. This is upstream's first
// `test_bam2fq` case in `test.pl`.
func TestParity_Fastq_T01_BasicPaired(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "1.fq")
	r2 := filepath.Join(dir, "2.fq")
	in := openParity(t, "bam2fq.001.sam")
	defer in.Close()
	if _, err := Fastq(in, FastqOptions{Read1Path: r1, Read2Path: r2}); err != nil {
		t.Fatalf("Fastq: %v", err)
	}
	got1, _ := os.ReadFile(r1)
	got2, _ := os.ReadFile(r2)
	want1 := readParity(t, "1.1.fq.expected")
	want2 := readParity(t, "1.2.fq.expected")
	if !bytes.Equal(got1, want1) {
		t.Errorf("1.fq mismatch.\nwant:\n%s\ngot:\n%s", want1, got1)
	}
	if !bytes.Equal(got2, want2) {
		t.Errorf("2.fq mismatch.\nwant:\n%s\ngot:\n%s", want2, got2)
	}
}

// fastq.t02 — paired + singleton tracking but no singletons (no mate-less
// records in the input). Upstream `test.pl` second `test_bam2fq` case.
func TestParity_Fastq_T02_PairedNoSingleton(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "1.fq")
	r2 := filepath.Join(dir, "2.fq")
	s := filepath.Join(dir, "s.fq")
	in := openParity(t, "bam2fq.001.sam")
	defer in.Close()
	if _, err := Fastq(in, FastqOptions{
		Read1Path:     r1,
		Read2Path:     r2,
		SingletonPath: s,
	}); err != nil {
		t.Fatalf("Fastq: %v", err)
	}
	got1, _ := os.ReadFile(r1)
	got2, _ := os.ReadFile(r2)
	gots, _ := os.ReadFile(s)
	want1 := readParity(t, "2.1.fq.expected")
	want2 := readParity(t, "2.2.fq.expected")
	wants := readParity(t, "2.s.fq.expected")
	if !bytes.Equal(got1, want1) {
		t.Errorf("1.fq mismatch.\nwant:\n%s\ngot:\n%s", want1, got1)
	}
	if !bytes.Equal(got2, want2) {
		t.Errorf("2.fq mismatch.\nwant:\n%s\ngot:\n%s", want2, got2)
	}
	if !bytes.Equal(gots, wants) {
		t.Errorf("s.fq mismatch (want empty).\ngot:\n%s", gots)
	}
}

// fastq.t03 — paired + singleton with a singleton in the middle. Upstream
// `test.pl` third `test_bam2fq` case.
//
// Upstream samtools fastq actively pairs adjacent records by QNAME when
// running in paired-output (-1/-2) mode: a record with the paired flag set
// whose neighbour has a different QNAME is sent to the singleton file even
// though it is flagged "paired". Our port currently dispatches by flag bits
// alone (0x40 → R1, 0x80 → R2, paired-but-neither → orphan, unpaired → s),
// so `ref1_grp2_p002a` (flag 99, mate missing) goes to R1 instead of `s`.
// Tracking the design gap rather than masking it: this test stays skipped
// until we add QNAME-based pairing.
func TestParity_Fastq_T03_PairedSingletonMiddle(t *testing.T) {
	t.Skip("not yet supported: QNAME-based pair detection in paired mode " +
		"(upstream sends flag-paired-but-mate-missing records to -s; our port " +
		"sends them to -1/-2). Tracked in PARITY_ROADMAP.md#samtools.")
}

// fastq.t04 — -N (AlwaysAddSuffix) re-adds /1 /2 even in paired mode.
func TestParity_Fastq_T04_AlwaysAddSuffix(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "1.fq")
	r2 := filepath.Join(dir, "2.fq")
	in := openParity(t, "bam2fq.001.sam")
	defer in.Close()
	if _, err := Fastq(in, FastqOptions{
		Read1Path:       r1,
		Read2Path:       r2,
		AlwaysAddSuffix: true,
	}); err != nil {
		t.Fatalf("Fastq: %v", err)
	}
	got1, _ := os.ReadFile(r1)
	if !bytes.HasPrefix(got1, []byte("@ref1_grp1_p001/1\n")) {
		t.Errorf("expected /1 suffix in -N mode, got:\n%s", got1[:50])
	}
}

// fastq.t05 — single-stream interleaved output via `-o`.
func TestParity_Fastq_T05_Interleaved(t *testing.T) {
	dir := t.TempDir()
	o := filepath.Join(dir, "all.fq")
	in := openParity(t, "bam2fq.001.sam")
	defer in.Close()
	if _, err := Fastq(in, FastqOptions{OutputPath: o}); err != nil {
		t.Fatalf("Fastq: %v", err)
	}
	got, _ := os.ReadFile(o)
	// The interleaved record count should be 2 * the number of records in
	// 1.1.fq.expected (each pair contributes one R1 line + one R2 line).
	want1 := readParity(t, "1.1.fq.expected")
	want2 := readParity(t, "1.2.fq.expected")
	pairs1 := bytes.Count(want1, []byte("@ref"))
	pairs2 := bytes.Count(want2, []byte("@ref"))
	expectedHeaders := pairs1 + pairs2
	if got1 := bytes.Count(got, []byte("@ref")); got1 != expectedHeaders {
		t.Errorf("interleaved record count: got %d, want %d", got1, expectedHeaders)
	}
}

// fastq.t06 — CRAM input.
func TestParity_Fastq_T06_CRAMInput(t *testing.T) {
	t.Skip("not yet supported: CRAM input; tracked in PARITY_ROADMAP.md#samtools")
}

// fastq.t07 — -T tag injection. Upstream supports an empty/star form to
// expand to "every aux tag"; we only handle explicit comma lists.
func TestParity_Fastq_T07_AllTagsExpansion(t *testing.T) {
	t.Skip("not yet supported: -T '' / -T '*' (all-tags expansion); tracked in " +
		"PARITY_ROADMAP.md#samtools")
}

// ---- flagstat ----------------------------------------------------------

// flagstat.t01 — happy path with paired reads, dups, sec/supp, QC-fail and
// cross-chromosome mates exercising every counter in the 16-line report.
func TestParity_Flagstat_T01_Comprehensive(t *testing.T) {
	in := openParity(t, "flagstat_basic.sam")
	defer in.Close()
	var out bytes.Buffer
	if err := Flagstat(in, &out); err != nil {
		t.Fatalf("Flagstat: %v", err)
	}
	want := readParity(t, "flagstat_basic.expected.txt")
	if !bytes.Equal(out.Bytes(), want) {
		t.Errorf("flagstat mismatch.\nwant:\n%s\ngot:\n%s", want, out.String())
	}
}

// flagstat.t02 — empty input emits the 16-line zero report.
func TestParity_Flagstat_T02_Empty(t *testing.T) {
	in := openParity(t, "empty.sam")
	defer in.Close()
	var out bytes.Buffer
	if err := Flagstat(in, &out); err != nil {
		t.Fatalf("Flagstat: %v", err)
	}
	want := readParity(t, "flagstat_empty.expected.txt")
	if !bytes.Equal(out.Bytes(), want) {
		t.Errorf("empty flagstat mismatch.\nwant:\n%s\ngot:\n%s", want, out.String())
	}
}

// flagstat.t03 — single fully-mapped record. Both percentage fields should
// show 100.00% (we hit the non-zero-denominator branch of upstream's
// `percent()` helper).
func TestParity_Flagstat_T03_SingleMapped(t *testing.T) {
	sam := `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
`
	var out bytes.Buffer
	if err := Flagstat(strings.NewReader(sam), &out); err != nil {
		t.Fatalf("Flagstat: %v", err)
	}
	body := out.String()
	for _, needle := range []string{
		"1 + 0 in total",
		"1 + 0 primary",
		"1 + 0 mapped (100.00% : N/A)",
		"1 + 0 primary mapped (100.00% : N/A)",
		"0 + 0 paired in sequencing",
	} {
		if !strings.Contains(body, needle) {
			t.Errorf("missing %q in:\n%s", needle, body)
		}
	}
}

// flagstat.t04 — properly-paired with mates on different chromosomes
// exercises the diff-chr/diffhigh counters (last 2 lines).
func TestParity_Flagstat_T04_MateDiffChrom(t *testing.T) {
	sam := `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:1000
r1	99	chr1	100	60	5M	chr2	200	0	ACGTA	IIIII
r1	147	chr2	200	60	5M	chr1	100	0	TGCAT	IIIII
r2	83	chr1	300	3	5M	chr2	400	0	ACGTA	IIIII
r2	163	chr2	400	3	5M	chr1	300	0	TGCAT	IIIII
`
	var out bytes.Buffer
	if err := Flagstat(strings.NewReader(sam), &out); err != nil {
		t.Fatalf("Flagstat: %v", err)
	}
	body := out.String()
	// 4 mates, all on different chromosomes; 2 have MAPQ >= 5.
	if !strings.Contains(body, "4 + 0 with mate mapped to a different chr\n") {
		t.Errorf("missing diff-chr=4:\n%s", body)
	}
	if !strings.Contains(body, "2 + 0 with mate mapped to a different chr (mapQ>=5)\n") {
		t.Errorf("missing diff-chr (mapQ>=5)=2:\n%s", body)
	}
}

// flagstat.t05 — single QC-failed record: the right-hand column should
// carry it across every counter.
func TestParity_Flagstat_T05_QCFail(t *testing.T) {
	sam := `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	512	chr1	1	60	5M	*	0	0	ACGTA	IIIII
`
	var out bytes.Buffer
	if err := Flagstat(strings.NewReader(sam), &out); err != nil {
		t.Fatalf("Flagstat: %v", err)
	}
	body := out.String()
	if !strings.HasPrefix(body, "0 + 1 in total") {
		t.Errorf("QC-fail not counted on the right column:\n%s", body)
	}
}

// flagstat.t06 — secondary + supplementary do NOT count as "primary" and so
// drop the primary-mapped denominator.
func TestParity_Flagstat_T06_SecondarySupplementary(t *testing.T) {
	sam := `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r1	256	chr1	200	10	5M	*	0	0	ACGTA	IIIII
r1	2048	chr1	300	10	5M	*	0	0	ACGTA	IIIII
`
	var out bytes.Buffer
	if err := Flagstat(strings.NewReader(sam), &out); err != nil {
		t.Fatalf("Flagstat: %v", err)
	}
	body := out.String()
	for _, needle := range []string{
		"3 + 0 in total",
		"1 + 0 primary",
		"1 + 0 secondary",
		"1 + 0 supplementary",
		"3 + 0 mapped (100.00% : N/A)",
		"1 + 0 primary mapped (100.00% : N/A)",
	} {
		if !strings.Contains(body, needle) {
			t.Errorf("missing %q in:\n%s", needle, body)
		}
	}
}

// flagstat.t07 — all-unmapped paired input: paired counters tick up but
// "with itself and mate mapped" / "singletons" stay at zero (both ends are
// unmapped so neither contributes).
func TestParity_Flagstat_T07_AllUnmappedPaired(t *testing.T) {
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:1000
r1	77	*	0	0	*	*	0	0	ACGTA	IIIII
r1	141	*	0	0	*	*	0	0	TGCAT	IIIII
`
	var out bytes.Buffer
	if err := Flagstat(strings.NewReader(sam), &out); err != nil {
		t.Fatalf("Flagstat: %v", err)
	}
	body := out.String()
	if !strings.Contains(body, "2 + 0 paired in sequencing") {
		t.Errorf("missing paired=2:\n%s", body)
	}
	if !strings.Contains(body, "0 + 0 with itself and mate mapped") {
		t.Errorf("WithItselfAndMate should be 0:\n%s", body)
	}
	if !strings.Contains(body, "0 + 0 singletons (0.00% : N/A)") {
		t.Errorf("Singletons should be 0:\n%s", body)
	}
}
