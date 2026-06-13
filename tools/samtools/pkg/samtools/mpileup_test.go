package samtools

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

// runMpileupOnSAM feeds one or more SAM-text inputs through Mpileup() with
// the supplied options and returns the emitted text. refFA and posFilter
// may be nil.
func runMpileupOnSAM(t *testing.T, sams []string, opts MpileupOptions, refFA *fasta.RandomAccess, pf *positionFilter) string {
	t.Helper()
	readers := make([]io.Reader, len(sams))
	for i, s := range sams {
		readers[i] = strings.NewReader(s)
	}
	var buf bytes.Buffer
	if err := Mpileup(readers, &buf, opts, refFA, pf); err != nil {
		t.Fatalf("Mpileup: %v", err)
	}
	return buf.String()
}

// minimalSAM is a SAM header + one read for the simplest possible test:
// 5M starting at chr1:10 against an N reference (no FASTA).
const minimalSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:200
r1	0	chr1	10	60	5M	*	0	0	ACGTA	IIIII
`

func TestMpileup_SingleRead_MatchesAtEachPosition(t *testing.T) {
	out := runMpileupOnSAM(t, []string{minimalSAM}, MpileupOptions{}, nil, nil)
	// Should emit lines for positions 10..14. ref column is N (no FASTA),
	// so every base is a mismatch ⇒ printed uppercase as-is. The first
	// position has the ^<charq> marker; the last has $.
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 5 {
		t.Fatalf("want 5 lines, got %d:\n%s", len(lines), out)
	}
	// Position 10: ^]A  (60 + 33 = 93 = ']')
	if !strings.HasPrefix(lines[0], "chr1\t10\tN\t1\t^]A\tI") {
		t.Errorf("line[0] = %q", lines[0])
	}
	// Middle position 12: just G.
	if lines[2] != "chr1\t12\tN\t1\tG\tI" {
		t.Errorf("line[2] = %q", lines[2])
	}
	// Last position 14: A$
	if !strings.HasSuffix(lines[4], "A$\tI") {
		t.Errorf("line[4] = %q", lines[4])
	}
}

func TestMpileup_ReferenceMatch_DotComma(t *testing.T) {
	// Two reads, one on each strand, with sequence that matches the
	// reference. Expected bases column: ".," .
	dir := t.TempDir()
	faPath := writeFasta(t, dir, "chr1", "AAAAAAAAAAACGTAAAAAAAAAAAAAAAA") // 30 bases
	ra, err := fasta.OpenRandomAccess(faPath)
	if err != nil {
		t.Fatalf("open ref: %v", err)
	}
	defer ra.Close()

	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
@RG	ID:x
r1	0	chr1	12	60	3M	*	0	0	CGT	III
r2	16	chr1	12	60	3M	*	0	0	CGT	III
`
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{MinBaseQ: 0}, ra, nil)
	// Position 12 is C (reference), both reads match. Forward = ".", reverse = ",".
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 3 {
		t.Fatalf("want 3 lines, got %d:\n%s", len(lines), out)
	}
	// Format: chr1\t12\tC\t2\t^]. ^],\tII (with start markers)
	if !strings.Contains(lines[0], "\tC\t2\t^].^],\t") {
		t.Errorf("line[0] missing dot/comma forward+reverse match: %q", lines[0])
	}
	// Middle: just .,
	if !strings.Contains(lines[1], "\tG\t2\t.,\t") {
		t.Errorf("line[1] = %q", lines[1])
	}
	// Last: end markers
	if !strings.Contains(lines[2], "\tT\t2\t.$,$\t") {
		t.Errorf("line[2] = %q", lines[2])
	}
}

func TestMpileup_Insertion(t *testing.T) {
	// Read with 2M2I2M at chr1:5: matches 2 ref bases, inserts 2, matches 2.
	// Expected at position 6: the base + "+2NN" annotation (no ref → 'NN'
	// uppercase since we have no FASTA → falls back to read bases).
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
r1	0	chr1	5	60	2M2I2M	*	0	0	ACGTAC	IIIIII
`
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{}, nil, nil)
	// Position 5 has the first M base 'A'. Position 6 has the second M
	// base 'C' followed by the +2GT insertion annotation.
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 4 {
		t.Fatalf("want 4 lines, got %d:\n%s", len(lines), out)
	}
	// Line for position 6 should embed "+2GT"
	if !strings.Contains(lines[1], "+2GT") {
		t.Errorf("line[1] missing +2GT: %q", lines[1])
	}
}

func TestMpileup_Deletion(t *testing.T) {
	// Read with 2M3D2M: matches 2 ref bases, deletes 3, matches 2.
	// Position 5: M-base, with "-3NNN" annotation. Position 6: M-base.
	// Positions 7-9: '*' deletion placeholders. Positions 10-11: M-bases.
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
r1	0	chr1	5	60	2M3D2M	*	0	0	ACGT	IIII
`
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{}, nil, nil)
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 7 {
		t.Fatalf("want 7 lines, got %d:\n%s", len(lines), out)
	}
	// Position 6 (second M-base) carries the -3NNN annotation.
	if !strings.Contains(lines[1], "-3NNN") {
		t.Errorf("line[1] missing -3NNN: %q", lines[1])
	}
	// Positions 7, 8, 9: '*' placeholders. Upstream renders these with
	// the quality of the *next* M-base in the read (here 'I' = 40 + 33),
	// not Phred 0 — see mp_D.out for the canonical example.
	for i, idx := range []int{2, 3, 4} {
		if !strings.Contains(lines[idx], "\t1\t*\tI") {
			t.Errorf("position %d: line[%d] should have deletion placeholder with next-base qual: %q", i+7, idx, lines[idx])
		}
	}
}

func TestMpileup_StartEndMarkers(t *testing.T) {
	// Two reads sharing the same position so we see both ^ and $ on
	// the same column when one ends and one starts.
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
r1	0	chr1	5	60	3M	*	0	0	ACG	III
r2	16	chr1	7	30	2M	*	0	0	GT	II
`
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{}, nil, nil)
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	// Position 7: r1's last (G$), r2's first (^?g — reverse-strand g, MAPQ 30+33='?').
	// Expected bases: G$^?g
	wantSub := "G$^?g"
	found := false
	for _, l := range lines {
		if strings.Contains(l, wantSub) {
			found = true
			break
		}
	}
	if !found {
		t.Errorf("missing %q in:\n%s", wantSub, out)
	}
}

func TestMpileup_MultiBAM(t *testing.T) {
	a := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
ra	0	chr1	10	60	3M	*	0	0	ACG	III
`
	b := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
rb	0	chr1	10	60	3M	*	0	0	TTT	III
`
	out := runMpileupOnSAM(t, []string{a, b}, MpileupOptions{}, nil, nil)
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 3 {
		t.Fatalf("want 3 lines, got %d:\n%s", len(lines), out)
	}
	// Each line should have 7 tab-separated columns: chrom, pos, ref,
	// d1, b1, q1, d2, b2, q2.
	for i, l := range lines {
		cols := strings.Split(l, "\t")
		if len(cols) != 9 {
			t.Errorf("line[%d] has %d cols, want 9: %q", i, len(cols), l)
		}
	}
}

func TestMpileup_RegionFilter(t *testing.T) {
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:200
r1	0	chr1	10	60	5M	*	0	0	ACGTA	IIIII
r2	0	chr1	100	60	5M	*	0	0	ACGTA	IIIII
`
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{Regions: []string{"chr1:100-104"}}, nil, nil)
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 5 {
		t.Fatalf("want 5 lines (positions 100..104), got %d:\n%s", len(lines), out)
	}
	for _, l := range lines {
		if !strings.HasPrefix(l, "chr1\t10") && !strings.HasPrefix(l, "chr1\t1") {
			continue
		}
		if strings.HasPrefix(l, "chr1\t10\t") || strings.HasPrefix(l, "chr1\t11\t") || strings.HasPrefix(l, "chr1\t12\t") {
			t.Errorf("region filter let through a position outside [100,104]: %q", l)
		}
	}
}

func TestMpileup_MinMAPQ_Drops(t *testing.T) {
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
r1	0	chr1	10	5	3M	*	0	0	ACG	III
`
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{MinMAPQ: 10}, nil, nil)
	if out != "" {
		t.Errorf("expected empty output (read filtered by MAPQ): %q", out)
	}
}

func TestMpileup_MinBaseQ_Drops(t *testing.T) {
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
r1	0	chr1	10	60	3M	*	0	0	ACG	!!!
`
	// Phred qualities "!!!" are 0 each; MinBaseQ 10 should drop them all,
	// yielding zero-depth output (which is omitted by default).
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{MinBaseQ: 10}, nil, nil)
	if out != "" {
		t.Errorf("expected empty output (bases all below MinBaseQ): %q", out)
	}
}

func TestMpileup_IgnoreOverlaps_DropsOneHalf(t *testing.T) {
	// Two paired reads (same QName) overlapping at chr1:10-12.
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
@RG	ID:x
pair	99	chr1	10	60	3M	=	10	5	ACG	III
pair	147	chr1	10	60	3M	=	10	-5	ACG	III
`
	// Without -x: depth = 2 everywhere.
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{}, nil, nil)
	if !strings.Contains(out, "\t2\t") {
		t.Errorf("without -x expected depth=2:\n%s", out)
	}
	// With -x: depth = 1 (one half of the pair is dropped).
	out2 := runMpileupOnSAM(t, []string{sam}, MpileupOptions{IgnoreOverlaps: true}, nil, nil)
	if !strings.Contains(out2, "\t1\t") {
		t.Errorf("with -x expected depth=1:\n%s", out2)
	}
}

func TestMpileup_AllPositions_InCoveredRange(t *testing.T) {
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
r1	0	chr1	10	60	3M	*	0	0	ACG	III
r2	0	chr1	15	60	3M	*	0	0	ACG	III
`
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{AllPositions: true}, nil, nil)
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	// Positions 10..17 inclusive (covered range), with 13,14 zero-depth.
	if len(lines) != 8 {
		t.Fatalf("want 8 lines (positions 10..17), got %d:\n%s", len(lines), out)
	}
	// Find the position 13 line: depth must be 0 and bases column '*'.
	found := false
	for _, l := range lines {
		if strings.HasPrefix(l, "chr1\t13\tN\t0\t*\t*") {
			found = true
		}
	}
	if !found {
		t.Errorf("missing zero-depth row at chr1:13:\n%s", out)
	}
}

func TestMpileup_AllPositionsAllChroms_EmitsEvenEmptyChrom(t *testing.T) {
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:5
@SQ	SN:chr2	LN:3
r1	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
`
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{AllPositionsAllChroms: true}, nil, nil)
	// Must include chr2 positions 1,2,3 with zero depth.
	for _, want := range []string{"chr2\t1\tN\t0\t*\t*", "chr2\t2\tN\t0\t*\t*", "chr2\t3\tN\t0\t*\t*"} {
		if !strings.Contains(out, want) {
			t.Errorf("missing %q in:\n%s", want, out)
		}
	}
}

func TestMpileup_OutputMapQAndBP(t *testing.T) {
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
r1	0	chr1	10	60	3M	*	0	0	ACG	III
`
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{OutputMapQ: true, OutputBP: true}, nil, nil)
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 3 {
		t.Fatalf("want 3 lines, got %d:\n%s", len(lines), out)
	}
	// 9 columns per line: chrom pos ref d b q mapq bp
	for i, l := range lines {
		cols := strings.Split(l, "\t")
		if len(cols) != 8 {
			t.Errorf("line[%d] has %d cols, want 8: %q", i, len(cols), l)
		}
	}
	// MAPQ column = ']' (60 + 33). BP column = read-base index 1..3.
	if !strings.Contains(lines[0], "\t]\t1") {
		t.Errorf("line[0] missing MAPQ ']' and BP '1': %q", lines[0])
	}
	if !strings.Contains(lines[2], "\t]\t3") {
		t.Errorf("line[2] missing BP '3': %q", lines[2])
	}
}

func TestMpileup_FastaRef_ProvidesRefBaseColumn(t *testing.T) {
	dir := t.TempDir()
	faPath := writeFasta(t, dir, "chr1", "AAAAAACGTAAAAAAA")
	ra, err := fasta.OpenRandomAccess(faPath)
	if err != nil {
		t.Fatalf("open ref: %v", err)
	}
	defer ra.Close()
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:16
r1	0	chr1	7	60	3M	*	0	0	CGT	III
`
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{}, ra, nil)
	if !strings.Contains(out, "chr1\t7\tC\t1\t^].\t") {
		t.Errorf("expected ref-base C and dot match at position 7:\n%s", out)
	}
}

// TestMpileup_BAQ_LowersQualities verifies that BAQ realignment runs by
// default when a reference is supplied (lowering the per-base qualities
// near the read ends), that -B/NoBAQ leaves the raw qualities untouched,
// and that -E/RedoBAQ is accepted rather than rejected.
func TestMpileup_BAQ_LowersQualities(t *testing.T) {
	dir := t.TempDir()
	// A 20 bp read that matches the reference exactly. BAQ tapers the
	// alignment-uncertainty cost in from each end, so the qualities near
	// the read edges drop below the flat input Phred 40 ('I').
	refSeq := "TTACGTACGTGGCCAATTGGACGTACGTAA"
	faPath := writeFasta(t, dir, "chr1", refSeq)
	ra, err := fasta.OpenRandomAccess(faPath)
	if err != nil {
		t.Fatalf("open ref: %v", err)
	}
	defer ra.Close()
	readSeq := refSeq[5:25] // 20 bp, exact match starting at 1-based pos 6
	sam := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:30\n" +
		"r1\t0\tchr1\t6\t60\t20M\t*\t0\t0\t" + readSeq + "\t" + strings.Repeat("I", 20) + "\n"

	// Default (BAQ on): at least one quality must drop below 'I'.
	withBAQ := runMpileupOnSAM(t, []string{sam}, MpileupOptions{}, ra, nil)
	if !mpileupHasLoweredQual(withBAQ) {
		t.Errorf("expected BAQ to lower some qualities below 'I':\n%s", withBAQ)
	}

	// -B (NoBAQ): every emitted base keeps its raw 'I' quality.
	noBAQ := runMpileupOnSAM(t, []string{sam}, MpileupOptions{NoBAQ: true}, ra, nil)
	if mpileupHasLoweredQual(noBAQ) {
		t.Errorf("with -B no quality should be lowered:\n%s", noBAQ)
	}

	// -E (RedoBAQ): accepted, and (like the default) applies BAQ.
	redoBAQ := runMpileupOnSAM(t, []string{sam}, MpileupOptions{RedoBAQ: true}, ra, nil)
	if !mpileupHasLoweredQual(redoBAQ) {
		t.Errorf("expected -E BAQ to lower some qualities:\n%s", redoBAQ)
	}
}

// mpileupHasLoweredQual reports whether any quality character in the 6th
// (quals) column of a text-mpileup output is below 'I' (Phred 40), i.e.
// BAQ lowered at least one base.
func mpileupHasLoweredQual(out string) bool {
	for _, ln := range strings.Split(strings.TrimRight(out, "\n"), "\n") {
		cols := strings.Split(ln, "\t")
		if len(cols) < 6 {
			continue
		}
		for i := 0; i < len(cols[5]); i++ {
			if cols[5][i] < 'I' {
				return true
			}
		}
	}
	return false
}

// writeFasta writes a tiny single-contig FASTA, builds its .fai sibling
// for random access, and returns the FASTA path.
func writeFasta(t *testing.T, dir, name, seq string) string {
	t.Helper()
	path := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(path, []byte(">"+name+"\n"+seq+"\n"), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	idx, err := fasta.BuildIndex(path)
	if err != nil {
		t.Fatalf("build index: %v", err)
	}
	if err := idx.Save(path + ".fai"); err != nil {
		t.Fatalf("save fai: %v", err)
	}
	return path
}

func TestMpileup_PositionsFile_Filters(t *testing.T) {
	dir := t.TempDir()
	posPath := filepath.Join(dir, "positions.tsv")
	if err := os.WriteFile(posPath, []byte("chr1\t11\nchr1\t13\n"), 0o644); err != nil {
		t.Fatalf("write positions: %v", err)
	}
	pf, err := loadPositionsFile(posPath)
	if err != nil {
		t.Fatalf("loadPositionsFile: %v", err)
	}
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
r1	0	chr1	10	60	5M	*	0	0	ACGTA	IIIII
`
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{}, nil, pf)
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 2 {
		t.Fatalf("positions filter: want 2 lines, got %d:\n%s", len(lines), out)
	}
	if !strings.HasPrefix(lines[0], "chr1\t11\t") || !strings.HasPrefix(lines[1], "chr1\t13\t") {
		t.Errorf("wrong positions emitted: %q %q", lines[0], lines[1])
	}
}

func TestMpileup_PositionsFile_BedFormat(t *testing.T) {
	dir := t.TempDir()
	bedPath := filepath.Join(dir, "regions.bed")
	if err := os.WriteFile(bedPath, []byte("chr1\t10\t12\n"), 0o644); err != nil {
		t.Fatalf("write bed: %v", err)
	}
	pf, err := loadPositionsFile(bedPath)
	if err != nil {
		t.Fatalf("loadPositionsFile: %v", err)
	}
	// BED 10..12 (half-open) => 1-based positions 11 and 12.
	if !pf.contains("chr1", 11) || !pf.contains("chr1", 12) {
		t.Errorf("expected positions 11 and 12 in filter")
	}
	if pf.contains("chr1", 10) || pf.contains("chr1", 13) {
		t.Errorf("filter should not include 10 or 13: %v", pf.byChrom)
	}
}

func TestMpileup_MaxDepth_Caps(t *testing.T) {
	// Build a SAM with 5 reads all starting at position 10.
	var b strings.Builder
	b.WriteString("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:30\n")
	for i := 0; i < 5; i++ {
		b.WriteString("r")
		b.WriteString(strings.Repeat(string(byte('0'+byte(i))), 1))
		b.WriteString("\t0\tchr1\t10\t60\t1M\t*\t0\t0\tA\tI\n")
	}
	out := runMpileupOnSAM(t, []string{b.String()}, MpileupOptions{MaxDepth: 3}, nil, nil)
	// Expect depth column to read 3, not 5.
	if !strings.Contains(out, "\tchr1\t10\tN\t3\t") && !strings.HasPrefix(out, "chr1\t10\tN\t3\t") {
		t.Errorf("expected max-depth cap to set depth=3: %q", out)
	}
}

func TestMpileupFile_EndToEnd(t *testing.T) {
	// File-driven entry point through MpileupFile, including the BAM-list
	// file path and the FASTA reference path.
	dir := t.TempDir()
	samPath := filepath.Join(dir, "in.sam")
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
r1	0	chr1	5	60	5M	*	0	0	ACGTA	IIIII
`
	if err := os.WriteFile(samPath, []byte(sam), 0o644); err != nil {
		t.Fatalf("write sam: %v", err)
	}
	faPath := writeFasta(t, dir, "chr1", strings.Repeat("A", 30))
	outPath := filepath.Join(dir, "pileup.tsv")
	outF, err := os.Create(outPath)
	if err != nil {
		t.Fatalf("create out: %v", err)
	}
	err = MpileupFile(MpileupOptions{
		Inputs:   []string{samPath},
		FastaRef: faPath,
		MinBaseQ: 0,
	}, outF)
	outF.Close()
	if err != nil {
		t.Fatalf("MpileupFile: %v", err)
	}
	body, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatalf("read out: %v", err)
	}
	if !strings.Contains(string(body), "chr1\t5\tA\t1\t^].\t") {
		t.Errorf("expected dot match at position 5 with ref A:\n%s", body)
	}
}

func TestMpileupFile_BamList(t *testing.T) {
	dir := t.TempDir()
	a := filepath.Join(dir, "a.sam")
	b := filepath.Join(dir, "b.sam")
	if err := os.WriteFile(a, []byte(`@HD	VN:1.6
@SQ	SN:chr1	LN:30
r1	0	chr1	10	60	3M	*	0	0	ACG	III
`), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(b, []byte(`@HD	VN:1.6
@SQ	SN:chr1	LN:30
r2	0	chr1	10	60	3M	*	0	0	TTT	III
`), 0o644); err != nil {
		t.Fatal(err)
	}
	list := filepath.Join(dir, "bams.list")
	if err := os.WriteFile(list, []byte("# a comment\n"+a+"\n"+b+"\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	var buf bytes.Buffer
	if err := MpileupFile(MpileupOptions{BamList: list}, &buf); err != nil {
		t.Fatalf("MpileupFile: %v", err)
	}
	if !strings.Contains(buf.String(), "\tchr1\t10\tN\t1\tA\tI\t1\tT\tI") &&
		!strings.HasPrefix(buf.String(), "chr1\t10\tN\t1\t^]A\tI\t1\t^]T\tI") {
		t.Errorf("expected two-BAM pileup at position 10:\n%s", buf.String())
	}
}

func TestMpileup_AnomalousPair_FilteredByDefault(t *testing.T) {
	// One read paired-but-mate-unmapped. Default behaviour: drop it.
	// With -A (CountOrphans): keep it.
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:30
r1	9	chr1	10	60	3M	*	0	0	ACG	III
`
	// Flag 9 = 0x1 (paired) | 0x8 (mate-unmapped) — anomalous.
	out := runMpileupOnSAM(t, []string{sam}, MpileupOptions{}, nil, nil)
	if out != "" {
		t.Errorf("expected anomalous pair dropped by default: %q", out)
	}
	out2 := runMpileupOnSAM(t, []string{sam}, MpileupOptions{CountOrphans: true}, nil, nil)
	if out2 == "" {
		t.Errorf("expected anomalous pair kept with -A")
	}
}

// multiContigGapSAM is a hand-built fixture exercising -aa zero-fill across:
//   - chr1 (LN=8): one read covering positions 3..5; positions 1,2,6,7,8 empty.
//   - chr2 (LN=4): completely uncovered.
//   - chr3 (LN=6): one read covering positions 2..3; positions 1,4,5,6 empty.
const multiContigGapSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:8
@SQ	SN:chr2	LN:4
@SQ	SN:chr3	LN:6
r1	0	chr1	3	60	3M	*	0	0	ACG	III
r2	0	chr3	2	60	2M	*	0	0	AC	II
`

// TestMpileup_AA_ZeroFillTableDriven asserts the exact set of rows emitted
// for -aa across multi-contig input with both covered and entirely-empty
// contigs. The default (no -a, no -aa) emits only covered positions; -a
// extends to the covered range per chrom; -aa extends to every position of
// every contig.
func TestMpileup_AA_ZeroFillTableDriven(t *testing.T) {
	// Build expected output for the -aa case explicitly so we catch any
	// drift in row formatting.
	want := strings.Join([]string{
		"chr1\t1\tN\t0\t*\t*",
		"chr1\t2\tN\t0\t*\t*",
		"chr1\t3\tN\t1\t^]A\tI",
		"chr1\t4\tN\t1\tC\tI",
		"chr1\t5\tN\t1\tG$\tI",
		"chr1\t6\tN\t0\t*\t*",
		"chr1\t7\tN\t0\t*\t*",
		"chr1\t8\tN\t0\t*\t*",
		"chr2\t1\tN\t0\t*\t*",
		"chr2\t2\tN\t0\t*\t*",
		"chr2\t3\tN\t0\t*\t*",
		"chr2\t4\tN\t0\t*\t*",
		"chr3\t1\tN\t0\t*\t*",
		"chr3\t2\tN\t1\t^]A\tI",
		"chr3\t3\tN\t1\tC$\tI",
		"chr3\t4\tN\t0\t*\t*",
		"chr3\t5\tN\t0\t*\t*",
		"chr3\t6\tN\t0\t*\t*",
	}, "\n") + "\n"

	tests := []struct {
		name string
		opts MpileupOptions
		want string
	}{
		{
			name: "default emits only covered positions",
			opts: MpileupOptions{},
			want: strings.Join([]string{
				"chr1\t3\tN\t1\t^]A\tI",
				"chr1\t4\tN\t1\tC\tI",
				"chr1\t5\tN\t1\tG$\tI",
				"chr3\t2\tN\t1\t^]A\tI",
				"chr3\t3\tN\t1\tC$\tI",
			}, "\n") + "\n",
		},
		{
			name: "a emits zero-depth inside covered range, skips empty chr2",
			opts: MpileupOptions{AllPositions: true},
			// chr1 covered range = [3,6), chr3 covered = [2,4).
			// -a does NOT extend to chr2 or to positions outside the
			// per-chrom covered range.
			want: strings.Join([]string{
				"chr1\t3\tN\t1\t^]A\tI",
				"chr1\t4\tN\t1\tC\tI",
				"chr1\t5\tN\t1\tG$\tI",
				"chr3\t2\tN\t1\t^]A\tI",
				"chr3\t3\tN\t1\tC$\tI",
			}, "\n") + "\n",
		},
		{
			name: "aa emits every position of every contig including empty chr2",
			opts: MpileupOptions{AllPositionsAllChroms: true},
			want: want,
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := runMpileupOnSAM(t, []string{multiContigGapSAM}, tc.opts, nil, nil)
			if got != tc.want {
				t.Errorf("output mismatch.\ngot:\n%s\nwant:\n%s", got, tc.want)
			}
		})
	}
}
