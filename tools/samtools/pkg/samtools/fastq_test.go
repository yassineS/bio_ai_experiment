package samtools

import (
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// fastqSAM is a hand-crafted SAM with name-sorted records covering every
// case the fastq tests need: paired (R1+R2), an orphan (only 0x40), a
// singleton (not paired), a reverse-strand read whose SEQ in the BAM is
// the reverse-complement, a record with auxiliary tags, secondary and
// supplementary alignments that the default exclude must drop.
// orph has flag 1 (paired) but neither 0x40 nor 0x80 set — that is the
// canonical orphan classification per upstream samtools fastq's `-0`.
const fastqSAM = `@HD	VN:1.6	SO:queryname
@SQ	SN:chr1	LN:1000
pa	99	chr1	100	60	5M	=	200	100	ACGTA	IIIII	NM:i:1
pa	147	chr1	200	60	5M	=	100	-100	TGCAT	IIIII	NM:i:0
pb	83	chr1	300	60	5M	=	400	100	ACGTA	IIIII
pb	163	chr1	400	60	5M	=	300	-100	TGCAT	IIIII
orph	1	chr1	500	60	5M	*	0	0	GGGGG	IIIII
sing	0	chr1	600	60	5M	*	0	0	CCCCC	IIIII
sec	256	chr1	700	60	5M	*	0	0	AAAAA	IIIII
sup	2048	chr1	800	60	5M	*	0	0	TTTTT	IIIII
`

// parseFastqRecords returns the set of records in a FASTQ string as a
// slice of [header, seq, +, qual] groups. It tolerates trailing newlines.
func parseFastqRecords(t *testing.T, s string) [][4]string {
	t.Helper()
	lines := strings.Split(strings.TrimRight(s, "\n"), "\n")
	if len(lines)%4 != 0 {
		t.Fatalf("FASTQ length %d not divisible by 4 — %q", len(lines), s)
	}
	out := make([][4]string, 0, len(lines)/4)
	for i := 0; i < len(lines); i += 4 {
		out = append(out, [4]string{lines[i], lines[i+1], lines[i+2], lines[i+3]})
	}
	return out
}

// readGzipFile returns the decompressed contents of a .gz file.
func readGzipFile(t *testing.T, path string) string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	gz, err := gzip.NewReader(f)
	if err != nil {
		t.Fatal(err)
	}
	defer gz.Close()
	b, err := io.ReadAll(gz)
	if err != nil {
		t.Fatal(err)
	}
	return string(b)
}

func TestFastqPairedSplit(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	orph := filepath.Join(dir, "orph.fq")
	sing := filepath.Join(dir, "sing.fq")
	counts, err := Fastq(strings.NewReader(fastqSAM), FastqOptions{
		Read1Path:     r1,
		Read2Path:     r2,
		OrphanPath:    orph,
		SingletonPath: sing,
	})
	if err != nil {
		t.Fatalf("Fastq: %v", err)
	}
	// Expected counts: 2 R1 (pa/1, pb/1), 2 R2 (pa/2, pb/2), 1 orphan
	// (orph), 1 singleton (sing). Secondary + supplementary excluded.
	if counts.Read1 != 2 {
		t.Errorf("Read1 count: got %d, want 2", counts.Read1)
	}
	if counts.Read2 != 2 {
		t.Errorf("Read2 count: got %d, want 2", counts.Read2)
	}
	if counts.Orphan != 1 {
		t.Errorf("Orphan count: got %d, want 1", counts.Orphan)
	}
	if counts.Singleton != 1 {
		t.Errorf("Singleton count: got %d, want 1", counts.Singleton)
	}

	r1Body, _ := os.ReadFile(r1)
	r1Recs := parseFastqRecords(t, string(r1Body))
	if len(r1Recs) != 2 {
		t.Fatalf("R1 has %d records, want 2", len(r1Recs))
	}
	// When -1 and -2 are both set (Read1Path + Read2Path), upstream
	// samtools drops the /1 /2 suffix unless -N (AlwaysAddSuffix) is on:
	// the separate output files already disambiguate mate identity.
	if r1Recs[0][0] != "@pa" {
		t.Errorf("R1[0] header: got %q, want @pa", r1Recs[0][0])
	}
	// pa R1 is forward strand (flag 99 = paired+proper+mate_rev+read1; no
	// 0x10), so SEQ unchanged.
	if r1Recs[0][1] != "ACGTA" {
		t.Errorf("R1[0] seq: got %q, want ACGTA", r1Recs[0][1])
	}

	r2Body, _ := os.ReadFile(r2)
	r2Recs := parseFastqRecords(t, string(r2Body))
	if len(r2Recs) != 2 {
		t.Fatalf("R2 has %d records, want 2", len(r2Recs))
	}
	// pa R2 is flag 147 = paired+proper+reverse+read2 → must reverse-
	// complement SEQ.
	if r2Recs[0][0] != "@pa" {
		t.Errorf("R2[0] header: got %q, want @pa (paired mode auto-drops suffix)", r2Recs[0][0])
	}
	if r2Recs[0][1] != "ATGCA" { // revcomp of TGCAT
		t.Errorf("R2[0] seq: got %q, want ATGCA (revcomp of TGCAT)", r2Recs[0][1])
	}
	// QUAL should be reversed too — "IIIII" reversed is still "IIIII"
	// since every char is identical, so just check it's non-empty.
	if r2Recs[0][3] != "IIIII" {
		t.Errorf("R2[0] qual: got %q, want IIIII", r2Recs[0][3])
	}

	orphBody, _ := os.ReadFile(orph)
	if !strings.Contains(string(orphBody), "@orph\n") {
		t.Errorf("orphan output missing @orph; got %q", orphBody)
	}

	singBody, _ := os.ReadFile(sing)
	if !strings.Contains(string(singBody), "@sing\n") {
		t.Errorf("singleton output missing @sing (no /1/2 suffix expected); got %q", singBody)
	}
}

func TestFastqInterleavedOutput(t *testing.T) {
	dir := t.TempDir()
	out := filepath.Join(dir, "all.fq")
	counts, err := Fastq(strings.NewReader(fastqSAM), FastqOptions{OutputPath: out})
	if err != nil {
		t.Fatal(err)
	}
	if counts.Output != 6 {
		t.Errorf("Output count: got %d, want 6 (4 paired + 1 orphan + 1 singleton)", counts.Output)
	}
	body, _ := os.ReadFile(out)
	recs := parseFastqRecords(t, string(body))
	if len(recs) != 6 {
		t.Fatalf("interleaved output: got %d records, want 6", len(recs))
	}
}

// TestFastqDefaultConsecutivePairOrder is the parity regression for the
// emission-order bug: upstream bam2fq groups CONSECUTIVE same-QNAME records
// and always flushes R1 before R2 (flush_rec writes best[1] then best[2]),
// even when the second mate appears FIRST in the file. Our default
// (interleaved-to-output) path must do the same.
func TestFastqDefaultConsecutivePairOrder(t *testing.T) {
	// pa's second mate (R2, 0x1|0x80=129) is written BEFORE its first mate
	// (R1, 0x1|0x40=65). Upstream still emits R1 (@pa/1) before R2 (@pa/2).
	const sam = `@HD	VN:1.6	SO:queryname
@SQ	SN:chr1	LN:1000
pa	129	chr1	200	60	5M	=	100	-100	TGCAT	IIIII
pa	65	chr1	100	60	5M	=	200	100	ACGTA	IIIII
`
	dir := t.TempDir()
	out := filepath.Join(dir, "all.fq")
	if _, err := Fastq(strings.NewReader(sam), FastqOptions{OutputPath: out}); err != nil {
		t.Fatal(err)
	}
	body, _ := os.ReadFile(out)
	recs := parseFastqRecords(t, string(body))
	if len(recs) != 2 {
		t.Fatalf("got %d records, want 2", len(recs))
	}
	if recs[0][0] != "@pa/1" || recs[1][0] != "@pa/2" {
		t.Errorf("consecutive pair order: got [%q, %q], want [@pa/1, @pa/2]",
			recs[0][0], recs[1][0])
	}
}

// TestFastqDefaultNonAdjacentMatesKeepFileOrder verifies that mates which are
// NOT consecutive are never joined: each emits as a singleton in its original
// file position (upstream groups only consecutive QNAMEs, so coordinate-sorted
// separated mates come out in file order, not paired/reordered).
func TestFastqDefaultNonAdjacentMatesKeepFileOrder(t *testing.T) {
	// File order: xa/1, xb/1, xa/2, xb/2 — mates separated by the other pair.
	const sam = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
xa	65	chr1	100	60	5M	=	500	100	ACGTA	IIIII
xb	65	chr1	150	60	5M	=	600	100	CCCCC	IIIII
xa	129	chr1	500	60	5M	=	100	-100	TGCAT	IIIII
xb	129	chr1	600	60	5M	=	150	-100	GGGGG	IIIII
`
	dir := t.TempDir()
	out := filepath.Join(dir, "all.fq")
	if _, err := Fastq(strings.NewReader(sam), FastqOptions{OutputPath: out}); err != nil {
		t.Fatal(err)
	}
	body, _ := os.ReadFile(out)
	recs := parseFastqRecords(t, string(body))
	want := []string{"@xa/1", "@xb/1", "@xa/2", "@xb/2"}
	if len(recs) != len(want) {
		t.Fatalf("got %d records, want %d", len(recs), len(want))
	}
	for i, w := range want {
		if recs[i][0] != w {
			t.Errorf("record %d: got %q, want %q", i, recs[i][0], w)
		}
	}
}

func TestFastqReverseStrandReverseComplement(t *testing.T) {
	// A single unpaired reverse-strand record. BAM stores SEQ as
	// revcomp(original), so emitting FASTQ must revcomp it back.
	samRev := `@HD	VN:1.6
@SQ	SN:chr1	LN:1000
r1	16	chr1	100	60	5M	*	0	0	ATGCA	HIJKL
`
	dir := t.TempDir()
	out := filepath.Join(dir, "out.fq")
	if _, err := Fastq(strings.NewReader(samRev), FastqOptions{OutputPath: out}); err != nil {
		t.Fatal(err)
	}
	body, _ := os.ReadFile(out)
	recs := parseFastqRecords(t, string(body))
	if len(recs) != 1 {
		t.Fatalf("got %d records, want 1", len(recs))
	}
	if recs[0][1] != "TGCAT" {
		t.Errorf("reverse strand seq: got %q, want TGCAT (revcomp of ATGCA)", recs[0][1])
	}
	if recs[0][3] != "LKJIH" {
		t.Errorf("reverse strand qual: got %q, want LKJIH (reversed)", recs[0][3])
	}
}

func TestFastqNoSuffix(t *testing.T) {
	dir := t.TempDir()
	out := filepath.Join(dir, "all.fq")
	if _, err := Fastq(strings.NewReader(fastqSAM), FastqOptions{
		OutputPath: out,
		NoSuffix:   true,
	}); err != nil {
		t.Fatal(err)
	}
	body, _ := os.ReadFile(out)
	if strings.Contains(string(body), "/1") || strings.Contains(string(body), "/2") {
		t.Errorf("NoSuffix output contains /1 or /2 suffix: %q", body)
	}
}

func TestFastqAlwaysAddSuffix(t *testing.T) {
	// If the qname already ends with /1, we should not double-suffix; but
	// with AlwaysAddSuffix=true we still append.
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:1000
r/1	64	chr1	100	60	5M	*	0	0	ACGTA	IIIII
`
	dir := t.TempDir()
	out := filepath.Join(dir, "all.fq")
	if _, err := Fastq(strings.NewReader(sam), FastqOptions{
		OutputPath:      out,
		AlwaysAddSuffix: true,
	}); err != nil {
		t.Fatal(err)
	}
	body, _ := os.ReadFile(out)
	if !strings.HasPrefix(string(body), "@r/1/1\n") {
		t.Errorf("AlwaysAddSuffix should double-suffix: got %q", body)
	}

	// Without AlwaysAddSuffix, the existing /1 suffix is preserved as-is.
	out2 := filepath.Join(dir, "all2.fq")
	if _, err := Fastq(strings.NewReader(sam), FastqOptions{OutputPath: out2}); err != nil {
		t.Fatal(err)
	}
	body2, _ := os.ReadFile(out2)
	if !strings.HasPrefix(string(body2), "@r/1\n") {
		t.Errorf("default mode should not double-suffix: got %q", body2)
	}
}

func TestFastqAddTags(t *testing.T) {
	dir := t.TempDir()
	out := filepath.Join(dir, "all.fq")
	if _, err := Fastq(strings.NewReader(fastqSAM), FastqOptions{
		OutputPath: out,
		AddTags:    []string{"NM"},
	}); err != nil {
		t.Fatal(err)
	}
	body, _ := os.ReadFile(out)
	// pa R1 has NM:i:1; header should contain it.
	if !strings.Contains(string(body), "@pa/1\tNM:i:1") {
		t.Errorf("AddTags NM should appear after qname: got %q", body)
	}
}

func TestFastqIncludeExcludeFlags(t *testing.T) {
	dir := t.TempDir()
	out := filepath.Join(dir, "all.fq")
	// Require 0x40 (read1) — should drop pa/2, pb/2, sing, etc.
	counts, err := Fastq(strings.NewReader(fastqSAM), FastqOptions{
		OutputPath:   out,
		IncludeFlags: 0x40,
	})
	if err != nil {
		t.Fatal(err)
	}
	// pa R1, pb R1 are the 0x40-set primaries. sec (256) is excluded by
	// default. So 2 records emitted.
	if counts.Output != 2 {
		t.Errorf("Output count with -f 0x40: got %d, want 2", counts.Output)
	}
}

func TestFastqGzipOutput(t *testing.T) {
	dir := t.TempDir()
	out := filepath.Join(dir, "all.fq.gz")
	if _, err := Fastq(strings.NewReader(fastqSAM), FastqOptions{
		OutputPath:    out,
		CompressLevel: 1,
	}); err != nil {
		t.Fatal(err)
	}
	content := readGzipFile(t, out)
	recs := parseFastqRecords(t, content)
	if len(recs) != 6 {
		t.Errorf("gzipped output: got %d records, want 6", len(recs))
	}
}

func TestFastqCoordinateSortedWarn(t *testing.T) {
	// Coordinate-sorted SAM in paired mode → counts.PairedCoordinateWarn
	// should be set and reads should fall through to the (here unset)
	// fallback. With no Output set, paired reads count as Dropped.
	sam := `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	99	chr1	100	60	5M	=	200	100	ACGTA	IIIII
r1	147	chr1	200	60	5M	=	100	-100	TGCAT	IIIII
`
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	counts, err := Fastq(strings.NewReader(sam), FastqOptions{
		Read1Path: r1,
		Read2Path: r2,
	})
	if err != nil {
		t.Fatal(err)
	}
	if !counts.PairedCoordinateWarn {
		t.Errorf("PairedCoordinateWarn should be set for coordinate-sorted paired input")
	}
}

func TestFastqUseOQ(t *testing.T) {
	// Record has Q="IIIII" but OQ:Z:!!!!!. With UseOQ, the output should
	// use the OQ value.
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:1000
r1	0	chr1	100	60	5M	*	0	0	ACGTA	IIIII	OQ:Z:!!!!!
`
	dir := t.TempDir()
	out := filepath.Join(dir, "out.fq")
	if _, err := Fastq(strings.NewReader(sam), FastqOptions{
		OutputPath: out,
		UseOQ:      true,
	}); err != nil {
		t.Fatal(err)
	}
	body, _ := os.ReadFile(out)
	recs := parseFastqRecords(t, string(body))
	if len(recs) != 1 {
		t.Fatalf("got %d records, want 1", len(recs))
	}
	if recs[0][3] != "!!!!!" {
		t.Errorf("UseOQ qual: got %q, want !!!!!", recs[0][3])
	}
}

func TestFastqMissingQual(t *testing.T) {
	// SAM with QUAL="*" → output should fall back to "!" per length.
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:1000
r1	0	chr1	100	60	5M	*	0	0	ACGTA	*
`
	dir := t.TempDir()
	out := filepath.Join(dir, "out.fq")
	if _, err := Fastq(strings.NewReader(sam), FastqOptions{OutputPath: out}); err != nil {
		t.Fatal(err)
	}
	body, _ := os.ReadFile(out)
	recs := parseFastqRecords(t, string(body))
	if len(recs) != 1 {
		t.Fatalf("got %d records, want 1", len(recs))
	}
	if recs[0][3] != "!!!!!" {
		t.Errorf("missing qual fallback: got %q, want !!!!!", recs[0][3])
	}
}

func TestParseAddTags(t *testing.T) {
	cases := []struct {
		in   string
		want []string
	}{
		{"", nil},
		{"NM", []string{"NM"}},
		{"NM,MD,RG", []string{"NM", "MD", "RG"}},
		{" NM , MD ", []string{"NM", "MD"}},
		{",,", nil},
	}
	for _, tc := range cases {
		got := ParseAddTags(tc.in)
		if len(got) != len(tc.want) {
			t.Errorf("ParseAddTags(%q): got %v, want %v", tc.in, got, tc.want)
			continue
		}
		for i := range got {
			if got[i] != tc.want[i] {
				t.Errorf("ParseAddTags(%q)[%d]: got %q, want %q", tc.in, i, got[i], tc.want[i])
			}
		}
	}
}

func TestParseCompressLevel(t *testing.T) {
	cases := []struct {
		in      string
		want    int
		wantErr bool
	}{
		{"", 0, false},
		{"0", 0, false},
		{"9", 9, false},
		{"10", 0, true},
		{"-1", 0, true},
		{"abc", 0, true},
	}
	for _, tc := range cases {
		got, err := ParseCompressLevel(tc.in)
		if (err != nil) != tc.wantErr {
			t.Errorf("ParseCompressLevel(%q): err=%v, wantErr=%v", tc.in, err, tc.wantErr)
		}
		if got != tc.want {
			t.Errorf("ParseCompressLevel(%q): got %d, want %d", tc.in, got, tc.want)
		}
	}
}

func TestReverseComplement(t *testing.T) {
	cases := []struct{ in, want string }{
		{"ACGT", "ACGT"},
		{"AAAA", "TTTT"},
		{"NNNN", "NNNN"},
		{"acgt", "acgt"},
		// RYKM → complement YRMK → reverse → KMRY.
		{"RYKM", "KMRY"},
		// BDHV → complement VHDB → reverse → BDHV.
		{"BDHV", "BDHV"},
		// SW → complement SW → reverse → WS.
		{"SW", "WS"},
		{"X", "X"}, // unknown passthrough
	}
	for _, tc := range cases {
		if got := reverseComplement(tc.in); got != tc.want {
			t.Errorf("reverseComplement(%q) = %q, want %q", tc.in, got, tc.want)
		}
	}
}

func TestErrNoFastqOutputSentinel(t *testing.T) {
	// Existence smoke test for the sentinel error (used by CLI driver).
	if ErrNoFastqOutput == nil {
		t.Fatal("ErrNoFastqOutput must be non-nil")
	}
}

// TestFastqPairedPartialFallback exercises the partial-output fallback
// paths: a paired-mode invocation with only -1 (no -2/-0/-s/-o) — R2,
// singletons, and orphans should fall through to "dropped".
func TestFastqPairedPartialFallback(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	counts, err := Fastq(strings.NewReader(fastqSAM), FastqOptions{Read1Path: r1})
	if err != nil {
		t.Fatal(err)
	}
	if counts.Read1 != 2 {
		t.Errorf("Read1 count: got %d, want 2", counts.Read1)
	}
	// R2 records (pa/2, pb/2), singleton (sing), orphan (orph) — all
	// have no configured sink → dropped.
	if counts.Dropped < 4 {
		t.Errorf("Dropped count: got %d, want at least 4", counts.Dropped)
	}
}

// TestFastqPairedWithOutputFallback exercises the path where the paired
// configuration has only -1 but also -o — R2/orphan/singleton fall through
// to the Output sink.
func TestFastqPairedWithOutputFallback(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	out := filepath.Join(dir, "out.fq")
	counts, err := Fastq(strings.NewReader(fastqSAM), FastqOptions{
		Read1Path:  r1,
		OutputPath: out,
	})
	if err != nil {
		t.Fatal(err)
	}
	if counts.Read1 != 2 {
		t.Errorf("Read1 count: got %d, want 2", counts.Read1)
	}
	// R2 (2) + orph (1) + sing (1) = 4 records go to output.
	if counts.Output != 4 {
		t.Errorf("Output count: got %d, want 4", counts.Output)
	}
}

// TestFastqExcludeFlagsAll covers the -G/UseExcludeAll branch.
func TestFastqExcludeFlagsAll(t *testing.T) {
	dir := t.TempDir()
	out := filepath.Join(dir, "out.fq")
	// Records r1 and r2 below both have FLAG=99 (paired+proper+mate_rev+read1).
	// Setting `-G 0x60` (mate_rev + read1) should drop both because the
	// combined bits are exactly set.
	sam := `@HD	VN:1.6	SO:queryname
@SQ	SN:chr1	LN:1000
r1	99	chr1	100	60	5M	*	0	0	ACGTA	IIIII
r2	0	chr1	200	60	5M	*	0	0	ACGTA	IIIII
`
	counts, err := Fastq(strings.NewReader(sam), FastqOptions{
		OutputPath:      out,
		ExcludeFlagsAll: 0x60,
		UseExcludeAll:   true,
	})
	if err != nil {
		t.Fatal(err)
	}
	if counts.Output != 1 {
		t.Errorf("with -G 0x60: got %d, want 1 (r2 only)", counts.Output)
	}
}

// TestFastqSingletonOnly covers the path where only -s is configured (no
// -1/-2/-o/-0): paired reads have nowhere to go and are dropped.
func TestFastqSingletonOnly(t *testing.T) {
	dir := t.TempDir()
	sing := filepath.Join(dir, "s.fq")
	counts, err := Fastq(strings.NewReader(fastqSAM), FastqOptions{SingletonPath: sing})
	if err != nil {
		t.Fatal(err)
	}
	// `sing` is the only unpaired record → Singleton count = 1, rest dropped.
	if counts.Singleton != 1 {
		t.Errorf("Singleton count: got %d, want 1", counts.Singleton)
	}
}

// loneMateSAM has two paired reads whose mates are absent from the file: a
// lone READ1 (flag 0x1|0x40=65) and a lone READ2 (flag 0x1|0x80=129). Both
// are forward strand so SEQ is written verbatim.
const loneMateSAM = `@HD	VN:1.6	SO:queryname
@SQ	SN:chr1	LN:1000
loneA	65	chr1	100	60	5M	*	0	0	ACGTA	IIIII
loneB	129	chr1	200	60	5M	*	0	0	CCCCC	IIIII
`

// TestFastqLoneMateToReadFiles is the parity test for the lone-mate routing
// fix: in the -1/-2 split path, a paired read whose mate is absent must be
// written to its R1/R2 file (upstream bam2fq flush_rec writes fpr[1]/fpr[2]
// when no -s file is given), not dropped into a nil sink.
func TestFastqLoneMateToReadFiles(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	counts, err := Fastq(strings.NewReader(loneMateSAM), FastqOptions{
		Read1Path: r1,
		Read2Path: r2,
	})
	if err != nil {
		t.Fatalf("Fastq: %v", err)
	}
	if counts.Read1 != 1 {
		t.Errorf("Read1 count: got %d, want 1 (lone R1 written to R1 file)", counts.Read1)
	}
	if counts.Read2 != 1 {
		t.Errorf("Read2 count: got %d, want 1 (lone R2 written to R2 file)", counts.Read2)
	}
	if counts.Dropped != 0 {
		t.Errorf("Dropped count: got %d, want 0 (lone mates must not be dropped)", counts.Dropped)
	}
	if body, _ := os.ReadFile(r1); !strings.Contains(string(body), "@loneA\nACGTA\n") {
		t.Errorf("R1 file missing lone READ1; got %q", body)
	}
	if body, _ := os.ReadFile(r2); !strings.Contains(string(body), "@loneB\nCCCCC\n") {
		t.Errorf("R2 file missing lone READ2; got %q", body)
	}
}

// TestFastqLoneMateToSingleton verifies that when a singleton (-s) file is
// configured, lone mates are routed there instead of the R1/R2 files —
// matching upstream bam2fq flush_rec's fpse branch.
func TestFastqLoneMateToSingleton(t *testing.T) {
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	sing := filepath.Join(dir, "s.fq")
	counts, err := Fastq(strings.NewReader(loneMateSAM), FastqOptions{
		Read1Path:     r1,
		Read2Path:     r2,
		SingletonPath: sing,
	})
	if err != nil {
		t.Fatalf("Fastq: %v", err)
	}
	if counts.Singleton != 2 {
		t.Errorf("Singleton count: got %d, want 2 (both lone mates to -s)", counts.Singleton)
	}
	if counts.Read1 != 0 || counts.Read2 != 0 {
		t.Errorf("Read1/Read2 counts: got %d/%d, want 0/0 (lone mates go to -s)", counts.Read1, counts.Read2)
	}
	body, _ := os.ReadFile(sing)
	if !strings.Contains(string(body), "@loneA\n") || !strings.Contains(string(body), "@loneB\n") {
		t.Errorf("singleton file missing lone mates; got %q", body)
	}
}

// TestFastqInputError exercises the error path when the reader fails on
// header parse.
func TestFastqInputError(t *testing.T) {
	// A SAM stream with a malformed body line — the header reads fine
	// but the body parse will surface an error.
	bad := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nthis-is-not-a-valid-sam-record\n"
	_, err := Fastq(strings.NewReader(bad), FastqOptions{OutputPath: "-"})
	if err == nil {
		t.Fatal("expected an error from malformed SAM body")
	}
}

// TestFastqOpenFastqOutputStdoutSink covers the "-" path that returns the
// stdout sink. We don't actually write to stdout in tests; this just
// ensures the constructor returns without error.
func TestFastqOpenFastqOutputStdoutSink(t *testing.T) {
	s, err := openFastqOutput("-", 0)
	if err != nil {
		t.Fatalf("openFastqOutput(-): %v", err)
	}
	if s.w == nil || s.bw == nil {
		t.Fatal("openFastqOutput(-): nil sink fields")
	}
	// Closer is nopCloseStdout — Close should return nil.
	if err := s.w.Close(); err != nil {
		t.Errorf("nopCloseStdout.Close: %v", err)
	}
}

// TestFastqOpenFastqOutputUnopenable covers the error path for a path
// that cannot be created.
func TestFastqOpenFastqOutputUnopenable(t *testing.T) {
	// A path that points into a non-existent directory will fail Create.
	_, err := openFastqOutput("/nonexistent-directory-xyz123/out.fq", 0)
	if err == nil {
		t.Fatal("expected error for unopenable path")
	}
}
