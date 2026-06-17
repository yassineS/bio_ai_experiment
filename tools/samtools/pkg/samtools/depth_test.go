package samtools

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// depthSAM is a small SAM file used by the depth tests. It places 5M reads
// at known coordinates so depth at each position is computable by hand.
//
//	r1: chr1:10..14 (5M)
//	r2: chr1:12..16 (5M) — overlaps r1 over chr1:12..14
//	r3: chr1:49..51 (3M run inside 1S3M1S) — soft clips don't consume ref.
//	r4: chr1:60..63 (2M2I2M) — insertion does not advance ref; covers
//	    chr1:60..61 and chr1:62..63 (positions 64+ NOT covered because
//	    the M3 run is only 2 long).
//	r5: chr2:5..7 (3M) — only contributes to chr2.
//	r6: chr1:200..204 (5M) but MAPQ=0 — for the MAPQ filter test.
//	r7: chr1:300..304 (5M) but FLAG=0x4 unmapped — for default exclude.
//	r8: chr1:400..404 (5M) with low quality scores ("!!!!!" = Phred 0)
//	    for the BaseQ filter test.
const depthSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:500
@SQ	SN:chr2	LN:20
r1	0	chr1	10	60	5M	*	0	0	ACGTA	IIIII
r2	0	chr1	12	60	5M	*	0	0	TGCAT	IIIII
r3	0	chr1	49	60	1S3M1S	*	0	0	NACGN	IIIII
r4	0	chr1	60	60	2M2I2M	*	0	0	ACTTGA	IIIIII
r5	0	chr2	5	60	3M	*	0	0	ACG	III
r6	0	chr1	200	0	5M	*	0	0	ACGTA	IIIII
r7	4	chr1	300	60	5M	*	0	0	ACGTA	IIIII
r8	0	chr1	400	60	5M	*	0	0	ACGTA	!!!!!
`

// parseDepthOut converts the tab-delimited depth output into a map keyed by
// "chrom:pos" for easy assertions.
func parseDepthOut(s string) map[string][]string {
	out := map[string][]string{}
	for _, line := range strings.Split(strings.TrimRight(s, "\n"), "\n") {
		if line == "" {
			continue
		}
		f := strings.Split(line, "\t")
		key := f[0] + ":" + f[1]
		out[key] = f[2:]
	}
	return out
}

func runDepth(t *testing.T, samInputs []string, opts DepthOptions) string {
	t.Helper()
	in := make([]io.Reader, len(samInputs))
	for i, s := range samInputs {
		in[i] = strings.NewReader(s)
	}
	var buf bytes.Buffer
	if err := Depth(in, &buf, opts); err != nil {
		t.Fatalf("Depth: %v", err)
	}
	return buf.String()
}

func TestDepthBasicCoverage(t *testing.T) {
	got := runDepth(t, []string{depthSAM}, DepthOptions{
		ExcludeFlags: 0x4, // default unmapped exclude
	})
	pos := parseDepthOut(got)

	// Spot-check coverage at chr1:10..16, chr1:50..52, chr1:60..63.
	cases := map[string]string{
		"chr1:10": "1",
		"chr1:11": "1",
		"chr1:12": "2",
		"chr1:13": "2",
		"chr1:14": "2",
		"chr1:15": "1",
		"chr1:16": "1",
		"chr1:49": "1",
		"chr1:50": "1",
		"chr1:51": "1",
		"chr1:60": "1",
		"chr1:61": "1",
		"chr1:62": "1",
		"chr1:63": "1",
		"chr2:5":  "1",
		"chr2:6":  "1",
		"chr2:7":  "1",
	}
	for k, want := range cases {
		gotDepths := pos[k]
		if len(gotDepths) != 1 || gotDepths[0] != want {
			t.Errorf("position %s: got %v, want depth %s", k, gotDepths, want)
		}
	}

	// Should NOT contain zero-depth positions by default.
	if _, ok := pos["chr1:1"]; ok {
		t.Errorf("default mode unexpectedly emitted chr1:1 (zero depth)")
	}
	// r6 has MAPQ=0 and ExcludeFlags=0x4 doesn't drop it — but we filter
	// it via MAPQ filter in another test.  Here, r6 should be counted.
	if d := pos["chr1:200"]; len(d) != 1 || d[0] != "1" {
		t.Errorf("chr1:200 (r6): got %v, want 1", d)
	}
	// r7 is unmapped (0x4) — should be excluded by default.
	if _, ok := pos["chr1:300"]; ok {
		t.Errorf("unmapped read r7 unexpectedly produced depth at chr1:300")
	}
}

func TestDepthCigarOpsExcludesInsertSoftHard(t *testing.T) {
	// CIGAR 1S3M1S — the 1S clips shouldn't count.
	// Read starts at SAM POS=49 (the M block is at chr1:50..52). Position
	// chr1:49 should NOT be in output because S doesn't consume ref.
	got := runDepth(t, []string{depthSAM}, DepthOptions{ExcludeFlags: 0x4})
	pos := parseDepthOut(got)
	// r3 has CIGAR 1S3M1S at POS=49 → covers chr1:49,50,51 (the M run).
	// The 1S clips before and after must NOT extend ref coverage.
	if d := pos["chr1:49"]; len(d) != 1 || d[0] != "1" {
		t.Errorf("chr1:49 (r3 M run): got %v, want depth 1", d)
	}
	if d := pos["chr1:51"]; len(d) != 1 || d[0] != "1" {
		t.Errorf("chr1:51 (r3 M run end): got %v, want depth 1", d)
	}
	// chr1:52 must NOT have depth — the trailing S clip should not advance
	// the ref position.
	if _, ok := pos["chr1:52"]; ok {
		t.Errorf("chr1:52 unexpectedly produced depth: trailing S clip should NOT count")
	}

	// r4: CIGAR 2M2I2M, POS=60. Insertion at ref position 62 should NOT
	// advance ref. So coverage is chr1:60..63 (4 positions).
	for p := 60; p <= 63; p++ {
		key := "chr1:" + itoa(p)
		if d := pos[key]; len(d) != 1 || d[0] != "1" {
			t.Errorf("%s (r4 2M2I2M): got %v, want depth 1", key, d)
		}
	}
	if _, ok := pos["chr1:64"]; ok {
		t.Errorf("chr1:64 should NOT have depth (2M2I2M only spans 4 ref bases)")
	}
}

func TestDepthAllPositions(t *testing.T) {
	// `-r` to restrict to a small range, with `-a` to emit zeros inside it.
	got := runDepth(t, []string{depthSAM}, DepthOptions{
		ExcludeFlags: 0x4,
		AllPositions: true,
		Regions:      []string{"chr1:8-14"},
	})
	pos := parseDepthOut(got)
	// Positions 8 and 9 are inside [8,14] but have zero depth — should still
	// be emitted because AllPositions=true.
	for _, p := range []int{8, 9} {
		key := "chr1:" + itoa(p)
		if d := pos[key]; len(d) != 1 || d[0] != "0" {
			t.Errorf("%s with -a: got %v, want depth 0", key, d)
		}
	}
	// Position 10 is the first read; should be depth 1.
	if d := pos["chr1:10"]; len(d) != 1 || d[0] != "1" {
		t.Errorf("chr1:10 with -a -r chr1:8-14: got %v, want 1", d)
	}
}

func TestDepthAllTrans(t *testing.T) {
	// `-A` emits every position of every reference, even where no reads
	// land. chr2 is 20bp, only 5..7 are covered.
	got := runDepth(t, []string{depthSAM}, DepthOptions{
		ExcludeFlags:      0x4,
		AllTransPositions: true,
	})
	pos := parseDepthOut(got)
	// chr2:1..4 should be 0; chr2:5..7 should be 1; chr2:8..20 should be 0.
	for _, p := range []int{1, 2, 3, 4, 8, 9, 10, 11, 20} {
		key := "chr2:" + itoa(p)
		if d := pos[key]; len(d) != 1 || d[0] != "0" {
			t.Errorf("%s with -A: got %v, want depth 0", key, d)
		}
	}
}

func TestDepthMinMAPQ(t *testing.T) {
	// r6 has MAPQ=0; filter it out with `-Q 1` (the mapping-quality knob).
	got := runDepth(t, []string{depthSAM}, DepthOptions{
		ExcludeFlags: 0x4,
		MinMAPQ:      1,
	})
	pos := parseDepthOut(got)
	if _, ok := pos["chr1:200"]; ok {
		t.Errorf("chr1:200 (r6 MAPQ=0) should be filtered by -q 1")
	}
	// r1 is still at MAPQ=60 → should still be there.
	if d := pos["chr1:10"]; len(d) != 1 || d[0] != "1" {
		t.Errorf("chr1:10: got %v, want 1", d)
	}
}

func TestDepthMinBaseQ(t *testing.T) {
	// r8 has Phred 0 quality everywhere (ASCII '!'). With `-q 1` (the base-
	// quality knob) every base fails the filter, so the depth at chr1:400..404
	// is 0. But those positions are still INSIDE r8's read span, so upstream
	// samtools depth emits a row for each of them with depth 0 (see
	// reference_code/samtools/bam2depth.c: the position prints while i is
	// within the read's [pos, bam_endpos) span; the qual filter only zeroes
	// the count). Our port must do the same.
	got := runDepth(t, []string{depthSAM}, DepthOptions{
		ExcludeFlags: 0x4,
		MinBaseQ:     1,
	})
	pos := parseDepthOut(got)
	for p := 400; p <= 404; p++ {
		key := "chr1:" + itoa(p)
		if d := pos[key]; len(d) != 1 || d[0] != "0" {
			t.Errorf("%s with -q 1 (r8 quality is 0): got %v, want depth 0 (in-span zero)", key, d)
		}
	}
	// Position 405 is one past r8's span (5M from POS 400 covers 400..404)
	// and is covered by no other read, so it must NOT be emitted.
	if _, ok := pos["chr1:405"]; ok {
		t.Errorf("chr1:405 is outside r8's span and should not be emitted")
	}
	// r1 (chr1:10..14, full-quality 'IIIII') is unaffected by the filter.
	if d := pos["chr1:10"]; len(d) != 1 || d[0] != "1" {
		t.Errorf("chr1:10: got %v, want 1", d)
	}
}

// TestDepthInteriorZeroBaseQFilter reproduces the verified upstream parity
// example for `samtools depth -q 10`: a position whose only covering base
// fails the base-quality filter is still emitted with depth 0 because it is
// interior to the read's covered span, bracketed by passing positions on
// either side. Upstream prints chr1\t8\t0 between chr1\t7\t1 and chr1\t9\t1;
// our port used to skip position 8 entirely.
func TestDepthInteriorZeroBaseQFilter(t *testing.T) {
	// One read at chr1:6..10 (5M). The base at the 3rd position (chr1:8) has
	// quality '+' = Phred 10; all others are 'I' = Phred 40. With -q 11 only
	// the chr1:8 base fails, leaving an interior depth-0 hole at chr1:8 that
	// is still inside the read span chr1:6..10.
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:100
r1	0	chr1	6	60	5M	*	0	0	ACGTA	II+II
`
	got := runDepth(t, []string{sam}, DepthOptions{
		ExcludeFlags: 0x4,
		MinBaseQ:     11,
	})
	pos := parseDepthOut(got)
	want := map[string]string{
		"chr1:6":  "1",
		"chr1:7":  "1",
		"chr1:8":  "0", // the only base here is Phred 10 < 11 -> filtered, but in span
		"chr1:9":  "1",
		"chr1:10": "1",
	}
	for k, w := range want {
		if d := pos[k]; len(d) != 1 || d[0] != w {
			t.Errorf("%s: got %v, want depth %s", k, d, w)
		}
	}
	// Leading/trailing uncovered positions are NOT emitted (no -a).
	if _, ok := pos["chr1:5"]; ok {
		t.Errorf("chr1:5 (before the read span) should not be emitted")
	}
	if _, ok := pos["chr1:11"]; ok {
		t.Errorf("chr1:11 (after the read span) should not be emitted")
	}
}

// TestDepthInteriorZeroDeletion verifies a deletion (CIGAR D) inside a read
// prints depth 0 at the deleted reference positions while flanking matches
// print their real depth — they are all inside the read's covered span.
func TestDepthInteriorZeroDeletion(t *testing.T) {
	// CIGAR 2M2D2M at POS=20 spans chr1:20..25 (2 match, 2 deleted, 2 match).
	// The deleted positions chr1:22,23 carry no base => depth 0, but they are
	// interior to the read span so they are emitted.
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:100
r1	0	chr1	20	60	2M2D2M	*	0	0	ACGT	IIII
`
	got := runDepth(t, []string{sam}, DepthOptions{ExcludeFlags: 0x4})
	pos := parseDepthOut(got)
	want := map[string]string{
		"chr1:20": "1",
		"chr1:21": "1",
		"chr1:22": "0", // deleted
		"chr1:23": "0", // deleted
		"chr1:24": "1",
		"chr1:25": "1",
	}
	for k, w := range want {
		if d := pos[k]; len(d) != 1 || d[0] != w {
			t.Errorf("%s: got %v, want depth %s", k, d, w)
		}
	}
}

func TestDepthMaxDepthCap(t *testing.T) {
	// Cap reported depth at 1; chr1:12..14 has two-read overlap so should
	// be reported as 1 not 2.
	got := runDepth(t, []string{depthSAM}, DepthOptions{
		ExcludeFlags: 0x4,
		MaxDepth:     1,
	})
	pos := parseDepthOut(got)
	for p := 12; p <= 14; p++ {
		key := "chr1:" + itoa(p)
		if d := pos[key]; len(d) != 1 || d[0] != "1" {
			t.Errorf("%s with -d 1: got %v, want capped 1", key, d)
		}
	}
}

func TestDepthMultiBAMParallel(t *testing.T) {
	// Two SAM inputs with disjoint reads — output should have two depth
	// columns and each column should equal the depth in its own input.
	samA := `@HD	VN:1.6
@SQ	SN:chr1	LN:100
ra	0	chr1	10	60	3M	*	0	0	ACG	III
`
	samB := `@HD	VN:1.6
@SQ	SN:chr1	LN:100
rb	0	chr1	11	60	3M	*	0	0	ACG	III
`
	got := runDepth(t, []string{samA, samB}, DepthOptions{ExcludeFlags: 0x4})
	pos := parseDepthOut(got)
	// chr1:10 should be (1, 0); chr1:11 should be (1, 1); chr1:13 should
	// be (0 — wait it would skip the line entirely if both zero). chr1:13
	// is depth 0 in samA, depth 1 in samB → emitted.
	if d := pos["chr1:10"]; len(d) != 2 || d[0] != "1" || d[1] != "0" {
		t.Errorf("chr1:10: got %v, want [1 0]", d)
	}
	if d := pos["chr1:11"]; len(d) != 2 || d[0] != "1" || d[1] != "1" {
		t.Errorf("chr1:11: got %v, want [1 1]", d)
	}
	if d := pos["chr1:12"]; len(d) != 2 || d[0] != "1" || d[1] != "1" {
		t.Errorf("chr1:12: got %v, want [1 1]", d)
	}
	if d := pos["chr1:13"]; len(d) != 2 || d[0] != "0" || d[1] != "1" {
		t.Errorf("chr1:13: got %v, want [0 1]", d)
	}
}

func TestDepthIncludeFlagsFilter(t *testing.T) {
	// Use `-f 1` (must be paired). r1..r5 don't have 0x1 set → should be
	// dropped, leaving zero coverage everywhere.
	got := runDepth(t, []string{depthSAM}, DepthOptions{
		ExcludeFlags: 0x4,
		IncludeFlags: 0x1,
	})
	pos := parseDepthOut(got)
	if len(pos) != 0 {
		t.Errorf("with `-f 1`, expected no output, got %d positions", len(pos))
	}
}

func TestDepthRegionsBED(t *testing.T) {
	dir := t.TempDir()
	bedPath := filepath.Join(dir, "regions.bed")
	if err := os.WriteFile(bedPath, []byte("chr1\t9\t13\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	got := runDepth(t, []string{depthSAM}, DepthOptions{
		ExcludeFlags: 0x4,
		BedPath:      bedPath,
		AllPositions: true,
	})
	pos := parseDepthOut(got)
	// BED region is chr1:9..13 (0-based half-open → 1-based 10..13).
	// chr1:10 depth 1, chr1:11 depth 1, chr1:12 depth 2, chr1:13 depth 2.
	for k, want := range map[string]string{
		"chr1:10": "1", "chr1:11": "1", "chr1:12": "2", "chr1:13": "2",
	} {
		if d := pos[k]; len(d) != 1 || d[0] != want {
			t.Errorf("%s: got %v, want %s", k, d, want)
		}
	}
	// Outside the BED region should NOT appear.
	if _, ok := pos["chr1:14"]; ok {
		t.Errorf("chr1:14 should not appear (outside BED)")
	}
}

func TestDepthMinReadLen(t *testing.T) {
	// r5 has query length 3; with -l 4 it should be filtered.
	got := runDepth(t, []string{depthSAM}, DepthOptions{
		ExcludeFlags: 0x4,
		MinReadLen:   4,
	})
	pos := parseDepthOut(got)
	for _, p := range []int{5, 6, 7} {
		key := "chr2:" + itoa(p)
		if _, ok := pos[key]; ok {
			t.Errorf("%s should be filtered by -l 4 (r5 is 3M)", key)
		}
	}
}

func TestDepthBadHeader(t *testing.T) {
	a := `@HD	VN:1.6
@SQ	SN:chr1	LN:100
r1	0	chr1	1	60	3M	*	0	0	ACG	III
`
	b := `@HD	VN:1.6
@SQ	SN:chr2	LN:100
r1	0	chr2	1	60	3M	*	0	0	ACG	III
`
	var buf bytes.Buffer
	err := Depth([]io.Reader{strings.NewReader(a), strings.NewReader(b)}, &buf, DepthOptions{})
	if err == nil || !strings.Contains(err.Error(), "@SQ ordering") {
		t.Fatalf("expected @SQ ordering error, got %v", err)
	}
}

func TestDepthNoInputs(t *testing.T) {
	var buf bytes.Buffer
	if err := Depth(nil, &buf, DepthOptions{}); err == nil {
		t.Fatal("expected error from empty input list")
	}
}

func TestMergeIntervals(t *testing.T) {
	in := [][2]int{{10, 20}, {15, 25}, {30, 40}, {35, 45}, {50, 55}}
	want := [][2]int{{10, 25}, {30, 45}, {50, 55}}
	got := mergeIntervals(in)
	if len(got) != len(want) {
		t.Fatalf("mergeIntervals: got %v, want %v", got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("interval %d: got %v, want %v", i, got[i], want[i])
		}
	}
}

func itoa(i int) string {
	if i == 0 {
		return "0"
	}
	neg := i < 0
	if neg {
		i = -i
	}
	var buf [20]byte
	n := len(buf)
	for i > 0 {
		n--
		buf[n] = byte('0' + i%10)
		i /= 10
	}
	if neg {
		n--
		buf[n] = '-'
	}
	return string(buf[n:])
}
