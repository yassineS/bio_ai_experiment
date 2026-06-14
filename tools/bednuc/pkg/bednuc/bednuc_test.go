package bednuc

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// writeTempFile writes content to a fresh file under dir and returns the path.
func writeTempFile(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", p, err)
	}
	return p
}

// runIt is the standard test harness: spool a FASTA + a BED to disk in a
// temp directory and call Run with the given options. Returns stdout and
// the warning buffer.
func runIt(t *testing.T, fasta, bed string, opts Options) (string, string) {
	t.Helper()
	dir := t.TempDir()
	faPath := writeTempFile(t, dir, "ref.fa", fasta)
	bedPath := writeTempFile(t, dir, "in.bed", bed)
	bedR, err := os.Open(bedPath)
	if err != nil {
		t.Fatalf("open bed: %v", err)
	}
	defer bedR.Close()
	var out, warn bytes.Buffer
	if _, err := Run(bedR, faPath, &out, &warn, opts); err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	return out.String(), warn.String()
}

// TestProfile_HandCounts validates the counter against hand-computed counts.
func TestProfile_HandCounts(t *testing.T) {
	// ACGTAN: A=2, C=1, G=1, T=1, N=1, oth=0.
	c := Profile([]byte("ACGTAN"), "", false)
	if c.A != 2 || c.C != 1 || c.G != 1 || c.T != 1 || c.N != 1 || c.Other != 0 || c.SeqLen != 6 {
		t.Errorf("unexpected counts: %+v", c)
	}
	// Lowercase + Us count too: 'a'=A, 'u'=T.
	c2 := Profile([]byte("aauu"), "", false)
	if c2.A != 2 || c2.T != 2 {
		t.Errorf("aauu counts: %+v", c2)
	}
	// 'X' is "other".
	c3 := Profile([]byte("XX"), "", false)
	if c3.Other != 2 || c3.SeqLen != 2 {
		t.Errorf("XX counts: %+v", c3)
	}
}

// TestCountPattern_Overlapping verifies upstream's overlapping substring
// semantics: 'AAA' in 'AAAAA' = 3 hits (positions 0,1,2).
func TestCountPattern_Overlapping(t *testing.T) {
	if n := countPattern([]byte("AAAAA"), "AAA", false); n != 3 {
		t.Errorf("AAA in AAAAA: want 3, got %d", n)
	}
	if n := countPattern([]byte("ACGTACGT"), "CGT", false); n != 2 {
		t.Errorf("CGT in ACGTACGT: want 2, got %d", n)
	}
	// Case-sensitive default: lowercase pattern misses uppercase haystack.
	if n := countPattern([]byte("AAAAA"), "aaa", false); n != 0 {
		t.Errorf("aaa in AAAAA case-sensitive: want 0, got %d", n)
	}
	// Case-insensitive flips that on.
	if n := countPattern([]byte("AAAAA"), "aaa", true); n != 3 {
		t.Errorf("aaa in AAAAA case-insensitive: want 3, got %d", n)
	}
}

// TestReverseComplement_Basic checks A<->T / C<->G and IUPAC fall-through.
func TestReverseComplement_Basic(t *testing.T) {
	got := string(ReverseComplement([]byte("ACGTN")))
	if got != "NACGT" {
		t.Errorf("rc(ACGTN): want NACGT, got %q", got)
	}
	// Lowercase preserved.
	got2 := string(ReverseComplement([]byte("acgt")))
	if got2 != "acgt" {
		t.Errorf("rc(acgt): want acgt, got %q", got2)
	}
	// Unknown char passes through.
	got3 := string(ReverseComplement([]byte("AZ")))
	if got3 != "ZT" {
		t.Errorf("rc(AZ): want ZT, got %q", got3)
	}
}

// TestComplement_IUPAC exercises every branch of the IUPAC table so the
// switch stays covered for downstream users that touch it directly.
func TestComplement_IUPAC(t *testing.T) {
	pairs := map[byte]byte{
		'A': 'T', 'a': 't', 'C': 'G', 'c': 'g',
		'G': 'C', 'g': 'c', 'T': 'A', 't': 'a',
		'U': 'A', 'u': 'a',
		'R': 'Y', 'r': 'y', 'Y': 'R', 'y': 'r',
		'K': 'M', 'k': 'm', 'M': 'K', 'm': 'k',
		'B': 'V', 'b': 'v', 'V': 'B', 'v': 'b',
		'D': 'H', 'd': 'h', 'H': 'D', 'h': 'd',
		'N': 'N', 'n': 'n', 'S': 'S', 's': 's',
		'W': 'W', 'w': 'w',
		'?': '?', '5': '5', // fall-through
	}
	for in, want := range pairs {
		if got := complement(in); got != want {
			t.Errorf("complement(%c)=%c, want %c", in, got, want)
		}
	}
}

// TestRun_DefaultColumns is a hand-computed end-to-end check.
//
// FASTA chr1 = ACGTACGTAC (10bp). Interval [0,5) -> ACGTA:
//
//	A=2, C=1, G=1, T=1, N=0, oth=0, len=5, %AT=3/5=0.600000, %GC=2/5=0.400000.
func TestRun_DefaultColumns(t *testing.T) {
	fasta := ">chr1\nACGTACGTAC\n"
	bed := "chr1\t0\t5\n"
	out, _ := runIt(t, fasta, bed, Options{})
	if !strings.HasPrefix(out, "#1_usercol\t2_usercol\t3_usercol\t4_pct_at") {
		t.Errorf("unexpected header: %q", strings.SplitN(out, "\n", 2)[0])
	}
	want := "chr1\t0\t5\t0.600000\t0.400000\t2\t1\t1\t1\t0\t0\t5\n"
	if !strings.Contains(out, want) {
		t.Errorf("expected row %q in output:\n%s", want, out)
	}
}

// TestRun_StrandRevComp confirms that -s reverse-complements before counting
// (counts themselves are symmetric for A/T and C/G so the only observable
// effect is on the -seq output, which we exercise here too).
func TestRun_StrandRevComp(t *testing.T) {
	fasta := ">chr1\nAAAACCCCGG\n" // 10 bp
	// Interval [0,4) = AAAA on '+' strand; on '-' strand RC = TTTT.
	bed := "chr1\t0\t4\tname\t0\t-\n"
	out, _ := runIt(t, fasta, bed, Options{ForceStrand: true, PrintSeq: true})
	// Counts: T=4 (because RC was applied), A=0.
	if !strings.Contains(out, "\t0\t0\t0\t4\t0\t0\t4\tTTTT\n") {
		t.Errorf("expected RC seq column TTTT, got:\n%s", out)
	}
	// Without -s, same record should report A=4.
	out2, _ := runIt(t, fasta, bed, Options{PrintSeq: true})
	if !strings.Contains(out2, "\t4\t0\t0\t0\t0\t0\t4\tAAAA\n") {
		t.Errorf("expected non-RC seq column AAAA, got:\n%s", out2)
	}
}

// TestRun_PatternColumn checks -pattern and -C interaction. Pattern "AC"
// inside ACGTACGTAC occurs at positions 0, 4 (case-sensitive) = 2 hits.
func TestRun_PatternColumn(t *testing.T) {
	fasta := ">chr1\nACGTACGTAC\n"
	bed := "chr1\t0\t10\n"
	out, _ := runIt(t, fasta, bed, Options{Pattern: "AC", HasPattern: true})
	if !strings.HasSuffix(strings.TrimSpace(out), "\t3") {
		// AC at 0, 4, 8 = 3 hits in ACGTACGTAC.
		t.Errorf("expected 3 pattern hits at end of line, got:\n%s", out)
	}
	// -C lowercases both; result should be unchanged.
	out2, _ := runIt(t, fasta, bed, Options{Pattern: "ac", HasPattern: true, IgnoreCase: true})
	if !strings.HasSuffix(strings.TrimSpace(out2), "\t3") {
		t.Errorf("expected 3 hits with -C lowercase pattern, got:\n%s", out2)
	}
	// Lowercase pattern without -C should miss (Fetch uppercases the seq).
	out3, _ := runIt(t, fasta, bed, Options{Pattern: "ac", HasPattern: true})
	if !strings.HasSuffix(strings.TrimSpace(out3), "\t0") {
		t.Errorf("expected 0 hits with case-sensitive lowercase pattern, got:\n%s", out3)
	}
}

// TestRun_SkipsMissingAndOOB exercises warning paths.
func TestRun_SkipsMissingAndOOB(t *testing.T) {
	fasta := ">chr1\nACGT\n"
	// Out-of-bounds end (5 > 4), missing contig chr2, zero-length [3,3).
	bed := "chr1\t0\t5\nchr2\t0\t3\nchr1\t3\t3\nchr1\t0\t4\n"
	out, warn := runIt(t, fasta, bed, Options{})
	// Only the last record (chr1\t0\t4) should produce a data row.
	dataLines := 0
	for _, line := range strings.Split(strings.TrimSpace(out), "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		dataLines++
	}
	if dataLines != 1 {
		t.Errorf("want 1 data row after skips, got %d:\n%s", dataLines, out)
	}
	if !strings.Contains(warn, "chr2") || !strings.Contains(warn, "beyond") {
		t.Errorf("expected warnings for chr2 + OOB, got: %q", warn)
	}
}

// TestRun_FullHeader documents the upstream `-fullHeader` behaviour. The htslib
// shipped with bedtools builds the `.fai` on the first whitespace token even
// with `-fullHeader`, so a BED whose chrom is the *full* multi-token header
// (`chr1 extra info`) is not found and produces no data row — upstream skips it
// with a "size (0 bp)" warning. A first-token chrom (`chr1`) still resolves,
// exactly as in the default mode.
func TestRun_FullHeader(t *testing.T) {
	dir := t.TempDir()
	fasta := ">chr1 extra info\nACGT\n"
	faPath := writeTempFile(t, dir, "ref.fa", fasta)

	// Full multi-token header: not found, no data row, "size (0 bp)" skip.
	full := runFullHeader(t, faPath, "chr1 extra info\t0\t4\n")
	if dataRows(full.out) != 0 {
		t.Errorf("full-token chrom should yield no data row, got:\n%s", full.out)
	}
	if !strings.Contains(full.warn, "beyond the length") || !strings.Contains(full.warn, "0 bp") {
		t.Errorf("expected a size-0 skip warning, got: %q", full.warn)
	}

	// First token: resolves exactly as default mode (the .fai keys on `chr1`).
	first := runFullHeader(t, faPath, "chr1\t0\t4\n")
	if dataRows(first.out) != 1 {
		t.Errorf("first-token chrom should resolve to 1 data row, got:\n%s", first.out)
	}
	if !strings.Contains(first.out, "\t1\t1\t1\t1\t0\t0\t4") {
		t.Errorf("unexpected counts for chr1:\n%s", first.out)
	}
}

type fullHeaderResult struct{ out, warn string }

// runFullHeader runs bednuc with -fullHeader over an inline BED against faPath.
func runFullHeader(t *testing.T, faPath, bed string) fullHeaderResult {
	t.Helper()
	var out, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), faPath, &out, &warn, Options{FullHeader: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	return fullHeaderResult{out: out.String(), warn: warn.String()}
}

// dataRows counts non-header, non-blank output lines.
func dataRows(out string) int {
	n := 0
	for _, line := range strings.Split(strings.TrimSpace(out), "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		n++
	}
	return n
}

// TestFormatHeader_PrintSeqAndPattern verifies column order for combined opts.
func TestFormatHeader_PrintSeqAndPattern(t *testing.T) {
	// bedType=4 (chrom,start,end,name) + seq + pattern.
	h := FormatHeader(4, true, true)
	want := "#1_usercol\t2_usercol\t3_usercol\t4_usercol\t5_pct_at\t6_pct_gc\t7_num_A\t8_num_C\t9_num_G\t10_num_T\t11_num_N\t12_num_oth\t13_seq_len\t14_seq\t15_user_patt_count\n"
	if h != want {
		t.Errorf("header mismatch.\nwant: %q\n got: %q", want, h)
	}
}
