package bedcoverage

import (
	"bytes"
	"fmt"
	"strings"
	"testing"
)

// samHeader is a minimal SAM header used by the alignment-input unit tests.
const samHeader = "@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"

// samRead builds a single mapped SAM alignment line. pos is 1-based.
func samRead(qname string, flag, pos int, cigar string) string {
	return fmt.Sprintf("%s\t%d\tchr1\t%d\t60\t%s\t*\t0\t0\t*\t*\n", qname, flag, pos, cigar)
}

// TestCoverage_SAMDatabase feeds a SAM text stream as the -b side and checks
// that coverage routes through the alnbed auto-detect path. The spliced read
// (10M10N10M) at chr1:1 yields a BED12 record with two blocks ([0,10),
// [20,30)); without -split the whole [0,30) span counts, with -split only the
// 20 bp of blocks count.
func TestCoverage_SAMDatabase(t *testing.T) {
	a := "chr1\t0\t30\n"
	b := samHeader + samRead("r1", 0, 1, "10M10N10M")

	// No -split: whole span [0,30) covers A entirely.
	got := runCoverage(t, a, b, Options{})
	want := "chr1\t0\t30\t1\t30\t30\t1.0000000\n"
	if got != want {
		t.Errorf("SAM -b no split:\nwant: %q\ngot:  %q", want, got)
	}

	// -split: only the two 10 bp blocks count -> 20 bp, fraction 20/30.
	got = runCoverage(t, a, b, Options{Split: true})
	want = "chr1\t0\t30\t2\t20\t30\t0.6666667\n"
	if got != want {
		t.Errorf("SAM -b split:\nwant: %q\ngot:  %q", want, got)
	}
}

// TestCoverage_SAMUnmappedSkipped verifies an unmapped read in the -b stream
// contributes nothing, matching upstream's BAM-to-BED conversion which only
// emits mapped alignments.
func TestCoverage_SAMUnmappedSkipped(t *testing.T) {
	a := "chr1\t0\t30\n"
	// flag 4 = unmapped.
	b := samHeader + samRead("r1", 4, 1, "30M")
	got := runCoverage(t, a, b, Options{})
	want := "chr1\t0\t30\t0\t0\t30\t0.0000000\n"
	if got != want {
		t.Errorf("unmapped skipped:\nwant: %q\ngot:  %q", want, got)
	}
}

// TestCoverage_BAMQuery exercises BAM/SAM on the -a (query) side: the
// alignment is converted to a BED12 record and coverage is computed against
// its whole span. The echoed A columns are the BED12 form of the read.
func TestCoverage_BAMQuery(t *testing.T) {
	// Query: a single 30 bp read at chr1:1 (one CIGAR block, [0,30)).
	a := samHeader + samRead("q1", 0, 1, "30M")
	b := "chr1\t0\t10\nchr1\t20\t30\n"
	got := runCoverage(t, a, b, Options{})
	// Two B features overlap, covering [0,10)+[20,30) = 20 bp of the 30 bp
	// read; count 2, covered 20, len 30, fraction 0.6666667. The A columns
	// are the BED12 echo of the read (block list without a trailing comma).
	want := "chr1\t0\t30\tq1\t60\t+\t0\t30\t0,0,0\t1\t30,\t0,\t2\t20\t30\t0.6666667\n"
	if got != want {
		t.Errorf("BAM -a query:\nwant: %q\ngot:  %q", want, got)
	}
}

// TestCoverage_SAMQuerySplitRejected confirms a spliced (blocked) alignment on
// the -a side under -split is rejected with a clear error rather than silently
// producing a wrong answer (upstream splits the query into blocks; we do not
// yet support that path).
func TestCoverage_SAMQuerySplitRejected(t *testing.T) {
	a := samHeader + samRead("q1", 0, 1, "10M10N10M")
	b := "chr1\t0\t30\n"
	var buf bytes.Buffer
	if _, err := Coverage(strings.NewReader(a), strings.NewReader(b), &buf, Options{Split: true}); err == nil {
		t.Fatal("expected error for spliced BAM/SAM query under -split")
	}
}

// TestSourceReader_AutoDetect checks the sniffing helper routes a leading-'@'
// SAM header to the alignment reader and plain BED text to the BED reader.
func TestSourceReader_AutoDetect(t *testing.T) {
	samSrc, err := sourceReader(strings.NewReader(samHeader + samRead("r1", 0, 1, "10M")))
	if err != nil {
		t.Fatalf("sourceReader(SAM): %v", err)
	}
	rec, err := samSrc.Read()
	if err != nil {
		t.Fatalf("read SAM record: %v", err)
	}
	if rec.Chrom != "chr1" || rec.ChromStart != 0 || rec.ChromEnd != 10 {
		t.Errorf("SAM record = %s:%d-%d, want chr1:0-10", rec.Chrom, rec.ChromStart, rec.ChromEnd)
	}

	bedSrc, err := sourceReader(strings.NewReader("chr1\t5\t15\n"))
	if err != nil {
		t.Fatalf("sourceReader(BED): %v", err)
	}
	rec, err = bedSrc.Read()
	if err != nil {
		t.Fatalf("read BED record: %v", err)
	}
	if rec.Chrom != "chr1" || rec.ChromStart != 5 || rec.ChromEnd != 15 {
		t.Errorf("BED record = %s:%d-%d, want chr1:5-15", rec.Chrom, rec.ChromStart, rec.ChromEnd)
	}
}
