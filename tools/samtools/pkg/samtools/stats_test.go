package samtools

import (
	"bytes"
	"hash/crc32"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// statsFixture returns the absolute path to a fixture under
// tools/samtools/testdata/parity/stat/.
func statsFixture(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join("..", "..", "testdata", "parity", "stat", name)
}

// extractSN returns only the body lines starting with "SN\t" from a
// stats-text blob, joined by '\n' and trailing-newline normalised.
func extractSN(blob string) string {
	var b strings.Builder
	for _, line := range strings.Split(blob, "\n") {
		if strings.HasPrefix(line, "SN\t") {
			b.WriteString(line)
			b.WriteByte('\n')
		}
	}
	return b.String()
}

// TestStatsSNParity compares our SN section against upstream's
// stats expected outputs. The other sections (FFQ/RL/MAPQ/IS/etc.) have
// version-dependent header text and ordering quirks that make
// byte-for-byte parity fragile, so we restrict the comparison to the
// section the entire downstream ecosystem actually consumes.
func TestStatsSNParity(t *testing.T) {
	cases := []struct {
		name   string
		input  string
		expect string
		opts   StatsOptions
	}{
		{"1_map_cigar", "1_map_cigar.sam", "1.stats.expected", StatsOptions{}},
		{"2_equal_cigar", "2_equal_cigar_full_seq.sam", "2.stats.expected", StatsOptions{}},
		{"5_insert_cigar", "5_insert_cigar.sam", "5.stats.expected", StatsOptions{}},
		{"7_supp", "7_supp.sam", "7.stats.expected", StatsOptions{}},
		{"8_secondary", "8_secondary.sam", "8.stats.expected", StatsOptions{}},
		{"10_map_cigar_unsorted", "10_map_cigar.sam", "10.stats.expected", StatsOptions{}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			in, err := os.Open(statsFixture(t, tc.input))
			if err != nil {
				t.Fatalf("open input: %v", err)
			}
			defer in.Close()
			var out bytes.Buffer
			if err := Stats(in, &out, tc.opts); err != nil {
				t.Fatalf("Stats: %v", err)
			}
			gotSN := extractSN(out.String())
			expected, err := os.ReadFile(statsFixture(t, tc.expect))
			if err != nil {
				t.Fatalf("read expected: %v", err)
			}
			wantSN := extractSN(string(expected))
			if gotSN != wantSN {
				t.Fatalf("SN section differs\n--- want\n%s\n--- got\n%s", wantSN, gotSN)
			}
		})
	}
}

// extractSection returns only the body lines starting with "<tag>\t" from a
// stats-text blob, joined by '\n'.
func extractSection(blob, tag string) string {
	prefix := tag + "\t"
	var b strings.Builder
	for _, line := range strings.Split(blob, "\n") {
		if strings.HasPrefix(line, prefix) {
			b.WriteString(line)
			b.WriteByte('\n')
		}
	}
	return b.String()
}

// TestStatsCycleSectionParity compares our per-cycle and base-content sections
// (FFQ/LFQ/GCF/GCL/GCC/GCT/IC/ID) against upstream's stats expected outputs.
// Only .sam fixtures are used: .bam fixtures would additionally require BGZF
// byte-parity which is unrelated to these sections.
func TestStatsCycleSectionParity(t *testing.T) {
	cases := []struct {
		name   string
		input  string
		expect string
	}{
		{"1_map_cigar", "1_map_cigar.sam", "1.stats.expected"},
		{"2_equal_cigar", "2_equal_cigar_full_seq.sam", "2.stats.expected"},
		{"5_insert_cigar", "5_insert_cigar.sam", "5.stats.expected"},
		{"7_supp", "7_supp.sam", "7.stats.expected"},
		{"8_secondary", "8_secondary.sam", "8.stats.expected"},
		{"10_map_cigar_unsorted", "10_map_cigar.sam", "10.stats.expected"},
	}
	sections := []string{"FFQ", "LFQ", "GCF", "GCL", "GCC", "GCT", "IC", "ID"}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			in, err := os.Open(statsFixture(t, tc.input))
			if err != nil {
				t.Fatalf("open input: %v", err)
			}
			defer in.Close()
			var out bytes.Buffer
			if err := Stats(in, &out, StatsOptions{}); err != nil {
				t.Fatalf("Stats: %v", err)
			}
			expected, err := os.ReadFile(statsFixture(t, tc.expect))
			if err != nil {
				t.Fatalf("read expected: %v", err)
			}
			got, want := out.String(), string(expected)
			for _, sec := range sections {
				if extractSection(got, sec) != extractSection(want, sec) {
					t.Errorf("%s section differs\n--- want\n%s\n--- got\n%s",
						sec, extractSection(want, sec), extractSection(got, sec))
				}
			}
		})
	}
}

// extractCommentHeader returns the first comment line (starting with "# ")
// whose remainder begins with prefix, trailing-newline normalised.
func extractCommentHeader(blob, prefix string) string {
	for _, line := range strings.Split(blob, "\n") {
		if strings.HasPrefix(line, "# "+prefix) {
			return line
		}
	}
	return ""
}

// TestStatsCHKParity compares our leading CHK CRC32 checksum block (both the
// two comment header lines and the CHK data row) against upstream's stats
// expected outputs. CHK byte-parity exercises the BAM 4-bit nibble packing
// used for the sequences checksum.
func TestStatsCHKParity(t *testing.T) {
	cases := []struct {
		name   string
		input  string
		expect string
	}{
		{"1_map_cigar", "1_map_cigar.sam", "1.stats.expected"},
		{"2_equal_cigar", "2_equal_cigar_full_seq.sam", "2.stats.expected"},
		{"5_insert_cigar", "5_insert_cigar.sam", "5.stats.expected"},
		{"7_supp", "7_supp.sam", "7.stats.expected"},
		{"8_secondary", "8_secondary.sam", "8.stats.expected"},
		{"10_map_cigar_unsorted", "10_map_cigar.sam", "10.stats.expected"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			in, err := os.Open(statsFixture(t, tc.input))
			if err != nil {
				t.Fatalf("open input: %v", err)
			}
			defer in.Close()
			var out bytes.Buffer
			if err := Stats(in, &out, StatsOptions{}); err != nil {
				t.Fatalf("Stats: %v", err)
			}
			expected, err := os.ReadFile(statsFixture(t, tc.expect))
			if err != nil {
				t.Fatalf("read expected: %v", err)
			}
			got, want := out.String(), string(expected)
			if extractSection(got, "CHK") != extractSection(want, "CHK") {
				t.Errorf("CHK row differs\n--- want\n%s--- got\n%s",
					extractSection(want, "CHK"), extractSection(got, "CHK"))
			}
			for _, prefix := range []string{"CHK, Checksum", "CHK, CRC32"} {
				if extractCommentHeader(got, prefix) != extractCommentHeader(want, prefix) {
					t.Errorf("CHK header %q differs\nwant: %q\ngot:  %q", prefix,
						extractCommentHeader(want, prefix), extractCommentHeader(got, prefix))
				}
			}
		})
	}
}

// TestStatsCOVParity compares our COV coverage-distribution histogram against
// upstream's stats expected outputs. Coordinate-sorted fixtures must emit COV;
// the unsorted fixture must omit the COV section entirely, matching upstream's
// is_sorted gating.
func TestStatsCOVParity(t *testing.T) {
	cases := []struct {
		name   string
		input  string
		expect string
		sorted bool
	}{
		{"1_map_cigar", "1_map_cigar.sam", "1.stats.expected", true},
		{"2_equal_cigar", "2_equal_cigar_full_seq.sam", "2.stats.expected", true},
		{"5_insert_cigar", "5_insert_cigar.sam", "5.stats.expected", true},
		{"7_supp", "7_supp.sam", "7.stats.expected", true},
		{"8_secondary", "8_secondary.sam", "8.stats.expected", true},
		{"10_map_cigar_unsorted", "10_map_cigar.sam", "10.stats.expected", false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			in, err := os.Open(statsFixture(t, tc.input))
			if err != nil {
				t.Fatalf("open input: %v", err)
			}
			defer in.Close()
			var out bytes.Buffer
			if err := Stats(in, &out, StatsOptions{}); err != nil {
				t.Fatalf("Stats: %v", err)
			}
			expected, err := os.ReadFile(statsFixture(t, tc.expect))
			if err != nil {
				t.Fatalf("read expected: %v", err)
			}
			got, want := out.String(), string(expected)
			if extractSection(got, "COV") != extractSection(want, "COV") {
				t.Errorf("COV rows differ\n--- want\n%s--- got\n%s",
					extractSection(want, "COV"), extractSection(got, "COV"))
			}
			gotHdr := extractCommentHeader(got, "Coverage distribution")
			wantHdr := extractCommentHeader(want, "Coverage distribution")
			if gotHdr != wantHdr {
				t.Errorf("COV header differs\nwant: %q\ngot:  %q", wantHdr, gotHdr)
			}
			if !tc.sorted && gotHdr != "" {
				t.Errorf("unsorted input must omit COV section, got header %q", gotHdr)
			}
		})
	}
}

// TestStatsCOVStreamingFlush confirms COV depth is accumulated in a bounded
// per-contig sliding window: positions strictly below an advancing record's
// start are flushed into the cov bin array, and a contig change flushes the
// previous contig's window in full. The covWindow map must never retain the
// whole genome — its size stays O(longest read span), not O(positions seen).
func TestStatsCOVStreamingFlush(t *testing.T) {
	mapped := func(rname string, pos int32, cig string) *sam.Record {
		ops, err := sam.ParseCigar(cig)
		if err != nil {
			t.Fatalf("ParseCigar(%q): %v", cig, err)
		}
		return &sam.Record{RName: rname, Pos: pos, Cigar: ops}
	}

	c := newStatsCounters()
	c.initCovBins(StatsOptions{})

	// First record on contig "a": 10bp at pos 0 fills positions 0..9.
	c.accumulateCoverage(mapped("a", 0, "10M"))
	if len(c.covWindow) != 10 {
		t.Fatalf("after first record: window size = %d, want 10", len(c.covWindow))
	}

	// A second record far downstream on the same contig finalizes every
	// earlier position. The window must shrink to just the new span.
	c.accumulateCoverage(mapped("a", 1_000_000, "5M"))
	if len(c.covWindow) != 5 {
		t.Fatalf("after distant record: window size = %d, want 5 "+
			"(earlier positions must be flushed, not retained)", len(c.covWindow))
	}

	// Switching contigs must flush "a"'s remaining window entirely.
	c.accumulateCoverage(mapped("b", 0, "3M"))
	if c.covContig != "b" {
		t.Fatalf("covContig = %q, want \"b\"", c.covContig)
	}
	if len(c.covWindow) != 3 {
		t.Fatalf("after contig change: window size = %d, want 3", len(c.covWindow))
	}

	// End-of-input flush empties the window; all depth is now binned.
	c.flushCoverageWindow(1 << 30)
	if len(c.covWindow) != 0 {
		t.Fatalf("after final flush: window size = %d, want 0", len(c.covWindow))
	}
	// 10 + 5 + 3 = 18 single-depth positions binned at depth 1.
	var total int64
	for _, v := range c.cov {
		total += v
	}
	if total != 18 {
		t.Fatalf("binned positions = %d, want 18", total)
	}
}

// TestStatsCHKMissingQual confirms a record with QUAL "*" feeds seq_len bytes
// of 0xff to the qualities checksum, matching upstream's bam_get_qual
// missing-qual buffer, while a SEQ-less record contributes only its name.
func TestStatsCHKMissingQual(t *testing.T) {
	c := newStatsCounters()
	c.updateChecksum(&sam.Record{QName: "r1", Seq: "ACGT"}) // QUAL "*"
	missing := []byte{0xff, 0xff, 0xff, 0xff}
	wantNames := crc32.ChecksumIEEE([]byte("r1"))
	wantReads := crc32.ChecksumIEEE((&sam.Record{Seq: "ACGT"}).PackedSeq())
	wantQuals := crc32.ChecksumIEEE(missing)
	if c.ChkNames != wantNames || c.ChkReads != wantReads || c.ChkQuals != wantQuals {
		t.Fatalf("missing-qual checksum mismatch: got (%08x %08x %08x) want (%08x %08x %08x)",
			c.ChkNames, c.ChkReads, c.ChkQuals, wantNames, wantReads, wantQuals)
	}
	// A record with SEQ "*" contributes only to ChkNames.
	c2 := newStatsCounters()
	c2.updateChecksum(&sam.Record{QName: "r2"})
	if c2.ChkNames != crc32.ChecksumIEEE([]byte("r2")) || c2.ChkReads != 0 || c2.ChkQuals != 0 {
		t.Fatalf("SEQ-less record should only checksum the name: got (%08x %08x %08x)",
			c2.ChkNames, c2.ChkReads, c2.ChkQuals)
	}
}

// TestStatsCycleSectionsSparseSuppressed confirms --sparse omits the per-cycle
// and base-content section bodies.
func TestStatsCycleSectionsSparseSuppressed(t *testing.T) {
	in, err := os.Open(statsFixture(t, "1_map_cigar.sam"))
	if err != nil {
		t.Fatalf("open input: %v", err)
	}
	defer in.Close()
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{Sparse: true}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	for _, sec := range []string{"FFQ", "LFQ", "GCF", "GCL", "GCC", "GCT", "IC", "ID"} {
		if extractSection(out.String(), sec) != "" {
			t.Errorf("--sparse should suppress %s section", sec)
		}
	}
}

// TestStatsCountersBasic exercises the StatsCounters surface independently
// of the SAM reader.
func TestStatsCountersBasic(t *testing.T) {
	c := newStatsCounters()
	if c.RL == nil || c.ISInw == nil {
		t.Fatal("maps not initialised")
	}
	// Empty counters should still emit a valid SN block.
	var out bytes.Buffer
	if err := c.Write(&out, StatsOptions{MaxInsertSize: 8000}); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if !strings.Contains(out.String(), "SN\traw total sequences:\t0") {
		t.Fatalf("expected zero-RawTotal line; got:\n%s", out.String())
	}
	if !strings.Contains(out.String(), "# Read lengths.") {
		t.Fatalf("missing RL section header")
	}
}

// TestStatsHelpers covers small internal helpers.
func TestStatsHelpers(t *testing.T) {
	t.Run("sortedKeys", func(t *testing.T) {
		m := map[int64]int64{3: 1, 1: 1, 2: 1}
		ks := sortedKeys(m)
		want := []int64{1, 2, 3}
		for i, v := range ks {
			if v != want[i] {
				t.Fatalf("sortedKeys[%d]: got %d, want %d", i, v, want[i])
			}
		}
	})
}

// TestStatsFilter exercises the -F / -l filter options.
func TestStatsFilter(t *testing.T) {
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:1000\n"
	body := "r1\t99\tc\t1\t40\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\n" +
		"r2\t99\tc\t1\t40\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\n" +
		"r3\t1024\tc\t1\t40\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
	in := strings.NewReader(hdr + body)
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{FilteringFlag: 0x400}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if !strings.Contains(out.String(), "SN\tfiltered sequences:\t1") {
		t.Fatalf("expected 1 filtered sequence; got:\n%s", out.String())
	}
}

// TestStatsRemoveDups verifies the --remove-dups option subtracts the
// duplicate-flagged record from `sequences`.
func TestStatsRemoveDups(t *testing.T) {
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:1000\n"
	body := "r1\t99\tc\t1\t40\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\n" +
		"r2\t1123\tc\t1\t40\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\n"
	in := strings.NewReader(hdr + body)
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{RemoveDups: true}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if !strings.Contains(out.String(), "SN\tsequences:\t1") {
		t.Fatalf("expected sequences:1 after --remove-dups; got:\n%s", out.String())
	}
}

// TestStatsMinMAPQ verifies that the -q filter routes low-MAPQ records to
// the "filtered" bucket rather than counting them.
func TestStatsMinMAPQ(t *testing.T) {
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:1000\n"
	body := "r1\t99\tc\t1\t10\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\n" +
		"r2\t99\tc\t1\t40\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\n"
	in := strings.NewReader(hdr + body)
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{MinMAPQ: 30}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if !strings.Contains(out.String(), "SN\tfiltered sequences:\t1") {
		t.Fatalf("expected 1 filtered (MAPQ<30); got:\n%s", out.String())
	}
}

// TestStatsRLAndMAPQ exercises the RL and MAPQ sections to keep coverage
// on the histogram writers.
func TestStatsRLAndMAPQ(t *testing.T) {
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:1000\n"
	body := "r1\t99\tc\t1\t40\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\n" +
		"r2\t99\tc\t1\t40\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\n"
	in := strings.NewReader(hdr + body)
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if !strings.Contains(out.String(), "RL\t10\t") {
		t.Fatalf("missing RL row")
	}
	if !strings.Contains(out.String(), "MAPQ\t40\t2") {
		t.Fatalf("missing MAPQ row")
	}
}

// TestStatsQCFailCountedInTotals verifies that QC-failed primary records
// are folded into RawTotal/Sequences/1st-2nd-fragment/TotalLength and the
// ReadsQCFailed counter — matching upstream stats.c
// (collect_orig_read_stats runs for every primary record regardless of
// QCFAIL).
func TestStatsQCFailCountedInTotals(t *testing.T) {
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:1000\n"
	// r1 = ordinary primary; r2 = QC-failed primary (flag 0x200 = 512 set,
	// added to 0x63=99 to make 0x263=611).
	body := "r1\t99\tc\t1\t40\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\n" +
		"r2\t611\tc\t100\t40\t10M\t=\t1\t-100\tACGTACGTAC\tIIIIIIIIII\n"
	in := strings.NewReader(hdr + body)
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	got := out.String()
	for _, want := range []string{
		"SN\traw total sequences:\t2",
		"SN\tsequences:\t2",
		"SN\treads QC failed:\t1",
		"SN\ttotal length:\t20",
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("missing %q in:\n%s", want, got)
		}
	}
}

// TestStatsSortedDemotionOnOutOfOrder verifies that a header advertising
// SO:coordinate but a body that is NOT coordinate-sorted gets demoted to
// "is sorted: 0" — matching upstream stats.c:2327+1347 which initialises
// is_sorted=1 then demotes on any out-of-order record.
func TestStatsSortedDemotionOnOutOfOrder(t *testing.T) {
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:1000\n"
	// r1 at pos=200, r2 at pos=100 on same ref → body is out of order
	// even though @HD claims SO:coordinate.
	body := "r1\t0\tc\t200\t40\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
		"r2\t0\tc\t100\t40\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
	in := strings.NewReader(hdr + body)
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if !strings.Contains(out.String(), "SN\tis sorted:\t0") {
		t.Fatalf("expected `is sorted:\t0` after out-of-order body; got:\n%s", out.String())
	}
}

// TestStatsInsertSizeFromTLEN verifies the IS histogram and insert-size
// SN counters come from |TLEN|, not from a position-derived computation.
// Setup: a single mate-pair where TLEN diverges from
// (right.end - left.pos + 1) because the right mate has a leading soft clip.
func TestStatsInsertSizeFromTLEN(t *testing.T) {
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:1000\n"
	// Mate 1: pos=1, M-only, TLEN=+50.
	// Mate 2: pos=41, "5S5M" (5bp soft clip + 5bp match), TLEN=-50.
	// Position-derived (right.end - left.pos + 1) = (41+5-1) - 1 + 1 = 45.
	// TLEN-derived = 50. We assert 50.
	body := "r1\t99\tc\t1\t40\t10M\t=\t41\t50\tACGTACGTAC\tIIIIIIIIII\n" +
		"r1\t147\tc\t41\t40\t5S5M\t=\t1\t-50\tACGTACGTAC\tIIIIIIIIII\n"
	in := strings.NewReader(hdr + body)
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	got := out.String()
	if !strings.Contains(got, "SN\tinsert size average:\t50.0") {
		t.Fatalf("expected insert size average 50 (from TLEN); got:\n%s", got)
	}
}

// TestStatsProperPairDenominator verifies the proper-pair percentage uses
// Sequences (= 1st+2nd+other) as the denominator, matching upstream
// stats.c:1606. Not ReadsPaired.
func TestStatsProperPairDenominator(t *testing.T) {
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:1000\n"
	// 2 proper pairs (4 records), 1 unpaired primary. Sequences=5 (1st=3,
	// 2nd=2), ReadsPaired=4, ReadsProperlyPaired=4. Upstream %: 4/5*100=80.
	// Old buggy %: 4/4*100=100.
	body := "r1\t99\tc\t1\t40\t10M\t=\t41\t50\tACGTACGTAC\tIIIIIIIIII\n" +
		"r1\t147\tc\t41\t40\t10M\t=\t1\t-50\tACGTACGTAC\tIIIIIIIIII\n" +
		"r2\t99\tc\t2\t40\t10M\t=\t42\t50\tACGTACGTAC\tIIIIIIIIII\n" +
		"r2\t147\tc\t42\t40\t10M\t=\t2\t-50\tACGTACGTAC\tIIIIIIIIII\n" +
		"r3\t0\tc\t100\t40\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
	in := strings.NewReader(hdr + body)
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if !strings.Contains(out.String(), "SN\tpercentage of properly paired reads (%):\t80.0") {
		t.Fatalf("expected proper-pair %% = 80.0 (4/5); got:\n%s", out.String())
	}
}

// TestStatsSparseOmitsHistograms confirms the -x/--sparse flag skips the
// RL/MAPQ/IS section blocks but keeps SN.
func TestStatsSparseOmitsHistograms(t *testing.T) {
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:1000\n"
	body := "r1\t99\tc\t1\t40\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\n"
	in := strings.NewReader(hdr + body)
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{Sparse: true}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if !strings.Contains(out.String(), "SN\t") {
		t.Fatalf("SN must still be emitted under --sparse")
	}
	if strings.Contains(out.String(), "# Read lengths.") {
		t.Fatalf("--sparse should suppress RL header")
	}
}

// TestStatsTargetRegionsParity validates --target-regions against upstream's
// stat/11 fixtures. The SN and COV sections are compared byte-for-byte to the
// golden outputs for both the default coverage threshold and -g 4. The full
// report is not compared because upstream emits version-dependent sections
// (FBC/FTC/LBC/LTC/GCD) and a dense IS table that this v1 does not implement.
func TestStatsTargetRegionsParity(t *testing.T) {
	cases := []struct {
		name   string
		expect string
		opts   StatsOptions
	}{
		{
			"default-threshold",
			"11.stats.expected",
			StatsOptions{TargetBED: statsFixture(t, "11.stats.targets")},
		},
		{
			"cov-threshold-4",
			"11.stats.g4.expected",
			StatsOptions{TargetBED: statsFixture(t, "11.stats.targets"), CovThreshold: 4},
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			in, err := os.Open(statsFixture(t, "11_target.sam"))
			if err != nil {
				t.Fatalf("open input: %v", err)
			}
			defer in.Close()
			var out bytes.Buffer
			if err := Stats(in, &out, tc.opts); err != nil {
				t.Fatalf("Stats: %v", err)
			}
			expected, err := os.ReadFile(statsFixture(t, tc.expect))
			if err != nil {
				t.Fatalf("read expected: %v", err)
			}
			got, want := out.String(), string(expected)
			if extractSN(got) != extractSN(want) {
				t.Errorf("SN section differs\n--- want\n%s\n--- got\n%s",
					extractSN(want), extractSN(got))
			}
			if extractSection(got, "COV") != extractSection(want, "COV") {
				t.Errorf("COV rows differ\n--- want\n%s--- got\n%s",
					extractSection(want, "COV"), extractSection(got, "COV"))
			}
		})
	}
}

// TestStatsTargetRegionsSNLines confirms the two target-specific SN lines are
// emitted only when a target file is given and that "bases inside the target"
// equals the merged-interval base total.
func TestStatsTargetRegionsSNLines(t *testing.T) {
	in, err := os.Open(statsFixture(t, "11_target.sam"))
	if err != nil {
		t.Fatalf("open input: %v", err)
	}
	defer in.Close()
	var out bytes.Buffer
	opts := StatsOptions{TargetBED: statsFixture(t, "11.stats.targets")}
	if err := Stats(in, &out, opts); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	s := out.String()
	if !strings.Contains(s, "SN\tbases inside the target:\t42") {
		t.Errorf("expected 42 target bases (intervals [10,24] + merged [30,56]); got:\n%s", s)
	}
	if !strings.Contains(s, "SN\tpercentage of target genome with coverage > 0 (%):\t100.00") {
		t.Errorf("expected 100.00%% target coverage; got:\n%s", s)
	}

	// Without a target file the two SN lines must be absent.
	in2, err := os.Open(statsFixture(t, "1_map_cigar.sam"))
	if err != nil {
		t.Fatalf("open input: %v", err)
	}
	defer in2.Close()
	var out2 bytes.Buffer
	if err := Stats(in2, &out2, StatsOptions{}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if strings.Contains(out2.String(), "bases inside the target") {
		t.Errorf("target SN lines must not appear without --target-regions")
	}
}

// TestStatsTargetRegionsFilter confirms reads on a reference absent from the
// target file are excluded entirely: stat/11 has 28 records (2 on contig
// "alpha"), and restricting to ref1 targets must drop those to 26 sequences.
func TestStatsTargetRegionsFilter(t *testing.T) {
	in, err := os.Open(statsFixture(t, "11_target.sam"))
	if err != nil {
		t.Fatalf("open input: %v", err)
	}
	defer in.Close()
	var out bytes.Buffer
	opts := StatsOptions{TargetBED: statsFixture(t, "11.stats.targets")}
	if err := Stats(in, &out, opts); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if !strings.Contains(out.String(), "SN\traw total sequences:\t26") {
		t.Errorf("off-target contig reads must be filtered (want 26 sequences); got:\n%s", out.String())
	}
}

// TestBwaTrimRead exercises the BWA-style trimming algorithm port, including
// the BWA_MIN_RDLEN early return and bwa's documented off-by-one.
func TestBwaTrimRead(t *testing.T) {
	mkquals := func(n int, q byte) []byte {
		b := make([]byte, n)
		for i := range b {
			b[i] = q
		}
		return b
	}
	cases := []struct {
		name     string
		trimQual int
		quals    []byte
		reverse  bool
		want     int
	}{
		{"too-short", 20, mkquals(34, 5), false, 0},
		{"high-quality-no-trim", 20, mkquals(40, 30), false, 0},
		{"low-quality-trims-to-min", 20, mkquals(40, 10), false, 5},
		{"low-quality-reverse", 20, mkquals(40, 10), true, 5},
		{"zero-threshold-no-trim", 0, mkquals(40, 10), false, 0},
		{"exact-min-length", 20, mkquals(35, 5), false, 0},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := bwaTrimRead(tc.trimQual, tc.quals, len(tc.quals), tc.reverse)
			if got != tc.want {
				t.Errorf("bwaTrimRead = %d, want %d", got, tc.want)
			}
		})
	}
}

// TestStatsTrimQualityCounter confirms -q/--trim-quality feeds the "bases
// trimmed" SN counter: a long low-quality read is trimmed while a short or
// high-quality read is not, and the default (TrimQuality 0) trims nothing.
func TestStatsTrimQualityCounter(t *testing.T) {
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:100000\n"
	// One 40bp read with uniformly low quality ('+' = Phred 10): with a
	// trim threshold of 20 bwa trims it down to BWA_MIN_RDLEN-1, i.e. 5
	// bases. A second high-quality read ('I' = Phred 40) trims nothing.
	lowQual := strings.Repeat("A", 40)
	lowQ := strings.Repeat("+", 40)
	hiQ := strings.Repeat("I", 40)
	body := "r1\t0\tc\t1\t40\t40M\t*\t0\t0\t" + lowQual + "\t" + lowQ + "\n" +
		"r2\t0\tc\t100\t40\t40M\t*\t0\t0\t" + lowQual + "\t" + hiQ + "\n"

	var trimmed bytes.Buffer
	if err := Stats(strings.NewReader(hdr+body), &trimmed, StatsOptions{TrimQuality: 20}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if !strings.Contains(trimmed.String(), "SN\tbases trimmed:\t5") {
		t.Errorf("expected 5 trimmed bases with -q 20; got:\n%s", extractSN(trimmed.String()))
	}

	var untrimmed bytes.Buffer
	if err := Stats(strings.NewReader(hdr+body), &untrimmed, StatsOptions{}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if !strings.Contains(untrimmed.String(), "SN\tbases trimmed:\t0") {
		t.Errorf("default (no -q) must trim nothing; got:\n%s", extractSN(untrimmed.String()))
	}
}

// TestStatsTargetRegionsUnsorted confirms --target-regions rejects an
// unsorted input, mirroring upstream stats.c's is_in_regions hard error.
func TestStatsTargetRegionsUnsorted(t *testing.T) {
	dir := t.TempDir()
	target := filepath.Join(dir, "t.regions")
	if err := os.WriteFile(target, []byte("c\t1 1000\n"), 0o644); err != nil {
		t.Fatalf("write target file: %v", err)
	}
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:100000\n"
	// Second record precedes the first on the same reference.
	body := "r1\t0\tc\t500\t40\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
		"r2\t0\tc\t100\t40\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
	var out bytes.Buffer
	err := Stats(strings.NewReader(hdr+body), &out, StatsOptions{TargetBED: target})
	if err == nil {
		t.Fatal("expected an error for unsorted input with --target-regions")
	}
	if !strings.Contains(err.Error(), "sorted") {
		t.Errorf("error should mention sorting; got: %v", err)
	}
}

// TestStatsTargetRegionsDeletionBoundary is a regression test for the on-target
// "bases mapped (cigar)" counter (addBasesMappedCigar) when a CIGAR D op
// precedes an M op that straddles a target-interval boundary. Upstream
// stats.c:1316 handles BAM_CDEL with only `readlen += ncig` and NEVER advances
// the reference cursor iref on a deletion. An earlier draft of this port did
// `iref += ncig` on the D op, mis-positioning every op after a deletion.
//
// Synthetic input: one mapped read, POS=8 (1-based), CIGAR 5M3D10M. The target
// file gives a single ref interval [10,20], so regFrom=10, regTo=20.
//
// Hand calculation following upstream (iref NOT advanced on D):
//
//	iref=8.
//	op 5M:  ncig=5; iref(8)<regFrom(10) -> ncig -= 10-8 => 3; add 3; iref += 5 => 13.
//	op 3D:  iref unchanged => 13.
//	op 10M: ncig=10; iref(13) not <10; iref+ncig-1 = 22 > regTo(20)
//	        -> ncig -= 22-20 => 8; add 8; iref += 10 => 23.
//	total bases mapped (cigar) = 3 + 8 = 11.
//
// With the old buggy `iref += ncig` on D, iref would be 16 at the 10M op:
// iref+ncig-1 = 25 > 20 -> ncig -= 5 => 5, giving total 3 + 5 = 8. The assertion
// on 11 therefore fails against the old behaviour.
func TestStatsTargetRegionsDeletionBoundary(t *testing.T) {
	dir := t.TempDir()
	target := filepath.Join(dir, "t.regions")
	if err := os.WriteFile(target, []byte("c\t10 20\n"), 0o644); err != nil {
		t.Fatalf("write target file: %v", err)
	}
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:100000\n"
	// 15 query bases (5M + 10M; the 3D consumes no query).
	body := "r1\t0\tc\t8\t40\t5M3D10M\t*\t0\t0\tACGTACGTACGTACG\tIIIIIIIIIIIIIII\n"
	var out bytes.Buffer
	if err := Stats(strings.NewReader(hdr+body), &out, StatsOptions{TargetBED: target}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if !strings.Contains(out.String(), "SN\tbases mapped (cigar):\t11\t") {
		t.Errorf("expected on-target bases mapped (cigar) = 11 (D must not advance iref); got:\n%s", out.String())
	}
}
