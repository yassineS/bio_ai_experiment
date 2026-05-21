package samtools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
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
