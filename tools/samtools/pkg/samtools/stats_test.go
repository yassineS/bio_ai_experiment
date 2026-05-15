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
	t.Run("sortInt64s", func(t *testing.T) {
		s := []int64{5, 1, 3, 2, 4}
		sortInt64s(s)
		want := []int64{1, 2, 3, 4, 5}
		for i, v := range s {
			if v != want[i] {
				t.Fatalf("sortInt64s[%d]: got %d, want %d", i, v, want[i])
			}
		}
	})
	t.Run("sortInt64s_short", func(t *testing.T) {
		// 1 element, 0 elements — should not panic.
		sortInt64s(nil)
		sortInt64s([]int64{42})
	})
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
