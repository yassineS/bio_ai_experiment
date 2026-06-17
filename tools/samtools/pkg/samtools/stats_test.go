package samtools

import (
	"bufio"
	"bytes"
	"fmt"
	"hash/crc32"
	"os"
	"os/exec"
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

// upstreamStats runs the live `samtools stats` binary on the given input
// fixture (plus any extra args, e.g. `-r ref.fa`) and returns its full
// text output. The upstream binary is built on demand; a build failure is
// fatal. This replaces reading committed `.stats.expected` golden files.
func upstreamStats(t *testing.T, inputPath string, extraArgs ...string) string {
	t.Helper()
	bin := upstreamSamtools(t)
	args := append([]string{"stats"}, extraArgs...)
	args = append(args, inputPath)
	cmd := exec.Command(bin, args...)
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream samtools stats %v: %v\n%s", args, err, errBuf.String())
	}
	return out.String()
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
		name  string
		input string
		opts  StatsOptions
	}{
		{"1_map_cigar", "1_map_cigar.sam", StatsOptions{}},
		{"2_equal_cigar", "2_equal_cigar_full_seq.sam", StatsOptions{}},
		{"5_insert_cigar", "5_insert_cigar.sam", StatsOptions{}},
		{"7_supp", "7_supp.sam", StatsOptions{}},
		{"8_secondary", "8_secondary.sam", StatsOptions{}},
		{"10_map_cigar_unsorted", "10_map_cigar.sam", StatsOptions{}},
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
			wantSN := extractSN(upstreamStats(t, statsFixture(t, tc.input)))
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
// (FFQ/LFQ/GCF/GCL/GCC/GCT/IC/ID) plus the IS insert-size section against
// upstream's stats expected outputs. Only .sam fixtures are used: .bam fixtures
// would additionally require BGZF byte-parity which is unrelated to these
// sections. The IS section exercises the full 0..ibulk-1 row range, including
// all-zero rows, in the default (non-sparse) mode.
func TestStatsCycleSectionParity(t *testing.T) {
	cases := []struct {
		name  string
		input string
	}{
		{"1_map_cigar", "1_map_cigar.sam"},
		{"2_equal_cigar", "2_equal_cigar_full_seq.sam"},
		{"5_insert_cigar", "5_insert_cigar.sam"},
		{"7_supp", "7_supp.sam"},
		{"8_secondary", "8_secondary.sam"},
		{"10_map_cigar_unsorted", "10_map_cigar.sam"},
	}
	sections := []string{"FFQ", "LFQ", "GCF", "GCL", "GCC", "GCT", "FBC", "FTC", "LBC", "LTC", "IC", "ID", "IS"}
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
			got, want := out.String(), upstreamStats(t, statsFixture(t, tc.input))
			for _, sec := range sections {
				if extractSection(got, sec) != extractSection(want, sec) {
					t.Errorf("%s section differs\n--- want\n%s\n--- got\n%s",
						sec, extractSection(want, sec), extractSection(got, sec))
				}
			}
		})
	}
}

// extractPrefixLines returns every line of a stats-text blob that starts with
// prefix, joined by '\n'. Unlike extractSection it does not append a tab, so
// it matches the barcode rows whose tag carries a trailing segment digit
// (e.g. "BCC1", "QTQ2").
func extractPrefixLines(blob, prefix string) string {
	var b strings.Builder
	for _, line := range strings.Split(blob, "\n") {
		if strings.HasPrefix(line, prefix) {
			b.WriteString(line)
			b.WriteByte('\n')
		}
	}
	return b.String()
}

// TestStatsBarcodeParity compares every section of our stats output against
// upstream's barcode golden files byte-for-byte. The two _ok fixtures exercise
// the per-barcode ACGT-content (<tag>C) and quality (<tag>Q) sections: the BC
// fixture has only BC/QT tags (BCC/QTQ), the OX fixture additionally carries
// OX/BZ tags (OXC/BZQ). The exact upstream invocations are test.pl:3325-3326
// `samtools stats <13_barcodes_ok[ _ox_bz].sam>`.
func TestStatsBarcodeParity(t *testing.T) {
	cases := []struct {
		name     string
		input    string
		sections []string
	}{
		{
			"bc", "13_barcodes_ok.sam",
			[]string{"CHK", "SN", "FFQ", "LFQ", "GCF", "GCL", "GCC", "GCT", "FBC", "FTC", "LBC", "LTC", "BCC", "QTQ", "IS", "RL", "FRL", "LRL", "MAPQ", "ID", "IC", "COV", "GCD"},
		},
		{
			"ox_bz", "13_barcodes_ok_ox_bz.sam",
			[]string{"CHK", "SN", "FFQ", "LFQ", "GCF", "GCL", "GCC", "GCT", "FBC", "FTC", "LBC", "LTC", "BCC", "QTQ", "OXC", "BZQ", "IS", "RL", "FRL", "LRL", "MAPQ", "ID", "IC", "COV", "GCD"},
		},
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
			got, want := out.String(), upstreamStats(t, statsFixture(t, tc.input))
			// The per-barcode sections carry a trailing segment digit on the
			// tag (BCC1/BCC2/QTQ1/...), so they need prefix matching.
			isBarcodeSec := map[string]bool{"BCC": true, "QTQ": true, "OXC": true, "BZQ": true}
			for _, sec := range tc.sections {
				var g, w string
				if isBarcodeSec[sec] {
					g, w = extractPrefixLines(got, sec), extractPrefixLines(want, sec)
				} else {
					g, w = extractSection(got, sec), extractSection(want, sec)
				}
				if g != w {
					t.Errorf("%s section differs\n--- want\n%s\n--- got\n%s", sec, w, g)
				}
			}
		})
	}
}

// TestStatsBarcodeFailInputs verifies the malformed-barcode fixtures (the
// inconsistent-length and misplaced-separator inputs from test.pl:3327-3329)
// run to completion: upstream warns on stderr and skips the offending record
// rather than aborting, so the run must succeed and still emit barcode
// sections. The exact byte output is intentionally NOT compared — test.pl
// marks these cases expect_fail, i.e. their output deliberately diverges from
// the clean baseline.
func TestStatsBarcodeFailInputs(t *testing.T) {
	for _, name := range []string{
		"13_barcodes_fail_bc_length.sam",
		"13_barcodes_fail_hyphen.sam",
		"13_barcodes_fail_qt_length.sam",
	} {
		t.Run(name, func(t *testing.T) {
			in, err := os.Open(statsFixture(t, name))
			if err != nil {
				t.Fatalf("open input: %v", err)
			}
			defer in.Close()
			var out bytes.Buffer
			if err := Stats(in, &out, StatsOptions{}); err != nil {
				t.Fatalf("Stats must not abort on a malformed barcode: %v", err)
			}
			if extractPrefixLines(out.String(), "BCC") == "" {
				t.Errorf("expected a BCC barcode section despite malformed records")
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
		name  string
		input string
	}{
		{"1_map_cigar", "1_map_cigar.sam"},
		{"2_equal_cigar", "2_equal_cigar_full_seq.sam"},
		{"5_insert_cigar", "5_insert_cigar.sam"},
		{"7_supp", "7_supp.sam"},
		{"8_secondary", "8_secondary.sam"},
		{"10_map_cigar_unsorted", "10_map_cigar.sam"},
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
			got, want := out.String(), upstreamStats(t, statsFixture(t, tc.input))
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
		sorted bool
	}{
		{"1_map_cigar", "1_map_cigar.sam", true},
		{"2_equal_cigar", "2_equal_cigar_full_seq.sam", true},
		{"5_insert_cigar", "5_insert_cigar.sam", true},
		{"7_supp", "7_supp.sam", true},
		{"8_secondary", "8_secondary.sam", true},
		{"10_map_cigar_unsorted", "10_map_cigar.sam", false},
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
			got, want := out.String(), upstreamStats(t, statsFixture(t, tc.input))
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

// TestStatsGCDParity compares our GCD GC-depth distribution against upstream's
// stats expected outputs in the no-reference path (GC approximated from the
// read sequences). Coordinate-sorted fixtures must emit one GCD row plus the
// section comment; the unsorted fixture must omit the section entirely,
// matching upstream's is_sorted gating. Every upstream fixture keeps its reads
// inside a single 20 kbp bin, so the finalised output is the empty placeholder
// row "GCD 0.0 100.000 0.000 0.000 0.000 0.000 0.000".
func TestStatsGCDParity(t *testing.T) {
	cases := []struct {
		name   string
		input  string
		sorted bool
	}{
		{"1_map_cigar", "1_map_cigar.sam", true},
		{"2_equal_cigar", "2_equal_cigar_full_seq.sam", true},
		{"5_insert_cigar", "5_insert_cigar.sam", true},
		{"7_supp", "7_supp.sam", true},
		{"8_secondary", "8_secondary.sam", true},
		{"10_map_cigar_unsorted", "10_map_cigar.sam", false},
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
			got, want := out.String(), upstreamStats(t, statsFixture(t, tc.input))
			if extractSection(got, "GCD") != extractSection(want, "GCD") {
				t.Errorf("GCD rows differ\n--- want\n%s--- got\n%s",
					extractSection(want, "GCD"), extractSection(got, "GCD"))
			}
			gotHdr := extractCommentHeader(got, "GC-depth")
			wantHdr := extractCommentHeader(want, "GC-depth")
			if gotHdr != wantHdr {
				t.Errorf("GCD header differs\nwant: %q\ngot:  %q", wantHdr, gotHdr)
			}
			if !tc.sorted && gotHdr != "" {
				t.Errorf("unsorted input must omit GCD section, got header %q", gotHdr)
			}
		})
	}
}

// TestStatsGCDFloat32Precision guards the float32 width of the GCD percentile
// scaling against a regression back to float64. Upstream's gcd_percentile is a
// `float` and the caller scales it by a `float avg_read_length` and the integer
// bin size in float, widening to double only at the printf. Computing this in
// float64 instead lands ~1 ULP off and flips the %.3f rounding on boundary
// values (the same class of bug as the errmod float-vs-double port). The
// depths/total/reads below were found by search to be exactly such a boundary:
// the float32 chain prints 0.122 (matching upstream) while a float64 chain
// would print 0.121. The multi-bin pipeline fixture surfaced the original miss;
// the in-tree GCD fixtures all fit one bin so they could not.
func TestStatsGCDFloat32Precision(t *testing.T) {
	grp := []gcDepth{{depth: 299}, {depth: 375}, {depth: 470}, {depth: 673},
		{depth: 854}, {depth: 859}, {depth: 916}}
	nbins := len(grp)
	// avg_read_length = (float)total_len / nreads, in float32 (stats.c:1586).
	avg := float32(482611) / float32(59383)
	binSize := float32(20000)
	got := fmt.Sprintf("%.3f", float64(gcdPercentile(grp, nbins, 10)*avg/binSize))
	if got != "0.122" {
		t.Errorf("GCD 10th-percentile scaled value = %s, want 0.122 "+
			"(float32 chain); a float64 chain would print 0.121 — the precision "+
			"width must match upstream", got)
	}
}

// TestStatsGCDReferenceParity exercises the --ref-seq GC-depth path: the
// upstream stat/1-8 fixtures were generated with `-r test.fa`, and because
// every read stays inside a single 20 kbp bin only the empty placeholder bin
// is finalised and printed — so the reference and no-reference paths emit
// byte-identical GCD output. This test feeds the reference FASTA and confirms
// the reference-derived GC path still matches the golden files exactly.
func TestStatsGCDReferenceParity(t *testing.T) {
	cases := []struct {
		input string
	}{
		{"1_map_cigar.sam"},
		{"2_equal_cigar_full_seq.sam"},
		{"5_insert_cigar.sam"},
		{"7_supp.sam"},
		{"8_secondary.sam"},
	}
	for _, tc := range cases {
		t.Run(tc.input, func(t *testing.T) {
			in, err := os.Open(statsFixture(t, tc.input))
			if err != nil {
				t.Fatalf("open input: %v", err)
			}
			defer in.Close()
			var out bytes.Buffer
			opts := StatsOptions{RefSeq: statsFixture(t, "test.fa")}
			if err := Stats(in, &out, opts); err != nil {
				t.Fatalf("Stats: %v", err)
			}
			got := out.String()
			want := upstreamStats(t, statsFixture(t, tc.input), "-r", statsFixture(t, "test.fa"))
			if extractSection(got, "GCD") != extractSection(want, "GCD") {
				t.Errorf("GCD rows differ (reference path)\n--- want\n%s--- got\n%s",
					extractSection(want, "GCD"), extractSection(got, "GCD"))
			}
		})
	}
}

// TestStatsGCDMultiBin exercises the multi-segment GC-depth algorithm that the
// single-bin upstream fixtures never reach. Synthetic reads are spread across
// three 20 kbp bins on one contig so several segments are finalised, sorted by
// GC and grouped — verifying the gcdPercentile interpolation, the
// unique-sequence-percentile formula and the upstream off-by-one whereby the
// gcd[0] placeholder and the last (gcd[gcdIdx]) segment are never finalised.
//
// With gcdIdx==3 the output loop finalises segments 0,1,2 but qsorts all four
// entries (0..3). Segment 3's raw, un-finalised GC value still participates in
// the sort and grouping — exactly as upstream stats.c does — so the printed
// rows derive from segments 0 (placeholder, GC 0), the un-finalised segment 3
// (raw GC 6.3 = 7*90/100), and finalised segment 1 (GC 10).
func TestStatsGCDMultiBin(t *testing.T) {
	c := newStatsCounters()
	c.gcdBinSize = 20000
	c.IsSorted = 1
	c.Sequences = 1
	c.TotalLength = 20000 // makes avg read length / bin size == 1.0

	mk := func(pos int32, gc, depth int) {
		for i := 0; i < depth; i++ {
			rec := &sam.Record{RName: "c", Pos: pos, Seq: strings.Repeat("N", 100)}
			c.accumulateGCD(rec, gc)
		}
	}
	// Bin A: pos 1, GC count 10 -> fraction 0.10, depth 3.
	// Bin B: pos 25001, GC count 50 -> fraction 0.50, depth 5.
	// Bin C: pos 50001, GC count 90 -> fraction 0.90, depth 7.
	mk(1, 10, 3)
	mk(25001, 50, 5)
	mk(50001, 90, 7)

	// Four segments exist: gcd[0] placeholder + three real bins, gcdIdx==3.
	if c.gcdIdx != 3 {
		t.Fatalf("expected gcdIdx 3, got %d", c.gcdIdx)
	}

	var out bytes.Buffer
	bw := bufio.NewWriter(&out)
	c.writeGCD(bw)
	bw.Flush()
	var gcd []string
	for _, l := range strings.Split(strings.TrimRight(out.String(), "\n"), "\n") {
		if strings.HasPrefix(l, "GCD\t") {
			gcd = append(gcd, l)
		}
	}
	// Three grouped rows: GC 0.0 (placeholder), 6.3 (un-finalised segment),
	// 10.0 (finalised bin A). Bin B (GC 50, segment 2) is finalised but sorts
	// last and is not reached by the igcd-bounded group loop.
	want := []string{
		"GCD\t0.0\t50.000\t0.000\t0.000\t0.000\t0.000\t0.000",
		"GCD\t6.3\t75.000\t7.000\t7.000\t7.000\t7.000\t7.000",
		"GCD\t10.0\t100.000\t3.000\t3.000\t3.000\t3.000\t3.000",
	}
	if len(gcd) != len(want) {
		t.Fatalf("expected %d GCD rows, got %d:\n%s", len(want), len(gcd), out.String())
	}
	for i := range want {
		if gcd[i] != want[i] {
			t.Errorf("row %d:\n want %q\n got  %q", i, want[i], gcd[i])
		}
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

// TestStatsSparseKeepsAllSections confirms -x/--sparse does NOT suppress whole
// sections: upstream stats.c consults `sparse` only at stats.c:1796, where it
// thins all-zero rows of the IS section. Every other section must be emitted
// identically with and without --sparse, and the non-IS bytes must be equal.
func TestStatsSparseKeepsAllSections(t *testing.T) {
	read := func(opts StatsOptions) string {
		in, err := os.Open(statsFixture(t, "1_map_cigar.sam"))
		if err != nil {
			t.Fatalf("open input: %v", err)
		}
		defer in.Close()
		var out bytes.Buffer
		if err := Stats(in, &out, opts); err != nil {
			t.Fatalf("Stats: %v", err)
		}
		return out.String()
	}
	dense := read(StatsOptions{})
	sparse := read(StatsOptions{Sparse: true})

	// Every non-IS section must survive --sparse unchanged.
	for _, sec := range []string{"CHK", "SN", "FFQ", "LFQ", "GCF", "GCL", "GCC", "GCT", "FBC", "FTC", "LBC", "LTC", "RL", "FRL", "LRL", "MAPQ", "IC", "ID", "COV", "GCD"} {
		if got := extractSection(sparse, sec); got == "" && extractSection(dense, sec) != "" {
			t.Errorf("--sparse must NOT suppress %s section", sec)
		}
		if extractSection(sparse, sec) != extractSection(dense, sec) {
			t.Errorf("%s section changed under --sparse", sec)
		}
	}

	// The IS section must still be present but thinned: fixture 1 has a
	// single pair at insert size 100, so dense emits 101 IS rows (0..100)
	// while sparse emits only the one non-zero row.
	denseIS := strings.Count(extractSection(dense, "IS"), "\n")
	sparseIS := strings.Count(extractSection(sparse, "IS"), "\n")
	if denseIS != 101 {
		t.Errorf("dense IS rows = %d, want 101", denseIS)
	}
	if sparseIS != 1 {
		t.Errorf("sparse IS rows = %d, want 1", sparseIS)
	}
	if !strings.Contains(extractSection(sparse, "IS"), "IS\t100\t1\t1\t0\t0\n") {
		t.Errorf("sparse IS missing the single non-zero row:\n%s", extractSection(sparse, "IS"))
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

// TestStatsSparseThinsISOnly confirms the -x/--sparse flag thins all-zero IS
// rows but keeps every section header, including RL — upstream's `sparse`
// touches the IS row loop only.
func TestStatsSparseThinsISOnly(t *testing.T) {
	hdr := "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:1000\n"
	// Both mates of one inward pair (insert size 100) — upstream classifies
	// each mate then halves, so a complete pair is needed to yield one row.
	body := "r1\t99\tc\t1\t40\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\n" +
		"r1\t147\tc\t100\t40\t10M\t=\t1\t-100\tACGTACGTAC\tIIIIIIIIII\n"
	in := strings.NewReader(hdr + body)
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{Sparse: true}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if !strings.Contains(out.String(), "SN\t") {
		t.Fatalf("SN must still be emitted under --sparse")
	}
	if !strings.Contains(out.String(), "# Read lengths.") {
		t.Fatalf("--sparse must NOT suppress the RL section")
	}
	if !strings.Contains(out.String(), "# Insert sizes.") {
		t.Fatalf("--sparse must NOT suppress the IS section header")
	}
	// Only the single non-zero IS row survives the sparse thinning.
	if n := strings.Count(extractSection(out.String(), "IS"), "\n"); n != 1 {
		t.Fatalf("sparse IS rows = %d, want 1", n)
	}
}

// TestStatsTargetRegionsParity validates --target-regions against upstream's
// stat/11 fixtures. The SN, COV and GCD sections are compared byte-for-byte to
// the golden outputs for both the default coverage threshold and -g 4. The
// full report is not compared because upstream emits version-dependent
// sections (FBC/FTC/LBC/LTC barcode tables) and a dense IS table that this v1
// does not implement.
func TestStatsTargetRegionsParity(t *testing.T) {
	cases := []struct {
		name      string
		extraArgs []string // upstream `samtools stats` flags before the input
		opts      StatsOptions
	}{
		{
			"default-threshold",
			[]string{"-t", statsFixture(t, "11.stats.targets")},
			StatsOptions{TargetBED: statsFixture(t, "11.stats.targets")},
		},
		{
			"cov-threshold-4",
			[]string{"-g", "4", "-t", statsFixture(t, "11.stats.targets")},
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
			got := out.String()
			want := upstreamStats(t, statsFixture(t, "11_target.sam"), tc.extraArgs...)
			if extractSN(got) != extractSN(want) {
				t.Errorf("SN section differs\n--- want\n%s\n--- got\n%s",
					extractSN(want), extractSN(got))
			}
			if extractSection(got, "COV") != extractSection(want, "COV") {
				t.Errorf("COV rows differ\n--- want\n%s--- got\n%s",
					extractSection(want, "COV"), extractSection(got, "COV"))
			}
			if extractSection(got, "GCD") != extractSection(want, "GCD") {
				t.Errorf("GCD rows differ\n--- want\n%s--- got\n%s",
					extractSection(want, "GCD"), extractSection(got, "GCD"))
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

// TestStatsMPCParity compares our MPC mismatches-per-cycle section against
// upstream's stats expected outputs. The stat/1-8 golden files were generated
// with `-r test.fa`, so the section comment lines and every cycle row
// (including the bumped trailing all-zero cycle) must match byte-for-byte.
// Fixture 7_supp exercises the only non-zero entry across the suite: a
// supplementary read whose first aligned base mismatches the reference, landing
// in the N column because the "*" quality string makes qual+1 wrap to 0.
func TestStatsMPCParity(t *testing.T) {
	cases := []struct {
		name  string
		input string
	}{
		{"1_map_cigar", "1_map_cigar.sam"},
		{"2_equal_cigar", "2_equal_cigar_full_seq.sam"},
		{"5_insert_cigar", "5_insert_cigar.sam"},
		{"7_supp", "7_supp.sam"},
		{"8_secondary", "8_secondary.sam"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			in, err := os.Open(statsFixture(t, tc.input))
			if err != nil {
				t.Fatalf("open input: %v", err)
			}
			defer in.Close()
			var out bytes.Buffer
			opts := StatsOptions{RefSeq: statsFixture(t, "test.fa")}
			if err := Stats(in, &out, opts); err != nil {
				t.Fatalf("Stats: %v", err)
			}
			got := out.String()
			want := upstreamStats(t, statsFixture(t, tc.input), "-r", statsFixture(t, "test.fa"))
			if extractSection(got, "MPC") != extractSection(want, "MPC") {
				t.Errorf("MPC rows differ\n--- want\n%s--- got\n%s",
					extractSection(want, "MPC"), extractSection(got, "MPC"))
			}
			for _, prefix := range []string{
				"Mismatches per cycle and quality",
				"Columns correspond to qualities, rows to cycles",
				"is the number of N's and the rest is the number of mismatches",
			} {
				if extractCommentHeader(got, prefix) != extractCommentHeader(want, prefix) {
					t.Errorf("MPC header %q differs\nwant: %q\ngot:  %q", prefix,
						extractCommentHeader(want, prefix), extractCommentHeader(got, prefix))
				}
			}
		})
	}
}

// TestStatsMPCGatedOnRefSeq confirms the MPC section is emitted only when
// --ref-seq is supplied, matching upstream's mpc_buf non-nil gate.
func TestStatsMPCGatedOnRefSeq(t *testing.T) {
	in, err := os.Open(statsFixture(t, "1_map_cigar.sam"))
	if err != nil {
		t.Fatalf("open input: %v", err)
	}
	defer in.Close()
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if extractSection(out.String(), "MPC") != "" {
		t.Errorf("MPC section must be omitted without --ref-seq")
	}
}

// TestStatsRFSParity compares our RFS reference-statistics section against
// upstream's stats expected outputs. The validated invocations are upstream's
// test.pl stats cases:
//
//	16: samtools stats --ref-stats 11_target.sam            (no reference FASTA)
//	17: samtools stats --ref-stats 11_target.sam -r test1.fa
//	19: samtools stats --ref-stats 11_target.sam -r test1.fa -t 11.stats.targets
//
// Upstream test 18 (a command-line region argument against a BAM) is not
// reproduced here: the Go CLI does not yet accept positional region arguments,
// which would need a region-iterator implementation outside the scope of the
// RFS section itself.
func TestStatsRFSParity(t *testing.T) {
	cases := []struct {
		name string
		args func(t *testing.T) []string // upstream `samtools stats` flags
		opts func(t *testing.T) StatsOptions
	}{
		{
			name: "no_reference",
			args: func(t *testing.T) []string { return []string{"--ref-stats"} },
			opts: func(t *testing.T) StatsOptions { return StatsOptions{RefStats: true} },
		},
		{
			name: "with_reference",
			args: func(t *testing.T) []string { return []string{"--ref-stats", "-r", statsFixture(t, "test1.fa")} },
			opts: func(t *testing.T) StatsOptions {
				return StatsOptions{RefStats: true, RefSeq: statsFixture(t, "test1.fa")}
			},
		},
		{
			name: "with_reference_and_targets",
			args: func(t *testing.T) []string {
				return []string{"--ref-stats", "-r", statsFixture(t, "test1.fa"), "-t", statsFixture(t, "11.stats.targets")}
			},
			opts: func(t *testing.T) StatsOptions {
				return StatsOptions{
					RefStats:  true,
					RefSeq:    statsFixture(t, "test1.fa"),
					TargetBED: statsFixture(t, "11.stats.targets"),
				}
			},
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
			if err := Stats(in, &out, tc.opts(t)); err != nil {
				t.Fatalf("Stats: %v", err)
			}
			want := extractSection(upstreamStats(t, statsFixture(t, "11_target.sam"), tc.args(t)...), "RFS")
			if extractSection(out.String(), "RFS") != want {
				t.Errorf("RFS rows differ\n--- want\n%s--- got\n%s",
					want, extractSection(out.String(), "RFS"))
			}
		})
	}
}

// TestStatsRFSHeaders confirms the three RFS comment header lines are emitted
// verbatim from upstream stats.c:1897-1899.
func TestStatsRFSHeaders(t *testing.T) {
	in, err := os.Open(statsFixture(t, "11_target.sam"))
	if err != nil {
		t.Fatalf("open input: %v", err)
	}
	defer in.Close()
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{RefStats: true}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	for _, prefix := range []string{
		"Reference statistics. Use `grep ^RFS",
		"Total count, Output count, Average GC",
		"Sequence name, Length, GC content, Unknown count",
	} {
		if extractCommentHeader(out.String(), prefix) == "" {
			t.Errorf("missing RFS header line %q", prefix)
		}
	}
}

// TestStatsRFSGatedOnRefStats confirms the RFS section is emitted only when
// --ref-stats is supplied, matching upstream's rstat non-nil gate.
func TestStatsRFSGatedOnRefStats(t *testing.T) {
	in, err := os.Open(statsFixture(t, "11_target.sam"))
	if err != nil {
		t.Fatalf("open input: %v", err)
	}
	defer in.Close()
	var out bytes.Buffer
	if err := Stats(in, &out, StatsOptions{}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if extractSection(out.String(), "RFS") != "" {
		t.Errorf("RFS section must be omitted without --ref-stats")
	}
}
