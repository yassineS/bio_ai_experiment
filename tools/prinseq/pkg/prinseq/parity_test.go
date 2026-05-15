package prinseq

// Byte-for-byte parity tests for prinseq's Go port against the upstream
// PRINSEQ-lite Perl reference (uwb-linux/prinseq @ 0.20.4). The upstream
// binary is `perl reference_code/prinseq/prinseq-lite.pl`. Each fixture
// under tools/prinseq/testdata/parity/ was generated once by running the
// upstream binary on a small representative corpus and capturing
// stdout / the `-out_good` file.
//
// Goal: every test that doesn't t.Skip should byte-match upstream's
// output. Where the upstream output format differs from ours
// (stats summary text vs structured numbers), we parse the upstream
// fixture and compare numbers rather than literal bytes.

import (
	"bufio"
	"bytes"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// readParity reads a fixture from tools/prinseq/testdata/parity/.
func readParity(t *testing.T, name string) []byte {
	t.Helper()
	path := filepath.Join("..", "..", "testdata", "parity", name)
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read parity fixture %s: %v", name, err)
	}
	return data
}

func mustEqualBytes(t *testing.T, label string, got, want []byte) {
	t.Helper()
	if !bytes.Equal(got, want) {
		t.Fatalf("%s mismatch.\nwant (%d bytes):\n%s\ngot (%d bytes):\n%s", label, len(want), want, len(got), got)
	}
}

// runFilter drives the Go Filter function with the given options on the
// named fixture and returns the resulting bytes.
func runFilter(t *testing.T, inputName string, isFastq bool, opts FilterOptions) []byte {
	t.Helper()
	in := readParity(t, inputName)
	var out bytes.Buffer
	if err := Filter(bytes.NewReader(in), &out, isFastq, opts); err != nil {
		t.Fatalf("Filter(%s): %v", inputName, err)
	}
	return out.Bytes()
}

// runFilterPaired drives the paired-end Filter and returns (r1, r2) bytes.
func runFilterPaired(t *testing.T, in1Name, in2Name string, isFastq bool, opts FilterOptions) ([]byte, []byte) {
	t.Helper()
	in1 := readParity(t, in1Name)
	in2 := readParity(t, in2Name)
	var o1, o2 bytes.Buffer
	if err := FilterPaired(bytes.NewReader(in1), bytes.NewReader(in2), &o1, &o2, isFastq, opts); err != nil {
		t.Fatalf("FilterPaired: %v", err)
	}
	return o1.Bytes(), o2.Bytes()
}

// parseUpstreamStats parses the upstream `prinseq-lite.pl -stats_info
// -stats_len` text output into a flat key->string map keyed on the
// section + name (e.g. "info.bases"). Sample input lines look like:
//
//	stats_info	bases	1150
//	stats_len	min	50
//
// The string value is returned verbatim so the caller can decide
// whether to compare as an integer or a float.
func parseUpstreamStats(t *testing.T, blob []byte) map[string]string {
	t.Helper()
	out := make(map[string]string)
	scanner := bufio.NewScanner(bytes.NewReader(blob))
	for scanner.Scan() {
		parts := strings.Split(scanner.Text(), "\t")
		if len(parts) < 3 {
			continue
		}
		section := strings.TrimPrefix(parts[0], "stats_")
		out[section+"."+parts[1]] = parts[2]
	}
	if err := scanner.Err(); err != nil {
		t.Fatalf("scan upstream stats: %v", err)
	}
	return out
}

// TestParity_Prinseq_Stats_FastaSmall verifies that our Go
// CalculateStats numbers match upstream `-stats_info -stats_len`
// for the small.fasta corpus.
func TestParity_Prinseq_Stats_FastaSmall(t *testing.T) {
	in := readParity(t, "small.fasta")
	stats, err := CalculateStats(bytes.NewReader(in), false)
	if err != nil {
		t.Fatalf("CalculateStats: %v", err)
	}
	want := parseUpstreamStats(t, readParity(t, "stats_info_len_fa.expected.txt"))

	mustAtoiEq(t, "info.bases", want, stats.TotalBases)
	mustAtoiEq(t, "info.reads", want, stats.NumReads)
	mustAtoiEq(t, "len.min", want, stats.MinLength)
	mustAtoiEq(t, "len.max", want, stats.MaxLength)
	// mean is a float formatted to 2 decimal places.
	if got := strconv.FormatFloat(stats.AvgLength, 'f', 2, 64); got != want["len.mean"] {
		t.Errorf("len.mean: got %s, want %s", got, want["len.mean"])
	}
}

// TestParity_Prinseq_Stats_FastqSmall verifies the same numbers for a
// FASTQ input (under default Phred+33).
func TestParity_Prinseq_Stats_FastqSmall(t *testing.T) {
	in := readParity(t, "small.fastq")
	stats, err := CalculateStats(bytes.NewReader(in), true)
	if err != nil {
		t.Fatalf("CalculateStats: %v", err)
	}
	want := parseUpstreamStats(t, readParity(t, "stats_info_len_fq.expected.txt"))

	mustAtoiEq(t, "info.bases", want, stats.TotalBases)
	mustAtoiEq(t, "info.reads", want, stats.NumReads)
	mustAtoiEq(t, "len.min", want, stats.MinLength)
	mustAtoiEq(t, "len.max", want, stats.MaxLength)
}

// mustAtoiEq compares an integer Go value against the same-keyed
// upstream stats value.
func mustAtoiEq(t *testing.T, key string, want map[string]string, got int) {
	t.Helper()
	w, err := strconv.Atoi(want[key])
	if err != nil {
		t.Fatalf("upstream stats key %s missing or non-integer: %q (%v)", key, want[key], err)
	}
	if w != got {
		t.Errorf("%s: got %d, want %d", key, got, w)
	}
}

// TestParity_Prinseq_Filter_MinLen verifies `-min_len 10` keeps only
// records >= 10 bases and discards the rest, with byte-for-byte FASTA
// output matching upstream `-line_width 0`.
func TestParity_Prinseq_Filter_MinLen(t *testing.T) {
	got := runFilter(t, "small.fasta", false, FilterOptions{MinLen: 10})
	want := readParity(t, "filt_min_len_10.expected.fasta")
	mustEqualBytes(t, "filt min_len 10", got, want)
}

// TestParity_Prinseq_Filter_MaxLen exercises the upper bound on length.
func TestParity_Prinseq_Filter_MaxLen(t *testing.T) {
	got := runFilter(t, "small.fasta", false, FilterOptions{MaxLen: 20})
	want := readParity(t, "filt_max_len_20.expected.fasta")
	mustEqualBytes(t, "filt max_len 20", got, want)
}

// TestParity_Prinseq_Filter_MinGC verifies `-min_gc 50`.
func TestParity_Prinseq_Filter_MinGC(t *testing.T) {
	got := runFilter(t, "small.fasta", false, FilterOptions{MinGC: 50})
	want := readParity(t, "filt_min_gc_50.expected.fasta")
	mustEqualBytes(t, "filt min_gc 50", got, want)
}

// TestParity_Prinseq_Filter_MaxGC verifies `-max_gc 60`.
func TestParity_Prinseq_Filter_MaxGC(t *testing.T) {
	got := runFilter(t, "small.fasta", false, FilterOptions{MaxGC: 60})
	want := readParity(t, "filt_max_gc_60.expected.fasta")
	mustEqualBytes(t, "filt max_gc 60", got, want)
}

// TestParity_Prinseq_Filter_NsMaxPercent verifies `-ns_max_p 5`.
func TestParity_Prinseq_Filter_NsMaxPercent(t *testing.T) {
	got := runFilter(t, "small.fasta", false, FilterOptions{MaxNsP: 5})
	want := readParity(t, "filt_ns_max_p_5.expected.fasta")
	mustEqualBytes(t, "filt ns_max_p 5", got, want)
}

// TestParity_Prinseq_Filter_NsMaxNumber verifies `-ns_max_n 2`.
func TestParity_Prinseq_Filter_NsMaxNumber(t *testing.T) {
	got := runFilter(t, "small.fasta", false, FilterOptions{MaxNsN: 2})
	want := readParity(t, "filt_ns_max_n_2.expected.fasta")
	mustEqualBytes(t, "filt ns_max_n 2", got, want)
}

// TestParity_Prinseq_Filter_MinQualMean verifies `-min_qual_mean 15`
// against a FASTQ input under default Phred+33.
func TestParity_Prinseq_Filter_MinQualMean(t *testing.T) {
	got := runFilter(t, "small.fastq", true, FilterOptions{MinQualMean: 15})
	want := readParity(t, "filt_min_qual_mean_15.expected.fastq")
	mustEqualBytes(t, "filt min_qual_mean 15", got, want)
}

// TestParity_Prinseq_Filter_Phred64MinQualMean verifies that mean
// quality is correctly decoded under Phred+64 (i.e. QualType:
// "illumina"). The fix for this is in this PR: the filter loop used
// to call calculateAvgQualityScore() which hard-coded offset 33,
// silently misclassifying Phred+64 reads. See
// tools/PARITY_VALIDATION.md.
func TestParity_Prinseq_Filter_Phred64MinQualMean(t *testing.T) {
	got := runFilter(t, "p64.fastq", true, FilterOptions{MinQualMean: 39, QualType: "illumina"})
	want := readParity(t, "filt_p64_min_qual_mean_39.expected.fastq")
	mustEqualBytes(t, "filt p64 min_qual_mean 39", got, want)
}

// TestParity_Prinseq_Filter_MultiCriteria exercises the union of a
// length window and a GC window on FASTA.
func TestParity_Prinseq_Filter_MultiCriteria(t *testing.T) {
	got := runFilter(t, "small.fasta", false, FilterOptions{MinLen: 10, MaxLen: 25, MinGC: 30, MaxGC: 60})
	want := readParity(t, "filt_multi.expected.fasta")
	mustEqualBytes(t, "filt multi", got, want)
}

// TestParity_Prinseq_Filter_TrimLeft verifies `-trim_left 5`.
func TestParity_Prinseq_Filter_TrimLeft(t *testing.T) {
	got := runFilter(t, "trim.fastq", true, FilterOptions{TrimLeft: 5})
	want := readParity(t, "trim_left_5.expected.fastq")
	mustEqualBytes(t, "trim_left 5", got, want)
}

// TestParity_Prinseq_Filter_TrimRight verifies `-trim_right 4`.
func TestParity_Prinseq_Filter_TrimRight(t *testing.T) {
	got := runFilter(t, "trim.fastq", true, FilterOptions{TrimRight: 4})
	want := readParity(t, "trim_right_4.expected.fastq")
	mustEqualBytes(t, "trim_right 4", got, want)
}

// TestParity_Prinseq_Filter_TrimQualLeft verifies `-trim_qual_left 20`
// which trims the 5' end while the per-base quality is below 20.
func TestParity_Prinseq_Filter_TrimQualLeft(t *testing.T) {
	got := runFilter(t, "trim.fastq", true, FilterOptions{TrimQualL: 20})
	want := readParity(t, "trim_qual_left_20.expected.fastq")
	mustEqualBytes(t, "trim_qual_left 20", got, want)
}

// TestParity_Prinseq_Filter_TrimQualRight verifies `-trim_qual_right 20`.
func TestParity_Prinseq_Filter_TrimQualRight(t *testing.T) {
	got := runFilter(t, "trim.fastq", true, FilterOptions{TrimQualR: 20})
	want := readParity(t, "trim_qual_right_20.expected.fastq")
	mustEqualBytes(t, "trim_qual_right 20", got, want)
}

// TestParity_Prinseq_Filter_TrimTailLeft exercises the poly-A/T head
// trim. Upstream only collapses a single homopolymer at a time (A or
// T, not both) — a quirk we previously got wrong and fixed in this PR.
func TestParity_Prinseq_Filter_TrimTailLeft(t *testing.T) {
	got := runFilter(t, "trim.fastq", true, FilterOptions{TrimTailLeft: 4})
	want := readParity(t, "trim_tail_left_4.expected.fastq")
	mustEqualBytes(t, "trim_tail_left 4", got, want)
}

// TestParity_Prinseq_Filter_TrimTailRight is the 3'-end counterpart.
func TestParity_Prinseq_Filter_TrimTailRight(t *testing.T) {
	got := runFilter(t, "trim.fastq", true, FilterOptions{TrimTailRight: 4})
	want := readParity(t, "trim_tail_right_4.expected.fastq")
	mustEqualBytes(t, "trim_tail_right 4", got, want)
}

// TestParity_Prinseq_Filter_Derep verifies that exact-duplicate removal
// (-derep 1) drops every record after the first occurrence of each
// unique sequence.
func TestParity_Prinseq_Filter_Derep(t *testing.T) {
	got := runFilter(t, "dups.fasta", false, FilterOptions{Derep: 1, DerepMin: 2})
	want := readParity(t, "derep_1.expected.fasta")
	mustEqualBytes(t, "derep 1", got, want)
}

// TestParity_Prinseq_FilterPaired_MinLen exercises the paired-end
// path: only records where both mates pass the filter survive, with
// the same ordering as upstream.
func TestParity_Prinseq_FilterPaired_MinLen(t *testing.T) {
	r1, r2 := runFilterPaired(t, "pe1.fastq", "pe2.fastq", true, FilterOptions{MinLen: 10})
	mustEqualBytes(t, "PE R1", r1, readParity(t, "pe_r1.expected.fastq"))
	mustEqualBytes(t, "PE R2", r2, readParity(t, "pe_r2.expected.fastq"))
}

// TestParity_Prinseq_Empty_NoCrash ensures both Filter and CalculateStats
// tolerate an empty input file. Upstream prints "Input sequences: 0"
// and exits cleanly; we should not panic, and we should produce empty
// output.
func TestParity_Prinseq_Empty_NoCrash(t *testing.T) {
	for _, fixture := range []string{"empty.fasta", "empty.fastq"} {
		fixture := fixture
		t.Run(fixture, func(t *testing.T) {
			in := readParity(t, fixture)
			isFQ := strings.HasSuffix(fixture, ".fastq")
			// Filter
			var out bytes.Buffer
			if err := Filter(bytes.NewReader(in), &out, isFQ, FilterOptions{MinLen: 1}); err != nil {
				t.Errorf("Filter on %s: %v", fixture, err)
			}
			if out.Len() != 0 {
				t.Errorf("Filter on empty input produced %d bytes, want 0", out.Len())
			}
			// CalculateStats
			if _, err := CalculateStats(bytes.NewReader(in), isFQ); err != nil {
				t.Errorf("CalculateStats on %s: %v", fixture, err)
			}
		})
	}
}
