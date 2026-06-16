package bedmerge

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"strings"
	"testing"
)

// errWriter is an io.Writer that always returns errFailingWrite on Write.
type errWriter struct {
	n int // bytes accepted before erroring; 0 means fail immediately
}

var errFailingWrite = errors.New("forced write failure")

func (e *errWriter) Write(p []byte) (int, error) {
	if e.n <= 0 {
		return 0, errFailingWrite
	}
	if len(p) <= e.n {
		e.n -= len(p)
		return len(p), nil
	}
	written := e.n
	e.n = 0
	return written, errFailingWrite
}

// errReader is an io.Reader that returns its err after returning its bytes.
type errReader struct {
	data []byte
	pos  int
	err  error
}

func (r *errReader) Read(p []byte) (int, error) {
	if r.pos >= len(r.data) {
		return 0, r.err
	}
	n := copy(p, r.data[r.pos:])
	r.pos += n
	return n, nil
}

// largeBedInput builds a BED text body with the given number of non-overlapping
// records on chromosome chr, so the output line count equals the input count.
func largeBedInput(chr string, n int) string {
	var sb strings.Builder
	for i := 0; i < n; i++ {
		fmt.Fprintf(&sb, "%s\t%d\t%d\n", chr, i*100, i*100+10)
	}
	return sb.String()
}

// TestMergeReadErrorMalformedRecord covers the chromStart parse-error path in
// the BED text reader by feeding a record with a non-numeric start.
func TestMergeReadErrorMalformedRecord(t *testing.T) {
	input := "chr1\tnotanumber\t200\tfoo\n" // 4 cols so BED detection can't apply
	var buf bytes.Buffer
	_, err := Merge(strings.NewReader(input), &buf, MergeOptions{})
	if err == nil {
		t.Fatalf("expected error for malformed record, got nil")
	}
}

// TestMergeSortTiebreakerByEnd covers the (chrom,start,end) sort comparator's
// end tiebreaker: two records share a start but differ in end.
func TestMergeSortTiebreakerByEnd(t *testing.T) {
	input := "chr1\t100\t150\nchr1\t100\t200\n"
	var buf bytes.Buffer
	count, err := Merge(strings.NewReader(input), &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	if count != 1 {
		t.Errorf("expected 1 merged interval, got %d", count)
	}
	if got := buf.String(); got != "chr1\t100\t200\n" {
		t.Errorf("Output mismatch: got %q", got)
	}
}

// TestMergeFlushError covers the bufio flush/write error in mergeRecords.
func TestMergeFlushError(t *testing.T) {
	input := largeBedInput("chr1", 1000)
	_, err := Merge(strings.NewReader(input), &errWriter{}, MergeOptions{})
	if err == nil {
		t.Fatalf("expected write/flush error, got nil")
	}
}

// TestPositionGroupsEmpty exercises the empty-slice early return.
func TestPositionGroupsEmpty(t *testing.T) {
	if got := positionGroups(nil, MergeOptions{}); got != nil {
		t.Errorf("positionGroups(nil) = %v, want nil", got)
	}
}

// TestBedGraphTakesColumnFour verifies -g emits column 4 (the bedGraph score)
// of the first record in the merged group.
func TestBedGraphTakesColumnFour(t *testing.T) {
	input := "chr1\t100\t200\t10\nchr1\t150\t250\t20\nchr2\t10\t20\t5\n"
	expected := "chr1\t100\t250\t10\nchr2\t10\t20\t5\n"
	var buf bytes.Buffer
	count, err := Merge(strings.NewReader(input), &buf, MergeOptions{
		OutputFields: OutputFields{BedGraph: true},
	})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	if count != 2 {
		t.Errorf("expected 2 intervals, got %d", count)
	}
	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected: %q\nGot: %q", expected, buf.String())
	}
}

// TestCountFlag verifies --count emits the merged group size as the 4th column.
func TestCountFlag(t *testing.T) {
	input := "chr1\t10\t20\nchr1\t15\t30\nchr1\t40\t50\n"
	expected := "chr1\t10\t30\t2\nchr1\t40\t50\t1\n"
	var buf bytes.Buffer
	if _, err := Merge(strings.NewReader(input), &buf, MergeOptions{
		OutputFields: OutputFields{Count: true},
	}); err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected: %q\nGot: %q", expected, buf.String())
	}
}

// TestMergeWithStatsColumnOps covers the ColumnOps branch in MergeWithStats.
func TestMergeWithStatsColumnOps(t *testing.T) {
	co, err := ParseColumnOps("5", "sum")
	if err != nil {
		t.Fatalf("ParseColumnOps failed: %v", err)
	}
	input := "chr1\t10\t20\ta\t5\nchr1\t15\t30\tb\t7\n"
	var buf bytes.Buffer
	stats, err := MergeWithStats(strings.NewReader(input), &buf, MergeOptions{ColumnOps: co})
	if err != nil {
		t.Fatalf("MergeWithStats failed: %v", err)
	}
	if stats.OutputIntervals != 1 {
		t.Errorf("expected 1 output interval, got %d", stats.OutputIntervals)
	}
	if got := buf.String(); got != "chr1\t10\t30\t12\n" {
		t.Errorf("unexpected output: %q", got)
	}

	// A non-numeric value in a sum group now yields the null value plus a
	// warning, not an error (parity with upstream merge.t23).
	badInput := "chr1\t10\t20\ta\tx\nchr1\t15\t30\tb\t7\n"
	buf.Reset()
	var warn bytes.Buffer
	_, err = MergeWithStats(strings.NewReader(badInput), &buf, MergeOptions{ColumnOps: co, Warn: &warn})
	if err != nil {
		t.Fatalf("expected non-numeric sum to warn, not error: %v", err)
	}
	if got := buf.String(); got != "chr1\t10\t30\t.\n" {
		t.Errorf("expected null value output, got %q", got)
	}
	if !strings.Contains(warn.String(), "Non numeric value x in 5") {
		t.Errorf("expected non-numeric warning, got %q", warn.String())
	}
}

// TestMergeWithStatsReadError covers the malformed-record path in MergeWithStats.
func TestMergeWithStatsReadError(t *testing.T) {
	var buf bytes.Buffer
	_, err := MergeWithStats(strings.NewReader("chr1\tnotanumber\t200\tx\n"), &buf, MergeOptions{})
	if err == nil {
		t.Fatalf("expected read error, got nil")
	}
}

// TestMergeWithStatsEmpty covers the empty-input path in MergeWithStats.
func TestMergeWithStatsEmpty(t *testing.T) {
	var buf bytes.Buffer
	stats, err := MergeWithStats(strings.NewReader(""), &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("MergeWithStats failed: %v", err)
	}
	if stats.InputIntervals != 0 || stats.OutputIntervals != 0 || stats.MergedCount != 0 {
		t.Errorf("expected all-zero stats, got %+v", stats)
	}
	if buf.Len() != 0 {
		t.Errorf("expected no output, got %q", buf.String())
	}
}

// TestMergeWithStatsTiebreakersAcrossChroms covers the chrom comparison and the
// end tiebreaker of the sort across multiple chromosomes.
func TestMergeWithStatsTiebreakersAcrossChroms(t *testing.T) {
	input := "chr2\t10\t20\nchr1\t100\t150\nchr1\t100\t200\n"
	expected := "chr1\t100\t200\nchr2\t10\t20\n"
	var buf bytes.Buffer
	stats, err := MergeWithStats(strings.NewReader(input), &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("MergeWithStats failed: %v", err)
	}
	if stats.InputIntervals != 3 || stats.OutputIntervals != 2 || stats.MergedCount != 1 {
		t.Errorf("unexpected stats: %+v", stats)
	}
	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected: %q\nGot: %q", expected, buf.String())
	}
}

// TestMergeWithStatsWriteError covers the write-error branch in MergeWithStats.
func TestMergeWithStatsWriteError(t *testing.T) {
	input := largeBedInput("chr1", 1000)
	_, err := MergeWithStats(strings.NewReader(input), &errWriter{}, MergeOptions{})
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
}

// --- input/colops coverage ---------------------------------------------------

// TestMergeColumnOpsSkipsCommentAndTrackLines covers the comment/track/browser/
// blank line skipping in the text reader.
func TestMergeColumnOpsSkipsCommentAndTrackLines(t *testing.T) {
	input := "" +
		"# leading comment\n" +
		"track name=foo\n" +
		"browser position chr1:1-1000\n" +
		"\n" +
		"chr1\t10\t20\ta\t5\n" +
		"   \n" +
		"chr1\t15\t30\tb\t7\n"
	co, err := ParseColumnOps("5", "sum")
	if err != nil {
		t.Fatalf("ParseColumnOps failed: %v", err)
	}
	var buf bytes.Buffer
	n, err := Merge(strings.NewReader(input), &buf, MergeOptions{ColumnOps: co})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	if n != 1 {
		t.Errorf("expected 1 merged interval, got %d", n)
	}
	if got := buf.String(); got != "chr1\t10\t30\t12\n" {
		t.Errorf("Output mismatch: got %q", got)
	}
}

// TestMergeColumnOpsStrandAssignment covers strand-aware merge (-s) using the
// BED strand column to keep opposite-strand overlapping records separate.
func TestMergeColumnOpsStrandAssignment(t *testing.T) {
	input := "" +
		"chr1\t10\t20\tn1\t1\t+\n" +
		"chr1\t15\t25\tn2\t1\t-\n"
	co, _ := ParseColumnOps("4", "distinct")
	var buf bytes.Buffer
	n, err := Merge(strings.NewReader(input), &buf, MergeOptions{ColumnOps: co, StrandSpec: true})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	if n != 2 {
		t.Errorf("expected 2 intervals (strand-separated), got %d", n)
	}
	expected := "chr1\t10\t20\tn1\nchr1\t15\t25\tn2\n"
	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected: %q\nGot: %q", expected, buf.String())
	}
}

// TestMergeColumnOpsMissingRequestedColumn verifies a -c column past the record
// width yields the null value (BAM-style), matching upstream's getColVal.
func TestMergeColumnOpsMissingRequestedColumn(t *testing.T) {
	input := "chr1\t10\t20\ta\n"
	co, _ := ParseColumnOps("6", "first")
	var buf bytes.Buffer
	if _, err := Merge(strings.NewReader(input), &buf, MergeOptions{ColumnOps: co}); err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	if got := buf.String(); got != "chr1\t10\t20\t\n" {
		t.Errorf("expected empty value for missing column, got %q", got)
	}
}

// TestMergeColumnOpsScannerError covers the scanner.Err() branch in readText.
func TestMergeColumnOpsScannerError(t *testing.T) {
	co, _ := ParseColumnOps("4", "first")
	r := &errReader{data: []byte("chr1\t10\t20\ta\n"), err: errors.New("boom")}
	var buf bytes.Buffer
	_, err := Merge(r, &buf, MergeOptions{ColumnOps: co})
	if err == nil {
		t.Fatalf("expected scanner error, got nil")
	}
}

// TestMergeColumnOpsEmptyInput covers the empty-input path.
func TestMergeColumnOpsEmptyInput(t *testing.T) {
	co, _ := ParseColumnOps("4", "first")
	var buf bytes.Buffer
	n, err := Merge(strings.NewReader(""), &buf, MergeOptions{ColumnOps: co})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	if n != 0 {
		t.Errorf("expected 0 intervals, got %d", n)
	}
	if buf.Len() != 0 {
		t.Errorf("expected no output, got %q", buf.String())
	}
}

// TestApplyOpMinMaxInteriorBranches covers the min/max inner-loop comparisons.
func TestApplyOpMinMaxInteriorBranches(t *testing.T) {
	if got, err := ApplyOp("min", 5, []string{"9", "3", "5"}); err != nil || got != "3" {
		t.Errorf("min: got %q, err %v; want 3", got, err)
	}
	if got, err := ApplyOp("max", 5, []string{"1", "8", "4"}); err != nil || got != "8" {
		t.Errorf("max: got %q, err %v; want 8", got, err)
	}
}

// TestApplyOpMedianEvenCount covers the even-length median branch.
func TestApplyOpMedianEvenCount(t *testing.T) {
	if got, err := ApplyOp("median", 5, []string{"2", "4", "6", "8"}); err != nil || got != "5" {
		t.Errorf("median(2,4,6,8) = %q, err %v; want 5", got, err)
	}
	if got, err := ApplyOp("median", 5, []string{"4", "7"}); err != nil || got != "5.5" {
		t.Errorf("median(4,7) = %q, err %v; want 5.5", got, err)
	}
}

// TestApplyOpNonNumericReturnsError confirms the standalone ApplyOp helper still
// errors on a non-numeric value (distinct from the merge path which warns).
func TestApplyOpNonNumericReturnsError(t *testing.T) {
	if _, err := ApplyOp("sum", 4, []string{"a"}); err == nil {
		t.Fatalf("expected non-numeric error from ApplyOp, got nil")
	}
}

// TestApplyOpUnsupportedOp covers the default fallthrough.
func TestApplyOpUnsupportedOp(t *testing.T) {
	if got, _ := ApplyOp("bogus", 4, []string{"a"}); got != "" {
		t.Errorf("expected empty result for unsupported op, got %q", got)
	}
}

// TestModeAntimodeNumeric covers mode/antimode tie handling on numeric input.
func TestModeAntimodeNumeric(t *testing.T) {
	if got := modeOrAntimodeNum([]float64{1, 1, 2}, false, DefaultPrecision); got != "2" {
		t.Errorf("antimode([1,1,2]) = %q, want 2", got)
	}
	if got := modeOrAntimodeNum([]float64{1, 1, 2}, true, DefaultPrecision); got != "1" {
		t.Errorf("mode([1,1,2]) = %q, want 1", got)
	}
}

// Compile-time sanity: assert errReader satisfies io.Reader.
var _ io.Reader = (*errReader)(nil)
