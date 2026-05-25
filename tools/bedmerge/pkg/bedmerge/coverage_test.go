package bedmerge

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
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

// errReader is an io.Reader that returns errFailingRead after returning n bytes.
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

// largeBedInput builds a BED text body with the given number of records on
// chromosome chr. Each record is non-overlapping so no merging occurs and the
// number of output lines matches the input.
func largeBedInput(chr string, n int) string {
	var sb strings.Builder
	for i := 0; i < n; i++ {
		fmt.Fprintf(&sb, "%s\t%d\t%d\n", chr, i*100, i*100+10)
	}
	return sb.String()
}

// TestMergeReadErrorMalformedRecord covers the bedReader.Read error path in
// Merge (bedmerge.go:55-57) by feeding a record with a non-numeric start.
func TestMergeReadErrorMalformedRecord(t *testing.T) {
	input := "chr1\tnotanumber\t200\n"
	var buf bytes.Buffer
	_, err := Merge(strings.NewReader(input), &buf, MergeOptions{})
	if err == nil {
		t.Fatalf("expected error for malformed BED record, got nil")
	}
	if !strings.Contains(err.Error(), "error reading BED record") {
		t.Errorf("expected wrapped read error, got: %v", err)
	}
}

// TestMergeSortTiebreakerByEnd covers the third-level sort comparator in Merge
// (bedmerge.go:84). With identical chrom and identical chromStart, the
// comparator falls through to comparing ChromEnd.
func TestMergeSortTiebreakerByEnd(t *testing.T) {
	// Two records with the same start but different ends. Output should be a
	// single merged interval spanning the wider record; the comparator code
	// path is exercised by the duplicate-start input.
	input := "chr1\t100\t150\nchr1\t100\t200\n"
	var buf bytes.Buffer
	count, err := Merge(strings.NewReader(input), &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	if count != 1 {
		t.Errorf("expected 1 merged interval, got %d", count)
	}
	expected := "chr1\t100\t200\n"
	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected: %q\nGot: %q", expected, buf.String())
	}
}

// TestMergeWriteError covers the writeIntervals error return in Merge
// (bedmerge.go:91-93). It uses a failing writer with enough records to
// overflow the internal bufio buffer, forcing bedWriter.Write to surface the
// error.
func TestMergeWriteError(t *testing.T) {
	// 1000 non-overlapping records, each about 20 bytes -> ~20kB output, well
	// over the default 4kB bufio buffer used inside bed.Writer.
	input := largeBedInput("chr1", 1000)
	_, err := Merge(strings.NewReader(input), &errWriter{}, MergeOptions{})
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
	if !strings.Contains(err.Error(), "error writing BED record") && !strings.Contains(err.Error(), "error flushing output") {
		t.Errorf("expected writeIntervals error, got: %v", err)
	}
}

// TestMergeIntervalsEmpty directly exercises the early-return branch in
// mergeIntervals (bedmerge.go:106-108) by passing a nil slice.
func TestMergeIntervalsEmpty(t *testing.T) {
	got := mergeIntervals(nil, MergeOptions{})
	if got != nil {
		t.Errorf("mergeIntervals(nil) = %v, want nil", got)
	}
}

// TestWriteIntervalsFlushError covers the flush-error branch in
// writeIntervals (bedmerge.go:204-206). The single record fits inside the
// bufio buffer so the underlying Write isn't called until Flush.
func TestWriteIntervalsFlushError(t *testing.T) {
	merged := []mergedInterval{{
		Record: &bed.Record{Chrom: "chr1", ChromStart: 1, ChromEnd: 2},
		count:  1,
	}}
	err := writeIntervals(&errWriter{}, merged, MergeOptions{})
	if err == nil {
		t.Fatalf("expected flush error, got nil")
	}
	if !strings.Contains(err.Error(), "error flushing output") {
		t.Errorf("expected flushing-output error, got: %v", err)
	}
}

// TestWriteBedGraphWriteError covers the writeBedGraph write-error branch
// (bedmerge.go:216-218). writeBedGraph writes directly to the underlying
// writer (no bufio), so a single record is enough.
func TestWriteBedGraphWriteError(t *testing.T) {
	merged := []mergedInterval{{
		Record: &bed.Record{Chrom: "chr1", ChromStart: 1, ChromEnd: 2, Score: 7},
		count:  1,
	}}
	err := writeBedGraph(&errWriter{}, merged)
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
	if !strings.Contains(err.Error(), "error writing bedGraph record") {
		t.Errorf("expected bedGraph write error, got: %v", err)
	}
}

// TestStreamingMergeEmpty covers the len(chromIntervals)==0 early return in
// the streamingMerge flushChrom closure (bedmerge.go:237-239), which is hit
// on a completely empty input via the final flush call (bedmerge.go:291).
func TestStreamingMergeEmpty(t *testing.T) {
	var buf bytes.Buffer
	count, err := Merge(strings.NewReader(""), &buf, MergeOptions{Streaming: true})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	if count != 0 {
		t.Errorf("expected 0 intervals, got %d", count)
	}
	if buf.Len() != 0 {
		t.Errorf("expected empty output, got %q", buf.String())
	}
}

// TestStreamingMergeSortTiebreakerByEnd covers the streaming sort comparator's
// ChromEnd tiebreaker (bedmerge.go:246) by giving two records on the same
// chromosome with identical starts but different ends.
func TestStreamingMergeSortTiebreakerByEnd(t *testing.T) {
	input := "chr1\t100\t150\nchr1\t100\t200\n"
	var buf bytes.Buffer
	count, err := Merge(strings.NewReader(input), &buf, MergeOptions{Streaming: true})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	if count != 1 {
		t.Errorf("expected 1 merged interval, got %d", count)
	}
	expected := "chr1\t100\t200\n"
	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected: %q\nGot: %q", expected, buf.String())
	}
}

// TestStreamingMergeWriteError covers the writeIntervals error propagation in
// flushChrom (bedmerge.go:251-253) and the chromosome-change flush error
// path (bedmerge.go:281-283) by switching chromosomes mid-stream while
// writing to a failing writer.
func TestStreamingMergeWriteError(t *testing.T) {
	// Enough records on chr1 so the first flushChrom (at the chr1 -> chr2
	// boundary) overflows bufio and surfaces the error. ~1000 records on chr1
	// produces well over 4kB.
	input := largeBedInput("chr1", 1000) + "chr2\t10\t20\n"
	_, err := Merge(strings.NewReader(input), &errWriter{}, MergeOptions{Streaming: true})
	if err == nil {
		t.Fatalf("expected streaming write error, got nil")
	}
}

// TestStreamingMergeFinalFlushError covers the final-flush error branch in
// streamingMerge (bedmerge.go:291-293) where the loop ends without a
// chromosome change. Many records on a single chromosome trigger bufio's
// downstream flush during the final write.
func TestStreamingMergeFinalFlushError(t *testing.T) {
	input := largeBedInput("chr1", 1000)
	_, err := Merge(strings.NewReader(input), &errWriter{}, MergeOptions{Streaming: true})
	if err == nil {
		t.Fatalf("expected final-flush write error, got nil")
	}
}

// TestStreamingMergeReadError covers the read-error path in streamingMerge
// (bedmerge.go:265-267) by passing a malformed BED record.
func TestStreamingMergeReadError(t *testing.T) {
	input := "chr1\tnotanumber\t200\n"
	var buf bytes.Buffer
	_, err := Merge(strings.NewReader(input), &buf, MergeOptions{Streaming: true})
	if err == nil {
		t.Fatalf("expected read error, got nil")
	}
	if !strings.Contains(err.Error(), "error reading BED record") {
		t.Errorf("expected wrapped read error, got: %v", err)
	}
}

// TestStreamingMergeBedGraphParsesName covers the bedGraph-name-as-score
// parsing block in streamingMerge (bedmerge.go:270-276) by feeding 4-column
// records whose name field is a number.
func TestStreamingMergeBedGraphParsesName(t *testing.T) {
	input := "chr1\t100\t200\t10\nchr1\t150\t250\t20\nchr2\t10\t20\t5\n"
	expected := "chr1\t100\t250\t10\nchr2\t10\t20\t5\n"
	var buf bytes.Buffer
	count, err := Merge(strings.NewReader(input), &buf, MergeOptions{
		Streaming:    true,
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

// TestMergeWithStatsColumnOpsError covers the ColumnOps branch in
// MergeWithStats (bedmerge.go:308-313): both the success return and the
// error-from-mergeWithColumnOps path.
func TestMergeWithStatsColumnOps(t *testing.T) {
	// Success path: stats only contain OutputIntervals.
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

	// Error path: non-numeric values for sum.
	badInput := "chr1\t10\t20\ta\tx\nchr1\t15\t30\tb\t7\n"
	buf.Reset()
	_, err = MergeWithStats(strings.NewReader(badInput), &buf, MergeOptions{ColumnOps: co})
	if err == nil {
		t.Fatalf("expected error for non-numeric sum value, got nil")
	}
}

// TestMergeWithStatsStreaming covers the streaming branch in MergeWithStats
// (bedmerge.go:317-324), including the error path returned from
// streamingMerge.
func TestMergeWithStatsStreaming(t *testing.T) {
	input := "chr1\t100\t200\nchr1\t150\t250\nchr2\t100\t200\n"
	var buf bytes.Buffer
	stats, err := MergeWithStats(strings.NewReader(input), &buf, MergeOptions{Streaming: true})
	if err != nil {
		t.Fatalf("MergeWithStats failed: %v", err)
	}
	if stats.OutputIntervals != 2 {
		t.Errorf("expected 2 output intervals, got %d", stats.OutputIntervals)
	}
	if got := buf.String(); got != "chr1\t100\t250\nchr2\t100\t200\n" {
		t.Errorf("unexpected output: %q", got)
	}

	// Error path: a malformed record reaches streamingMerge.
	buf.Reset()
	_, err = MergeWithStats(strings.NewReader("chr1\tnotanumber\t200\n"), &buf, MergeOptions{Streaming: true})
	if err == nil {
		t.Fatalf("expected streaming error, got nil")
	}
}

// TestMergeWithStatsReadError covers the malformed-record read-error path in
// MergeWithStats (bedmerge.go:336-338).
func TestMergeWithStatsReadError(t *testing.T) {
	var buf bytes.Buffer
	_, err := MergeWithStats(strings.NewReader("chr1\tnotanumber\t200\n"), &buf, MergeOptions{})
	if err == nil {
		t.Fatalf("expected read error, got nil")
	}
	if !strings.Contains(err.Error(), "error reading BED record") {
		t.Errorf("expected wrapped read error, got: %v", err)
	}
}

// TestMergeWithStatsBedGraphParsesName covers the bedGraph-name-as-score
// parsing block in MergeWithStats (bedmerge.go:341-347).
func TestMergeWithStatsBedGraphParsesName(t *testing.T) {
	input := "chr1\t100\t200\t10\nchr1\t150\t250\t20\n"
	expected := "chr1\t100\t250\t10\n"
	var buf bytes.Buffer
	stats, err := MergeWithStats(strings.NewReader(input), &buf, MergeOptions{
		OutputFields: OutputFields{BedGraph: true},
	})
	if err != nil {
		t.Fatalf("MergeWithStats failed: %v", err)
	}
	if stats.OutputIntervals != 1 {
		t.Errorf("expected 1 output interval, got %d", stats.OutputIntervals)
	}
	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected: %q\nGot: %q", expected, buf.String())
	}
}

// TestMergeWithStatsEmpty covers the empty-input early-return branch in
// MergeWithStats (bedmerge.go:357-359).
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

// TestMergeWithStatsSortTiebreakersAcrossChroms covers both the
// chrom-comparison and ChromEnd-tiebreaker branches of the MergeWithStats
// sort (bedmerge.go:363-365 and 369). Records span multiple chromosomes and
// include duplicate starts within one chromosome.
func TestMergeWithStatsSortTiebreakersAcrossChroms(t *testing.T) {
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

// TestMergeWithStatsWriteError covers the writeIntervals error branch in
// MergeWithStats (bedmerge.go:378-380).
func TestMergeWithStatsWriteError(t *testing.T) {
	input := largeBedInput("chr1", 1000)
	_, err := MergeWithStats(strings.NewReader(input), &errWriter{}, MergeOptions{})
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
}

// --- colops.go coverage ------------------------------------------------------

// TestMergeColumnOpsSkipsCommentAndTrackLines covers the skip-empty/comment/
// track/browser branch (colops.go:126-127).
func TestMergeColumnOpsSkipsCommentAndTrackLines(t *testing.T) {
	input := "" +
		"# leading comment\n" +
		"track name=foo\n" +
		"browser position chr1:1-1000\n" +
		"\n" +
		"chr1\t10\t20\ta\t5\n" +
		"   \n" + // whitespace-only line
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

// TestMergeColumnOpsTooFewFields covers the <3-fields error path
// (colops.go:130-132).
func TestMergeColumnOpsTooFewFields(t *testing.T) {
	input := "chr1\t10\n"
	co, _ := ParseColumnOps("4", "first")
	var buf bytes.Buffer
	_, err := Merge(strings.NewReader(input), &buf, MergeOptions{ColumnOps: co})
	if err == nil {
		t.Fatalf("expected error for too-few-fields, got nil")
	}
	if !strings.Contains(err.Error(), "at least 3 fields") {
		t.Errorf("expected field-count error, got: %v", err)
	}
}

// TestMergeColumnOpsBadStart covers the chromStart parse error
// (colops.go:134-136).
func TestMergeColumnOpsBadStart(t *testing.T) {
	input := "chr1\tNaN\t20\ta\t5\n"
	co, _ := ParseColumnOps("5", "sum")
	var buf bytes.Buffer
	_, err := Merge(strings.NewReader(input), &buf, MergeOptions{ColumnOps: co})
	if err == nil {
		t.Fatalf("expected chromStart parse error, got nil")
	}
	if !strings.Contains(err.Error(), "invalid chromStart") {
		t.Errorf("expected chromStart error, got: %v", err)
	}
}

// TestMergeColumnOpsBadEnd covers the chromEnd parse error
// (colops.go:138-140).
func TestMergeColumnOpsBadEnd(t *testing.T) {
	input := "chr1\t10\tNaN\ta\t5\n"
	co, _ := ParseColumnOps("5", "sum")
	var buf bytes.Buffer
	_, err := Merge(strings.NewReader(input), &buf, MergeOptions{ColumnOps: co})
	if err == nil {
		t.Fatalf("expected chromEnd parse error, got nil")
	}
	if !strings.Contains(err.Error(), "invalid chromEnd") {
		t.Errorf("expected chromEnd error, got: %v", err)
	}
}

// TestMergeColumnOpsStrandAssignment covers the strand-assignment branch
// (colops.go:142-144) by using >5 columns plus strand-specific merging that
// keeps two records separate.
func TestMergeColumnOpsStrandAssignment(t *testing.T) {
	// Two overlapping intervals with different concrete strands; with
	// StrandSpec they must not merge, and the strand column (#6) must drive
	// the decision -- which only happens if colops.go:142-144 ran.
	input := "" +
		"chr1\t10\t20\tn1\t1\t+\n" +
		"chr1\t15\t25\tn2\t1\t-\n"
	co, _ := ParseColumnOps("4", "distinct")
	var buf bytes.Buffer
	n, err := Merge(strings.NewReader(input), &buf, MergeOptions{
		ColumnOps:  co,
		StrandSpec: true,
	})
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

// TestMergeColumnOpsMissingRequestedColumn covers the "requested column N but
// input has only M columns" error (colops.go:147-149).
func TestMergeColumnOpsMissingRequestedColumn(t *testing.T) {
	// Request column 6 but input is only 4 columns wide.
	input := "chr1\t10\t20\ta\n"
	co, _ := ParseColumnOps("6", "first")
	var buf bytes.Buffer
	_, err := Merge(strings.NewReader(input), &buf, MergeOptions{ColumnOps: co})
	if err == nil {
		t.Fatalf("expected missing-column error, got nil")
	}
	if !strings.Contains(err.Error(), "requested column 6") {
		t.Errorf("expected requested-column error, got: %v", err)
	}
}

// TestMergeColumnOpsScannerError covers the scanner.Err() branch
// (colops.go:159-161). A short input followed by a Read returning a non-EOF
// error makes bufio.Scanner.Err() non-nil.
func TestMergeColumnOpsScannerError(t *testing.T) {
	co, _ := ParseColumnOps("4", "first")
	r := &errReader{
		data: []byte("chr1\t10\t20\ta\n"),
		err:  errors.New("boom"),
	}
	var buf bytes.Buffer
	_, err := Merge(r, &buf, MergeOptions{ColumnOps: co})
	if err == nil {
		t.Fatalf("expected scanner error, got nil")
	}
	if !strings.Contains(err.Error(), "error reading BED input") {
		t.Errorf("expected wrapped scanner error, got: %v", err)
	}
}

// TestMergeColumnOpsEmptyInput covers the empty-intervals early return
// (colops.go:163-165).
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

// TestMergeColumnOpsSortTiebreakerByEnd covers the ChromEnd-tiebreaker branch
// in the colops sort (colops.go:174). Two records share chrom + start but
// differ in end.
func TestMergeColumnOpsSortTiebreakerByEnd(t *testing.T) {
	input := "chr1\t100\t150\ta\t1\nchr1\t100\t200\tb\t2\n"
	co, _ := ParseColumnOps("4", "distinct")
	var buf bytes.Buffer
	n, err := Merge(strings.NewReader(input), &buf, MergeOptions{ColumnOps: co})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	if n != 1 {
		t.Errorf("expected 1 merged interval, got %d", n)
	}
	if got := buf.String(); got != "chr1\t100\t200\ta,b\n" {
		t.Errorf("Output mismatch: got %q", got)
	}
}

// TestMergeColumnOpsFlushGroupWriteError covers the write-error branch inside
// flushGroup (colops.go:205-207) by writing a large number of records to a
// failing writer. The bufio.Writer used by mergeWithColumnOps surfaces the
// underlying-writer error once its buffer fills.
func TestMergeColumnOpsFlushGroupWriteError(t *testing.T) {
	// Non-overlapping records so each forms its own output line; >4kB total
	// output causes bufio to flush downstream and surface the error.
	var sb strings.Builder
	for i := 0; i < 1000; i++ {
		fmt.Fprintf(&sb, "chr1\t%d\t%d\ta%d\t%d\n", i*100, i*100+10, i, i)
	}
	co, _ := ParseColumnOps("4", "distinct")
	_, err := Merge(strings.NewReader(sb.String()), &errWriter{}, MergeOptions{ColumnOps: co})
	if err == nil {
		t.Fatalf("expected colops write error, got nil")
	}
}

// TestMergeColumnOpsApplyOpError covers the error-propagation branch inside
// flushGroup when applyOp returns an error (colops.go:234-236). A
// non-numeric value in the second interval of a merged group makes "sum" fail
// only when the group is flushed.
func TestMergeColumnOpsApplyOpErrorMidStream(t *testing.T) {
	// Two intervals on chr1 form a merged group (good + bad value), then a
	// third on chr2 forces the chr1 group to flush -- hitting the error path
	// inside the loop rather than the final flush.
	// Note: "NaN" parses as a valid float, so we use a clearly non-numeric
	// token instead.
	input := "chr1\t10\t20\ta\t5\nchr1\t15\t30\tb\tnope\nchr2\t1\t2\tc\t1\n"
	co, _ := ParseColumnOps("5", "sum")
	var buf bytes.Buffer
	_, err := Merge(strings.NewReader(input), &buf, MergeOptions{ColumnOps: co})
	if err == nil {
		t.Fatalf("expected applyOp error during flushGroup, got nil")
	}
	if !strings.Contains(err.Error(), "non-numeric") {
		t.Errorf("expected non-numeric error, got: %v", err)
	}
}

// TestMergeColumnOpsFlushFinalError covers the final w.Flush error branch
// (colops.go:247-249). Output is small enough that the failure only surfaces
// during the explicit flush at the end of mergeWithColumnOps.
func TestMergeColumnOpsFlushFinalError(t *testing.T) {
	input := "chr1\t10\t20\ta\t5\nchr1\t15\t30\tb\t7\n"
	co, _ := ParseColumnOps("5", "sum")
	_, err := Merge(strings.NewReader(input), &errWriter{}, MergeOptions{ColumnOps: co})
	if err == nil {
		t.Fatalf("expected final flush error, got nil")
	}
	if !strings.Contains(err.Error(), "error flushing output") {
		t.Errorf("expected flushing-output error, got: %v", err)
	}
}

// TestApplyOpMinMaxInteriorBranches covers the inner-loop comparison branches
// of "min" and "max" (colops.go:277-279 and the symmetric "max" branch). The
// first value must NOT be the extremum so the if-body actually runs.
func TestApplyOpMinMaxInteriorBranches(t *testing.T) {
	// Min: first value (9) larger than later one (3) -> inner-if taken.
	if got, err := applyOp("min", 5, []string{"9", "3", "5"}); err != nil || got != "3" {
		t.Errorf("min: got %q, err %v; want 3", got, err)
	}
	// Max: first value (1) smaller than later one (8) -> inner-if taken.
	if got, err := applyOp("max", 5, []string{"1", "8", "4"}); err != nil || got != "8" {
		t.Errorf("max: got %q, err %v; want 8", got, err)
	}
}

// TestApplyOpMedianEvenCount covers the even-length branch of "median"
// (colops.go:303-305). Three input intervals yield three values; merge them
// into a single group with the median across an *even* number of values by
// using applyOp directly.
func TestApplyOpMedianEvenCount(t *testing.T) {
	// Even-length: median = (4+6)/2 = 5.
	if got, err := applyOp("median", 5, []string{"2", "4", "6", "8"}); err != nil || got != "5" {
		t.Errorf("median(2,4,6,8) = %q, err %v; want 5", got, err)
	}
	// Even-length producing a non-integer: median = (4+7)/2 = 5.5.
	if got, err := applyOp("median", 5, []string{"4", "7"}); err != nil || got != "5.5" {
		t.Errorf("median(4,7) = %q, err %v; want 5.5", got, err)
	}
}

// TestApplyOpUnsupportedOpDefensive covers the defensive fallthrough error
// at the bottom of applyOp (colops.go:344). ParseColumnOps prevents this
// path under normal use, so we call applyOp directly.
func TestApplyOpUnsupportedOpDefensive(t *testing.T) {
	_, err := applyOp("bogus", 4, []string{"a"})
	if err == nil {
		t.Fatalf("expected unsupported-op error, got nil")
	}
	if !strings.Contains(err.Error(), "unsupported operation") {
		t.Errorf("expected unsupported-op error message, got: %v", err)
	}
}

// TestModeOrAntimodeDecrementBranch covers the antimode branch where a later
// value's count is strictly less than the current best (colops.go:367-369).
// The first value must therefore be more common than a later, distinct one.
func TestModeOrAntimodeDecrementBranch(t *testing.T) {
	// vals: "x" appears twice, "y" once. Order = [x, y]. best=x, bestCount=2.
	// When the loop visits "y" (count 1), the antimode branch must update
	// best to "y" because 1 < 2.
	if got := modeOrAntimode([]string{"x", "x", "y"}, false); got != "y" {
		t.Errorf("antimode([x,x,y]) = %q, want y", got)
	}
	// Sanity: same data with mode=true returns "x".
	if got := modeOrAntimode([]string{"x", "x", "y"}, true); got != "x" {
		t.Errorf("mode([x,x,y]) = %q, want x", got)
	}
}

// TestFlushColumnGroupEmpty covers the defensive empty-group early return in
// flushColumnGroup (colops.go). The main control flow in mergeWithColumnOps
// never passes an empty group, but the helper still guards against it.
func TestFlushColumnGroupEmpty(t *testing.T) {
	co := &ColumnOps{Columns: []int{4}, Ops: []string{"first"}}
	var buf bytes.Buffer
	n, err := flushColumnGroup(&buf, co, nil, "")
	if err != nil {
		t.Fatalf("flushColumnGroup(nil) err = %v", err)
	}
	if n != 0 {
		t.Errorf("flushColumnGroup(nil) = %d, want 0", n)
	}
	if buf.Len() != 0 {
		t.Errorf("expected no output, got %q", buf.String())
	}
}

// Compile-time sanity: assert errReader satisfies io.Reader.
var _ io.Reader = (*errReader)(nil)
