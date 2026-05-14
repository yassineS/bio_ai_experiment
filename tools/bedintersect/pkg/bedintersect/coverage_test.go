package bedintersect

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// errWriter is an io.Writer that always returns errFailingWrite after n bytes
// have been accepted. If n is 0, it fails on the first Write.
type errWriter struct {
	n int
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

// Compile-time check that errWriter satisfies io.Writer.
var _ io.Writer = (*errWriter)(nil)

// largeBedInput builds a BED text body with n non-overlapping records on
// chromosome chr. Used to overflow the bufio.Writer inside bed.Writer and
// trigger downstream Write errors.
func largeBedInput(chr string, n int) string {
	var sb strings.Builder
	for i := 0; i < n; i++ {
		fmt.Fprintf(&sb, "%s\t%d\t%d\n", chr, i*100, i*100+10)
	}
	return sb.String()
}

// --- Intersect: error paths ----------------------------------------------

// TestIntersectReadErrorB covers the "error reading B intervals" branch
// (bedintersect.go:40) by feeding a malformed B record.
func TestIntersectReadErrorB(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr1\tnotanumber\t200\n"
	var buf bytes.Buffer
	_, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1})
	if err == nil {
		t.Fatalf("expected error for malformed B record, got nil")
	}
	if !strings.Contains(err.Error(), "error reading B intervals") {
		t.Errorf("expected wrapped B-read error, got: %v", err)
	}
}

// TestIntersectReadErrorA covers the "error reading A intervals" branch
// (bedintersect.go:79) by feeding a malformed A record.
func TestIntersectReadErrorA(t *testing.T) {
	fileA := "chr1\tnotanumber\t200\n"
	fileB := "chr1\t100\t200\n"
	var buf bytes.Buffer
	_, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1})
	if err == nil {
		t.Fatalf("expected error for malformed A record, got nil")
	}
	if !strings.Contains(err.Error(), "error reading A intervals") {
		t.Errorf("expected wrapped A-read error, got: %v", err)
	}
}

// TestIntersectWriteErrorIntersection covers the "error writing result" path
// for the default intersection output by overflowing bufio with many records.
func TestIntersectWriteErrorIntersection(t *testing.T) {
	// Many A records that each overlap the single big B record, generating
	// enough output to overflow the 4 KiB bufio buffer inside bed.Writer.
	var a strings.Builder
	for i := 0; i < 1000; i++ {
		fmt.Fprintf(&a, "chr1\t%d\t%d\n", i*100, i*100+50)
	}
	fileB := "chr1\t0\t10000000\n"
	_, err := Intersect(strings.NewReader(a.String()), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1})
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
	if !strings.Contains(err.Error(), "error writing result") && !strings.Contains(err.Error(), "error flushing output") {
		t.Errorf("expected write/flush error, got: %v", err)
	}
}

// TestIntersectFlushError covers the "error flushing output" branch
// (bedintersect.go:179). One small overlap fits in bufio, so the error
// surfaces during Flush, not on Write.
func TestIntersectFlushError(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr1\t150\t250\n"
	_, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1})
	if err == nil {
		t.Fatalf("expected flush error, got nil")
	}
	if !strings.Contains(err.Error(), "error flushing output") && !strings.Contains(err.Error(), "error writing result") {
		t.Errorf("expected flushing-output error, got: %v", err)
	}
}

// TestIntersectWriteErrorNoOverlap covers the write-error branch in the
// NoOverlap output path (bedintersect.go:138).
func TestIntersectWriteErrorNoOverlap(t *testing.T) {
	fileA := largeBedInput("chr1", 1000)
	fileB := "chr2\t100\t200\n" // different chr, so every A has no overlap
	_, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1, NoOverlap: true})
	if err == nil {
		t.Fatalf("expected write error in NoOverlap mode, got nil")
	}
}

// TestIntersectWriteErrorCount covers the write-error branch in the Count
// output path (bedintersect.go:151).
func TestIntersectWriteErrorCount(t *testing.T) {
	fileA := largeBedInput("chr1", 1000)
	fileB := "chr1\t0\t100000000\n"
	_, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1, Count: true})
	if err == nil {
		t.Fatalf("expected write error in Count mode, got nil")
	}
}

// TestIntersectWriteErrorDistance covers the write-error branch in the
// Distance output path (bedintersect.go:94 and 114).
func TestIntersectWriteErrorDistance(t *testing.T) {
	fileA := largeBedInput("chr1", 1000)
	fileB := "chr1\t10000000\t10000010\n"
	_, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1, Distance: true})
	if err == nil {
		t.Fatalf("expected write error in Distance mode, got nil")
	}
}

// TestIntersectWriteErrorDistanceNoB covers the write-error branch in the
// Distance-no-B-on-chromosome path (bedintersect.go:114) — A records on a
// chromosome with no B intervals trigger the "-1" branch.
func TestIntersectWriteErrorDistanceNoB(t *testing.T) {
	fileA := largeBedInput("chr1", 1000)
	fileB := "chr2\t100\t200\n"
	_, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1, Distance: true})
	if err == nil {
		t.Fatalf("expected write error in Distance-no-B mode, got nil")
	}
}

// TestIntersectDistanceNoBOnChromosome covers the non-error half of the
// "no B intervals on this chromosome" path with Distance (writes "-1").
func TestIntersectDistanceNoBOnChromosome(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr2\t100\t200\n"
	var buf bytes.Buffer
	count, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, Distance: true})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}
	if count != 1 {
		t.Errorf("expected 1 record, got %d", count)
	}
	if !strings.Contains(buf.String(), "-1") {
		t.Errorf("expected -1 distance for A with no B on same chromosome, got %q", buf.String())
	}
}

// TestIntersectClosestNoBOnChromosome covers the Closest mode when there are
// no B intervals on A's chromosome — the closest==nil branch with no Distance
// flag, which should not write anything (count stays 0).
func TestIntersectClosestNoBOnChromosome(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr2\t100\t200\n"
	var buf bytes.Buffer
	count, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, Closest: true})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}
	if count != 0 {
		t.Errorf("expected 0 closest output when no B on chromosome, got %d", count)
	}
	if buf.Len() != 0 {
		t.Errorf("expected empty output, got %q", buf.String())
	}
}

// TestIntersectClosestWriteError covers the write-error branch in the Closest
// output path (bedintersect.go:101).
func TestIntersectClosestWriteError(t *testing.T) {
	// Each A finds a closest B (same chrom).
	fileA := largeBedInput("chr1", 1000)
	fileB := "chr1\t0\t100000000\n"
	_, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1, Closest: true})
	if err == nil {
		t.Fatalf("expected write error in Closest mode, got nil")
	}
}

// TestIntersectNoOverlapWithSomeOverlapping covers the suppression branch in
// NoOverlap mode: an A interval with overlaps must NOT be written. Some
// records have overlaps (skipped) and some don't (written), which exercises
// both halves of the if-len(overlaps)==0 branch.
func TestIntersectNoOverlapWithSomeOverlapping(t *testing.T) {
	fileA := "chr1\t100\t200\nchr1\t300\t400\nchr1\t500\t600\n"
	fileB := "chr1\t150\t250\nchr1\t550\t650\n" // overlaps with A#1 and A#3
	var buf bytes.Buffer
	count, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, NoOverlap: true})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}
	if count != 1 {
		t.Errorf("expected 1 non-overlapping record, got %d", count)
	}
	if !strings.Contains(buf.String(), "300\t400") {
		t.Errorf("expected the middle record (300-400) in output, got %q", buf.String())
	}
}

// TestIntersectUseTreeNoMatchingChromosome covers the UseTree code path where
// A's chromosome has no entry in chromTrees (no B intervals at all on that
// chrom). The lookup misses and overlaps stays nil.
func TestIntersectUseTreeNoMatchingChromosome(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr2\t100\t200\n"
	var buf bytes.Buffer
	count, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, UseTree: true})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}
	if count != 0 {
		t.Errorf("expected 0 overlaps, got %d", count)
	}
}

// --- IntersectWithStats: every uncovered branch ---------------------------

// TestIntersectWithStatsReadErrorB covers the malformed-B branch.
func TestIntersectWithStatsReadErrorB(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr1\tnot\t200\n"
	var buf bytes.Buffer
	_, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1})
	if err == nil || !strings.Contains(err.Error(), "error reading B intervals") {
		t.Fatalf("expected B-read error, got %v", err)
	}
}

// TestIntersectWithStatsReadErrorA covers the malformed-A branch.
func TestIntersectWithStatsReadErrorA(t *testing.T) {
	fileA := "chr1\tnot\t200\n"
	fileB := "chr1\t100\t200\n"
	var buf bytes.Buffer
	_, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1})
	if err == nil || !strings.Contains(err.Error(), "error reading A intervals") {
		t.Fatalf("expected A-read error, got %v", err)
	}
}

// TestIntersectWithStatsSortDifferentChromosomes exercises the chromosome
// branch of the sort comparator (bedintersect.go:327-329) where the two B
// records have different chromosome names and the comparator returns the
// chrom-based result.
func TestIntersectWithStatsSortDifferentChromosomes(t *testing.T) {
	// chr2 first, chr1 second — sort must reorder by chrom.
	fileA := "chr1\t100\t200\n"
	fileB := "chr2\t10\t20\nchr1\t150\t250\n"
	var buf bytes.Buffer
	stats, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1})
	if err != nil {
		t.Fatalf("IntersectWithStats failed: %v", err)
	}
	if stats.IntervalsB != 2 || stats.Overlaps != 1 {
		t.Errorf("unexpected stats: %+v", stats)
	}
}

// TestIntersectWithStatsUseTree exercises the UseTree build path and lookup
// in IntersectWithStats (bedintersect.go:341-345, 405-409).
func TestIntersectWithStatsUseTree(t *testing.T) {
	fileA := "chr1\t100\t200\nchr2\t100\t200\nchr3\t100\t200\n"
	fileB := "chr1\t150\t250\nchr2\t150\t250\n"
	var buf bytes.Buffer
	stats, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, UseTree: true})
	if err != nil {
		t.Fatalf("IntersectWithStats failed: %v", err)
	}
	if stats.IntervalsAHit != 2 {
		t.Errorf("expected 2 hits, got %d", stats.IntervalsAHit)
	}
	if stats.IntervalsAMiss != 1 {
		t.Errorf("expected 1 miss, got %d", stats.IntervalsAMiss)
	}
}

// TestIntersectWithStatsDistance covers the Distance branch (closest!=nil
// and the Distance write path) in IntersectWithStats.
func TestIntersectWithStatsDistance(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr1\t300\t400\n"
	var buf bytes.Buffer
	stats, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, Distance: true})
	if err != nil {
		t.Fatalf("IntersectWithStats failed: %v", err)
	}
	if stats.IntervalsAHit != 1 {
		t.Errorf("expected 1 hit, got %d", stats.IntervalsAHit)
	}
	if !strings.Contains(buf.String(), "100\t200\t100") {
		t.Errorf("expected distance 100 in output, got %q", buf.String())
	}
}

// TestIntersectWithStatsDistanceNoBOnChromosome covers the closest==nil with
// Distance branch (writes "-1").
func TestIntersectWithStatsDistanceNoBOnChromosome(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr2\t100\t200\n"
	var buf bytes.Buffer
	stats, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, Distance: true})
	if err != nil {
		t.Fatalf("IntersectWithStats failed: %v", err)
	}
	if stats.IntervalsAMiss != 1 {
		t.Errorf("expected 1 miss, got %d", stats.IntervalsAMiss)
	}
	if !strings.Contains(buf.String(), "-1") {
		t.Errorf("expected -1 distance in output, got %q", buf.String())
	}
}

// TestIntersectWithStatsClosest covers the Closest mode write path in
// IntersectWithStats.
func TestIntersectWithStatsClosest(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr1\t300\t400\n"
	var buf bytes.Buffer
	stats, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, Closest: true})
	if err != nil {
		t.Fatalf("IntersectWithStats failed: %v", err)
	}
	if stats.IntervalsAHit != 1 {
		t.Errorf("expected 1 hit, got %d", stats.IntervalsAHit)
	}
	if !strings.Contains(buf.String(), "300\t400") {
		t.Errorf("expected B record in output, got %q", buf.String())
	}
}

// TestIntersectWithStatsClosestNoBOnChromosome covers the closest==nil with
// no Distance flag set — increments IntervalsAMiss without writing.
func TestIntersectWithStatsClosestNoBOnChromosome(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr2\t100\t200\n"
	var buf bytes.Buffer
	stats, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, Closest: true})
	if err != nil {
		t.Fatalf("IntersectWithStats failed: %v", err)
	}
	if stats.IntervalsAMiss != 1 {
		t.Errorf("expected 1 miss, got %d", stats.IntervalsAMiss)
	}
	if buf.Len() != 0 {
		t.Errorf("expected empty output, got %q", buf.String())
	}
}

// TestIntersectWithStatsNoOverlap covers the NoOverlap (-v) write path.
func TestIntersectWithStatsNoOverlap(t *testing.T) {
	fileA := "chr1\t100\t200\nchr1\t300\t400\n"
	fileB := "chr1\t150\t250\n"
	var buf bytes.Buffer
	stats, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, NoOverlap: true})
	if err != nil {
		t.Fatalf("IntersectWithStats failed: %v", err)
	}
	if stats.IntervalsAHit != 1 || stats.IntervalsAMiss != 1 {
		t.Errorf("unexpected stats: %+v", stats)
	}
	if !strings.Contains(buf.String(), "300\t400") || strings.Contains(buf.String(), "100\t200") {
		t.Errorf("expected only the non-overlapping record, got %q", buf.String())
	}
}

// TestIntersectWithStatsCount covers the Count write path.
func TestIntersectWithStatsCount(t *testing.T) {
	fileA := "chr1\t100\t500\n"
	fileB := "chr1\t150\t200\nchr1\t250\t300\n"
	var buf bytes.Buffer
	stats, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, Count: true})
	if err != nil {
		t.Fatalf("IntersectWithStats failed: %v", err)
	}
	if stats.Overlaps != 2 {
		t.Errorf("expected 2 overlaps, got %d", stats.Overlaps)
	}
	if !strings.Contains(buf.String(), "\t2\n") {
		t.Errorf("expected count=2 in output, got %q", buf.String())
	}
}

// TestIntersectWithStatsWriteA covers the WriteA branch in IntersectWithStats.
func TestIntersectWithStatsWriteA(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr1\t150\t250\n"
	var buf bytes.Buffer
	_, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, WriteA: true})
	if err != nil {
		t.Fatalf("IntersectWithStats failed: %v", err)
	}
	if strings.TrimSpace(buf.String()) != "chr1\t100\t200" {
		t.Errorf("expected A record, got %q", buf.String())
	}
}

// TestIntersectWithStatsWriteB covers the WriteB branch in IntersectWithStats.
func TestIntersectWithStatsWriteB(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr1\t150\t250\n"
	var buf bytes.Buffer
	_, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, WriteB: true})
	if err != nil {
		t.Fatalf("IntersectWithStats failed: %v", err)
	}
	if strings.TrimSpace(buf.String()) != "chr1\t150\t250" {
		t.Errorf("expected B record, got %q", buf.String())
	}
}

// TestIntersectWithStatsWriteErrorDefault covers the write-error branch in
// the default intersection write path.
func TestIntersectWithStatsWriteErrorDefault(t *testing.T) {
	var a strings.Builder
	for i := 0; i < 1000; i++ {
		fmt.Fprintf(&a, "chr1\t%d\t%d\n", i*100, i*100+50)
	}
	fileB := "chr1\t0\t100000000\n"
	_, err := IntersectWithStats(strings.NewReader(a.String()), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1})
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
}

// TestIntersectWithStatsFlushError covers the flush-error branch.
func TestIntersectWithStatsFlushError(t *testing.T) {
	fileA := "chr1\t100\t200\n"
	fileB := "chr1\t150\t250\n"
	_, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1})
	if err == nil {
		t.Fatalf("expected flush/write error, got nil")
	}
}

// TestIntersectWithStatsWriteErrorNoOverlap covers the write-error branch in
// IntersectWithStats NoOverlap mode.
func TestIntersectWithStatsWriteErrorNoOverlap(t *testing.T) {
	fileA := largeBedInput("chr1", 1000)
	fileB := "chr2\t100\t200\n"
	_, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1, NoOverlap: true})
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
}

// TestIntersectWithStatsWriteErrorCount covers the write-error branch in
// IntersectWithStats Count mode.
func TestIntersectWithStatsWriteErrorCount(t *testing.T) {
	fileA := largeBedInput("chr1", 1000)
	fileB := "chr1\t0\t100000000\n"
	_, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1, Count: true})
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
}

// TestIntersectWithStatsWriteErrorDistance covers the write-error branch in
// IntersectWithStats Distance mode (hit path).
func TestIntersectWithStatsWriteErrorDistance(t *testing.T) {
	fileA := largeBedInput("chr1", 1000)
	fileB := "chr1\t10000000\t10000010\n"
	_, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1, Distance: true})
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
}

// TestIntersectWithStatsWriteErrorDistanceNoB covers the write-error branch
// in IntersectWithStats Distance mode when there's no B on the chromosome
// (the "-1" branch).
func TestIntersectWithStatsWriteErrorDistanceNoB(t *testing.T) {
	fileA := largeBedInput("chr1", 1000)
	fileB := "chr2\t100\t200\n"
	_, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1, Distance: true})
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
}

// TestIntersectWithStatsWriteErrorClosest covers the write-error branch in
// IntersectWithStats Closest mode.
func TestIntersectWithStatsWriteErrorClosest(t *testing.T) {
	fileA := largeBedInput("chr1", 1000)
	fileB := "chr1\t0\t100000000\n"
	_, err := IntersectWithStats(strings.NewReader(fileA), strings.NewReader(fileB), &errWriter{}, IntersectOptions{MinOverlap: 1, Closest: true})
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
}

// --- findClosest: uncovered branches --------------------------------------

// TestFindClosestOverlapping covers the "Overlapping" branch where the A and
// B intervals overlap, producing dist = 0.
func TestFindClosestOverlapping(t *testing.T) {
	a := &bed.Record{Chrom: "chr1", ChromStart: 100, ChromEnd: 200}
	bs := []*bed.Record{
		{Chrom: "chr1", ChromStart: 150, ChromEnd: 250},
	}
	closest, dist := findClosest(a, bs, IntersectOptions{})
	if closest == nil || dist != 0 {
		t.Errorf("expected overlapping with dist=0, got closest=%v dist=%d", closest, dist)
	}
}

// TestFindClosestStrandMismatch covers the strand-mismatch skip branch in
// findClosest. With StrandSpec set and A's strand != B's strand, the B record
// is skipped entirely; the remaining (non-matching-chrom) B is also skipped,
// leaving the result nil/-1.
func TestFindClosestStrandMismatch(t *testing.T) {
	a := &bed.Record{Chrom: "chr1", ChromStart: 100, ChromEnd: 200, Strand: "+"}
	bs := []*bed.Record{
		{Chrom: "chr1", ChromStart: 300, ChromEnd: 400, Strand: "-"},
	}
	closest, dist := findClosest(a, bs, IntersectOptions{StrandSpec: true})
	if closest != nil {
		t.Errorf("expected no closest with strand mismatch, got %+v (dist=%d)", closest, dist)
	}
}

// TestFindClosestStrandMatch covers the StrandSpec branch with a matching
// strand — execution falls through the strand check and computes distance.
func TestFindClosestStrandMatch(t *testing.T) {
	a := &bed.Record{Chrom: "chr1", ChromStart: 100, ChromEnd: 200, Strand: "+"}
	bs := []*bed.Record{
		{Chrom: "chr1", ChromStart: 300, ChromEnd: 400, Strand: "+"},
	}
	closest, dist := findClosest(a, bs, IntersectOptions{StrandSpec: true})
	if closest == nil || dist != 100 {
		t.Errorf("expected closest with dist=100, got %+v (dist=%d)", closest, dist)
	}
}

// --- findOverlaps: uncovered branches -------------------------------------

// TestFindOverlapsStrandMismatch covers the strand-mismatch skip in
// findOverlaps.
func TestFindOverlapsStrandMismatch(t *testing.T) {
	fileA := "chr1\t100\t200\tnameA\t0\t+\n"
	fileB := "chr1\t150\t250\tnameB\t0\t-\n"
	var buf bytes.Buffer
	count, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, StrandSpec: true})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}
	if count != 0 {
		t.Errorf("expected 0 overlaps with strand mismatch, got %d", count)
	}
}

// TestFindOverlapsStrandMatch covers the StrandSpec branch with matching
// strands — should produce an overlap.
func TestFindOverlapsStrandMatch(t *testing.T) {
	fileA := "chr1\t100\t200\tnameA\t0\t+\n"
	fileB := "chr1\t150\t250\tnameB\t0\t+\n"
	var buf bytes.Buffer
	count, err := Intersect(strings.NewReader(fileA), strings.NewReader(fileB), &buf, IntersectOptions{MinOverlap: 1, StrandSpec: true})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}
	if count != 1 {
		t.Errorf("expected 1 overlap with matching strand, got %d", count)
	}
}

// --- IntervalTree: uncovered branches -------------------------------------

// TestBuildTreeLeftSubtreeMaxLarger covers the "node.Max = node.Left.Max"
// branch in buildTree by giving the left subtree an interval whose ChromEnd
// is larger than the root's. The Max propagated up should be the left's.
func TestBuildTreeLeftSubtreeMaxLarger(t *testing.T) {
	// Sorted by ChromStart so buildTree builds a balanced tree:
	//   mid = (0+2)/2 = 1 -> root = intervals[1], with ChromEnd = 150.
	//   Left = intervals[0] (ChromStart 50, ChromEnd 1000) — large end.
	//   Right = intervals[2].
	// The root must update Max = max(150, 1000, 250) = 1000 via Left.
	intervals := []*bed.Record{
		{Chrom: "chr1", ChromStart: 50, ChromEnd: 1000},
		{Chrom: "chr1", ChromStart: 100, ChromEnd: 150},
		{Chrom: "chr1", ChromStart: 200, ChromEnd: 250},
	}
	tree := NewIntervalTree(intervals)
	if tree.Root.Max != 1000 {
		t.Errorf("expected root.Max=1000 (from left subtree), got %d", tree.Root.Max)
	}
	// Also exercise a query that prunes via the (query.ChromStart >= node.Max)
	// check (bedintersect.go:77) at the root level.
	res := tree.Query(&bed.Record{Chrom: "chr1", ChromStart: 5000, ChromEnd: 6000})
	if len(res) != 0 {
		t.Errorf("expected 0 results from out-of-range query, got %d", len(res))
	}
}

// TestQueryNodeNilGuard covers the nil-node guard in queryNode by invoking
// Query on a tree whose Root has nil children. While buildTree wouldn't
// normally pass a nil node (the callers check), the explicit guard is reached
// when the tree has a single node and the recursive descent reaches a leaf.
// We also assert tree.Query handles a Root==nil tree (empty intervals).
func TestQueryNodeNilGuard(t *testing.T) {
	// Empty tree -> Root nil -> Query returns nil before entering queryNode.
	empty := NewIntervalTree(nil)
	if got := empty.Query(&bed.Record{Chrom: "chr1", ChromStart: 1, ChromEnd: 2}); got != nil {
		t.Errorf("empty Query: expected nil, got %v", got)
	}
	// Single-node tree: queryNode will be invoked with a non-nil node, then
	// recurse into nil Left/Right and hit the guard.
	tree := NewIntervalTree([]*bed.Record{{Chrom: "chr1", ChromStart: 100, ChromEnd: 200}})
	res := tree.Query(&bed.Record{Chrom: "chr1", ChromStart: 150, ChromEnd: 175})
	if len(res) != 1 {
		t.Errorf("expected 1 result, got %d", len(res))
	}

	// Note: the explicit nil-guard branch inside the internal queryNode
	// recursion is now covered by the shared bed.IntervalTree tests in
	// pkg/bioformats/bed/. We can no longer call it directly from this
	// package because the method moved out of bedintersect when the tree was
	// lifted into pkg/bioformats/bed.
}

// TestFindOverlapsChromMismatch directly calls findOverlaps with a B slice
// that contains a record on a different chromosome to A. Normal callers feed
// the per-chromosome chromIndex, so this defensive skip branch
// (bedintersect.go:198-200) is otherwise unreachable.
func TestFindOverlapsChromMismatch(t *testing.T) {
	a := &bed.Record{Chrom: "chr1", ChromStart: 100, ChromEnd: 200}
	bs := []*bed.Record{
		{Chrom: "chr2", ChromStart: 100, ChromEnd: 300}, // skipped: wrong chrom
		{Chrom: "chr1", ChromStart: 150, ChromEnd: 250}, // overlaps
	}
	overlaps := findOverlaps(a, bs, IntersectOptions{MinOverlap: 1})
	if len(overlaps) != 1 {
		t.Fatalf("expected 1 overlap (chr1 only), got %d", len(overlaps))
	}
	if overlaps[0].B.Chrom != "chr1" {
		t.Errorf("expected overlap with chr1 B, got %s", overlaps[0].B.Chrom)
	}
}

// TestFindClosestChromMismatch directly calls findClosest with a B slice
// containing a wrong-chromosome record. The defensive chrom-mismatch skip
// (bedintersect.go:260-262) is otherwise unreachable through normal callers.
func TestFindClosestChromMismatch(t *testing.T) {
	a := &bed.Record{Chrom: "chr1", ChromStart: 100, ChromEnd: 200}
	bs := []*bed.Record{
		{Chrom: "chr2", ChromStart: 100, ChromEnd: 300}, // skipped: wrong chrom
		{Chrom: "chr1", ChromStart: 300, ChromEnd: 400}, // closest
	}
	closest, dist := findClosest(a, bs, IntersectOptions{})
	if closest == nil || closest.Chrom != "chr1" || dist != 100 {
		t.Errorf("expected closest to be chr1 at dist=100, got %+v dist=%d", closest, dist)
	}
}
