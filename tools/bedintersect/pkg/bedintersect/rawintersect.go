// Raw, column-preserving implementation of every upstream-parity intersect
// output mode (default intersection, -wa, -wb, -wo, -wao, -loj, -c, -C, -u, -v).
// Unlike the legacy typed bed.Record path it echoes A's (and B's, where the mode
// prints it) original input columns verbatim and re-encodes clipped coordinates
// per the source format, matching `bedtools intersect` byte-for-byte. It also
// supports BED/VCF/GFF/BAM/CRAM inputs and multiple B files via readAllB.
package bedintersect

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
)

// inFinder returns the B records overlapping A. It abstracts over the two B
// lookup strategies: a per-chromosome linear scan, and the augmented
// interval-tree index used under -sortedtree (opts.UseTree) for large B files.
type inFinder interface {
	overlaps(a *inRecord, opts IntersectOptions) []rawHit
}

// linearFinder scans the chromosome's B slice in file order. O(A*B) per
// chromosome but allocation-free in setup; the default for modest inputs.
type linearFinder struct {
	byChrom map[string][]*inRecord
}

func (f linearFinder) overlaps(a *inRecord, opts IntersectOptions) []rawHit {
	return rawOverlaps(a, f.byChrom[a.chrom], opts)
}

// treeFinder answers overlap queries with one augmented interval tree per
// chromosome in O(log n + k). Tree queries return hits out of file order, so
// the candidate slice is re-sorted by each record's original position before
// the fraction/strand filters run, preserving upstream's B-file output order.
type treeFinder struct {
	trees map[string]*inIntervalTree
}

func newTreeFinder(byChrom map[string][]*inRecord) treeFinder {
	trees := make(map[string]*inIntervalTree, len(byChrom))
	for chrom, recs := range byChrom {
		trees[chrom] = newInIntervalTree(recs)
	}
	return treeFinder{trees: trees}
}

func (f treeFinder) overlaps(a *inRecord, opts IntersectOptions) []rawHit {
	tree, ok := f.trees[a.chrom]
	if !ok {
		return nil
	}
	cands := tree.query(a)
	sort.Slice(cands, func(i, j int) bool { return cands[i].order < cands[j].order })
	return rawOverlaps(a, cands, opts)
}

// newFinder picks the tree- or linear-scan B index based on opts.UseTree, and
// stamps each B record with its in-chromosome order so the tree path can restore
// file order. The records are kept in their readAllB order — grouped by B file,
// then by position within file — so multi-database hits print exactly as
// upstream does.
func newFinder(bRecords []*inRecord, opts IntersectOptions) inFinder {
	byChrom := make(map[string][]*inRecord)
	for _, b := range bRecords {
		b.order = len(byChrom[b.chrom])
		byChrom[b.chrom] = append(byChrom[b.chrom], b)
	}
	if opts.UseTree {
		return newTreeFinder(byChrom)
	}
	return linearFinder{byChrom: byChrom}
}

// inIntervalNode is a node in the augmented interval tree over inRecords. max is
// the largest end across the subtree, used to prune queries.
type inIntervalNode struct {
	rec         *inRecord
	max         int
	left, right *inIntervalNode
}

// inIntervalTree is a balanced augmented BST over a single chromosome's B
// records, mirroring pkg/htsgo/bed.IntervalTree but keyed on inRecord so the
// raw column-preserving path can use it for BED/VCF/GFF/BAM inputs alike.
type inIntervalTree struct {
	root *inIntervalNode
}

// newInIntervalTree builds a balanced tree from one chromosome's records. The
// slice is sorted by start so the median split yields a balanced tree.
func newInIntervalTree(recs []*inRecord) *inIntervalTree {
	if len(recs) == 0 {
		return &inIntervalTree{}
	}
	sorted := append([]*inRecord(nil), recs...)
	sort.Slice(sorted, func(i, j int) bool { return sorted[i].start < sorted[j].start })
	return &inIntervalTree{root: buildInTree(sorted, 0, len(sorted)-1)}
}

func buildInTree(recs []*inRecord, lo, hi int) *inIntervalNode {
	if lo > hi {
		return nil
	}
	mid := (lo + hi) / 2
	node := &inIntervalNode{rec: recs[mid], max: recs[mid].end}
	node.left = buildInTree(recs, lo, mid-1)
	node.right = buildInTree(recs, mid+1, hi)
	if node.left != nil && node.left.max > node.max {
		node.max = node.left.max
	}
	if node.right != nil && node.right.max > node.max {
		node.max = node.right.max
	}
	return node
}

// query returns every record whose span overlaps a's span (half-open). The
// chromosome is implied by the tree, so it is not re-compared here.
func (t *inIntervalTree) query(a *inRecord) []*inRecord {
	if t.root == nil {
		return nil
	}
	var out []*inRecord
	queryInNode(t.root, a, &out)
	return out
}

func queryInNode(node *inIntervalNode, a *inRecord, out *[]*inRecord) {
	if node == nil || a.start >= node.max {
		return
	}
	queryInNode(node.left, a, out)
	if a.start < node.rec.end && a.end > node.rec.start {
		*out = append(*out, node.rec)
	}
	if a.end > node.rec.start {
		queryInNode(node.right, a, out)
	}
}

// intersectRaw runs every upstream-parity output mode over the raw inRecord
// model, across one or more B files. Returns the number of output lines written.
func intersectRaw(readerA io.Reader, readersB []io.Reader, writer io.Writer, opts IntersectOptions) (int, error) {
	bRecords, dbType, dbFields, err := readAllB(readersB)
	if err != nil {
		return 0, err
	}
	finder := newFinder(bRecords, opts)

	aRecords, hdr, err := readInRecordsWithHeader(readerA, opts.Header)
	if err != nil {
		return 0, fmt.Errorf("error reading A intervals: %w", err)
	}
	if err := checkSorted(opts, aRecords, bRecords); err != nil {
		return 0, err
	}
	emitNameWarning(opts, aRecords, bRecords)

	bw := bufio.NewWriter(writer)
	if hdr != "" {
		if _, err := bw.WriteString(hdr); err != nil {
			return 0, fmt.Errorf("error writing header: %w", err)
		}
	}
	em := &emitter{bw: bw, dbType: dbType, dbFields: dbFields, opts: opts}
	for _, a := range aRecords {
		hits := finder.overlaps(a, opts)
		if !opts.SortOut {
			orderHitsByBin(hits)
		}
		if err := em.emit(a, hits); err != nil {
			return em.count, err
		}
	}
	if err := bw.Flush(); err != nil {
		return em.count, fmt.Errorf("error flushing output: %w", err)
	}
	return em.count, nil
}

// intersectRawWithStats is intersectRaw plus per-A hit accounting for the
// bedintersect-only --stats summary. It supports the same output modes over the
// raw inRecord model and one or more B files.
func intersectRawWithStats(readerA io.Reader, readersB []io.Reader, writer io.Writer, opts IntersectOptions) (*Stats, error) {
	bRecords, dbType, dbFields, err := readAllB(readersB)
	if err != nil {
		return nil, err
	}
	finder := newFinder(bRecords, opts)
	aRecords, hdr, err := readInRecordsWithHeader(readerA, opts.Header)
	if err != nil {
		return nil, fmt.Errorf("error reading A intervals: %w", err)
	}

	if err := checkSorted(opts, aRecords, bRecords); err != nil {
		return nil, err
	}
	emitNameWarning(opts, aRecords, bRecords)

	stats := &Stats{IntervalsB: len(bRecords)}
	bw := bufio.NewWriter(writer)
	if hdr != "" {
		if _, err := bw.WriteString(hdr); err != nil {
			return stats, fmt.Errorf("error writing header: %w", err)
		}
	}
	em := &emitter{bw: bw, dbType: dbType, dbFields: dbFields, opts: opts}
	for _, a := range aRecords {
		stats.IntervalsA++
		hits := finder.overlaps(a, opts)
		if !opts.SortOut {
			orderHitsByBin(hits)
		}
		if len(hits) > 0 {
			stats.IntervalsAHit++
			stats.Overlaps += len(hits)
		} else {
			stats.IntervalsAMiss++
		}
		if err := em.emit(a, hits); err != nil {
			return stats, err
		}
	}
	if err := bw.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output: %w", err)
	}
	return stats, nil
}

// emitter renders the output line(s) for each A record and its B hits, holding
// the buffered writer, the B-file classification (for null-B placeholders) and
// the resolved options. It is shared by intersectRaw and the stats path so the
// two stay byte-for-byte identical.
type emitter struct {
	bw       *bufio.Writer
	dbType   dbRecordType
	dbFields int
	opts     IntersectOptions
	count    int
}

func (e *emitter) line(s string) error {
	if _, err := e.bw.WriteString(s); err != nil {
		return fmt.Errorf("error writing result: %w", err)
	}
	if err := e.bw.WriteByte('\n'); err != nil {
		return fmt.Errorf("error writing result: %w", err)
	}
	e.count++
	return nil
}

// emit writes the output for one A record under the active output mode.
func (e *emitter) emit(a *inRecord, hits []rawHit) error {
	opts := e.opts
	switch {
	case opts.CountEach:
		return e.emitCountEach(a, hits)
	case opts.Count:
		return e.line(a.line + "\t" + strconv.Itoa(len(hits)))
	case opts.NoOverlap:
		if len(hits) == 0 {
			return e.line(a.line)
		}
		return nil
	case opts.Unique:
		if len(hits) > 0 {
			return e.line(a.line)
		}
		return nil
	case opts.usesJoinMode():
		return e.emitJoin(a, hits)
	default:
		// Default / -wa / -wb (without -wa).
		for _, h := range hits {
			var s string
			switch {
			case opts.WriteB:
				// -wb alone: A clipped to the overlap, then the DB-id column (if
				// multiple B files) and the full original B.
				s = a.clippedLine(h.start, h.end) + "\t" + e.bWithLabel(h.b)
			case opts.WriteA:
				s = a.line
			default:
				s = a.clippedLine(h.start, h.end)
			}
			if err := e.line(s); err != nil {
				return err
			}
		}
		return nil
	}
}

// emitCountEach implements -C: one line per B file reporting that file's overlap
// count with A (0 included), in file order. With a single B file the DB-id
// column is omitted, matching upstream.
func (e *emitter) emitCountEach(a *inRecord, hits []rawHit) error {
	opts := e.opts
	nFiles := len(opts.FilePaths)
	if nFiles < 1 {
		nFiles = 1
	}
	counts := make([]int, nFiles)
	for _, h := range hits {
		if h.b.dbID < nFiles {
			counts[h.b.dbID]++
		}
	}
	if !opts.multiDB() {
		return e.line(a.line + "\t" + strconv.Itoa(counts[0]))
	}
	for i := 0; i < nFiles; i++ {
		if err := e.line(a.line + "\t" + opts.dbLabel(i) + "\t" + strconv.Itoa(counts[i])); err != nil {
			return err
		}
	}
	return nil
}

// emitJoin implements the -loj / -wo / -wao / (-wa -wb) output modes, echoing A
// and B columns verbatim. Hits are grouped by B file (or sorted by position
// under -sortout) and each is prefixed with the DB-id column when multiple B
// files are present.
func (e *emitter) emitJoin(a *inRecord, hits []rawHit) error {
	opts := e.opts
	if len(hits) == 0 {
		nullB := e.bWithNullLabel()
		switch {
		case opts.WriteAllOverlap:
			return e.line(a.line + "\t" + nullB + "\t0")
		case opts.LeftJoin:
			return e.line(a.line + "\t" + nullB)
		}
		return nil
	}
	if opts.SortOut {
		hits = append([]rawHit(nil), hits...)
		sort.SliceStable(hits, func(i, j int) bool {
			if hits[i].b.chrom != hits[j].b.chrom {
				return hits[i].b.chrom < hits[j].b.chrom
			}
			if hits[i].b.start != hits[j].b.start {
				return hits[i].b.start < hits[j].b.start
			}
			return hits[i].b.end < hits[j].b.end
		})
	}
	for _, h := range hits {
		s := a.line + "\t" + e.bWithLabel(h.b)
		if opts.WriteOverlap || opts.WriteAllOverlap {
			s += "\t" + strconv.Itoa(h.overlapBases)
		}
		if err := e.line(s); err != nil {
			return err
		}
	}
	return nil
}

// bWithLabel renders a B record's columns, prefixing the DB-id column when more
// than one B file is present (matching upstream's fileId/name prefix).
func (e *emitter) bWithLabel(b *inRecord) string {
	if e.opts.multiDB() {
		return e.opts.dbLabel(b.dbID) + "\t" + b.line
	}
	return b.line
}

// bWithNullLabel renders the null-B placeholder, prefixed with the DB-id column
// when multiple B files are present. Upstream prints the null record once for an
// A with no hits in any database, using the FIRST B file's classification for
// the placeholder shape and "." for the DB-id column.
func (e *emitter) bWithNullLabel() string {
	null := nullDBString(e.dbType, e.dbFields)
	if e.opts.multiDB() {
		return ".\t" + null
	}
	return null
}

// rawHit records a B record overlapping A, with the 0-based overlap span (used
// for the default-mode clip) and the overlapping base count (used by -wo/-wao).
type rawHit struct {
	b            *inRecord
	start        int
	end          int
	overlapBases int
}

// rawOverlaps returns the B records overlapping A (in B order), applying the
// strand filter and the -f/-F/-r/-e fraction tests with split-aware block math
// when opts.Split is set. The overlap span carried back is the whole-record
// intersection (used to clip A in default mode); the overlapping base count is
// the split-aware non-redundant count under -split, else the whole-span overlap,
// with the zero-length correction undone so -wo/-wao report upstream's value.
func rawOverlaps(a *inRecord, bRecords []*inRecord, opts IntersectOptions) []rawHit {
	if a.unmapped {
		return nil // an unmapped alignment never overlaps anything
	}
	if opts.Split {
		return rawOverlapsSplit(a, bRecords, opts)
	}
	var hits []rawHit
	aLen := a.end - a.start
	aStart, aEnd, aZero := effectiveBounds(a.start, a.end)
	for _, b := range bRecords {
		if b.unmapped || a.chrom != b.chrom {
			continue
		}
		if !strandOK(opts, a.strand, b.strand) {
			continue
		}
		bStart, bEnd, bZero := effectiveBounds(b.start, b.end)
		overlapStart := max(aStart, bStart)
		overlapEnd := min(aEnd, bEnd)
		overlapLen := overlapEnd - overlapStart
		if overlapLen <= 0 {
			continue
		}
		if overlapLen < opts.MinOverlap {
			continue
		}
		// The clip span printed in default mode uses the ORIGINAL (unexpanded)
		// coordinates: a zero-length A inside B clips to A's own [p,p] span, not
		// the [p-1,p+1] detection window. Mirrors upstream RecordOutputMgr.
		clipStart := max(a.start, b.start)
		clipEnd := min(a.end, b.end)
		// Undo the 1bp zero-length expansion when reporting overlap bases
		// (upstream RecordOutputMgr::reportOverlapDetail does maxStart++/minEnd--).
		overlapBases := overlapLen
		if aZero || bZero {
			overlapBases = (overlapEnd - 1) - (overlapStart + 1)
			if overlapBases < 0 {
				overlapBases = 0
			}
		}
		if !fractionOK(overlapBases, aLen, b.end-b.start, opts) {
			continue
		}
		hits = append(hits, rawHit{
			b:            b,
			start:        clipStart,
			end:          clipEnd,
			overlapBases: overlapBases,
		})
	}
	return hits
}

// rawOverlapsSplit is rawOverlaps with -split block math, mirroring
// BlockMgr::findBlockedOverlaps. A B record is a candidate if any of its blocks
// overlaps any A block; the per-hit overlap count is the total block overlap for
// that B. The -f/-F/-r fraction tests are applied ONCE across ALL hits combined
// (non-redundant overlap over A's block-sum and the summed block lengths of the
// hit B records), so either every hit for this A passes or all are dropped.
func rawOverlapsSplit(a *inRecord, bRecords []*inRecord, opts IntersectOptions) []rawHit {
	aBlocks := blocksOf(a, true)
	aBlockSum := blockSum(aBlocks)
	var hits []rawHit
	var allOverlaps []block
	hitBlockSum := 0
	for _, b := range bRecords {
		if b.unmapped || a.chrom != b.chrom {
			continue
		}
		if !strandOK(opts, a.strand, b.strand) {
			continue
		}
		bBlocks := blocksOf(b, true)
		overlapBases := 0
		var clipStart, clipEnd int
		first := true
		for _, hb := range bBlocks {
			for _, kb := range aBlocks {
				s := max(kb.start, hb.start)
				e := min(kb.end, hb.end)
				if e > s {
					overlapBases += e - s
					allOverlaps = append(allOverlaps, block{s, e})
					if first || s < clipStart {
						clipStart = s
					}
					if first || e > clipEnd {
						clipEnd = e
					}
					first = false
				}
			}
		}
		if overlapBases <= 0 {
			continue
		}
		hitBlockSum += blockSum(bBlocks)
		hits = append(hits, rawHit{
			b:            b,
			start:        clipStart,
			end:          clipEnd,
			overlapBases: overlapBases,
		})
	}
	// Apply the fraction tests once across the combined overlap, exactly as
	// BlockMgr::findBlockedOverlaps does with its totalUniqueOverlap.
	if len(hits) > 0 && (opts.FractionA > 0 || opts.FractionB > 0 || opts.Reciprocal) {
		uniq := nonRedundantOverlap(allOverlaps)
		if !fractionOK(uniq, aBlockSum, hitBlockSum, opts) {
			return nil
		}
	}
	return hits
}

// fractionOK applies the -f/-F/-r/-e fraction tests for one overlap, mirroring
// Record::sameChromIntersects. By default both the -f (fraction of A) and -F
// (fraction of B) thresholds must hold; -e requires only one of them; -r forces
// the B threshold up to the -f value (reciprocal). With no fraction requested at
// all, any positive overlap passes.
func fractionOK(overlap, lenA, lenB int, opts IntersectOptions) bool {
	fA := opts.FractionA
	fB := opts.FractionB
	if opts.Reciprocal && fB <= 0 {
		// -r without an explicit -F mirrors -F equal to -f.
		fB = fA
	}
	if fA <= 0 && fB <= 0 {
		return true
	}
	aPass := fA <= 0 || fraction(overlap, lenA) >= fA
	bPass := fB <= 0 || fraction(overlap, lenB) >= fB
	if opts.EitherFraction && fA > 0 && fB > 0 {
		// -e: either threshold suffices (only meaningful when both are set).
		return aPass || bPass
	}
	return aPass && bPass
}
