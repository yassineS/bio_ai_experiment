// Raw, column-preserving implementation of the upstream-parity intersect
// output modes (default intersection, -wa, -wb, -c, -v). Unlike the legacy
// typed bed.Record path, this echoes A's (and B's, for -wb) original input
// columns verbatim and re-encodes clipped coordinates per the source format,
// matching `bedtools intersect` byte-for-byte. It also supports BED/VCF/GFF/BAM
// inputs via readInRecords.
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
// stamps each B record with its in-chromosome order so the tree path can
// restore file order.
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

// intersectRaw runs the default / -wa / -wb / -c / -v output modes over the
// raw inRecord model. Returns the number of output lines written.
func intersectRaw(readerA, readerB io.Reader, writer io.Writer, opts IntersectOptions) (int, error) {
	bRecords, err := readInRecords(readerB)
	if err != nil {
		return 0, fmt.Errorf("error reading B intervals: %w", err)
	}
	// Index B by chromosome. Either path preserves upstream's B-file output
	// order: the linear scan walks the slice in file order, and the tree path
	// re-sorts each query's candidates back into it.
	finder := newFinder(bRecords, opts)

	aRecords, err := readInRecords(readerA)
	if err != nil {
		return 0, fmt.Errorf("error reading A intervals: %w", err)
	}

	bw := bufio.NewWriter(writer)
	count := 0
	for _, a := range aRecords {
		hits := finder.overlaps(a, opts)
		if err := emitRawOutput(bw, a, hits, opts, &count); err != nil {
			return count, err
		}
	}
	if err := bw.Flush(); err != nil {
		return count, fmt.Errorf("error flushing output: %w", err)
	}
	return count, nil
}

// emitRawOutput writes the output line(s) for one A record and its B hits under
// the default / -wa / -wb / -c / -v output modes, exactly as upstream `bedtools
// intersect` prints them, and bumps *count once per emitted line. It is shared
// by intersectRaw and intersectRawWithStats so the two stay byte-for-byte
// identical, and writes straight to the buffer (no per-record slice).
func emitRawOutput(bw *bufio.Writer, a *inRecord, hits []rawHit, opts IntersectOptions, count *int) error {
	emit := func(s string) error {
		if _, err := bw.WriteString(s); err != nil {
			return fmt.Errorf("error writing result: %w", err)
		}
		if err := bw.WriteByte('\n'); err != nil {
			return fmt.Errorf("error writing result: %w", err)
		}
		*count++
		return nil
	}
	switch {
	case opts.NoOverlap: // -v: emit A (verbatim) only when there are no hits
		if len(hits) == 0 {
			return emit(a.line)
		}
		return nil
	case opts.Count: // -c: A (verbatim) + the overlap count as a trailing column
		return emit(a.line + "\t" + strconv.Itoa(len(hits)))
	default:
		for _, h := range hits {
			var line string
			switch {
			case opts.WriteA && opts.WriteB:
				line = a.line + "\t" + h.b.line
			case opts.WriteB:
				// -wb alone: A clipped to the overlap, then full original B.
				line = a.clippedLine(h.start, h.end) + "\t" + h.b.line
			case opts.WriteA:
				line = a.line
			default:
				// Default: A clipped to the overlap span, columns verbatim.
				line = a.clippedLine(h.start, h.end)
			}
			if err := emit(line); err != nil {
				return err
			}
		}
		return nil
	}
}

// intersectRawWithStats is intersectRaw plus per-A hit accounting for the -S
// stats summary. It supports the same output modes (default / -wa / -wb / -c /
// -v) over the raw inRecord model.
func intersectRawWithStats(readerA, readerB io.Reader, writer io.Writer, opts IntersectOptions) (*Stats, error) {
	bRecords, err := readInRecords(readerB)
	if err != nil {
		return nil, fmt.Errorf("error reading B intervals: %w", err)
	}
	finder := newFinder(bRecords, opts)
	aRecords, err := readInRecords(readerA)
	if err != nil {
		return nil, fmt.Errorf("error reading A intervals: %w", err)
	}

	stats := &Stats{IntervalsB: len(bRecords)}
	bw := bufio.NewWriter(writer)
	discard := 0
	for _, a := range aRecords {
		stats.IntervalsA++
		hits := finder.overlaps(a, opts)
		if len(hits) > 0 {
			stats.IntervalsAHit++
			stats.Overlaps += len(hits)
		} else {
			stats.IntervalsAMiss++
		}
		if err := emitRawOutput(bw, a, hits, opts, &discard); err != nil {
			return stats, err
		}
	}
	if err := bw.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output: %w", err)
	}
	return stats, nil
}

// rawHit records a B record overlapping A, with the 0-based overlap span used
// for the default-mode clip.
type rawHit struct {
	b     *inRecord
	start int
	end   int
}

// rawOverlaps returns the B records overlapping A (in B position order),
// applying the -s strand filter and the -f/-F/-r fraction tests with
// split-aware block math when opts.Split is set. The overlap span carried back
// is the whole-record intersection (used to clip A in default mode); under
// -split the fraction tests use the non-redundant block overlap but the clip
// span is still the whole-span intersection, matching upstream.
func rawOverlaps(a *inRecord, bRecords []*inRecord, opts IntersectOptions) []rawHit {
	var hits []rawHit
	aLen := a.end - a.start
	for _, b := range bRecords {
		if a.chrom != b.chrom {
			continue
		}
		if opts.StrandSpec && !sameStrandMatch(a.strand, b.strand) {
			continue
		}
		overlapStart := max(a.start, b.start)
		overlapEnd := min(a.end, b.end)
		overlapLen := overlapEnd - overlapStart
		if overlapLen <= 0 {
			continue
		}
		if overlapLen < opts.MinOverlap {
			continue
		}

		fracOverlap := overlapLen
		lenA := aLen
		lenB := b.end - b.start
		if opts.Split {
			blockOverlap, aSum, bSum := splitBlockOverlapIn(a, b)
			if blockOverlap <= 0 {
				continue
			}
			fracOverlap = blockOverlap
			lenA = aSum
			lenB = bSum
		}
		if opts.FractionA > 0 && fraction(fracOverlap, lenA) < opts.FractionA {
			continue
		}
		if opts.FractionB > 0 && fraction(fracOverlap, lenB) < opts.FractionB {
			continue
		}
		hits = append(hits, rawHit{b: b, start: overlapStart, end: overlapEnd})
	}
	return hits
}

// splitBlockOverlapIn mirrors splitBlockOverlap (the typed-record version) but
// operates on inRecords, expanding BED12 / BAM block columns into absolute
// intervals and returning the non-redundant block overlap and the summed block
// lengths used by the fraction tests.
func splitBlockOverlapIn(a, b *inRecord) (overlap, aSum, bSum int) {
	aBlocks := inBlocks(a)
	bBlocks := inBlocks(b)
	var ovs []block
	for _, hb := range bBlocks {
		for _, kb := range aBlocks {
			s := max(kb.start, hb.start)
			e := min(kb.end, hb.end)
			if e > s {
				ovs = append(ovs, block{s, e})
			}
		}
	}
	return nonRedundantOverlap(ovs), blockSum(aBlocks), blockSum(bBlocks)
}

// inBlocks expands an inRecord into absolute block intervals. BED12 and BAM
// records carry blockCount / blockSizes / blockStarts columns (fields 9..11);
// any other record (or a malformed BED12) is treated as a single whole-span
// block.
func inBlocks(r *inRecord) []block {
	switch r.format {
	case fmtBED, fmtBAM:
		if blks, ok := bed12BlocksFromFields(r.fields, r.start); ok {
			return blks
		}
	}
	// Single whole-span block, with the zero-length [p,p] -> [p-1,p+1] expansion
	// so a zero-length record can still intersect under -split (matching
	// blocksOf in the join path and upstream's adjustZeroLength).
	s, e, _ := effectiveBounds(r.start, r.end)
	return []block{{s, e}}
}

// bed12BlocksFromFields expands a BED12 record's blocks into absolute
// [start,end) ranges from its raw fields, returning ok=false when the record is
// not parseable as BED12 (fewer than 12 columns or malformed block columns).
func bed12BlocksFromFields(fields []string, recStart int) ([]block, bool) {
	if len(fields) < 12 {
		return nil, false
	}
	blockCount, err := strconv.Atoi(fields[9])
	if err != nil || blockCount <= 0 {
		return nil, false
	}
	sizes := splitCSV(fields[10])
	starts := splitCSV(fields[11])
	if len(sizes) != blockCount || len(starts) != blockCount {
		return nil, false
	}
	blks := make([]block, 0, blockCount)
	for i := 0; i < blockCount; i++ {
		off, err := strconv.Atoi(starts[i])
		if err != nil {
			return nil, false
		}
		size, err := strconv.Atoi(sizes[i])
		if err != nil {
			return nil, false
		}
		s := recStart + off
		blks = append(blks, block{s, s + size})
	}
	return blks, true
}
