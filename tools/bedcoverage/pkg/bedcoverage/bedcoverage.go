// Package bedcoverage computes per-A coverage statistics from a B set of
// intervals, mirroring upstream `bedtools coverage`.
//
// For each record in A the tool reports, in order:
//
//   - count: number of B features that overlap A
//   - bp covered: total number of bases in A covered by at least one B feature
//   - length of A
//   - fraction (bp covered / length of A)
//
// Modes:
//   - default: append the four numbers to A's existing columns
//   - -counts: append only the count
//   - -d: emit one line per base in A: A's columns, 1-based offset within A,
//     depth at that base
//   - -hist: append per-depth-bucket histogram lines for each A, plus a
//     final "all" summary line aggregated across all A records
//   - numeric ops (-mean / -median / -min / -max / -sum): collapse the
//     per-base-depth vector with the requested op and append the single number
//
// Optional filters:
//   - -s / -S: same-strand / opposite-strand only
//   - -f / -F: minimum fraction of A / B that must overlap before B contributes
package bedcoverage

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnbed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// recordSource is the minimal record-stream interface Coverage consumes.
// Both *bed.Reader (BED text input) and *alnbed.Reader (BAM/SAM input)
// satisfy it.
type recordSource interface {
	Read() (*bed.Record, error)
}

// sourceReader auto-detects whether r is a SAM/BAM alignment stream or a BED
// text stream and returns the matching record source. Upstream
// `bedtools coverage` accepts BAM on both -a (-abam) and -b; a BAM alignment
// becomes a BED12 record (its CIGAR blocks as BED12 blocks), so the -split
// block-awareness already in Coverage composes for free on the -b side.
func sourceReader(r io.Reader) (recordSource, error) {
	br := bufio.NewReader(r)
	head, _ := br.Peek(16)
	if alnbed.LooksLikeAlignment(head) {
		return alnbed.NewReader(br)
	}
	return bed.NewReader(br), nil
}

// Mode selects which output shape is emitted.
type Mode int

const (
	// ModeDefault emits A + count + covered_bp + len_A + fraction.
	ModeDefault Mode = iota
	// ModeCounts emits A + count only.
	ModeCounts
	// ModeDepth emits one line per base position inside A:
	// A + 1-based position within A + depth.
	ModeDepth
	// ModeHist emits per-depth-bucket histogram lines per A plus an "all"
	// aggregate trailer.
	ModeHist
	// ModeMean / ModeMedian / ModeMin / ModeMax / ModeSum collapse the
	// per-base depth vector via the requested aggregation. Result is appended
	// as a single column after A's columns.
	ModeMean
	ModeMedian
	ModeMin
	ModeMax
	ModeSum
)

// Options controls Coverage behaviour.
type Options struct {
	Mode Mode

	// SameStrand: require A.Strand == B.Strand (skip if A or B has empty
	// strand). Mirrors `bedtools coverage -s`.
	SameStrand bool
	// OppositeStrand: require A.Strand != B.Strand (and both non-empty).
	// Mirrors `bedtools coverage -S`.
	OppositeStrand bool

	// FractionA: minimum fraction of A that must overlap a single B record
	// before that B counts. 0 disables the check.
	FractionA float64
	// FractionB: minimum fraction of B that must overlap A.
	FractionB float64
	// Reciprocal ("-r"): require the overlap fraction be reciprocal for A AND
	// B — B must overlap FractionA of A and A must also overlap FractionA of B
	// (FractionB is taken to equal FractionA). Mirrors `bedtools coverage -r`.
	Reciprocal bool

	// Either ("-e"): require the minimum fraction be satisfied for A OR B,
	// rather than the default AND across the supplied -f/-F thresholds. Mirrors
	// `bedtools coverage -e` (e.g. with -f 0.9 -F 0.1, count when 90% of A OR
	// 10% of B is covered).
	Either bool

	// Split ("-split") makes coverage block-aware. On the database (-b) side it
	// expands BED12 records into their blocks before indexing, so coverage is
	// counted against each block rather than the whole record span. On the
	// query (-a) side a blocked record (a BED12 line or a spliced/N-CIGAR BAM
	// alignment) is split into its sub-blocks: overlap is computed only against
	// those blocks (introns/gaps are excluded), while the reported length-of-A
	// and the per-base depth vector still span the record's full [start,end) —
	// matching upstream bedtools coverage -split (coverageFile.cpp).
	Split bool
}

// Coverage runs the coverage calculation streaming records from readerA,
// indexing readerB into an interval tree first, and writing the result to
// writer. It returns the number of A records processed.
func Coverage(readerA, readerB io.Reader, writer io.Writer, opts Options) (int, error) {
	// Read and index B. The B side auto-detects BAM/SAM vs BED; a BAM
	// alignment arrives as a BED12 record, so -split's block expansion below
	// works for BAM input too.
	srcB, err := sourceReader(readerB)
	if err != nil {
		return 0, fmt.Errorf("error reading B intervals: %w", err)
	}
	bRecords, err := readAll(srcB)
	if err != nil {
		return 0, fmt.Errorf("error reading B intervals: %w", err)
	}
	// -split: expand BED12 database records into their constituent blocks
	// so coverage is counted per block instead of per whole-record span.
	if opts.Split {
		bRecords = expandBlocks(bRecords)
	}
	trees := indexB(bRecords)

	// Stream A (also auto-detecting BAM/SAM vs BED).
	bedReaderA, err := sourceReader(readerA)
	if err != nil {
		return 0, fmt.Errorf("error reading A intervals: %w", err)
	}
	bw := bufio.NewWriter(writer)
	defer bw.Flush()

	// Histogram mode aggregates an "all" footer across all A records.
	allDepthCounts := map[int]int{}
	allLen := 0

	// st holds per-call reusable scratch buffers so the hot loop allocates
	// nothing per A record on the steady-state path (mirrors the bedmerge /
	// bedintersect playbook).
	st := &covState{}

	count := 0
	for {
		recA, err := bedReaderA.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return count, fmt.Errorf("error reading A intervals: %w", err)
		}
		count++

		bMatches := st.selectOverlapping(recA, trees[recA.Chrom], opts)

		// hitCount is the number of contributing B records; depths is the
		// per-base depth vector spanning A's full [start,end). With -split over
		// a blocked query (a BED12 record or a spliced/N-CIGAR BAM alignment),
		// overlap is restricted to A's sub-blocks (introns stay depth 0) and
		// only B records overlapping a block contribute — but the depth vector
		// (and the reported length-of-A) still spans the full record, matching
		// upstream coverageFile.cpp where _queryLen = endPos - startPos.
		var hitCount int
		var depths []int
		if opts.Split && recA.BlockCount > 0 && len(recA.BlockSizes) > 0 {
			hitCount, depths = st.splitDepth(recA, bMatches)
		} else {
			hitCount = len(bMatches)
			depths = st.perBaseDepth(recA, bMatches)
		}

		// prefix holds A's original columns rendered once. For ModeDepth it is
		// reused across every per-base output line, so the (potentially large)
		// column reconstruction happens once per A record rather than once per
		// emitted line.
		prefix := st.recordPrefix(recA)

		switch opts.Mode {
		case ModeCounts:
			st.line = st.line[:0]
			st.line = append(st.line, prefix...)
			st.line = appendIntCol(st.line, hitCount)
			st.line = append(st.line, '\n')
			if _, err := bw.Write(st.line); err != nil {
				return count, err
			}
		case ModeDepth:
			for i, d := range depths {
				st.line = st.line[:0]
				st.line = append(st.line, prefix...)
				st.line = appendIntCol(st.line, i+1)
				st.line = appendIntCol(st.line, d)
				st.line = append(st.line, '\n')
				if _, err := bw.Write(st.line); err != nil {
					return count, err
				}
			}
		case ModeHist:
			counts := map[int]int{}
			for _, d := range depths {
				counts[d]++
				allDepthCounts[d]++
			}
			allLen += len(depths)
			keys := sortedKeys(counts)
			lenA := recA.ChromEnd - recA.ChromStart
			for _, d := range keys {
				bp := counts[d]
				frac := float64(bp) / float64(lenA)
				if lenA == 0 {
					frac = 0
				}
				st.line = st.line[:0]
				st.line = append(st.line, prefix...)
				st.line = appendIntCol(st.line, d)
				st.line = appendIntCol(st.line, bp)
				st.line = appendIntCol(st.line, lenA)
				st.line = appendFracCol(st.line, frac)
				st.line = append(st.line, '\n')
				if _, err := bw.Write(st.line); err != nil {
					return count, err
				}
			}
		case ModeMean, ModeMedian, ModeMin, ModeMax, ModeSum:
			val, ok := depthOp(opts.Mode, depths)
			st.line = st.line[:0]
			st.line = append(st.line, prefix...)
			st.line = append(st.line, '\t')
			switch {
			case !ok:
				st.line = append(st.line, '.')
			case opts.Mode == ModeMean:
				// Upstream `bedtools coverage -mean` accumulates the mean as a
				// 32-bit float and prints it with 7 decimals, so the output
				// carries float32 rounding noise (e.g. 1.3200001). Reproduce it
				// by narrowing to float32 before formatting.
				st.line = strconv.AppendFloat(st.line, float64(float32(val)), 'f', 7, 64)
			default:
				st.line = appendFloatLoose(st.line, val)
			}
			st.line = append(st.line, '\n')
			if _, err := bw.Write(st.line); err != nil {
				return count, err
			}
		default: // ModeDefault
			covered := coveredFromDepths(depths)
			lenA := recA.ChromEnd - recA.ChromStart
			frac := 0.0
			if lenA > 0 {
				frac = float64(covered) / float64(lenA)
			}
			st.line = st.line[:0]
			st.line = append(st.line, prefix...)
			st.line = appendIntCol(st.line, hitCount)
			st.line = appendIntCol(st.line, covered)
			st.line = appendIntCol(st.line, lenA)
			st.line = appendFracCol(st.line, frac)
			st.line = append(st.line, '\n')
			if _, err := bw.Write(st.line); err != nil {
				return count, err
			}
		}
	}

	// Emit "all" footer for hist mode.
	if opts.Mode == ModeHist && allLen > 0 {
		keys := sortedKeys(allDepthCounts)
		for _, d := range keys {
			bp := allDepthCounts[d]
			frac := float64(bp) / float64(allLen)
			if _, err := fmt.Fprintf(bw, "all\t%d\t%d\t%d\t%s\n", d, bp, allLen, formatFraction(frac)); err != nil {
				return count, err
			}
		}
	}

	return count, nil
}

// readAll drains every record from src into a slice, returning io.EOF as a
// clean end (not an error). It is the recordSource equivalent of
// bed.Reader.ReadAll, used so the B side can be either a BED or a BAM/SAM
// stream.
func readAll(src recordSource) ([]*bed.Record, error) {
	var out []*bed.Record
	for {
		rec, err := src.Read()
		if err == io.EOF {
			return out, nil
		}
		if err != nil {
			return out, err
		}
		out = append(out, rec)
	}
}

// expandBlocks replaces each BED12 record carrying block information with one
// record per block (the block interval [ChromStart+BlockStarts[i],
// +BlockSizes[i])); records without blocks pass through unchanged. Used to
// implement the database (-b) side of coverage -split.
func expandBlocks(records []*bed.Record) []*bed.Record {
	out := make([]*bed.Record, 0, len(records))
	for _, r := range records {
		if r.BlockCount <= 0 || len(r.BlockSizes) == 0 {
			out = append(out, r)
			continue
		}
		for i := range r.BlockSizes {
			s := r.ChromStart
			if i < len(r.BlockStarts) {
				s += r.BlockStarts[i]
			}
			block := *r
			block.ChromStart = s
			block.ChromEnd = s + r.BlockSizes[i]
			block.BlockCount = 0
			block.BlockSizes = nil
			block.BlockStarts = nil
			out = append(out, &block)
		}
	}
	return out
}

// indexB returns a per-chromosome interval tree for B. Records are sorted by
// (start, end) within each chromosome so the tree is balanced.
func indexB(records []*bed.Record) map[string]*bed.IntervalTree {
	if len(records) == 0 {
		return nil
	}
	byChrom := map[string][]*bed.Record{}
	for _, r := range records {
		byChrom[r.Chrom] = append(byChrom[r.Chrom], r)
	}
	for chrom := range byChrom {
		recs := byChrom[chrom]
		sort.SliceStable(recs, func(i, j int) bool {
			if recs[i].ChromStart != recs[j].ChromStart {
				return recs[i].ChromStart < recs[j].ChromStart
			}
			return recs[i].ChromEnd < recs[j].ChromEnd
		})
		byChrom[chrom] = recs
	}
	trees := map[string]*bed.IntervalTree{}
	for chrom, recs := range byChrom {
		trees[chrom] = bed.NewIntervalTree(recs)
	}
	return trees
}

// covState carries the reusable scratch buffers a single Coverage call threads
// through its hot loop so the steady-state per-A-record path allocates nothing.
// All buffers are reset (sliced to length 0) before each reuse; their backing
// arrays grow once and are then recycled. A covState is not safe for concurrent
// use — Coverage owns exactly one.
type covState struct {
	matches []*bed.Record // selectOverlapping scratch (filtered B hits)
	depths  []int         // per-base depth vector scratch
	prefix  []byte        // A's original columns rendered as bytes
	line    []byte        // one output line under construction
}

// selectOverlapping is a stateless wrapper over covState.selectOverlapping that
// allocates a fresh result slice. The hot path uses the covState method (which
// recycles scratch); this free function exists for callers/tests that want an
// independently owned slice.
func selectOverlapping(recA *bed.Record, tree *bed.IntervalTree, opts Options) []*bed.Record {
	var st covState
	got := st.selectOverlapping(recA, tree, opts)
	if got == nil {
		return nil
	}
	return append([]*bed.Record(nil), got...)
}

// selectOverlapping returns B records that overlap recA AND pass the
// strand / fraction filters. The returned slice aliases st.matches and is only
// valid until the next call; callers consume it within the same loop iteration.
func (st *covState) selectOverlapping(recA *bed.Record, tree *bed.IntervalTree, opts Options) []*bed.Record {
	if tree == nil {
		return nil
	}
	candidates := tree.Query(recA)
	if len(candidates) == 0 {
		return nil
	}
	out := st.matches[:0]
	for _, b := range candidates {
		if !strandPass(recA, b, opts) {
			continue
		}
		// Under -split, upstream `bedtools coverage` does NOT apply the -f / -F
		// (and hence -r / -e) overlap-fraction thresholds at all: its blocked
		// path (coverageFile.cpp::checkSplits) keeps the BlockMgr *overlapSet*,
		// which is populated for every block intersection regardless of the
		// fraction tests, instead of the fraction-filtered *resultSet* that the
		// plain intersect path uses. So any B that overlaps any A block is
		// counted, whatever -f / -F were given. We therefore skip fractionPass
		// when Split is set (verified empirically against bedtools 2.31.1: even
		// `-f 1.0` / `-F 1.0` / `-r` leave the count unchanged under -split).
		if !opts.Split && !fractionPass(recA, b, opts) {
			continue
		}
		out = append(out, b)
	}
	st.matches = out
	return out
}

// strandPass checks SameStrand / OppositeStrand filters.
func strandPass(a, b *bed.Record, opts Options) bool {
	if opts.SameStrand {
		if a.Strand == "" || b.Strand == "" {
			return false
		}
		if a.Strand != b.Strand {
			return false
		}
	}
	if opts.OppositeStrand {
		if a.Strand == "" || b.Strand == "" {
			return false
		}
		if a.Strand == b.Strand {
			return false
		}
	}
	return true
}

// fractionPass enforces -f / -F (and -r). Returns true if the B record should
// be counted against A.
func fractionPass(a, b *bed.Record, opts Options) bool {
	if opts.FractionA == 0 && opts.FractionB == 0 {
		return true
	}
	overlapStart := a.ChromStart
	if b.ChromStart > overlapStart {
		overlapStart = b.ChromStart
	}
	overlapEnd := a.ChromEnd
	if b.ChromEnd < overlapEnd {
		overlapEnd = b.ChromEnd
	}
	ov := overlapEnd - overlapStart
	if ov <= 0 {
		return false
	}
	lenA := a.ChromEnd - a.ChromStart
	lenB := b.ChromEnd - b.ChromStart
	passA := opts.FractionA == 0 || (lenA > 0 && float64(ov)/float64(lenA) >= opts.FractionA)
	passB := opts.FractionB == 0 || (lenB > 0 && float64(ov)/float64(lenB) >= opts.FractionB)
	if opts.Reciprocal {
		// -r: the A fraction must also hold on the B side (FractionB == FractionA).
		passB = lenB > 0 && float64(ov)/float64(lenB) >= opts.FractionA
		return passA && passB
	}
	if opts.Either {
		// -e: satisfy A OR B (only the thresholds actually supplied count; a
		// 0 threshold means "no constraint", so it must not make OR trivially
		// true). When neither is supplied, fall through to the default.
		if opts.FractionA > 0 || opts.FractionB > 0 {
			eitherA := opts.FractionA > 0 && passA
			eitherB := opts.FractionB > 0 && passB
			return eitherA || eitherB
		}
	}
	// Default `bedtools coverage` semantics: when both -f and -F are given,
	// BOTH must hold (AND across the supplied thresholds).
	return passA && passB
}

// coveredFromDepths returns the number of bases with depth >= 1 in a per-base
// depth vector (the covered-bp column shared by the default mode).
func coveredFromDepths(depths []int) int {
	n := 0
	for _, d := range depths {
		if d > 0 {
			n++
		}
	}
	return n
}

// depthBuf returns a zeroed []int of length n drawn from st.depths, growing the
// backing array only when n exceeds its current capacity. The returned slice
// aliases st.depths and is valid until the next depthBuf call.
func (st *covState) depthBuf(n int) []int {
	if cap(st.depths) < n {
		st.depths = make([]int, n)
	} else {
		st.depths = st.depths[:n]
		for i := range st.depths {
			st.depths[i] = 0
		}
	}
	return st.depths
}

// perBaseDepth returns the per-base depth vector for A's interval given the
// matching B records. The returned slice aliases st.depths (see depthBuf).
func (st *covState) perBaseDepth(a *bed.Record, bs []*bed.Record) []int {
	lenA := a.ChromEnd - a.ChromStart
	if lenA <= 0 {
		return nil
	}
	d := st.depthBuf(lenA)
	for _, b := range bs {
		start := b.ChromStart - a.ChromStart
		end := b.ChromEnd - a.ChromStart
		if start < 0 {
			start = 0
		}
		if end > lenA {
			end = lenA
		}
		for i := start; i < end; i++ {
			d[i]++
		}
	}
	return d
}

// block is a half-open sub-interval [start,end) of a blocked query record, in
// absolute (chromosome) coordinates.
type block struct {
	start int
	end   int
}

// queryBlocks expands a blocked (BED12 or spliced-BAM-derived) record into its
// constituent sub-blocks in absolute coordinates. Each block is
// [ChromStart+BlockStarts[i], +BlockSizes[i]). Records without block info yield
// a single block spanning the whole record. Mirrors upstream GetBedBlocks /
// GetBamBlocks (BlockedIntervals.cpp): M/=/X consume and N skips have already
// been resolved into BlockStarts/BlockSizes by the BED12 parser and by
// pkg/htsgo/alnbed for spliced BAM.
func queryBlocks(a *bed.Record) []block {
	if a.BlockCount <= 0 || len(a.BlockSizes) == 0 {
		return []block{{start: a.ChromStart, end: a.ChromEnd}}
	}
	blocks := make([]block, 0, len(a.BlockSizes))
	for i := range a.BlockSizes {
		s := a.ChromStart
		if i < len(a.BlockStarts) {
			s += a.BlockStarts[i]
		}
		blocks = append(blocks, block{start: s, end: s + a.BlockSizes[i]})
	}
	return blocks
}

// splitDepth computes coverage for a blocked query record under -split. It
// returns (hitCount, depths) where:
//
//   - depths spans A's full [ChromStart,ChromEnd) — intronic bases between
//     blocks stay at depth 0 and still count toward the reported length-of-A,
//     matching upstream coverageFile.cpp (_queryLen = endPos - startPos).
//   - per-base depth is only incremented over the intersection of each B record
//     with A's sub-blocks (gaps/introns are never counted).
//   - hitCount is the number of distinct B records overlapping at least one
//     block, matching upstream's _hitCount after findBlockedOverlaps swaps the
//     hit set for the blocked-overlap set.
//
// splitDepth is a stateless wrapper over covState.splitDepth that returns an
// independently owned depth slice. The hot path uses the covState method (which
// recycles scratch); this free function exists for tests.
func splitDepth(a *bed.Record, bs []*bed.Record) (int, []int) {
	var st covState
	hc, d := st.splitDepth(a, bs)
	if d == nil {
		return hc, nil
	}
	return hc, append([]int(nil), d...)
}

func (st *covState) splitDepth(a *bed.Record, bs []*bed.Record) (int, []int) {
	lenA := a.ChromEnd - a.ChromStart
	if lenA <= 0 {
		return 0, nil
	}
	blocks := queryBlocks(a)
	d := st.depthBuf(lenA)
	hitCount := 0
	for _, b := range bs {
		for _, blk := range blocks {
			start := b.ChromStart
			if blk.start > start {
				start = blk.start
			}
			end := b.ChromEnd
			if blk.end < end {
				end = blk.end
			}
			if end <= start {
				continue
			}
			// Upstream findBlockedOverlaps pushes one overlap sub-interval per
			// (query-block x hit-block) intersection, and makeDepthCount counts
			// _hitCount over those swapped entries. So a single B record that
			// straddles an intron and overlaps two query blocks is counted
			// twice — once per block it touches. Match that by incrementing
			// hitCount per overlapping block, not per B record. (The B side is
			// already expanded to one record per block by expandBlocks under
			// -split, so each b here is a single block.)
			hitCount++
			for i := start - a.ChromStart; i < end-a.ChromStart; i++ {
				d[i]++
			}
		}
	}
	return hitCount, d
}

// depthOp applies a numeric op to the per-base depth vector. Returns
// (value, ok=false) when depths is empty (the only legitimate "no data" case).
func depthOp(mode Mode, depths []int) (float64, bool) {
	if len(depths) == 0 {
		return 0, false
	}
	switch mode {
	case ModeMean:
		sum := 0
		for _, d := range depths {
			sum += d
		}
		return float64(sum) / float64(len(depths)), true
	case ModeMedian:
		sorted := append([]int(nil), depths...)
		sort.Ints(sorted)
		n := len(sorted)
		if n%2 == 1 {
			return float64(sorted[n/2]), true
		}
		return float64(sorted[n/2-1]+sorted[n/2]) / 2, true
	case ModeMin:
		m := depths[0]
		for _, d := range depths[1:] {
			if d < m {
				m = d
			}
		}
		return float64(m), true
	case ModeMax:
		m := depths[0]
		for _, d := range depths[1:] {
			if d > m {
				m = d
			}
		}
		return float64(m), true
	case ModeSum:
		s := 0
		for _, d := range depths {
			s += d
		}
		return float64(s), true
	}
	return 0, false
}

// sortedKeys returns the int keys of m sorted ascending.
func sortedKeys(m map[int]int) []int {
	keys := make([]int, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	return keys
}

// recordPrefix renders A's original columns into st.prefix as a tab-separated
// byte run (no leading or trailing tab) and returns the slice. Output columns
// are appended by the caller via appendIntCol / appendFracCol, each of which
// emits a leading '\t'. The returned slice aliases st.prefix and is valid until
// the next recordPrefix call. This is the byte-buffer equivalent of the former
// recordColumns + strings.Join, avoiding the per-record []string and Itoa heap
// allocations; the emitted bytes are identical to the previous tab-joined form.
//
// Fields beyond chrom/start/end are emitted only when they were populated,
// mirroring the conservative bed.Writer round-trip behaviour.
func (st *covState) recordPrefix(r *bed.Record) []byte {
	out := st.prefix[:0]
	out = append(out, r.Chrom...)
	out = appendIntCol(out, r.ChromStart)
	out = appendIntCol(out, r.ChromEnd)
	// A BED12 record (block information present, e.g. from a BED12 file or a
	// BAM alignment) is echoed as the full 12 columns — including thickStart/
	// thickEnd/itemRgb even when zero. The block lists are echoed VERBATIM:
	// upstream bedtools preserves whatever text was read for the blockSizes /
	// blockStarts columns (a trailing comma is kept if present, omitted if
	// absent), so a record read from BED text round-trips exactly. BAM-derived
	// records carry no raw block text, so they fall back to the trailing-comma
	// form upstream emits for `-abam` (e.g. "50,50,").
	if r.BlockCount != 0 || len(r.BlockSizes) > 0 {
		rgb := r.ItemRGB
		if rgb == "" {
			rgb = "0"
		}
		out = appendStrCol(out, r.Name)
		out = appendIntCol(out, r.Score)
		out = appendStrCol(out, r.Strand)
		out = appendIntCol(out, r.ThickStart)
		out = appendIntCol(out, r.ThickEnd)
		out = appendStrCol(out, rgb)
		out = appendIntCol(out, r.BlockCount)
		out = append(out, '\t')
		out = appendBlockField(out, r.RawBlockSizes, r.BlockSizes)
		out = append(out, '\t')
		out = appendBlockField(out, r.RawBlockStarts, r.BlockStarts)
		for _, ef := range r.ExtraFields {
			out = appendStrCol(out, ef)
		}
		st.prefix = out
		return out
	}
	// The Name/Score/Strand chain only fires once Name is non-empty, matching
	// the conservative BED-aware emit logic in pkg/htsgo/bed.
	if r.Name == "" && r.Score == 0 && r.Strand == "" && len(r.ExtraFields) == 0 {
		st.prefix = out
		return out
	}
	out = appendStrCol(out, r.Name)
	if r.Score != 0 || r.Strand != "" {
		out = appendIntCol(out, r.Score)
	}
	if r.Strand != "" {
		out = appendStrCol(out, r.Strand)
	}
	if r.ThickStart != 0 || r.ThickEnd != 0 {
		out = appendIntCol(out, r.ThickStart)
		out = appendIntCol(out, r.ThickEnd)
	}
	if r.ItemRGB != "" {
		out = appendStrCol(out, r.ItemRGB)
	}
	for _, ef := range r.ExtraFields {
		out = appendStrCol(out, ef)
	}
	st.prefix = out
	return out
}

// recordColumns reconstructs the original column list from a parsed Record as a
// []string. It is the structural twin of recordPrefix and is kept for tests and
// any caller that needs the columns split out; the hot output path uses
// recordPrefix (byte buffer) instead. Both must stay in lockstep.
func recordColumns(r *bed.Record) []string {
	var st covState
	prefix := string(st.recordPrefix(r))
	return strings.Split(prefix, "\t")
}

// appendIntCol appends "\t<n>" to b using strconv.AppendInt (no heap alloc).
func appendIntCol(b []byte, n int) []byte {
	b = append(b, '\t')
	return strconv.AppendInt(b, int64(n), 10)
}

// appendStrCol appends "\t<s>" to b.
func appendStrCol(b []byte, s string) []byte {
	b = append(b, '\t')
	return append(b, s...)
}

// appendFracCol appends "\t<fraction>" using the upstream-faithful 7-decimal
// float32-narrowed formatting (see formatFraction).
func appendFracCol(b []byte, v float64) []byte {
	b = append(b, '\t')
	return strconv.AppendFloat(b, float64(float32(v)), 'f', 7, 64)
}

// appendBlockField renders one BED12 block column (blockSizes or blockStarts)
// into b. See blockField for the raw-vs-synthesized rule.
func appendBlockField(b []byte, raw string, vs []int) []byte {
	if raw != "" {
		return append(b, raw...)
	}
	for _, v := range vs {
		b = strconv.AppendInt(b, int64(v), 10)
		b = append(b, ',')
	}
	return b
}

// appendFloatLoose appends v to b with up to 7 significant digits, trimming
// trailing zeros — the byte-buffer equivalent of formatFloatLoose.
func appendFloatLoose(b []byte, v float64) []byte {
	if v == float64(int64(v)) {
		return strconv.AppendInt(b, int64(v), 10)
	}
	return strconv.AppendFloat(b, v, 'g', -1, 64)
}

// blockField renders one BED12 block column (blockSizes or blockStarts) for
// echo. When raw is non-empty — i.e. the record came from BED text and the
// reader retained the exact column text — it is returned verbatim, so a
// trailing comma is preserved if (and only if) the input had one, matching
// upstream bedtools which echoes the block columns unchanged. When raw is empty
// (e.g. a BAM-derived record, which has no source text), it falls back to the
// trailing-comma form upstream emits for synthesized BED12 records.
func blockField(raw string, vs []int) string {
	if raw != "" {
		return raw
	}
	return joinTrailingComma(vs)
}

// joinTrailingComma renders a block-size/block-start list as
// "v0,v1,...,vN," — the UCSC BED12 form with a trailing comma that bedtools
// echoes verbatim.
func joinTrailingComma(vs []int) string {
	var sb strings.Builder
	for _, v := range vs {
		sb.WriteString(strconv.Itoa(v))
		sb.WriteByte(',')
	}
	return sb.String()
}

// formatFraction prints the fraction column using 7 fixed decimals, matching
// upstream `bedtools coverage` (e.g. "1.0000000", "0.7600000").
//
// Upstream computes the covered-fraction as a 32-bit float (the
// numerator/denominator division happens in float arithmetic in
// coverageFile.cpp / RecordOutputMgr) and prints it with 7 decimals, so the
// last digit carries float32 rounding. For example 7/19 prints as
// "0.3684210", not the float64-rounded "0.3684211". Narrow to float32 before
// formatting to reproduce upstream byte-for-byte.
func formatFraction(v float64) string {
	return strconv.FormatFloat(float64(float32(v)), 'f', 7, 64)
}

// formatFloatLoose prints a number with up to 7 significant digits, trimming
// trailing zeros. Used by the numeric-op output columns where upstream uses
// %g-style formatting.
func formatFloatLoose(v float64) string {
	if v == float64(int64(v)) {
		return strconv.FormatInt(int64(v), 10)
	}
	return strconv.FormatFloat(v, 'g', -1, 64)
}
