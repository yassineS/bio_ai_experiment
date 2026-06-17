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

		bMatches := selectOverlapping(recA, trees[recA.Chrom], opts)

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
			hitCount, depths = splitDepth(recA, bMatches)
		} else {
			hitCount = len(bMatches)
			depths = perBaseDepth(recA, bMatches)
		}

		switch opts.Mode {
		case ModeCounts:
			if err := writeWithExtra(bw, recA, strconv.Itoa(hitCount)); err != nil {
				return count, err
			}
		case ModeDepth:
			for i, d := range depths {
				if err := writeWithExtra(bw, recA, strconv.Itoa(i+1), strconv.Itoa(d)); err != nil {
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
				if err := writeWithExtra(bw, recA,
					strconv.Itoa(d),
					strconv.Itoa(bp),
					strconv.Itoa(lenA),
					formatFraction(frac),
				); err != nil {
					return count, err
				}
			}
		case ModeMean, ModeMedian, ModeMin, ModeMax, ModeSum:
			val, ok := depthOp(opts.Mode, depths)
			var s string
			switch {
			case !ok:
				s = "."
			case opts.Mode == ModeMean:
				// Upstream `bedtools coverage -mean` accumulates the mean as a
				// 32-bit float and prints it with 7 decimals, so the output
				// carries float32 rounding noise (e.g. 1.3200001). Reproduce it
				// by narrowing to float32 before formatting.
				s = strconv.FormatFloat(float64(float32(val)), 'f', 7, 64)
			default:
				s = formatFloatLoose(val)
			}
			if err := writeWithExtra(bw, recA, s); err != nil {
				return count, err
			}
		default: // ModeDefault
			covered := coveredFromDepths(depths)
			lenA := recA.ChromEnd - recA.ChromStart
			frac := 0.0
			if lenA > 0 {
				frac = float64(covered) / float64(lenA)
			}
			if err := writeWithExtra(bw, recA,
				strconv.Itoa(hitCount),
				strconv.Itoa(covered),
				strconv.Itoa(lenA),
				formatFraction(frac),
			); err != nil {
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

// selectOverlapping returns B records that overlap recA AND pass the
// strand / fraction filters.
func selectOverlapping(recA *bed.Record, tree *bed.IntervalTree, opts Options) []*bed.Record {
	if tree == nil {
		return nil
	}
	candidates := tree.Query(recA)
	if len(candidates) == 0 {
		return nil
	}
	out := candidates[:0:0] // new slice; never reuse the tree-owned one
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

// perBaseDepth returns the per-base depth vector for A's interval given the
// matching B records.
func perBaseDepth(a *bed.Record, bs []*bed.Record) []int {
	lenA := a.ChromEnd - a.ChromStart
	if lenA <= 0 {
		return nil
	}
	d := make([]int, lenA)
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
func splitDepth(a *bed.Record, bs []*bed.Record) (int, []int) {
	lenA := a.ChromEnd - a.ChromStart
	if lenA <= 0 {
		return 0, nil
	}
	blocks := queryBlocks(a)
	d := make([]int, lenA)
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

// writeWithExtra emits a tab-separated line: A's original columns followed by
// the extra columns. Uses the record's parsed fields rather than re-stringing
// every BED field because the BED writer is conservative about how many
// optional columns it emits (a 6-field input must round-trip as 6 columns).
func writeWithExtra(w io.Writer, r *bed.Record, extra ...string) error {
	cols := recordColumns(r)
	all := append(cols, extra...)
	_, err := fmt.Fprintln(w, strings.Join(all, "\t"))
	return err
}

// recordColumns reconstructs the original column list from a parsed Record.
// Fields beyond chrom/start/end are emitted only when they were populated.
// This mirrors the bed.Writer behaviour, but keeps the columns as a []string
// so we can append further columns cleanly.
func recordColumns(r *bed.Record) []string {
	out := []string{r.Chrom, strconv.Itoa(r.ChromStart), strconv.Itoa(r.ChromEnd)}
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
		out = append(out,
			r.Name, strconv.Itoa(r.Score), r.Strand,
			strconv.Itoa(r.ThickStart), strconv.Itoa(r.ThickEnd), rgb,
			strconv.Itoa(r.BlockCount),
			blockField(r.RawBlockSizes, r.BlockSizes),
			blockField(r.RawBlockStarts, r.BlockStarts),
		)
		out = append(out, r.ExtraFields...)
		return out
	}
	// The Name/Score/Strand chain only fires once Name is non-empty, matching
	// the conservative BED-aware emit logic in pkg/htsgo/bed.
	if r.Name == "" && r.Score == 0 && r.Strand == "" && len(r.ExtraFields) == 0 {
		return out
	}
	out = append(out, r.Name)
	if r.Score != 0 || r.Strand != "" {
		out = append(out, strconv.Itoa(r.Score))
	}
	if r.Strand != "" {
		out = append(out, r.Strand)
	}
	if r.ThickStart != 0 || r.ThickEnd != 0 {
		out = append(out, strconv.Itoa(r.ThickStart), strconv.Itoa(r.ThickEnd))
	}
	if r.ItemRGB != "" {
		out = append(out, r.ItemRGB)
	}
	if len(r.ExtraFields) > 0 {
		out = append(out, r.ExtraFields...)
	}
	return out
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
