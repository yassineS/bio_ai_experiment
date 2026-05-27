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

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

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
	// Reciprocal: when true, both FractionA and FractionB must be satisfied
	// (default behaviour is "either satisfies if non-zero").
	Reciprocal bool
}

// Coverage runs the coverage calculation streaming records from readerA,
// indexing readerB into an interval tree first, and writing the result to
// writer. It returns the number of A records processed.
func Coverage(readerA, readerB io.Reader, writer io.Writer, opts Options) (int, error) {
	// Read and index B.
	bedReaderB := bed.NewReader(readerB)
	bRecords, err := bedReaderB.ReadAll()
	if err != nil {
		return 0, fmt.Errorf("error reading B intervals: %w", err)
	}
	trees := indexB(bRecords)

	// Stream A.
	bedReaderA := bed.NewReader(readerA)
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
		switch opts.Mode {
		case ModeCounts:
			if err := writeWithExtra(bw, recA, strconv.Itoa(len(bMatches))); err != nil {
				return count, err
			}
		case ModeDepth:
			depths := perBaseDepth(recA, bMatches)
			for i, d := range depths {
				if err := writeWithExtra(bw, recA, strconv.Itoa(i+1), strconv.Itoa(d)); err != nil {
					return count, err
				}
			}
		case ModeHist:
			depths := perBaseDepth(recA, bMatches)
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
			depths := perBaseDepth(recA, bMatches)
			val, ok := depthOp(opts.Mode, depths)
			var s string
			if !ok {
				s = "."
			} else {
				s = formatFloatLoose(val)
			}
			if err := writeWithExtra(bw, recA, s); err != nil {
				return count, err
			}
		default: // ModeDefault
			covered := coveredBases(recA, bMatches)
			lenA := recA.ChromEnd - recA.ChromStart
			frac := 0.0
			if lenA > 0 {
				frac = float64(covered) / float64(lenA)
			}
			if err := writeWithExtra(bw, recA,
				strconv.Itoa(len(bMatches)),
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
		if !fractionPass(recA, b, opts) {
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
		return passA && passB
	}
	// Match `bedtools coverage` semantics: when both -f and -F are given
	// without -r, BOTH must still hold (the default is AND across the
	// supplied thresholds). This matches upstream too — `-r` only matters
	// when one of -f/-F is supplied by itself.
	return passA && passB
}

// coveredBases returns the number of bases in A covered by at least one of
// the matching B records (depth >= 1).
func coveredBases(a *bed.Record, bs []*bed.Record) int {
	lenA := a.ChromEnd - a.ChromStart
	if lenA <= 0 || len(bs) == 0 {
		return 0
	}
	covered := make([]bool, lenA)
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
			covered[i] = true
		}
	}
	n := 0
	for _, c := range covered {
		if c {
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
	// The Name/Score/Strand chain only fires once Name is non-empty, matching
	// the conservative BED12-aware emit logic in pkg/htsgo/bed.
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
	if r.BlockCount != 0 {
		out = append(out, strconv.Itoa(r.BlockCount))
	}
	if len(r.BlockSizes) > 0 {
		sizes := make([]string, len(r.BlockSizes))
		for i, s := range r.BlockSizes {
			sizes[i] = strconv.Itoa(s)
		}
		// BED12 convention preserves the trailing comma on block-size /
		// block-start columns; upstream bedtools coverage emits with it.
		out = append(out, strings.Join(sizes, ",")+",")
	}
	if len(r.BlockStarts) > 0 {
		starts := make([]string, len(r.BlockStarts))
		for i, s := range r.BlockStarts {
			starts[i] = strconv.Itoa(s)
		}
		out = append(out, strings.Join(starts, ",")+",")
	}
	if len(r.ExtraFields) > 0 {
		out = append(out, r.ExtraFields...)
	}
	return out
}

// formatFraction prints the fraction column using 7 fixed decimals, matching
// upstream `bedtools coverage` (e.g. "1.0000000", "0.7600000").
func formatFraction(v float64) string {
	return strconv.FormatFloat(v, 'f', 7, 64)
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
