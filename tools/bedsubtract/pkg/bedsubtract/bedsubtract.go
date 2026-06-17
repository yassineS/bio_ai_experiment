package bedsubtract

import (
	"bufio"
	"fmt"
	"io"
	"sort"
)

// Options controls how Subtract processes intervals.
type Options struct {
	// RemoveEntire ("-A") drops any A interval that overlaps B at all,
	// rather than splitting it.
	RemoveEntire bool
	// MinFraction ("-f", 0..1). When > 0, an overlap with an individual B
	// interval is only subtracted if it covers at least this fraction of A.
	//
	// With RemoveSum ("-N") set, MinFraction instead becomes the union
	// coverage threshold above which the whole A feature is dropped.
	MinFraction float64
	// Reciprocal ("-r"): require the per-B overlap to cover at least
	// MinFraction of BOTH A and B (upstream sets overlapFractionB =
	// overlapFractionA when -f and -r are given). Must be used solely with -f
	// (upstream rejects -r combined with -F).
	Reciprocal bool
	// RemoveSum ("-N") sums coverage across ALL overlapping B intervals
	// (their union with A) and, if that union covers strictly more than
	// MinFraction of A, drops the entire A feature; otherwise A is emitted
	// unchanged (never split into partial segments). This mirrors
	// upstream bedtools subtract -N (subtractFile.cpp): -f is reused as the
	// union-coverage threshold and per-B fraction filtering is disabled.
	RemoveSum bool
	// SameStrand ("-s"): only consider B intervals on the same strand as A.
	SameStrand bool
	// OppositeStrand ("-S"): only consider B intervals on the opposite
	// strand of A.
	OppositeStrand bool
}

// Validate reports configuration errors.
func (o Options) Validate() error {
	if o.SameStrand && o.OppositeStrand {
		return fmt.Errorf("-s and -S are mutually exclusive")
	}
	if o.MinFraction < 0 || o.MinFraction > 1 {
		return fmt.Errorf("min-fraction must be between 0 and 1, got %v", o.MinFraction)
	}
	// Upstream bedtools subtract -N folds the -f value into the union
	// coverage threshold and then validates it, so -N requires an explicit
	// -f in the open-low/closed-high range (0.0, 1.0].
	if o.RemoveSum && (o.MinFraction <= 0 || o.MinFraction > 1) {
		return fmt.Errorf("-f must be in the range (0.0, 1.0].")
	}
	return nil
}

// Subtract reads BED records from readerA and readerB, subtracts B intervals
// from each A interval, and writes the resulting segments to writer. Output
// preserves the column count of A. Returns the number of output rows.
func Subtract(readerA, readerB io.Reader, writer io.Writer, opts Options) (int, error) {
	if err := opts.Validate(); err != nil {
		return 0, err
	}
	aRows, err := readRows(readerA)
	if err != nil {
		return 0, fmt.Errorf("error reading A: %w", err)
	}
	bRows, err := readRows(readerB)
	if err != nil {
		return 0, fmt.Errorf("error reading B: %w", err)
	}

	// Bucket B by chromosome and sort each bucket by start for predictable
	// processing.
	bByChrom := make(map[string][]*row, 16)
	for _, b := range bRows {
		bByChrom[b.chrom] = append(bByChrom[b.chrom], b)
	}
	for _, list := range bByChrom {
		sort.Slice(list, func(i, j int) bool { return list[i].start < list[j].start })
	}

	bw := bufio.NewWriter(writer)
	count := 0
	for _, a := range aRows {
		segs, drop := subtractOne(a, bByChrom[a.chrom], opts)
		if drop {
			continue
		}
		for _, s := range segs {
			if err := writeRow(bw, s); err != nil {
				return count, fmt.Errorf("error writing output: %w", err)
			}
			count++
		}
	}
	if err := bw.Flush(); err != nil {
		return count, fmt.Errorf("error flushing output: %w", err)
	}
	return count, nil
}

// subtractOne computes the result of subtracting all eligible B intervals
// from a. It returns the resulting segments in input/start order. If drop is
// true the entire A interval is suppressed (used by -A); otherwise segs may
// be empty (A fully consumed by B).
func subtractOne(a *row, bs []*row, opts Options) (segs []*row, drop bool) {
	// Collect eligible B intervals (chromosome already matches).
	var eligible []*row
	aLen := a.length()
	for _, b := range bs {
		if !strandMatch(a, b, opts) {
			continue
		}
		if b.end <= a.start || b.start >= a.end {
			continue
		}
		// Compute overlap.
		ovStart := a.start
		if b.start > ovStart {
			ovStart = b.start
		}
		ovEnd := a.end
		if b.end < ovEnd {
			ovEnd = b.end
		}
		ovLen := ovEnd - ovStart
		if ovLen <= 0 {
			continue
		}
		// In -N (RemoveSum) mode every overlap is collected for the union
		// coverage calculation; the per-B fraction filter does not apply
		// (upstream sets the per-B overlap fraction to 1E-9).
		if !opts.RemoveSum && opts.MinFraction > 0 {
			if aLen > 0 && float64(ovLen)/float64(aLen) < opts.MinFraction {
				continue
			}
			// -r (reciprocal): the overlap must additionally cover at least
			// MinFraction of B. Upstream sets overlapFractionB =
			// overlapFractionA when -f and -r are both given, so the same
			// threshold applies to the B side.
			if opts.Reciprocal {
				bLen := b.length()
				if bLen > 0 && float64(ovLen)/float64(bLen) < opts.MinFraction {
					continue
				}
			}
		}
		eligible = append(eligible, b)
	}

	if len(eligible) == 0 {
		return []*row{a}, false
	}

	// -N / RemoveSum: drop A entirely iff the union of all B overlaps covers
	// strictly more than MinFraction of A; otherwise emit A unchanged.
	if opts.RemoveSum {
		covered := unionCoverage(a, eligible)
		if aLen > 0 && float64(covered)/float64(aLen) > opts.MinFraction {
			return nil, true
		}
		return []*row{a}, false
	}

	if opts.RemoveEntire {
		return nil, true
	}

	// Subtract eligible Bs from [a.start, a.end). Sort by start and process
	// to produce remaining segments.
	sort.Slice(eligible, func(i, j int) bool { return eligible[i].start < eligible[j].start })
	cur := a.start
	for _, b := range eligible {
		bs := b.start
		be := b.end
		if bs < a.start {
			bs = a.start
		}
		if be > a.end {
			be = a.end
		}
		if be <= cur {
			continue
		}
		if bs > cur {
			segs = append(segs, a.withSpan(cur, bs))
		}
		if be > cur {
			cur = be
		}
	}
	if cur < a.end {
		segs = append(segs, a.withSpan(cur, a.end))
	}
	return segs, false
}

// unionCoverage returns the number of bases of a covered by the union of
// the eligible B intervals (each clamped to a's span). It mirrors upstream
// bedtools subtract -N, which paints a per-base covered/uncovered bitmap of
// the query and counts the covered bases.
func unionCoverage(a *row, eligible []*row) int {
	type span struct{ start, end int }
	spans := make([]span, 0, len(eligible))
	for _, b := range eligible {
		s := b.start
		if s < a.start {
			s = a.start
		}
		e := b.end
		if e > a.end {
			e = a.end
		}
		if e > s {
			spans = append(spans, span{s, e})
		}
	}
	if len(spans) == 0 {
		return 0
	}
	sort.Slice(spans, func(i, j int) bool { return spans[i].start < spans[j].start })
	covered := 0
	curStart, curEnd := spans[0].start, spans[0].end
	for _, sp := range spans[1:] {
		if sp.start > curEnd {
			covered += curEnd - curStart
			curStart, curEnd = sp.start, sp.end
			continue
		}
		if sp.end > curEnd {
			curEnd = sp.end
		}
	}
	covered += curEnd - curStart
	return covered
}

// strandMatch returns true if b should be considered for subtraction from a
// under the strand options.
func strandMatch(a, b *row, opts Options) bool {
	if !opts.SameStrand && !opts.OppositeStrand {
		return true
	}
	if a.strand == "" || a.strand == "." || b.strand == "" || b.strand == "." {
		// With strand filtering enabled, missing or unknown strand on
		// either side means "cannot determine same/opposite": exclude
		// the B interval (matches bedtools subtract behaviour).
		return false
	}
	if opts.SameStrand {
		return a.strand == b.strand
	}
	return a.strand != b.strand
}
