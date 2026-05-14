// Package bedjaccard computes the Jaccard similarity between two sorted BED
// files, mirroring `bedtools jaccard`.
//
// Given two BED files A and B that are pre-sorted by (chrom, start), it
// computes:
//
//   - intersection: total bases shared between A and B
//   - union:        |A| + |B| - intersection
//   - jaccard:      intersection / union (0 if union is zero)
//   - n:            number of (A, B) interval pairs that overlap
//
// The algorithm performs a single linear sweep, holding only a small "active"
// window of B intervals at any time, so memory use is independent of file
// size beyond the sort prerequisite.
package bedjaccard

import (
	"errors"
	"fmt"
	"io"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// Options configures the Jaccard computation.
type Options struct {
	// SameStrand causes only same-strand B intervals to be considered for an A.
	SameStrand bool
	// OppositeStrand causes only opposite-strand B intervals to be considered.
	OppositeStrand bool
	// FractionA: require at least this fraction of A to overlap B for the
	// pair to count (0..1). Zero disables the check.
	FractionA float64
	// FractionB: require at least this fraction of B to overlap A.
	FractionB float64
}

// Result is the one-line summary written by Run.
type Result struct {
	Intersection int
	Union        int
	Jaccard      float64
	N            int
}

// Run reads sorted BED records from a and b, computes the Jaccard summary, and
// writes a two-line tab-separated table (header then values) to w. The Result
// is also returned for programmatic use.
func Run(a, b io.Reader, w io.Writer, opts Options) (*Result, error) {
	if opts.SameStrand && opts.OppositeStrand {
		return nil, errors.New("-s and -S are mutually exclusive")
	}
	if opts.FractionA < 0 || opts.FractionA > 1 {
		return nil, fmt.Errorf("-f must be in [0,1], got %v", opts.FractionA)
	}
	if opts.FractionB < 0 || opts.FractionB > 1 {
		return nil, fmt.Errorf("-F must be in [0,1], got %v", opts.FractionB)
	}

	res, err := jaccard(a, b, opts)
	if err != nil {
		return nil, err
	}

	if _, err := fmt.Fprintln(w, "intersection\tunion\tjaccard\tn_intersections"); err != nil {
		return nil, err
	}
	if _, err := fmt.Fprintf(w, "%d\t%d\t%s\t%d\n",
		res.Intersection, res.Union, formatJaccard(res.Jaccard), res.N); err != nil {
		return nil, err
	}
	return res, nil
}

// formatJaccard renders the ratio with C++ ostream's default precision
// (6 significant digits with %g-style trimming), which is what upstream
// `bedtools jaccard` uses when it prints the ratio via `cout`.
func formatJaccard(j float64) string {
	return strconv.FormatFloat(j, 'g', 6, 64)
}

// jaccard does the streaming sweep.
func jaccard(aReader, bReader io.Reader, opts Options) (*Result, error) {
	ra := bed.NewReader(aReader)
	rb := bed.NewReader(bReader)

	var active []*bed.Record
	var (
		lastA, lastB *bed.Record
		bExhausted   bool

		totalA, totalB, totalIntersect, totalPairs int
	)

	for {
		recA, err := ra.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("reading A: %w", err)
		}
		if lastA != nil && !sortedAfter(lastA, recA) {
			return nil, fmt.Errorf("input A is not sorted: %s:%d..%d before %s:%d..%d",
				lastA.Chrom, lastA.ChromStart, lastA.ChromEnd,
				recA.Chrom, recA.ChromStart, recA.ChromEnd)
		}
		lastA = recA
		totalA += recA.ChromEnd - recA.ChromStart

		// Pull more B until the active window covers recA's territory.
		for !bExhausted && needMoreB(active, recA) {
			recB, err := rb.Read()
			if err == io.EOF {
				bExhausted = true
				break
			}
			if err != nil {
				return nil, fmt.Errorf("reading B: %w", err)
			}
			if lastB != nil && !sortedAfter(lastB, recB) {
				return nil, fmt.Errorf("input B is not sorted: %s:%d..%d before %s:%d..%d",
					lastB.Chrom, lastB.ChromStart, lastB.ChromEnd,
					recB.Chrom, recB.ChromStart, recB.ChromEnd)
			}
			lastB = recB
			totalB += recB.ChromEnd - recB.ChromStart
			active = append(active, recB)
		}

		// Drop B intervals that can no longer match recA or any later A.
		active = pruneActive(active, recA)

		// Score recA against the current active window.
		for _, b := range active {
			if b.Chrom != recA.Chrom {
				continue
			}
			if !strandOK(recA, b, opts) {
				continue
			}
			start := recA.ChromStart
			if b.ChromStart > start {
				start = b.ChromStart
			}
			end := recA.ChromEnd
			if b.ChromEnd < end {
				end = b.ChromEnd
			}
			if end <= start {
				continue
			}
			ov := end - start
			if !fractionOK(recA, b, ov, opts) {
				continue
			}
			totalIntersect += ov
			totalPairs++
		}
	}

	// Drain any remaining B records so we account for |B|.
	for !bExhausted {
		recB, err := rb.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("reading B: %w", err)
		}
		if lastB != nil && !sortedAfter(lastB, recB) {
			return nil, fmt.Errorf("input B is not sorted: %s:%d..%d before %s:%d..%d",
				lastB.Chrom, lastB.ChromStart, lastB.ChromEnd,
				recB.Chrom, recB.ChromStart, recB.ChromEnd)
		}
		lastB = recB
		totalB += recB.ChromEnd - recB.ChromStart
	}

	union := totalA + totalB - totalIntersect
	jacc := 0.0
	if union > 0 {
		jacc = float64(totalIntersect) / float64(union)
	}
	return &Result{
		Intersection: totalIntersect,
		Union:        union,
		Jaccard:      jacc,
		N:            totalPairs,
	}, nil
}

// needMoreB returns true if more B records should be pulled to fully cover
// candidates for recA: when active is empty, the last B is on an earlier chrom,
// or the last B starts before recA ends on the same chrom.
func needMoreB(active []*bed.Record, recA *bed.Record) bool {
	if len(active) == 0 {
		return true
	}
	last := active[len(active)-1]
	if last.Chrom != recA.Chrom {
		return last.Chrom < recA.Chrom
	}
	return last.ChromStart < recA.ChromEnd
}

// pruneActive removes B intervals that can no longer overlap recA or any
// subsequent A (A is sorted).
func pruneActive(active []*bed.Record, recA *bed.Record) []*bed.Record {
	out := active[:0]
	for _, b := range active {
		if b.Chrom < recA.Chrom {
			continue
		}
		if b.Chrom == recA.Chrom && b.ChromEnd <= recA.ChromStart {
			continue
		}
		out = append(out, b)
	}
	return out
}

// sortedAfter checks whether next comes at or after prev in (chrom, start) order.
func sortedAfter(prev, next *bed.Record) bool {
	if prev.Chrom != next.Chrom {
		return prev.Chrom < next.Chrom
	}
	return next.ChromStart >= prev.ChromStart
}

// strandOK applies -s/-S filtering. If neither flag is set, any strand
// combination passes. Empty/dot strands never match an explicit filter
// (matching bedtools' behaviour: BED6 is required for -s/-S).
func strandOK(a, b *bed.Record, opts Options) bool {
	switch {
	case opts.SameStrand:
		if a.Strand == "" || a.Strand == "." || b.Strand == "" || b.Strand == "." {
			return false
		}
		return a.Strand == b.Strand
	case opts.OppositeStrand:
		if a.Strand == "" || a.Strand == "." || b.Strand == "" || b.Strand == "." {
			return false
		}
		return (a.Strand == "+" && b.Strand == "-") || (a.Strand == "-" && b.Strand == "+")
	}
	return true
}

// fractionOK applies -f / -F overlap-fraction filters.
func fractionOK(a, b *bed.Record, overlap int, opts Options) bool {
	if opts.FractionA > 0 {
		lenA := a.ChromEnd - a.ChromStart
		if lenA == 0 || float64(overlap)/float64(lenA) < opts.FractionA {
			return false
		}
	}
	if opts.FractionB > 0 {
		lenB := b.ChromEnd - b.ChromStart
		if lenB == 0 || float64(overlap)/float64(lenB) < opts.FractionB {
			return false
		}
	}
	return true
}
