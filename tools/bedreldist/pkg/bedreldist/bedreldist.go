// Package bedreldist computes the distribution of "relative distances"
// between intervals in two BED files, mirroring `bedtools reldist`.
//
// For each interval a in A, we find the two database (B) interval midpoints
// that flank a's midpoint on the same chromosome. The relative distance is
//
//	min(|m - left|, |m - right|) / (right - left)
//
// where m is a's midpoint and left/right are the nearest B-midpoints to its
// left and right. By construction this value lies in [0.0, 0.5].
//
// The default output is a histogram with bins of width 0.01:
//
//	reldist\tcount\ttotal\tfraction
//	0.00\t<count>\t<total>\t<fraction>
//	...
//
// With Detail=true, each A interval gets a per-row line of the form
//
//	chrom\tstart\tend[\t...other A fields...]\t<reldist>
//
// where <reldist> is formatted with %.3f, matching upstream.
package bedreldist

import (
	"fmt"
	"io"
	"math"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// Options controls bedreldist behaviour.
type Options struct {
	// Detail, when true, switches output from the histogram to a per-A row
	// containing the input fields followed by the relative distance.
	Detail bool
}

// HistogramBin holds the count for a single 0.01-wide reldist bucket.
type HistogramBin struct {
	// Bin is the lower edge of the bucket (0.00 .. 0.50 inclusive).
	Bin float64
	// Count is the number of A intervals that fell into this bucket.
	Count int
}

// Result captures the histogram bins and the number of A intervals that
// contributed to the summary. Intervals whose chromosome is absent from B,
// or that fall before the first B-midpoint or after the last B-midpoint, do
// not contribute (matching upstream).
type Result struct {
	Bins  []HistogramBin
	Total int
}

// Run reads BED records from a (queries) and b (database), and writes either
// the histogram (default) or per-A detail rows to w. It returns the
// histogram-mode summary regardless of Detail; callers that only need a
// programmatic result can ignore the write side.
func Run(a, b io.Reader, w io.Writer, opts Options) (*Result, error) {
	// Load and sort all B midpoints by chromosome.
	mids, err := loadMidpoints(b)
	if err != nil {
		return nil, fmt.Errorf("reading B: %w", err)
	}

	ra := bed.NewReader(a)
	counts := map[int]int{} // key: floor(rel_dist * 100), 0..50
	total := 0

	for {
		rec, err := ra.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("reading A: %w", err)
		}

		chromMids, ok := mids[rec.Chrom]
		if !ok || len(chromMids) == 0 {
			continue
		}
		midpoint := (rec.ChromStart + rec.ChromEnd) / 2

		// lower_bound: first index whose value >= midpoint.
		idx := sort.Search(len(chromMids), func(i int) bool {
			return chromMids[i] >= midpoint
		})
		// Upstream: low_idx = idx if idx == 0 else idx-1; high_idx = low_idx+1.
		var lowIdx int
		if idx == 0 {
			lowIdx = 0
		} else {
			lowIdx = idx - 1
		}
		highIdx := lowIdx + 1
		// "make sure we don't run off the boundaries"
		if lowIdx == len(chromMids)-1 {
			continue
		}
		left := chromMids[lowIdx]
		right := chromMids[highIdx]
		// Upstream guards against `left > midpoint` (happens at idx==0 when the
		// query sits to the left of every B midpoint).
		if left > midpoint {
			continue
		}
		leftDist := abs(midpoint - left)
		rightDist := abs(midpoint - right)
		var rel float64
		minDist := leftDist
		if rightDist < minDist {
			minDist = rightDist
		}
		if minDist == 0 {
			rel = 0.0
		} else {
			// right-left is guaranteed > 0 because chromMids is sorted and
			// lowIdx < highIdx (we already skipped the lowIdx == last case).
			// If chromMids has a run of duplicates with right == left, treat
			// as zero relative distance, since minDist is also 0 in that case
			// (otherwise we'd divide by zero).
			if right == left {
				rel = 0.0
			} else {
				rel = float64(minDist) / float64(right-left)
			}
		}
		if opts.Detail {
			if err := writeDetail(w, rec, rel); err != nil {
				return nil, err
			}
		}
		// Round down to two decimals, matching `floor(rel*100)/100`.
		bin := int(math.Floor(rel * 100))
		if bin < 0 {
			bin = 0
		}
		if bin > 50 {
			bin = 50
		}
		counts[bin]++
		total++
	}

	res := &Result{Total: total}
	// Emit bins in ascending order. Upstream uses `map<double,size_t>` which
	// iterates in key-sorted order and ONLY emits non-empty bins; replicate
	// that exactly so histogram output is byte-for-byte identical.
	keys := make([]int, 0, len(counts))
	for k := range counts {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	for _, k := range keys {
		res.Bins = append(res.Bins, HistogramBin{Bin: float64(k) / 100.0, Count: counts[k]})
	}

	if !opts.Detail {
		if err := writeHistogram(w, res); err != nil {
			return nil, err
		}
	}
	return res, nil
}

func writeHistogram(w io.Writer, res *Result) error {
	if _, err := fmt.Fprint(w, "reldist\tcount\ttotal\tfraction\n"); err != nil {
		return err
	}
	if res.Total == 0 {
		return nil
	}
	for _, b := range res.Bins {
		frac := float64(b.Count) / float64(res.Total)
		if _, err := fmt.Fprintf(w, "%.2f\t%d\t%d\t%.3f\n", b.Bin, b.Count, res.Total, frac); err != nil {
			return err
		}
	}
	return nil
}

// writeDetail mirrors upstream `_bedA->reportBedTab(bed); printf("%.3lf\n", rel_dist);`
// We emit chrom, start, end (always), then any of name/score/strand that are
// non-default, then the relative distance. The BED reader does not preserve
// the original raw line, so the emitted prefix may differ from the input on
// records that mix populated and skipped optional fields; for typical BED3-6
// inputs this matches.
func writeDetail(w io.Writer, rec *bed.Record, rel float64) error {
	if _, err := fmt.Fprintf(w, "%s\t%d\t%d", rec.Chrom, rec.ChromStart, rec.ChromEnd); err != nil {
		return err
	}
	if rec.Name != "" {
		if _, err := fmt.Fprintf(w, "\t%s", rec.Name); err != nil {
			return err
		}
		if rec.Strand != "" || rec.Score != 0 {
			if _, err := fmt.Fprintf(w, "\t%d", rec.Score); err != nil {
				return err
			}
			if rec.Strand != "" {
				if _, err := fmt.Fprintf(w, "\t%s", rec.Strand); err != nil {
					return err
				}
			}
		}
	}
	_, err := fmt.Fprintf(w, "\t%.3f\n", rel)
	return err
}

// loadMidpoints reads every record in b and returns a per-chromosome,
// ascending-sorted slice of midpoints.
func loadMidpoints(r io.Reader) (map[string][]int, error) {
	rb := bed.NewReader(r)
	mids := map[string][]int{}
	for {
		rec, err := rb.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		mp := (rec.ChromStart + rec.ChromEnd) / 2
		mids[rec.Chrom] = append(mids[rec.Chrom], mp)
	}
	for k := range mids {
		sort.Ints(mids[k])
	}
	return mids, nil
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}
