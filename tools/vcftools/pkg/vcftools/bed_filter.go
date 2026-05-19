// BED-based site filtering for vcftools.
//
// `--bed FILE` keeps only sites whose 1-based VCF POS falls inside any
// interval listed in the BED file. `--exclude-bed FILE` is the inverse.
//
// BED files are 0-based half-open ([start, end)) per the UCSC spec. A VCF
// site at position P (1-based) is considered "inside" interval [s, e) when
//
//	s < P <= e         (equivalently, s <= P-1 < e)
//
// This mirrors upstream vcftools, which treats the VCF position as a single
// 1-bp feature for membership testing and does not consider the reference
// allele length. Tools that want overlap semantics on indels can post-
// process the recoded VCF.
//
// We re-use `pkg/htsgo/bed` for parsing (so `.gz` files work via
// `pkg/htsgo/iohelper`). Intervals are stored as a sorted slice per
// chromosome and queried with binary search, giving O(log N) membership
// tests after an O(N log N) load. Overlapping or out-of-order intervals in
// the input are tolerated — we merge during load so each chromosome holds a
// disjoint sorted set.
package vcftools

import (
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// bedInterval is a single half-open [start, end) interval in BED coordinates.
type bedInterval struct {
	start int
	end   int
}

// bedRegions holds per-chromosome sorted, merged intervals.
type bedRegions struct {
	byChrom map[string][]bedInterval
}

// loadBedRegions parses a BED file (optionally `.gz`) into per-chromosome
// merged-and-sorted intervals. Empty files produce an empty (non-nil) result.
func loadBedRegions(filename string) (*bedRegions, error) {
	f, err := iohelper.OpenReader(filename)
	if err != nil {
		return nil, fmt.Errorf("opening BED file %s: %w", filename, err)
	}
	defer f.Close()

	reader := bed.NewReader(f)
	raw := make(map[string][]bedInterval)
	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("parsing BED file %s: %w", filename, err)
		}
		if rec.ChromEnd <= rec.ChromStart {
			// Skip degenerate / zero-length intervals; vcftools also drops
			// these silently.
			continue
		}
		raw[rec.Chrom] = append(raw[rec.Chrom], bedInterval{
			start: rec.ChromStart,
			end:   rec.ChromEnd,
		})
	}

	r := &bedRegions{byChrom: make(map[string][]bedInterval, len(raw))}
	for chrom, ivs := range raw {
		r.byChrom[chrom] = mergeIntervals(ivs)
	}
	return r, nil
}

// mergeIntervals sorts and merges overlapping intervals so each chromosome's
// list is disjoint and sorted by start.
func mergeIntervals(ivs []bedInterval) []bedInterval {
	if len(ivs) == 0 {
		return ivs
	}
	sort.Slice(ivs, func(i, j int) bool {
		if ivs[i].start != ivs[j].start {
			return ivs[i].start < ivs[j].start
		}
		return ivs[i].end < ivs[j].end
	})
	merged := ivs[:0]
	cur := ivs[0]
	for _, iv := range ivs[1:] {
		if iv.start <= cur.end {
			if iv.end > cur.end {
				cur.end = iv.end
			}
			continue
		}
		merged = append(merged, cur)
		cur = iv
	}
	merged = append(merged, cur)
	return merged
}

// containsVCFPos reports whether the 1-based VCF position pos is inside any
// interval on chromosome chrom. Translating to BED coordinates: the position
// is "inside" [s, e) iff s <= pos-1 < e, i.e. s < pos <= e.
func (r *bedRegions) containsVCFPos(chrom string, pos int) bool {
	if r == nil {
		return false
	}
	ivs, ok := r.byChrom[chrom]
	if !ok {
		return false
	}
	zero := pos - 1
	// Find the first interval whose end is strictly greater than zero
	// (i.e. could contain zero). Since intervals are disjoint and sorted by
	// start, the candidate is at this index.
	i := sort.Search(len(ivs), func(i int) bool {
		return ivs[i].end > zero
	})
	if i == len(ivs) {
		return false
	}
	return ivs[i].start <= zero && zero < ivs[i].end
}
