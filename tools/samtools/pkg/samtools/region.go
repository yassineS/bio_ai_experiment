package samtools

import (
	"fmt"
	"sort"
	"strconv"
	"strings"
)

// Region is a parsed region-query specifier: a chromosome name and a
// half-open 1-based coordinate range. End == 0 means "to the end of the
// chromosome".
type Region struct {
	Chrom string
	// Beg is 1-based inclusive (matching SAM conventions for the CLI).
	Beg int
	// End is 1-based inclusive. End == 0 means open-ended.
	End int
}

// ParseRegion parses a region specifier in the conventional samtools form:
//
//	chrom
//	chrom:start
//	chrom:start-end
//	chrom:start-
//
// Underscores and digit-separators in coordinates (e.g. "1,000") are
// stripped to match upstream samtools. start defaults to 1; an absent or
// trailing-dash end defaults to "end of chromosome".
func ParseRegion(spec string) (Region, error) {
	if spec == "" {
		return Region{}, fmt.Errorf("samtools: empty region")
	}
	// Split on the last ':' to allow ':'-containing chrom names. Real SAM
	// chrom names rarely contain ':' but we keep the rule consistent with
	// htslib's hts_parse_reg.
	colon := strings.LastIndex(spec, ":")
	if colon < 0 {
		return Region{Chrom: spec, Beg: 1, End: 0}, nil
	}
	chrom := spec[:colon]
	rest := spec[colon+1:]
	if chrom == "" {
		return Region{}, fmt.Errorf("samtools: empty chromosome in region %q", spec)
	}
	// rest must be "start" or "start-" or "start-end".
	beg := 1
	end := 0
	if rest == "" {
		// "chr:" — treated as the whole chromosome.
		return Region{Chrom: chrom, Beg: 1, End: 0}, nil
	}
	dash := strings.Index(rest, "-")
	startStr := rest
	endStr := ""
	if dash >= 0 {
		startStr = rest[:dash]
		endStr = rest[dash+1:]
	}
	if startStr != "" {
		n, err := parseCoord(startStr)
		if err != nil {
			return Region{}, fmt.Errorf("samtools: bad region start %q: %w", startStr, err)
		}
		beg = n
	}
	if endStr != "" {
		n, err := parseCoord(endStr)
		if err != nil {
			return Region{}, fmt.Errorf("samtools: bad region end %q: %w", endStr, err)
		}
		end = n
	}
	if beg < 1 {
		beg = 1
	}
	if end > 0 && end < beg {
		return Region{}, fmt.Errorf("samtools: region end %d < beg %d", end, beg)
	}
	return Region{Chrom: chrom, Beg: beg, End: end}, nil
}

// parseCoord strips digit-separators (commas, underscores) and parses an
// unsigned integer.
func parseCoord(s string) (int, error) {
	if s == "" {
		return 0, fmt.Errorf("empty coordinate")
	}
	clean := strings.NewReplacer(",", "", "_", "").Replace(s)
	n, err := strconv.Atoi(clean)
	if err != nil {
		return 0, err
	}
	return n, nil
}

// Overlaps reports whether [pos, pos+refLen) on 0-based 'refID' overlaps the
// region's half-open 0-based [Beg-1, End) (with End == 0 meaning "open").
func (r Region) OverlapsRef(refID int, refIDForChrom int, pos0 int, refLen int) bool {
	if refID != refIDForChrom {
		return false
	}
	// Compute the half-open interval in 0-based form.
	regBeg0 := r.Beg - 1
	if regBeg0 < 0 {
		regBeg0 = 0
	}
	regEnd0 := r.End
	if regEnd0 <= 0 {
		regEnd0 = 1 << 30 // effectively infinity for SAM coordinates.
	}
	recBeg := pos0
	recEnd := pos0 + refLen
	if recEnd <= recBeg {
		recEnd = recBeg + 1
	}
	return recBeg < regEnd0 && recEnd > regBeg0
}

// RegionChunks returns the union of BAI chunks across regions (deduplicated
// and merged). It also returns the parallel list of (refID, 0-based half-
// open beg, end) tuples for caller-side overlap filtering after a chunk has
// been read.
type ResolvedRegion struct {
	Region Region
	RefID  int
	Beg0   int
	End0   int // exclusive, half-open
}

// ResolveRegions parses every spec and resolves each chrom name to a refID
// via the provided lookup. Regions referencing unknown chromosomes are
// dropped (returned in `unknown`); they do not abort the resolve. ResolveRegions
// preserves the input order in `resolved` to keep downstream filtering
// stable.
func ResolveRegions(specs []string, lookup func(chrom string) int) (resolved []ResolvedRegion, unknown []string, err error) {
	for _, spec := range specs {
		reg, perr := ParseRegion(spec)
		if perr != nil {
			return nil, nil, perr
		}
		rid := lookup(reg.Chrom)
		if rid < 0 {
			unknown = append(unknown, reg.Chrom)
			continue
		}
		beg0 := reg.Beg - 1
		if beg0 < 0 {
			beg0 = 0
		}
		end0 := reg.End
		if end0 <= 0 {
			end0 = 1 << 30
		}
		resolved = append(resolved, ResolvedRegion{Region: reg, RefID: rid, Beg0: beg0, End0: end0})
	}
	return resolved, unknown, nil
}

// UnionChunks combines the BAI region-chunks across every ResolvedRegion
// into a single sorted, merged []BAIChunk. It is the entry point used by
// the seek-and-scan region-query path in samtools view.
func UnionChunks(idx *BAIIndex, regions []ResolvedRegion) []BAIChunk {
	var all []BAIChunk
	for _, r := range regions {
		all = append(all, idx.RegionChunks(r.RefID, r.Beg0, r.End0)...)
	}
	if len(all) == 0 {
		return nil
	}
	sort.Slice(all, func(a, b int) bool { return all[a].Beg < all[b].Beg })
	merged := all[:0]
	cur := all[0]
	for i := 1; i < len(all); i++ {
		if all[i].Beg <= cur.End {
			if all[i].End > cur.End {
				cur.End = all[i].End
			}
		} else {
			merged = append(merged, cur)
			cur = all[i]
		}
	}
	merged = append(merged, cur)
	return merged
}
