// Package bedintersect provides functionality to find intersecting intervals between two BED files.
package bedintersect

import (
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// IntersectOptions contains options for finding intersections.
type IntersectOptions struct {
	MinOverlap int     // Minimum overlap required (default: 1bp)
	FractionA  float64 // Minimum fraction of A that must overlap (0.0-1.0)
	FractionB  float64 // Minimum fraction of B that must overlap (0.0-1.0)
	StrandSpec bool    // Only report hits on same strand
	NoOverlap  bool    // Report A entries with no overlap with B
	WriteA     bool    // Write original A entry (default: write intersection)
	WriteB     bool    // Write B entry instead of A
	Count      bool    // For each A, report count of B overlaps
	Reciprocal bool    // Require reciprocal overlap (both -f and -F must be satisfied)
	Distance   bool    // Report distance to nearest B feature
	Closest    bool    // Report closest B feature for each A
	UseTree    bool    // Use interval tree for large B files

	// LeftJoin enables -loj: report every A record, appending the overlapping
	// B record (or a null B placeholder when there are no overlaps).
	LeftJoin bool
	// WriteOverlap enables -wo: report A and B for each overlap, followed by
	// the number of overlapping bases.
	WriteOverlap bool
	// WriteAllOverlap enables -wao: like WriteOverlap, but also report A with a
	// null B and an overlap count of 0 when A has no overlaps.
	WriteAllOverlap bool
	// Split enables -split: treat BED12 blocks as separate intervals when
	// determining overlaps and counting overlapping bases.
	Split bool
}

// Intersect finds intervals in A that overlap with intervals in B.
func Intersect(readerA, readerB io.Reader, writer io.Writer, opts IntersectOptions) (int, error) {
	// The join/overlap output modes (-loj, -wo, -wao, -wa -wb) echo A and B
	// columns verbatim and in B-file order, so they use a separate raw,
	// line-preserving code path.
	if opts.usesJoinMode() {
		return intersectJoin(readerA, readerB, writer, opts)
	}
	// Read all B intervals (database to search against)
	bedReaderB := bed.NewReader(readerB)
	var intervalsB []*bed.Record

	for {
		record, err := bedReaderB.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, fmt.Errorf("error reading B intervals: %w", err)
		}
		intervalsB = append(intervalsB, record)
	}

	// Sort B intervals for efficient searching
	sort.Slice(intervalsB, func(i, j int) bool {
		if intervalsB[i].Chrom != intervalsB[j].Chrom {
			return intervalsB[i].Chrom < intervalsB[j].Chrom
		}
		return intervalsB[i].ChromStart < intervalsB[j].ChromStart
	})

	// Create chromosome index and optionally interval trees
	chromIndex := make(map[string][]*bed.Record)
	chromTrees := make(map[string]*IntervalTree)

	for _, interval := range intervalsB {
		chromIndex[interval.Chrom] = append(chromIndex[interval.Chrom], interval)
	}

	// Build interval trees for each chromosome if UseTree is enabled
	if opts.UseTree {
		for chrom, intervals := range chromIndex {
			chromTrees[chrom] = NewIntervalTree(intervals)
		}
	}

	// Process A intervals
	bedReaderA := bed.NewReader(readerA)
	bedWriter := bed.NewWriter(writer)
	count := 0

	for {
		recordA, err := bedReaderA.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, fmt.Errorf("error reading A intervals: %w", err)
		}

		// Handle distance and closest modes
		if opts.Distance || opts.Closest {
			closest, dist := findClosest(recordA, chromIndex[recordA.Chrom], opts)
			if closest != nil {
				if opts.Distance {
					// Write A interval with distance in name field
					result := &bed.Record{
						Chrom:      recordA.Chrom,
						ChromStart: recordA.ChromStart,
						ChromEnd:   recordA.ChromEnd,
						Name:       fmt.Sprintf("%d", dist),
					}
					if err := bedWriter.Write(result); err != nil {
						return 0, fmt.Errorf("error writing result: %w", err)
					}
					count++
				} else if opts.Closest {
					// Write the closest B interval
					if err := bedWriter.Write(closest); err != nil {
						return 0, fmt.Errorf("error writing result: %w", err)
					}
					count++
				}
			} else {
				// No B intervals on this chromosome
				if opts.Distance {
					result := &bed.Record{
						Chrom:      recordA.Chrom,
						ChromStart: recordA.ChromStart,
						ChromEnd:   recordA.ChromEnd,
						Name:       "-1",
					}
					if err := bedWriter.Write(result); err != nil {
						return 0, fmt.Errorf("error writing result: %w", err)
					}
					count++
				}
			}
			continue
		}

		// Find overlaps using interval tree or linear search
		var overlaps []*Overlap
		if opts.UseTree {
			if tree, ok := chromTrees[recordA.Chrom]; ok {
				candidates := tree.Query(recordA)
				overlaps = findOverlaps(recordA, candidates, opts)
			}
		} else {
			overlaps = findOverlaps(recordA, chromIndex[recordA.Chrom], opts)
		}

		if opts.NoOverlap {
			// Report if no overlaps found
			if len(overlaps) == 0 {
				if err := bedWriter.Write(recordA); err != nil {
					return 0, fmt.Errorf("error writing result: %w", err)
				}
				count++
			}
		} else if opts.Count {
			// Append the overlap count as the final column, preserving all
			// of A's columns (matches `bedtools intersect -c`).
			result := *recordA
			result.ExtraFields = append(append([]string(nil), recordA.ExtraFields...), fmt.Sprintf("%d", len(overlaps)))
			if err := bedWriter.Write(&result); err != nil {
				return 0, fmt.Errorf("error writing result: %w", err)
			}
			count++
		} else {
			// Report each overlap
			for _, overlap := range overlaps {
				var result *bed.Record
				if opts.WriteB {
					result = overlap.B
				} else if opts.WriteA {
					result = recordA
				} else {
					// Write intersection: clip A to the overlap range but
					// preserve A's name/score/strand and any extra columns
					// (matches `bedtools intersect` default behaviour).
					clip := *recordA
					clip.ChromStart = max(recordA.ChromStart, overlap.B.ChromStart)
					clip.ChromEnd = min(recordA.ChromEnd, overlap.B.ChromEnd)
					result = &clip
				}
				if err := bedWriter.Write(result); err != nil {
					return 0, fmt.Errorf("error writing result: %w", err)
				}
				count++
			}
		}
	}

	if err := bedWriter.Flush(); err != nil {
		return 0, fmt.Errorf("error flushing output: %w", err)
	}

	return count, nil
}

// Overlap represents an overlapping interval pair.
type Overlap struct {
	A          *bed.Record
	B          *bed.Record
	OverlapLen int
}

// findOverlaps finds all intervals in B that overlap with A.
func findOverlaps(a *bed.Record, bIntervals []*bed.Record, opts IntersectOptions) []*Overlap {
	var overlaps []*Overlap

	for _, b := range bIntervals {
		// Skip if chromosomes don't match
		if a.Chrom != b.Chrom {
			continue
		}

		// Check strand if required
		if opts.StrandSpec && !sameStrandMatch(a.Strand, b.Strand) {
			continue
		}

		// Calculate overlap
		overlapStart := max(a.ChromStart, b.ChromStart)
		overlapEnd := min(a.ChromEnd, b.ChromEnd)
		overlapLen := overlapEnd - overlapStart

		// Check if there's an overlap
		if overlapLen <= 0 {
			continue
		}

		// Lengths used for the -f/-F fraction tests. With -split these are the
		// summed block lengths and the overlap is the non-redundant block
		// overlap (mirroring BlockMgr::findBlockedOverlaps); otherwise they are
		// the whole-span lengths.
		lenA := a.ChromEnd - a.ChromStart
		lenB := b.ChromEnd - b.ChromStart
		fracOverlap := overlapLen
		if opts.Split {
			blockOverlap, aSum, bSum := splitBlockOverlap(a, b)
			if blockOverlap <= 0 {
				continue
			}
			fracOverlap = blockOverlap
			lenA = aSum
			lenB = bSum
		}

		// Check minimum overlap
		if overlapLen < opts.MinOverlap {
			continue
		}

		// Check fraction of A that overlaps
		fracA := float64(fracOverlap) / float64(lenA)
		if opts.FractionA > 0 && fracA < opts.FractionA {
			continue
		}

		// Check fraction of B that overlaps
		fracB := float64(fracOverlap) / float64(lenB)
		if opts.FractionB > 0 && fracB < opts.FractionB {
			continue
		}

		// Note: Reciprocal mode is satisfied by the two checks above when
		// callers supply both FractionA and FractionB, which is how the CLI
		// wires it (-r combined with -f and -F). No additional check needed.

		overlaps = append(overlaps, &Overlap{
			A:          a,
			B:          b,
			OverlapLen: overlapLen,
		})
	}

	return overlaps
}

// findClosest finds the closest interval in B to A and returns the distance.
// Returns the closest interval and distance (0 if overlapping, positive if upstream/downstream).
func findClosest(a *bed.Record, bIntervals []*bed.Record, opts IntersectOptions) (*bed.Record, int) {
	var closest *bed.Record
	minDist := -1

	for _, b := range bIntervals {
		// Skip if chromosomes don't match
		if a.Chrom != b.Chrom {
			continue
		}

		// Check strand if required
		if opts.StrandSpec && !sameStrandMatch(a.Strand, b.Strand) {
			continue
		}

		// Calculate distance
		var dist int
		if a.ChromEnd <= b.ChromStart {
			// A is upstream of B
			dist = b.ChromStart - a.ChromEnd
		} else if b.ChromEnd <= a.ChromStart {
			// B is upstream of A
			dist = a.ChromStart - b.ChromEnd
		} else {
			// Overlapping
			dist = 0
		}

		// Update closest if this is closer
		if minDist == -1 || dist < minDist {
			minDist = dist
			closest = b
		}
	}

	return closest, minDist
}

// Stats contains statistics about the intersect operation.
type Stats struct {
	IntervalsA     int
	IntervalsB     int
	Overlaps       int
	IntervalsAHit  int
	IntervalsAMiss int
}

// IntersectWithStats performs intersection and returns detailed statistics.
func IntersectWithStats(readerA, readerB io.Reader, writer io.Writer, opts IntersectOptions) (*Stats, error) {
	// Read all B intervals
	bedReaderB := bed.NewReader(readerB)
	var intervalsB []*bed.Record

	for {
		record, err := bedReaderB.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("error reading B intervals: %w", err)
		}
		intervalsB = append(intervalsB, record)
	}

	sort.Slice(intervalsB, func(i, j int) bool {
		if intervalsB[i].Chrom != intervalsB[j].Chrom {
			return intervalsB[i].Chrom < intervalsB[j].Chrom
		}
		return intervalsB[i].ChromStart < intervalsB[j].ChromStart
	})

	chromIndex := make(map[string][]*bed.Record)
	chromTrees := make(map[string]*IntervalTree)

	for _, interval := range intervalsB {
		chromIndex[interval.Chrom] = append(chromIndex[interval.Chrom], interval)
	}

	// Build interval trees for each chromosome if UseTree is enabled
	if opts.UseTree {
		for chrom, intervals := range chromIndex {
			chromTrees[chrom] = NewIntervalTree(intervals)
		}
	}

	stats := &Stats{
		IntervalsB: len(intervalsB),
	}

	// Process A intervals
	bedReaderA := bed.NewReader(readerA)
	bedWriter := bed.NewWriter(writer)

	for {
		recordA, err := bedReaderA.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("error reading A intervals: %w", err)
		}

		stats.IntervalsA++

		// Handle distance and closest modes
		if opts.Distance || opts.Closest {
			closest, dist := findClosest(recordA, chromIndex[recordA.Chrom], opts)
			if closest != nil {
				stats.IntervalsAHit++
				if opts.Distance {
					result := &bed.Record{
						Chrom:      recordA.Chrom,
						ChromStart: recordA.ChromStart,
						ChromEnd:   recordA.ChromEnd,
						Name:       fmt.Sprintf("%d", dist),
					}
					if err := bedWriter.Write(result); err != nil {
						return nil, fmt.Errorf("error writing result: %w", err)
					}
				} else if opts.Closest {
					if err := bedWriter.Write(closest); err != nil {
						return nil, fmt.Errorf("error writing result: %w", err)
					}
				}
			} else {
				stats.IntervalsAMiss++
				if opts.Distance {
					result := &bed.Record{
						Chrom:      recordA.Chrom,
						ChromStart: recordA.ChromStart,
						ChromEnd:   recordA.ChromEnd,
						Name:       "-1",
					}
					if err := bedWriter.Write(result); err != nil {
						return nil, fmt.Errorf("error writing result: %w", err)
					}
				}
			}
			continue
		}

		// Find overlaps using interval tree or linear search
		var overlaps []*Overlap
		if opts.UseTree {
			if tree, ok := chromTrees[recordA.Chrom]; ok {
				candidates := tree.Query(recordA)
				overlaps = findOverlaps(recordA, candidates, opts)
			}
		} else {
			overlaps = findOverlaps(recordA, chromIndex[recordA.Chrom], opts)
		}

		if len(overlaps) > 0 {
			stats.IntervalsAHit++
			stats.Overlaps += len(overlaps)
		} else {
			stats.IntervalsAMiss++
		}

		// Write output (same logic as Intersect)
		if opts.NoOverlap {
			if len(overlaps) == 0 {
				if err := bedWriter.Write(recordA); err != nil {
					return nil, fmt.Errorf("error writing result: %w", err)
				}
			}
		} else if opts.Count {
			result := &bed.Record{
				Chrom:      recordA.Chrom,
				ChromStart: recordA.ChromStart,
				ChromEnd:   recordA.ChromEnd,
				Name:       fmt.Sprintf("%d", len(overlaps)),
			}
			if err := bedWriter.Write(result); err != nil {
				return nil, fmt.Errorf("error writing result: %w", err)
			}
		} else {
			for _, overlap := range overlaps {
				var result *bed.Record
				if opts.WriteB {
					result = overlap.B
				} else if opts.WriteA {
					result = recordA
				} else {
					result = &bed.Record{
						Chrom:      recordA.Chrom,
						ChromStart: max(recordA.ChromStart, overlap.B.ChromStart),
						ChromEnd:   min(recordA.ChromEnd, overlap.B.ChromEnd),
					}
				}
				if err := bedWriter.Write(result); err != nil {
					return nil, fmt.Errorf("error writing result: %w", err)
				}
			}
		}
	}

	if err := bedWriter.Flush(); err != nil {
		return nil, fmt.Errorf("error flushing output: %w", err)
	}

	return stats, nil
}

// splitBlockOverlap returns the non-redundant overlapping bases between A and B
// when treating BED12 records as their constituent blocks, along with the
// summed A-block and B-block lengths used by the fraction tests. Records that
// are not valid BED12 fall back to their whole span (a single block).
func splitBlockOverlap(a, b *bed.Record) (overlap, aSum, bSum int) {
	aBlocks := typedBlocks(a)
	bBlocks := typedBlocks(b)
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

// typedBlocks expands a typed bed.Record into absolute block intervals, or
// returns the whole span when the record has no BED12 blocks.
func typedBlocks(r *bed.Record) []block {
	if r.BlockCount <= 0 || len(r.BlockSizes) != r.BlockCount || len(r.BlockStarts) != r.BlockCount {
		return []block{{r.ChromStart, r.ChromEnd}}
	}
	blks := make([]block, 0, r.BlockCount)
	for i := 0; i < r.BlockCount; i++ {
		s := r.ChromStart + r.BlockStarts[i]
		blks = append(blks, block{s, s + r.BlockSizes[i]})
	}
	return blks
}

// sameStrandMatch reports whether two strand columns count as the same strand
// under `-s`, matching upstream Record::sameChromIntersects. Only "+" and "-"
// are real strands; ".", "*", and a missing column are all UNKNOWN and can
// never satisfy a same-strand requirement, so a hit is reported only when both
// strands are a known value and equal.
func sameStrandMatch(a, b string) bool {
	if a != "+" && a != "-" {
		return false
	}
	return a == b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
