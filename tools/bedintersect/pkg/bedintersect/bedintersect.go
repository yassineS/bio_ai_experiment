// Package bedintersect provides functionality to find intersecting intervals between two BED files.
package bedintersect

import (
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// IntersectOptions contains options for finding intersections.
type IntersectOptions struct {
	MinOverlap   int     // Minimum overlap required (default: 1bp)
	FractionA    float64 // Minimum fraction of A that must overlap (0.0-1.0)
	FractionB    float64 // Minimum fraction of B that must overlap (0.0-1.0)
	StrandSpec   bool    // Only report hits on same strand
	NoOverlap    bool    // Report A entries with no overlap with B
	WriteA       bool    // Write original A entry (default: write intersection)
	WriteB       bool    // Write B entry instead of A
	Count        bool    // For each A, report count of B overlaps
}

// Intersect finds intervals in A that overlap with intervals in B.
func Intersect(readerA, readerB io.Reader, writer io.Writer, opts IntersectOptions) (int, error) {
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

	// Create interval tree index for each chromosome
	chromIndex := make(map[string][]*bed.Record)
	for _, interval := range intervalsB {
		chromIndex[interval.Chrom] = append(chromIndex[interval.Chrom], interval)
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

		// Find overlaps
		overlaps := findOverlaps(recordA, chromIndex[recordA.Chrom], opts)

		if opts.NoOverlap {
			// Report if no overlaps found
			if len(overlaps) == 0 {
				if err := bedWriter.Write(recordA); err != nil {
					return 0, fmt.Errorf("error writing result: %w", err)
				}
				count++
			}
		} else if opts.Count {
			// Report count of overlaps
			result := &bed.Record{
				Chrom:      recordA.Chrom,
				ChromStart: recordA.ChromStart,
				ChromEnd:   recordA.ChromEnd,
				Name:       fmt.Sprintf("%d", len(overlaps)),
			}
			if err := bedWriter.Write(result); err != nil {
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
					// Write intersection
					result = &bed.Record{
						Chrom:      recordA.Chrom,
						ChromStart: max(recordA.ChromStart, overlap.B.ChromStart),
						ChromEnd:   min(recordA.ChromEnd, overlap.B.ChromEnd),
					}
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
		if opts.StrandSpec {
			if a.Strand != "" && b.Strand != "" && a.Strand != b.Strand {
				continue
			}
		}

		// Calculate overlap
		overlapStart := max(a.ChromStart, b.ChromStart)
		overlapEnd := min(a.ChromEnd, b.ChromEnd)
		overlapLen := overlapEnd - overlapStart

		// Check if there's an overlap
		if overlapLen <= 0 {
			continue
		}

		// Check minimum overlap
		if overlapLen < opts.MinOverlap {
			continue
		}

		// Check fraction of A that overlaps
		if opts.FractionA > 0 {
			lenA := a.ChromEnd - a.ChromStart
			fracA := float64(overlapLen) / float64(lenA)
			if fracA < opts.FractionA {
				continue
			}
		}

		// Check fraction of B that overlaps
		if opts.FractionB > 0 {
			lenB := b.ChromEnd - b.ChromStart
			fracB := float64(overlapLen) / float64(lenB)
			if fracB < opts.FractionB {
				continue
			}
		}

		overlaps = append(overlaps, &Overlap{
			A:          a,
			B:          b,
			OverlapLen: overlapLen,
		})
	}

	return overlaps
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
	for _, interval := range intervalsB {
		chromIndex[interval.Chrom] = append(chromIndex[interval.Chrom], interval)
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
		overlaps := findOverlaps(recordA, chromIndex[recordA.Chrom], opts)

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
