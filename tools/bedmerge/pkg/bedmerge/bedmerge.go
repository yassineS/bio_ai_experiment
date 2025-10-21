// Package bedmerge provides functionality to merge overlapping or adjacent BED intervals.
package bedmerge

import (
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// MergeOptions contains options for merging BED intervals.
type MergeOptions struct {
	MaxDistance int  // Maximum distance between intervals to merge (default: 0)
	StrandSpec  bool // Merge only intervals on the same strand
}

// Merge reads BED intervals, sorts them, and merges overlapping/adjacent intervals.
// Returns the number of merged intervals.
func Merge(reader io.Reader, writer io.Writer, opts MergeOptions) (int, error) {
	// Read all intervals
	bedReader := bed.NewReader(reader)
	var intervals []*bed.Record
	
	for {
		record, err := bedReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, fmt.Errorf("error reading BED record: %w", err)
		}
		intervals = append(intervals, record)
	}

	if len(intervals) == 0 {
		return 0, nil
	}

	// Sort intervals by chromosome, then start position
	sort.Slice(intervals, func(i, j int) bool {
		if intervals[i].Chrom != intervals[j].Chrom {
			return intervals[i].Chrom < intervals[j].Chrom
		}
		if intervals[i].ChromStart != intervals[j].ChromStart {
			return intervals[i].ChromStart < intervals[j].ChromStart
		}
		return intervals[i].ChromEnd < intervals[j].ChromEnd
	})

	// Merge overlapping/adjacent intervals
	merged := mergeIntervals(intervals, opts)

	// Write merged intervals
	bedWriter := bed.NewWriter(writer)
	for _, record := range merged {
		if err := bedWriter.Write(record); err != nil {
			return 0, fmt.Errorf("error writing BED record: %w", err)
		}
	}
	if err := bedWriter.Flush(); err != nil {
		return 0, fmt.Errorf("error flushing output: %w", err)
	}

	return len(merged), nil
}

// mergeIntervals performs the actual merging of sorted intervals.
func mergeIntervals(intervals []*bed.Record, opts MergeOptions) []*bed.Record {
	if len(intervals) == 0 {
		return nil
	}

	merged := []*bed.Record{}
	current := &bed.Record{
		Chrom:      intervals[0].Chrom,
		ChromStart: intervals[0].ChromStart,
		ChromEnd:   intervals[0].ChromEnd,
		Strand:     intervals[0].Strand,
	}

	for i := 1; i < len(intervals); i++ {
		interval := intervals[i]

		// Check if we can merge with current
		canMerge := false
		
		// Same chromosome?
		if interval.Chrom == current.Chrom {
			// Check strand if needed
			if !opts.StrandSpec || interval.Strand == current.Strand || current.Strand == "." || interval.Strand == "." {
				// Check if overlapping or within max distance
				if interval.ChromStart <= current.ChromEnd+opts.MaxDistance {
					canMerge = true
				}
			}
		}

		if canMerge {
			// Extend current interval
			if interval.ChromEnd > current.ChromEnd {
				current.ChromEnd = interval.ChromEnd
			}
		} else {
			// Save current and start new interval
			merged = append(merged, current)
			current = &bed.Record{
				Chrom:      interval.Chrom,
				ChromStart: interval.ChromStart,
				ChromEnd:   interval.ChromEnd,
				Strand:     interval.Strand,
			}
		}
	}

	// Add the last interval
	merged = append(merged, current)

	return merged
}

// Stats contains statistics about the merge operation.
type Stats struct {
	InputIntervals  int
	OutputIntervals int
	MergedCount     int
}

// MergeWithStats performs merge and returns detailed statistics.
func MergeWithStats(reader io.Reader, writer io.Writer, opts MergeOptions) (*Stats, error) {
	// Read all intervals
	bedReader := bed.NewReader(reader)
	var intervals []*bed.Record
	
	for {
		record, err := bedReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("error reading BED record: %w", err)
		}
		intervals = append(intervals, record)
	}

	stats := &Stats{
		InputIntervals: len(intervals),
	}

	if len(intervals) == 0 {
		return stats, nil
	}

	// Sort intervals
	sort.Slice(intervals, func(i, j int) bool {
		if intervals[i].Chrom != intervals[j].Chrom {
			return intervals[i].Chrom < intervals[j].Chrom
		}
		if intervals[i].ChromStart != intervals[j].ChromStart {
			return intervals[i].ChromStart < intervals[j].ChromStart
		}
		return intervals[i].ChromEnd < intervals[j].ChromEnd
	})

	// Merge overlapping/adjacent intervals
	merged := mergeIntervals(intervals, opts)
	stats.OutputIntervals = len(merged)
	stats.MergedCount = stats.InputIntervals - stats.OutputIntervals

	// Write merged intervals
	bedWriter := bed.NewWriter(writer)
	for _, record := range merged {
		if err := bedWriter.Write(record); err != nil {
			return nil, fmt.Errorf("error writing BED record: %w", err)
		}
	}
	if err := bedWriter.Flush(); err != nil {
		return nil, fmt.Errorf("error flushing output: %w", err)
	}

	return stats, nil
}
