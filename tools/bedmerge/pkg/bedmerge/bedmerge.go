// Package bedmerge provides functionality to merge overlapping or adjacent BED intervals.
package bedmerge

import (
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// OutputFields specifies which fields to include in the output.
type OutputFields struct {
	Name     bool // Include name field
	Score    bool // Include score field
	Strand   bool // Include strand field
	Count    bool // Include count of merged intervals as name field
	BedGraph bool // Output in bedGraph format (chrom, start, end, score)
}

// MergeOptions contains options for merging BED intervals.
type MergeOptions struct {
	MaxDistance  int          // Maximum distance between intervals to merge (default: 0)
	StrandSpec   bool         // Merge only intervals on the same strand
	OutputFields OutputFields // Fields to include in output
	Streaming    bool         // Use streaming mode for large files
	// ColumnOps, when non-nil, requests bedtools-merge-style aggregation of
	// input columns (the -c/-o options). When set, output columns are
	// chrom, start, end followed by one aggregated value per requested column.
	ColumnOps *ColumnOps
}

// Merge reads BED intervals, sorts them, and merges overlapping/adjacent intervals.
// Returns the number of merged intervals.
func Merge(reader io.Reader, writer io.Writer, opts MergeOptions) (int, error) {
	// Column-aggregation mode (bedtools merge -c/-o style).
	if opts.ColumnOps != nil {
		return mergeWithColumnOps(reader, writer, opts)
	}

	// Use streaming mode if requested
	if opts.Streaming {
		return streamingMerge(reader, writer, opts)
	}

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

		// If bedGraph input mode, parse the name field as score
		if opts.OutputFields.BedGraph && record.Name != "" && record.Score == 0 {
			// Try to parse name as score
			var score int
			if _, err := fmt.Sscanf(record.Name, "%d", &score); err == nil {
				record.Score = score
				record.Name = ""
			}
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
	if err := writeIntervals(writer, merged, opts); err != nil {
		return 0, err
	}

	return len(merged), nil
}

// mergedInterval represents a merged interval with metadata.
type mergedInterval struct {
	*bed.Record
	count int // Number of original intervals merged into this one
}

// mergeIntervals performs the actual merging of sorted intervals.
func mergeIntervals(intervals []*bed.Record, opts MergeOptions) []mergedInterval {
	if len(intervals) == 0 {
		return nil
	}

	merged := []mergedInterval{}
	current := mergedInterval{
		Record: &bed.Record{
			Chrom:      intervals[0].Chrom,
			ChromStart: intervals[0].ChromStart,
			ChromEnd:   intervals[0].ChromEnd,
			Strand:     intervals[0].Strand,
			Name:       intervals[0].Name,
			Score:      intervals[0].Score,
		},
		count: 1,
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
			current.count++

			// For bedGraph, we might want to average or sum scores
			// For now, we keep the first score
		} else {
			// Save current and start new interval
			merged = append(merged, current)
			current = mergedInterval{
				Record: &bed.Record{
					Chrom:      interval.Chrom,
					ChromStart: interval.ChromStart,
					ChromEnd:   interval.ChromEnd,
					Strand:     interval.Strand,
					Name:       interval.Name,
					Score:      interval.Score,
				},
				count: 1,
			}
		}
	}

	// Add the last interval
	merged = append(merged, current)

	return merged
}

// writeIntervals writes merged intervals according to output options.
func writeIntervals(writer io.Writer, merged []mergedInterval, opts MergeOptions) error {
	// Handle bedGraph format separately (requires score without name)
	if opts.OutputFields.BedGraph {
		return writeBedGraph(writer, merged)
	}

	bedWriter := bed.NewWriter(writer)

	for _, m := range merged {
		record := m.Record

		// Handle count field
		if opts.OutputFields.Count {
			record.Name = fmt.Sprintf("%d", m.count)
		} else if !opts.OutputFields.Name {
			record.Name = ""
		}

		if !opts.OutputFields.Score {
			record.Score = 0
		}

		if !opts.OutputFields.Strand {
			record.Strand = ""
		}

		if err := bedWriter.Write(record); err != nil {
			return fmt.Errorf("error writing BED record: %w", err)
		}
	}

	if err := bedWriter.Flush(); err != nil {
		return fmt.Errorf("error flushing output: %w", err)
	}

	return nil
}

// writeBedGraph writes intervals in bedGraph format (chrom, start, end, score).
func writeBedGraph(writer io.Writer, merged []mergedInterval) error {
	for _, m := range merged {
		line := fmt.Sprintf("%s\t%d\t%d\t%d\n",
			m.Chrom, m.ChromStart, m.ChromEnd, m.Score)
		if _, err := writer.Write([]byte(line)); err != nil {
			return fmt.Errorf("error writing bedGraph record: %w", err)
		}
	}
	return nil
}

// streamingMerge performs streaming merge for very large files.
// This uses a sliding window approach to avoid loading all intervals into memory.
func streamingMerge(reader io.Reader, writer io.Writer, opts MergeOptions) (int, error) {
	bedReader := bed.NewReader(reader)

	// For streaming mode, we need to read intervals in chunks per chromosome
	// We'll collect all intervals for one chromosome, merge them, write output,
	// then move to the next chromosome. This assumes input is roughly sorted by chromosome.

	var currentChrom string
	var chromIntervals []*bed.Record
	outputCount := 0

	flushChrom := func() error {
		if len(chromIntervals) == 0 {
			return nil
		}

		// Sort intervals for this chromosome
		sort.Slice(chromIntervals, func(i, j int) bool {
			if chromIntervals[i].ChromStart != chromIntervals[j].ChromStart {
				return chromIntervals[i].ChromStart < chromIntervals[j].ChromStart
			}
			return chromIntervals[i].ChromEnd < chromIntervals[j].ChromEnd
		})

		// Merge and write
		merged := mergeIntervals(chromIntervals, opts)
		if err := writeIntervals(writer, merged, opts); err != nil {
			return err
		}

		outputCount += len(merged)
		chromIntervals = nil
		return nil
	}

	for {
		record, err := bedReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, fmt.Errorf("error reading BED record: %w", err)
		}

		// If bedGraph input mode, parse the name field as score
		if opts.OutputFields.BedGraph && record.Name != "" && record.Score == 0 {
			// Try to parse name as score
			var score int
			if _, err := fmt.Sscanf(record.Name, "%d", &score); err == nil {
				record.Score = score
				record.Name = ""
			}
		}

		// If chromosome changed, flush previous chromosome
		if record.Chrom != currentChrom && currentChrom != "" {
			if err := flushChrom(); err != nil {
				return 0, err
			}
		}

		currentChrom = record.Chrom
		chromIntervals = append(chromIntervals, record)
	}

	// Flush last chromosome
	if err := flushChrom(); err != nil {
		return 0, err
	}

	return outputCount, nil
}

// Stats contains statistics about the merge operation.
type Stats struct {
	InputIntervals  int
	OutputIntervals int
	MergedCount     int
}

// MergeWithStats performs merge and returns detailed statistics.
func MergeWithStats(reader io.Reader, writer io.Writer, opts MergeOptions) (*Stats, error) {
	// Column-aggregation mode: report only the output count.
	if opts.ColumnOps != nil {
		count, err := mergeWithColumnOps(reader, writer, opts)
		if err != nil {
			return nil, err
		}
		return &Stats{OutputIntervals: count}, nil
	}

	// Streaming mode doesn't support detailed stats tracking currently
	if opts.Streaming {
		count, err := streamingMerge(reader, writer, opts)
		if err != nil {
			return nil, err
		}
		return &Stats{
			OutputIntervals: count,
		}, nil
	}

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

		// If bedGraph input mode, parse the name field as score
		if opts.OutputFields.BedGraph && record.Name != "" && record.Score == 0 {
			// Try to parse name as score
			var score int
			if _, err := fmt.Sscanf(record.Name, "%d", &score); err == nil {
				record.Score = score
				record.Name = ""
			}
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
	if err := writeIntervals(writer, merged, opts); err != nil {
		return nil, err
	}

	return stats, nil
}
