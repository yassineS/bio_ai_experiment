// bedmerge merges overlapping or adjacent BED intervals.
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedmerge/pkg/bedmerge"
)

const usage = `bedmerge - Merge overlapping or adjacent BED intervals

Usage:
  bedmerge [options] <input.bed>

Description:
  Reads BED intervals, sorts them by chromosome and position, and merges
  overlapping or adjacent intervals. Intervals are merged if they overlap
  or are within the specified maximum distance.

Options:
  -d, --distance INT    Maximum distance between intervals to merge (default: 0)
  -s, --strand          Merge only intervals on the same strand
  -i, --input FILE      Input BED file (default: stdin)
  -o, --output FILE     Output BED file (default: stdout)
  -S, --stats           Print merge statistics to stderr
  -h, --help            Show this help message

Examples:
  # Merge overlapping intervals
  bedmerge input.bed > merged.bed

  # Merge intervals within 100bp
  bedmerge -d 100 input.bed > merged.bed

  # Merge strand-specific intervals
  bedmerge -s input.bed > merged.bed

  # Show statistics
  bedmerge -S input.bed > merged.bed

  # Use stdin/stdout
  cat input.bed | bedmerge > merged.bed

Format:
  Input: BED format (tab-delimited, minimum 3 columns: chrom, start, end)
  Output: BED3 format (chrom, start, end)

Notes:
  - Coordinates are 0-based, half-open [start, end)
  - Input does not need to be sorted
  - Adjacent intervals (touching but not overlapping) are merged by default
  - Use -d to merge nearby intervals within specified distance
`

func main() {
	// Define flags
	distance := flag.Int("d", 0, "Maximum distance between intervals to merge")
	flag.IntVar(distance, "distance", 0, "Maximum distance between intervals to merge")
	
	strandSpec := flag.Bool("s", false, "Merge only intervals on the same strand")
	flag.BoolVar(strandSpec, "strand", false, "Merge only intervals on the same strand")
	
	inputFile := flag.String("i", "", "Input BED file (default: stdin)")
	flag.StringVar(inputFile, "input", "", "Input BED file (default: stdin)")
	
	outputFile := flag.String("o", "", "Output BED file (default: stdout)")
	flag.StringVar(outputFile, "output", "", "Output BED file (default: stdout)")
	
	showStats := flag.Bool("S", false, "Print merge statistics to stderr")
	flag.BoolVar(showStats, "stats", false, "Print merge statistics to stderr")
	
	help := flag.Bool("h", false, "Show help message")
	flag.BoolVar(help, "help", false, "Show help message")

	flag.Parse()

	if *help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}

	// Determine input file
	input := *inputFile
	if input == "" && flag.NArg() > 0 {
		input = flag.Arg(0)
	}

	// Open input
	inputReader, err := iohelper.OpenReader(input)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input: %v\n", err)
		os.Exit(1)
	}
	defer inputReader.Close()

	// Open output
	outputWriter, err := iohelper.OpenWriter(*outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		os.Exit(1)
	}
	defer outputWriter.Close()

	// Set merge options
	opts := bedmerge.MergeOptions{
		MaxDistance: *distance,
		StrandSpec:  *strandSpec,
	}

	// Perform merge
	if *showStats {
		stats, err := bedmerge.MergeWithStats(inputReader, outputWriter, opts)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error merging intervals: %v\n", err)
			os.Exit(1)
		}
		fmt.Fprintf(os.Stderr, "Input intervals:  %d\n", stats.InputIntervals)
		fmt.Fprintf(os.Stderr, "Output intervals: %d\n", stats.OutputIntervals)
		fmt.Fprintf(os.Stderr, "Merged: %d intervals\n", stats.MergedCount)
	} else {
		count, err := bedmerge.Merge(inputReader, outputWriter, opts)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error merging intervals: %v\n", err)
			os.Exit(1)
		}
		_ = count // Silently succeed
	}
}
