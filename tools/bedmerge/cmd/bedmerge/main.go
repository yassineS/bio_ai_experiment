// bedmerge merges overlapping or adjacent BED intervals.
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedmerge/pkg/bedmerge"
)

const version = "1.0.0"

const usage = `bedmerge - Merge overlapping or adjacent BED intervals

Usage:
  bedmerge [options] <input.bed>

Description:
  Reads BED intervals, sorts them by chromosome and position, and merges
  overlapping or adjacent intervals. Intervals are merged if they overlap
  or are within the specified maximum distance.

Options:
  -d, --distance INT     Maximum distance between intervals to merge (default: 0)
  -s, --strand           Merge only intervals on the same strand
  -S <+|->               Merge only intervals on the given strand
  -i, --input FILE       Input BED file (default: stdin)
      --output FILE      Output BED file (default: stdout)
      --stats            Print merge statistics to stderr
      --count            Output count of merged intervals as name field
  -g, --bedgraph         Input/output in bedGraph format (chrom, start, end, score)
  -c, --columns LIST     Comma-separated 1-based input columns to aggregate
                         (bedtools merge -c style); requires -o
  -o, --operations LIST  Comma-separated operations, one per -c column or a
                         single op applied to all. Supported: sum, min, max,
                         mean, median, count, count_distinct, distinct,
                         collapse, first, last, mode, antimode
      --delim CHAR       Delimiter joining collapse/distinct/freq list values
                         (default ",")
      --streaming        Use streaming mode for very large files
  -h, --help             Show this help message
  -v, --version          Show version information and exit

Examples:
  # Merge overlapping intervals
  bedmerge input.bed > merged.bed

  # Merge intervals within 100bp
  bedmerge -d 100 input.bed > merged.bed

  # Merge strand-specific intervals
  bedmerge -s input.bed > merged.bed

  # Show statistics
  bedmerge -S input.bed > merged.bed

  # Output with merge count
  bedmerge --count input.bed > merged.bed

  # Aggregate columns over merged groups (bedtools merge -c/-o style)
  bedmerge -c 4,5 -o distinct,sum input.bed > merged.bed
  bedmerge -c 5,6 -o mean input.bed > merged.bed

  # Merge bedGraph files
  bedmerge -g input.bedgraph > merged.bedgraph

  # Use streaming mode for large files
  bedmerge --streaming large.bed > merged.bed

  # Use stdin/stdout
  cat input.bed | bedmerge > merged.bed

Format:
  Input: BED format (tab-delimited, minimum 3 columns: chrom, start, end)
         or bedGraph format with -g flag (4 columns: chrom, start, end, score)
  Output: BED3 format (chrom, start, end) by default; with -c/-o the output is
          chrom, start, end followed by one aggregated value per requested column

Notes:
  - Coordinates are 0-based, half-open [start, end)
  - Input does not need to be sorted
  - Adjacent intervals (touching but not overlapping) are merged by default
  - Use -d to merge nearby intervals within specified distance
  - Streaming mode processes by chromosome and is more memory efficient
`

func main() {
	fs := flag.CommandLine

	var distance int
	cliflag.IntVar(fs, &distance, "d", "distance", 0, "Maximum distance between intervals to merge")

	var strandSpec bool
	cliflag.BoolVar(fs, &strandSpec, "s", "strand", false, "Merge only intervals on the same strand")

	var strandFilter string
	fs.StringVar(&strandFilter, "S", "", "Merge only intervals on the given strand (+ or -)")

	var inputFile string
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input BED file (default: stdin)")

	var outputFile string
	cliflag.StringVar(fs, &outputFile, "", "output", "", "Output BED file (default: stdout)")

	var showStats bool
	cliflag.BoolVar(fs, &showStats, "", "stats", false, "Print merge statistics to stderr")

	var showCount bool
	cliflag.BoolVar(fs, &showCount, "", "count", false, "Output count of merged intervals as name field")

	var bedGraph bool
	cliflag.BoolVar(fs, &bedGraph, "g", "bedgraph", false, "Input/output in bedGraph format")

	var columns string
	cliflag.StringVar(fs, &columns, "c", "columns", "", "Comma-separated 1-based input columns to aggregate")

	var operations string
	cliflag.StringVar(fs, &operations, "o", "operations", "", "Comma-separated operations, one per -c column or one applied to all")

	var delim string
	fs.StringVar(&delim, "delim", ",", "Delimiter for collapse/distinct/freq list joins")

	var streaming bool
	cliflag.BoolVar(fs, &streaming, "", "streaming", false, "Use streaming mode for very large files")

	var help bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help message")

	var showVersion bool
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version information")

	flag.Parse()

	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}

	if showVersion {
		fmt.Printf("bedmerge version %s\n", version)
		os.Exit(0)
	}

	columnOps, err := bedmerge.ParseColumnOps(columns, operations)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
	if columnOps != nil {
		columnOps.Delim = delim
	}

	// Determine input file
	input := inputFile
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
	outputWriter, err := iohelper.OpenWriter(outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		os.Exit(1)
	}
	defer outputWriter.Close()

	// Set merge options
	opts := bedmerge.MergeOptions{
		MaxDistance:  distance,
		StrandSpec:   strandSpec,
		StrandFilter: strandFilter,
		Streaming:    streaming,
		ColumnOps:    columnOps,
		OutputFields: bedmerge.OutputFields{
			Count:    showCount,
			BedGraph: bedGraph,
		},
	}

	// Perform merge
	if showStats {
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
