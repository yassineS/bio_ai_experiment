// bedintersect finds intersecting intervals between two BED files.
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedintersect/pkg/bedintersect"
)

const usage = `bedintersect - Find intersecting intervals between two BED files

Usage:
  bedintersect [options] -a <fileA.bed> -b <fileB.bed>

Description:
  Report intervals in A that overlap intervals in B. By default, reports
  the intersection (overlapping portion) of each pair. Various options
  control what gets reported and the minimum overlap required.

Options:
  -a, --input-a FILE    Input BED file A (required)
  -b, --input-b FILE    Input BED file B (required)
  -o, --output FILE     Output file (default: stdout)
  -m, --min-overlap INT Minimum overlap required (default: 1)
  -f, --fraction-a NUM  Minimum fraction of A that must overlap (0.0-1.0)
  -F, --fraction-b NUM  Minimum fraction of B that must overlap (0.0-1.0)
  -r, --reciprocal      Require reciprocal overlap (both -f and -F)
  -s, --strand          Only report hits on same strand
  -v, --invert          Report A entries with NO overlap with B
  -wa, --write-a        Write original A entry (default: write intersection)
  -wb, --write-b        Write B entry instead of A
  -c, --count           Report count of B overlaps for each A
  -d, --distance        Report distance to nearest B feature
  -k, --closest         Output closest B feature for each A
  -t, --tree            Use interval tree for large B files
  -S, --stats           Print statistics to stderr
  -h, --help            Show this help message

Examples:
  # Find overlapping regions
  bedintersect -a genes.bed -b peaks.bed

  # Require 50bp minimum overlap
  bedintersect -a genes.bed -b peaks.bed -m 50

  # Require 80% of gene overlaps peak
  bedintersect -a genes.bed -b peaks.bed -f 0.8

  # Require reciprocal 50% overlap
  bedintersect -a genes.bed -b peaks.bed -f 0.5 -F 0.5 -r

  # Report genes that don't overlap peaks
  bedintersect -a genes.bed -b peaks.bed -v

  # Count overlaps per gene
  bedintersect -a genes.bed -b peaks.bed -c

  # Get original A entries
  bedintersect -a genes.bed -b peaks.bed -wa

  # Get B entries that overlap A
  bedintersect -a genes.bed -b peaks.bed -wb

  # Report distance to nearest peak
  bedintersect -a genes.bed -b peaks.bed -d

  # Report closest peak for each gene
  bedintersect -a genes.bed -b peaks.bed -k

  # Use interval tree for large files
  bedintersect -a genes.bed -b large_features.bed -t

  # Strand-specific intersection
  bedintersect -a genes.bed -b peaks.bed -s

Format:
  Input: BED format (tab-delimited, minimum 3 columns)
  Output: Depends on options (-wa, -wb, -c, -d, -k)
  Default: Intersection coordinates (chrom, start, end)

Notes:
  - Coordinates are 0-based, half-open [start, end)
  - Files do not need to be sorted
  - Multiple B hits generate multiple output lines (unless -c)
`

func main() {
	// Define flags
	inputA := flag.String("a", "", "Input BED file A (required)")
	flag.StringVar(inputA, "input-a", "", "Input BED file A (required)")

	inputB := flag.String("b", "", "Input BED file B (required)")
	flag.StringVar(inputB, "input-b", "", "Input BED file B (required)")

	output := flag.String("o", "", "Output file (default: stdout)")
	flag.StringVar(output, "output", "", "Output file (default: stdout)")

	minOverlap := flag.Int("m", 1, "Minimum overlap required")
	flag.IntVar(minOverlap, "min-overlap", 1, "Minimum overlap required")

	fractionA := flag.Float64("f", 0.0, "Minimum fraction of A that must overlap")
	flag.Float64Var(fractionA, "fraction-a", 0.0, "Minimum fraction of A that must overlap")

	fractionB := flag.Float64("F", 0.0, "Minimum fraction of B that must overlap")
	flag.Float64Var(fractionB, "fraction-b", 0.0, "Minimum fraction of B that must overlap")

	strandSpec := flag.Bool("s", false, "Only report hits on same strand")
	flag.BoolVar(strandSpec, "strand", false, "Only report hits on same strand")

	invert := flag.Bool("v", false, "Report A entries with NO overlap with B")
	flag.BoolVar(invert, "invert", false, "Report A entries with NO overlap with B")

	writeA := flag.Bool("wa", false, "Write original A entry")
	flag.BoolVar(writeA, "write-a", false, "Write original A entry")

	writeB := flag.Bool("wb", false, "Write B entry instead of A")
	flag.BoolVar(writeB, "write-b", false, "Write B entry instead of A")

	writeOverlap := flag.Bool("wo", false, "Write A and B plus the overlap length per hit")
	flag.BoolVar(writeOverlap, "write-overlap", false, "Write A and B plus the overlap length per hit")

	writeAllOverlap := flag.Bool("wao", false, "Like -wo, but emit every A (with B columns = '.', '-1' for misses)")
	flag.BoolVar(writeAllOverlap, "write-all-overlap", false, "Like -wo, but emit every A (with B columns = '.', '-1' for misses)")

	count := flag.Bool("c", false, "Report count of B overlaps for each A")
	flag.BoolVar(count, "count", false, "Report count of B overlaps for each A")

	showStats := flag.Bool("S", false, "Print statistics to stderr")
	flag.BoolVar(showStats, "stats", false, "Print statistics to stderr")

	reciprocal := flag.Bool("r", false, "Require reciprocal overlap (both -f and -F)")
	flag.BoolVar(reciprocal, "reciprocal", false, "Require reciprocal overlap (both -f and -F)")

	distance := flag.Bool("d", false, "Report distance to nearest B feature")
	flag.BoolVar(distance, "distance", false, "Report distance to nearest B feature")

	closest := flag.Bool("k", false, "Output closest B feature for each A")
	flag.BoolVar(closest, "closest", false, "Output closest B feature for each A")

	useTree := flag.Bool("t", false, "Use interval tree for large B files")
	flag.BoolVar(useTree, "tree", false, "Use interval tree for large B files")

	help := flag.Bool("h", false, "Show help message")
	flag.BoolVar(help, "help", false, "Show help message")

	flag.Parse()

	if *help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}

	// Validate required inputs
	if *inputA == "" || *inputB == "" {
		fmt.Fprintln(os.Stderr, "Error: Both -a and -b are required")
		fmt.Fprint(os.Stderr, "\nUsage: bedintersect -a fileA.bed -b fileB.bed [options]\n")
		fmt.Fprintln(os.Stderr, "Use -h for help")
		os.Exit(1)
	}

	// Open inputs
	readerA, err := iohelper.OpenReader(*inputA)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening file A: %v\n", err)
		os.Exit(1)
	}
	defer readerA.Close()

	readerB, err := iohelper.OpenReader(*inputB)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening file B: %v\n", err)
		os.Exit(1)
	}
	defer readerB.Close()

	// Open output
	writer, err := iohelper.OpenWriter(*output)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		os.Exit(1)
	}
	defer writer.Close()

	// Set options
	opts := bedintersect.IntersectOptions{
		MinOverlap:      *minOverlap,
		FractionA:       *fractionA,
		FractionB:       *fractionB,
		StrandSpec:      *strandSpec,
		NoOverlap:       *invert,
		WriteA:          *writeA,
		WriteB:          *writeB,
		WriteOverlap:    *writeOverlap,
		WriteAllOverlap: *writeAllOverlap,
		Count:           *count,
		Reciprocal:      *reciprocal,
		Distance:        *distance,
		Closest:         *closest,
		UseTree:         *useTree,
	}

	// Perform intersection
	if *showStats {
		stats, err := bedintersect.IntersectWithStats(readerA, readerB, writer, opts)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error finding intersections: %v\n", err)
			os.Exit(1)
		}
		fmt.Fprintf(os.Stderr, "Intervals in A: %d\n", stats.IntervalsA)
		fmt.Fprintf(os.Stderr, "Intervals in B: %d\n", stats.IntervalsB)
		fmt.Fprintf(os.Stderr, "A intervals with hits: %d\n", stats.IntervalsAHit)
		fmt.Fprintf(os.Stderr, "A intervals with no hits: %d\n", stats.IntervalsAMiss)
		fmt.Fprintf(os.Stderr, "Total overlaps: %d\n", stats.Overlaps)
	} else {
		_, err := bedintersect.Intersect(readerA, readerB, writer, opts)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error finding intersections: %v\n", err)
			os.Exit(1)
		}
	}
}
