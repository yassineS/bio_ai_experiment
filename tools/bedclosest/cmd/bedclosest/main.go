// bedclosest finds, for each A interval, the closest B interval (mirrors
// `bedtools closest`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedclosest/pkg/bedclosest"
)

const version = "1.0.0"

const usage = `bedclosest - Find the closest B interval for each A interval

Usage:
  bedclosest -a <fileA.bed> -b <fileB.bed> [options]

Description:
  For each interval in A (sorted), find the closest interval in B (also
  sorted) on the same chromosome and report A's columns + B's columns +
  the signed distance. Distance is 0 when A and B overlap. For tied
  distances the default is to emit one row per tied B (-t all).

  Both inputs MUST be sorted on (chrom, start). bedclosest errors out
  clearly when they are not.

  Note: unlike bedtools closest, the distance column is printed BY DEFAULT.
  Use -d=false to suppress it.

Options:
  -a, --a FILE             Input BED file A (sorted; use '-' for stdin)
  -b, --b FILE             Input BED file B (sorted; use '-' for stdin)
  -o, --output FILE        Output BED file ('-' for stdout, default: stdout)
  -d, --distance           Print signed distance column (default: true)
  -D MODE                  Strandedness of the distance sign: ref (default),
                           a (relative to A's strand), or b (relative to B's).
  -N                       Require strict overlap; non-overlapping B intervals
                           are treated as infinitely far away.
  -t MODE                  Tie-break mode for equally close B's:
                             all   - emit one row per tied B (default)
                             first - emit only the first (in B's input order)
                             last  - emit only the last
  -h, --help               Show this help message and exit
  -v, --version            Show version information and exit

Examples:
  # Closest peak for each gene (both sorted)
  bedclosest -a genes.sorted.bed -b peaks.sorted.bed > out.bed

  # Suppress the distance column
  bedclosest -a a.bed -b b.bed --distance=false > out.bed

  # Only report when A overlaps a B
  bedclosest -a a.bed -b b.bed -N > out.bed

  # Single hit per A (first in B input order on ties)
  bedclosest -a a.bed -b b.bed -t first > out.bed

Format:
  Input: BED format (tab-delimited, minimum 3 columns), sorted on
         (chrom, start).
  Output: A's columns, then B's columns, then signed distance (if -d).
`

func main() {
	fs := flag.CommandLine
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	var aFile, bFile, outputFile string
	// "-a" and "--a" are treated identically by Go's flag package; register once.
	fs.StringVar(&aFile, "a", "", "Input BED file A (sorted)")
	fs.StringVar(&bFile, "b", "", "Input BED file B (sorted)")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output BED file")

	// Distance: default true. Use BoolVar so users can write --distance=false.
	var printDist bool
	fs.BoolVar(&printDist, "d", true, "")
	fs.BoolVar(&printDist, "distance", true, "Print signed distance column (default true)")

	var distMode string
	fs.StringVar(&distMode, "D", "ref", "Distance sign mode: ref|a|b")

	var requireOverlap bool
	fs.BoolVar(&requireOverlap, "N", false, "Require strict overlap")

	var tieMode string
	fs.StringVar(&tieMode, "t", "all", "Tie-break: all|first|last")

	var help, showVersion bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help message")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version information")

	flag.Parse()

	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if showVersion {
		fmt.Printf("bedclosest version %s\n", version)
		os.Exit(0)
	}

	if aFile == "" || bFile == "" {
		fmt.Fprintln(os.Stderr, "Error: -a and -b are required")
		os.Exit(2)
	}
	if aFile == "-" && bFile == "-" {
		fmt.Fprintln(os.Stderr, "Error: -a and -b cannot both be '-' (stdin)")
		os.Exit(2)
	}

	dm, err := parseDistanceMode(distMode)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(2)
	}
	tb, err := parseTieBreak(tieMode)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(2)
	}

	readerA, err := iohelper.OpenReader(aFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening A: %v\n", err)
		os.Exit(1)
	}
	defer readerA.Close()

	readerB, err := iohelper.OpenReader(bFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening B: %v\n", err)
		os.Exit(1)
	}
	defer readerB.Close()

	writer, err := iohelper.OpenWriter(outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		os.Exit(1)
	}
	defer writer.Close()

	opts := bedclosest.Options{
		PrintDistance:  printDist,
		DistanceMode:   dm,
		RequireOverlap: requireOverlap,
		TieBreak:       tb,
	}
	if _, err := bedclosest.Closest(readerA, readerB, writer, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

// parseDistanceMode converts the -D flag string into a bedclosest.DistanceMode.
func parseDistanceMode(s string) (bedclosest.DistanceMode, error) {
	switch s {
	case "ref":
		return bedclosest.DistanceRef, nil
	case "a":
		return bedclosest.DistanceA, nil
	case "b":
		return bedclosest.DistanceB, nil
	default:
		return 0, fmt.Errorf("invalid -D value %q (expected ref|a|b)", s)
	}
}

// parseTieBreak converts the -t flag string into a bedclosest.TieBreak.
func parseTieBreak(s string) (bedclosest.TieBreak, error) {
	switch s {
	case "all":
		return bedclosest.TieAll, nil
	case "first":
		return bedclosest.TieFirst, nil
	case "last":
		return bedclosest.TieLast, nil
	default:
		return 0, fmt.Errorf("invalid -t value %q (expected all|first|last)", s)
	}
}
