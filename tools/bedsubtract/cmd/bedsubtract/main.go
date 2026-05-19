// bedsubtract subtracts overlapping regions of B intervals from A intervals
// (mirrors `bedtools subtract`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedsubtract/pkg/bedsubtract"
)

const version = "1.0.0"

const usage = `bedsubtract - Subtract B intervals from A intervals

Usage:
  bedsubtract -a <fileA.bed> -b <fileB.bed> [options]

Description:
  For every interval in A, subtract any overlap with intervals in B. The
  remaining segments of A are emitted (in input order). When an interval
  in B punches a hole in the middle of an A interval, the A interval is
  split into multiple output rows. A's columns beyond chrom/start/end are
  preserved through the split (BED3 in -> BED3 out, BED6 in -> BED6 out,
  etc.).

Options:
  -a, --a FILE             Input BED file A (use '-' for stdin)
  -b, --b FILE             Input BED file B (use '-' for stdin)
  -o, --output FILE        Output BED file ('-' for stdout, default: stdout)
  -A                       If any part of A overlaps B, drop the entire A
                           interval (do not split it).
  -N, --min-fraction NUM   Only subtract B from A when the overlap covers at
                           least NUM (0..1) of A.
  -s, --strand             Only subtract B intervals on the same strand as A
                           (BED6+).
  -S                       Only subtract B intervals on the opposite strand
                           from A (BED6+).
  -h, --help               Show this help message and exit
  -v, --version            Show version information and exit

Examples:
  # Subtract peaks from genes
  bedsubtract -a genes.bed -b peaks.bed > genes_minus_peaks.bed

  # Drop any A interval that touches B
  bedsubtract -a genes.bed -b blacklist.bed -A > clean_genes.bed

  # Strand-aware subtraction
  bedsubtract -a a.bed -b b.bed -s > out.bed

  # Stream A from stdin
  cat a.bed | bedsubtract -a - -b b.bed > out.bed

Notes:
  - Coordinates are 0-based, half-open [start, end).
  - Only one of -a and -b may be '-' (stdin) per invocation.
  - With -s or -S, B intervals lacking a strand are skipped.
`

func main() {
	fs := flag.CommandLine
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	var aFile, bFile, outputFile string
	// Both "-a" and "--a" should accept the same value. Go's flag package
	// treats "-a" and "--a" identically, so registering one name covers both.
	fs.StringVar(&aFile, "a", "", "Input BED file A")
	fs.StringVar(&bFile, "b", "", "Input BED file B")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output BED file")

	var removeEntire, sameStrand, oppositeStrand bool
	fs.BoolVar(&removeEntire, "A", false, "Drop A on any overlap")
	cliflag.BoolVar(fs, &sameStrand, "s", "strand", false, "Same-strand only")
	fs.BoolVar(&oppositeStrand, "S", false, "Opposite-strand only")

	var minFraction float64
	cliflag.Float64Var(fs, &minFraction, "N", "min-fraction", 0, "Min overlap fraction of A")

	var help, showVersion bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help message")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version information")

	flag.Parse()

	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if showVersion {
		fmt.Printf("bedsubtract version %s\n", version)
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

	opts := bedsubtract.Options{
		RemoveEntire:   removeEntire,
		MinFraction:    minFraction,
		SameStrand:     sameStrand,
		OppositeStrand: oppositeStrand,
	}
	if _, err := bedsubtract.Subtract(readerA, readerB, writer, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}
