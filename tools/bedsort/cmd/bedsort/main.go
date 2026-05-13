// bedsort sorts BED intervals using various sort modes (mirrors `bedtools sort`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedsort/pkg/bedsort"
)

const version = "1.0.0"

const usage = `bedsort - Sort BED intervals

Usage:
  bedsort [options] [<input.bed>]

Description:
  Reads a BED file (or stdin with '-') and writes the records out in sorted
  order. The default sort is by chromosome (lexicographic), then by chromStart
  ascending, then by chromEnd ascending. The full set of input columns is
  preserved through the sort (BED3, BED6, BED12, or any number of extra
  columns).

Options:
  -i, --input FILE         Input BED file ('-' for stdin, default: stdin).
                           Plain or .gz; transparently decompressed.
  -o, --output FILE        Output BED file ('-' for stdout, default: stdout)
      --sizeA              Sort by interval size ascending
      --sizeD              Sort by interval size descending
      --chrThenSizeA       Sort by chromosome, then by interval size ascending
      --chrThenSizeD       Sort by chromosome, then by interval size descending
      --chrThenScoreA      Sort by chromosome, then by score (col 5) ascending
      --chrThenScoreD      Sort by chromosome, then by score (col 5) descending
      --faidx FILE         Order chromosomes by their appearance in FILE
                           (a .fai or chrom-sizes file). Chromosomes that are
                           not listed in FILE sort after the listed ones, in
                           lexicographic order.
  -g, --genome FILE        Alias for --faidx
  -h, --help               Show this help message and exit
  -v, --version            Show version information and exit

Examples:
  # Default sort (chrom + start + end)
  bedsort -i input.bed -o sorted.bed

  # Read from stdin, write to stdout
  cat input.bed | bedsort > sorted.bed

  # Sort by interval size, descending
  bedsort --sizeD input.bed

  # Order chromosomes by the order they appear in a .fai file
  bedsort --faidx hg38.fa.fai input.bed
`

func main() {
	fs := flag.CommandLine
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	var inputFile string
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input BED file ('-' for stdin)")

	var outputFile string
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output BED file ('-' for stdout)")

	var sizeA, sizeD, chrSizeA, chrSizeD, chrScoreA, chrScoreD bool
	cliflag.BoolVar(fs, &sizeA, "", "sizeA", false, "Sort by interval size ascending")
	cliflag.BoolVar(fs, &sizeD, "", "sizeD", false, "Sort by interval size descending")
	cliflag.BoolVar(fs, &chrSizeA, "", "chrThenSizeA", false, "Sort by chromosome then size ascending")
	cliflag.BoolVar(fs, &chrSizeD, "", "chrThenSizeD", false, "Sort by chromosome then size descending")
	cliflag.BoolVar(fs, &chrScoreA, "", "chrThenScoreA", false, "Sort by chromosome then score ascending")
	cliflag.BoolVar(fs, &chrScoreD, "", "chrThenScoreD", false, "Sort by chromosome then score descending")

	var faidxFile string
	cliflag.StringVar(fs, &faidxFile, "", "faidx", "", "Chromosome order file (.fai or chrom-sizes)")
	cliflag.StringVar(fs, &faidxFile, "g", "genome", "", "Alias for --faidx")

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
		fmt.Printf("bedsort version %s\n", version)
		os.Exit(0)
	}

	// Resolve sort mode (at most one mode flag may be set).
	mode, err := resolveMode(sizeA, sizeD, chrSizeA, chrSizeD, chrScoreA, chrScoreD)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(2)
	}

	// Positional input file (e.g. `bedsort sorted.bed`) overrides default stdin
	// when -i wasn't supplied. This matches the bedtools convention.
	input := inputFile
	if input == "" && flag.NArg() > 0 {
		input = flag.Arg(0)
	}

	// Optional chromosome ordering.
	var chromOrder []string
	if faidxFile != "" {
		fr, ferr := iohelper.OpenReader(faidxFile)
		if ferr != nil {
			fmt.Fprintf(os.Stderr, "Error opening %s: %v\n", faidxFile, ferr)
			os.Exit(1)
		}
		chromOrder, ferr = bedsort.ReadFaidx(fr)
		fr.Close()
		if ferr != nil {
			fmt.Fprintf(os.Stderr, "Error reading %s: %v\n", faidxFile, ferr)
			os.Exit(1)
		}
	}

	in, err := iohelper.OpenReader(input)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input: %v\n", err)
		os.Exit(1)
	}
	defer in.Close()

	out, err := iohelper.OpenWriter(outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if err := bedsort.Run(in, out, bedsort.Options{Mode: mode, ChromOrder: chromOrder}); err != nil {
		fmt.Fprintf(os.Stderr, "Error sorting: %v\n", err)
		os.Exit(1)
	}
}

// resolveMode picks the single sort mode requested on the command line.
//
// It is an error to set more than one of the mutually exclusive --sizeA /
// --sizeD / --chrThenSize{A,D} / --chrThenScore{A,D} flags.
func resolveMode(sizeA, sizeD, chrSizeA, chrSizeD, chrScoreA, chrScoreD bool) (bedsort.SortMode, error) {
	flags := []struct {
		set  bool
		mode bedsort.SortMode
		name string
	}{
		{sizeA, bedsort.ModeSizeA, "--sizeA"},
		{sizeD, bedsort.ModeSizeD, "--sizeD"},
		{chrSizeA, bedsort.ModeChrThenSizeA, "--chrThenSizeA"},
		{chrSizeD, bedsort.ModeChrThenSizeD, "--chrThenSizeD"},
		{chrScoreA, bedsort.ModeChrThenScoreA, "--chrThenScoreA"},
		{chrScoreD, bedsort.ModeChrThenScoreD, "--chrThenScoreD"},
	}
	var chosen bedsort.SortMode = bedsort.ModeChrom
	var chosenName string
	for _, f := range flags {
		if !f.set {
			continue
		}
		if chosenName != "" {
			return 0, fmt.Errorf("cannot combine %s with %s", chosenName, f.name)
		}
		chosen = f.mode
		chosenName = f.name
	}
	return chosen, nil
}
