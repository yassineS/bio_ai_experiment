// bedunionbedg combines multiple BEDGRAPH files into a single matrix, emitting
// the value of each input over every interval boundary (mirrors
// `bedtools unionbedg`, aka unionBedGraphs).
package main

import (
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedunionbedg/pkg/bedunionbedg"
)

const version = "1.0.0"

const usage = `bedunionbedg - Combine multiple BedGraph files into one matrix

Usage:
  bedunionbedg [options] -i FILE1 FILE2 .. FILEn
  Assumes each BedGraph file is sorted by chrom/start and that the intervals in
  each are non-overlapping.

Options:
  -i FILE1 FILE2 ..  Input BedGraph files (two or more). Required.
  -o, --output FILE  Output file ('-' for stdout, default: stdout)
  -header            Print a header line (chrom/start/end + names of each file).
  -names NAME ..     A list of names (one per file) printed in the header.
  -g FILE            Genome (chrom-sizes) file, used to calculate empty regions.
  -empty             Report empty regions (requires -g).
  -filler TEXT       Text to use for intervals with no value (default '0').
  -h, --help         Show this help message and exit
  -v, --version      Show version information and exit
`

func main() {
	args := os.Args[1:]

	var (
		inputFiles  []string
		names       []string
		genomeFile  string
		outputFile  string
		filler      = "0"
		printHeader bool
		printEmpty  bool
		help        bool
		showVersion bool
	)

	// Manual parse: -i and -names are variadic (consume args until the next
	// token starting with '-'), matching upstream's hand-rolled parser.
	for i := 0; i < len(args); i++ {
		a := args[i]
		switch a {
		case "-h", "--help":
			help = true
		case "-v", "--version":
			showVersion = true
		case "-header":
			printHeader = true
		case "-empty":
			printEmpty = true
		case "-g", "--genome":
			if i+1 < len(args) {
				i++
				genomeFile = args[i]
			}
		case "-o", "--output":
			if i+1 < len(args) {
				i++
				outputFile = args[i]
			}
		case "-filler":
			if i+1 < len(args) {
				i++
				filler = args[i]
			}
		case "-i", "--input":
			for i+1 < len(args) && !startsWithDash(args[i+1]) {
				i++
				inputFiles = append(inputFiles, args[i])
			}
		case "-names":
			for i+1 < len(args) && !startsWithDash(args[i+1]) {
				i++
				names = append(names, args[i])
			}
		default:
			fmt.Fprintf(os.Stderr, "Error: unrecognized parameter: %s\n", a)
			fmt.Fprint(os.Stderr, usage)
			os.Exit(1)
		}
	}

	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if showVersion {
		fmt.Printf("bedunionbedg version %s\n", version)
		os.Exit(0)
	}

	if len(inputFiles) == 0 {
		fmt.Fprintln(os.Stderr, "Error: missing BedGraph file names (-i) to combine.")
		os.Exit(1)
	}
	if len(inputFiles) == 1 {
		fmt.Fprintln(os.Stderr, "Error: Only a single BedGraph file was specified. Nothing to combine, exiting.")
		os.Exit(1)
	}
	if printEmpty && genomeFile == "" {
		fmt.Fprintln(os.Stderr, "Error: when using -empty, the genome sizes file (-g) must be specified using '-g FILE'.")
		os.Exit(1)
	}
	if len(names) > 0 && len(names) != len(inputFiles) {
		fmt.Fprintln(os.Stderr, "Error: The number of file titles (-names) does not match the number of files (-i).")
		os.Exit(1)
	}

	opts := bedunionbedg.Options{
		PrintHeader: printHeader,
		Names:       names,
		PrintEmpty:  printEmpty,
		Filler:      filler,
	}

	if genomeFile != "" {
		gr, err := iohelper.OpenReader(genomeFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening genome file: %v\n", err)
			os.Exit(1)
		}
		opts.Sizes, err = bedunionbedg.ReadChromSizes(gr)
		gr.Close()
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading genome file: %v\n", err)
			os.Exit(1)
		}
	}

	readers := make([]io.Reader, 0, len(inputFiles))
	closers := make([]io.Closer, 0, len(inputFiles))
	for _, f := range inputFiles {
		rc, err := iohelper.OpenReader(f)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening input %q: %v\n", f, err)
			os.Exit(1)
		}
		readers = append(readers, rc)
		closers = append(closers, rc)
	}
	defer func() {
		for _, c := range closers {
			c.Close()
		}
	}()

	out, err := iohelper.OpenWriter(outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if err := bedunionbedg.Union(readers, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

// startsWithDash reports whether s begins with '-', marking the end of a
// variadic flag's argument list.
func startsWithDash(s string) bool {
	return len(s) > 0 && s[0] == '-'
}
