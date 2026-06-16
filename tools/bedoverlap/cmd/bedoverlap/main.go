// bedoverlap computes the overlap (positive) or distance (negative) between two
// intervals described by four columns of each input line, appending the result
// as a new trailing column (mirrors `bedtools overlap`, aka getOverlap).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedoverlap/pkg/bedoverlap"
)

const version = "1.0.0"

const usage = `bedoverlap - Compute overlap/distance between two intervals on a line

Usage:
  bedoverlap -i <file> -cols s1,e1,s2,e2

Description:
  Reads each input line, computes the overlap (positive values) or distance
  (negative values) between the two intervals given by the four 1-based,
  comma-separated columns, and reports the result as a new final column on the
  same line. The overlap is min(end1,end2) - max(start1,start2).

Options:
  -i, --input FILE   Input file ('-' or 'stdin' for stdin, default: stdin)
  -o, --output FILE  Output file ('-' for stdout, default: stdout)
      --cols COLS    Comma-separated 1-based columns: start1,end1,start2,end2.
                     Required.
  -h, --help         Show this help message and exit
  -v, --version      Show version information and exit

Example:
  $ bedtools window -a A.bed -b B.bed -w 10 | bedoverlap -i stdin -cols 2,3,6,7
  chr1  10  20  A   chr1  15  25  B   5
  chr1  10  20  C   chr1  25  35  D   -5
`

func main() {
	fs := flag.CommandLine
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	var inputFile, outputFile, colsSpec string
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input file")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output file")
	// "cols" is conventionally a single-dash flag in bedtools. Go's flag
	// package treats `-cols` and `--cols` identically.
	fs.StringVar(&colsSpec, "cols", "", "Comma-separated columns: start1,end1,start2,end2")

	var help, showVersion bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help message")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version information")

	flag.Parse()

	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if showVersion {
		fmt.Printf("bedoverlap version %s\n", version)
		os.Exit(0)
	}

	cols, err := bedoverlap.ParseCols(colsSpec)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(2)
	}

	input := inputFile
	if input == "" && flag.NArg() > 0 {
		input = flag.Arg(0)
	}
	if input == "stdin" {
		input = "-"
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

	if err := bedoverlap.Overlap(in, out, cols); err != nil {
		fmt.Fprintf(os.Stderr, "%v\n", err)
		os.Exit(1)
	}
}
