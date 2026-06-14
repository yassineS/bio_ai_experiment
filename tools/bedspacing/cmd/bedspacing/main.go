// Command bedspacing is a pure-Go reimplementation of `bedtools spacing`. It
// appends a "distance to previous interval on the same chromosome" column.
// See pkg/bedspacing for behaviour.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedspacing/pkg/bedspacing"
)

const version = "0.1.0"

const usage = `bedspacing - distance between adjacent intervals on the same chromosome.

Usage:
  bedspacing [-i BED] [-o OUT]

I/O:
  -i, --input FILE        Input BED or SAM/BAM ('-' or empty = stdin).
                          SAM/BAM is auto-detected and emitted as BED12.
                          Transparent gzip.
  -o, --output FILE       Output file (default stdout). '-' = stdout.

Standard:
  -h, --help              Show this help.
  -v, --version           Show version.

Each output row is the input row plus a new tab-separated column:
  "."   for the first interval on its chromosome
  "-1"  if the interval overlaps the previous one on its chromosome
  "0"   if it exactly abuts (this.start == prev.end)
  N>0   otherwise: this.start - prev.end (gap in bases)

Note: bedspacing does NOT sort. Pipe a sorted BED in (e.g. via 'bedsort') to
get the conventional genome-sorted spacing report.
`

func main() {
	fs := flag.NewFlagSet("bedspacing", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		inFile      string
		outFile     string
		showHelp    bool
		showVersion bool
	)
	cliflag.StringVar(fs, &inFile, "i", "input", "", "Input BED file")
	cliflag.StringVar(fs, &outFile, "o", "output", "", "Output file")
	fs.BoolVar(&showHelp, "h", false, "Help")
	fs.BoolVar(&showHelp, "help", false, "Help")
	fs.BoolVar(&showVersion, "v", false, "Version")
	fs.BoolVar(&showVersion, "version", false, "Version")

	if err := fs.Parse(os.Args[1:]); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, usage)
		os.Exit(2)
	}
	if showHelp {
		fmt.Print(usage)
		return
	}
	if showVersion {
		fmt.Println(version)
		return
	}

	in, err := iohelper.OpenReader(inFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bedspacing: open input: %v\n", err)
		os.Exit(1)
	}
	defer in.Close()

	out, err := iohelper.OpenWriter(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bedspacing: open output: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if _, err := bedspacing.Spacing(in, out); err != nil {
		fmt.Fprintf(os.Stderr, "bedspacing: %v\n", err)
		os.Exit(1)
	}
}
