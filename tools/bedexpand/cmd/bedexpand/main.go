// Command bedexpand is a pure-Go reimplementation of `bedtools expand`. It
// reads a tab-delimited file and emits one output row per element of one or
// more comma-separated list columns. See pkg/bedexpand for behaviour.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedexpand/pkg/bedexpand"
)

const version = "0.1.0"

const usage = `bedexpand - expand comma-separated list columns into separate rows.

Usage:
  bedexpand -c COLS [-i FILE] [-o OUT]

Required:
  -c, --columns COLS      1-based comma-separated columns to expand (e.g. '4'
                          or '4,5'). Multiple columns are zipped in lock-step;
                          the user-specified order decides which list supplies
                          the value at each expanded slot (matches upstream).

I/O:
  -i, --input FILE        Input file ('-' or empty = stdin). Transparent gzip.
  -o, --output FILE       Output file (default stdout). '-' = stdout.

Standard:
  -h, --help              Show this help.
  -v, --version           Show version.
`

func main() {
	fs := flag.NewFlagSet("bedexpand", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		inFile      string
		outFile     string
		colsArg     string
		showHelp    bool
		showVersion bool
	)
	cliflag.StringVar(fs, &inFile, "i", "input", "", "Input file")
	cliflag.StringVar(fs, &outFile, "o", "output", "", "Output file")
	cliflag.StringVar(fs, &colsArg, "c", "columns", "", "Columns to expand (1-based, comma-separated)")
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
	if colsArg == "" {
		fmt.Fprintln(os.Stderr, "bedexpand: -c COLS is required")
		fmt.Fprint(os.Stderr, usage)
		os.Exit(2)
	}
	cols, err := bedexpand.ParseColumns(colsArg)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bedexpand: %v\n", err)
		os.Exit(2)
	}

	in, err := iohelper.OpenReader(inFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bedexpand: open input: %v\n", err)
		os.Exit(1)
	}
	defer in.Close()

	out, err := iohelper.OpenWriter(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bedexpand: open output: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if _, err := bedexpand.Expand(in, out, bedexpand.Options{Columns: cols}); err != nil {
		fmt.Fprintf(os.Stderr, "bedexpand: %v\n", err)
		os.Exit(1)
	}
}
