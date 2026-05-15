// Command bedlinks is a pure-Go reimplementation of `bedtools links`. It
// reads BED intervals and emits an HTML page of UCSC Genome Browser links,
// one row per record. See pkg/bedlinks for the output shape and parity
// matrix.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedlinks/pkg/bedlinks"
)

const version = "0.1.0"

const usage = `bedlinks - emit a UCSC-link HTML page from BED intervals.

Usage:
  bedlinks [OPTIONS] -i <bed> > out.html

I/O:
  -i, --input  FILE     BED input ('-' or empty = stdin). Transparent gzip.
  -o, --output FILE     Output file (default stdout). '-' = stdout.

UCSC options (match upstream 'bedtools links'):
  -base URL             UCSC base URL.  Default: http://genome.ucsc.edu
  -org STRING           Organism token. Default: human
  -db STRING            UCSC build/db.  Default: hg18

Standard:
  -h, --help            Show this help.
  -v, --version         Show version.
`

func main() {
	os.Exit(run(os.Args[1:], os.Stdout, os.Stderr))
}

func run(args []string, stdout, stderr io.Writer) int {
	fs := flag.NewFlagSet("bedlinks", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		inFile, outFile       string
		base, org, db         string
		showHelp, showVersion bool
	)
	cliflag.StringVar(fs, &inFile, "i", "input", "", "Input BED")
	cliflag.StringVar(fs, &outFile, "o", "output", "", "Output file")
	fs.StringVar(&base, "base", bedlinks.DefaultBase, "UCSC base URL")
	fs.StringVar(&org, "org", bedlinks.DefaultOrg, "UCSC organism")
	fs.StringVar(&db, "db", bedlinks.DefaultDB, "UCSC db / build")
	fs.BoolVar(&showHelp, "h", false, "Help")
	fs.BoolVar(&showHelp, "help", false, "Help")
	fs.BoolVar(&showVersion, "v", false, "Version")
	fs.BoolVar(&showVersion, "version", false, "Version")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(stderr, err)
		fmt.Fprint(stderr, usage)
		return 2
	}
	if showHelp {
		fmt.Fprint(stdout, usage)
		return 0
	}
	if showVersion {
		fmt.Fprintln(stdout, version)
		return 0
	}

	in, err := iohelper.OpenReader(inFile)
	if err != nil {
		fmt.Fprintf(stderr, "bedlinks: open input: %v\n", err)
		return 1
	}
	defer in.Close()

	out, err := iohelper.OpenWriter(outFile)
	if err != nil {
		fmt.Fprintf(stderr, "bedlinks: open output: %v\n", err)
		return 1
	}
	defer out.Close()

	bedFile := inFile
	if bedFile == "" || bedFile == "-" {
		bedFile = "stdin"
	}

	if _, err := bedlinks.Run(in, out, bedlinks.Options{
		Base:    base,
		Org:     org,
		DB:      db,
		BedFile: bedFile,
	}); err != nil {
		fmt.Fprintf(stderr, "bedlinks: %v\n", err)
		return 1
	}
	return 0
}
