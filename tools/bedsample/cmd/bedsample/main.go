// Command bedsample is a pure-Go reimplementation of `bedtools sample`. It
// draws N random records (without replacement, reservoir-style) from a BED
// input. See pkg/bedsample for behaviour.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedsample/pkg/bedsample"
)

const version = "0.1.0"

const usage = `bedsample - random subsample of N BED records (without replacement).

Usage:
  bedsample -n N [-seed SEED] [-header] [-i FILE] [-o OUT]

Required:
  -n, --number N          Number of records to draw (must be > 0).

Optional:
  -seed, --seed SEED      PRNG seed for deterministic output. 0 = time-based
                          (default).
  -header, --header       Forward '#', 'track', and 'browser' lines verbatim
                          to the output before the sampled records.

I/O:
  -i, --input FILE        Input BED ('-' or empty = stdin). Transparent gzip.
  -o, --output FILE       Output file (default stdout). '-' = stdout.

Standard:
  -h, --help              Show this help.
  -v, --version           Show version.
`

func main() {
	fs := flag.NewFlagSet("bedsample", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		inFile      string
		outFile     string
		n           int
		seed        int64
		header      bool
		showHelp    bool
		showVersion bool
	)
	cliflag.StringVar(fs, &inFile, "i", "input", "", "Input BED file")
	cliflag.StringVar(fs, &outFile, "o", "output", "", "Output file")
	cliflag.IntVar(fs, &n, "n", "number", 0, "Number of records to draw")
	// math/rand seed is int64, but cliflag only exposes int. We register
	// the long form ourselves to match.
	fs.Int64Var(&seed, "seed", 0, "PRNG seed (0 = time-based)")
	fs.BoolVar(&header, "header", false, "Forward header lines")
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
	if n <= 0 {
		fmt.Fprintln(os.Stderr, "bedsample: -n N is required and must be > 0")
		fmt.Fprint(os.Stderr, usage)
		os.Exit(2)
	}

	in, err := iohelper.OpenReader(inFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bedsample: open input: %v\n", err)
		os.Exit(1)
	}
	defer in.Close()

	out, err := iohelper.OpenWriter(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bedsample: open output: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if _, err := bedsample.Sample(in, out, bedsample.Options{N: n, Seed: seed, Header: header}); err != nil {
		fmt.Fprintf(os.Stderr, "bedsample: %v\n", err)
		os.Exit(1)
	}
}
