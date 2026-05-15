// Command bedpairtobed is a pure-Go reimplementation of `bedtools pairtobed`.
// It reports overlaps between a BEDPE A file and a regular BED B file.
// See pkg/bedpairtobed for the algorithm and parity matrix.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedpairtobed/pkg/bedpairtobed"
)

const version = "0.1.0"

const usage = `bedpairtobed - report overlaps between BEDPE pairs and a BED file.

Usage:
  bedpairtobed -a <BEDPE> -b <BED> [OPTIONS]

I/O:
  -a FILE                 BEDPE A input ('-' = stdin). Transparent gzip.
  -b FILE                 BED B input.
  -o, --output FILE       Output file (default stdout). '-' = stdout.

Reporting:
  -type {either|both|notboth|neither|xor|notxor}
                          Default: either.
                          - either : >=1 end of A overlaps any B record.
                          - both   : BOTH ends of A overlap some B record.
                          - notboth: NOT both ends overlap a B record.
                          - neither: NO end overlaps any B record.
                          - xor    : exactly one end overlaps any B record.
                          - notxor : either both or neither end overlaps
                                     (deviation: not in upstream).

Filters:
  -f FRAC                 Min fraction of A end length covered by B
                          (default 1e-9 -> 1bp).
  -s, --same-strand       Require matching strands between A end and B hit.
  -S, --opposite-strand   Require opposite strands.
  -is, --ignore-strand    Ignore strands entirely (overrides -s/-S).

Standard:
  -h, --help              Show this help.
  -v, --version           Show version.
`

func main() { os.Exit(run(os.Args[1:], os.Stdout, os.Stderr)) }

func run(args []string, stdout, stderr io.Writer) int {
	fs := flag.NewFlagSet("bedpairtobed", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		aFile, bFile, outFile string
		typeStr               string
		fracA                 float64
		same, opp, ignore     bool
		showHelp, showVersion bool
	)
	fs.StringVar(&aFile, "a", "", "")
	fs.StringVar(&bFile, "b", "", "")
	cliflag.StringVar(fs, &outFile, "o", "output", "", "Output file")
	fs.StringVar(&typeStr, "type", "either", "")
	fs.Float64Var(&fracA, "f", 1e-9, "")
	cliflag.BoolVar(fs, &same, "s", "same-strand", false, "Require same strand")
	cliflag.BoolVar(fs, &opp, "S", "opposite-strand", false, "Require opposite strand")
	cliflag.BoolVar(fs, &ignore, "is", "ignore-strand", false, "Ignore strand")
	cliflag.BoolVar(fs, &showHelp, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version")

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
	if aFile == "" || bFile == "" {
		fmt.Fprintln(stderr, "bedpairtobed: -a and -b are required")
		return 2
	}

	a, err := iohelper.OpenReader(aFile)
	if err != nil {
		fmt.Fprintf(stderr, "bedpairtobed: open -a: %v\n", err)
		return 1
	}
	defer a.Close()
	b, err := iohelper.OpenReader(bFile)
	if err != nil {
		fmt.Fprintf(stderr, "bedpairtobed: open -b: %v\n", err)
		return 1
	}
	defer b.Close()
	out, err := iohelper.OpenWriter(outFile)
	if err != nil {
		fmt.Fprintf(stderr, "bedpairtobed: open output: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bedpairtobed.Run(a, b, out, bedpairtobed.Options{
		Type:           bedpairtobed.Type(typeStr),
		MinFractionA:   fracA,
		SameStrand:     same,
		OppositeStrand: opp,
		IgnoreStrand:   ignore,
	}); err != nil {
		fmt.Fprintf(stderr, "bedpairtobed: %v\n", err)
		return 1
	}
	return 0
}
