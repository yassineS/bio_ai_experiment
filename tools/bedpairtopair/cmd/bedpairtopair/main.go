// Command bedpairtopair is a pure-Go reimplementation of `bedtools pairtopair`.
// It reports overlaps between two BEDPE files. See pkg/bedpairtopair for the
// algorithm and parity matrix.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedpairtopair/pkg/bedpairtopair"
)

const version = "0.1.0"

const usage = `bedpairtopair - report overlaps between two BEDPE files.

Usage:
  bedpairtopair -a <BEDPE> -b <BEDPE> [OPTIONS]

I/O:
  -a FILE                 BEDPE A input ('-' = stdin). Transparent gzip.
  -b FILE                 BEDPE B input.
  -o, --output FILE       Output file (default stdout). '-' = stdout.

Reporting:
  -type {both|notboth|either|neither}
                          Default: both.
                          - both    : both ends of A overlap the SAME B pair.
                          - notboth : NOT both ends of A overlap the same B.
                          - either  : at least one end of A overlaps any B.
                          - neither : NO end of A overlaps any B.

Filters:
  -f FRAC                 Min fraction of A end length covered by B end
                          (default 1e-9 -> 1bp).
  -slop N                 Add N bp of slop to each end of A before search.
  -ss                     Make -slop strand-aware (extend in strand dir only).
  -is, --ignore-strand    Ignore strand entirely (overrides -s/-S).
  -s, --same-strand       Require matching strands.
  -S, --opposite-strand   Require opposite strands.
  -rdn                    Require A and B pairs to have different names
                          (avoid self-hits between datasets that share names).

Standard:
  -h, --help              Show this help.
  -v, --version           Show version.
`

func main() { os.Exit(run(os.Args[1:], os.Stdout, os.Stderr)) }

func run(args []string, stdout, stderr io.Writer) int {
	fs := flag.NewFlagSet("bedpairtopair", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		aFile, bFile, outFile string
		typeStr               string
		frac                  float64
		slop                  int
		stranded              bool
		same, opp, ignore     bool
		rdn                   bool
		showHelp, showVersion bool
	)
	fs.StringVar(&aFile, "a", "", "")
	fs.StringVar(&bFile, "b", "", "")
	cliflag.StringVar(fs, &outFile, "o", "output", "", "Output file")
	fs.StringVar(&typeStr, "type", "both", "")
	fs.Float64Var(&frac, "f", 1e-9, "")
	fs.IntVar(&slop, "slop", 0, "")
	fs.BoolVar(&stranded, "ss", false, "")
	cliflag.BoolVar(fs, &same, "s", "same-strand", false, "Require same strand")
	cliflag.BoolVar(fs, &opp, "S", "opposite-strand", false, "Require opposite strand")
	cliflag.BoolVar(fs, &ignore, "is", "ignore-strand", false, "Ignore strand")
	fs.BoolVar(&rdn, "rdn", false, "")
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
		fmt.Fprintln(stderr, "bedpairtopair: -a and -b are required")
		return 2
	}

	a, err := iohelper.OpenReader(aFile)
	if err != nil {
		fmt.Fprintf(stderr, "bedpairtopair: open -a: %v\n", err)
		return 1
	}
	defer a.Close()
	b, err := iohelper.OpenReader(bFile)
	if err != nil {
		fmt.Fprintf(stderr, "bedpairtopair: open -b: %v\n", err)
		return 1
	}
	defer b.Close()
	out, err := iohelper.OpenWriter(outFile)
	if err != nil {
		fmt.Fprintf(stderr, "bedpairtopair: open output: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bedpairtopair.Run(a, b, out, bedpairtopair.Options{
		Type:                  bedpairtopair.Type(typeStr),
		MinFraction:           frac,
		Slop:                  slop,
		StrandedSlop:          stranded,
		IgnoreStrand:          ignore,
		SameStrand:            same,
		OppositeStrand:        opp,
		RequireDifferentNames: rdn,
	}); err != nil {
		fmt.Fprintf(stderr, "bedpairtopair: %v\n", err)
		return 1
	}
	return 0
}
