// bedslop extends BED intervals by a fixed or fractional amount, clipping to
// chromosome boundaries (mirrors `bedtools slop`).
package main

import (
	"flag"
	"fmt"
	"math"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedslop/pkg/bedslop"
)

const version = "1.0.0"

const usage = `bedslop - Extend BED intervals by a fixed or fractional amount

Usage:
  bedslop [options] -g <genome.sizes>

Description:
  Reads a BED file and extends each interval by a fixed number of bases (or a
  fraction of the interval length, with --pct) on the left, right, or both
  sides. Extensions are clipped to [0, chromSize] using the supplied chrom-
  sizes file. Intervals whose extended length is non-positive are dropped with
  a warning written to stderr.

Options:
  -i, --input FILE         Input BED file ('-' for stdin, default: stdin)
  -o, --output FILE        Output BED file ('-' for stdout, default: stdout)
  -g, --genome FILE        Chromosome sizes file: 'chrom<TAB>size' per line
                           (also accepts samtools .fai). Required.
  -b NUM                   Extend by NUM bases (or fraction with --pct) on
                           both sides. Mutually exclusive with -l/-r.
  -l NUM                   Extend by NUM on the left ("upstream") side.
                           Requires -r.
  -r NUM                   Extend by NUM on the right ("downstream") side.
                           Requires -l.
  -s, --strand             Respect strand: swap left/right semantics for
                           records on the '-' strand.
  -pct, --percentage       Treat -b/-l/-r as fractions (0..1) of the interval
                           length rather than absolute base counts.
  -h, --help               Show this help message and exit
  -v, --version            Show version information and exit

Examples:
  # Extend every interval by 50bp on each side
  bedslop -i input.bed -g hg38.sizes -b 50 > extended.bed

  # Asymmetric: 100bp upstream, 25bp downstream
  bedslop -i input.bed -g hg38.sizes -l 100 -r 25 > extended.bed

  # Strand-aware: upstream means 5' relative to transcription
  bedslop -i input.bed -g hg38.sizes -l 100 -r 25 -s > extended.bed

  # Extend by 25% of the interval length on each side
  bedslop -i input.bed -g hg38.sizes -b 0.25 --pct > extended.bed
`

func main() {
	fs := flag.CommandLine
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	var inputFile, outputFile, genomeFile string
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input BED file")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output BED file")
	cliflag.StringVar(fs, &genomeFile, "g", "genome", "", "Chrom-sizes file (required)")

	// Use sentinel NaN to detect whether each flag was set.
	bothSet := math.NaN()
	leftSet := math.NaN()
	rightSet := math.NaN()
	var bothVal, leftVal, rightVal float64
	bothVal = bothSet
	leftVal = leftSet
	rightVal = rightSet
	fs.Float64Var(&bothVal, "b", math.NaN(), "Extend NUM bases on both sides")
	fs.Float64Var(&leftVal, "l", math.NaN(), "Extend NUM bases on the left")
	fs.Float64Var(&rightVal, "r", math.NaN(), "Extend NUM bases on the right")

	var strandSpec bool
	cliflag.BoolVar(fs, &strandSpec, "s", "strand", false, "Respect strand")

	var pct bool
	// "pct" is conventionally a single-dash short flag in bedtools. Go's flag
	// package treats `-pct` and `--pct` identically, so registering both
	// `pct` and `percentage` on the FlagSet covers `-pct`, `--pct`, and
	// `--percentage`.
	fs.BoolVar(&pct, "pct", false, "")
	fs.BoolVar(&pct, "percentage", false, "Treat -b/-l/-r as fractions")

	var help, showVersion bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help message")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version information")

	flag.Parse()

	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if showVersion {
		fmt.Printf("bedslop version %s\n", version)
		os.Exit(0)
	}

	if genomeFile == "" {
		fmt.Fprintln(os.Stderr, "Error: -g/--genome is required")
		os.Exit(2)
	}

	opts, err := buildOptions(bothVal, leftVal, rightVal, strandSpec, pct)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(2)
	}

	// Load chromosome sizes.
	gr, err := iohelper.OpenReader(genomeFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening genome file: %v\n", err)
		os.Exit(1)
	}
	sizes, err := bedslop.ReadChromSizes(gr)
	gr.Close()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error reading genome file: %v\n", err)
		os.Exit(1)
	}

	// Determine input path: prefer -i, fall back to positional.
	input := inputFile
	if input == "" && flag.NArg() > 0 {
		input = flag.Arg(0)
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

	if _, err := bedslop.Slop(in, out, os.Stderr, sizes, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

// buildOptions resolves the set of -b / -l / -r flags into bedslop.Options or
// returns an error explaining how they conflict.
func buildOptions(bothVal, leftVal, rightVal float64, strand, pct bool) (bedslop.Options, error) {
	bothSet := !math.IsNaN(bothVal)
	leftSet := !math.IsNaN(leftVal)
	rightSet := !math.IsNaN(rightVal)
	if bothSet && (leftSet || rightSet) {
		return bedslop.Options{}, fmt.Errorf("-b is mutually exclusive with -l/-r")
	}
	if !bothSet && (leftSet != rightSet) {
		return bedslop.Options{}, fmt.Errorf("-l and -r must be used together")
	}
	if !bothSet && !leftSet && !rightSet {
		return bedslop.Options{}, fmt.Errorf("must specify -b or both -l and -r")
	}
	opts := bedslop.Options{StrandSpec: strand, Pct: pct}
	if bothSet {
		opts.Both = true
		opts.BothAdd = bothVal
	} else {
		opts.LeftAdd = leftVal
		opts.RightAdd = rightVal
	}
	return opts, nil
}
