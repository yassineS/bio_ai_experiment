// bedflank emits the flanking regions of each BED interval (mirrors
// `bedtools flank`).
package main

import (
	"flag"
	"fmt"
	"math"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedflank/pkg/bedflank"
)

const version = "1.0.0"

const usage = `bedflank - Emit the flanking regions of each BED interval

Usage:
  bedflank -i <input.bed> -g <genome.sizes> [-b N | -l L -r R] [options]

Description:
  Like bedtools flank: for each interval [s, e) on a chromosome of size C,
  emit the left flank [max(0, s-l), s) and the right flank [e, min(C, e+r)).
  The interval itself is NOT emitted (use bedslop for that). Empty flanks
  are skipped. With -s/--strand the flanks are interpreted relative to the
  transcribed strand (l and r swap on '-'-strand entries).

Options:
  -i, --input FILE         Input BED file ('-' for stdin, default: stdin)
  -o, --output FILE        Output BED file ('-' for stdout, default: stdout)
  -g, --genome FILE        Chromosome sizes file: 'chrom<TAB>size' per line
                           (also accepts samtools .fai). Required.
  -b NUM                   Symmetric flank size (NUM bases per side, or a
                           fraction with --pct). Mutually exclusive with -l/-r.
  -l NUM                   Left ("upstream") flank size. Requires -r.
  -r NUM                   Right ("downstream") flank size. Requires -l.
  -s, --strand             Respect strand: swap l/r on '-'-strand records.
  -pct, --percentage       Treat -b/-l/-r as fractions (0..1) of the
                           interval length rather than absolute base counts.
  -h, --help               Show this help message and exit
  -v, --version            Show version information and exit

Examples:
  # 50bp flank on each side
  bedflank -i input.bed -g hg38.sizes -b 50 > flanks.bed

  # 100bp upstream, 25bp downstream
  bedflank -i input.bed -g hg38.sizes -l 100 -r 25 > flanks.bed

  # Strand-aware promoter-like flanks
  bedflank -i tss.bed -g hg38.sizes -l 1000 -r 100 -s > regions.bed

  # 10% of each interval as flank
  bedflank -i input.bed -g hg38.sizes -b 0.1 --pct > flanks.bed
`

func main() {
	fs := flag.CommandLine
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	var inputFile, outputFile, genomeFile string
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input BED file")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output BED file")
	cliflag.StringVar(fs, &genomeFile, "g", "genome", "", "Chrom-sizes file (required)")

	// Use sentinel NaN so we can tell whether each was set on the command line.
	var bothVal, leftVal, rightVal float64
	bothVal = math.NaN()
	leftVal = math.NaN()
	rightVal = math.NaN()
	fs.Float64Var(&bothVal, "b", math.NaN(), "Symmetric flank size")
	fs.Float64Var(&leftVal, "l", math.NaN(), "Left flank size (requires -r)")
	fs.Float64Var(&rightVal, "r", math.NaN(), "Right flank size (requires -l)")

	var strandSpec bool
	cliflag.BoolVar(fs, &strandSpec, "s", "strand", false, "Respect strand")

	var pct bool
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
		fmt.Printf("bedflank version %s\n", version)
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

	gr, err := iohelper.OpenReader(genomeFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening genome file: %v\n", err)
		os.Exit(1)
	}
	sizes, err := bedflank.ReadChromSizes(gr)
	gr.Close()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error reading genome file: %v\n", err)
		os.Exit(1)
	}

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

	if _, err := bedflank.Flank(in, out, os.Stderr, sizes, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

// buildOptions resolves the set of -b / -l / -r flags into bedflank.Options or
// returns an error explaining how they conflict.
func buildOptions(bothVal, leftVal, rightVal float64, strand, pct bool) (bedflank.Options, error) {
	bothSet := !math.IsNaN(bothVal)
	leftSet := !math.IsNaN(leftVal)
	rightSet := !math.IsNaN(rightVal)
	if bothSet && (leftSet || rightSet) {
		return bedflank.Options{}, fmt.Errorf("-b is mutually exclusive with -l/-r")
	}
	if !bothSet && (leftSet != rightSet) {
		return bedflank.Options{}, fmt.Errorf("-l and -r must be used together")
	}
	if !bothSet && !leftSet && !rightSet {
		return bedflank.Options{}, fmt.Errorf("must specify -b or both -l and -r")
	}
	opts := bedflank.Options{StrandSpec: strand, Pct: pct}
	if bothSet {
		opts.Both = true
		opts.BothAdd = bothVal
	} else {
		opts.LeftAdd = leftVal
		opts.RightAdd = rightVal
	}
	return opts, nil
}
