// bedgenomecov computes per-base or summarised coverage for BED intervals,
// mirroring `bedtools genomecov`.
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedgenomecov/pkg/bedgenomecov"
)

const version = "1.0.0"

const usage = `bedgenomecov - Coverage histogram, bedGraph or per-base depth for BED intervals

Usage:
  bedgenomecov -i <intervals.bed> -g <chrom.sizes> [options]

Description:
  Reads BED/bedGraph intervals and a chromosome-sizes genome file and produces
  coverage information. By default a histogram is written; switch to
  --bedGraph/-bg, -bga, -d or -dz for alternative output styles.

Options:
  -i, --input FILE        Input BED file (default: stdin, '-' for stdin)
  -g, --genome FILE       Chromosome sizes file (chrom<TAB>size per line, required)
      --output FILE       Output file (default: stdout)
  -bg, --bedGraph         Emit non-zero runs of constant depth as bedGraph
  -bga                    Emit every run of constant depth (includes zero)
  -d,  --per-base         Per-base depth (1-based positions)
  -dz, --per-base-nonzero Per-base depth, skip positions with depth 0
      --strand +|-        Count only intervals on the given strand (BED6 col 6)
      --max N             Cap histogram depth at N (--max-depth alias)
      --scale FLOAT       Multiply every depth by FLOAT (default 1.0)
  -5,  --five-prime       Count only the 5'-most base of each interval
  -3,  --three-prime      Count only the 3'-most base of each interval
      --trackline         Prepend a UCSC trackline to -bg/-bga output
      --trackopts STR     Extra trackline options appended after "track"
  -h,  --help             Show this help message
  -v,  --version          Show version information

Examples:
  bedgenomecov -i reads.bed -g hg38.sizes                  # histogram (default)
  bedgenomecov -i reads.bed -g hg38.sizes -bg              # bedGraph (non-zero)
  bedgenomecov -i reads.bed -g hg38.sizes -bga             # bedGraph (with zero)
  bedgenomecov -i reads.bed -g hg38.sizes -d               # per-base
  bedgenomecov -i reads.bed -g hg38.sizes -dz              # per-base, non-zero only
  bedgenomecov -i reads.bed.gz -g hg38.sizes --scale=0.5   # normalised

Notes:
  - Coordinates are 0-based, half-open; per-base output uses 1-based positions
    for compatibility with bedtools.
  - The genome file establishes chromosome ordering for all output modes.
  - Use '-' for stdin/stdout and '--' to end option parsing.
`

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, "Error:", err)
		os.Exit(1)
	}
}

// run is the testable entry point — kept tiny so most logic lives in pkg/.
func run(argv []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bedgenomecov", flag.ContinueOnError)
	fs.SetOutput(stderr)

	var (
		inputFile  string
		genomeFile string
		outputFile string

		bedGraph    bool
		bedGraphAll bool
		perBase     bool
		perBaseNZ   bool

		strand    string
		maxDepth  int
		scale     float64
		fivePrime bool
		threePrm  bool
		trackLine bool
		trackOpts string

		help    bool
		showVer bool
	)

	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input BED file (default: stdin)")
	cliflag.StringVar(fs, &genomeFile, "g", "genome", "", "Chromosome sizes file (required)")
	cliflag.StringVar(fs, &outputFile, "", "output", "", "Output file (default: stdout)")

	cliflag.BoolVar(fs, &bedGraph, "bg", "bedGraph", false, "Emit non-zero bedGraph runs")
	cliflag.BoolVar(fs, &bedGraphAll, "bga", "bedgraph-all", false, "Emit all bedGraph runs (incl. zero)")
	cliflag.BoolVar(fs, &perBase, "d", "per-base", false, "Per-base depth (1-based)")
	cliflag.BoolVar(fs, &perBaseNZ, "dz", "per-base-nonzero", false, "Per-base depth, skip zero")

	cliflag.StringVar(fs, &strand, "", "strand", "", "Strand filter (+ or -)")
	cliflag.IntVar(fs, &maxDepth, "", "max", 0, "Cap histogram depth")
	fs.IntVar(&maxDepth, "max-depth", 0, "")
	cliflag.Float64Var(fs, &scale, "", "scale", 1.0, "Depth multiplier")

	cliflag.BoolVar(fs, &fivePrime, "5", "five-prime", false, "Count only 5' end")
	cliflag.BoolVar(fs, &threePrm, "3", "three-prime", false, "Count only 3' end")
	cliflag.BoolVar(fs, &trackLine, "", "trackline", false, "Emit UCSC trackline header")
	cliflag.StringVar(fs, &trackOpts, "", "trackopts", "", "Extra trackline options")

	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVer, "v", "version", false, "Show version")

	fs.Usage = func() { fmt.Fprint(stderr, usage) }

	if err := fs.Parse(argv); err != nil {
		return err
	}
	if help {
		fmt.Fprint(stderr, usage)
		return nil
	}
	if showVer {
		fmt.Fprintf(stdout, "bedgenomecov version %s\n", version)
		return nil
	}

	if genomeFile == "" {
		return fmt.Errorf("missing -g/--genome (use -h for help)")
	}

	// Pick output mode (mutually exclusive).
	mode := bedgenomecov.ModeHistogram
	modeCount := 0
	if bedGraph {
		mode = bedgenomecov.ModeBedGraph
		modeCount++
	}
	if bedGraphAll {
		mode = bedgenomecov.ModeBedGraphAll
		modeCount++
	}
	if perBase {
		mode = bedgenomecov.ModePerBase
		modeCount++
	}
	if perBaseNZ {
		mode = bedgenomecov.ModePerBaseNonZero
		modeCount++
	}
	if modeCount > 1 {
		return fmt.Errorf("the output-mode flags -bg/-bga/-d/-dz are mutually exclusive")
	}

	// Resolve input file. Positional arg may be used if -i not given.
	if inputFile == "" && fs.NArg() > 0 {
		inputFile = fs.Arg(0)
	}

	genomeR, err := iohelper.OpenReader(genomeFile)
	if err != nil {
		return fmt.Errorf("opening genome file: %w", err)
	}
	defer genomeR.Close()
	genome, err := bedgenomecov.ReadGenome(genomeR)
	if err != nil {
		return err
	}

	inR, err := iohelper.OpenReader(inputFile)
	if err != nil {
		return fmt.Errorf("opening input: %w", err)
	}
	defer inR.Close()

	outW, err := iohelper.OpenWriter(outputFile)
	if err != nil {
		return fmt.Errorf("opening output: %w", err)
	}
	defer outW.Close()

	opts := bedgenomecov.Options{
		Mode:       mode,
		Strand:     strand,
		MaxDepth:   maxDepth,
		Scale:      scale,
		FivePrime:  fivePrime,
		ThreePrime: threePrm,
		TrackLine:  trackLine,
		TrackOpts:  trackOpts,
	}
	return bedgenomecov.Run(inR, genome, outW, opts)
}
