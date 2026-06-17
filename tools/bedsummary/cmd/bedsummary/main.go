// bedsummary reports per-chromosome interval summary statistics for a
// BED/GFF/VCF file against a genome (Go port of `bedtools summary`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedsummary/pkg/bedsummary"
)

const usage = `bedsummary - Per-chromosome interval summary statistics

Usage:
  bedsummary -i FILE.bed -g GENOME [options]

Options:
  -i, --input FILE        Input BED/GFF/VCF file (required, '-' for stdin)
  -g, --genome FILE       Genome (chrom-sizes) file (required)
  -o, --output FILE       Output file (default: stdout)
  --no-header             Suppress the column-header line
  -h, --help              Show this help
  -v, --version           Show version and exit

Output (TSV):
  chrom  chrom_length  num_ivls  total_ivl_bp  chrom_frac_genome
  frac_all_ivls  frac_all_bp  min  max  mean

Notes:
  - Chromosomes are reported in the order they appear in the genome file.
  - Chromosomes with no intervals are reported with -1 min/max/mean.
  - A final "all" row aggregates over the whole input.
`

const version = "1.0.0"

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func run(args []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bedsummary", flag.ContinueOnError)
	fs.SetOutput(stderr)
	fs.Usage = func() { fmt.Fprint(stderr, usage) }

	var (
		input    string
		genome   string
		output   string
		noHeader bool
		showHelp bool
		showVer  bool
	)
	cliflag.StringVar(fs, &input, "i", "input", "", "Input BED file")
	cliflag.StringVar(fs, &genome, "g", "genome", "", "Genome (chrom-sizes) file")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file")
	fs.BoolVar(&noHeader, "no-header", false, "Suppress header row")
	cliflag.BoolVar(fs, &showHelp, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVer, "v", "version", false, "Show version")

	if err := fs.Parse(args); err != nil {
		return err
	}
	if showHelp {
		fmt.Fprint(stdout, usage)
		return nil
	}
	if showVer {
		fmt.Fprintln(stdout, version)
		return nil
	}
	if input == "" {
		fmt.Fprint(stderr, usage)
		return fmt.Errorf("error: -i/--input is required")
	}
	if genome == "" {
		fmt.Fprint(stderr, usage)
		return fmt.Errorf("error: -g/--genome is required")
	}

	gr, err := iohelper.OpenReader(genome)
	if err != nil {
		return fmt.Errorf("opening genome: %w", err)
	}
	defer gr.Close()
	g, err := bedsummary.ParseGenome(gr)
	if err != nil {
		return err
	}

	r, err := iohelper.OpenReader(input)
	if err != nil {
		return fmt.Errorf("opening input: %w", err)
	}
	defer r.Close()

	w, err := iohelper.OpenWriter(output)
	if err != nil {
		return fmt.Errorf("opening output: %w", err)
	}
	defer w.Close()

	return bedsummary.Run(r, g, w, bedsummary.Options{NoHeader: noHeader})
}
