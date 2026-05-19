// bedsummary computes per-chromosome interval-length summary statistics for
// a BED file (Go port of `bedtools summary`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedsummary/pkg/bedsummary"
)

const usage = `bedsummary - Per-chromosome interval-length summary stats

Usage:
  bedsummary -i FILE.bed [options]

Options:
  -i, --input FILE        Input BED file (required, '-' for stdin)
  -o, --output FILE       Output file (default: stdout)
  --no-header             Suppress the column-header line
  --skip-all              Suppress the trailing "all" aggregate row
  -h, --help              Show this help
  -v, --version           Show version and exit

Output (TSV):
  chrom  num_ivls  total_ivl_bp  min_ivl_bp  max_ivl_bp  mean_ivl_bp  median_ivl_bp

Notes:
  - Chromosomes are emitted in the order they first appear in input.
  - Mean / median are emitted as integers when integer-valued, otherwise
    formatted with 3-decimal precision.
`

const version = "0.1.0"

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
		output   string
		noHeader bool
		skipAll  bool
		showHelp bool
		showVer  bool
	)
	cliflag.StringVar(fs, &input, "i", "input", "", "Input BED file")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file")
	fs.BoolVar(&noHeader, "no-header", false, "Suppress header row")
	fs.BoolVar(&skipAll, "skip-all", false, "Suppress aggregate 'all' row")
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

	return bedsummary.Run(r, w, bedsummary.Options{NoHeader: noHeader, SkipAll: skipAll})
}
