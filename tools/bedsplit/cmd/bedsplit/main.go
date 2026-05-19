// bedsplit splits a BED file into N approximately-equal-sized output files
// (Go port of `bedtools split`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedsplit/pkg/bedsplit"
)

const usage = `bedsplit - Split a BED file into N approximately equal-sized files

Usage:
  bedsplit -i FILE.bed -p PREFIX -n N [-a {simple|size}]

Options:
  -i, --input FILE        Input BED file (required, '-' for stdin)
  -p, --prefix STRING     Output filename prefix (required)
  -n, --number INT        Number of output files (required, >=1)
  -a, --algorithm STRING  Partitioning algorithm: simple | size (default: size)
                            simple = equal record counts per file
                            size   = balance total bp per file (LPT)
  -h, --help              Show help
  -v, --version           Show version

Output:
  Each shard is written to <prefix>.NNNNN.bed (1-based, zero-padded to 5
  digits). A manifest TSV (filename<TAB>total_bp<TAB>num_records) is
  written to stdout.

Notes:
  - If N exceeds the record count, only one file per record is created.
  - Header lines (#, track, browser) and blank lines are dropped from input.
`

const version = "0.1.0"

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func run(args []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bedsplit", flag.ContinueOnError)
	fs.SetOutput(stderr)
	fs.Usage = func() { fmt.Fprint(stderr, usage) }

	var (
		input    string
		prefix   string
		n        int
		alg      string
		showHelp bool
		showVer  bool
	)
	cliflag.StringVar(fs, &input, "i", "input", "", "Input BED file")
	cliflag.StringVar(fs, &prefix, "p", "prefix", "", "Output prefix")
	cliflag.IntVar(fs, &n, "n", "number", 0, "Number of files")
	cliflag.StringVar(fs, &alg, "a", "algorithm", "size", "Algorithm")
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
	if input == "" || prefix == "" || n < 1 {
		fmt.Fprint(stderr, usage)
		return fmt.Errorf("error: -i, -p, and -n (>=1) are required")
	}

	algorithm, err := bedsplit.ParseAlgorithm(alg)
	if err != nil {
		return err
	}

	r, err := iohelper.OpenReader(input)
	if err != nil {
		return fmt.Errorf("opening input: %w", err)
	}
	defer r.Close()

	_, err = bedsplit.Split(r, stdout, bedsplit.Options{N: n, Prefix: prefix, Algorithm: algorithm})
	return err
}
