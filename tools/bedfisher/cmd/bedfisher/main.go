// bedfisher computes Fisher's exact test of overlap enrichment between two
// BED files over a genome (Go port of `bedtools fisher`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedfisher/pkg/bedfisher"
)

const version = "1.0.0"

const usage = `bedfisher - Fisher's exact test of overlap enrichment

Usage:
  bedfisher -a <A.bed> -b <B.bed> -g <genome> [options]

Description:
  Builds the 2x2 contingency table of overlap counts between A and B,
  estimates the number of "possible" intervals from the mean interval
  length and genome size, and runs Fisher's exact test (two-tailed) on
  the resulting table. Output matches 'bedtools fisher'.

Options:
  -a, --a FILE         BED file A (queries; required, '-' for stdin)
  -b, --b FILE         BED file B (database; required)
  -g, --g FILE         Chrom-sizes / genome file (required)
      --output FILE    Output file (default: stdout)
  -f FRACTION          Minimum fraction of A overlapped
  -F FRACTION          Minimum fraction of B overlapped
  -r, --reciprocal     Require -f to apply to BOTH A and B (same threshold)
  -s, --strand         Same-strand overlaps only
  -S, --opposite-strand Opposite-strand overlaps only
  -m, --merge          Merge overlapping A records before the test
  -h, --help           Show this help message
  -v, --version        Show version information
`

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, "Error:", err)
		os.Exit(1)
	}
}

func run(argv []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bedfisher", flag.ContinueOnError)
	fs.SetOutput(stderr)

	var (
		fileA, fileB, fileG, output string
		fractionA, fractionB        float64
		reciprocal, same, opposite  bool
		mergeA                      bool
		help, showVer               bool
	)
	cliflag.StringVar(fs, &fileA, "a", "", "", "BED file A")
	cliflag.StringVar(fs, &fileB, "b", "", "", "BED file B")
	cliflag.StringVar(fs, &fileG, "g", "", "", "Genome file")
	cliflag.StringVar(fs, &output, "", "output", "", "Output file (default: stdout)")
	cliflag.Float64Var(fs, &fractionA, "f", "fraction-a", 0.0, "Fraction of A overlapped")
	cliflag.Float64Var(fs, &fractionB, "F", "fraction-b", 0.0, "Fraction of B overlapped")
	cliflag.BoolVar(fs, &reciprocal, "r", "reciprocal", false, "Reciprocal -f")
	cliflag.BoolVar(fs, &same, "s", "strand", false, "Same-strand overlaps only")
	cliflag.BoolVar(fs, &opposite, "S", "opposite-strand", false, "Opposite-strand overlaps only")
	cliflag.BoolVar(fs, &mergeA, "m", "merge", false, "Pre-merge A")
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
		fmt.Fprintf(stdout, "bedfisher version %s\n", version)
		return nil
	}
	if fileA == "" || fileB == "" || fileG == "" {
		return fmt.Errorf("-a, -b, and -g are all required (use -h for help)")
	}

	rA, err := iohelper.OpenReader(fileA)
	if err != nil {
		return fmt.Errorf("opening A: %w", err)
	}
	defer rA.Close()
	rB, err := iohelper.OpenReader(fileB)
	if err != nil {
		return fmt.Errorf("opening B: %w", err)
	}
	defer rB.Close()
	// The genome file is a small chrom-sizes file; open it as a normal file
	// (not via iohelper, since '-' for stdin would conflict with -a).
	gF, err := os.Open(fileG)
	if err != nil {
		return fmt.Errorf("opening -g: %w", err)
	}
	defer gF.Close()
	w, err := iohelper.OpenWriter(output)
	if err != nil {
		return fmt.Errorf("opening output: %w", err)
	}
	defer w.Close()

	opts := bedfisher.Options{
		FractionA:      fractionA,
		FractionB:      fractionB,
		Reciprocal:     reciprocal,
		SameStrand:     same,
		OppositeStrand: opposite,
		MergeInputs:    mergeA,
	}
	if _, err := bedfisher.Run(rA, rB, gF, w, opts); err != nil {
		return err
	}
	return nil
}
