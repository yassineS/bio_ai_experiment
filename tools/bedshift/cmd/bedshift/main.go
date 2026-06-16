// bedshift shifts each feature of a BED/GFF/VCF file by a requested number of
// base pairs, clamping to chromosome bounds (mirrors `bedtools shift`).
package main

import (
	"flag"
	"fmt"
	"os"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedshift/pkg/bedshift"
)

const version = "1.0.0"

const usage = `bedshift - Shift each feature by a requested number of base pairs

Usage:
  bedshift [options] -i <bed/gff/vcf> -g <genome> [-s <int> | (-p and -m)]

Description:
  Reads a BED/GFF/VCF file and shifts each feature by a number of base pairs.
  Use -s to shift every feature by the same amount, or -p and -m together to
  shift '+' and '-' strand features by different amounts. With -pct the shift is
  a fraction of each feature's length. Starts are clamped to >= 0 and ends to
  the chromosome length.

Options:
  -i, --input FILE   Input BED/GFF/VCF file ('-' for stdin, default: stdin)
  -o, --output FILE  Output file ('-' for stdout, default: stdout)
  -g, --genome FILE  Chromosome sizes file: 'chrom<TAB>size' per line
                     (also accepts a samtools .fai). Required.
  -s NUM             Shift every feature by NUM bp (integer, or float with -pct).
                     Mutually exclusive with -p/-m.
  -p NUM             Shift '+' strand features by NUM bp. Requires -m.
  -m NUM             Shift '-' strand features by NUM bp. Requires -p.
  -pct, --percentage Treat -s/-p/-m as a fraction of the feature's length.
  -header            Print the header from the input file prior to results.
  -h, --help         Show this help message and exit
  -v, --version      Show version information and exit
`

func main() {
	fs := flag.CommandLine
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	var inputFile, outputFile, genomeFile string
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input file")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output file")
	cliflag.StringVar(fs, &genomeFile, "g", "genome", "", "Chrom-sizes file (required)")

	// -s/-p/-m accept integers or floats; capture as strings so we can detect
	// which were supplied and parse with float32 precision like upstream.
	var sStr, pStr, mStr string
	fs.StringVar(&sStr, "s", "", "Shift all features by NUM bp")
	fs.StringVar(&pStr, "p", "", "Shift '+' strand features by NUM bp")
	fs.StringVar(&mStr, "m", "", "Shift '-' strand features by NUM bp")

	var pct bool
	fs.BoolVar(&pct, "pct", false, "")
	fs.BoolVar(&pct, "percentage", false, "Treat -s/-p/-m as a fraction")

	var printHeader bool
	fs.BoolVar(&printHeader, "header", false, "Print header before results")

	var help, showVersion bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help message")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version information")

	flag.Parse()

	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if showVersion {
		fmt.Printf("bedshift version %s\n", version)
		os.Exit(0)
	}

	opts, err := buildOptions(sStr, pStr, mStr, pct)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(2)
	}

	if genomeFile == "" {
		fmt.Fprintln(os.Stderr, "Error: Need both a BED (-i) and a genome (-g) file.")
		os.Exit(2)
	}

	gr, err := iohelper.OpenReader(genomeFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening genome file: %v\n", err)
		os.Exit(1)
	}
	sizes, err := bedshift.ReadChromSizes(gr)
	gr.Close()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error reading genome file: %v\n", err)
		os.Exit(1)
	}

	input := inputFile
	if input == "" && flag.NArg() > 0 {
		input = flag.Arg(0)
	}
	if input == "stdin" {
		input = "-"
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

	if _, err := bedshift.Shift(in, out, sizes, opts, printHeader); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

// buildOptions resolves the -s / -p / -m / -pct flags into bedshift.Options,
// reproducing upstream's validation: exactly one of (-s) or (-p and -m) must be
// supplied, and -p/-m must be used together.
func buildOptions(sStr, pStr, mStr string, pct bool) (bedshift.Options, error) {
	haveAll := sStr != ""
	havePlus := pStr != ""
	haveMinus := mStr != ""

	// Upstream: if ((haveMinus && havePlus) == haveAll) -> error.
	if (haveMinus && havePlus) == haveAll {
		return bedshift.Options{}, fmt.Errorf("Need -m and -p together or -s alone.")
	}
	// Upstream: if ((!havePlus && haveMinus) || (havePlus && !haveMinus)) -> error.
	if (!havePlus && haveMinus) || (havePlus && !haveMinus) {
		return bedshift.Options{}, fmt.Errorf("Need both -m and -p.")
	}

	opts := bedshift.Options{Fractional: pct}
	if haveAll {
		v, err := parseFloat32(sStr, "-s")
		if err != nil {
			return bedshift.Options{}, err
		}
		opts.ShiftPlus = v
		opts.ShiftMinus = v
	} else {
		p, err := parseFloat32(pStr, "-p")
		if err != nil {
			return bedshift.Options{}, err
		}
		m, err := parseFloat32(mStr, "-m")
		if err != nil {
			return bedshift.Options{}, err
		}
		opts.ShiftPlus = p
		opts.ShiftMinus = m
	}
	return opts, nil
}

// parseFloat32 parses a shift amount with the same 32-bit float precision
// upstream uses (its shift fields are C floats parsed via atof).
func parseFloat32(s, name string) (float32, error) {
	v, err := strconv.ParseFloat(s, 32)
	if err != nil {
		return 0, fmt.Errorf("invalid value %q for %s: %v", s, name, err)
	}
	return float32(v), nil
}
