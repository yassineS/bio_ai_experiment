// bedcomplement emits the genomic regions NOT covered by a sorted BED file
// (mirrors `bedtools complement`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedcomplement/pkg/bedcomplement"
)

const version = "1.0.0"

const usage = `bedcomplement - Emit genomic regions NOT covered by a BED file

Usage:
  bedcomplement [options] -g <genome.sizes>

Description:
  Reads a sorted BED file and emits the complementary intervals (the gaps
  between consecutive intervals, plus the leading gap [0, first.start) and
  the trailing gap [last.end, chromSize)) as BED3.

  The input must be sorted on (chrom, start). bedcomplement detects
  out-of-order input on the fly and aborts with a clear error. Chromosomes
  that appear only in the chrom-sizes file (with no intervals) produce a
  single full-length record 'chrom<TAB>0<TAB>chromSize'.

Options:
  -i, --input FILE         Input BED file ('-' for stdin, default: stdin)
  -o, --output FILE        Output BED file ('-' for stdout, default: stdout)
  -g, --genome FILE        Chromosome sizes file (required). One
                           'chrom<TAB>size' per line (or samtools .fai).
  -h, --help               Show this help message and exit
  -v, --version            Show version information and exit

Examples:
  # Print all regions not covered by genes.bed
  bedcomplement -i genes.sorted.bed -g hg38.sizes > intergenic.bed

  # Read sorted BED from stdin
  cat sorted.bed | bedcomplement -g hg38.sizes > complement.bed
`

func main() {
	fs := flag.CommandLine
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	var inputFile, outputFile, genomeFile string
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input BED file")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output BED file")
	cliflag.StringVar(fs, &genomeFile, "g", "genome", "", "Chrom-sizes file (required)")

	var help, showVersion bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help message")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version information")

	flag.Parse()

	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if showVersion {
		fmt.Printf("bedcomplement version %s\n", version)
		os.Exit(0)
	}
	if genomeFile == "" {
		fmt.Fprintln(os.Stderr, "Error: -g/--genome is required")
		os.Exit(2)
	}

	gr, err := iohelper.OpenReader(genomeFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening genome file: %v\n", err)
		os.Exit(1)
	}
	sizes, order, err := bedcomplement.ReadChromSizes(gr)
	gr.Close()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error reading genome file: %v\n", err)
		os.Exit(1)
	}

	// Positional input file overrides default stdin.
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

	if _, err := bedcomplement.Complement(in, out, os.Stderr, sizes, order); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}
