// bed12tobed6 splits each BED12 record into its individual blocks, emitting
// one BED6 record per block (mirrors `bedtools bed12tobed6`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bed12tobed6/pkg/bed12tobed6"
)

const version = "1.0.0"

const usage = `bed12tobed6 - Split BED12 records into per-block BED6 records

Usage:
  bed12tobed6 [options]

Description:
  Reads a BED12 file and writes one BED6 record per block. Each block becomes
  a 6-column record (chrom, start, end, name, score, strand). With -n the
  score column is populated with the 1-based block index; on '-' strand
  records the indices are reversed so the first emitted block carries the
  highest index (matches upstream bedtools).

Options:
  -i, --input FILE         Input BED12 file ('-' for stdin, default: stdin)
  -o, --output FILE        Output BED6 file ('-' for stdout, default: stdout)
  -n, --number             Number the blocks (1-based) into the score column.
  -h, --help               Show this help message and exit
  -v, --version            Show version information and exit
`

func main() {
	fs := flag.CommandLine
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	var inputFile, outputFile string
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input BED12 file")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output BED6 file")
	var numberBlocks bool
	cliflag.BoolVar(fs, &numberBlocks, "n", "number", false, "Number the blocks")
	var help, showVersion bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help message")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version information")

	flag.Parse()
	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if showVersion {
		fmt.Printf("bed12tobed6 version %s\n", version)
		os.Exit(0)
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

	if _, err := bed12tobed6.Convert(in, out, bed12tobed6.Options{NumberBlocks: numberBlocks}); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}
