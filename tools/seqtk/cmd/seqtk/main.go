// seqtk: A fast FASTA/Q processor implemented in Go
//
// This is a Go reimplementation of the popular seqtk tool with enhanced performance
// and better error handling. It provides common operations on FASTA and FASTQ files.
//
// Usage:
//   seqtk <command> [options] <input>
//
// Commands:
//   seq        Transform sequences (reverse complement, etc.)
//   sample     Subsample sequences
//   trimfq     Trim FASTQ sequences based on quality
//   fq2fa      Convert FASTQ to FASTA
//   comp       Get sequence composition statistics
//
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
	"github.com/yassineS/bio_ai_experiment/tools/seqtk/pkg/seqtk"
)

const version = "1.0.0"

func main() {
	if len(os.Args) < 2 {
		printUsage()
		os.Exit(1)
	}

	command := os.Args[1]

	switch command {
	case "seq":
		seqCommand()
	case "sample":
		sampleCommand()
	case "trimfq":
		trimfqCommand()
	case "fq2fa":
		fq2faCommand()
	case "comp":
		compCommand()
	case "version", "-v", "--version":
		fmt.Printf("seqtk version %s\n", version)
	case "help", "-h", "--help":
		printUsage()
	default:
		fmt.Fprintf(os.Stderr, "Unknown command: %s\n", command)
		printUsage()
		os.Exit(1)
	}
}

func printUsage() {
	fmt.Fprintf(os.Stderr, `seqtk - A fast FASTA/Q processor (Go implementation)

Usage: seqtk <command> [options] <input>

Commands:
  seq        Transform sequences (reverse complement, etc.)
  sample     Subsample sequences
  trimfq     Trim FASTQ sequences based on quality
  fq2fa      Convert FASTQ to FASTA
  comp       Get sequence composition statistics
  version    Show version information
  help       Show this help message

Use 'seqtk <command> -h' for help on a specific command.

Examples:
  seqtk comp sequences.fasta
  seqtk fq2fa reads.fastq > reads.fasta
  seqtk seq -r sequences.fasta > rev_comp.fasta
  seqtk sample reads.fastq 0.1 > sample.fastq
  seqtk trimfq -q 20 reads.fastq > trimmed.fastq

`)
}

func seqCommand() {
	fs := flag.NewFlagSet("seq", flag.ExitOnError)
	revComp := fs.Bool("r", false, "reverse complement")
	phred64 := fs.Bool("6", false, "use Phred+64 quality encoding (default: Phred+33)")
	output := fs.String("o", "", "output file (default: stdout)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk seq [options] <input>

Transform sequences.

Options:
  -r         Reverse complement
  -6         Use Phred+64 quality encoding (default: Phred+33)
  -o FILE    Output file (default: stdout)

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 1 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)
	input, err := os.Open(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	// Determine output
	var out *os.File
	if *output != "" {
		out, err = os.Create(*output)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
		defer out.Close()
	} else {
		out = os.Stdout
	}

	// Detect file type
	isFastq, err := seqtk.GetFileType(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error detecting file type: %v\n", err)
		os.Exit(1)
	}

	encoding := fastq.Phred33
	if *phred64 {
		encoding = fastq.Phred64
	}

	// Reopen file for processing
	input.Close()
	input, err = os.Open(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}

	if *revComp {
		if err := seqtk.ReverseComplement(input, out, isFastq, encoding); err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v\n", err)
			os.Exit(1)
		}
	} else {
		// Just copy through (could add more transformations here)
		fmt.Fprintf(os.Stderr, "No transformation specified. Use -r for reverse complement.\n")
		os.Exit(1)
	}
}

func sampleCommand() {
	fs := flag.NewFlagSet("sample", flag.ExitOnError)
	phred64 := fs.Bool("6", false, "use Phred+64 quality encoding (default: Phred+33)")
	output := fs.String("o", "", "output file (default: stdout)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk sample [options] <input> <fraction>

Subsample sequences randomly.

Arguments:
  <input>      Input FASTA/FASTQ file
  <fraction>   Fraction of sequences to sample (0.0-1.0)

Options:
  -6         Use Phred+64 quality encoding (default: Phred+33)
  -o FILE    Output file (default: stdout)

Example:
  seqtk sample reads.fastq 0.1 > sample.fastq

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 2 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)
	var fraction float64
	if _, err := fmt.Sscanf(fs.Arg(1), "%f", &fraction); err != nil {
		fmt.Fprintf(os.Stderr, "Invalid fraction: %v\n", err)
		os.Exit(1)
	}

	input, err := os.Open(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	var out *os.File
	if *output != "" {
		out, err = os.Create(*output)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
		defer out.Close()
	} else {
		out = os.Stdout
	}

	isFastq, err := seqtk.GetFileType(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error detecting file type: %v\n", err)
		os.Exit(1)
	}

	encoding := fastq.Phred33
	if *phred64 {
		encoding = fastq.Phred64
	}

	// Reopen for processing
	input.Close()
	input, err = os.Open(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}

	if err := seqtk.Sample(input, out, fraction, isFastq, encoding); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func trimfqCommand() {
	fs := flag.NewFlagSet("trimfq", flag.ExitOnError)
	quality := fs.Int("q", 20, "minimum quality threshold")
	phred64 := fs.Bool("6", false, "use Phred+64 quality encoding (default: Phred+33)")
	output := fs.String("o", "", "output file (default: stdout)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk trimfq [options] <input>

Trim FASTQ sequences based on quality.

Options:
  -q INT     Minimum quality threshold (default: 20)
  -6         Use Phred+64 quality encoding (default: Phred+33)
  -o FILE    Output file (default: stdout)

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 1 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)
	input, err := os.Open(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	var out *os.File
	if *output != "" {
		out, err = os.Create(*output)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
		defer out.Close()
	} else {
		out = os.Stdout
	}

	encoding := fastq.Phred33
	if *phred64 {
		encoding = fastq.Phred64
	}

	if err := seqtk.TrimQuality(input, out, *quality, encoding); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func fq2faCommand() {
	fs := flag.NewFlagSet("fq2fa", flag.ExitOnError)
	phred64 := fs.Bool("6", false, "use Phred+64 quality encoding (default: Phred+33)")
	output := fs.String("o", "", "output file (default: stdout)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk fq2fa [options] <input>

Convert FASTQ to FASTA.

Options:
  -6         Use Phred+64 quality encoding (default: Phred+33)
  -o FILE    Output file (default: stdout)

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 1 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)
	input, err := os.Open(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	var out *os.File
	if *output != "" {
		out, err = os.Create(*output)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
		defer out.Close()
	} else {
		out = os.Stdout
	}

	encoding := fastq.Phred33
	if *phred64 {
		encoding = fastq.Phred64
	}

	if err := seqtk.ConvertFastqToFasta(input, out, encoding); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func compCommand() {
	fs := flag.NewFlagSet("comp", flag.ExitOnError)
	phred64 := fs.Bool("6", false, "use Phred+64 quality encoding (default: Phred+33)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk comp [options] <input>

Get sequence composition statistics.

Options:
  -6         Use Phred+64 quality encoding for FASTQ (default: Phred+33)

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 1 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)

	isFastq, err := seqtk.GetFileType(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error detecting file type: %v\n", err)
		os.Exit(1)
	}

	// Open file fresh for reading
	input, err := os.Open(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	var stats *seqtk.Stats
	if isFastq {
		encoding := fastq.Phred33
		if *phred64 {
			encoding = fastq.Phred64
		}
		stats, err = seqtk.CalculateFastqStats(input, encoding)
	} else {
		stats, err = seqtk.CalculateFastaStats(input)
	}

	if err != nil {
		fmt.Fprintf(os.Stderr, "Error calculating statistics: %v\n", err)
		os.Exit(1)
	}

	fmt.Printf("Number of sequences: %d\n", stats.NumSequences)
	fmt.Printf("Total bases: %d\n", stats.TotalBases)
	fmt.Printf("Min length: %d\n", stats.MinLength)
	fmt.Printf("Max length: %d\n", stats.MaxLength)
	fmt.Printf("Average length: %.2f\n", stats.AvgLength)
	fmt.Printf("GC content: %.2f%%\n", stats.GCContent)
	if isFastq {
		fmt.Printf("Average quality: %.2f\n", stats.AvgQuality)
	}
}
