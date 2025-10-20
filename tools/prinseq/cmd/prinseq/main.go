package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/tools/prinseq/pkg/prinseq"
)

const version = "1.0.0"

func main() {
	if len(os.Args) < 2 {
		printUsage()
		os.Exit(1)
	}

	command := os.Args[1]

	switch command {
	case "stats":
		runStats(os.Args[2:])
	case "filter":
		runFilter(os.Args[2:])
	case "version", "-v", "--version":
		fmt.Printf("prinseq version %s\n", version)
	case "help", "-h", "--help":
		printUsage()
	default:
		fmt.Fprintf(os.Stderr, "Unknown command: %s\n\n", command)
		printUsage()
		os.Exit(1)
	}
}

func printUsage() {
	fmt.Println(`prinseq - Sequence quality control and preprocessing tool

Usage:
  prinseq <command> [options]

Commands:
  stats     Calculate sequence statistics
  filter    Filter sequences based on quality criteria
  version   Print version information
  help      Print this help message

Use "prinseq <command> -h" for more information about a command.`)
}

func runStats(args []string) {
	fs := flag.NewFlagSet("stats", flag.ExitOnError)
	fastq := fs.String("fastq", "", "Input FASTQ file (use '-' for stdin)")
	fasta := fs.String("fasta", "", "Input FASTA file (use '-' for stdin)")

	fs.Usage = func() {
		fmt.Println(`Usage: prinseq stats [options]

Calculate sequence statistics for FASTA or FASTQ files.

Options:`)
		fs.PrintDefaults()
	}

	if err := fs.Parse(args); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	// Determine input file and format
	var inputFile string
	var isFastq bool

	if *fastq != "" {
		inputFile = *fastq
		isFastq = true
	} else if *fasta != "" {
		inputFile = *fasta
		isFastq = false
	} else {
		fmt.Fprintln(os.Stderr, "Error: Either -fastq or -fasta must be specified")
		fs.Usage()
		os.Exit(1)
	}

	// Open input file
	var reader *os.File
	if inputFile == "-" {
		reader = os.Stdin
	} else {
		var err error
		reader, err = os.Open(inputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening file: %v\n", err)
			os.Exit(1)
		}
		defer reader.Close()
	}

	// Calculate statistics
	stats, err := prinseq.CalculateStats(reader, isFastq)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error calculating statistics: %v\n", err)
		os.Exit(1)
	}

	// Print statistics
	fmt.Printf("Number of reads: %d\n", stats.NumReads)
	fmt.Printf("Total bases: %d\n", stats.TotalBases)
	fmt.Printf("Min length: %d\n", stats.MinLength)
	fmt.Printf("Max length: %d\n", stats.MaxLength)
	fmt.Printf("Average length: %.2f\n", stats.AvgLength)
	fmt.Printf("GC content: %.2f%%\n", stats.GCContent)
	fmt.Printf("Number of Ns: %d\n", stats.NumNs)
	if isFastq {
		fmt.Printf("Average quality: %.2f\n", stats.AvgQuality)
	}
}

func runFilter(args []string) {
	fs := flag.NewFlagSet("filter", flag.ExitOnError)
	
	// Input/output options
	fastq := fs.String("fastq", "", "Input FASTQ file (use '-' for stdin)")
	fasta := fs.String("fasta", "", "Input FASTA file (use '-' for stdin)")
	outGood := fs.String("out_good", "", "Output file for filtered sequences (default: stdout)")

	// Filter options
	minLen := fs.Int("min_len", 0, "Minimum sequence length")
	maxLen := fs.Int("max_len", 0, "Maximum sequence length")
	minGC := fs.Float64("min_gc", 0, "Minimum GC content percentage (0-100)")
	maxGC := fs.Float64("max_gc", 0, "Maximum GC content percentage (0-100)")
	minQualMean := fs.Float64("min_qual_mean", 0, "Minimum mean quality score")
	maxQualMean := fs.Float64("max_qual_mean", 0, "Maximum mean quality score")
	maxNsP := fs.Float64("ns_max_p", 0, "Maximum percentage of Ns allowed")
	maxNsN := fs.Int("ns_max_n", 0, "Maximum number of Ns allowed")

	fs.Usage = func() {
		fmt.Println(`Usage: prinseq filter [options]

Filter sequences based on quality criteria.

Options:`)
		fs.PrintDefaults()
	}

	if err := fs.Parse(args); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	// Determine input file and format
	var inputFile string
	var isFastq bool

	if *fastq != "" {
		inputFile = *fastq
		isFastq = true
	} else if *fasta != "" {
		inputFile = *fasta
		isFastq = false
	} else {
		fmt.Fprintln(os.Stderr, "Error: Either -fastq or -fasta must be specified")
		fs.Usage()
		os.Exit(1)
	}

	// Open input file
	var reader *os.File
	if inputFile == "-" {
		reader = os.Stdin
	} else {
		var err error
		reader, err = os.Open(inputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening file: %v\n", err)
			os.Exit(1)
		}
		defer reader.Close()
	}

	// Open output file
	var writer *os.File
	if *outGood == "" {
		writer = os.Stdout
	} else {
		var err error
		writer, err = os.Create(*outGood)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
		defer writer.Close()
	}

	// Set filter options
	opts := prinseq.FilterOptions{
		MinLen:      *minLen,
		MaxLen:      *maxLen,
		MinGC:       *minGC,
		MaxGC:       *maxGC,
		MinQualMean: *minQualMean,
		MaxQualMean: *maxQualMean,
		MaxNsP:      *maxNsP,
		MaxNsN:      *maxNsN,
	}

	// Filter sequences
	if err := prinseq.Filter(reader, writer, isFastq, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error filtering sequences: %v\n", err)
		os.Exit(1)
	}
}
