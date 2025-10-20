package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/tools/prinseq/pkg/prinseq"
)

const version = "1.0.0"

// openInput opens the input file or returns stdin if filename is "-"
func openInput(filename string) (io.ReadCloser, error) {
	if filename == "-" {
		return io.NopCloser(os.Stdin), nil
	}
	return os.Open(filename)
}

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
	stats, err2 := prinseq.CalculateStats(reader, isFastq)
	if err2 != nil {
		fmt.Fprintf(os.Stderr, "Error calculating statistics: %v\n", err2)
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
	fastq2 := fs.String("fastq2", "", "Input paired-end FASTQ file 2")
	fasta2 := fs.String("fasta2", "", "Input paired-end FASTA file 2")
	outGood := fs.String("out_good", "", "Output file for filtered sequences (default: stdout)")
	outGood2 := fs.String("out_good2", "", "Output file for paired-end file 2 (required with fastq2/fasta2)")

	// Filter options
	minLen := fs.Int("min_len", 0, "Minimum sequence length")
	maxLen := fs.Int("max_len", 0, "Maximum sequence length")
	minGC := fs.Float64("min_gc", 0, "Minimum GC content percentage (0-100)")
	maxGC := fs.Float64("max_gc", 0, "Maximum GC content percentage (0-100)")
	minQualMean := fs.Float64("min_qual_mean", 0, "Minimum mean quality score")
	maxQualMean := fs.Float64("max_qual_mean", 0, "Maximum mean quality score")
	maxNsP := fs.Float64("ns_max_p", 0, "Maximum percentage of Ns allowed")
	maxNsN := fs.Int("ns_max_n", 0, "Maximum number of Ns allowed")

	// Trimming options
	trimLeft := fs.Int("trim_left", 0, "Trim bases from 5' end")
	trimRight := fs.Int("trim_right", 0, "Trim bases from 3' end")
	trimLeftP := fs.Int("trim_left_p", 0, "Trim percentage from 5' end")
	trimRightP := fs.Int("trim_right_p", 0, "Trim percentage from 3' end")
	trimQualL := fs.Int("trim_qual_left", 0, "Trim 5' end by quality threshold")
	trimQualR := fs.Int("trim_qual_right", 0, "Trim 3' end by quality threshold")
	trimNsLeft := fs.Int("trim_ns_left", 0, "Trim poly-N tail from 5' end (min length)")
	trimNsRight := fs.Int("trim_ns_right", 0, "Trim poly-N tail from 3' end (min length)")
	trimTailLeft := fs.Int("trim_tail_left", 0, "Trim poly-A/T tail from 5' end (min length)")
	trimTailRight := fs.Int("trim_tail_right", 0, "Trim poly-A/T tail from 3' end (min length)")

	// Duplicate removal options
	derep := fs.Int("derep", 0, "Remove duplicates: 1=exact, 4=reverse complement, 5=both")
	derepMin := fs.Int("derep_min", 2, "Minimum number of duplicates to keep")

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

	// Check for paired-end mode
	isPaired := (*fastq2 != "" || *fasta2 != "")

	// Determine input file and format
	var inputFile, inputFile2 string
	var isFastq bool

	if *fastq != "" {
		inputFile = *fastq
		inputFile2 = *fastq2
		isFastq = true
	} else if *fasta != "" {
		inputFile = *fasta
		inputFile2 = *fasta2
		isFastq = false
	} else {
		fmt.Fprintln(os.Stderr, "Error: Either -fastq or -fasta must be specified")
		fs.Usage()
		os.Exit(1)
	}

	// Validate paired-end requirements
	if isPaired && *outGood2 == "" {
		fmt.Fprintln(os.Stderr, "Error: -out_good2 required when using paired-end input")
		os.Exit(1)
	}

	// Set filter options
	opts := prinseq.FilterOptions{
		MinLen:        *minLen,
		MaxLen:        *maxLen,
		MinGC:         *minGC,
		MaxGC:         *maxGC,
		MinQualMean:   *minQualMean,
		MaxQualMean:   *maxQualMean,
		MaxNsP:        *maxNsP,
		MaxNsN:        *maxNsN,
		TrimLeft:      *trimLeft,
		TrimRight:     *trimRight,
		TrimLeftP:     *trimLeftP,
		TrimRightP:    *trimRightP,
		TrimQualL:     *trimQualL,
		TrimQualR:     *trimQualR,
		TrimNsLeft:    *trimNsLeft,
		TrimNsRight:   *trimNsRight,
		TrimTailLeft:  *trimTailLeft,
		TrimTailRight: *trimTailRight,
		Derep:         *derep,
		DerepMin:      *derepMin,
	}

	if isPaired {
		// Paired-end mode
		reader1, err := openInput(inputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening file 1: %v\n", err)
			os.Exit(1)
		}
		defer reader1.Close()

		reader2, err := openInput(inputFile2)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening file 2: %v\n", err)
			os.Exit(1)
		}
		defer reader2.Close()

		var writer1, writer2 io.WriteCloser
		if *outGood == "" {
			writer1 = os.Stdout
		} else {
			writer1, err = os.Create(*outGood)
			if err != nil {
				fmt.Fprintf(os.Stderr, "Error creating output file 1: %v\n", err)
				os.Exit(1)
			}
			defer writer1.Close()
		}

		writer2, err = os.Create(*outGood2)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file 2: %v\n", err)
			os.Exit(1)
		}
		defer writer2.Close()

		// Filter paired sequences
		if err := prinseq.FilterPaired(reader1, reader2, writer1, writer2, isFastq, opts); err != nil {
			fmt.Fprintf(os.Stderr, "Error filtering paired sequences: %v\n", err)
			os.Exit(1)
		}
	} else {
		// Single-end mode
		reader, err := openInput(inputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening file: %v\n", err)
			os.Exit(1)
		}
		defer reader.Close()

		var writer io.WriteCloser
		if *outGood == "" {
			writer = os.Stdout
		} else {
			writer, err = os.Create(*outGood)
			if err != nil {
				fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
				os.Exit(1)
			}
			defer writer.Close()
		}

		// Filter sequences
		if err2 := prinseq.Filter(reader, writer, isFastq, opts); err2 != nil {
			fmt.Fprintf(os.Stderr, "Error filtering sequences: %v\n", err2)
			os.Exit(1)
		}
	}
}
