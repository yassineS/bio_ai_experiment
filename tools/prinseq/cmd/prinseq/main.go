package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
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
	
	var fastq, fasta string
	cliflag.StringVar(fs, &fastq, "", "fastq", "", "Input FASTQ file (use '-' for stdin)")
	cliflag.StringVar(fs, &fasta, "", "fasta", "", "Input FASTA file (use '-' for stdin)")

	fs.Usage = func() {
		fmt.Print(`Usage: prinseq stats [options]

Calculate sequence statistics for FASTA or FASTQ files.

Options:
  --fasta FILE              Input FASTA file (use '-' for stdin)
  --fastq FILE              Input FASTQ file (use '-' for stdin)
`)
	}

	if err := fs.Parse(args); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	// Determine input file and format
	var inputFile string
	var isFastq bool

	if fastq != "" {
		inputFile = fastq
		isFastq = true
	} else if fasta != "" {
		inputFile = fasta
		isFastq = false
	} else {
		fmt.Fprintln(os.Stderr, "Error: Either --fastq or --fasta must be specified")
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
	var input1, input2, output1, output2 string
	var fasta, fastq bool
	
	cliflag.StringVar(fs, &input1, "i", "input", "", "Primary input file (use '-' for stdin)")
	cliflag.StringVar(fs, &input2, "", "input2", "", "Paired-end input file 2")
	cliflag.StringVar(fs, &output1, "o", "output", "", "Output file for filtered sequences (default: stdout)")
	cliflag.StringVar(fs, &output2, "", "output2", "", "Output file for paired-end file 2")
	cliflag.BoolVar(fs, &fasta, "", "fasta", false, "Input is FASTA format")
	cliflag.BoolVar(fs, &fastq, "", "fastq", false, "Input is FASTQ format")

	// Filter options
	var minLen, maxLen, maxNsN int
	var minGC, maxGC, minQualMean, maxQualMean, maxNsP float64
	
	cliflag.IntVar(fs, &minLen, "l", "min-length", 0, "Minimum sequence length")
	cliflag.IntVar(fs, &maxLen, "L", "max-length", 0, "Maximum sequence length")
	cliflag.Float64Var(fs, &minGC, "g", "min-gc", 0, "Minimum GC content percentage (0-100)")
	cliflag.Float64Var(fs, &maxGC, "G", "max-gc", 0, "Maximum GC content percentage (0-100)")
	cliflag.Float64Var(fs, &minQualMean, "q", "min-quality", 0, "Minimum mean quality score")
	cliflag.Float64Var(fs, &maxQualMean, "Q", "max-quality", 0, "Maximum mean quality score")
	cliflag.Float64Var(fs, &maxNsP, "N", "max-ns-percent", 0, "Maximum percentage of Ns allowed")
	cliflag.IntVar(fs, &maxNsN, "n", "max-ns", 0, "Maximum number of Ns allowed")

	// Trimming options
	var trimLeft, trimRight, trimLeftP, trimRightP int
	var trimQualL, trimQualR, trimNsLeft, trimNsRight int
	var trimTailLeft, trimTailRight int
	
	cliflag.IntVar(fs, &trimLeft, "", "trim-left", 0, "Trim bases from 5' end")
	cliflag.IntVar(fs, &trimRight, "", "trim-right", 0, "Trim bases from 3' end")
	cliflag.IntVar(fs, &trimLeftP, "", "trim-left-p", 0, "Trim percentage from 5' end")
	cliflag.IntVar(fs, &trimRightP, "", "trim-right-p", 0, "Trim percentage from 3' end")
	cliflag.IntVar(fs, &trimQualL, "", "trim-qual-left", 0, "Quality threshold for 5' trimming")
	cliflag.IntVar(fs, &trimQualR, "", "trim-qual-right", 0, "Quality threshold for 3' trimming")
	cliflag.IntVar(fs, &trimNsLeft, "", "trim-n-left", 0, "Trim poly-N from 5' end (min length)")
	cliflag.IntVar(fs, &trimNsRight, "", "trim-n-right", 0, "Trim poly-N from 3' end (min length)")
	cliflag.IntVar(fs, &trimTailLeft, "", "trim-tail-left", 0, "Trim poly-A/T from 5' end (min length)")
	cliflag.IntVar(fs, &trimTailRight, "", "trim-tail-right", 0, "Trim poly-A/T from 3' end (min length)")

	// Duplicate removal options
	var derep, derepMin int
	cliflag.IntVar(fs, &derep, "d", "derep", 0, "Remove duplicates: 1=exact, 4=revcomp, 5=both")
	cliflag.IntVar(fs, &derepMin, "", "derep-min", 2, "Minimum occurrences to keep")

	fs.Usage = func() {
		fmt.Print(`Usage: prinseq filter [options]

Filter sequences based on quality criteria.

Input/Output Options:
  -i, --input FILE          Primary input file (use '-' for stdin)
  --input2 FILE             Paired-end input file 2
  -o, --output FILE         Output file (default: stdout)
  --output2 FILE            Output file for paired-end file 2
  --fasta                   Input is FASTA format
  --fastq                   Input is FASTQ format

Filter Options:
  -l, --min-length INT      Minimum sequence length
  -L, --max-length INT      Maximum sequence length
  -g, --min-gc FLOAT        Minimum GC content (0-100)
  -G, --max-gc FLOAT        Maximum GC content (0-100)
  -q, --min-quality FLOAT   Minimum mean quality score
  -Q, --max-quality FLOAT   Maximum mean quality score
  -n, --max-ns INT          Maximum number of Ns
  -N, --max-ns-percent FLOAT Maximum percentage of Ns (0-100)

Trimming Options:
  --trim-left INT           Trim bases from 5' end
  --trim-right INT          Trim bases from 3' end
  --trim-left-p INT         Trim percentage from 5' end
  --trim-right-p INT        Trim percentage from 3' end
  --trim-qual-left INT      Quality threshold for 5' trimming
  --trim-qual-right INT     Quality threshold for 3' trimming
  --trim-n-left INT         Trim poly-N from 5' end
  --trim-n-right INT        Trim poly-N from 3' end
  --trim-tail-left INT      Trim poly-A/T from 5' end
  --trim-tail-right INT     Trim poly-A/T from 3' end

Duplicate Removal Options:
  -d, --derep MODE          Remove duplicates (1=exact, 4=revcomp, 5=both)
  --derep-min INT           Minimum occurrences to keep (default: 2)

Examples:
  # Filter by length using short options
  prinseq filter -i reads.fastq -o filtered.fastq -l 100 -L 500

  # Filter with quality and GC using long options
  prinseq filter --input reads.fastq --min-quality 20 --min-gc 40 --max-gc 60

  # Trim and filter
  prinseq filter -i reads.fastq --trim-qual-left 20 --trim-qual-right 20 -l 100

  # Paired-end filtering
  prinseq filter -i R1.fastq --input2 R2.fastq -o out_R1.fastq --output2 out_R2.fastq

  # Remove duplicates
  prinseq filter -i seqs.fasta -d 1 --derep-min 2 -o unique.fasta
`)
	}

	if err := fs.Parse(args); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	// Determine format
	var isFastq bool
	if fastq && fasta {
		fmt.Fprintln(os.Stderr, "Error: Cannot specify both --fasta and --fastq")
		os.Exit(1)
	} else if fastq {
		isFastq = true
	} else if fasta {
		isFastq = false
	} else {
		// Auto-detect from input filename
		if input1 != "" && input1 != "-" {
			if hasSuffix(input1, ".fastq", ".fq", ".fastq.gz", ".fq.gz") {
				isFastq = true
			} else if hasSuffix(input1, ".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz") {
				isFastq = false
			} else {
				fmt.Fprintln(os.Stderr, "Error: Cannot auto-detect format. Use --fasta or --fastq")
				os.Exit(1)
			}
		} else {
			fmt.Fprintln(os.Stderr, "Error: Must specify --fasta or --fastq")
			os.Exit(1)
		}
	}

	// Check for paired-end mode
	isPaired := (input2 != "")

	// Validate inputs
	if input1 == "" {
		fmt.Fprintln(os.Stderr, "Error: Input file required (-i or --input)")
		fs.Usage()
		os.Exit(1)
	}

	// Validate paired-end requirements
	if isPaired && output2 == "" {
		fmt.Fprintln(os.Stderr, "Error: --output2 required when using paired-end input")
		os.Exit(1)
	}

	// Set filter options
	opts := prinseq.FilterOptions{
		MinLen:        minLen,
		MaxLen:        maxLen,
		MinGC:         minGC,
		MaxGC:         maxGC,
		MinQualMean:   minQualMean,
		MaxQualMean:   maxQualMean,
		MaxNsP:        maxNsP,
		MaxNsN:        maxNsN,
		TrimLeft:      trimLeft,
		TrimRight:     trimRight,
		TrimLeftP:     trimLeftP,
		TrimRightP:    trimRightP,
		TrimQualL:     trimQualL,
		TrimQualR:     trimQualR,
		TrimNsLeft:    trimNsLeft,
		TrimNsRight:   trimNsRight,
		TrimTailLeft:  trimTailLeft,
		TrimTailRight: trimTailRight,
		Derep:         derep,
		DerepMin:      derepMin,
	}

	if isPaired {
		// Paired-end mode
		reader1, err := openInput(input1)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening file 1: %v\n", err)
			os.Exit(1)
		}
		defer reader1.Close()

		reader2, err := openInput(input2)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening file 2: %v\n", err)
			os.Exit(1)
		}
		defer reader2.Close()

		var writer1, writer2 io.WriteCloser
		if output1 == "" {
			writer1 = os.Stdout
		} else {
			writer1, err = os.Create(output1)
			if err != nil {
				fmt.Fprintf(os.Stderr, "Error creating output file 1: %v\n", err)
				os.Exit(1)
			}
			defer writer1.Close()
		}

		writer2, err = os.Create(output2)
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
		reader, err := openInput(input1)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening file: %v\n", err)
			os.Exit(1)
		}
		defer reader.Close()

		var writer io.WriteCloser
		if output1 == "" {
			writer = os.Stdout
		} else {
			writer, err = os.Create(output1)
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

// hasSuffix checks if a string ends with any of the given suffixes
func hasSuffix(s string, suffixes ...string) bool {
	for _, suffix := range suffixes {
		if len(s) >= len(suffix) && s[len(s)-len(suffix):] == suffix {
			return true
		}
	}
	return false
}
