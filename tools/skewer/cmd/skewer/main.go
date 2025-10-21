package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/skewer/pkg/skewer"
)

const usage = `skewer - Fast adapter trimming tool for FASTQ files

Usage:
  skewer <command> [options]

Commands:
  se    Trim single-end reads
  pe    Trim paired-end reads

For command-specific help:
  skewer se -h
  skewer pe -h

Version: 1.0.0 (Go implementation)
`

func main() {
	if len(os.Args) < 2 {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(1)
	}

	command := os.Args[1]
	switch command {
	case "se":
		runSingleEnd()
	case "pe":
		runPairedEnd()
	case "-h", "--help", "help":
		fmt.Print(usage)
		os.Exit(0)
	case "-v", "--version", "version":
		fmt.Println("skewer version 1.0.0 (Go implementation)")
		os.Exit(0)
	default:
		fmt.Fprintf(os.Stderr, "Error: unknown command %q\n\n", command)
		fmt.Fprint(os.Stderr, usage)
		os.Exit(1)
	}
}

func runSingleEnd() {
	fs := flag.NewFlagSet("skewer se", flag.ExitOnError)
	
	var (
		inputFile      string
		outputFile     string
		adapter3       string
		adapter5       string
		qualType       string
		minLength      int
		qualThreshold  int
		minOverlap     int
		errorRate      float64
		quiet          bool
	)
	
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input FASTQ file (required)")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output trimmed file (default: stdout)")
	cliflag.StringVar(fs, &adapter3, "x", "adapter3", "", "3' adapter sequence")
	cliflag.StringVar(fs, &adapter5, "y", "adapter5", "", "5' adapter sequence")
	cliflag.StringVar(fs, &qualType, "t", "qual-type", "sanger", "Quality type: sanger, illumina (default: sanger)")
	cliflag.IntVar(fs, &minLength, "l", "min-length", 18, "Minimum read length (default: 18)")
	cliflag.IntVar(fs, &qualThreshold, "q", "qual-threshold", 0, "Quality threshold for trimming (default: 0)")
	cliflag.IntVar(fs, &minOverlap, "m", "min-overlap", 3, "Minimum overlap for adapter detection (default: 3)")
	cliflag.Float64Var(fs, &errorRate, "r", "error-rate", 0.1, "Maximum error rate (default: 0.1)")
	cliflag.BoolVar(fs, &quiet, "", "quiet", false, "Don't print statistics")
	
	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: skewer se [options]\n\n")
		fmt.Fprintf(os.Stderr, "Options:\n")
		fmt.Fprintf(os.Stderr, "  -i, --input FILE          Input FASTQ file (required)\n")
		fmt.Fprintf(os.Stderr, "  -o, --output FILE         Output trimmed file (default: stdout)\n")
		fmt.Fprintf(os.Stderr, "  -x, --adapter3 SEQ        3' adapter sequence\n")
		fmt.Fprintf(os.Stderr, "  -y, --adapter5 SEQ        5' adapter sequence\n")
		fmt.Fprintf(os.Stderr, "  -t, --qual-type TYPE      Quality type: sanger, illumina (default: sanger)\n")
		fmt.Fprintf(os.Stderr, "  -l, --min-length INT      Minimum read length (default: 18)\n")
		fmt.Fprintf(os.Stderr, "  -q, --qual-threshold INT  Quality threshold for trimming (default: 0)\n")
		fmt.Fprintf(os.Stderr, "  -m, --min-overlap INT     Minimum overlap for adapter detection (default: 3)\n")
		fmt.Fprintf(os.Stderr, "  -r, --error-rate FLOAT    Maximum error rate (default: 0.1)\n")
		fmt.Fprintf(os.Stderr, "  --quiet                   Don't print statistics\n")
		fmt.Fprintf(os.Stderr, "\nExample:\n")
		fmt.Fprintf(os.Stderr, "  skewer se -i input.fastq -o output.fastq -x AGATCGGAAGAGC\n")
		fmt.Fprintf(os.Stderr, "  skewer se -i input.fastq.gz -o output.fastq.gz -x AGATCGGAAGAGC -l 20\n")
	}
	
	fs.Parse(os.Args[2:])
	
	// Validate required arguments
	if inputFile == "" {
		fmt.Fprintln(os.Stderr, "Error: -i/--input is required")
		fs.Usage()
		os.Exit(1)
	}
	
	// Determine quality encoding
	encoding := getQualityEncoding(qualType)
	
	// Open input file (with automatic gzip support)
	input, err := iohelper.OpenReader(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()
	
	// Open output file (with automatic gzip support)
	outFileName := outputFile
	if outFileName == "" {
		outFileName = "-"
	}
	output, err := iohelper.OpenWriter(outFileName)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer output.Close()
	
	// Set up trim options
	opts := skewer.TrimOptions{
		Adapter3:      adapter3,
		Adapter5:      adapter5,
		MinLength:     minLength,
		QualThreshold: qualThreshold,
		MinOverlap:    minOverlap,
		ErrorRate:     errorRate,
	}
	
	// Perform trimming
	stats, err := skewer.TrimSingleEnd(input, output, encoding, opts)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error during trimming: %v\n", err)
		os.Exit(1)
	}
	
	// Print statistics
	if !quiet {
		printStats(stats, "SE")
	}
}

func runPairedEnd() {
	fs := flag.NewFlagSet("skewer pe", flag.ExitOnError)
	
	var (
		inputFile1     string
		inputFile2     string
		outputFile1    string
		outputFile2    string
		outputSingle   string
		adapter3       string
		adapter5       string
		qualType       string
		minLength      int
		qualThreshold  int
		minOverlap     int
		errorRate      float64
		quiet          bool
	)
	
	cliflag.StringVar(fs, &inputFile1, "i", "input1", "", "First input FASTQ file (required)")
	cliflag.StringVar(fs, &inputFile2, "j", "input2", "", "Second input FASTQ file (required)")
	cliflag.StringVar(fs, &outputFile1, "o", "output1", "", "First output trimmed file (required)")
	cliflag.StringVar(fs, &outputFile2, "p", "output2", "", "Second output trimmed file (required)")
	cliflag.StringVar(fs, &outputSingle, "s", "single", "", "Output single-end reads (optional)")
	cliflag.StringVar(fs, &adapter3, "x", "adapter3", "", "3' adapter sequence")
	cliflag.StringVar(fs, &adapter5, "y", "adapter5", "", "5' adapter sequence")
	cliflag.StringVar(fs, &qualType, "t", "qual-type", "sanger", "Quality type: sanger, illumina (default: sanger)")
	cliflag.IntVar(fs, &minLength, "l", "min-length", 18, "Minimum read length (default: 18)")
	cliflag.IntVar(fs, &qualThreshold, "q", "qual-threshold", 0, "Quality threshold for trimming (default: 0)")
	cliflag.IntVar(fs, &minOverlap, "m", "min-overlap", 3, "Minimum overlap for adapter detection (default: 3)")
	cliflag.Float64Var(fs, &errorRate, "r", "error-rate", 0.1, "Maximum error rate (default: 0.1)")
	cliflag.BoolVar(fs, &quiet, "", "quiet", false, "Don't print statistics")
	
	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: skewer pe [options]\n\n")
		fmt.Fprintf(os.Stderr, "Options:\n")
		fmt.Fprintf(os.Stderr, "  -i, --input1 FILE         First input FASTQ file (required)\n")
		fmt.Fprintf(os.Stderr, "  -j, --input2 FILE         Second input FASTQ file (required)\n")
		fmt.Fprintf(os.Stderr, "  -o, --output1 FILE        First output trimmed file (required)\n")
		fmt.Fprintf(os.Stderr, "  -p, --output2 FILE        Second output trimmed file (required)\n")
		fmt.Fprintf(os.Stderr, "  -s, --single FILE         Output single-end reads (optional)\n")
		fmt.Fprintf(os.Stderr, "  -x, --adapter3 SEQ        3' adapter sequence\n")
		fmt.Fprintf(os.Stderr, "  -y, --adapter5 SEQ        5' adapter sequence\n")
		fmt.Fprintf(os.Stderr, "  -t, --qual-type TYPE      Quality type: sanger, illumina (default: sanger)\n")
		fmt.Fprintf(os.Stderr, "  -l, --min-length INT      Minimum read length (default: 18)\n")
		fmt.Fprintf(os.Stderr, "  -q, --qual-threshold INT  Quality threshold for trimming (default: 0)\n")
		fmt.Fprintf(os.Stderr, "  -m, --min-overlap INT     Minimum overlap for adapter detection (default: 3)\n")
		fmt.Fprintf(os.Stderr, "  -r, --error-rate FLOAT    Maximum error rate (default: 0.1)\n")
		fmt.Fprintf(os.Stderr, "  --quiet                   Don't print statistics\n")
		fmt.Fprintf(os.Stderr, "\nExample:\n")
		fmt.Fprintf(os.Stderr, "  skewer pe -i R1.fastq -j R2.fastq -o R1_trim.fastq -p R2_trim.fastq -x AGATCGGAAGAGC\n")
	}
	
	fs.Parse(os.Args[2:])
	
	// Validate required arguments
	if inputFile1 == "" || inputFile2 == "" {
		fmt.Fprintln(os.Stderr, "Error: both -i/--input1 and -j/--input2 are required")
		fs.Usage()
		os.Exit(1)
	}
	if outputFile1 == "" || outputFile2 == "" {
		fmt.Fprintln(os.Stderr, "Error: both -o/--output1 and -p/--output2 are required")
		fs.Usage()
		os.Exit(1)
	}
	
	// Determine quality encoding
	encoding := getQualityEncoding(qualType)
	
	// Open input files (with automatic gzip support)
	input1, err := iohelper.OpenReader(inputFile1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening first input file: %v\n", err)
		os.Exit(1)
	}
	defer input1.Close()
	
	input2, err := iohelper.OpenReader(inputFile2)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening second input file: %v\n", err)
		os.Exit(1)
	}
	defer input2.Close()
	
	// Open output files (with automatic gzip support)
	output1, err := iohelper.OpenWriter(outputFile1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating first output file: %v\n", err)
		os.Exit(1)
	}
	defer output1.Close()
	
	output2, err := iohelper.OpenWriter(outputFile2)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating second output file: %v\n", err)
		os.Exit(1)
	}
	defer output2.Close()
	
	// Open optional single output file (with automatic gzip support)
	var outSingle io.Writer
	if outputSingle != "" {
		f, err := iohelper.OpenWriter(outputSingle)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating single output file: %v\n", err)
			os.Exit(1)
		}
		defer f.Close()
		outSingle = f
	}
	
	// Set up trim options
	opts := skewer.TrimOptions{
		Adapter3:      adapter3,
		Adapter5:      adapter5,
		MinLength:     minLength,
		QualThreshold: qualThreshold,
		MinOverlap:    minOverlap,
		ErrorRate:     errorRate,
	}
	
	// Perform trimming
	stats, err := skewer.TrimPairedEnd(input1, input2, output1, output2, outSingle, encoding, opts)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error during trimming: %v\n", err)
		os.Exit(1)
	}
	
	// Print statistics
	if !quiet {
		printStats(stats, "PE")
	}
}

func getQualityEncoding(qualType string) fastq.QualityEncoding {
	switch qualType {
	case "sanger", "phred33":
		return fastq.Phred33
	case "illumina", "phred64":
		return fastq.Phred64
	default:
		fmt.Fprintf(os.Stderr, "Warning: unknown quality type %q, using sanger\n", qualType)
		return fastq.Phred33
	}
}

func printStats(stats *skewer.TrimStats, mode string) {
	fmt.Fprintf(os.Stderr, "\n%s Adapter Trimming Stats:\n", mode)
	fmt.Fprintf(os.Stderr, "  Total reads:        %d\n", stats.TotalReads)
	fmt.Fprintf(os.Stderr, "  Trimmed reads:      %d (%.2f%%)\n", 
		stats.TrimmedReads, 
		100.0*float64(stats.TrimmedReads)/float64(stats.TotalReads))
	fmt.Fprintf(os.Stderr, "  3' adapters found:  %d (%.2f%%)\n", 
		stats.AdapterFound3, 
		100.0*float64(stats.AdapterFound3)/float64(stats.TotalReads))
	fmt.Fprintf(os.Stderr, "  5' adapters found:  %d (%.2f%%)\n", 
		stats.AdapterFound5, 
		100.0*float64(stats.AdapterFound5)/float64(stats.TotalReads))
	fmt.Fprintf(os.Stderr, "  Discarded reads:    %d (%.2f%%)\n", 
		stats.DiscardedReads, 
		100.0*float64(stats.DiscardedReads)/float64(stats.TotalReads))
	fmt.Fprintf(os.Stderr, "  Kept reads:         %d (%.2f%%)\n", 
		stats.TotalReads-stats.DiscardedReads, 
		100.0*float64(stats.TotalReads-stats.DiscardedReads)/float64(stats.TotalReads))
	fmt.Fprintf(os.Stderr, "  Total bases:        %d\n", stats.TotalBases)
	fmt.Fprintf(os.Stderr, "  Trimmed bases:      %d (%.2f%%)\n", 
		stats.TrimmedBases, 
		100.0*float64(stats.TrimmedBases)/float64(stats.TotalBases))
	fmt.Fprintln(os.Stderr)
}
