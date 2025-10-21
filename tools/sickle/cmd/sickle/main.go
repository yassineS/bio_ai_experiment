package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/sickle/pkg/sickle"
)

const usage = `sickle - A windowed adaptive trimming tool for FASTQ files using quality scores

Usage:
  sickle <command> [options]

Commands:
  se    Trim single-end reads
  pe    Trim paired-end reads

For command-specific help:
  sickle se -h
  sickle pe -h

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
		fmt.Println("sickle version 1.0.0 (Go implementation)")
		os.Exit(0)
	default:
		fmt.Fprintf(os.Stderr, "Error: unknown command %q\n\n", command)
		fmt.Fprint(os.Stderr, usage)
		os.Exit(1)
	}
}

func runSingleEnd() {
	fs := flag.NewFlagSet("sickle se", flag.ExitOnError)
	
	var (
		fastqFile      string
		outputFile     string
		qualType       string
		qualThreshold  int
		lengthThreshold int
		noFivePrime    bool
		truncateN      bool
		quiet          bool
	)
	
	cliflag.StringVar(fs, &fastqFile, "f", "fastq-file", "", "Input FASTQ file (required)")
	cliflag.StringVar(fs, &outputFile, "o", "output-file", "", "Output trimmed file (default: stdout)")
	cliflag.StringVar(fs, &qualType, "t", "qual-type", "sanger", "Quality type: sanger, illumina, solexa (default: sanger)")
	cliflag.IntVar(fs, &qualThreshold, "q", "qual-threshold", 20, "Threshold for trimming (default: 20)")
	cliflag.IntVar(fs, &lengthThreshold, "l", "length-threshold", 20, "Minimum length to keep (default: 20)")
	cliflag.BoolVar(fs, &noFivePrime, "x", "no-fiveprime", false, "Don't trim 5' end")
	cliflag.BoolVar(fs, &truncateN, "n", "trunc-n", false, "Truncate sequences at position of first N")
	cliflag.BoolVar(fs, &quiet, "", "quiet", false, "Don't print statistics")
	
	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: sickle se [options]\n\n")
		fmt.Fprintf(os.Stderr, "Options:\n")
		fmt.Fprintf(os.Stderr, "  -f, --fastq-file FILE       Input FASTQ file (required)\n")
		fmt.Fprintf(os.Stderr, "  -o, --output-file FILE      Output trimmed file (default: stdout)\n")
		fmt.Fprintf(os.Stderr, "  -t, --qual-type TYPE        Quality type: sanger, illumina, solexa (default: sanger)\n")
		fmt.Fprintf(os.Stderr, "  -q, --qual-threshold INT    Threshold for trimming (default: 20)\n")
		fmt.Fprintf(os.Stderr, "  -l, --length-threshold INT  Minimum length to keep (default: 20)\n")
		fmt.Fprintf(os.Stderr, "  -x, --no-fiveprime          Don't trim 5' end\n")
		fmt.Fprintf(os.Stderr, "  -n, --trunc-n               Truncate sequences at position of first N\n")
		fmt.Fprintf(os.Stderr, "  --quiet                     Don't print statistics\n")
		fmt.Fprintf(os.Stderr, "\nExample:\n")
		fmt.Fprintf(os.Stderr, "  sickle se -f input.fastq -o output.fastq -q 20 -l 20\n")
	}
	
	fs.Parse(os.Args[2:])
	
	// Validate required arguments
	if fastqFile == "" {
		fmt.Fprintln(os.Stderr, "Error: -f/--fastq-file is required")
		fs.Usage()
		os.Exit(1)
	}
	
	// Determine quality encoding
	encoding := getQualityEncoding(qualType)
	
	// Open input file (with automatic gzip support)
	inputFile, err := iohelper.OpenReader(fastqFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer inputFile.Close()
	
	// Open output file (with automatic gzip support)
	outFileName := outputFile
	if outFileName == "" {
		outFileName = "-"
	}
	outFile, err := iohelper.OpenWriter(outFileName)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer outFile.Close()
	
	// Set up trim options
	opts := sickle.TrimOptions{
		QualThreshold:   qualThreshold,
		LengthThreshold: lengthThreshold,
		NoFivePrime:     noFivePrime,
		TruncateN:       truncateN,
		WindowSize:      10,
	}
	
	// Perform trimming
	stats, err := sickle.TrimSingleEnd(inputFile, outFile, encoding, opts)
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
	fs := flag.NewFlagSet("sickle pe", flag.ExitOnError)
	
	var (
		fastqFile1      string
		fastqFile2      string
		outputFile1     string
		outputFile2     string
		outputSingle    string
		qualType        string
		qualThreshold   int
		lengthThreshold int
		noFivePrime     bool
		truncateN       bool
		quiet           bool
	)
	
	cliflag.StringVar(fs, &fastqFile1, "f", "fastq-file", "", "First input FASTQ file (required)")
	cliflag.StringVar(fs, &fastqFile2, "r", "reverse-file", "", "Second input FASTQ file (required)")
	cliflag.StringVar(fs, &outputFile1, "o", "output-file", "", "First output trimmed file (required)")
	cliflag.StringVar(fs, &outputFile2, "p", "output-paired", "", "Second output trimmed file (required)")
	cliflag.StringVar(fs, &outputSingle, "s", "output-single", "", "Output single-end reads (optional)")
	cliflag.StringVar(fs, &qualType, "t", "qual-type", "sanger", "Quality type: sanger, illumina, solexa (default: sanger)")
	cliflag.IntVar(fs, &qualThreshold, "q", "qual-threshold", 20, "Threshold for trimming (default: 20)")
	cliflag.IntVar(fs, &lengthThreshold, "l", "length-threshold", 20, "Minimum length to keep (default: 20)")
	cliflag.BoolVar(fs, &noFivePrime, "x", "no-fiveprime", false, "Don't trim 5' end")
	cliflag.BoolVar(fs, &truncateN, "n", "trunc-n", false, "Truncate sequences at position of first N")
	cliflag.BoolVar(fs, &quiet, "", "quiet", false, "Don't print statistics")
	
	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: sickle pe [options]\n\n")
		fmt.Fprintf(os.Stderr, "Options:\n")
		fmt.Fprintf(os.Stderr, "  -f, --fastq-file FILE       First input FASTQ file (required)\n")
		fmt.Fprintf(os.Stderr, "  -r, --reverse-file FILE     Second input FASTQ file (required)\n")
		fmt.Fprintf(os.Stderr, "  -o, --output-file FILE      First output trimmed file (required)\n")
		fmt.Fprintf(os.Stderr, "  -p, --output-paired FILE    Second output trimmed file (required)\n")
		fmt.Fprintf(os.Stderr, "  -s, --output-single FILE    Output single-end reads (optional)\n")
		fmt.Fprintf(os.Stderr, "  -t, --qual-type TYPE        Quality type: sanger, illumina, solexa (default: sanger)\n")
		fmt.Fprintf(os.Stderr, "  -q, --qual-threshold INT    Threshold for trimming (default: 20)\n")
		fmt.Fprintf(os.Stderr, "  -l, --length-threshold INT  Minimum length to keep (default: 20)\n")
		fmt.Fprintf(os.Stderr, "  -x, --no-fiveprime          Don't trim 5' end\n")
		fmt.Fprintf(os.Stderr, "  -n, --trunc-n               Truncate sequences at position of first N\n")
		fmt.Fprintf(os.Stderr, "  --quiet                     Don't print statistics\n")
		fmt.Fprintf(os.Stderr, "\nExample:\n")
		fmt.Fprintf(os.Stderr, "  sickle pe -f input1.fastq -r input2.fastq -o output1.fastq -p output2.fastq -s singles.fastq\n")
	}
	
	fs.Parse(os.Args[2:])
	
	// Validate required arguments
	if fastqFile1 == "" || fastqFile2 == "" {
		fmt.Fprintln(os.Stderr, "Error: both -f/--fastq-file and -r/--reverse-file are required")
		fs.Usage()
		os.Exit(1)
	}
	if outputFile1 == "" || outputFile2 == "" {
		fmt.Fprintln(os.Stderr, "Error: both -o/--output-file and -p/--output-paired are required")
		fs.Usage()
		os.Exit(1)
	}
	
	// Determine quality encoding
	encoding := getQualityEncoding(qualType)
	
	// Open input files (with automatic gzip support)
	f1, err := iohelper.OpenReader(fastqFile1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening first input file: %v\n", err)
		os.Exit(1)
	}
	defer f1.Close()
	
	f2, err := iohelper.OpenReader(fastqFile2)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening second input file: %v\n", err)
		os.Exit(1)
	}
	defer f2.Close()
	
	// Open output files (with automatic gzip support)
	out1, err := iohelper.OpenWriter(outputFile1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating first output file: %v\n", err)
		os.Exit(1)
	}
	defer out1.Close()
	
	out2, err := iohelper.OpenWriter(outputFile2)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating second output file: %v\n", err)
		os.Exit(1)
	}
	defer out2.Close()
	
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
	opts := sickle.TrimOptions{
		QualThreshold:   qualThreshold,
		LengthThreshold: lengthThreshold,
		NoFivePrime:     noFivePrime,
		TruncateN:       truncateN,
		WindowSize:      10,
	}
	
	// Perform trimming
	stats, err := sickle.TrimPairedEnd(f1, f2, out1, out2, outSingle, encoding, opts)
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
	case "solexa":
		// Solexa uses a different formula, but we'll approximate with Phred64
		return fastq.Phred64
	default:
		fmt.Fprintf(os.Stderr, "Warning: unknown quality type %q, using sanger\n", qualType)
		return fastq.Phred33
	}
}

func printStats(stats *sickle.TrimStats, mode string) {
	fmt.Fprintf(os.Stderr, "\n%s Trimming Stats:\n", mode)
	fmt.Fprintf(os.Stderr, "  Total reads:     %d\n", stats.TotalReads)
	fmt.Fprintf(os.Stderr, "  Trimmed reads:   %d (%.2f%%)\n", 
		stats.TrimmedReads, 
		100.0*float64(stats.TrimmedReads)/float64(stats.TotalReads))
	fmt.Fprintf(os.Stderr, "  Discarded reads: %d (%.2f%%)\n", 
		stats.DiscardedReads, 
		100.0*float64(stats.DiscardedReads)/float64(stats.TotalReads))
	fmt.Fprintf(os.Stderr, "  Kept reads:      %d (%.2f%%)\n", 
		stats.TotalReads-stats.DiscardedReads, 
		100.0*float64(stats.TotalReads-stats.DiscardedReads)/float64(stats.TotalReads))
	fmt.Fprintf(os.Stderr, "  Total bases:     %d\n", stats.TotalBases)
	fmt.Fprintf(os.Stderr, "  Trimmed bases:   %d (%.2f%%)\n", 
		stats.TrimmedBases, 
		100.0*float64(stats.TrimmedBases)/float64(stats.TotalBases))
	fmt.Fprintln(os.Stderr)
}
