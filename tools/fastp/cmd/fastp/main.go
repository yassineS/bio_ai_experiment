package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/fastp/pkg/fastp"
)

const usage = `fastp - All-in-one FASTQ preprocessor

Usage:
  fastp -i input.fastq -o output.fastq [options]

Options:
  Input/Output:
    -i, --input FILE          Input FASTQ file (required)
    -o, --output FILE         Output FASTQ file (required)
  
  Adapter Trimming:
    -x, --adapter3 SEQ        3' adapter sequence
    -y, --adapter5 SEQ        5' adapter sequence
  
  Quality Filtering:
    -q, --qual-threshold INT  Quality threshold (default: 15)
    --qual-percent INT        Percent of bases meeting quality (default: 40)
  
  Length Filtering:
    -l, --min-length INT      Minimum read length (default: 15)
    --max-length INT          Maximum read length (0 = no limit)
  
  Content Filtering:
    --max-n-count INT         Maximum N count (default: 5)
    --max-n-percent FLOAT     Maximum N percentage (default: 20.0)
  
  Poly-tail Trimming:
    --trim-poly-g             Enable poly-G tail trimming
    --trim-poly-x             Enable poly-X tail trimming
    --poly-g-min-len INT      Minimum poly-G length (default: 10)
  
  Complexity Filtering:
    --low-complexity          Enable complexity filtering
    --complexity-threshold FLOAT  Complexity threshold (default: 0.3)
  
  Other:
    -t, --qual-type TYPE      Quality type: sanger, illumina (default: sanger)
    --quiet                   Don't print statistics

Examples:
  # Basic adapter trimming and filtering
  fastp -i input.fastq -o output.fastq -x AGATCGGAAGAGC
  
  # Comprehensive preprocessing
  fastp -i input.fastq -o output.fastq \
    -x AGATCGGAAGAGC -q 20 -l 30 \
    --trim-poly-g --max-n-count 3
  
  # With gzip files
  fastp -i input.fastq.gz -o output.fastq.gz -x AGATCGGAAGAGC

Version: 1.0.0 (Go implementation)
`

func main() {
	fs := flag.NewFlagSet("fastp", flag.ExitOnError)
	
	var (
		inputFile           string
		outputFile          string
		adapter3            string
		adapter5            string
		qualType            string
		qualThreshold       int
		qualPercent         int
		minLength           int
		maxLength           int
		maxNCount           int
		maxNPercent         float64
		trimPolyG           bool
		trimPolyX           bool
		polyGMinLen         int
		lowComplexity       bool
		complexityThreshold float64
		quiet               bool
	)
	
	// Input/Output
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input FASTQ file (required)")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output FASTQ file (required)")
	
	// Adapter trimming
	cliflag.StringVar(fs, &adapter3, "x", "adapter3", "", "3' adapter sequence")
	cliflag.StringVar(fs, &adapter5, "y", "adapter5", "", "5' adapter sequence")
	
	// Quality filtering
	cliflag.IntVar(fs, &qualThreshold, "q", "qual-threshold", 15, "Quality threshold (default: 15)")
	cliflag.IntVar(fs, &qualPercent, "", "qual-percent", 40, "Percent of bases meeting quality (default: 40)")
	
	// Length filtering
	cliflag.IntVar(fs, &minLength, "l", "min-length", 15, "Minimum read length (default: 15)")
	cliflag.IntVar(fs, &maxLength, "", "max-length", 0, "Maximum read length (0 = no limit)")
	
	// Content filtering
	cliflag.IntVar(fs, &maxNCount, "", "max-n-count", 5, "Maximum N count (default: 5)")
	cliflag.Float64Var(fs, &maxNPercent, "", "max-n-percent", 20.0, "Maximum N percentage (default: 20.0)")
	
	// Poly-tail trimming
	cliflag.BoolVar(fs, &trimPolyG, "", "trim-poly-g", false, "Enable poly-G tail trimming")
	cliflag.BoolVar(fs, &trimPolyX, "", "trim-poly-x", false, "Enable poly-X tail trimming")
	cliflag.IntVar(fs, &polyGMinLen, "", "poly-g-min-len", 10, "Minimum poly-G length (default: 10)")
	
	// Complexity filtering
	cliflag.BoolVar(fs, &lowComplexity, "", "low-complexity", false, "Enable complexity filtering")
	cliflag.Float64Var(fs, &complexityThreshold, "", "complexity-threshold", 0.3, "Complexity threshold (default: 0.3)")
	
	// Other
	cliflag.StringVar(fs, &qualType, "t", "qual-type", "sanger", "Quality type: sanger, illumina (default: sanger)")
	cliflag.BoolVar(fs, &quiet, "", "quiet", false, "Don't print statistics")
	
	fs.Usage = func() {
		fmt.Print(usage)
	}
	
	if len(os.Args) < 2 {
		fs.Usage()
		os.Exit(1)
	}
	
	fs.Parse(os.Args[1:])
	
	// Validate required arguments
	if inputFile == "" || outputFile == "" {
		fmt.Fprintln(os.Stderr, "Error: both -i/--input and -o/--output are required")
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
	output, err := iohelper.OpenWriter(outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer output.Close()
	
	// Set up processing options
	opts := fastp.ProcessOptions{
		Adapter3:            adapter3,
		Adapter5:            adapter5,
		QualThreshold:       qualThreshold,
		MinLength:           minLength,
		MaxLength:           maxLength,
		QualPercent:         qualPercent,
		LowComplexity:       lowComplexity,
		ComplexityThreshold: complexityThreshold,
		TrimPolyG:           trimPolyG,
		TrimPolyX:           trimPolyX,
		PolyGMinLen:         polyGMinLen,
		MaxNCount:           maxNCount,
		MaxNPercent:         maxNPercent,
		LengthRequired:      minLength,
		LengthLimit:         maxLength,
	}
	
	// Perform processing
	stats, err := fastp.ProcessSingleEnd(input, output, encoding, opts)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error during processing: %v\n", err)
		os.Exit(1)
	}
	
	// Print statistics
	if !quiet {
		printStats(stats)
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

func printStats(stats *fastp.ProcessStats) {
	fmt.Fprintf(os.Stderr, "\nFastp Processing Statistics:\n")
	fmt.Fprintf(os.Stderr, "  Total reads:           %d\n", stats.TotalReads)
	fmt.Fprintf(os.Stderr, "  Total bases:           %d\n", stats.TotalBases)
	fmt.Fprintf(os.Stderr, "  Clean reads:           %d (%.2f%%)\n",
		stats.CleanReads,
		100.0*float64(stats.CleanReads)/float64(stats.TotalReads))
	fmt.Fprintf(os.Stderr, "  Clean bases:           %d (%.2f%%)\n",
		stats.CleanBases,
		100.0*float64(stats.CleanBases)/float64(stats.TotalBases))
	
	if stats.AdapterTrimmedReads > 0 {
		fmt.Fprintf(os.Stderr, "  Adapter trimmed:       %d (%.2f%%)\n",
			stats.AdapterTrimmedReads,
			100.0*float64(stats.AdapterTrimmedReads)/float64(stats.TotalReads))
		fmt.Fprintf(os.Stderr, "  Adapter bases removed: %d\n", stats.AdapterTrimmedBases)
	}
	
	if stats.PolyGTrimmedReads > 0 {
		fmt.Fprintf(os.Stderr, "  Poly-G trimmed:        %d (%.2f%%)\n",
			stats.PolyGTrimmedReads,
			100.0*float64(stats.PolyGTrimmedReads)/float64(stats.TotalReads))
		fmt.Fprintf(os.Stderr, "  Poly-G bases removed:  %d\n", stats.PolyGTrimmedBases)
	}
	
	if stats.LowQualityReads > 0 {
		fmt.Fprintf(os.Stderr, "  Low quality filtered:  %d (%.2f%%)\n",
			stats.LowQualityReads,
			100.0*float64(stats.LowQualityReads)/float64(stats.TotalReads))
	}
	
	if stats.TooShortReads > 0 {
		fmt.Fprintf(os.Stderr, "  Too short filtered:    %d (%.2f%%)\n",
			stats.TooShortReads,
			100.0*float64(stats.TooShortReads)/float64(stats.TotalReads))
	}
	
	if stats.TooLongReads > 0 {
		fmt.Fprintf(os.Stderr, "  Too long filtered:     %d (%.2f%%)\n",
			stats.TooLongReads,
			100.0*float64(stats.TooLongReads)/float64(stats.TotalReads))
	}
	
	if stats.TooManyNReads > 0 {
		fmt.Fprintf(os.Stderr, "  Too many N filtered:   %d (%.2f%%)\n",
			stats.TooManyNReads,
			100.0*float64(stats.TooManyNReads)/float64(stats.TotalReads))
	}
	
	fmt.Fprintln(os.Stderr)
}
