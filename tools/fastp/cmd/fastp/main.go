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
  fastp -I read1.fastq -O out1.fastq --in2 read2.fastq --out2 out2.fastq [options]

Options:
  Input/Output:
    -i, --input FILE          Input FASTQ file (single-end, required)
    -o, --output FILE         Output FASTQ file (single-end, required)
    -I, --in1 FILE            Input FASTQ file read 1 (paired-end)
    --in2 FILE                Input FASTQ file read 2 (paired-end)
    -O, --out1 FILE           Output FASTQ file read 1 (paired-end)
    --out2 FILE               Output FASTQ file read 2 (paired-end)
  
  Adapter Trimming:
    -x, --adapter3 SEQ        3' adapter sequence
    -y, --adapter5 SEQ        5' adapter sequence
    --detect-adapter          Auto-detect adapter sequences
  
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

  Sliding-window Quality Trimming:
    -5, --cut-front           Drop low-quality bases from the 5' end (sliding window)
    -3, --cut-tail            Drop low-quality bases from the 3' end (sliding window)
    -r, --cut-right           Cut from the first low-quality window onward (5'->3')
    -W, --cut-window-size INT Window size for --cut-front/--cut-tail/--cut-right (default: 4)
    -M, --cut-mean-quality INT  Mean quality threshold for the sliding window (default: 20)

  Complexity Filtering:
    --low-complexity          Enable complexity filtering
    --complexity-threshold FLOAT  Complexity threshold (default: 0.3)
  
  UMI/Barcode Processing:
    --umi                     Enable UMI processing (use with --umi_loc / --umi_len)
    --umi_loc STRING          UMI location: read1, read2, per_read, index1,
                              index2, per_index (default: read1 SE / per_read PE)
    --umi_len INT             UMI length in bases (read1/read2/per_read)
    --umi_prefix STRING       Prefix prepended to UMI in the read name
    --umi_skip INT            Bases to skip after the UMI bases (default: 0)
    --umi-length INT          [legacy] alias of --umi_len
    --umi-location STRING     [legacy] alias of --umi_loc
    --umi-skip INT            [legacy] alias of --umi_skip

  Duplication Evaluation:
    --dup_calc_accuracy INT   Duplication accuracy bucket (1-6; 0 disables)
    --dedup                   Drop duplicate reads from the output stream
  
  Base Correction:
    --base-correction         Enable base correction
    --correction-threshold INT  Base correction quality threshold (default: 20)
  
  Overlap Analysis (Paired-end):
    --merge-overlap           Merge overlapping paired-end reads
    --min-overlap INT         Minimum overlap length (default: 30)
    --max-mismatch INT        Maximum mismatches in overlap (default: 5)
  
  Performance:
    -w, --threads INT         Number of threads (default: 1)
  
  Reporting:
    --html FILE               Self-contained HTML report (no JS, no CDN)
    --json FILE               JSON report (fastp-compatible schema)
    --detect_adapter_for_pe   Overlap-based adapter detection for PE reads

  Other:
    -t, --qual-type TYPE      Quality type: sanger, illumina (default: sanger)
    -h, --help                Show this help and exit
    --quiet                   Don't print statistics

Examples:
  # Basic adapter trimming and filtering
  fastp -i input.fastq -o output.fastq -x AGATCGGAAGAGC
  
  # Auto-detect adapter
  fastp -i input.fastq -o output.fastq --detect-adapter
  
  # With UMI extraction (legacy flags still work)
  fastp -i input.fastq -o output.fastq --umi --umi_loc read1 --umi_len 8

  # Duplication evaluation + drop duplicates
  fastp -i input.fastq -o output.fastq --dup_calc_accuracy 3 --dedup
  
  # Sliding-window quality trimming from both ends
  fastp -i input.fastq -o output.fastq --cut-front --cut-tail -W 4 -M 20

  # Base correction
  fastp -i input.fastq -o output.fastq --base-correction
  
  # Merge overlapping paired-end reads
  fastp -I R1.fastq -O out1.fastq --in2 R2.fastq --out2 out2.fastq --merge-overlap
  
  # Multi-threaded with HTML + JSON reports
  fastp -i input.fastq -o output.fastq -w 4 --html report.html --json report.json

  # PE adapter detection via overlap analysis
  fastp -I r1.fq -O r1.out.fq --in2 r2.fq --out2 r2.out.fq --detect_adapter_for_pe

  # Comprehensive preprocessing
  fastp -i input.fastq -o output.fastq \
    -x AGATCGGAAGAGC -q 20 -l 30 \
    --trim-poly-g --max-n-count 3 \
    --base-correction -w 4 --html report.html --json report.json

Version: 1.0.0 (Go implementation)
`

func main() {
	fs := flag.NewFlagSet("fastp", flag.ExitOnError)

	var (
		inputFile           string
		outputFile          string
		in1File             string
		in2File             string
		out1File            string
		out2File            string
		adapter3            string
		adapter5            string
		detectAdapter       bool
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
		cutFront            bool
		cutTail             bool
		cutRight            bool
		cutWindowSize       int
		cutMeanQuality      int
		lowComplexity       bool
		complexityThreshold float64
		quiet               bool
		// New features
		umiLength           int
		umiLocation         string
		umiSkip             int
		baseCorrection      bool
		correctionThreshold int
		mergeOverlap        bool
		minOverlap          int
		maxMismatch         int
		threads             int
		htmlReport          string
		jsonReport          string
		detectAdapterForPE  bool
		showHelp            bool
		// Duplication evaluation
		dupCalcAccuracy int
		dedup           bool
		// UMI (fastp-aligned flag names)
		umiEnable bool
		umiLoc    string
		umiLen    int
		umiPrefix string
	)

	// Input/Output
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input FASTQ file (single-end)")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output FASTQ file (single-end)")
	cliflag.StringVar(fs, &in1File, "I", "in1", "", "Input FASTQ file read 1 (paired-end)")
	cliflag.StringVar(fs, &in2File, "", "in2", "", "Input FASTQ file read 2 (paired-end)")
	cliflag.StringVar(fs, &out1File, "O", "out1", "", "Output FASTQ file read 1 (paired-end)")
	cliflag.StringVar(fs, &out2File, "", "out2", "", "Output FASTQ file read 2 (paired-end)")

	// Adapter trimming
	cliflag.StringVar(fs, &adapter3, "x", "adapter3", "", "3' adapter sequence")
	cliflag.StringVar(fs, &adapter5, "y", "adapter5", "", "5' adapter sequence")
	cliflag.BoolVar(fs, &detectAdapter, "", "detect-adapter", false, "Auto-detect adapter sequences")

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

	// Sliding-window quality trimming
	cliflag.BoolVar(fs, &cutFront, "5", "cut-front", false, "Drop low-quality bases from the 5' end (sliding window)")
	cliflag.BoolVar(fs, &cutTail, "3", "cut-tail", false, "Drop low-quality bases from the 3' end (sliding window)")
	cliflag.BoolVar(fs, &cutRight, "r", "cut-right", false, "Cut from the first low-quality window onward (5'->3')")
	cliflag.IntVar(fs, &cutWindowSize, "W", "cut-window-size", 4, "Window size for sliding-window trimming (default: 4)")
	cliflag.IntVar(fs, &cutMeanQuality, "M", "cut-mean-quality", 20, "Mean quality threshold for the sliding window (default: 20)")

	// Complexity filtering
	cliflag.BoolVar(fs, &lowComplexity, "", "low-complexity", false, "Enable complexity filtering")
	cliflag.Float64Var(fs, &complexityThreshold, "", "complexity-threshold", 0.3, "Complexity threshold (default: 0.3)")

	// UMI/barcode processing — legacy flag names.
	cliflag.IntVar(fs, &umiLength, "", "umi-length", 0, "UMI length (legacy alias of --umi_len)")
	cliflag.StringVar(fs, &umiLocation, "", "umi-location", "read1", "UMI location (legacy alias of --umi_loc)")
	cliflag.IntVar(fs, &umiSkip, "", "umi-skip", 0, "Bases to skip after UMI (legacy alias of --umi_skip)")
	// UMI/barcode processing — fastp-aligned flag names.
	cliflag.BoolVar(fs, &umiEnable, "", "umi", false, "Enable UMI processing")
	cliflag.StringVar(fs, &umiLoc, "", "umi_loc", "", "UMI location: read1|read2|per_read|index1|index2|per_index")
	cliflag.IntVar(fs, &umiLen, "", "umi_len", 0, "UMI length in bases (read1/read2/per_read modes)")
	cliflag.StringVar(fs, &umiPrefix, "", "umi_prefix", "", "Prefix prepended to UMI in the read name")
	cliflag.IntVar(fs, &umiSkip, "", "umi_skip", 0, "Bases to skip after UMI bases")

	// Duplication evaluation
	cliflag.IntVar(fs, &dupCalcAccuracy, "", "dup_calc_accuracy", 0, "Duplication accuracy bucket (1-6; 0 = disabled)")
	cliflag.BoolVar(fs, &dedup, "", "dedup", false, "Drop duplicate reads from the output stream")

	// Base correction
	cliflag.BoolVar(fs, &baseCorrection, "", "base-correction", false, "Enable base correction")
	cliflag.IntVar(fs, &correctionThreshold, "", "correction-threshold", 20, "Base correction quality threshold")

	// Overlap analysis (paired-end)
	cliflag.BoolVar(fs, &mergeOverlap, "", "merge-overlap", false, "Merge overlapping paired-end reads")
	cliflag.IntVar(fs, &minOverlap, "", "min-overlap", 30, "Minimum overlap length")
	cliflag.IntVar(fs, &maxMismatch, "", "max-mismatch", 5, "Maximum mismatches in overlap")

	// Multi-threading
	cliflag.IntVar(fs, &threads, "w", "threads", 1, "Number of threads (default: 1)")

	// Reporting outputs. --html is long-only (upstream fastp uses -h for
	// help; we keep -h reserved for help to avoid colliding with that
	// muscle memory). --json is also long-only for symmetry.
	cliflag.StringVar(fs, &htmlReport, "", "html", "", "HTML report output file")
	cliflag.StringVar(fs, &jsonReport, "", "json", "", "JSON report output file")
	cliflag.BoolVar(fs, &showHelp, "h", "help", false, "Show usage and exit")

	// Automatic adapter detection
	cliflag.BoolVar(fs, &detectAdapterForPE, "", "detect_adapter_for_pe", false, "Enable overlap-based adapter detection for paired-end")

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

	if showHelp {
		fs.Usage()
		os.Exit(0)
	}

	// Determine mode: paired-end or single-end
	isPaired := (in1File != "" && in2File != "" && out1File != "" && out2File != "")
	isSingle := (inputFile != "" && outputFile != "")

	if !isPaired && !isSingle {
		fmt.Fprintln(os.Stderr, "Error: must specify either:")
		fmt.Fprintln(os.Stderr, "  Single-end: -i/--input and -o/--output")
		fmt.Fprintln(os.Stderr, "  Paired-end: -I/--in1, --in2, -O/--out1, --out2")
		fs.Usage()
		os.Exit(1)
	}

	if isPaired && isSingle {
		fmt.Fprintln(os.Stderr, "Error: cannot specify both single-end and paired-end options")
		os.Exit(1)
	}

	// Determine quality encoding
	encoding := getQualityEncoding(qualType)

	// Set up processing options
	opts := fastp.ProcessOptions{
		Adapter3:            adapter3,
		Adapter5:            adapter5,
		DetectAdapter:       detectAdapter,
		QualThreshold:       qualThreshold,
		MinLength:           minLength,
		MaxLength:           maxLength,
		QualPercent:         qualPercent,
		LowComplexity:       lowComplexity,
		ComplexityThreshold: complexityThreshold,
		TrimPolyG:           trimPolyG,
		TrimPolyX:           trimPolyX,
		PolyGMinLen:         polyGMinLen,
		CutFront:            cutFront,
		CutTail:             cutTail,
		CutRight:            cutRight,
		CutWindowSize:       cutWindowSize,
		CutMeanQuality:      cutMeanQuality,
		MaxNCount:           maxNCount,
		MaxNPercent:         maxNPercent,
		LengthRequired:      minLength,
		LengthLimit:         maxLength,
		UMILength:           umiLength,
		UMILocation:         umiLocation,
		UMI:                 umiEnable,
		UMILoc:              umiLoc,
		UMILen:              umiLen,
		UMIPrefix:           umiPrefix,
		UMISkip:             umiSkip,
		DupCalcAccuracy:     dupCalcAccuracy,
		Dedup:               dedup,
		BaseCorrection:      baseCorrection,
		CorrectionThreshold: correctionThreshold,
		MergeOverlap:        mergeOverlap,
		MinOverlap:          minOverlap,
		MaxMismatch:         maxMismatch,
		Threads:             threads,
		HTMLReport:          htmlReport,
		JSONReport:          jsonReport,
		DetectAdapterPE:     detectAdapterForPE,
		DetectAdapterSE:     detectAdapter, // legacy --detect-adapter triggers SE detection
	}

	var stats *fastp.ProcessStats
	var err error

	if isPaired {
		// Paired-end mode
		input1, err := iohelper.OpenReader(in1File)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening input file 1: %v\n", err)
			os.Exit(1)
		}
		defer input1.Close()

		input2, err := iohelper.OpenReader(in2File)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening input file 2: %v\n", err)
			os.Exit(1)
		}
		defer input2.Close()

		output1, err := iohelper.OpenWriter(out1File)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file 1: %v\n", err)
			os.Exit(1)
		}
		defer output1.Close()

		output2, err := iohelper.OpenWriter(out2File)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file 2: %v\n", err)
			os.Exit(1)
		}
		defer output2.Close()

		stats, err = fastp.ProcessPairedEnd(input1, input2, output1, output2, encoding, opts)
	} else {
		// Single-end mode
		input, err := iohelper.OpenReader(inputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
			os.Exit(1)
		}
		defer input.Close()

		output, err := iohelper.OpenWriter(outputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
		defer output.Close()

		stats, err = fastp.ProcessSingleEnd(input, output, encoding, opts)
	}

	if err != nil {
		fmt.Fprintf(os.Stderr, "Error during processing: %v\n", err)
		os.Exit(1)
	}

	// Generate HTML report if requested
	if htmlReport != "" {
		if err := fastp.WriteHTMLReport(htmlReport, stats); err != nil {
			fmt.Fprintf(os.Stderr, "Error generating HTML report: %v\n", err)
			os.Exit(1)
		}
		if !quiet {
			fmt.Fprintf(os.Stderr, "HTML report written to: %s\n", htmlReport)
		}
	}

	// Generate JSON report if requested
	if jsonReport != "" {
		if err := fastp.WriteJSONReport(jsonReport, stats); err != nil {
			fmt.Fprintf(os.Stderr, "Error generating JSON report: %v\n", err)
			os.Exit(1)
		}
		if !quiet {
			fmt.Fprintf(os.Stderr, "JSON report written to: %s\n", jsonReport)
		}
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

	if stats.DetectedAdapter != "" {
		fmt.Fprintf(os.Stderr, "  Detected adapter:      %s\n", stats.DetectedAdapter)
	}

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

	if stats.QualityCutReads > 0 {
		fmt.Fprintf(os.Stderr, "  Sliding-window cut:    %d (%.2f%%)\n",
			stats.QualityCutReads,
			100.0*float64(stats.QualityCutReads)/float64(stats.TotalReads))
		fmt.Fprintf(os.Stderr, "  Sliding-window bases:  %d\n", stats.QualityCutBases)
	}

	if stats.UMIExtracted > 0 || stats.UMIProcessed > 0 {
		processed := stats.UMIProcessed
		if processed == 0 {
			processed = stats.UMIExtracted
		}
		fmt.Fprintf(os.Stderr, "  UMIs processed:        %d\n", processed)
	}

	if stats.DupTotal > 0 {
		fmt.Fprintf(os.Stderr, "  Duplication rate:      %.2f%% (%d / %d)\n",
			100.0*stats.DupRate, int64(stats.DupRate*float64(stats.DupTotal)+0.5), stats.DupTotal)
		if stats.DedupDropped > 0 {
			fmt.Fprintf(os.Stderr, "  Dedup dropped:         %d\n", stats.DedupDropped)
		}
	}

	if stats.BasesCorrected > 0 {
		fmt.Fprintf(os.Stderr, "  Bases corrected:       %d\n", stats.BasesCorrected)
	}

	if stats.MergedReads > 0 {
		fmt.Fprintf(os.Stderr, "  Overlapping merged:    %d (%.2f%%)\n",
			stats.MergedReads,
			100.0*float64(stats.MergedReads)/float64(stats.TotalReads/2))
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
