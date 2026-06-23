package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/fastp/pkg/fastp"
)

// version is the fastp port's version, reported by -v/--version.
const version = "1.0.0"

const usage = `fastp - All-in-one FASTQ preprocessor

Usage:
  fastp -i input.fastq -o output.fastq [options]
  fastp -i read1.fastq -I read2.fastq -o out1.fastq -O out2.fastq [options]

Options:
  Input/Output (upstream-exact short flags):
    -i, --in1 FILE            Read 1 input (or single-end input; --input alias)
    -I, --in2 FILE            Read 2 input (paired-end)
    -o, --out1 FILE           Read 1 output (or single-end output; --output alias)
    -O, --out2 FILE           Read 2 output (paired-end)

  Adapter Trimming:
    -x, --adapter3 SEQ        3' adapter sequence
    -y, --adapter5 SEQ        5' adapter sequence
    --detect-adapter          Auto-detect adapter sequences
    -A, --disable_adapter_trimming  Disable all adapter trimming
    --adapter_fasta FILE      Trim reads by all sequences in this FASTA file

  Quality Filtering:
    -q, --qualified_quality_phred INT  Qualified base quality threshold (default: 15)
    -u, --unqualified_percent_limit INT  Max percent of bases below threshold (default: 40)
    -Q, --disable_quality_filtering  Disable quality filtering (N-base, low-quality-percent)

  Length Filtering:
    -l, --length_required INT Minimum read length (default: 15; --min-length alias)
    --max-length INT          Maximum read length (0 = no limit)
    -L, --disable_length_filtering   Disable length filtering (length_required / length_limit)

  Content Filtering:
    -n, --n_base_limit INT    Maximum N count (default: 5; --max-n-count alias)
    --max-n-percent FLOAT     Maximum N percentage (default: 20.0)

  Poly-tail Trimming:
    -g, --trim_poly_g         Enable poly-G tail trimming
    -G, --disable_trim_poly_g Disable poly-G tail trimming (overrides --trim_poly_g)
    --trim-poly-x             Enable poly-X tail trimming
    --poly-g-min-len INT      Minimum poly-G length (default: 10)
    --poly_x_min_len INT      Minimum poly-X length (default: 10)

  Overlap Merge (Paired-end):
    -m, --merge               Merge overlapping read pairs into single reads
    --merged_out FILE         File to store merged reads (or use stdout)
    --include_unmerged        Write unmerged/unpaired reads to the merge stream

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
    -c, --correction          Enable overlap-based base correction (PE only)
    --overlap_len_require INT     Min length to detect PE overlap (default: 30)
    --overlap_diff_limit INT      Max mismatched bases in PE overlap (default: 5)
    --overlap_diff_percent_limit INT  Max mismatch percent in PE overlap (default: 20)

  Overrepresentation Analysis:
    -p, --overrepresentation_analysis  Enable overrepresented sequence analysis
    -P, --overrepresentation_sampling INT  1-in-N sampling rate (default: 20)

  Output Splitting:
    -s, --split INT           Split output into INT files (2-999)
    -S, --split_by_lines INT  Split output by lines per file (>=1000, mult. of 4)
    -d, --split_prefix_digits INT  Zero-pad width for split prefix (default: 4)

  Performance:
    -w, --threads INT         Number of threads (default: 1)
  
  Reporting:
    -h, --html FILE           Self-contained HTML report (no JS, no CDN)
    -j, --json FILE           JSON report (fastp-compatible schema)
    -2, --detect_adapter_for_pe   Overlap-based adapter detection for PE reads

  Other:
    -t, --qual-type TYPE      Quality type: sanger, illumina (default: sanger)
    -?, --help                Show this help and exit (-h is the HTML report)
    -v, --version             Show version information and exit
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
  fastp -i R1.fastq -I R2.fastq -o out1.fastq -O out2.fastq --merge-overlap

  # Multi-threaded with HTML + JSON reports
  fastp -i input.fastq -o output.fastq -w 4 -h report.html -j report.json

  # PE adapter detection via overlap analysis
  fastp -i r1.fq -I r2.fq -o r1.out.fq -O r2.out.fq --detect_adapter_for_pe

  # Comprehensive preprocessing
  fastp -i input.fastq -o output.fastq \
    -x AGATCGGAAGAGC -q 20 -l 30 \
    --trim-poly-g -n 3 \
    --base-correction -w 4 -h report.html -j report.json

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
		disableAdapter      bool
		disableQualityFilt  bool
		disableLengthFilt   bool
		disablePolyG        bool
		adapterFasta        string
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
		polyXMinLen         int
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
		merge               bool
		mergedOut           string
		includeUnmerged     bool
		threads             int
		htmlReport          string
		jsonReport          string
		detectAdapterForPE  bool
		showHelp            bool
		showVersion         bool
		// Duplication evaluation
		dupCalcAccuracy int
		dedup           bool
		// UMI (fastp-aligned flag names)
		umiEnable bool
		umiLoc    string
		umiLen    int
		umiPrefix string
		// Overlap-based correction (fastp-aligned flag names)
		correction              bool
		overlapLenRequire       int
		overlapDiffLimit        int
		overlapDiffPercentLimit int
		// Overrepresentation analysis
		overrepAnalysis bool
		overrepSampling int
		// Output splitting
		splitNumber       int
		splitByLines      int
		splitPrefixDigits int
	)

	// Input/Output.
	//
	// Short flags are upstream-exact so a stock fastp command line is
	// drop-in: -i=read1 (or the single-end input), -I=read2, -o=read1
	// output (or the single-end output), -O=read2 output. read1/read1-out
	// are stored in in1File/out1File; single-end mode is detected later by
	// the absence of in2File. The legacy long names (--input/--output) are
	// retained as aliases of --in1/--out1 so existing single-end command
	// lines keep working, and --in2/--out2 remain as the GNU long forms.
	cliflag.StringVar(fs, &in1File, "i", "in1", "", "Input FASTQ file read 1 (or single-end input)")
	cliflag.StringVar(fs, &in1File, "", "input", "", "Single-end input FASTQ (alias of --in1)")
	cliflag.StringVar(fs, &in2File, "I", "in2", "", "Input FASTQ file read 2 (paired-end)")
	cliflag.StringVar(fs, &out1File, "o", "out1", "", "Output FASTQ file read 1 (or single-end output)")
	cliflag.StringVar(fs, &out1File, "", "output", "", "Single-end output FASTQ (alias of --out1)")
	cliflag.StringVar(fs, &out2File, "O", "out2", "", "Output FASTQ file read 2 (paired-end)")

	// Adapter trimming
	cliflag.StringVar(fs, &adapter3, "x", "adapter3", "", "3' adapter sequence")
	cliflag.StringVar(fs, &adapter5, "y", "adapter5", "", "5' adapter sequence")
	cliflag.BoolVar(fs, &detectAdapter, "", "detect-adapter", false, "Auto-detect adapter sequences")
	cliflag.BoolVar(fs, &disableAdapter, "A", "disable_adapter_trimming", false, "Disable all adapter trimming")
	cliflag.StringVar(fs, &adapterFasta, "", "adapter_fasta", "", "Trim reads by all sequences in this FASTA file")

	// Quality filtering. Short flags and upstream long aliases are
	// upstream-exact: -q/--qualified_quality_phred sets the per-base
	// qualified threshold, and -u/--unqualified_percent_limit sets the
	// maximum percentage of sub-threshold bases tolerated before a read is
	// discarded (NOT the percent that must pass). --qual-threshold and
	// --qual-percent are retained as legacy aliases.
	cliflag.IntVar(fs, &qualThreshold, "q", "qual-threshold", 15, "Qualified base quality threshold (default: 15)")
	cliflag.IntVar(fs, &qualThreshold, "", "qualified_quality_phred", 15, "Qualified base quality threshold (upstream alias of -q)")
	cliflag.IntVar(fs, &qualPercent, "u", "qual-percent", 40, "Max percent of bases allowed below the quality threshold (default: 40)")
	cliflag.IntVar(fs, &qualPercent, "", "unqualified_percent_limit", 40, "Max percent of unqualified bases (upstream alias of -u)")
	cliflag.BoolVar(fs, &disableQualityFilt, "Q", "disable_quality_filtering", false, "Disable quality filtering (N-base, low-quality-percent, avg-quality)")

	// Length filtering. -l/--length_required is upstream-exact;
	// --min-length is a legacy alias.
	cliflag.IntVar(fs, &minLength, "l", "min-length", 15, "Minimum read length (default: 15)")
	cliflag.IntVar(fs, &minLength, "", "length_required", 15, "Minimum read length (upstream alias of -l)")
	cliflag.IntVar(fs, &maxLength, "", "max-length", 0, "Maximum read length (0 = no limit)")
	cliflag.BoolVar(fs, &disableLengthFilt, "L", "disable_length_filtering", false, "Disable length filtering (length_required / length_limit discard)")

	// Content filtering. -n/--n_base_limit is upstream-exact; --max-n-count
	// is a legacy alias. (--max-n-percent has no upstream equivalent.)
	cliflag.IntVar(fs, &maxNCount, "n", "max-n-count", 5, "Maximum N count (default: 5)")
	cliflag.IntVar(fs, &maxNCount, "", "n_base_limit", 5, "Maximum N count (upstream alias of -n)")
	cliflag.Float64Var(fs, &maxNPercent, "", "max-n-percent", 20.0, "Maximum N percentage (default: 20.0)")

	// Poly-tail trimming
	cliflag.BoolVar(fs, &trimPolyG, "", "trim-poly-g", false, "Enable poly-G tail trimming")
	cliflag.BoolVar(fs, &trimPolyG, "g", "trim_poly_g", false, "Enable poly-G tail trimming (upstream spelling)")
	cliflag.BoolVar(fs, &disablePolyG, "G", "disable_trim_poly_g", false, "Disable poly-G tail trimming (overrides --trim_poly_g)")
	cliflag.BoolVar(fs, &trimPolyX, "", "trim-poly-x", false, "Enable poly-X tail trimming")
	cliflag.IntVar(fs, &polyGMinLen, "", "poly-g-min-len", 10, "Minimum poly-G length (default: 10)")
	cliflag.IntVar(fs, &polyGMinLen, "", "poly_g_min_len", 10, "Minimum poly-G length (alias of --poly-g-min-len)")
	cliflag.IntVar(fs, &polyXMinLen, "", "poly_x_min_len", 10, "Minimum poly-X length (default: 10)")

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
	cliflag.BoolVar(fs, &mergeOverlap, "", "merge-overlap", false, "Merge overlapping paired-end reads (legacy)")
	cliflag.IntVar(fs, &minOverlap, "", "min-overlap", 30, "Minimum overlap length")
	cliflag.IntVar(fs, &maxMismatch, "", "max-mismatch", 5, "Maximum mismatches in overlap")

	// Overlap-driven merge writer (upstream -m/--merge family).
	cliflag.BoolVar(fs, &merge, "m", "merge", false, "Merge overlapping read pairs into single reads")
	cliflag.StringVar(fs, &mergedOut, "", "merged_out", "", "File to store merged reads")
	cliflag.BoolVar(fs, &includeUnmerged, "", "include_unmerged", false, "Write unmerged/unpaired reads to the merge stream")

	// Overlap-based base correction (paired-end), upstream flag names.
	cliflag.BoolVar(fs, &correction, "c", "correction", false, "Enable base correction in overlapped regions (PE only)")
	cliflag.IntVar(fs, &overlapLenRequire, "", "overlap_len_require", 30, "Minimum length to detect PE overlap (default 30)")
	cliflag.IntVar(fs, &overlapDiffLimit, "", "overlap_diff_limit", 5, "Maximum mismatched bases in PE overlap (default 5)")
	cliflag.IntVar(fs, &overlapDiffPercentLimit, "", "overlap_diff_percent_limit", 20, "Maximum mismatch percentage in PE overlap (default 20)")

	// Overrepresentation analysis.
	cliflag.BoolVar(fs, &overrepAnalysis, "p", "overrepresentation_analysis", false, "Enable overrepresented sequence analysis")
	cliflag.IntVar(fs, &overrepSampling, "P", "overrepresentation_sampling", 20, "One in N reads sampled for ORA (1-10000, default 20)")

	// Output splitting.
	cliflag.IntVar(fs, &splitNumber, "s", "split", 0, "Split output into this many files (2-999)")
	cliflag.IntVar(fs, &splitByLines, "S", "split_by_lines", 0, "Split output by max lines per file (>=1000, multiple of 4)")
	cliflag.IntVar(fs, &splitPrefixDigits, "d", "split_prefix_digits", 4, "Zero-pad width for split file prefixes (1-10, default 4)")

	// Multi-threading
	cliflag.IntVar(fs, &threads, "w", "threads", 1, "Number of threads (default: 1)")

	// Reporting outputs. Short flags are upstream-exact: -h=html and -j=json
	// (upstream reserves -? / --help for help, NOT -h). We keep --help and
	// add -? so the help path is still reachable without colliding with the
	// html report flag.
	cliflag.StringVar(fs, &htmlReport, "h", "html", "", "HTML report output file")
	cliflag.StringVar(fs, &jsonReport, "j", "json", "", "JSON report output file")
	cliflag.BoolVar(fs, &showHelp, "?", "help", false, "Show usage and exit")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version information and exit")

	// Automatic adapter detection. -2 is upstream's short flag for
	// detect_adapter_for_pe.
	cliflag.BoolVar(fs, &detectAdapterForPE, "2", "detect_adapter_for_pe", false, "Enable overlap-based adapter detection for paired-end")

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

	if err := cliflag.Parse(fs, os.Args[1:]); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(2)
	}

	if showHelp {
		fs.Usage()
		os.Exit(0)
	}

	if showVersion {
		fmt.Printf("fastp version %s\n", version)
		os.Exit(0)
	}

	// Determine mode: paired-end or single-end. As in upstream fastp,
	// paired-end is selected purely by the presence of both -i/--in1 and
	// -I/--in2; the single-end input/output reuse in1File/out1File (which
	// -i/-o feed). In merge mode the merged output is written to
	// --merged_out, so out1/out2 are optional (they only carry unmerged
	// pairs when --include_unmerged is NOT used).
	inputFile = in1File
	outputFile = out1File
	isPaired := (in1File != "" && in2File != "" && out1File != "" && out2File != "")
	if merge && in1File != "" && in2File != "" && (mergedOut != "" || (out1File != "" && out2File != "")) {
		isPaired = true
	}
	isSingle := (in1File != "" && in2File == "" && out1File != "")

	if !isPaired && !isSingle {
		fmt.Fprintln(os.Stderr, "Error: must specify either:")
		fmt.Fprintln(os.Stderr, "  Single-end: -i/--in1 (or --input) and -o/--out1 (or --output)")
		fmt.Fprintln(os.Stderr, "  Paired-end: -i/--in1, -I/--in2, -o/--out1, -O/--out2")
		fs.Usage()
		os.Exit(1)
	}

	if isPaired && isSingle {
		fmt.Fprintln(os.Stderr, "Error: cannot specify both single-end and paired-end options")
		os.Exit(1)
	}

	// Determine quality encoding
	encoding := getQualityEncoding(qualType)

	// Load --adapter_fasta sequences if requested. Sequences shorter than
	// 6 bp are skipped with a stderr warning, mirroring upstream.
	var adapterFastaSeqs []string
	if adapterFasta != "" {
		fa, oerr := iohelper.OpenReader(adapterFasta)
		if oerr != nil {
			fmt.Fprintf(os.Stderr, "Error opening --adapter_fasta file: %v\n", oerr)
			os.Exit(1)
		}
		seqs, skipped := fastp.LoadAdapterFasta(fa)
		fa.Close()
		for _, s := range skipped {
			fmt.Fprintf(os.Stderr, "skip too short adapter sequence in %s (6bp required): %s\n", adapterFasta, s)
		}
		adapterFastaSeqs = seqs
	}

	// Set up processing options
	opts := fastp.ProcessOptions{
		Adapter3:                adapter3,
		Adapter5:                adapter5,
		DetectAdapter:           detectAdapter,
		DisableAdapterTrimming:  disableAdapter,
		DisableQualityFiltering: disableQualityFilt,
		DisableLengthFiltering:  disableLengthFilt,
		DisableTrimPolyG:        disablePolyG,
		AdapterFasta:            adapterFastaSeqs,
		PolyXMinLen:             polyXMinLen,
		Merge:                   merge,
		IncludeUnmerged:         includeUnmerged,
		QualThreshold:           qualThreshold,
		MinLength:               minLength,
		MaxLength:               maxLength,
		QualPercent:             qualPercent,
		LowComplexity:           lowComplexity,
		ComplexityThreshold:     complexityThreshold,
		TrimPolyG:               trimPolyG,
		TrimPolyX:               trimPolyX,
		PolyGMinLen:             polyGMinLen,
		CutFront:                cutFront,
		CutTail:                 cutTail,
		CutRight:                cutRight,
		CutWindowSize:           cutWindowSize,
		CutMeanQuality:          cutMeanQuality,
		MaxNCount:               maxNCount,
		MaxNPercent:             maxNPercent,
		LengthRequired:          minLength,
		LengthLimit:             maxLength,
		UMILength:               umiLength,
		UMILocation:             umiLocation,
		UMI:                     umiEnable,
		UMILoc:                  umiLoc,
		UMILen:                  umiLen,
		UMIPrefix:               umiPrefix,
		UMISkip:                 umiSkip,
		DupCalcAccuracy:         dupCalcAccuracy,
		Dedup:                   dedup,
		BaseCorrection:          baseCorrection,
		CorrectionThreshold:     correctionThreshold,
		Correction:              correction,
		OverlapRequire:          overlapLenRequire,
		OverlapDiffLimit:        overlapDiffLimit,
		OverlapDiffPercentLimit: overlapDiffPercentLimit,
		MergeOverlap:            mergeOverlap,
		MinOverlap:              minOverlap,
		MaxMismatch:             maxMismatch,
		OverrepAnalysis:         overrepAnalysis,
		OverrepSampling:         overrepSampling,
		SplitNumber:             splitNumber,
		SplitByLines:            splitByLines,
		SplitPrefixDigits:       splitPrefixDigits,
		Threads:                 threads,
		HTMLReport:              htmlReport,
		JSONReport:              jsonReport,
		DetectAdapterPE:         detectAdapterForPE,
		DetectAdapterSE:         detectAdapter, // legacy --detect-adapter triggers SE detection
	}

	// Splitting mode: --split / --split_by_lines route output across
	// numbered files. Validate the knobs the way upstream does.
	splitEnabled := splitNumber > 0 || splitByLines > 0
	if splitEnabled {
		if err := validateSplit(splitNumber, splitByLines, splitPrefixDigits, mergeOverlap); err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v\n", err)
			os.Exit(1)
		}
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

		if splitEnabled {
			stats, err = fastp.ProcessPairedEndSplit(input1, input2, out1File, out2File, encoding, opts)
		} else if merge {
			// Merge mode: merged reads go to --merged_out (or stdout / out1
			// fallback); out1/out2 receive unmerged pairs when given.
			var mergeOutW io.WriteCloser
			if mergedOut != "" {
				mergeOutW, err = iohelper.OpenWriter(mergedOut)
				if err != nil {
					fmt.Fprintf(os.Stderr, "Error creating merged output file: %v\n", err)
					os.Exit(1)
				}
				defer mergeOutW.Close()
			}
			out1W := openOptionalWriter(out1File)
			out2W := openOptionalWriter(out2File)
			defer out1W.Close()
			defer out2W.Close()
			var mw io.Writer
			if mergeOutW != nil {
				mw = mergeOutW
			}
			stats, err = fastp.ProcessPairedEndMerge(input1, input2, out1W, out2W, mw, encoding, opts)
		} else {
			output1, oerr := iohelper.OpenWriter(out1File)
			if oerr != nil {
				fmt.Fprintf(os.Stderr, "Error creating output file 1: %v\n", oerr)
				os.Exit(1)
			}
			defer output1.Close()

			output2, oerr := iohelper.OpenWriter(out2File)
			if oerr != nil {
				fmt.Fprintf(os.Stderr, "Error creating output file 2: %v\n", oerr)
				os.Exit(1)
			}
			defer output2.Close()

			stats, err = fastp.ProcessPairedEnd(input1, input2, output1, output2, encoding, opts)
		}
	} else {
		// Single-end mode
		input, err := iohelper.OpenReader(inputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
			os.Exit(1)
		}
		defer input.Close()

		if splitEnabled {
			stats, err = fastp.ProcessSingleEndSplit(input, outputFile, encoding, opts)
		} else {
			output, oerr := iohelper.OpenWriter(outputFile)
			if oerr != nil {
				fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", oerr)
				os.Exit(1)
			}
			defer output.Close()

			stats, err = fastp.ProcessSingleEnd(input, output, encoding, opts)
		}
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

// validateSplit checks the --split / --split_by_lines / --split_prefix_digits
// inputs the way upstream options.cpp does and reports the first problem.
func validateSplit(splitNumber, splitByLines, digits int, merge bool) error {
	if merge {
		return fmt.Errorf("splitting mode cannot work with merging mode")
	}
	if splitNumber > 0 && splitByLines > 0 {
		return fmt.Errorf("you cannot set both --split and --split_by_lines, please choose either")
	}
	if digits < 0 || digits > 10 {
		return fmt.Errorf("--split_prefix_digits should be 0~10")
	}
	if splitNumber > 0 && (splitNumber < 2 || splitNumber >= 1000) {
		return fmt.Errorf("--split (number of files) should be 2~999")
	}
	if splitByLines > 0 {
		if splitByLines%4 != 0 {
			return fmt.Errorf("--split_by_lines should be a multiple of 4")
		}
		if splitByLines < 1000 {
			return fmt.Errorf("--split_by_lines should be >= 1000")
		}
	}
	return nil
}

// openOptionalWriter opens path for writing, or returns a no-op writer
// when path is empty. Used in merge mode where --out1/--out2 are optional
// (they only carry unmerged pairs).
func openOptionalWriter(path string) io.WriteCloser {
	if path == "" {
		return nopWriteCloser{io.Discard}
	}
	w, err := iohelper.OpenWriter(path)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file %q: %v\n", path, err)
		os.Exit(1)
	}
	return w
}

// nopWriteCloser adds a no-op Close to an io.Writer.
type nopWriteCloser struct{ io.Writer }

func (nopWriteCloser) Close() error { return nil }

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
