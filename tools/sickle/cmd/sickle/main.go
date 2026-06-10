package main

import (
	"bufio"
	"compress/gzip"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/sickle/pkg/sickle"
)

const usage = `sickle - A windowed adaptive trimming tool for FASTQ files using quality scores

Usage:
  sickle <command> [options]

Commands:
  se    Trim single-end reads
  pe    Trim paired-end reads
  batch Trim multiple files in parallel

For command-specific help:
  sickle se -h
  sickle pe -h
  sickle batch -h

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
	case "batch":
		runBatch()
	case "-h", "--help", "help":
		fmt.Print(usage)
		os.Exit(0)
	case "-v", "--version", "version":
		fmt.Println(sickleVersion)
		os.Exit(0)
	default:
		fmt.Fprintf(os.Stderr, "Error: unknown command %q\n\n", command)
		fmt.Fprint(os.Stderr, usage)
		os.Exit(1)
	}
}

// gzipWriteCloser wraps a destination WriteCloser so that data written to it is
// gzip-compressed. Closing it flushes and closes the gzip stream and then the
// underlying destination, mirroring upstream sickle's -g/--gzip-output.
type gzipWriteCloser struct {
	gz   *gzip.Writer
	dest io.WriteCloser
}

func (g *gzipWriteCloser) Write(p []byte) (int, error) { return g.gz.Write(p) }

func (g *gzipWriteCloser) Close() error {
	if err := g.gz.Close(); err != nil {
		g.dest.Close()
		return err
	}
	return g.dest.Close()
}

// maybeGzip wraps w in a gzip compressor when on is true; otherwise it returns w
// unchanged. It lets the -g/--gzip-output compatibility flag force gzip output
// regardless of the output filename's extension, as upstream sickle does.
func maybeGzip(w io.WriteCloser, on bool) io.WriteCloser {
	if !on {
		return w
	}
	return &gzipWriteCloser{gz: gzip.NewWriter(w), dest: w}
}

// mustOpenWriter opens filename for writing (with transparent gzip-by-extension
// via iohelper) and exits with an error message tagged by label on failure.
func mustOpenWriter(filename, label string) io.WriteCloser {
	w, err := iohelper.OpenWriter(filename)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating %s file: %v\n", label, err)
		os.Exit(1)
	}
	return w
}

func runSingleEnd() {
	fs := flag.NewFlagSet("sickle se", flag.ExitOnError)
	var hv subFlags
	hv.register(fs)

	var (
		fastqFile       string
		outputFile      string
		qualType        string
		qualThreshold   int
		lengthThreshold int
		windowSize      int
		noFivePrime     bool
		truncateN       bool
		quiet           bool
		gzipOutput      bool
		debugCompat     bool
		jsonOutput      string
		htmlReport      string
		progress        bool
		autoDetect      bool
		recalibrate     bool
	)

	cliflag.StringVar(fs, &fastqFile, "f", "fastq-file", "", "Input FASTQ file (required)")
	cliflag.StringVar(fs, &outputFile, "o", "output-file", "", "Output trimmed file (default: stdout)")
	cliflag.StringVar(fs, &qualType, "t", "qual-type", "auto", "Quality type: auto, sanger, illumina, solexa (default: auto)")
	cliflag.IntVar(fs, &qualThreshold, "q", "qual-threshold", 20, "Threshold for trimming (default: 20)")
	cliflag.IntVar(fs, &lengthThreshold, "l", "length-threshold", 20, "Minimum length to keep (default: 20)")
	cliflag.IntVar(fs, &windowSize, "w", "window-size", 10, "Window size for quality assessment (default: 10)")
	cliflag.BoolVar(fs, &noFivePrime, "x", "no-fiveprime", false, "Don't trim 5' end")
	cliflag.BoolVar(fs, &truncateN, "n", "trunc-n", false, "Truncate sequences at position of first N")
	cliflag.BoolVar(fs, &quiet, "z", "quiet", false, "Don't print statistics")
	cliflag.BoolVar(fs, &gzipOutput, "g", "gzip-output", false, "Gzip-compress the output stream (upstream -g)")
	cliflag.BoolVar(fs, &debugCompat, "d", "debug", false, "Upstream debug flag, accepted for compatibility (no effect)")
	cliflag.StringVar(fs, &jsonOutput, "", "json", "", "Output statistics in JSON format to file")
	cliflag.StringVar(fs, &htmlReport, "", "html", "", "Generate HTML report to file")
	cliflag.BoolVar(fs, &progress, "", "progress", false, "Show progress reporting")
	cliflag.BoolVar(fs, &autoDetect, "", "auto-detect", false, "Force auto-detection of quality encoding (same as -t auto)")
	cliflag.BoolVar(fs, &recalibrate, "", "recalibrate", false, "Recalibrate quality scores")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: sickle se [options]\n\n")
		fmt.Fprintf(os.Stderr, "Options:\n")
		fmt.Fprintf(os.Stderr, "  -f, --fastq-file FILE       Input FASTQ file (required)\n")
		fmt.Fprintf(os.Stderr, "  -o, --output-file FILE      Output trimmed file (default: stdout)\n")
		fmt.Fprintf(os.Stderr, "  -t, --qual-type TYPE        Quality type: auto, sanger, illumina, solexa (default: auto)\n")
		fmt.Fprintf(os.Stderr, "  -q, --qual-threshold INT    Threshold for trimming (default: 20)\n")
		fmt.Fprintf(os.Stderr, "  -l, --length-threshold INT  Minimum length to keep (default: 20)\n")
		fmt.Fprintf(os.Stderr, "  -w, --window-size INT       Window size for quality assessment (default: 10)\n")
		fmt.Fprintf(os.Stderr, "  -x, --no-fiveprime          Don't trim 5' end\n")
		fmt.Fprintf(os.Stderr, "  -n, --trunc-n               Truncate sequences at position of first N\n")
		fmt.Fprintf(os.Stderr, "  -z, --quiet                 Don't print statistics\n")
		fmt.Fprintf(os.Stderr, "  -g, --gzip-output           Gzip-compress the output stream (upstream -g)\n")
		fmt.Fprintf(os.Stderr, "  -d, --debug                 Upstream debug flag, accepted for compatibility (no effect)\n")
		fmt.Fprintf(os.Stderr, "  --json FILE                 Output statistics in JSON format to file\n")
		fmt.Fprintf(os.Stderr, "  --html FILE                 Generate HTML report to file\n")
		fmt.Fprintf(os.Stderr, "  --progress                  Show progress reporting\n")
		fmt.Fprintf(os.Stderr, "  --auto-detect               Force auto-detection of quality encoding (same as -t auto)\n")
		fmt.Fprintf(os.Stderr, "  --recalibrate               Recalibrate quality scores\n")
		fmt.Fprintf(os.Stderr, "\nExample:\n")
		fmt.Fprintf(os.Stderr, "  sickle se -f input.fastq -o output.fastq -q 20 -l 20\n")
	}

	if err := cliflag.Parse(fs, os.Args[2:]); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(2)
	}

	hv.handle(fs)

	// Validate required arguments
	if fastqFile == "" {
		fmt.Fprintln(os.Stderr, "Error: -f/--fastq-file is required")
		fs.Usage()
		os.Exit(1)
	}

	// Open input file (with automatic gzip support).
	inputFile, err := iohelper.OpenReader(fastqFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer inputFile.Close()

	// Wrap in a bufio.Reader large enough that DetectEncoding's Peek can find
	// ~10000 quality characters without consuming the underlying stream. We
	// keep using this *bufio.Reader for the actual trimming, so any bytes we
	// peeked at remain available for reading.
	bufInput := bufio.NewReaderSize(inputFile, 256*1024)

	// Determine quality encoding. The CLI accepts an explicit name (sanger,
	// illumina, solexa, phred33, phred64) or "auto" to detect from the input.
	// The legacy --auto-detect boolean is honored as an override of -t.
	encoding, err := resolveEncoding(qualType, autoDetect, bufInput, quiet)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	// Open output file (with automatic gzip support)
	outFileName := outputFile
	if outFileName == "" {
		outFileName = "-"
	}
	rawOut, err := iohelper.OpenWriter(outFileName)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	outFile := maybeGzip(rawOut, gzipOutput)
	defer outFile.Close()

	// The upstream -d/--debug flag only printed internal cut coordinates; it is
	// accepted for command-line compatibility but has no effect here.
	_ = debugCompat

	// Set up trim options
	opts := sickle.TrimOptions{
		QualThreshold:   qualThreshold,
		LengthThreshold: lengthThreshold,
		NoFivePrime:     noFivePrime,
		TruncateN:       truncateN,
		WindowSize:      windowSize,
		Progress:        progress,
		Recalibrate:     recalibrate,
	}

	// Perform trimming
	stats, err := sickle.TrimSingleEnd(bufInput, outFile, encoding, opts)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error during trimming: %v\n", err)
		os.Exit(1)
	}

	// Print statistics
	if !quiet {
		printStats(stats, "SE")
	}

	// Save JSON statistics
	if jsonOutput != "" {
		if err := saveJSONStats(stats, jsonOutput); err != nil {
			fmt.Fprintf(os.Stderr, "Error saving JSON statistics: %v\n", err)
			os.Exit(1)
		}
	}

	// Generate HTML report
	if htmlReport != "" {
		if err := generateHTMLReport(stats, htmlReport, "SE"); err != nil {
			fmt.Fprintf(os.Stderr, "Error generating HTML report: %v\n", err)
			os.Exit(1)
		}
	}
}

func runPairedEnd() {
	fs := flag.NewFlagSet("sickle pe", flag.ExitOnError)
	var hv subFlags
	hv.register(fs)

	var (
		fastqFile1      string
		fastqFile2      string
		outputFile1     string
		outputFile2     string
		outputSingle    string
		qualType        string
		qualThreshold   int
		lengthThreshold int
		windowSize      int
		noFivePrime     bool
		truncateN       bool
		quiet           bool
		gzipOutput      bool
		debugCompat     bool
		jsonOutput      string
		htmlReport      string
		progress        bool
		autoDetect      bool
		recalibrate     bool
	)

	cliflag.StringVar(fs, &fastqFile1, "f", "fastq-file", "", "First input FASTQ file (required)")
	cliflag.StringVar(fs, &fastqFile2, "r", "reverse-file", "", "Second input FASTQ file (required)")
	cliflag.StringVar(fs, &outputFile1, "o", "output-file", "", "First output trimmed file (required)")
	cliflag.StringVar(fs, &outputFile2, "p", "output-paired", "", "Second output trimmed file (required)")
	cliflag.StringVar(fs, &outputSingle, "s", "output-single", "", "Output single-end reads (optional)")
	cliflag.StringVar(fs, &qualType, "t", "qual-type", "auto", "Quality type: auto, sanger, illumina, solexa (default: auto)")
	cliflag.IntVar(fs, &qualThreshold, "q", "qual-threshold", 20, "Threshold for trimming (default: 20)")
	cliflag.IntVar(fs, &lengthThreshold, "l", "length-threshold", 20, "Minimum length to keep (default: 20)")
	cliflag.IntVar(fs, &windowSize, "w", "window-size", 10, "Window size for quality assessment (default: 10)")
	cliflag.BoolVar(fs, &noFivePrime, "x", "no-fiveprime", false, "Don't trim 5' end")
	cliflag.BoolVar(fs, &truncateN, "n", "trunc-n", false, "Truncate sequences at position of first N")
	cliflag.BoolVar(fs, &quiet, "z", "quiet", false, "Don't print statistics")
	cliflag.BoolVar(fs, &gzipOutput, "g", "gzip-output", false, "Gzip-compress the output streams (upstream -g)")
	cliflag.BoolVar(fs, &debugCompat, "d", "debug", false, "Upstream debug flag, accepted for compatibility (no effect)")
	cliflag.StringVar(fs, &jsonOutput, "", "json", "", "Output statistics in JSON format to file")
	cliflag.StringVar(fs, &htmlReport, "", "html", "", "Generate HTML report to file")
	cliflag.BoolVar(fs, &progress, "", "progress", false, "Show progress reporting")
	cliflag.BoolVar(fs, &autoDetect, "", "auto-detect", false, "Force auto-detection of quality encoding (same as -t auto)")
	cliflag.BoolVar(fs, &recalibrate, "", "recalibrate", false, "Recalibrate quality scores")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: sickle pe [options]\n\n")
		fmt.Fprintf(os.Stderr, "Options:\n")
		fmt.Fprintf(os.Stderr, "  -f, --fastq-file FILE       First input FASTQ file (required)\n")
		fmt.Fprintf(os.Stderr, "  -r, --reverse-file FILE     Second input FASTQ file (required)\n")
		fmt.Fprintf(os.Stderr, "  -o, --output-file FILE      First output trimmed file (required)\n")
		fmt.Fprintf(os.Stderr, "  -p, --output-paired FILE    Second output trimmed file (required)\n")
		fmt.Fprintf(os.Stderr, "  -s, --output-single FILE    Output single-end reads (optional)\n")
		fmt.Fprintf(os.Stderr, "  -t, --qual-type TYPE        Quality type: auto, sanger, illumina, solexa (default: auto)\n")
		fmt.Fprintf(os.Stderr, "  -q, --qual-threshold INT    Threshold for trimming (default: 20)\n")
		fmt.Fprintf(os.Stderr, "  -l, --length-threshold INT  Minimum length to keep (default: 20)\n")
		fmt.Fprintf(os.Stderr, "  -w, --window-size INT       Window size for quality assessment (default: 10)\n")
		fmt.Fprintf(os.Stderr, "  -x, --no-fiveprime          Don't trim 5' end\n")
		fmt.Fprintf(os.Stderr, "  -n, --trunc-n               Truncate sequences at position of first N\n")
		fmt.Fprintf(os.Stderr, "  -z, --quiet                 Don't print statistics\n")
		fmt.Fprintf(os.Stderr, "  -g, --gzip-output           Gzip-compress the output streams (upstream -g)\n")
		fmt.Fprintf(os.Stderr, "  -d, --debug                 Upstream debug flag, accepted for compatibility (no effect)\n")
		fmt.Fprintf(os.Stderr, "  --json FILE                 Output statistics in JSON format to file\n")
		fmt.Fprintf(os.Stderr, "  --html FILE                 Generate HTML report to file\n")
		fmt.Fprintf(os.Stderr, "  --progress                  Show progress reporting\n")
		fmt.Fprintf(os.Stderr, "  --auto-detect               Force auto-detection of quality encoding (same as -t auto)\n")
		fmt.Fprintf(os.Stderr, "  --recalibrate               Recalibrate quality scores\n")
		fmt.Fprintf(os.Stderr, "\nExample:\n")
		fmt.Fprintf(os.Stderr, "  sickle pe -f input1.fastq -r input2.fastq -o output1.fastq -p output2.fastq -s singles.fastq\n")
		fmt.Fprintf(os.Stderr, "\nNote: when --qual-type is auto, the encoding is detected from R1 and applied\n")
		fmt.Fprintf(os.Stderr, "to both R1 and R2. Mismatched encodings between mates are not handled separately.\n")
	}

	if err := cliflag.Parse(fs, os.Args[2:]); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(2)
	}

	hv.handle(fs)

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

	// Wrap both inputs in bufio.Reader so DetectEncoding can Peek at R1 without
	// consuming bytes that the trimmer still needs to read. R2 is wrapped for
	// symmetry — the encoding detected on R1 is applied to both files.
	bufF1 := bufio.NewReaderSize(f1, 256*1024)
	bufF2 := bufio.NewReaderSize(f2, 256*1024)

	// Determine quality encoding. The CLI accepts an explicit name (sanger,
	// illumina, solexa, phred33, phred64) or "auto" to detect from R1. The
	// legacy --auto-detect boolean is honored as an override of -t.
	encoding, err := resolveEncoding(qualType, autoDetect, bufF1, quiet)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	// The upstream -d/--debug flag only printed internal cut coordinates; it is
	// accepted for command-line compatibility but has no effect here.
	_ = debugCompat

	// Open output files (with automatic gzip support)
	raw1, err := iohelper.OpenWriter(outputFile1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating first output file: %v\n", err)
		os.Exit(1)
	}
	out1 := maybeGzip(raw1, gzipOutput)
	defer out1.Close()

	raw2, err := iohelper.OpenWriter(outputFile2)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating second output file: %v\n", err)
		os.Exit(1)
	}
	out2 := maybeGzip(raw2, gzipOutput)
	defer out2.Close()

	// Open optional single output file (with automatic gzip support)
	var outSingle io.Writer
	if outputSingle != "" {
		f := maybeGzip(mustOpenWriter(outputSingle, "single output"), gzipOutput)
		defer f.Close()
		outSingle = f
	}

	// Set up trim options
	opts := sickle.TrimOptions{
		QualThreshold:   qualThreshold,
		LengthThreshold: lengthThreshold,
		NoFivePrime:     noFivePrime,
		TruncateN:       truncateN,
		WindowSize:      windowSize,
		Progress:        progress,
		Recalibrate:     recalibrate,
	}

	// Perform trimming
	stats, err := sickle.TrimPairedEnd(bufF1, bufF2, out1, out2, outSingle, encoding, opts)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error during trimming: %v\n", err)
		os.Exit(1)
	}

	// Print statistics
	if !quiet {
		printStats(stats, "PE")
	}

	// Save JSON statistics
	if jsonOutput != "" {
		if err := saveJSONStats(stats, jsonOutput); err != nil {
			fmt.Fprintf(os.Stderr, "Error saving JSON statistics: %v\n", err)
			os.Exit(1)
		}
	}

	// Generate HTML report
	if htmlReport != "" {
		if err := generateHTMLReport(stats, htmlReport, "PE"); err != nil {
			fmt.Fprintf(os.Stderr, "Error generating HTML report: %v\n", err)
			os.Exit(1)
		}
	}
}

func runBatch() {
	fs := flag.NewFlagSet("sickle batch", flag.ExitOnError)
	var hv subFlags
	hv.register(fs)

	var (
		fileList        string
		outputDir       string
		qualType        string
		qualThreshold   int
		lengthThreshold int
		windowSize      int
		noFivePrime     bool
		truncateN       bool
		quiet           bool
		jsonOutput      bool
		htmlReport      bool
		progress        bool
		autoDetect      bool
		recalibrate     bool
		numWorkers      int
	)

	cliflag.StringVar(fs, &fileList, "i", "input-list", "", "File containing list of input FASTQ files (one per line)")
	cliflag.StringVar(fs, &outputDir, "o", "output-dir", ".", "Output directory for trimmed files")
	cliflag.StringVar(fs, &qualType, "t", "qual-type", "auto", "Quality type: auto, sanger, illumina, solexa (default: auto)")
	cliflag.IntVar(fs, &qualThreshold, "q", "qual-threshold", 20, "Threshold for trimming (default: 20)")
	cliflag.IntVar(fs, &lengthThreshold, "l", "length-threshold", 20, "Minimum length to keep (default: 20)")
	cliflag.IntVar(fs, &windowSize, "w", "window-size", 10, "Window size for quality assessment (default: 10)")
	cliflag.IntVar(fs, &numWorkers, "j", "jobs", 4, "Number of parallel workers (default: 4)")
	cliflag.BoolVar(fs, &noFivePrime, "x", "no-fiveprime", false, "Don't trim 5' end")
	cliflag.BoolVar(fs, &truncateN, "n", "trunc-n", false, "Truncate sequences at position of first N")
	cliflag.BoolVar(fs, &quiet, "z", "quiet", false, "Don't print statistics")
	cliflag.BoolVar(fs, &jsonOutput, "", "json", false, "Output statistics in JSON format for each file")
	cliflag.BoolVar(fs, &htmlReport, "", "html", false, "Generate HTML report for each file")
	cliflag.BoolVar(fs, &progress, "", "progress", false, "Show progress reporting")
	cliflag.BoolVar(fs, &autoDetect, "", "auto-detect", false, "Force auto-detection of quality encoding (same as -t auto)")
	cliflag.BoolVar(fs, &recalibrate, "", "recalibrate", false, "Recalibrate quality scores")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: sickle batch [options]\n\n")
		fmt.Fprintf(os.Stderr, "Batch mode processes multiple FASTQ files in parallel.\n\n")
		fmt.Fprintf(os.Stderr, "Options:\n")
		fmt.Fprintf(os.Stderr, "  -i, --input-list FILE       File containing list of input FASTQ files (required)\n")
		fmt.Fprintf(os.Stderr, "  -o, --output-dir DIR        Output directory for trimmed files (default: .)\n")
		fmt.Fprintf(os.Stderr, "  -t, --qual-type TYPE        Quality type: auto, sanger, illumina, solexa (default: auto)\n")
		fmt.Fprintf(os.Stderr, "  -q, --qual-threshold INT    Threshold for trimming (default: 20)\n")
		fmt.Fprintf(os.Stderr, "  -l, --length-threshold INT  Minimum length to keep (default: 20)\n")
		fmt.Fprintf(os.Stderr, "  -w, --window-size INT       Window size for quality assessment (default: 10)\n")
		fmt.Fprintf(os.Stderr, "  -j, --jobs INT              Number of parallel workers (default: 4)\n")
		fmt.Fprintf(os.Stderr, "  -x, --no-fiveprime          Don't trim 5' end\n")
		fmt.Fprintf(os.Stderr, "  -n, --trunc-n               Truncate sequences at position of first N\n")
		fmt.Fprintf(os.Stderr, "  --quiet                     Don't print statistics\n")
		fmt.Fprintf(os.Stderr, "  --json                      Output statistics in JSON format for each file\n")
		fmt.Fprintf(os.Stderr, "  --html                      Generate HTML report for each file\n")
		fmt.Fprintf(os.Stderr, "  --progress                  Show progress reporting\n")
		fmt.Fprintf(os.Stderr, "  --auto-detect               Force auto-detection of quality encoding (same as -t auto)\n")
		fmt.Fprintf(os.Stderr, "  --recalibrate               Recalibrate quality scores\n")
		fmt.Fprintf(os.Stderr, "\nExample:\n")
		fmt.Fprintf(os.Stderr, "  sickle batch -i files.txt -o trimmed_output -j 8\n")
		fmt.Fprintf(os.Stderr, "\nThe input list file should contain one FASTQ file path per line.\n")
	}

	if err := cliflag.Parse(fs, os.Args[2:]); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(2)
	}

	hv.handle(fs)

	// Validate required arguments
	if fileList == "" {
		fmt.Fprintln(os.Stderr, "Error: -i/--input-list is required")
		fs.Usage()
		os.Exit(1)
	}

	// Create output directory if it doesn't exist
	if err := os.MkdirAll(outputDir, 0755); err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output directory: %v\n", err)
		os.Exit(1)
	}

	// Read file list
	f, err := os.Open(fileList)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening file list: %v\n", err)
		os.Exit(1)
	}
	defer f.Close()

	var files []string
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line != "" && !strings.HasPrefix(line, "#") {
			files = append(files, line)
		}
	}
	if err := scanner.Err(); err != nil {
		fmt.Fprintf(os.Stderr, "Error reading file list: %v\n", err)
		os.Exit(1)
	}

	if len(files) == 0 {
		fmt.Fprintln(os.Stderr, "Error: no files found in input list")
		os.Exit(1)
	}

	fmt.Fprintf(os.Stderr, "Processing %d files with %d workers...\n", len(files), numWorkers)

	// Set up worker pool
	type job struct {
		inputFile string
		index     int
	}

	type result struct {
		inputFile  string
		outputFile string
		stats      *sickle.TrimStats
		err        error
	}

	jobs := make(chan job, len(files))
	results := make(chan result, len(files))

	var wg sync.WaitGroup

	// Start workers
	for w := 0; w < numWorkers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j := range jobs {
				// Generate output filename
				baseName := filepath.Base(j.inputFile)
				outputFile := filepath.Join(outputDir, "trimmed_"+baseName)

				// Open files
				input, err := iohelper.OpenReader(j.inputFile)
				if err != nil {
					results <- result{inputFile: j.inputFile, err: err}
					continue
				}

				// Wrap in a bufio.Reader so we can Peek at quality bytes for
				// auto-detection without consuming them — the trimmer reads
				// from the same buffered reader and re-uses the peeked bytes.
				bufInput := bufio.NewReaderSize(input, 256*1024)

				// Determine quality encoding (per-file, since different files
				// in a batch may have different encodings).
				encoding, encErr := resolveEncoding(qualType, autoDetect, bufInput, true /* quiet inside workers */)
				if encErr != nil {
					input.Close()
					results <- result{inputFile: j.inputFile, err: encErr}
					continue
				}

				output, err := iohelper.OpenWriter(outputFile)
				if err != nil {
					input.Close()
					results <- result{inputFile: j.inputFile, err: err}
					continue
				}

				// Set up trim options
				opts := sickle.TrimOptions{
					QualThreshold:   qualThreshold,
					LengthThreshold: lengthThreshold,
					NoFivePrime:     noFivePrime,
					TruncateN:       truncateN,
					WindowSize:      windowSize,
					Progress:        progress,
					Recalibrate:     recalibrate,
				}

				// Perform trimming
				stats, err := sickle.TrimSingleEnd(bufInput, output, encoding, opts)
				input.Close()
				output.Close()

				if err != nil {
					results <- result{inputFile: j.inputFile, outputFile: outputFile, err: err}
					continue
				}

				// Save JSON if requested
				if jsonOutput {
					jsonFile := outputFile + ".json"
					if err := saveJSONStats(stats, jsonFile); err != nil {
						fmt.Fprintf(os.Stderr, "Warning: failed to save JSON for %s: %v\n", j.inputFile, err)
					}
				}

				// Generate HTML if requested
				if htmlReport {
					htmlFile := outputFile + ".html"
					if err := generateHTMLReport(stats, htmlFile, "SE"); err != nil {
						fmt.Fprintf(os.Stderr, "Warning: failed to generate HTML for %s: %v\n", j.inputFile, err)
					}
				}

				results <- result{inputFile: j.inputFile, outputFile: outputFile, stats: stats, err: nil}
			}
		}()
	}

	// Send jobs
	for i, file := range files {
		jobs <- job{inputFile: file, index: i}
	}
	close(jobs)

	// Wait for all workers to finish
	go func() {
		wg.Wait()
		close(results)
	}()

	// Collect results
	successCount := 0
	failCount := 0
	for r := range results {
		if r.err != nil {
			fmt.Fprintf(os.Stderr, "Error processing %s: %v\n", r.inputFile, r.err)
			failCount++
		} else {
			if !quiet {
				fmt.Fprintf(os.Stderr, "✓ %s -> %s (kept %d/%d reads)\n",
					r.inputFile, r.outputFile,
					r.stats.TotalReads-r.stats.DiscardedReads, r.stats.TotalReads)
			}
			successCount++
		}
	}

	fmt.Fprintf(os.Stderr, "\nBatch processing complete: %d succeeded, %d failed\n", successCount, failCount)

	if failCount > 0 {
		os.Exit(1)
	}
}

// resolveEncoding turns the user-supplied -t/--qual-type value into a concrete
// fastq.QualityEncoding. When qualType is "auto" (the default) — or the legacy
// --auto-detect boolean flag is set — it peeks at br via sickle.DetectEncoding
// to infer the encoding, then logs a one-line stderr message describing what
// it picked (unless quiet is true). For any other value it normalizes through
// sickle.EncodingFromName.
//
// The bufio.Reader br is only used when detection actually runs; callers must
// continue reading from the same br afterwards because Peek does not consume.
func resolveEncoding(qualType string, autoDetect bool, br *bufio.Reader, quiet bool) (fastq.QualityEncoding, error) {
	if autoDetect || qualType == "" || qualType == "auto" {
		res, err := sickle.DetectEncoding(br)
		if err != nil {
			return fastq.Phred33, fmt.Errorf("auto-detecting quality encoding: %w", err)
		}
		if !quiet {
			label := "Phred+33"
			if res.Encoding == fastq.Phred64 {
				label = "Phred+64"
			}
			suffix := ""
			if res.Ambiguous {
				suffix = " [ambiguous: byte range did not match a single encoding cleanly; defaulted to sanger]"
			}
			fmt.Fprintf(os.Stderr, "sickle: auto-detected quality encoding: %s (%s)%s\n", res.Name, label, suffix)
		}
		return res.Encoding, nil
	}
	enc, err := sickle.EncodingFromName(qualType)
	if err != nil {
		return fastq.Phred33, err
	}
	return enc, nil
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

// saveJSONStats saves trimming statistics to a JSON file
func saveJSONStats(stats *sickle.TrimStats, filename string) error {
	type JSONStats struct {
		TotalReads          int     `json:"total_reads"`
		TrimmedReads        int     `json:"trimmed_reads"`
		TrimmedPercent      float64 `json:"trimmed_percent"`
		DiscardedReads      int     `json:"discarded_reads"`
		DiscardedPercent    float64 `json:"discarded_percent"`
		KeptReads           int     `json:"kept_reads"`
		KeptPercent         float64 `json:"kept_percent"`
		TotalBases          int64   `json:"total_bases"`
		TrimmedBases        int64   `json:"trimmed_bases"`
		TrimmedBasesPercent float64 `json:"trimmed_bases_percent"`
	}

	jsonStats := JSONStats{
		TotalReads:          stats.TotalReads,
		TrimmedReads:        stats.TrimmedReads,
		TrimmedPercent:      100.0 * float64(stats.TrimmedReads) / float64(stats.TotalReads),
		DiscardedReads:      stats.DiscardedReads,
		DiscardedPercent:    100.0 * float64(stats.DiscardedReads) / float64(stats.TotalReads),
		KeptReads:           stats.TotalReads - stats.DiscardedReads,
		KeptPercent:         100.0 * float64(stats.TotalReads-stats.DiscardedReads) / float64(stats.TotalReads),
		TotalBases:          stats.TotalBases,
		TrimmedBases:        stats.TrimmedBases,
		TrimmedBasesPercent: 100.0 * float64(stats.TrimmedBases) / float64(stats.TotalBases),
	}

	data, err := json.MarshalIndent(jsonStats, "", "  ")
	if err != nil {
		return err
	}

	return os.WriteFile(filename, data, 0644)
}

// generateHTMLReport generates an HTML report with trimming statistics
func generateHTMLReport(stats *sickle.TrimStats, filename, mode string) error {
	htmlTemplate := `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Sickle Trimming Report</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            max-width: 900px;
            margin: 40px auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .container {
            background-color: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        h1 {
            color: #333;
            border-bottom: 3px solid #4CAF50;
            padding-bottom: 10px;
        }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 20px;
            margin-top: 20px;
        }
        .stat-box {
            background-color: #f9f9f9;
            padding: 20px;
            border-radius: 4px;
            border-left: 4px solid #4CAF50;
        }
        .stat-label {
            font-size: 14px;
            color: #666;
            margin-bottom: 5px;
        }
        .stat-value {
            font-size: 28px;
            font-weight: bold;
            color: #333;
        }
        .stat-percent {
            font-size: 16px;
            color: #4CAF50;
            margin-left: 10px;
        }
        .progress-bar {
            width: 100%;
            height: 30px;
            background-color: #e0e0e0;
            border-radius: 15px;
            overflow: hidden;
            margin-top: 10px;
        }
        .progress-fill {
            height: 100%;
            background: linear-gradient(90deg, #4CAF50 0%, #45a049 100%);
            transition: width 0.3s ease;
        }
        .footer {
            margin-top: 30px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #666;
            font-size: 12px;
            text-align: center;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>Sickle %s Trimming Report</h1>
        
        <div class="stats-grid">
            <div class="stat-box">
                <div class="stat-label">Total Reads</div>
                <div class="stat-value">%d</div>
            </div>
            
            <div class="stat-box">
                <div class="stat-label">Kept Reads</div>
                <div class="stat-value">%d<span class="stat-percent">(%.2f%%)</span></div>
                <div class="progress-bar">
                    <div class="progress-fill" style="width: %.2f%%"></div>
                </div>
            </div>
            
            <div class="stat-box">
                <div class="stat-label">Trimmed Reads</div>
                <div class="stat-value">%d<span class="stat-percent">(%.2f%%)</span></div>
                <div class="progress-bar">
                    <div class="progress-fill" style="width: %.2f%%"></div>
                </div>
            </div>
            
            <div class="stat-box">
                <div class="stat-label">Discarded Reads</div>
                <div class="stat-value">%d<span class="stat-percent">(%.2f%%)</span></div>
                <div class="progress-bar">
                    <div class="progress-fill" style="width: %.2f%%; background: linear-gradient(90deg, #f44336 0%, #da190b 100%);"></div>
                </div>
            </div>
            
            <div class="stat-box">
                <div class="stat-label">Total Bases</div>
                <div class="stat-value">%d</div>
            </div>
            
            <div class="stat-box">
                <div class="stat-label">Trimmed Bases</div>
                <div class="stat-value">%d<span class="stat-percent">(%.2f%%)</span></div>
                <div class="progress-bar">
                    <div class="progress-fill" style="width: %.2f%%"></div>
                </div>
            </div>
        </div>
        
        <div class="footer">
            Generated by Sickle v1.0.0 (Go implementation)
        </div>
    </div>
</body>
</html>`

	keptReads := stats.TotalReads - stats.DiscardedReads
	keptPercent := 100.0 * float64(keptReads) / float64(stats.TotalReads)
	trimmedPercent := 100.0 * float64(stats.TrimmedReads) / float64(stats.TotalReads)
	discardedPercent := 100.0 * float64(stats.DiscardedReads) / float64(stats.TotalReads)
	trimmedBasesPercent := 100.0 * float64(stats.TrimmedBases) / float64(stats.TotalBases)

	html := fmt.Sprintf(htmlTemplate,
		mode,
		stats.TotalReads,
		keptReads, keptPercent, keptPercent,
		stats.TrimmedReads, trimmedPercent, trimmedPercent,
		stats.DiscardedReads, discardedPercent, discardedPercent,
		stats.TotalBases,
		stats.TrimmedBases, trimmedBasesPercent, trimmedBasesPercent,
	)

	return os.WriteFile(filename, []byte(html), 0644)
}
