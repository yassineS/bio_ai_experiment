package main

import (
	"encoding/json"
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
		windowSize     int
		noFivePrime    bool
		truncateN      bool
		quiet          bool
		jsonOutput     string
		htmlReport     string
		progress       bool
		autoDetect     bool
		recalibrate    bool
	)
	
	cliflag.StringVar(fs, &fastqFile, "f", "fastq-file", "", "Input FASTQ file (required)")
	cliflag.StringVar(fs, &outputFile, "o", "output-file", "", "Output trimmed file (default: stdout)")
	cliflag.StringVar(fs, &qualType, "t", "qual-type", "sanger", "Quality type: sanger, illumina, solexa (default: sanger)")
	cliflag.IntVar(fs, &qualThreshold, "q", "qual-threshold", 20, "Threshold for trimming (default: 20)")
	cliflag.IntVar(fs, &lengthThreshold, "l", "length-threshold", 20, "Minimum length to keep (default: 20)")
	cliflag.IntVar(fs, &windowSize, "w", "window-size", 10, "Window size for quality assessment (default: 10)")
	cliflag.BoolVar(fs, &noFivePrime, "x", "no-fiveprime", false, "Don't trim 5' end")
	cliflag.BoolVar(fs, &truncateN, "n", "trunc-n", false, "Truncate sequences at position of first N")
	cliflag.BoolVar(fs, &quiet, "", "quiet", false, "Don't print statistics")
	cliflag.StringVar(fs, &jsonOutput, "", "json", "", "Output statistics in JSON format to file")
	cliflag.StringVar(fs, &htmlReport, "", "html", "", "Generate HTML report to file")
	cliflag.BoolVar(fs, &progress, "", "progress", false, "Show progress reporting")
	cliflag.BoolVar(fs, &autoDetect, "", "auto-detect", false, "Auto-detect quality encoding")
	cliflag.BoolVar(fs, &recalibrate, "", "recalibrate", false, "Recalibrate quality scores")
	
	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: sickle se [options]\n\n")
		fmt.Fprintf(os.Stderr, "Options:\n")
		fmt.Fprintf(os.Stderr, "  -f, --fastq-file FILE       Input FASTQ file (required)\n")
		fmt.Fprintf(os.Stderr, "  -o, --output-file FILE      Output trimmed file (default: stdout)\n")
		fmt.Fprintf(os.Stderr, "  -t, --qual-type TYPE        Quality type: sanger, illumina, solexa (default: sanger)\n")
		fmt.Fprintf(os.Stderr, "  -q, --qual-threshold INT    Threshold for trimming (default: 20)\n")
		fmt.Fprintf(os.Stderr, "  -l, --length-threshold INT  Minimum length to keep (default: 20)\n")
		fmt.Fprintf(os.Stderr, "  -w, --window-size INT       Window size for quality assessment (default: 10)\n")
		fmt.Fprintf(os.Stderr, "  -x, --no-fiveprime          Don't trim 5' end\n")
		fmt.Fprintf(os.Stderr, "  -n, --trunc-n               Truncate sequences at position of first N\n")
		fmt.Fprintf(os.Stderr, "  --quiet                     Don't print statistics\n")
		fmt.Fprintf(os.Stderr, "  --json FILE                 Output statistics in JSON format to file\n")
		fmt.Fprintf(os.Stderr, "  --html FILE                 Generate HTML report to file\n")
		fmt.Fprintf(os.Stderr, "  --progress                  Show progress reporting\n")
		fmt.Fprintf(os.Stderr, "  --auto-detect               Auto-detect quality encoding\n")
		fmt.Fprintf(os.Stderr, "  --recalibrate               Recalibrate quality scores\n")
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
	var encoding fastq.QualityEncoding
	if autoDetect {
		detected, err := detectQualityEncoding(fastqFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error auto-detecting quality encoding: %v\n", err)
			os.Exit(1)
		}
		encoding = detected
		if !quiet {
			encodingName := "Phred+33 (Sanger)"
			if detected == fastq.Phred64 {
				encodingName = "Phred+64 (Illumina)"
			}
			fmt.Fprintf(os.Stderr, "Auto-detected quality encoding: %s\n", encodingName)
		}
	} else {
		encoding = getQualityEncoding(qualType)
	}
	
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
		WindowSize:      windowSize,
		Progress:        progress,
		Recalibrate:     recalibrate,
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
	cliflag.StringVar(fs, &qualType, "t", "qual-type", "sanger", "Quality type: sanger, illumina, solexa (default: sanger)")
	cliflag.IntVar(fs, &qualThreshold, "q", "qual-threshold", 20, "Threshold for trimming (default: 20)")
	cliflag.IntVar(fs, &lengthThreshold, "l", "length-threshold", 20, "Minimum length to keep (default: 20)")
	cliflag.IntVar(fs, &windowSize, "w", "window-size", 10, "Window size for quality assessment (default: 10)")
	cliflag.BoolVar(fs, &noFivePrime, "x", "no-fiveprime", false, "Don't trim 5' end")
	cliflag.BoolVar(fs, &truncateN, "n", "trunc-n", false, "Truncate sequences at position of first N")
	cliflag.BoolVar(fs, &quiet, "", "quiet", false, "Don't print statistics")
	cliflag.StringVar(fs, &jsonOutput, "", "json", "", "Output statistics in JSON format to file")
	cliflag.StringVar(fs, &htmlReport, "", "html", "", "Generate HTML report to file")
	cliflag.BoolVar(fs, &progress, "", "progress", false, "Show progress reporting")
	cliflag.BoolVar(fs, &autoDetect, "", "auto-detect", false, "Auto-detect quality encoding")
	cliflag.BoolVar(fs, &recalibrate, "", "recalibrate", false, "Recalibrate quality scores")
	
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
		fmt.Fprintf(os.Stderr, "  -w, --window-size INT       Window size for quality assessment (default: 10)\n")
		fmt.Fprintf(os.Stderr, "  -x, --no-fiveprime          Don't trim 5' end\n")
		fmt.Fprintf(os.Stderr, "  -n, --trunc-n               Truncate sequences at position of first N\n")
		fmt.Fprintf(os.Stderr, "  --quiet                     Don't print statistics\n")
		fmt.Fprintf(os.Stderr, "  --json FILE                 Output statistics in JSON format to file\n")
		fmt.Fprintf(os.Stderr, "  --html FILE                 Generate HTML report to file\n")
		fmt.Fprintf(os.Stderr, "  --progress                  Show progress reporting\n")
		fmt.Fprintf(os.Stderr, "  --auto-detect               Auto-detect quality encoding\n")
		fmt.Fprintf(os.Stderr, "  --recalibrate               Recalibrate quality scores\n")
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
	var encoding fastq.QualityEncoding
	if autoDetect {
		detected, err := detectQualityEncoding(fastqFile1)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error auto-detecting quality encoding: %v\n", err)
			os.Exit(1)
		}
		encoding = detected
		if !quiet {
			encodingName := "Phred+33 (Sanger)"
			if detected == fastq.Phred64 {
				encodingName = "Phred+64 (Illumina)"
			}
			fmt.Fprintf(os.Stderr, "Auto-detected quality encoding: %s\n", encodingName)
		}
	} else {
		encoding = getQualityEncoding(qualType)
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
		WindowSize:      windowSize,
		Progress:        progress,
		Recalibrate:     recalibrate,
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

// detectQualityEncoding automatically detects the quality encoding from a FASTQ file
func detectQualityEncoding(filename string) (fastq.QualityEncoding, error) {
	f, err := iohelper.OpenReader(filename)
	if err != nil {
		return fastq.Phred33, err
	}
	defer f.Close()
	
	// Read a sample of records to detect encoding
	reader := fastq.NewReader(f, fastq.Phred33) // Start with Phred33 for reading
	minQual := 255
	maxQual := 0
	samplesRead := 0
	maxSamples := 10000 // Sample first 10k reads
	
	for samplesRead < maxSamples {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fastq.Phred33, err
		}
		
		for _, q := range record.Quality {
			if int(q) < minQual {
				minQual = int(q)
			}
			if int(q) > maxQual {
				maxQual = int(q)
			}
		}
		samplesRead++
	}
	
	// Phred+33 range: 33-126 (quality 0-93)
	// Phred+64 range: 64-126 (quality 0-62)
	// If we see chars < 59 (which would be quality -5 in Phred+64), it's Phred+33
	// If all chars >= 64, it's likely Phred+64
	if minQual < 59 {
		return fastq.Phred33, nil
	}
	return fastq.Phred64, nil
}

// saveJSONStats saves trimming statistics to a JSON file
func saveJSONStats(stats *sickle.TrimStats, filename string) error {
	type JSONStats struct {
		TotalReads     int     `json:"total_reads"`
		TrimmedReads   int     `json:"trimmed_reads"`
		TrimmedPercent float64 `json:"trimmed_percent"`
		DiscardedReads int     `json:"discarded_reads"`
		DiscardedPercent float64 `json:"discarded_percent"`
		KeptReads      int     `json:"kept_reads"`
		KeptPercent    float64 `json:"kept_percent"`
		TotalBases     int64   `json:"total_bases"`
		TrimmedBases   int64   `json:"trimmed_bases"`
		TrimmedBasesPercent float64 `json:"trimmed_bases_percent"`
	}
	
	jsonStats := JSONStats{
		TotalReads:     stats.TotalReads,
		TrimmedReads:   stats.TrimmedReads,
		TrimmedPercent: 100.0 * float64(stats.TrimmedReads) / float64(stats.TotalReads),
		DiscardedReads: stats.DiscardedReads,
		DiscardedPercent: 100.0 * float64(stats.DiscardedReads) / float64(stats.TotalReads),
		KeptReads:      stats.TotalReads - stats.DiscardedReads,
		KeptPercent:    100.0 * float64(stats.TotalReads-stats.DiscardedReads) / float64(stats.TotalReads),
		TotalBases:     stats.TotalBases,
		TrimmedBases:   stats.TrimmedBases,
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
