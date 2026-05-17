package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"os"
	"strings"
	"time"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/prinseq/pkg/prinseq"
)

// nowFunc is overridable in tests.
var nowFunc = time.Now

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
	case "graph":
		runGraph(os.Args[2:])
	case "graph_data", "graph-data", "graphdata":
		runGraphData(os.Args[2:])
	case "report":
		runReport(os.Args[2:])
	case "benchmark":
		runBenchmark(os.Args[2:])
	case "api":
		runAPI(os.Args[2:])
	case "batch":
		runBatch(os.Args[2:])
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
  stats       Calculate sequence statistics
  filter      Filter sequences based on quality criteria
  graph       Generate quality graphs
  graph_data  Emit upstream prinseq-lite .gd graph-data JSON
  report      Generate HTML quality report
  benchmark   Run performance benchmarks
  api         Start REST API server
  batch       Process multiple files in parallel
  version     Print version information
  help        Print this help message

Use "prinseq <command> -h" for more information about a command.`)
}

func runStats(args []string) {
	fs := flag.NewFlagSet("stats", flag.ExitOnError)

	var fastq, fasta string
	var jsonOutput, enhanced bool
	cliflag.StringVar(fs, &fastq, "", "fastq", "", "Input FASTQ file (use '-' for stdin)")
	cliflag.StringVar(fs, &fasta, "", "fasta", "", "Input FASTA file (use '-' for stdin)")
	cliflag.BoolVar(fs, &jsonOutput, "j", "json", false, "Output statistics in JSON format")
	cliflag.BoolVar(fs, &enhanced, "e", "enhanced", false, "Calculate enhanced statistics (distributions, dinucleotides)")

	fs.Usage = func() {
		fmt.Print(`Usage: prinseq stats [options]

Calculate sequence statistics for FASTA or FASTQ files.

Options:
  --fasta FILE              Input FASTA file (use '-' for stdin)
  --fastq FILE              Input FASTQ file (use '-' for stdin)
  -j, --json                Output statistics in JSON format
  -e, --enhanced            Calculate enhanced statistics
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
	var stats *prinseq.Stats
	var err2 error
	if enhanced {
		stats, err2 = prinseq.CalculateEnhancedStats(reader, isFastq)
	} else {
		stats, err2 = prinseq.CalculateStats(reader, isFastq)
	}

	if err2 != nil {
		fmt.Fprintf(os.Stderr, "Error calculating statistics: %v\n", err2)
		os.Exit(1)
	}

	// Print statistics
	if jsonOutput {
		encoder := json.NewEncoder(os.Stdout)
		encoder.SetIndent("", "  ")
		if err := encoder.Encode(stats); err != nil {
			fmt.Fprintf(os.Stderr, "Error encoding JSON: %v\n", err)
			os.Exit(1)
		}
	} else {
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
}

func runFilter(args []string) {
	fs := flag.NewFlagSet("filter", flag.ExitOnError)

	// Input/output options
	var input1, input2, output1, output2, outBad string
	var fasta, fastq bool

	cliflag.StringVar(fs, &input1, "i", "input", "", "Primary input file (use '-' for stdin)")
	cliflag.StringVar(fs, &input2, "", "input2", "", "Paired-end input file 2")
	cliflag.StringVar(fs, &output1, "o", "output", "", "Output file for filtered sequences (default: stdout)")
	cliflag.StringVar(fs, &output2, "", "output2", "", "Output file for paired-end file 2")
	cliflag.StringVar(fs, &outBad, "", "out-bad", "", "Output file for rejected sequences")
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

	// Upstream-named aliases for the N-filter knobs (single-dash forms
	// from prinseq-lite.pl). Go's `flag` package treats `-foo` and
	// `--foo` identically, so registering `ns_max_p` covers both
	// `-ns_max_p` (upstream POSIX-ish form) and the long `--ns_max_p`.
	// They alias the same destination as `--max-ns-percent` etc.
	fs.Float64Var(&maxNsP, "ns_max_p", 0, "")
	fs.IntVar(&maxNsN, "ns_max_n", 0, "")

	// --noniupac strict-IUPAC filter (upstream prinseq-lite.pl:3478,
	// `[^ACGTN]/o`).
	var noniupac bool
	fs.BoolVar(&noniupac, "noniupac", false, "")

	// --seq_id <prefix> renames headers to "<prefix><counter>" on
	// records that pass all filters (prinseq-lite.pl:3640-3648).
	var seqID string
	fs.StringVar(&seqID, "seq_id", "", "")

	// --seq_id_mappings <file> writes a TSV of original-to-new ids
	// (prinseq-lite.pl:3646; requires -seq_id).
	var seqIDMappings string
	fs.StringVar(&seqIDMappings, "seq_id_mappings", "", "")

	// --out_format INT (1=FASTA, 2=FASTA+QUAL, 3=FASTQ,
	// 4=FASTQ+FASTA, 5=FASTQ+FASTA+QUAL); see prinseq-lite.pl:242-247
	// and the mode 1-5 branches at lines 769-789, 1302-1348, 3711-3714.
	var outFormat int
	fs.IntVar(&outFormat, "out_format", 0, "")

	// --phred64 is an input-encoding toggle, equivalent to passing
	// --qual-type illumina (prinseq-lite.pl:230-232, 760-764). We
	// register it as a bool here for upstream compatibility; it is
	// translated into opts.QualType below.
	var phred64 bool
	fs.BoolVar(&phred64, "phred64", false, "")

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

	// Quality encoding option
	var qualType string
	cliflag.StringVar(fs, &qualType, "t", "qual-type", "sanger", "Quality type: sanger (Phred+33) or illumina (Phred+64)")

	// Complexity filtering options
	var lcMethod string
	var lcThreshold float64
	cliflag.StringVar(fs, &lcMethod, "", "lc-method", "", "Low complexity filter method: dust or entropy")
	cliflag.Float64Var(fs, &lcThreshold, "", "lc-threshold", 0, "Low complexity threshold (default: 7 for dust, 70 for entropy)")

	fs.Usage = func() {
		fmt.Print(`Usage: prinseq filter [options]

Filter sequences based on quality criteria.

Input/Output Options:
  -i, --input FILE          Primary input file (use '-' for stdin)
  --input2 FILE             Paired-end input file 2
  -o, --output FILE         Output file (default: stdout)
  --output2 FILE            Output file for paired-end file 2
  --out-bad FILE            Output file for rejected sequences
  --fasta                   Input is FASTA format
  --fastq                   Input is FASTQ format

Filter Options:
  -l, --min-length INT      Minimum sequence length
  -L, --max-length INT      Maximum sequence length
  -g, --min-gc FLOAT        Minimum GC content (0-100)
  -G, --max-gc FLOAT        Maximum GC content (0-100)
  -q, --min-quality FLOAT   Minimum mean quality score
  -Q, --max-quality FLOAT   Maximum mean quality score
  -n, --max-ns INT          Maximum number of Ns (alias: -ns_max_n)
  -N, --max-ns-percent FLOAT Maximum percentage of Ns (0-100) (alias: -ns_max_p)
  --noniupac                Filter sequences with bases outside ACGTN

Identifier / Output-format Options:
  --seq_id PREFIX           Rename passing records to "<PREFIX><N>"
  --seq_id_mappings FILE    Write "<orig>\t<new>" TSV (requires --seq_id)
  --out_format INT          1=FASTA, 2=FASTA+QUAL, 3=FASTQ,
                            4=FASTQ+FASTA, 5=FASTQ+FASTA+QUAL
  --phred64                 Input FASTQ uses Phred+64 (alias for
                            -t illumina)

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

Quality Encoding:
  -t, --qual-type TYPE      Quality encoding: sanger (Phred+33, default) or illumina (Phred+64)

Complexity Filtering:
  --lc-method METHOD        Low complexity method: dust or entropy
  --lc-threshold FLOAT      Low complexity threshold (default: 7 for dust, 70 for entropy)

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
  
  # Use Phred+64 encoding (Illumina 1.3-1.7)
  prinseq filter -i reads.fastq -o filtered.fastq -t illumina -l 100
  
  # Filter low complexity sequences
  prinseq filter -i seqs.fasta -o filtered.fasta --lc-method dust --lc-threshold 7
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

	// Set default thresholds for complexity filtering if method is specified
	if lcMethod != "" && lcThreshold == 0 {
		if lcMethod == "dust" {
			lcThreshold = 7
		} else if lcMethod == "entropy" {
			lcThreshold = 70
		}
	}

	// --phred64 is an alias for --qual-type illumina (matches upstream
	// prinseq-lite.pl:230-232). We honour --qual-type as authoritative
	// if both are set, but otherwise flip the encoding here.
	if phred64 && qualType == "sanger" {
		qualType = "illumina"
	}

	// Validate --seq_id_mappings / --seq_id coupling. Upstream
	// (prinseq-lite.pl:945-946) prints "option -seq_id_mappings
	// requires option -seq_id" and exits.
	if seqIDMappings != "" && seqID == "" {
		fmt.Fprintln(os.Stderr, "Error: option --seq_id_mappings requires option --seq_id")
		os.Exit(1)
	}

	// Validate --out_format range.
	if outFormat < 0 || outFormat > 5 {
		fmt.Fprintln(os.Stderr, "Error: --out_format must be an integer between 1 and 5")
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
		QualType:      qualType,
		LcMethod:      lcMethod,
		LcThreshold:   lcThreshold,
		NonIUPAC:      noniupac,
		SeqID:         seqID,
		OutFormat:     outFormat,
	}

	// Open the seq_id_mappings TSV writer (when requested). Upstream
	// truncates the file on each run (open mode ">"), so we use
	// os.Create here.
	if seqIDMappings != "" {
		mapW, err := os.Create(seqIDMappings)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating seq_id_mappings file: %v\n", err)
			os.Exit(1)
		}
		defer mapW.Close()
		opts.SeqIDMap = mapW
	}

	// For multi-stream --out_format modes (2, 4, 5) we need additional
	// output files. Upstream derives them from `-out_good <prefix>`
	// using `.fasta` and `.qual` suffixes; in this Go port we use the
	// `--output` value as that prefix when set. If `--output` is unset
	// (i.e. primary stream goes to stdout), require an explicit
	// prefix via an extra positional path or fall back to "out".
	if outFormat == 2 || outFormat == 4 || outFormat == 5 {
		prefix := output1
		if prefix == "" {
			// Upstream raises an error here (line 801-802): you cannot
			// write multi-file output to STDOUT.
			fmt.Fprintln(os.Stderr, "Error: --out_format 2/4/5 require --output to provide a filename prefix (cannot stream multiple files to stdout)")
			os.Exit(1)
		}
		if outFormat == 4 || outFormat == 5 {
			f, err := os.Create(prefix + ".fasta")
			if err != nil {
				fmt.Fprintf(os.Stderr, "Error creating .fasta output: %v\n", err)
				os.Exit(1)
			}
			defer f.Close()
			opts.FastaOut = f
		}
		if outFormat == 2 || outFormat == 5 {
			f, err := os.Create(prefix + ".qual")
			if err != nil {
				fmt.Fprintf(os.Stderr, "Error creating .qual output: %v\n", err)
				os.Exit(1)
			}
			defer f.Close()
			opts.QualOut = f
		}
	}

	// Open bad output file if specified
	if outBad != "" {
		badWriter, err := os.Create(outBad)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating bad output file: %v\n", err)
			os.Exit(1)
		}
		defer badWriter.Close()
		opts.OutBad = badWriter
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

func runGraph(args []string) {
	fs := flag.NewFlagSet("graph", flag.ExitOnError)

	var fastq, fasta, graphType, output string
	var svg bool
	cliflag.StringVar(fs, &fastq, "", "fastq", "", "Input FASTQ file")
	cliflag.StringVar(fs, &fasta, "", "fasta", "", "Input FASTA file")
	cliflag.StringVar(fs, &graphType, "t", "type", "length", "Graph type: length, gc, quality, dinucleotides, positional_quality")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout)")
	cliflag.BoolVar(fs, &svg, "", "svg", false, "Generate SVG output")

	fs.Usage = func() {
		fmt.Print(`Usage: prinseq graph [options]

Generate quality graphs from sequence statistics.

Options:
  --fasta FILE       Input FASTA file
  --fastq FILE       Input FASTQ file
  -t, --type TYPE    Graph type (length, gc, quality, dinucleotides, positional_quality)
  -o, --output FILE  Output file (default: stdout)
  --svg              Generate SVG output (default: ASCII)
`)
	}

	if err := fs.Parse(args); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

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

	reader, err := os.Open(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening file: %v\n", err)
		os.Exit(1)
	}
	defer reader.Close()

	stats, err := prinseq.CalculateEnhancedStats(reader, isFastq)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error calculating stats: %v\n", err)
		os.Exit(1)
	}

	var writer io.WriteCloser = os.Stdout
	if output != "" {
		writer, err = os.Create(output)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
		defer writer.Close()
	}

	if svg {
		err = prinseq.GenerateSVG(stats, writer)
	} else {
		err = prinseq.GenerateGraph(stats, prinseq.GraphType(graphType), writer)
	}

	if err != nil {
		fmt.Fprintf(os.Stderr, "Error generating graph: %v\n", err)
		os.Exit(1)
	}
}

// runGraphData implements the upstream `--graph_data` flag from
// `prinseq-lite.pl`. It collects the full graph-data stat tables
// for one input file (FASTA or FASTQ) and writes the upstream
// `.gd` JSON to either the path given on the flag or the
// upstream-default `<input>__.gd`.
//
// Compared to upstream, key order in every map is sorted
// lexicographically (intentional deviation; see graphdata.go and
// docs/UPSTREAM_BUGS.md). The two leading #-comment lines are
// emitted with the real timestamp and argv when --no-header is
// not specified.
func runGraphData(args []string) {
	fs := flag.NewFlagSet("graph_data", flag.ExitOnError)

	var fastq, fasta, graphData, graphStats string
	var phred64, qualNoScale, noHeader bool
	cliflag.StringVar(fs, &fastq, "", "fastq", "", "Input FASTQ file (use '-' for stdin)")
	cliflag.StringVar(fs, &fasta, "", "fasta", "", "Input FASTA file (use '-' for stdin)")
	cliflag.StringVar(fs, &graphData, "", "graph_data", "", "Output .gd path (default: <input>__.gd)")
	cliflag.StringVar(fs, &graphStats, "", "graph_stats", "", "Comma-separated stat selector: ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn")
	cliflag.BoolVar(fs, &phred64, "", "phred64", false, "Input FASTQ uses Phred+64")
	cliflag.BoolVar(fs, &qualNoScale, "", "qual_noscale", false, "Disable per-position quality relative-bin table")
	cliflag.BoolVar(fs, &noHeader, "", "no-header", false, "Omit the two leading #-comment lines (test mode)")

	fs.Usage = func() {
		fmt.Print(`Usage: prinseq graph_data --fastq FILE [options]

Emit the upstream prinseq-lite .gd graph-data JSON for a single
input file.

Options:
  --fastq FILE          Input FASTQ file
  --fasta FILE          Input FASTA file
  --graph_data FILE     Output .gd path (default: <input>__.gd)
  --graph_stats CODES   Comma-separated stat selector (ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn)
  --phred64             Input FASTQ uses Phred+64 encoding
  --qual_noscale        Disable the relative (100-bin) quality table
  --no-header           Omit the leading #-comment lines (testing only)
`)
	}

	if err := fs.Parse(args); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	var inputFile string
	var isFastq bool
	switch {
	case fastq != "":
		inputFile = fastq
		isFastq = true
	case fasta != "":
		inputFile = fasta
		isFastq = false
	default:
		fmt.Fprintln(os.Stderr, "Error: Either --fastq or --fasta must be specified")
		fs.Usage()
		os.Exit(1)
	}

	opts := prinseq.DefaultGraphDataOptions()
	if graphStats != "" {
		if err := prinseq.ParseGraphStatsCSV(graphStats, &opts); err != nil {
			fmt.Fprintf(os.Stderr, "Error parsing --graph_stats: %v\n", err)
			os.Exit(1)
		}
	}
	opts.Phred64 = phred64
	opts.QualNoscale = qualNoScale
	opts.Filename1 = filenameHex(inputFile)

	reader, err := openInput(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening file: %v\n", err)
		os.Exit(1)
	}
	defer reader.Close()

	g, err := prinseq.CollectGraphData(reader, isFastq, opts)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error collecting graph data: %v\n", err)
		os.Exit(1)
	}

	outPath := prinseq.ResolveGraphDataPath(graphData, inputFile)
	out, err := os.Create(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	var header prinseq.GDHeader
	if !noHeader {
		header = prinseq.GDHeader{
			Version:   "0.20.4",
			Timestamp: nowUpstreamFmt(),
			Command:   strings.Join(append([]string{"prinseq", "graph_data"}, args...), " "),
		}
	}
	if err := g.EmitGD(out, header); err != nil {
		fmt.Fprintf(os.Stderr, "Error writing graph data: %v\n", err)
		os.Exit(1)
	}
}

// filenameHex hex-encodes the basename of a path, matching upstream's
// convertStringToInt (prinseq-lite.pl:4851-4855). When the path is
// empty or "-" (stdin), we return "stdin" so downstream renderers can
// distinguish piped input from a missing filename.
func filenameHex(path string) string {
	if path == "" || path == "-" {
		return "stdin"
	}
	// Take basename.
	base := path
	if i := strings.LastIndex(path, "/"); i >= 0 {
		base = path[i+1:]
	}
	var b strings.Builder
	for i := 0; i < len(base); i++ {
		fmt.Fprintf(&b, "%x", base[i])
	}
	return b.String()
}

// nowUpstreamFmt returns the current local time in the upstream
// graph-data header timestamp format: "MM/DD/YYYY HH:MM:SS"
// (prinseq-lite.pl:2277).
func nowUpstreamFmt() string {
	t := nowFunc()
	return fmt.Sprintf("%02d/%02d/%04d %02d:%02d:%02d",
		t.Month(), t.Day(), t.Year(), t.Hour(), t.Minute(), t.Second())
}

func runReport(args []string) {
	fs := flag.NewFlagSet("report", flag.ExitOnError)

	var fastq, fasta, output string
	cliflag.StringVar(fs, &fastq, "", "fastq", "", "Input FASTQ file")
	cliflag.StringVar(fs, &fasta, "", "fasta", "", "Input FASTA file")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output HTML file (default: stdout)")

	fs.Usage = func() {
		fmt.Print(`Usage: prinseq report [options]

Generate an HTML quality report with embedded graphs.

Options:
  --fasta FILE       Input FASTA file
  --fastq FILE       Input FASTQ file
  -o, --output FILE  Output HTML file (default: stdout)
`)
	}

	if err := fs.Parse(args); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

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

	reader, err := os.Open(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening file: %v\n", err)
		os.Exit(1)
	}
	defer reader.Close()

	stats, err := prinseq.CalculateEnhancedStats(reader, isFastq)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error calculating stats: %v\n", err)
		os.Exit(1)
	}

	var writer io.WriteCloser = os.Stdout
	if output != "" {
		writer, err = os.Create(output)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
		defer writer.Close()
	}

	if err := prinseq.GenerateHTMLReport(stats, writer); err != nil {
		fmt.Fprintf(os.Stderr, "Error generating report: %v\n", err)
		os.Exit(1)
	}
}

func runBenchmark(args []string) {
	fs := flag.NewFlagSet("benchmark", flag.ExitOnError)

	var fastq, fasta string
	var jsonOutput bool
	cliflag.StringVar(fs, &fastq, "", "fastq", "", "Input FASTQ file")
	cliflag.StringVar(fs, &fasta, "", "fasta", "", "Input FASTA file")
	cliflag.BoolVar(fs, &jsonOutput, "j", "json", false, "Output results in JSON format")

	fs.Usage = func() {
		fmt.Print(`Usage: prinseq benchmark [options]

Run performance benchmarks on sequence processing operations.

Options:
  --fasta FILE       Input FASTA file
  --fastq FILE       Input FASTQ file
  -j, --json         Output results in JSON format
`)
	}

	if err := fs.Parse(args); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

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

	reader, err := os.Open(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening file: %v\n", err)
		os.Exit(1)
	}
	defer reader.Close()

	results, err := prinseq.RunBenchmarkSuite(reader, isFastq)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error running benchmark: %v\n", err)
		os.Exit(1)
	}

	if jsonOutput {
		encoder := json.NewEncoder(os.Stdout)
		encoder.SetIndent("", "  ")
		if err := encoder.Encode(results); err != nil {
			fmt.Fprintf(os.Stderr, "Error encoding JSON: %v\n", err)
			os.Exit(1)
		}
	} else {
		fmt.Print(prinseq.FormatBenchmarkResults(results))
	}
}

func runAPI(args []string) {
	fs := flag.NewFlagSet("api", flag.ExitOnError)

	var addr string
	cliflag.StringVar(fs, &addr, "a", "addr", ":8080", "Server address")

	fs.Usage = func() {
		fmt.Print(`Usage: prinseq api [options]

Start a REST API server for PRINSEQ operations.

Options:
  -a, --addr ADDR    Server address (default: :8080)

Examples:
  prinseq api --addr :8080
  prinseq api --addr localhost:9000
`)
	}

	if err := fs.Parse(args); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	server := prinseq.NewAPIServer(addr)
	if err := server.Start(); err != nil {
		fmt.Fprintf(os.Stderr, "Error starting server: %v\n", err)
		os.Exit(1)
	}
}

func runBatch(args []string) {
	fs := flag.NewFlagSet("batch", flag.ExitOnError)

	var outputDir string
	var workers int
	var generateReport bool
	var isFastq bool
	var minLen int
	var minGC, maxGC float64

	cliflag.StringVar(fs, &outputDir, "o", "output", "", "Output directory")
	cliflag.IntVar(fs, &workers, "w", "workers", 4, "Number of parallel workers")
	cliflag.BoolVar(fs, &generateReport, "r", "report", false, "Generate HTML reports")
	cliflag.BoolVar(fs, &isFastq, "", "fastq", false, "Input files are FASTQ format")
	cliflag.IntVar(fs, &minLen, "l", "min-length", 0, "Minimum sequence length for filtering")
	cliflag.Float64Var(fs, &minGC, "g", "min-gc", 0, "Minimum GC content")
	cliflag.Float64Var(fs, &maxGC, "G", "max-gc", 0, "Maximum GC content")

	fs.Usage = func() {
		fmt.Print(`Usage: prinseq batch [options] <input_files...>

Process multiple files in parallel.

Options:
  -o, --output DIR      Output directory
  -w, --workers N       Number of parallel workers (default: 4)
  -r, --report          Generate HTML reports
  --fastq               Input files are FASTQ format
  -l, --min-length INT  Minimum sequence length for filtering
  -g, --min-gc FLOAT    Minimum GC content
  -G, --max-gc FLOAT    Maximum GC content

Examples:
  prinseq batch --fastq -o output -w 8 *.fastq
  prinseq batch --fastq -o output -r *.fastq
`)
	}

	if err := fs.Parse(args); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	inputFiles := fs.Args()
	if len(inputFiles) == 0 {
		fmt.Fprintln(os.Stderr, "Error: No input files specified")
		fs.Usage()
		os.Exit(1)
	}

	config := prinseq.BatchProcessConfig{
		InputFiles:     inputFiles,
		OutputDir:      outputDir,
		IsFastq:        isFastq,
		Workers:        workers,
		GenerateReport: generateReport,
		FilterOpts: prinseq.FilterOptions{
			MinLen: minLen,
			MinGC:  minGC,
			MaxGC:  maxGC,
		},
	}

	fmt.Printf("Processing %d files with %d workers...\n", len(inputFiles), workers)

	results, err := prinseq.BatchProcess(config)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error in batch processing: %v\n", err)
		os.Exit(1)
	}

	// Print summary
	fmt.Printf("\nProcessing complete!\n")
	fmt.Printf("Successfully processed: %d files\n", len(results))

	for _, result := range results {
		if result.Error != nil {
			fmt.Printf("  ✗ %s: %v\n", result.Filename, result.Error)
		} else {
			fmt.Printf("  ✓ %s: %d reads, %.2f avg length\n",
				result.Filename, result.Stats.NumReads, result.Stats.AvgLength)
		}
	}
}
