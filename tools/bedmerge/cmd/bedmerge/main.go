// bedmerge merges overlapping or adjacent BED/GFF/VCF/BAM intervals, a drop-in
// re-implementation of `bedtools merge`.
package main

import (
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedmerge/pkg/bedmerge"
)

const version = "1.0.0"

const usage = `bedmerge - Merge overlapping or adjacent BED/GFF/VCF/BAM intervals

Usage:
  bedmerge [options] -i <bed/gff/vcf/bam>

Description:
  Reads intervals (BED, GFF, VCF, or BAM; auto-detected), sorts them by
  chromosome and position, and merges overlapping or adjacent intervals.
  Intervals are merged if they overlap or are within the -d distance.

Options:
  -i, --input FILE       Input file (BED/GFF/VCF/BAM; default: stdin)
  -s, --strand           Force strandedness: only merge same-strand features
  -S <+|->               Force merge of one specific strand only
  -d, --distance INT     Maximum distance between features to merge (default: 0;
                         negative values require that many bp of overlap)
  -c, --columns LIST     Comma-separated 1-based input columns to aggregate
  -o, --operations LIST  Comma-separated operations, one per -c column or a
                         single op applied to all. Supported: sum, min, max,
                         absmin, absmax, mean, median, mode, antimode, stdev,
                         sstdev, collapse, distinct, distinct_sort_num,
                         distinct_sort_num_desc, distinct_only, count,
                         count_distinct, first, last (default: sum)
      --delim CHAR       Delimiter for collapse/distinct list values (default ",")
      --prec INT         Decimal precision for output (default: 10)
      --header           Print the header from the input file before results
      --bed              If using BAM input, write output as BED
      --nobuf            Disable buffered output
      --iobuf SIZE       Input buffer size (suffixes K/M/G allowed)
      --output FILE      Output file (default: stdout)
      --stats            Print merge statistics to stderr
      --count            Output count of merged intervals as name field
  -g, --bedgraph         Treat input column 5 as a bedGraph score
  -h, --help             Show this help message
  -v, --version          Show version information and exit

Examples:
  bedmerge -i input.bed > merged.bed
  bedmerge -i input.bed -d 100 > merged.bed
  bedmerge -i input.bed -s > merged.bed
  bedmerge -i input.bed -c 4,5 -o distinct,sum > merged.bed
  bedmerge -i input.bam -c 1 -o collapse > merged.bed
  cat input.bed | bedmerge > merged.bed

Notes:
  - Coordinates are 0-based, half-open [start, end)
  - Input does not need to be pre-sorted; bedmerge sorts internally
  - Adjacent (book-ended) intervals are merged by default
`

func main() {
	os.Exit(run(os.Args[1:]))
}

// run parses args, performs the merge, and returns the process exit code. It is
// separated from main so the behaviour can be unit-tested.
func run(args []string) int {
	// -iobuf must be validated even when given without a value (the flag parser
	// would otherwise consume the next token), so pre-scan it to emit the exact
	// upstream messages (merge.t37-t40).
	if code, exit := checkIobuf(args); exit {
		return code
	}

	fs := flag.NewFlagSet("bedmerge", flag.ContinueOnError)
	fs.SetOutput(os.Stderr)
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	var distance int
	cliflag.IntVar(fs, &distance, "d", "distance", 0, "Maximum distance between intervals to merge")

	var strandSpec bool
	cliflag.BoolVar(fs, &strandSpec, "s", "strand", false, "Merge only intervals on the same strand")

	var strandFilter string
	fs.StringVar(&strandFilter, "S", "", "Merge only intervals on the given strand (+ or -)")

	var inputFile string
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input file (default: stdin)")

	var outputFile string
	cliflag.StringVar(fs, &outputFile, "", "output", "", "Output file (default: stdout)")

	var showStats bool
	cliflag.BoolVar(fs, &showStats, "", "stats", false, "Print merge statistics to stderr")

	var showCount bool
	cliflag.BoolVar(fs, &showCount, "", "count", false, "Output count of merged intervals as name field")

	var bedGraph bool
	cliflag.BoolVar(fs, &bedGraph, "g", "bedgraph", false, "Treat input column 5 as a bedGraph score")

	var columns string
	cliflag.StringVar(fs, &columns, "c", "columns", "", "Comma-separated 1-based input columns to aggregate")

	var operations string
	cliflag.StringVar(fs, &operations, "o", "operations", "", "Comma-separated operations, one per -c column or one applied to all")

	var delim string
	fs.StringVar(&delim, "delim", ",", "Delimiter for collapse/distinct/freq list joins")

	var precision int
	cliflag.IntVar(fs, &precision, "", "prec", bedmerge.DefaultPrecision, "Decimal precision for output")

	var header bool
	cliflag.BoolVar(fs, &header, "", "header", false, "Print the input file's header before results")

	var bedOut bool
	cliflag.BoolVar(fs, &bedOut, "", "bed", false, "Write BAM output as BED (always BED here)")
	_ = bedOut // Output is always BED; the flag is accepted for compatibility.

	var nobuf bool
	cliflag.BoolVar(fs, &nobuf, "", "nobuf", false, "Disable buffered output")
	_ = nobuf // Accepted for compatibility; output correctness is unaffected.

	var iobuf string
	fs.StringVar(&iobuf, "iobuf", "", "Input buffer size (suffixes K/M/G allowed)")

	// Deprecated upstream flags: -n and -nms are no longer supported and must
	// produce the documented error (merge.t2 / merge.t4).
	var deprecatedN, deprecatedNMS bool
	fs.BoolVar(&deprecatedN, "n", false, "(deprecated)")
	fs.BoolVar(&deprecatedNMS, "nms", false, "(deprecated)")

	var streaming bool
	cliflag.BoolVar(fs, &streaming, "", "streaming", false, "Use streaming mode for very large files")

	var help bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help message")

	var showVersion bool
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version information")

	if err := fs.Parse(args); err != nil {
		return 2 // flag already printed the error and usage (conventional exit 2).
	}

	switch {
	case help:
		fmt.Fprint(os.Stderr, usage)
		return 0
	case showVersion:
		fmt.Printf("bedmerge version %s\n", version)
		return 0
	case deprecatedN:
		fmt.Fprintln(os.Stderr, "***** ERROR: -n option is deprecated. Please see the documentation for the -c and -o column operation options. *****")
		return 1
	case deprecatedNMS:
		fmt.Fprintln(os.Stderr, "***** ERROR: -nms option is deprecated. Please see the documentation for the -c and -o column operation options. *****")
		return 1
	}

	// -S must be exactly + or - (merge.t18).
	if strandFilter != "" && strandFilter != "+" && strandFilter != "-" {
		fmt.Fprintln(os.Stderr, "***** ERROR: -S option must be followed by + or -. *****")
		return 1
	}

	columnOps, err := bedmerge.ParseColumnOps(columns, operations)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		return 1
	}
	if columnOps != nil {
		columnOps.Delim = delim
	}

	// Determine input file (positional arg also accepted).
	input := inputFile
	if input == "" && fs.NArg() > 0 {
		input = fs.Arg(0)
	}

	inputReader, err := iohelper.OpenReader(input)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input: %v\n", err)
		return 1
	}
	defer inputReader.Close()

	outputWriter, err := iohelper.OpenWriter(outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		return 1
	}
	defer outputWriter.Close()

	if header {
		if err := printHeader(input, outputWriter); err != nil {
			fmt.Fprintf(os.Stderr, "Error reading header: %v\n", err)
			return 1
		}
	}

	opts := bedmerge.MergeOptions{
		MaxDistance:  distance,
		StrandSpec:   strandSpec,
		StrandFilter: strandFilter,
		Streaming:    streaming,
		ColumnOps:    columnOps,
		Precision:    precision,
		Warn:         os.Stderr,
		OutputFields: bedmerge.OutputFields{
			Count:    showCount,
			BedGraph: bedGraph,
		},
	}

	if showStats {
		stats, err := bedmerge.MergeWithStats(inputReader, outputWriter, opts)
		if err != nil {
			reportMergeError(err, input)
			return 1
		}
		fmt.Fprintf(os.Stderr, "Input intervals:  %d\n", stats.InputIntervals)
		fmt.Fprintf(os.Stderr, "Output intervals: %d\n", stats.OutputIntervals)
		fmt.Fprintf(os.Stderr, "Merged: %d intervals\n", stats.MergedCount)
		return 0
	}

	if _, err := bedmerge.Merge(inputReader, outputWriter, opts); err != nil {
		reportMergeError(err, input)
		return 1
	}
	return 0
}

// reportMergeError prints the upstream-formatted message for the sentinel errors
// (stranded-VCF, bad -S argument) and a generic message otherwise.
func reportMergeError(err error, input string) {
	var bamErr *bedmerge.BAMColumnError
	switch {
	case errors.As(err, &bamErr):
		if bamErr.Flags {
			fmt.Fprintln(os.Stderr, "***** ERROR: Requested column 2 of a BAM file, which is the Flags field.")
		} else {
			name := input
			if name == "" || name == "-" {
				name = "stdin"
			}
			fmt.Fprintf(os.Stderr, "***** ERROR: Requested column %d, but database file %s only has fields 1 - 11.\n", bamErr.Column, name)
		}
	case errors.Is(err, bedmerge.ErrStrandedVCF):
		name := input
		if name == "" || name == "-" {
			name = "stdin"
		}
		fmt.Fprintf(os.Stderr, "***** ERROR: stranded merge not supported for VCF file %s. *****\n", name)
	case errors.Is(err, bedmerge.ErrBadStrandArg):
		fmt.Fprintln(os.Stderr, "***** ERROR: -S option must be followed by + or -. *****")
	default:
		fmt.Fprintf(os.Stderr, "Error merging intervals: %v\n", err)
	}
}

// printHeader copies leading comment/track/browser lines from the input file to
// w, mirroring upstream `bedtools merge -header`. It opens the file separately
// (the merge reader is consumed independently) and stops at the first data line.
func printHeader(input string, w io.Writer) error {
	r, err := iohelper.OpenReader(input)
	if err != nil {
		return err
	}
	defer r.Close()
	return bedmerge.WriteHeader(r, w)
}

// checkIobuf validates a -iobuf argument exactly as upstream does
// (merge.t37-t40). It returns (exitCode, true) when bedmerge should stop with
// that code, or (0, false) to continue. A valid -iobuf is left for the flag
// parser to consume.
func checkIobuf(args []string) (int, bool) {
	for i := 0; i < len(args); i++ {
		a := args[i]
		if a != "-iobuf" && a != "--iobuf" {
			continue
		}
		if i+1 >= len(args) {
			fmt.Fprintln(os.Stderr, "***** ERROR: -iobuf option given, but size of input buffer not specified. *****")
			return 1, true
		}
		if code, ok := validateIobuf(args[i+1]); !ok {
			return code, true
		}
		return 0, false
	}
	return 0, false
}

// validateIobuf checks a single -iobuf value, emitting the matching upstream
// message and returning ok=false (with the exit code) when malformed.
func validateIobuf(val string) (int, bool) {
	suffix := byte(0)
	numStr := val
	if len(val) > 0 {
		last := val[len(val)-1]
		if last < '0' || last > '9' {
			suffix = last
			numStr = val[:len(val)-1]
		}
	}
	if suffix != 0 && suffix != 'K' && suffix != 'M' && suffix != 'G' {
		fmt.Fprintf(os.Stderr, "***** ERROR: Unrecognized memory buffer size suffix '%c' given. *****\n", suffix)
		return 1, false
	}
	n, err := strconv.Atoi(numStr)
	if err != nil {
		fmt.Fprintln(os.Stderr, "***** ERROR: argument passed to -iobuf is not numeric. *****")
		return 1, false
	}
	mult := 1
	switch suffix {
	case 'K':
		mult = 1 << 10
	case 'M':
		mult = 1 << 20
	case 'G':
		mult = 1 << 30
	}
	if n*mult < 8 {
		fmt.Fprintln(os.Stderr, "***** ERROR: specified buffer size is too small. *****")
		return 1, false
	}
	return 0, true
}
