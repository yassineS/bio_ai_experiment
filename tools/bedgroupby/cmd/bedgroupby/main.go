// bedgroupby groups TSV/BED records by one or more columns and applies
// aggregation operations to each group, mirroring `bedtools groupby`.
package main

import (
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedgroupby/pkg/bedgroupby"
)

const usage = `bedgroupby - Group TSV/BED records and aggregate columns

Usage:
  bedgroupby [options] -i <input> -c <cols> [-o <ops>]

Description:
  Reads a TSV/BED-like file and groups consecutive records sharing identical
  values in the grouping columns (default: columns 1,2,3). For each group the
  requested aggregation operations are applied to the requested columns, and
  one row per group is emitted.

Options:
  -i, --input FILE       Input file (default: stdin; use '-' for stdin)
      --output FILE      Output file (default: stdout)
  -g, --grp LIST         Grouping columns (1-based, ranges like "2-4" allowed;
                         default: 1,2,3)
  -c, --opCols LIST      Comma-separated columns to aggregate (required)
  -o, --ops LIST         Operations: one per column, or one applied to all
                         (default: sum). Supported: sum, min, max, absmin,
                         absmax, mean, median, stdev, sstdev, count,
                         count_distinct, distinct, collapse, first, last,
                         mode, antimode.
      --full             Emit every column of the first record in each group
                         before the aggregated values.
      --ignorecase       Case-insensitive grouping.
      --inheader         Treat the first non-empty line as a header even if
                         it does not start with one of the recognised marker
                         prefixes (#, track, browser, @).
      --outheader        Print the detected/synthesised header on output.
      --header           Combination of --inheader and --outheader.
  -h, --help             Show this help message.
  -v, --version          Show version and exit.

Examples:
  # Default: group by chrom/start/end, sum column 5
  bedgroupby -i values.bed -c 5

  # Group by columns 2-4 of a non-BED file, mean of column 6
  bedgroupby -g 2-4 -i values.bed -c 6 -o mean

  # Multiple ops, one per column
  bedgroupby -c 2,3 -o distinct,min -i bug569.txt

Notes:
  - Input must already be sorted on the grouping columns (matches upstream).
  - '-' for input means standard input.
  - Gzip and BGZF inputs are decompressed transparently.
`

const version = "bedgroupby 0.1.0"

func main() {
	fs := flag.CommandLine

	var inputFile, outputFile string
	cliflag.StringVar(fs, &inputFile, "i", "input", "", "Input file (default: stdin)")
	cliflag.StringVar(fs, &outputFile, "", "output", "", "Output file (default: stdout)")

	var groupSpec string
	cliflag.StringVar(fs, &groupSpec, "g", "grp", "", "Grouping columns")

	var aggCols string
	cliflag.StringVar(fs, &aggCols, "c", "opCols", "", "Aggregation columns")

	var ops string
	cliflag.StringVar(fs, &ops, "o", "ops", "", "Aggregation operations")

	var full bool
	cliflag.BoolVar(fs, &full, "", "full", false, "Emit full first-record columns before aggregates")

	var ignoreCase bool
	cliflag.BoolVar(fs, &ignoreCase, "", "ignorecase", false, "Case-insensitive grouping")

	var inHeader bool
	cliflag.BoolVar(fs, &inHeader, "", "inheader", false, "Treat first line as header")

	var outHeader bool
	cliflag.BoolVar(fs, &outHeader, "", "outheader", false, "Print header on output")

	var header bool
	cliflag.BoolVar(fs, &header, "", "header", false, "Equivalent to --inheader --outheader")

	var help bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help message")

	var showVersion bool
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version")

	flag.Parse()

	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if showVersion {
		fmt.Println(version)
		os.Exit(0)
	}

	groupCols, err := bedgroupby.ParseGroupSpec(groupSpec)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	if strings.TrimSpace(aggCols) == "" {
		fmt.Fprintln(os.Stderr, "Error: -c/--opCols is required")
		os.Exit(1)
	}
	var aggColsList []int
	for _, p := range strings.Split(aggCols, ",") {
		p = strings.TrimSpace(p)
		if p == "" {
			continue
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error: invalid column %q in -c/--opCols: %v\n", p, err)
			os.Exit(1)
		}
		aggColsList = append(aggColsList, n)
	}

	var opList []string
	if strings.TrimSpace(ops) != "" {
		for _, p := range strings.Split(ops, ",") {
			opList = append(opList, strings.TrimSpace(p))
		}
	}

	input := inputFile
	if input == "" && flag.NArg() > 0 {
		input = flag.Arg(0)
	}
	in, err := iohelper.OpenReader(input)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input: %v\n", err)
		os.Exit(1)
	}
	defer in.Close()

	out, err := iohelper.OpenWriter(outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if _, err := bedgroupby.Group(in, out, bedgroupby.Options{
		GroupCols:  groupCols,
		AggCols:    aggColsList,
		Ops:        opList,
		Full:       full,
		IgnoreCase: ignoreCase,
		InHeader:   inHeader,
		OutHeader:  outHeader,
		Header:     header,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}
