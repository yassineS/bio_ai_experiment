// bedmap applies column aggregation ops to A's intervals using B as the
// value source. Mirrors `bedtools map`.
package main

import (
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedmap/pkg/bedmap"
)

const usage = `bedmap - Apply column aggregation ops to A's intervals using B

Usage:
  bedmap [options] -a <A.bed> -b <B.bed>

Description:
  For each interval in A, find overlapping B intervals, take values from a
  single B column (default: column 5), and apply an aggregation operation
  (default: sum). When no B intervals overlap A, emit the configured null
  placeholder (default: ".") for each requested column.

Options:
  -a, --input-a FILE    A intervals (required; '-' for stdin)
  -b, --input-b FILE    B intervals (required)
      --output FILE     Output file (default: stdout)
  -c, --columns LIST    Comma-separated 1-based B column indices (default: 5)
  -o, --operations LIST One op per column or a single op for all (default: sum).
                        Supported: sum, min, max, mean, median, count,
                        count_distinct, distinct, collapse, first, last, mode,
                        antimode
      --null STRING     Placeholder for no overlap (default: ".")
      --delim STRING    Separator for collapse / distinct (default: ",")
  -s, --strand          Only consider B on the same strand as A
  -S, --opposite        Only consider B on the opposite strand
  -f NUM,--fraction-a   Min fraction of A that must overlap a B record
  -F NUM,--fraction-b   Min fraction of B that must overlap A
  -r,    --reciprocal   Require both -f and -F (default already AND)
  -h,    --help         Show this help
  -v,    --version      Show version

Examples:
  # Sum values in B's column 5 (default) for each A:
  bedmap -a A.bed -b B.bed

  # Collapse names from B's column 4:
  bedmap -a A.bed -b B.bed -c 4 -o collapse

  # Mean of B's column 5, missing → "NA":
  bedmap -a A.bed -b B.bed -o mean --null NA

  # Multiple ops on one column:
  bedmap -a A.bed -b B.bed -c 5 -o min,max
`

const version = "bedmap 0.1.0"

func main() {
	fs := flag.CommandLine

	var inputA, inputB, output string
	cliflag.StringVar(fs, &inputA, "a", "input-a", "", "A intervals")
	cliflag.StringVar(fs, &inputB, "b", "input-b", "", "B intervals")
	cliflag.StringVar(fs, &output, "", "output", "", "Output (default stdout)")

	var cols, ops, null, delim string
	cliflag.StringVar(fs, &cols, "c", "columns", "", "B columns to aggregate (default 5)")
	cliflag.StringVar(fs, &ops, "o", "operations", "", "Aggregation ops (default sum)")
	cliflag.StringVar(fs, &null, "", "null", ".", "Placeholder for no overlap")
	cliflag.StringVar(fs, &delim, "", "delim", ",", "Separator for collapse / distinct")

	var sameStrand, oppStrand bool
	cliflag.BoolVar(fs, &sameStrand, "s", "strand", false, "Same strand only")
	cliflag.BoolVar(fs, &oppStrand, "S", "opposite", false, "Opposite strand only")

	var fracA, fracB float64
	cliflag.Float64Var(fs, &fracA, "f", "fraction-a", 0, "Min fraction of A")
	cliflag.Float64Var(fs, &fracB, "F", "fraction-b", 0, "Min fraction of B")

	var reciprocal bool
	cliflag.BoolVar(fs, &reciprocal, "r", "reciprocal", false, "Require both -f and -F")

	var help, showVersion bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help")
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

	if inputA == "" || inputB == "" {
		fmt.Fprintln(os.Stderr, "Error: -a and -b are both required")
		fmt.Fprintln(os.Stderr, "Use -h for help.")
		os.Exit(1)
	}

	colList, err := parseCols(cols)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
	opList := parseOps(ops)

	rA, err := iohelper.OpenReader(inputA)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening A: %v\n", err)
		os.Exit(1)
	}
	defer rA.Close()
	rB, err := iohelper.OpenReader(inputB)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening B: %v\n", err)
		os.Exit(1)
	}
	defer rB.Close()
	w, err := iohelper.OpenWriter(output)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		os.Exit(1)
	}
	defer w.Close()

	if _, err := bedmap.Map(rA, rB, w, bedmap.Options{
		Columns:        colList,
		Ops:            opList,
		Null:           null,
		Delim:          delim,
		SameStrand:     sameStrand,
		OppositeStrand: oppStrand,
		FractionA:      fracA,
		FractionB:      fracB,
		Reciprocal:     reciprocal,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

// parseCols parses a comma-separated 1-based column-number list.
func parseCols(s string) ([]int, error) {
	s = strings.TrimSpace(s)
	if s == "" {
		return nil, nil
	}
	var out []int
	for _, p := range strings.Split(s, ",") {
		p = strings.TrimSpace(p)
		if p == "" {
			continue
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			return nil, fmt.Errorf("invalid column %q in -c: %v", p, err)
		}
		out = append(out, n)
	}
	return out, nil
}

// parseOps splits a comma-separated op list.
func parseOps(s string) []string {
	s = strings.TrimSpace(s)
	if s == "" {
		return nil
	}
	parts := strings.Split(s, ",")
	for i := range parts {
		parts[i] = strings.TrimSpace(parts[i])
	}
	return parts
}
