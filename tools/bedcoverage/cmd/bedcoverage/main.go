// bedcoverage reports, for each interval in A, the number of features in B
// that overlap, the total bp covered, the length of A, and the fraction of A
// covered. Mirrors `bedtools coverage`.
package main

import (
	"errors"
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedcoverage/pkg/bedcoverage"
)

const usage = `bedcoverage - Per-A coverage statistics from B

Usage:
  bedcoverage [options] -a <A.bed> -b <B.bed>

Description:
  For each interval in A, report how much it is covered by B. The default
  output appends four columns to each A line:

      <count> <bp_covered> <length_A> <fraction>

  Various option flags switch to different summary shapes.

Options:
  -a, --input-a FILE     A intervals (required; '-' for stdin)
  -b, --input-b FILE     B intervals (required)
      --output FILE      Output file (default: stdout)
  -counts, --counts      Just append the overlap count to A
  -d,      --depth       Per-base depth: emit one line per base in A
  -hist,   --hist        Depth histogram per A plus an "all" footer
  -mean                  Append mean depth across A
  -median                Append median depth across A
  -min                   Append min depth across A
  -max                   Append max depth across A
  -sum                   Append sum of per-base depths across A
  -s,      --strand      Only count B features on the same strand as A
  -S,      --opposite    Only count B features on the opposite strand
  -f NUM, --fraction-a   Min fraction of A that must overlap a B feature
  -F NUM, --fraction-b   Min fraction of B that must overlap A
  -r,     --reciprocal   Require reciprocal -f overlap (A AND B)
  -e,     --either       Require -f OR -F (instead of the default AND)
  -h,     --help         Show this help
  -v,     --version      Show version
`

const version = "bedcoverage 0.1.0"

func main() {
	fs := flag.CommandLine

	var inputA, inputB, output string
	cliflag.StringVar(fs, &inputA, "a", "input-a", "", "A intervals")
	cliflag.StringVar(fs, &inputB, "b", "input-b", "", "B intervals")
	cliflag.StringVar(fs, &output, "", "output", "", "Output file (default stdout)")

	// Legacy upstream aliases for supplying BAM input. `bedtools coverage`
	// historically used -abam/-ibam for the query (A) and -bbam for the
	// database (B); the modern -a/-b already auto-detect BAM, so these alias
	// straight onto inputA/inputB.
	var abam, ibam, bbam string
	fs.StringVar(&abam, "abam", "", "BAM query alias for -a")
	fs.StringVar(&ibam, "ibam", "", "BAM query alias for -a")
	fs.StringVar(&bbam, "bbam", "", "BAM database alias for -b")

	// -sorted selects upstream's chromsweep (linear-merge) algorithm. Our
	// interval-tree path produces identical output, so we accept the flag for
	// drop-in compatibility and otherwise treat it as a no-op.
	var sorted bool
	cliflag.BoolVar(fs, &sorted, "", "sorted", false, "Inputs are position-sorted (accepted; no-op)")
	// -g supplies a genome file in upstream's -sorted mode (chrom ordering).
	// We don't need it, but accept it so sorted-mode command lines parse.
	var genome string
	cliflag.StringVar(fs, &genome, "g", "genome", "", "Genome file for -sorted (accepted; unused)")

	var counts, depth, hist bool
	cliflag.BoolVar(fs, &counts, "", "counts", false, "Just append the count")
	cliflag.BoolVar(fs, &depth, "d", "depth", false, "Per-base depth")
	cliflag.BoolVar(fs, &hist, "", "hist", false, "Depth histogram")

	var doMean, doMedian, doMin, doMax, doSum bool
	cliflag.BoolVar(fs, &doMean, "", "mean", false, "Mean depth")
	cliflag.BoolVar(fs, &doMedian, "", "median", false, "Median depth")
	cliflag.BoolVar(fs, &doMin, "", "min", false, "Min depth")
	cliflag.BoolVar(fs, &doMax, "", "max", false, "Max depth")
	cliflag.BoolVar(fs, &doSum, "", "sum", false, "Sum of per-base depths")

	var sameStrand, oppStrand bool
	cliflag.BoolVar(fs, &sameStrand, "s", "strand", false, "Same strand")
	cliflag.BoolVar(fs, &oppStrand, "S", "opposite", false, "Opposite strand")

	var fracA, fracB float64
	cliflag.Float64Var(fs, &fracA, "f", "fraction-a", 0, "Min fraction of A")
	cliflag.Float64Var(fs, &fracB, "F", "fraction-b", 0, "Min fraction of B")

	var reciprocal bool
	cliflag.BoolVar(fs, &reciprocal, "r", "reciprocal", false, "Require reciprocal -f overlap (A AND B)")
	var either bool
	cliflag.BoolVar(fs, &either, "e", "either", false, "Require -f OR -F (instead of AND)")

	var split bool
	cliflag.BoolVar(fs, &split, "", "split", false, "Treat BED12 -b records as their blocks")

	var help, showVersion bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version")

	flag.Parse()

	// Fold the legacy BAM aliases onto -a / -b. A given -abam/-ibam/-bbam
	// supplies the corresponding input when the modern flag was not used.
	if inputA == "" {
		if abam != "" {
			inputA = abam
		} else if ibam != "" {
			inputA = ibam
		}
	}
	if inputB == "" && bbam != "" {
		inputB = bbam
	}

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

	// Resolve mode. Numeric ops, hist, depth, and counts are mutually
	// exclusive — last-one-wins on the CLI is too surprising, so reject
	// combinations explicitly.
	mode, err := resolveMode(counts, depth, hist, doMean, doMedian, doMin, doMax, doSum)
	if err != nil {
		// Match upstream bedtools' exact stderr text for this error so the
		// message is drop-in compatible (the upstream coverage tool prints a
		// single fixed line regardless of which pair was combined).
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	if sameStrand && oppStrand {
		fmt.Fprintln(os.Stderr, "Error: -s and -S are mutually exclusive")
		os.Exit(1)
	}

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

	opts := bedcoverage.Options{
		Mode:           mode,
		SameStrand:     sameStrand,
		OppositeStrand: oppStrand,
		FractionA:      fracA,
		FractionB:      fracB,
		Reciprocal:     reciprocal,
		Either:         either,
		Split:          split,
	}
	if _, err := bedcoverage.Coverage(rA, rB, w, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

// errMutuallyExclusiveModes is the fixed error upstream `bedtools coverage`
// prints (verbatim, including the surrounding asterisks) when more than one of
// -counts / -d / -mean / -hist is requested. It is reproduced byte-for-byte so
// the CLI is drop-in compatible with scripts that match upstream's stderr.
var errMutuallyExclusiveModes = errors.New("***** ERROR: -counts, -d, -mean, and -hist are all mutually exclusive options. *****")

// resolveMode picks at most one output mode from the mutually exclusive set.
func resolveMode(counts, depth, hist, mean, median, mn, mx, sum bool) (bedcoverage.Mode, error) {
	flags := []struct {
		on   bool
		mode bedcoverage.Mode
		name string
	}{
		{counts, bedcoverage.ModeCounts, "-counts"},
		{depth, bedcoverage.ModeDepth, "-d"},
		{hist, bedcoverage.ModeHist, "-hist"},
		{mean, bedcoverage.ModeMean, "-mean"},
		{median, bedcoverage.ModeMedian, "-median"},
		{mn, bedcoverage.ModeMin, "-min"},
		{mx, bedcoverage.ModeMax, "-max"},
		{sum, bedcoverage.ModeSum, "-sum"},
	}
	var chosen []string
	var mode bedcoverage.Mode = bedcoverage.ModeDefault
	for _, f := range flags {
		if f.on {
			chosen = append(chosen, f.name)
			mode = f.mode
		}
	}
	if len(chosen) > 1 {
		return 0, errMutuallyExclusiveModes
	}
	return mode, nil
}
