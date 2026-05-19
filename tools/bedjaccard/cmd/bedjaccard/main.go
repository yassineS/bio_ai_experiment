// bedjaccard computes the Jaccard similarity between two sorted BED files.
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedjaccard/pkg/bedjaccard"
)

const version = "1.0.0"

const usage = `bedjaccard - Compute the Jaccard similarity of two sorted BED files

Usage:
  bedjaccard -a <fileA.bed> -b <fileB.bed> [options]

Description:
  Reads two BED files A and B (both pre-sorted by chrom then start) and writes
  a two-line tab-separated summary to stdout:

      intersection<TAB>union<TAB>jaccard<TAB>n_intersections
      <int>        <int>  <float>  <int>

  Where intersection is the total number of bases covered by both A and B,
  union = |A| + |B| - intersection, jaccard = intersection / union (0 if union
  is 0), and n_intersections is the count of overlapping interval pairs.

Options:
  -a, --a FILE          First sorted BED file (required, '-' for stdin)
  -b, --b FILE          Second sorted BED file (required)
      --output FILE     Output file (default: stdout)
  -s, --strand          Same-strand overlaps only (BED6 strand column)
  -S, --opposite-strand Opposite-strand overlaps only (BED6 strand column)
  -f FRACTION           Require >= FRACTION of A overlapped by B (0..1)
  -F FRACTION           Require >= FRACTION of B overlapped by A (0..1)
  -h, --help            Show this help message
  -v, --version         Show version information

Examples:
  bedjaccard -a peaksA.bed -b peaksB.bed
  bedjaccard --a a.bed --b b.bed -f 0.5
  bedjaccard -a a.bed -b b.bed -s

Notes:
  - Both inputs MUST be sorted on (chrom, start); the tool errors out otherwise.
  - Coordinates are 0-based, half-open.
  - Use '-' for stdin and '--' to end option parsing.
`

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, "Error:", err)
		os.Exit(1)
	}
}

// run is the testable entry point.
func run(argv []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bedjaccard", flag.ContinueOnError)
	fs.SetOutput(stderr)

	var (
		fileA, fileB, output string
		same, opposite       bool
		fractionA, fractionB float64
		help, showVer        bool
	)

	// `-a/--a` and `-b/--b` use a single-letter name; Go's `flag` accepts both
	// `-a` and `--a` for any flag, so we register it once.
	cliflag.StringVar(fs, &fileA, "a", "", "", "First sorted BED file (required)")
	cliflag.StringVar(fs, &fileB, "b", "", "", "Second sorted BED file (required)")
	cliflag.StringVar(fs, &output, "", "output", "", "Output file (default: stdout)")

	cliflag.BoolVar(fs, &same, "s", "strand", false, "Same-strand overlaps only")
	cliflag.BoolVar(fs, &opposite, "S", "opposite-strand", false, "Opposite-strand overlaps only")
	cliflag.Float64Var(fs, &fractionA, "f", "fraction-a", 0.0, "Fraction of A overlapped (0..1)")
	cliflag.Float64Var(fs, &fractionB, "F", "fraction-b", 0.0, "Fraction of B overlapped (0..1)")

	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVer, "v", "version", false, "Show version")

	fs.Usage = func() { fmt.Fprint(stderr, usage) }

	if err := fs.Parse(argv); err != nil {
		return err
	}
	if help {
		fmt.Fprint(stderr, usage)
		return nil
	}
	if showVer {
		fmt.Fprintf(stdout, "bedjaccard version %s\n", version)
		return nil
	}

	if fileA == "" || fileB == "" {
		return fmt.Errorf("both -a/--a and -b/--b are required (use -h for help)")
	}

	rA, err := iohelper.OpenReader(fileA)
	if err != nil {
		return fmt.Errorf("opening A: %w", err)
	}
	defer rA.Close()
	rB, err := iohelper.OpenReader(fileB)
	if err != nil {
		return fmt.Errorf("opening B: %w", err)
	}
	defer rB.Close()
	w, err := iohelper.OpenWriter(output)
	if err != nil {
		return fmt.Errorf("opening output: %w", err)
	}
	defer w.Close()

	opts := bedjaccard.Options{
		SameStrand:     same,
		OppositeStrand: opposite,
		FractionA:      fractionA,
		FractionB:      fractionB,
	}
	if _, err := bedjaccard.Run(rA, rB, w, opts); err != nil {
		return err
	}
	return nil
}
