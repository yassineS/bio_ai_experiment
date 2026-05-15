// bedmulticov reports the per-input overlap count for each interval in a
// primary BED file (Go port of `bedtools multicov`).
package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedmulticov/pkg/bedmulticov"
)

const version = "1.0.0"

const usage = `bedmulticov - Per-interval overlap counts against N BED/BAM files

Usage:
  bedmulticov -bed <A.bed> -files <B1.bed> [<B2.bed> ...] [options]
  bedmulticov -bed <A.bed> -bams  <B1.bam> [<B2.bam> ...] [options]

Options:
  -bed FILE              A intervals (required, '-' for stdin)
  -files FILE..          One or more B files (BED). Variadic.
  -bams FILE..           Alias for -files; reserved for BAM inputs (not yet
                         implemented — see README).
  -q,  --mapq N          Minimum MAPQ for BAM inputs (default 0; ignored for
                         BED inputs).
  -D,  --max-depth N     Cap per-position depth for BAM inputs (default 64000;
                         ignored for BED inputs).
  -s,  --strand          Same-strand overlaps only.
  -S,  --opposite        Opposite-strand overlaps only.
  -f FRACTION            Minimum fraction of A overlapped (0,1].
  -F FRACTION            Minimum fraction of B overlapped (0,1].
  -r,  --reciprocal      Apply -f to BOTH A and B.
  -o,  --output FILE     Output file (default: stdout).
  -h,  --help            Show this help.
  -v,  --version         Show version.

Notes:
  - BAM inputs are not yet implemented; supplying a .bam path returns an
    explicit error.
  - The output preserves A's columns verbatim and appends one integer
    count column per input file, matching upstream's column ordering.

Examples:
  bedmulticov -bed a.bed -files b1.bed b2.bed
  bedmulticov -bed a.bed -files b1.bed b2.bed -s -f 0.5
`

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, "Error:", err)
		os.Exit(1)
	}
}

func run(argv []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bedmulticov", flag.ContinueOnError)
	fs.SetOutput(stderr)

	var (
		bedPath, output       string
		mapq, maxDepth        int
		sameStrand, oppStrand bool
		reciprocal            bool
		fractionA, fractionB  float64
		help, showVer         bool
	)

	cliflag.StringVar(fs, &bedPath, "", "bed", "", "A BED file (required)")
	cliflag.IntVar(fs, &mapq, "q", "mapq", 0, "Min MAPQ for BAM inputs")
	cliflag.IntVar(fs, &maxDepth, "D", "max-depth", 64000, "Max per-position depth for BAM inputs")
	cliflag.BoolVar(fs, &sameStrand, "s", "strand", false, "Same-strand overlaps only")
	cliflag.BoolVar(fs, &oppStrand, "S", "opposite", false, "Opposite-strand overlaps only")
	cliflag.Float64Var(fs, &fractionA, "f", "fraction-a", 0.0, "Fraction of A overlapped")
	cliflag.Float64Var(fs, &fractionB, "F", "fraction-b", 0.0, "Fraction of B overlapped")
	cliflag.BoolVar(fs, &reciprocal, "r", "reciprocal", false, "Apply -f to both A and B")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout)")
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVer, "v", "version", false, "Show version")

	fs.Usage = func() { fmt.Fprint(stderr, usage) }

	// Pre-extract variadic -files / -bams paths. Everything between a
	// trigger flag and the next dash-prefixed token (or EOF) is a file.
	filesPaths, argv1 := extractVarArg(argv, []string{"-files", "--files"})
	bamsPaths, argvRest := extractVarArg(argv1, []string{"-bams", "--bams"})
	all := append([]string{}, filesPaths...)
	all = append(all, bamsPaths...)

	if err := fs.Parse(argvRest); err != nil {
		return err
	}
	if help {
		fmt.Fprint(stderr, usage)
		return nil
	}
	if showVer {
		fmt.Fprintf(stdout, "bedmulticov version %s\n", version)
		return nil
	}
	if bedPath == "" {
		return fmt.Errorf("-bed is required (use -h for help)")
	}
	if len(all) == 0 {
		return fmt.Errorf("at least one -files or -bams <FILE> is required (use -h for help)")
	}
	// Surface BAM paths explicitly — they are not yet supported.
	for _, p := range all {
		lp := strings.ToLower(p)
		if strings.HasSuffix(lp, ".bam") || strings.HasSuffix(lp, ".cram") {
			return fmt.Errorf("BAM/CRAM input not yet implemented: %q (see README)", p)
		}
	}
	// `-q` and `-D` only meaningful for BAM; tolerated for BED so the CLI
	// matches upstream's accepted flag set.
	_ = mapq
	_ = maxDepth

	aR, err := iohelper.OpenReader(bedPath)
	if err != nil {
		return fmt.Errorf("opening A: %w", err)
	}
	defer aR.Close()

	bRs := make([]io.Reader, 0, len(all))
	closers := make([]io.Closer, 0, len(all))
	defer func() {
		for _, c := range closers {
			c.Close()
		}
	}()
	for _, p := range all {
		f, err := iohelper.OpenReader(p)
		if err != nil {
			return fmt.Errorf("opening input %q: %w", p, err)
		}
		closers = append(closers, f)
		bRs = append(bRs, f)
	}

	w, err := iohelper.OpenWriter(output)
	if err != nil {
		return fmt.Errorf("opening output: %w", err)
	}
	defer w.Close()

	opts := bedmulticov.Options{
		FractionA:      fractionA,
		FractionB:      fractionB,
		Reciprocal:     reciprocal,
		SameStrand:     sameStrand,
		OppositeStrand: oppStrand,
	}
	if _, err := bedmulticov.Run(aR, bRs, w, opts); err != nil {
		return err
	}
	return nil
}

// extractVarArg pulls the values that follow any of the given trigger
// flags. Values continue until the next argument that starts with '-'
// or EOF. The unused tail is returned for flag.Parse to consume.
func extractVarArg(argv, triggers []string) (values, rest []string) {
	rest = make([]string, 0, len(argv))
	for i := 0; i < len(argv); i++ {
		a := argv[i]
		matched := false
		for _, t := range triggers {
			if a == t {
				matched = true
				break
			}
		}
		if !matched {
			rest = append(rest, a)
			continue
		}
		i++
		for ; i < len(argv); i++ {
			v := argv[i]
			if len(v) > 0 && v[0] == '-' {
				i--
				break
			}
			values = append(values, v)
		}
	}
	return values, rest
}
