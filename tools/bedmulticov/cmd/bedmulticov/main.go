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
  -files FILE..          One or more B files (BED or BAM, auto-detected by
                         extension). Variadic.
  -bams FILE..           Alias for -files (kept for upstream compatibility).
  -q,  --mapq N          Minimum MAPQ for BAM inputs (default 0; ignored for
                         BED inputs).
  -D,  --max-depth N     Cap per-A-interval depth count for BAM inputs
                         (default 64000; ignored for BED inputs).
  -s,  --strand          Same-strand overlaps only.
  -S,  --opposite        Opposite-strand overlaps only.
  -f FRACTION            Minimum fraction of A overlapped (0,1].
  -F FRACTION            Minimum fraction of B overlapped (0,1].
  -r,  --reciprocal      Apply -f to BOTH A and B.
       --split           BAM CIGAR N-op block-aware coverage: count only
                         contiguous reference-consuming op runs (M/=/X/D),
                         skipping spanning N gaps. Ignored for BED inputs.
  -o,  --output FILE     Output file (default: stdout).
  -h,  --help            Show this help.
  -v,  --version         Show version.

Notes:
  - BAM input is supported (BGZF-wrapped, decoded via pkg/bioformats/sam).
  - CRAM input is NOT yet supported — a clear error is surfaced for any
    .cram path; see docs/CRAM_DESIGN.md.
  - The output preserves A's columns verbatim and appends one integer
    count column per input file, matching upstream's column ordering.

Examples:
  bedmulticov -bed a.bed -files b1.bed b2.bed
  bedmulticov -bed a.bed -bams b1.bam b2.bam -q 20
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
		split                 bool
		fractionA, fractionB  float64
		help, showVer         bool
	)

	cliflag.StringVar(fs, &bedPath, "", "bed", "", "A BED file (required)")
	cliflag.IntVar(fs, &mapq, "q", "mapq", 0, "Min MAPQ for BAM inputs")
	cliflag.IntVar(fs, &maxDepth, "D", "max-depth", 64000, "Cap per-A-interval depth count for BAM inputs (0 disables)")
	cliflag.BoolVar(fs, &sameStrand, "s", "strand", false, "Same-strand overlaps only")
	cliflag.BoolVar(fs, &oppStrand, "S", "opposite", false, "Opposite-strand overlaps only")
	cliflag.Float64Var(fs, &fractionA, "f", "fraction-a", 0.0, "Fraction of A overlapped")
	cliflag.Float64Var(fs, &fractionB, "F", "fraction-b", 0.0, "Fraction of B overlapped")
	cliflag.BoolVar(fs, &reciprocal, "r", "reciprocal", false, "Apply -f to both A and B")
	cliflag.BoolVar(fs, &split, "", "split", false, "BAM CIGAR N-op block-aware coverage")
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
	// Reject CRAM explicitly — we don't have a CRAM reader yet.
	for _, p := range all {
		if strings.HasSuffix(strings.ToLower(p), ".cram") {
			return fmt.Errorf("CRAM input not yet supported: %q (see docs/CRAM_DESIGN.md)", p)
		}
	}

	aR, err := iohelper.OpenReader(bedPath)
	if err != nil {
		return fmt.Errorf("opening A: %w", err)
	}
	defer aR.Close()

	sources := make([]bedmulticov.Source, 0, len(all))
	closers := make([]io.Closer, 0, len(all))
	defer func() {
		for _, c := range closers {
			c.Close()
		}
	}()
	for _, p := range all {
		if strings.HasSuffix(strings.ToLower(p), ".bam") {
			// BAM is BGZF-wrapped; iohelper would auto-decode the BGZF
			// layer and break sam.NewBAMReader. Open raw.
			f, err := os.Open(p)
			if err != nil {
				return fmt.Errorf("opening input %q: %w", p, err)
			}
			closers = append(closers, f)
			sources = append(sources, bedmulticov.Source{Reader: f, Kind: bedmulticov.SourceBAM})
			continue
		}
		f, err := iohelper.OpenReader(p)
		if err != nil {
			return fmt.Errorf("opening input %q: %w", p, err)
		}
		closers = append(closers, f)
		sources = append(sources, bedmulticov.Source{Reader: f, Kind: bedmulticov.SourceBED})
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
		MinMAPQ:        mapq,
		MaxDepth:       maxDepth,
		Split:          split,
	}
	if _, err := bedmulticov.RunSources(aR, sources, w, opts); err != nil {
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
