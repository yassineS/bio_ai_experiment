// bedannotate annotates A's intervals with overlap stats from N input
// BED files (Go port of `bedtools annotate`).
package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedannotate/pkg/bedannotate"
)

const version = "1.0.0"

const usage = `bedannotate - Annotate A intervals with overlap stats from N BED files

Usage:
  bedannotate -i <A.bed> -files <B1.bed> [<B2.bed> ...] [options]

Options:
  -i,  --input FILE        A intervals (required, '-' for stdin)
       --files FILE..      One or more B files to compute coverage against
       --names N1 N2 ..    Header labels (variadic; a single comma-separated
                           token is also accepted). Triggers a header line.
       --counts            Emit per-B count of overlapping records
       --both              Emit count and coverage fraction per B (interleaved)
  -s,  --strand            Restrict overlaps to same-strand pairs
  -S,  --opposite          Restrict overlaps to opposite-strand pairs
  -o,  --output FILE       Output file (default: stdout)
  -h,  --help              Show this help
  -v,  --version           Show version

Notes:
  - Default appends one fraction (in [0,1], %f-formatted) per B file.
  - -counts replaces fractions with integer overlap counts.
  - -both interleaves count + fraction per B (2N columns total).
  - A '#' header line is emitted ONLY when -names is given (matching upstream).
  - Records are reported grouped by chromosome then by UCSC bin.

Examples:
  bedannotate -i a.bed --files b1.bed b2.bed
  bedannotate -i a.bed --files b1.bed b2.bed --names exons introns --both
`

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, "Error:", err)
		os.Exit(1)
	}
}

func run(argv []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bedannotate", flag.ContinueOnError)
	fs.SetOutput(stderr)

	var (
		aPath, output      string
		counts, both       bool
		sameStrand, oppStr bool
		help, showVer      bool
	)

	cliflag.StringVar(fs, &aPath, "i", "input", "", "A BED file (required)")
	cliflag.BoolVar(fs, &counts, "", "counts", false, "Emit overlap counts")
	cliflag.BoolVar(fs, &both, "", "both", false, "Emit count + fraction per B")
	cliflag.BoolVar(fs, &sameStrand, "s", "strand", false, "Same-strand overlaps only")
	cliflag.BoolVar(fs, &oppStr, "S", "opposite", false, "Opposite-strand overlaps only")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout)")
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVer, "v", "version", false, "Show version")

	fs.Usage = func() { fmt.Fprintf(stderr, "%s", usage) }

	// Pre-extract the variadic `-files <FILE...>` and `-names <NAME...>` args,
	// which Go's flag package cannot model: everything between the trigger flag
	// and the next `-`-prefixed token (or EOF) is a value. `namesGiven` records
	// whether -names was present at all, since an empty -names still differs
	// from omitting it (a header is emitted only when -names is given).
	filesPaths, argvRest := extractVarArg(argv, []string{"-files", "--files"})
	nameValues, namesGiven, argvRest := extractVarArgFlagged(argvRest, []string{"-names", "--names"})

	if err := fs.Parse(argvRest); err != nil {
		return err
	}
	if help {
		fmt.Fprintf(stderr, "%s", usage)
		return nil
	}
	if showVer {
		fmt.Fprintf(stdout, "bedannotate version %s\n", version)
		return nil
	}
	if aPath == "" {
		return fmt.Errorf("-i / --input is required (use -h for help)")
	}
	if len(filesPaths) == 0 {
		return fmt.Errorf("at least one --files <FILE> is required (use -h for help)")
	}

	// Mode resolution.
	mode := bedannotate.ModeFraction
	if both {
		mode = bedannotate.ModeBoth
	} else if counts {
		mode = bedannotate.ModeCounts
	}

	// Names: a header line is emitted ONLY when -names is explicitly given
	// (matching upstream — file basenames do NOT trigger a header). Leave the
	// slice nil otherwise so no header is printed. Accept both the upstream
	// space-separated form (-names b1 b2) and a single comma-separated token
	// (--names b1,b2) for convenience.
	var nameSlice []string
	if namesGiven {
		if len(nameValues) == 1 && strings.Contains(nameValues[0], ",") {
			nameSlice = strings.Split(nameValues[0], ",")
		} else {
			nameSlice = nameValues
		}
		if len(nameSlice) == 0 {
			// -names given with no values: emit a header with empty labels,
			// one per file (matching upstream's empty-title behaviour).
			nameSlice = make([]string, len(filesPaths))
		}
	}

	aR, err := iohelper.OpenReader(aPath)
	if err != nil {
		return fmt.Errorf("opening A: %w", err)
	}
	defer aR.Close()

	bRs := make([]io.Reader, 0, len(filesPaths))
	closers := make([]io.Closer, 0, len(filesPaths))
	defer func() {
		for _, c := range closers {
			c.Close()
		}
	}()
	for _, p := range filesPaths {
		f, err := iohelper.OpenReader(p)
		if err != nil {
			return fmt.Errorf("opening B file %q: %w", p, err)
		}
		closers = append(closers, f)
		bRs = append(bRs, f)
	}

	w, err := iohelper.OpenWriter(output)
	if err != nil {
		return fmt.Errorf("opening output: %w", err)
	}
	defer w.Close()

	opts := bedannotate.Options{
		Mode:           mode,
		Names:          nameSlice,
		SameStrand:     sameStrand,
		OppositeStrand: oppStr,
	}
	if _, err := bedannotate.Run(aR, bRs, w, opts); err != nil {
		return err
	}
	return nil
}

// extractVarArgFlagged is like extractVarArg but also reports whether any of
// the trigger flags appeared at all (so an empty value list can be
// distinguished from the flag being absent).
func extractVarArgFlagged(argv, triggers []string) (values []string, found bool, rest []string) {
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
		found = true
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
	return values, found, rest
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
