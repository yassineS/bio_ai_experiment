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
       --names N1,N2,..    Comma-separated header labels (defaults: file basenames)
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
  - When -names is supplied (or files are supplied via --files), a header
    row prefixed with '#' is emitted before the data.

Examples:
  bedannotate -i a.bed --files b1.bed b2.bed
  bedannotate -i a.bed --files b1.bed b2.bed --names exons,introns --both
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
		names              string
		counts, both       bool
		sameStrand, oppStr bool
		help, showVer      bool
	)

	cliflag.StringVar(fs, &aPath, "i", "input", "", "A BED file (required)")
	cliflag.StringVar(fs, &names, "", "names", "", "Comma-separated header labels")
	cliflag.BoolVar(fs, &counts, "", "counts", false, "Emit overlap counts")
	cliflag.BoolVar(fs, &both, "", "both", false, "Emit count + fraction per B")
	cliflag.BoolVar(fs, &sameStrand, "s", "strand", false, "Same-strand overlaps only")
	cliflag.BoolVar(fs, &oppStr, "S", "opposite", false, "Opposite-strand overlaps only")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout)")
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVer, "v", "version", false, "Show version")

	fs.Usage = func() { fmt.Fprintf(stderr, "%s", usage) }

	// Pre-extract --files / -files <FILE...> from argv: it's variadic,
	// which Go's flag package can't model directly. Everything between
	// `-files` (or `--files`) and the next `-`-prefixed token (or EOF)
	// is a B-file path.
	filesPaths, argvRest := extractVarArg(argv, []string{"-files", "--files"})

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

	// Names: explicit -names overrides; otherwise derive from basenames.
	var nameSlice []string
	if names != "" {
		nameSlice = strings.Split(names, ",")
	} else {
		nameSlice = bedannotate.DefaultNames(filesPaths)
	}
	if len(nameSlice) != len(filesPaths) {
		return fmt.Errorf("--names supplies %d labels but --files has %d entries",
			len(nameSlice), len(filesPaths))
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
