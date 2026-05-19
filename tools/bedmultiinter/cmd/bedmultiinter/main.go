// bedmultiinter performs a multi-way intersection across N BED files
// (Go port of `bedtools multiinter`).
package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedmultiinter/pkg/bedmultiinter"
)

const version = "1.0.0"

const usage = `bedmultiinter - Multi-way intersection across N BED files

Usage:
  bedmultiinter -i <FILE1> <FILE2> [<FILE3> ...] [options]

Options:
  -i FILE..              N input BED files (variadic; required, '-' for stdin).
  -names CSV             Comma-separated labels (one per input). Default:
                         per-file 1-based index in the 'list' column and the
                         filename in the header row.
  -empty                 Emit '0'-count regions at chrom heads, tails, and
                         gaps between covered spans. Requires -g.
  -g FILE                Chrom-sizes genome file (required with -empty).
  -cluster               Collapse adjacent same-active-set segments.
  -header                Emit a column-header row before the data.
  -filler TEXT           Indicator for "file not contributing" cells.
                         Default: '0'.
  -o, --output FILE      Output file (default: stdout).
  -h, --help             Show this help.
  -v, --version          Show version.

Output columns:
  chrom  start  end  num  list  <per-file 0/1>

Examples:
  bedmultiinter -i a.bed b.bed c.bed
  bedmultiinter -i a.bed b.bed c.bed -header -names A,B,C
  bedmultiinter -i a.bed b.bed c.bed -empty -g sizes.txt -header
`

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, "Error:", err)
		os.Exit(1)
	}
}

func run(argv []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bedmultiinter", flag.ContinueOnError)
	fs.SetOutput(stderr)

	var (
		names         string
		genome        string
		filler        string
		output        string
		empty         bool
		cluster       bool
		header        bool
		help, showVer bool
	)

	cliflag.StringVar(fs, &names, "", "names", "", "Comma-separated labels")
	cliflag.StringVar(fs, &genome, "g", "genome", "", "Chrom-sizes file")
	cliflag.StringVar(fs, &filler, "", "filler", "0", "Indicator for absent files")
	cliflag.BoolVar(fs, &empty, "", "empty", false, "Emit 0-count regions")
	cliflag.BoolVar(fs, &cluster, "", "cluster", false, "Collapse same-set runs")
	cliflag.BoolVar(fs, &header, "", "header", false, "Emit header row")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout)")
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVer, "v", "version", false, "Show version")

	fs.Usage = func() { fmt.Fprint(stderr, usage) }

	// Variadic -i pre-extracted before the rest of the flags are parsed.
	filesPaths, argvRest := extractVarArg(argv, []string{"-i", "--i"})

	if err := fs.Parse(argvRest); err != nil {
		return err
	}
	if help {
		fmt.Fprint(stderr, usage)
		return nil
	}
	if showVer {
		fmt.Fprintf(stdout, "bedmultiinter version %s\n", version)
		return nil
	}
	if len(filesPaths) == 0 {
		return fmt.Errorf("at least one -i <FILE> is required (use -h for help)")
	}
	if len(filesPaths) < 2 {
		return fmt.Errorf("multiinter needs >=2 input files; got %d", len(filesPaths))
	}
	if empty && genome == "" {
		return fmt.Errorf("-empty requires -g <genome-file>")
	}

	// Resolve names: explicit -names overrides; otherwise no per-row
	// labels (numeric indices in the list column) but the header (if
	// requested) shows raw filenames as supplied.
	var nameSlice []string
	if names != "" {
		nameSlice = strings.Split(names, ",")
		if len(nameSlice) != len(filesPaths) {
			return fmt.Errorf("-names supplies %d labels but -i has %d files",
				len(nameSlice), len(filesPaths))
		}
	}

	// Optional chrom-size map for -empty.
	var sizes map[string]int
	if empty {
		gR, err := iohelper.OpenReader(genome)
		if err != nil {
			return fmt.Errorf("opening genome file: %w", err)
		}
		defer gR.Close()
		sizes, err = bedmultiinter.ReadGenomeSizes(gR)
		if err != nil {
			return fmt.Errorf("parsing genome file: %w", err)
		}
	}

	// Open the input files.
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

	opts := bedmultiinter.Options{
		Names:      nameSlice,
		Filenames:  filesPaths,
		Empty:      empty,
		ChromSizes: sizes,
		Cluster:    cluster,
		Header:     header,
		Filler:     filler,
	}
	if _, err := bedmultiinter.Run(bRs, w, opts); err != nil {
		return err
	}
	return nil
}

// extractVarArg pulls values following any of the given trigger flags
// until the next dash-prefixed token (or EOF).
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
