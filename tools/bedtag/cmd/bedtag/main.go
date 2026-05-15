// bedtag annotates each interval in A with the name (or another column) of
// any overlapping interval in B (Go port of `bedtools tag`).
package main

import (
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedtag/pkg/bedtag"
)

const usage = `bedtag - Annotate A intervals with B's name column

Usage:
  bedtag -a A.bed -b B.bed[,B2.bed,...] [options]

Options:
  -a, --input-a FILE      Input BED file A (required)
  -b, --input-b LIST      One or more comma-separated B BED files (required)
  -o, --output FILE       Output file (default: stdout)
  -i, --tag INT           1-based column from B to use as the tag (default 4)
  --names LIST            Comma-separated names; replaces B's tag column.
                          Length must equal the number of B files.
  --labels                Prefix each tag with "<bfile>=".
  -s, --strand            Only same-strand B records contribute tags
  -S, --inverse-strand    Only opposite-strand B records contribute tags
  -m, --min-overlap INT   Minimum bp overlap to consider (default 1)
  -f, --fraction-a NUM    Minimum fraction of A that must overlap (0-1)
  -F, --fraction-b NUM    Minimum fraction of B that must overlap (0-1)
  -h, --help              Show help
  -v, --version           Show version

Output:
  Each A line, verbatim, with one extra TSV column listing the comma-separated
  tags from overlapping B records (empty when none).
`

const version = "0.1.0"

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func run(args []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bedtag", flag.ContinueOnError)
	fs.SetOutput(stderr)
	fs.Usage = func() { fmt.Fprint(stderr, usage) }

	var (
		inputA     string
		inputB     string
		output     string
		tagCol     int
		namesStr   string
		labels     bool
		strand     bool
		inverseStr bool
		minOverlap int
		fractionA  float64
		fractionB  float64
		showHelp   bool
		showVer    bool
	)
	cliflag.StringVar(fs, &inputA, "a", "input-a", "", "BED A")
	cliflag.StringVar(fs, &inputB, "b", "input-b", "", "BED B (comma list)")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output")
	cliflag.IntVar(fs, &tagCol, "i", "tag", 4, "1-based column for tag value in B")
	fs.StringVar(&namesStr, "names", "", "Replacement names")
	fs.BoolVar(&labels, "labels", false, "Prefix tags with B file name")
	cliflag.BoolVar(fs, &strand, "s", "strand", false, "Same-strand only")
	cliflag.BoolVar(fs, &inverseStr, "S", "inverse-strand", false, "Opposite-strand only")
	cliflag.IntVar(fs, &minOverlap, "m", "min-overlap", 1, "Min overlap bp")
	cliflag.Float64Var(fs, &fractionA, "f", "fraction-a", 0, "Min fraction A")
	cliflag.Float64Var(fs, &fractionB, "F", "fraction-b", 0, "Min fraction B")
	cliflag.BoolVar(fs, &showHelp, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVer, "v", "version", false, "Show version")

	if err := fs.Parse(args); err != nil {
		return err
	}
	if showHelp {
		fmt.Fprint(stdout, usage)
		return nil
	}
	if showVer {
		fmt.Fprintln(stdout, version)
		return nil
	}
	if inputA == "" || inputB == "" {
		fmt.Fprint(stderr, usage)
		return fmt.Errorf("error: -a and -b are required")
	}

	bFiles := strings.Split(inputB, ",")
	for i := range bFiles {
		bFiles[i] = strings.TrimSpace(bFiles[i])
	}
	var names []string
	if namesStr != "" {
		names = strings.Split(namesStr, ",")
		for i := range names {
			names[i] = strings.TrimSpace(names[i])
		}
		if len(names) != len(bFiles) {
			return fmt.Errorf("--names has %d entries but %d B files supplied", len(names), len(bFiles))
		}
	}

	aR, err := iohelper.OpenReader(inputA)
	if err != nil {
		return err
	}
	defer aR.Close()

	sources := make([]bedtag.Source, len(bFiles))
	for i, p := range bFiles {
		bR, err := iohelper.OpenReader(p)
		if err != nil {
			return fmt.Errorf("opening B file %s: %w", p, err)
		}
		defer bR.Close()
		sources[i] = bedtag.Source{Name: p, Reader: bR}
	}

	w, err := iohelper.OpenWriter(output)
	if err != nil {
		return err
	}
	defer w.Close()

	_, err = bedtag.Tag(aR, sources, w, bedtag.Options{
		TagColumn:     tagCol,
		Names:         names,
		Labels:        labels,
		StrandSpec:    strand,
		InverseStrand: inverseStr,
		MinOverlap:    minOverlap,
		FractionA:     fractionA,
		FractionB:     fractionB,
	})
	return err
}
