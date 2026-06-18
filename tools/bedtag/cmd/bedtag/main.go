// bedtag annotates each interval in A with the name (or another column) of
// any overlapping interval in B (Go port of `bedtools tag`).
package main

import (
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedtag/pkg/bedtag"
)

const usage = `bedtag - Annotate alignments/intervals from overlaps (Go port of bedtools tag).

Upstream tagBam mode (BAM in, BAM out — selected by -files):
  bedtag -i <BAM> -files F1 [F2 ...] (-labels L1 [L2 ...] | -names | -scores) [options]
    -i FILE        Input BAM.
    -files F...    Annotation BED/GFF/VCF files (space-separated).
    -labels L...   Per-file label written to the tag on any overlap (default mode).
    -names         Tag with the overlapping records' name fields instead.
    -scores        Tag with the overlapping records' score fields instead.
    -tag XY        Two-character aux tag to write (default YB).
    -s | -S        Require same / opposite strand overlaps.
    -f FLOAT       Min overlap as a fraction of the alignment (default 1e-9).

BED-in/BED-out extension mode:
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
	args := os.Args[1:]
	// The upstream `bedtools tag` model — `-i <BAM> -files F1 .. -labels L1 ..`
	// — is selected by the presence of -files (its variadic flags can't be
	// expressed with the standard flag package, so it has a dedicated parser).
	// Without -files we fall back to the BED-in/BED-out extension below.
	if hasFlag(args, "-files") {
		if err := runTagBAM(args, os.Stdout, os.Stderr); err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		return
	}
	if err := run(args, os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

// hasFlag reports whether token appears as an argument.
func hasFlag(args []string, token string) bool {
	for _, a := range args {
		if a == token {
			return true
		}
	}
	return false
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

// runTagBAM implements the upstream `bedtools tag` (tagBam) model: tag a BAM's
// alignments from overlaps with annotation files. Its -files / -labels flags
// are variadic (space-separated), so the argv is scanned manually, mirroring
// tagBamMain.cpp.
func runTagBAM(args []string, stdout, stderr *os.File) error {
	var (
		bamFile   string
		files     []string
		labels    []string
		tag       = "YB"
		fraction  float64
		useNames  bool
		useScores bool
		sameStr   bool
		diffStr   bool
	)
	for i := 0; i < len(args); i++ {
		switch args[i] {
		case "-i":
			if i+1 < len(args) {
				bamFile = args[i+1]
				i++
			}
		case "-files":
			for i+1 < len(args) && !strings.HasPrefix(args[i+1], "-") {
				files = append(files, args[i+1])
				i++
			}
		case "-labels":
			for i+1 < len(args) && !strings.HasPrefix(args[i+1], "-") {
				labels = append(labels, args[i+1])
				i++
			}
		case "-names":
			useNames = true
		case "-scores":
			useScores = true
		case "-s":
			sameStr = true
		case "-S":
			diffStr = true
		case "-f":
			if i+1 < len(args) {
				f, err := strconv.ParseFloat(args[i+1], 64)
				if err != nil {
					return fmt.Errorf("invalid -f %q: %w", args[i+1], err)
				}
				fraction = f
				i++
			}
		case "-tag":
			if i+1 < len(args) {
				tag = args[i+1]
				i++
			}
		case "-h", "--help":
			fmt.Fprint(stdout, usage)
			return nil
		case "-v", "--version":
			fmt.Fprintln(stdout, version)
			return nil
		}
	}
	if bamFile == "" {
		return fmt.Errorf("error: bedtools tag requires -i <BAM>")
	}
	if len(files) == 0 {
		return fmt.Errorf("error: bedtools tag requires -files")
	}
	if !useNames && !useScores && len(labels) == 0 {
		return fmt.Errorf("error: need -labels or -names or -scores")
	}

	mode := bedtag.TagModeLabels
	switch {
	case useNames:
		mode = bedtag.TagModeNames
	case useScores:
		mode = bedtag.TagModeScores
	}

	annoFiles := make([][]*bed.Record, len(files))
	for i, p := range files {
		br, err := iohelper.OpenReader(p)
		if err != nil {
			return fmt.Errorf("opening annotation file %s: %w", p, err)
		}
		recs, err := bed.NewReader(br).ReadAll()
		_ = br.Close()
		if err != nil {
			return fmt.Errorf("reading annotation file %s: %w", p, err)
		}
		annoFiles[i] = recs
	}

	// BAM is BGZF; sam.NewBAMReader handles the block decompression itself, so
	// open the file raw rather than through iohelper (which would gunzip it and
	// hand the reader already-inflated bytes).
	in, err := os.Open(bamFile)
	if err != nil {
		return err
	}
	defer in.Close()

	_, err = bedtag.TagBAM(in, annoFiles, stdout, bedtag.TagBAMOptions{
		Tag:            tag,
		Mode:           mode,
		Labels:         labels,
		SameStrand:     sameStr,
		OppositeStrand: diffStr,
		MinFraction:    fraction,
	})
	return err
}
