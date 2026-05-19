// Command bedmakewindows is a pure-Go reimplementation of
// `bedtools makewindows`. It partitions either chromosomes from a genome
// file or arbitrary BED intervals into fixed-width or fixed-count windows.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedmakewindows/pkg/bedmakewindows"
)

const version = "0.1.0"

const usage = `bedmakewindows - partition intervals into windows.

Usage:
  bedmakewindows -g GENOME -w WIDTH [-s STEP] [-n COUNT] [-i NAMING] [-reverse] [-o OUT]
  bedmakewindows -b BED    -w WIDTH [-s STEP] [-n COUNT] [-i NAMING] [-reverse] [-o OUT]

Inputs (exactly one of -g/-b required):
  -g, --genome FILE       Genome chrom-sizes file (CHROM\tSIZE per line).
  -b, --bed FILE          BED file of source intervals (use '-' for stdin).

Window strategy (exactly one of -w/-n required):
  -w, --window-size N     Window width in bases.
  -s, --step-size N       Slide between windows (defaults to -w when omitted).
  -n, --count N           Partition each interval into N equal windows.

Naming:
  -i, --id-name TYPE      One of: 'srcwinnum' | 'winnum' | 'src' | 'none' (default 'none').
      --reverse           Reverse per-interval window numbering (last window = 1).

Output:
  -o, --output FILE       Output file (default stdout). '-' for stdout.

Standard:
  -h, --help              Show this help.
  -v, --version           Show version.
`

func main() {
	fs := flag.NewFlagSet("bedmakewindows", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		genome      string
		bed         string
		width       int
		step        int
		count       int
		idName      string
		reverse     bool
		outFile     string
		showHelp    bool
		showVersion bool
	)
	cliflag.StringVar(fs, &genome, "g", "genome", "", "Genome chrom-sizes file")
	cliflag.StringVar(fs, &bed, "b", "bed", "", "Source BED file")
	cliflag.IntVar(fs, &width, "w", "window-size", 0, "Window size")
	cliflag.IntVar(fs, &step, "s", "step-size", 0, "Step")
	cliflag.IntVar(fs, &count, "n", "count", 0, "Window count")
	cliflag.StringVar(fs, &idName, "i", "id-name", "none", "Naming strategy")
	fs.BoolVar(&reverse, "reverse", false, "Reverse numbering")
	cliflag.StringVar(fs, &outFile, "o", "output", "", "Output path")
	fs.BoolVar(&showHelp, "h", false, "Help")
	fs.BoolVar(&showHelp, "help", false, "Help")
	fs.BoolVar(&showVersion, "v", false, "Version")
	fs.BoolVar(&showVersion, "version", false, "Version")

	if err := fs.Parse(os.Args[1:]); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, usage)
		os.Exit(2)
	}
	if showHelp {
		fmt.Print(usage)
		return
	}
	if showVersion {
		fmt.Println(version)
		return
	}
	if (genome == "" && bed == "") || (genome != "" && bed != "") {
		fmt.Fprintln(os.Stderr, "bedmakewindows: exactly one of -g / -b required")
		fmt.Fprint(os.Stderr, usage)
		os.Exit(2)
	}
	if (width == 0 && count == 0) || (width > 0 && count > 0) {
		fmt.Fprintln(os.Stderr, "bedmakewindows: exactly one of -w / -n required")
		os.Exit(2)
	}

	naming, err := bedmakewindows.ParseNaming(idName)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bedmakewindows: %v\n", err)
		os.Exit(2)
	}

	var intervals []bedmakewindows.Interval
	switch {
	case genome != "":
		f, err := iohelper.OpenReader(genome)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bedmakewindows: open genome: %v\n", err)
			os.Exit(1)
		}
		defer f.Close()
		intervals, err = bedmakewindows.FromGenome(f)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bedmakewindows: parse genome: %v\n", err)
			os.Exit(1)
		}
	case bed != "":
		f, err := iohelper.OpenReader(bed)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bedmakewindows: open bed: %v\n", err)
			os.Exit(1)
		}
		defer f.Close()
		intervals, err = bedmakewindows.FromBED(f)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bedmakewindows: parse bed: %v\n", err)
			os.Exit(1)
		}
	}

	out, err := iohelper.OpenWriter(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bedmakewindows: open output: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if _, err := bedmakewindows.MakeWindows(intervals, out, os.Stderr, bedmakewindows.Options{
		Width:   width,
		Step:    step,
		Count:   count,
		Reverse: reverse,
		Naming:  naming,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "bedmakewindows: %v\n", err)
		os.Exit(1)
	}
}
