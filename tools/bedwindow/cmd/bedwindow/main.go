// bedwindow finds A intervals overlapping B intervals after expanding B by
// a window (Go port of `bedtools window`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedwindow/pkg/bedwindow"
)

const usage = `bedwindow - Overlap A with B after expanding B by a window

Usage:
  bedwindow -a A.bed -b B.bed [options]

Options:
  -a, --input-a FILE      Input BED file A (required)
  -b, --input-b FILE      Input BED file B (required)
  -o, --output FILE       Output file (default: stdout)
  -w, --window INT        Extend B intervals by N bp on both sides (default 0)
  -l, --left INT          Extend B intervals to the left only (overrides -w)
  -r, --right INT         Extend B intervals to the right only (overrides -w)
  -sm                     Same-strand only
  -sw                     Opposite-strand only
  -wa, --write-a          Write the original A entry only
  -wb, --write-b          Write the original B entry only
  -u, --unique            (Same as -wa: write each A at most once if it has any hit)
  -c, --count             Append B-hit count to A
  -v, --invert            Report only A entries with no B overlap
  -m, --min-overlap INT   Minimum bp overlap to count (default 1)
  -h, --help              Show help
  -v, --version           Show version

Default output:
  A<TAB>B for each (A, B) overlap pair (matches upstream).
`

const version = "0.1.0"

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func run(args []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bedwindow", flag.ContinueOnError)
	fs.SetOutput(stderr)
	fs.Usage = func() { fmt.Fprint(stderr, usage) }

	var (
		inputA   string
		inputB   string
		output   string
		window   int
		left     int
		right    int
		sm       bool
		sw       bool
		writeA   bool
		writeB   bool
		count    bool
		invert   bool
		minOL    int
		showHelp bool
		showVer  bool
	)
	cliflag.StringVar(fs, &inputA, "a", "input-a", "", "BED A")
	cliflag.StringVar(fs, &inputB, "b", "input-b", "", "BED B")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output")
	cliflag.IntVar(fs, &window, "w", "window", 0, "Window bp")
	cliflag.IntVar(fs, &left, "l", "left", -1, "Left extension bp (-1 = use -w)")
	cliflag.IntVar(fs, &right, "r", "right", -1, "Right extension bp (-1 = use -w)")
	fs.BoolVar(&sm, "sm", false, "Same-strand")
	fs.BoolVar(&sw, "sw", false, "Opposite-strand")
	cliflag.BoolVar(fs, &writeA, "wa", "write-a", false, "Write A only")
	cliflag.BoolVar(fs, &writeB, "wb", "write-b", false, "Write B only")
	cliflag.BoolVar(fs, &count, "c", "count", false, "Count hits per A")
	cliflag.BoolVar(fs, &invert, "v", "invert", false, "Report A with no overlap")
	cliflag.IntVar(fs, &minOL, "m", "min-overlap", 1, "Min overlap bp")
	cliflag.BoolVar(fs, &showHelp, "h", "help", false, "Show help")
	fs.BoolVar(&showVer, "version", false, "Show version")

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

	// -l / -r override -w on whichever side they are set.
	leftBP := window
	rightBP := window
	if left >= 0 {
		leftBP = left
	}
	if right >= 0 {
		rightBP = right
	}

	aR, err := iohelper.OpenReader(inputA)
	if err != nil {
		return err
	}
	defer aR.Close()
	bR, err := iohelper.OpenReader(inputB)
	if err != nil {
		return err
	}
	defer bR.Close()
	w, err := iohelper.OpenWriter(output)
	if err != nil {
		return err
	}
	defer w.Close()

	_, err = bedwindow.Window(aR, bR, w, bedwindow.Options{
		Left:          leftBP,
		Right:         rightBP,
		StrandSpec:    sm,
		InverseStrand: sw,
		WriteA:        writeA,
		WriteB:        writeB,
		Count:         count,
		Invert:        invert,
		MinOverlap:    minOL,
	})
	return err
}
