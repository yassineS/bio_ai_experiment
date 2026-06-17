// bedwindow examines a window around each feature in A and reports features in
// B that overlap that window (Go port of `bedtools window`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedwindow/pkg/bedwindow"
)

const usage = `bedwindow - Examine a window around each A feature and report overlapping B features

Usage:
  bedwindow -a A.bed -b B.bed [options]

Options:
  -a, --input-a FILE      Input BED file A (required)
  -b, --input-b FILE      Input BED file B (required)
  -o, --output FILE       Output file (default: stdout)
  -w, --window INT        Bp added upstream and downstream of each A entry (default 1000)
  -l, --left INT          Bp added upstream (left of) each A entry (default 1000)
  -r, --right INT         Bp added downstream (right of) each A entry (default 1000)
  -sw                     Define -l/-r based on A's strand (swap for '-' strand)
  -sm                     Only report B hits on the SAME strand as A
  -Sm                     Only report B hits on the OPPOSITE strand to A
  -u, --unique            Write each A entry once if it has any B overlap
  -c, --count             Append the B-hit count to each A entry (0 included)
  -v, --invert            Report only A entries with NO B overlap
  -wa, --write-a          Write the original A entry only
  -wb, --write-b          Write the original B entry only
  -h, --help              Show help
  --version               Show version

Default output:
  A<TAB>B for each (A, B) overlap pair (matches upstream).
`

const version = "0.2.0"

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
		sw       bool
		sm       bool
		bigSm    bool
		writeA   bool
		writeB   bool
		anyHit   bool
		count    bool
		invert   bool
		showHelp bool
		showVer  bool
	)
	// Sentinel -1 means "not set"; the default window of 1000 is applied below.
	const unset = -1
	cliflag.StringVar(fs, &inputA, "a", "input-a", "", "BED A")
	cliflag.StringVar(fs, &inputB, "b", "input-b", "", "BED B")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output")
	cliflag.IntVar(fs, &window, "w", "window", unset, "Window bp (default 1000)")
	cliflag.IntVar(fs, &left, "l", "left", unset, "Left extension bp (default 1000)")
	cliflag.IntVar(fs, &right, "r", "right", unset, "Right extension bp (default 1000)")
	fs.BoolVar(&sw, "sw", false, "Strand windows (define -l/-r by A strand)")
	fs.BoolVar(&sm, "sm", false, "Same-strand B hits only")
	fs.BoolVar(&bigSm, "Sm", false, "Opposite-strand B hits only")
	cliflag.BoolVar(fs, &writeA, "wa", "write-a", false, "Write A only")
	cliflag.BoolVar(fs, &writeB, "wb", "write-b", false, "Write B only")
	cliflag.BoolVar(fs, &anyHit, "u", "unique", false, "Write A once if any hit")
	cliflag.BoolVar(fs, &count, "c", "count", false, "Count hits per A")
	cliflag.BoolVar(fs, &invert, "v", "invert", false, "Report A with no overlap")
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

	// Validate the slop combination the way upstream does.
	haveW := window != unset
	haveL := left != unset
	haveR := right != unset
	if haveW && (haveL || haveR) {
		return fmt.Errorf("error: cannot combine -w with -l or -r")
	}
	if (haveL && !haveR) || (haveR && !haveL) {
		return fmt.Errorf("error: please specify both -l and -r")
	}
	if anyHit && invert {
		return fmt.Errorf("error: request either -u or -v, not both")
	}
	if anyHit && count {
		return fmt.Errorf("error: request either -u or -c, not both")
	}
	if sm && bigSm {
		return fmt.Errorf("error: use either -sm or -Sm, not both")
	}

	// Resolve the window. Upstream defaults both slops to 1000; -w sets both,
	// -l/-r set each side.
	leftBP, rightBP := 1000, 1000
	if haveW {
		leftBP, rightBP = window, window
	}
	if haveL {
		leftBP = left
	}
	if haveR {
		rightBP = right
	}
	if leftBP < 0 {
		return fmt.Errorf("error: upstream window (-l) must be positive")
	}
	if rightBP < 0 {
		return fmt.Errorf("error: downstream window (-r) must be positive")
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
		StrandWindows: sw,
		StrandSpec:    sm,
		InverseStrand: bigSm,
		WriteA:        writeA,
		WriteB:        writeB,
		AnyHit:        anyHit,
		Count:         count,
		Invert:        invert,
	})
	return err
}
