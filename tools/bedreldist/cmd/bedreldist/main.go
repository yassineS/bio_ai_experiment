// bedreldist computes the distribution of relative distances between
// intervals in two BED files (Go port of `bedtools reldist`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedreldist/pkg/bedreldist"
)

const version = "1.0.0"

const usage = `bedreldist - Relative-distance distribution between two BED files

Usage:
  bedreldist -a <fileA.bed> -b <fileB.bed> [options]

Description:
  For each interval in A, find the two nearest B-interval midpoints that
  flank A's midpoint on the same chromosome and compute the relative
  distance (in [0, 0.5]). The default output is a tab-separated histogram
  with 0.01-wide bins:

      reldist  count  total  fraction

  With -detail / --detail, each A interval is emitted on its own line
  followed by the computed relative distance.

Options:
  -a, --a FILE         BED file A (queries; required, '-' for stdin)
  -b, --b FILE         BED file B (database; required)
      --output FILE    Output file (default: stdout)
  -detail, --detail    Emit per-A relative distances instead of the histogram
  -h, --help           Show this help message
  -v, --version        Show version information

Notes:
  - Coordinates are 0-based, half-open.
  - Intervals whose chromosome is absent from B, or that fall before the
    first B-midpoint, do not contribute (matches upstream).
  - Histogram bins are emitted in ascending order, omitting empty bins.
`

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, "Error:", err)
		os.Exit(1)
	}
}

func run(argv []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bedreldist", flag.ContinueOnError)
	fs.SetOutput(stderr)

	var (
		fileA, fileB, output string
		detail               bool
		help, showVer        bool
	)

	cliflag.StringVar(fs, &fileA, "a", "", "", "BED file A (required)")
	cliflag.StringVar(fs, &fileB, "b", "", "", "BED file B (required)")
	cliflag.StringVar(fs, &output, "", "output", "", "Output file (default: stdout)")
	// Upstream uses a single-dash multi-letter flag `-detail`; Go's flag
	// package accepts both `-detail` and `--detail` for any flag, so
	// registering it under the long name covers both spellings.
	cliflag.BoolVar(fs, &detail, "", "detail", false, "Per-A detail rows")
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVer, "v", "version", false, "Show version")

	fs.Usage = func() { fmt.Fprint(stderr, usage) }

	if err := fs.Parse(argv); err != nil {
		return err
	}
	if help {
		fmt.Fprint(stderr, usage)
		return nil
	}
	if showVer {
		fmt.Fprintf(stdout, "bedreldist version %s\n", version)
		return nil
	}

	if fileA == "" || fileB == "" {
		return fmt.Errorf("both -a/--a and -b/--b are required (use -h for help)")
	}

	rA, err := iohelper.OpenReader(fileA)
	if err != nil {
		return fmt.Errorf("opening A: %w", err)
	}
	defer rA.Close()
	rB, err := iohelper.OpenReader(fileB)
	if err != nil {
		return fmt.Errorf("opening B: %w", err)
	}
	defer rB.Close()
	w, err := iohelper.OpenWriter(output)
	if err != nil {
		return fmt.Errorf("opening output: %w", err)
	}
	defer w.Close()

	if _, err := bedreldist.Run(rA, rB, w, bedreldist.Options{Detail: detail}); err != nil {
		return err
	}
	return nil
}
