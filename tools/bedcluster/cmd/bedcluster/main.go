// bedcluster groups overlapping intervals into clusters and tags each input
// record with a cluster ID (Go port of `bedtools cluster`).
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bedcluster/pkg/bedcluster"
)

const usage = `bedcluster - Cluster overlapping intervals (port of bedtools cluster)

Usage:
  bedcluster -i FILE.bed [options]

Options:
  -i, --input FILE     Input BED file (required, '-' for stdin)
  -o, --output FILE    Output file (default: stdout)
  -d, --distance INT   Max gap (bp) between intervals to still cluster (default 0)
  -s, --strand         Cluster only intervals on the same strand
  -h, --help           Show help
  -v, --version        Show version and exit

Output:
  Each input record is emitted with an additional final column holding its
  1-based cluster ID. Records sharing a cluster ID are within MaxDistance
  bp of each other.

Notes:
  - With -d 0 (default), book-ended intervals (one ending exactly where the
    next starts) DO cluster together — matches upstream bedtools.
  - Cluster IDs are global / monotonic and never reset across chromosomes.
`

const version = "0.1.0"

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func run(args []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bedcluster", flag.ContinueOnError)
	fs.SetOutput(stderr)
	fs.Usage = func() { fmt.Fprint(stderr, usage) }

	var (
		input    string
		output   string
		distance int
		strand   bool
		showHelp bool
		showVer  bool
	)
	cliflag.StringVar(fs, &input, "i", "input", "", "Input BED file")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file")
	cliflag.IntVar(fs, &distance, "d", "distance", 0, "Max gap (bp)")
	cliflag.BoolVar(fs, &strand, "s", "strand", false, "Cluster only same-strand")
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
	if input == "" {
		fmt.Fprint(stderr, usage)
		return fmt.Errorf("error: -i/--input is required")
	}

	r, err := iohelper.OpenReader(input)
	if err != nil {
		return fmt.Errorf("opening input: %w", err)
	}
	defer r.Close()

	w, err := iohelper.OpenWriter(output)
	if err != nil {
		return fmt.Errorf("opening output: %w", err)
	}
	defer w.Close()

	_, err = bedcluster.Cluster(r, w, bedcluster.Options{MaxDistance: distance, StrandSpec: strand})
	return err
}
