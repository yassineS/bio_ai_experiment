// Command bedigv is a pure-Go reimplementation of `bedtools igv`. It reads a
// BED file and emits an IGV batch-mode script that takes one snapshot per
// interval. See pkg/bedigv for behaviour and the parity matrix.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedigv/pkg/bedigv"
)

const version = "0.1.0"

const usage = `bedigv - emit an IGV batch script from BED intervals.

Usage:
  bedigv [OPTIONS] -i <bed>

I/O:
  -i, --input  FILE     BED input ('-' or empty = stdin). Transparent gzip.
  -o, --output FILE     Output file (default stdout). '-' = stdout.

IGV options (match upstream 'bedtools igv'):
  -path PATH            Directory IGV writes snapshots into. Default './'.
  -sess FILE            Existing IGV session to 'load' before snapshots.
  -sort TYPE            Sort BAM reads per snapshot. One of:
                          base, position, strand, quality, sample, readGroup.
                        Default: no sort.
  -clps                 Collapse aligned reads before snapshot.
  -name                 Append the BED record's name column to each
                        snapshot's filename. Errors when name is empty.
  -slop N               Flank each interval by N bp in the 'goto' locus.
                        Default 0; must be >= 0. The snapshot filename
                        keeps the original coordinates.
  -img TYPE             Snapshot extension. One of: png, eps, svg, jpg.
                        Default: png.

Standard:
  -h, --help            Show this help.
  -v, --version         Show version.
`

func main() {
	os.Exit(run(os.Args[1:], os.Stdout, os.Stderr))
}

func run(args []string, stdout, stderr io.Writer) int {
	fs := flag.NewFlagSet("bedigv", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		inFile, outFile       string
		path, sess, sortType  string
		imgType               string
		clps, useNames        bool
		slop                  int
		showHelp, showVersion bool
	)
	cliflag.StringVar(fs, &inFile, "i", "input", "", "Input BED")
	cliflag.StringVar(fs, &outFile, "o", "output", "", "Output file")
	fs.StringVar(&path, "path", "./", "Snapshot directory")
	fs.StringVar(&sess, "sess", "", "IGV session file to load")
	fs.StringVar(&sortType, "sort", "", "BAM sort directive")
	fs.BoolVar(&clps, "clps", false, "Emit a `collapse` line per record")
	fs.BoolVar(&useNames, "name", false, "Append BED name column to filenames")
	fs.IntVar(&slop, "slop", 0, "Flank each interval (bp) in the goto locus")
	fs.StringVar(&imgType, "img", "png", "Snapshot extension (png|eps|svg|jpg)")
	cliflag.BoolVar(fs, &showHelp, "h", "help", false, "Help")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Version")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(stderr, err)
		fmt.Fprint(stderr, usage)
		return 2
	}
	if showHelp {
		fmt.Fprint(stdout, usage)
		return 0
	}
	if showVersion {
		fmt.Fprintln(stdout, version)
		return 0
	}

	opts := bedigv.Options{
		Path:      path,
		Session:   sess,
		Sort:      bedigv.SortType(sortType),
		Collapse:  clps,
		UseNames:  useNames,
		Slop:      slop,
		ImageType: bedigv.ImageType(imgType),
	}

	switch opts.ImageType {
	case bedigv.ImagePNG, bedigv.ImageEPS, bedigv.ImageSVG, bedigv.ImageJPG:
	default:
		fmt.Fprintf(stderr, "bedigv: invalid -img %q (must be png|eps|svg|jpg)\n", imgType)
		return 2
	}
	if opts.Sort != bedigv.SortNone && !bedigv.IsValidSort(opts.Sort) {
		fmt.Fprintf(stderr, "bedigv: invalid -sort %q\n", sortType)
		return 2
	}
	if opts.Slop < 0 {
		fmt.Fprintf(stderr, "bedigv: -slop must be >= 0 (got %d)\n", opts.Slop)
		return 2
	}

	in, err := iohelper.OpenReader(inFile)
	if err != nil {
		fmt.Fprintf(stderr, "bedigv: open input: %v\n", err)
		return 1
	}
	defer in.Close()

	out, err := iohelper.OpenWriter(outFile)
	if err != nil {
		fmt.Fprintf(stderr, "bedigv: open output: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bedigv.Run(in, out, opts); err != nil {
		fmt.Fprintf(stderr, "bedigv: %v\n", err)
		return 1
	}
	return 0
}
