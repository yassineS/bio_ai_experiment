// Command bedbamtobed is a pure-Go reimplementation of `bedtools bamtobed`
// (aka bamToBed). It converts BAM/SAM alignments to BED6, blocked BED12, or
// BEDPE. See pkg/bedbamtobed for the conversion logic and parity matrix.
package main

import (
	"bufio"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedbamtobed/pkg/bedbamtobed"
)

const version = "0.1.0"

const usage = `bedbamtobed - convert BAM alignments to BED6, BED12, or BEDPE.

Usage:
  bedbamtobed [OPTIONS] -i <bam>

Options:
  -i FILE     BAM/SAM input ('-' = stdin, the default).
  -bedpe      Write BEDPE format (requires BAM grouped/sorted by query).
  -mate1      With -bedpe, always report mate one as the first block.
  -bed12      Write "blocked" BED12 format. Forces -split.
  -split      Report "split" alignments as separate BED entries (splits on N).
  -splitD     Split alignments on N and D CIGAR ops. Forces -split.
  -ed         Use BAM edit distance (NM tag) as the BED score.
  -tag TAG    Use another NUMERIC BAM tag for the BED score (disallowed -bedpe).
  -color R,G,B  itemRgb for BED12 output. Default 255,0,0.
  -cigar      Append the CIGAR string as a trailing column.
  -h, --help     Show this help.
  -v, --version  Show version.
`

func main() { os.Exit(run(os.Args[1:], os.Stdout, os.Stderr)) }

func run(args []string, stdout, stderr io.Writer) int {
	var (
		bamFile = "-"
		opts    bedbamtobed.Options
	)
	opts.Color = "255,0,0"

	for i := 0; i < len(args); i++ {
		a := args[i]
		switch a {
		case "-h", "--help":
			fmt.Fprint(stderr, usage)
			return 0
		case "-v", "--version":
			fmt.Fprintf(stdout, "bedbamtobed %s\n", version)
			return 0
		case "-i":
			if i+1 < len(args) {
				bamFile = args[i+1]
				i++
			}
		case "-bedpe":
			opts.WriteBedPE = true
		case "-bed12":
			opts.WriteBed12 = true
		case "-split":
			opts.ObeySplits = true
		case "-splitD":
			opts.ObeySplits = true
			opts.SplitOnDeletions = true
		case "-ed", "-edit":
			opts.UseEditDistance = true
			opts.Tag = "NM"
		case "-cigar":
			opts.UseCigar = true
		case "-mate1":
			opts.Mate1First = true
		case "-color":
			if i+1 < len(args) {
				opts.Color = args[i+1]
				i++
			}
		case "-tag":
			if i+1 < len(args) {
				opts.Tag = args[i+1]
				i++
			}
		default:
			fmt.Fprintf(stderr, "\n*****ERROR: Unrecognized parameter: %s *****\n\n", a)
			fmt.Fprint(stderr, usage)
			return 1
		}
	}

	if err := opts.Validate(); err != nil {
		fmt.Fprintf(stderr, "\n*****\n*****ERROR: %s\n*****\n", err)
		return 1
	}

	in, err := iohelper.OpenReader(bamFile)
	if err != nil {
		fmt.Fprintf(stderr, "Failed to open BAM file %s\n", bamFile)
		return 1
	}
	defer in.Close()

	bw := bufio.NewWriter(stdout)
	defer bw.Flush()

	if _, err := bedbamtobed.Run(in, bw, opts); err != nil {
		bw.Flush()
		fmt.Fprintf(stderr, "%s\n", err)
		return 1
	}
	return 0
}
