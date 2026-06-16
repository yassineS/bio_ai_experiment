// Command bedtobam is a pure-Go reimplementation of `bedtools bedtobam` (aka
// bedToBam). It converts BED (or BED12) records into BAM alignments against a
// header derived from a genome (chrom sizes) file. See pkg/bedtobam for the
// conversion logic and parity matrix.
package main

import (
	"bufio"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedtobam/pkg/bedtobam"
)

const version = "0.1.0"

const usage = `bedtobam - convert BED/BED12 features to BAM.

Usage:
  bedtobam [OPTIONS] -i <bed> -g <genome>

Options:
  -i FILE     BED/GFF/VCF input ('-' = stdin, the default).
  -g FILE     Genome (chrom sizes) file. Required.
  -mapq INT   MAPQ for the BAM records. Default 255.
  -bed12      Treat the input as BED12; the CIGAR reflects BED "blocks".
  -ubam       Write uncompressed BAM (default writes compressed BAM).
  -h, --help     Show this help.
  -v, --version  Show version.

Notes:
  (1) BED files must be at least BED4 to create BAM (needs a name field).
`

func main() { os.Exit(run(os.Args[1:], os.Stdout, os.Stderr)) }

func run(args []string, stdout, stderr io.Writer) int {
	var (
		bedFile    = "-"
		genomeFile string
		opts       = bedtobam.Options{MapQ: 255}
	)

	for i := 0; i < len(args); i++ {
		a := args[i]
		switch a {
		case "-h", "--help":
			fmt.Fprint(stderr, usage)
			return 0
		case "-v", "--version":
			fmt.Fprintf(stdout, "bedtobam %s\n", version)
			return 0
		case "-i":
			if i+1 < len(args) {
				bedFile = args[i+1]
				i++
			}
		case "-g":
			if i+1 < len(args) {
				genomeFile = args[i+1]
				i++
			}
		case "-mapq":
			if i+1 < len(args) {
				v, err := parseInt(args[i+1])
				if err != nil {
					fmt.Fprintf(stderr, "\n*****\n*****ERROR: bad -mapq value %q.\n*****\n", args[i+1])
					return 1
				}
				opts.MapQ = v
				i++
			}
		case "-bed12":
			opts.BED12 = true
		case "-ubam":
			opts.Uncompressed = true
		default:
			fmt.Fprintf(stderr, "\n*****ERROR: Unrecognized parameter: %s *****\n\n", a)
			fmt.Fprint(stderr, usage)
			return 1
		}
	}

	if genomeFile == "" {
		fmt.Fprint(stderr, "\n*****\n*****ERROR: Need -g (genome) file. \n*****\n")
		return 1
	}
	if opts.MapQ < 0 || opts.MapQ > 255 {
		fmt.Fprint(stderr, "\n*****\n*****ERROR: MAPQ must be in range [0,255]. \n*****\n")
		return 1
	}
	opts.GenomeFileName = genomeFile

	gf, err := iohelper.OpenReader(genomeFile)
	if err != nil {
		fmt.Fprintf(stderr, "Error: Unable to open genome file %s.  Exiting!\n", genomeFile)
		return 1
	}
	genome, err := bedtobam.ReadGenome(gf)
	gf.Close()
	if err != nil {
		fmt.Fprintf(stderr, "%s\n", err)
		return 1
	}

	in, err := iohelper.OpenReader(bedFile)
	if err != nil {
		fmt.Fprintf(stderr, "Error: Unable to open BED file %s.  Exiting!\n", bedFile)
		return 1
	}
	defer in.Close()

	bw := bufio.NewWriter(stdout)
	if _, err := bedtobam.Run(in, bw, genome, opts); err != nil {
		bw.Flush()
		fmt.Fprintf(stderr, "%s\n", err)
		return 1
	}
	return bwFlush(bw, stderr)
}

// bwFlush flushes the buffered writer, reporting any error as a non-zero exit.
func bwFlush(bw *bufio.Writer, stderr io.Writer) int {
	if err := bw.Flush(); err != nil {
		fmt.Fprintf(stderr, "%s\n", err)
		return 1
	}
	return 0
}

// parseInt parses a base-10 integer, matching upstream's atoi-style handling.
func parseInt(s string) (int, error) {
	var n int
	_, err := fmt.Sscanf(s, "%d", &n)
	return n, err
}
