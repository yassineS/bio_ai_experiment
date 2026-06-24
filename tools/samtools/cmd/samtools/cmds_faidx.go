package main

// CLI wiring for the faidx and fqidx subcommands. These build a FASTA/FASTQ
// index when called with no regions, or extract regions to FASTA/FASTQ when
// regions (or a -r region file) are supplied. The flag surface mirrors
// `samtools faidx --help` / `samtools fqidx --help`.

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/samtools/pkg/samtools"
)

// faidxUsage / fqidxUsage reproduce `samtools {faidx,fqidx} --help` byte for
// byte (including the "Option: " trailing space and the global-options
// footer). The footer's --output-fmt-option / -@/--threads / --write-index
// lines are listed for parity even though those flags are no-ops in this port.
const faidxUsage = `Usage: samtools faidx <file.fa|file.fa.gz> [<reg> [...]]
Option: 
  -o, --output FILE        Write FASTA to file.
  -n, --length INT         Length of FASTA sequence line. [60]
  -c, --continue           Continue after trying to retrieve missing region.
  -r, --region-file FILE   File of regions.  Format is chr:from-to. One per line.
  -i, --reverse-complement Reverse complement sequences.
      --mark-strand TYPE   Add strand indicator to sequence name
                           TYPE = rc   for /rc on negative strand (default)
                                  no   for no strand indicator
                                  sign for (+) / (-)
                                  custom,<pos>,<neg> for custom indicator
      --fai-idx      FILE  name of the index file (default file.fa.fai).
      --gzi-idx      FILE  name of compressed file index (default file.fa.gz.gzi).
  -f, --fastq              File and index in FASTQ format.
  -h, --help               This message.
      --output-fmt-option OPT[=VAL]
               Specify a single output file format option in the form
               of OPTION or OPTION=VALUE
  -@, --threads INT
               Number of additional threads to use [0]
      --write-index
               Automatically index the output files [off]

See https://www.htslib.org/doc/samtools.html#GLOBAL_COMMAND_OPTIONS
for more details.
`

const fqidxUsage = `Usage: samtools fqidx <file.fq|file.fq.gz> [<reg> [...]]
Option: 
  -o, --output FILE        Write FASTQ to file.
  -n, --length INT         Length of FASTQ sequence line. [60]
  -c, --continue           Continue after trying to retrieve missing region.
  -r, --region-file FILE   File of regions.  Format is chr:from-to. One per line.
  -i, --reverse-complement Reverse complement sequences.
      --mark-strand TYPE   Add strand indicator to sequence name
                           TYPE = rc   for /rc on negative strand (default)
                                  no   for no strand indicator
                                  sign for (+) / (-)
                                  custom,<pos>,<neg> for custom indicator
      --fai-idx      FILE  name of the index file (default file.fq.fai).
      --gzi-idx      FILE  name of compressed file index (default file.fq.gz.gzi).
  -h, --help               This message.
      --output-fmt-option OPT[=VAL]
               Specify a single output file format option in the form
               of OPTION or OPTION=VALUE
  -@, --threads INT
               Number of additional threads to use [0]
      --write-index
               Automatically index the output files [off]

See https://www.htslib.org/doc/samtools.html#GLOBAL_COMMAND_OPTIONS
for more details.
`

func runFaidx(args []string) int { return runFaidxCore(args, samtools.FaidxFASTA, faidxUsage) }

func runFqidx(args []string) int { return runFaidxCore(args, samtools.FaidxFASTQ, fqidxUsage) }

func runFaidxCore(args []string, format samtools.FaidxFormat, usage string) int {
	name := "faidx"
	if format == samtools.FaidxFASTQ {
		name = "fqidx"
	}
	fs := flag.NewFlagSet("samtools "+name, flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		output     string
		lineLenStr string
		cont       bool
		regionFile string
		rev        bool
		markStrand string
		faiName    string
		gziName    string
		asFastq    bool
		showHelp   bool
		showVer    bool
		threads    int // accepted, no-op (single-threaded)
	)
	// -n carries the upstream default sentinel of -60 ("same as input data").
	cliflag.StringVar(fs, &output, "o", "output", "", "Write to file.")
	cliflag.StringVar(fs, &lineLenStr, "n", "length", "-60", "Length of sequence line.")
	cliflag.BoolVar(fs, &cont, "c", "continue", false, "Continue past missing regions.")
	cliflag.StringVar(fs, &regionFile, "r", "region-file", "", "File of regions.")
	cliflag.BoolVar(fs, &rev, "i", "reverse-complement", false, "Reverse complement.")
	fs.StringVar(&markStrand, "mark-strand", "rc", "")
	fs.StringVar(&faiName, "fai-idx", "", "")
	fs.StringVar(&gziName, "gzi-idx", "", "")
	cliflag.BoolVar(fs, &asFastq, "f", "fastq", false, "FASTQ format.")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Additional threads (no-op).")

	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprint(os.Stderr, usage)
		return 1
	}
	_ = threads
	if showHelp {
		fmt.Print(usage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	// No positional argument: upstream prints usage to stdout and exits 0.
	if fs.NArg() == 0 {
		fmt.Print(usage)
		return 0
	}

	// -f promotes a faidx invocation to FASTQ format (only valid for faidx).
	if asFastq {
		format = samtools.FaidxFASTQ
	}

	opts := samtools.DefaultFaidxOptions(format)
	opts.Output = output
	opts.Continue = cont
	opts.RegionFile = regionFile
	opts.ReverseComplement = rev
	opts.FaiName = faiName
	opts.GziName = gziName

	// Resolve -n. Default sentinel -60 means "same as input"; a user-supplied
	// negative value warns and falls back to 60 (matching upstream).
	if lineLenStr != "-60" {
		v, perr := parseFaidxLineLen(lineLenStr)
		if perr != nil {
			fmt.Fprint(os.Stderr, usage)
			return 1
		}
		if v < 0 {
			fmt.Fprintf(os.Stderr, "[faidx] bad line length '%s', using default:%d\n", lineLenStr, 60)
			v = 60
		}
		opts.LineLen = v
	}

	if err := samtools.ParseMarkStrand(&opts, markStrand); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, usage)
		return 1
	}

	path := fs.Arg(0)
	regions := fs.Args()[1:]

	// Resolve the output sink. "-" or empty means stdout.
	out, oerr := openOut(output)
	if oerr != nil {
		fmt.Fprintf(os.Stderr, "[faidx] Cannot open \"%s\" for writing.\n", output)
		return 1
	}
	defer out.Close()

	return samtools.Faidx(path, regions, opts, out, os.Stderr)
}

// parseFaidxLineLen parses the -n value as a base-10 integer (strtol-style:
// leading sign permitted, trailing garbage ignored).
func parseFaidxLineLen(s string) (int, error) {
	i := 0
	neg := false
	if i < len(s) && (s[i] == '+' || s[i] == '-') {
		neg = s[i] == '-'
		i++
	}
	n := 0
	seen := false
	for i < len(s) && s[i] >= '0' && s[i] <= '9' {
		n = n*10 + int(s[i]-'0')
		i++
		seen = true
	}
	if !seen {
		return 0, fmt.Errorf("invalid line length %q", s)
	}
	if neg {
		n = -n
	}
	return n, nil
}
