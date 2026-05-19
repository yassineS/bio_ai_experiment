// bednuc profiles the nucleotide content of each interval in a BED file
// against an indexed FASTA reference (Go port of `bedtools nuc`).
package main

import (
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bednuc/pkg/bednuc"
)

const version = "1.0.0"

const usage = `bednuc - Nucleotide content of BED intervals

Usage:
  bednuc -fi <fasta> -bed <bed> [options]

Options:
  -fi  FASTA           Indexed FASTA reference (required)
  -bed BED             BED/GFF/VCF intervals ('-' for stdin)
  -s,  --strand        Reverse-complement '-' strand intervals before counting
  -seq, --seq          Emit the extracted sequence as an extra column
  -pattern STR         Count occurrences of substring STR per interval
  -C,  --ignorecase    Ignore case when matching -pattern (upstream default
                       is case-sensitive)
  -fullHeader          Match contigs by their full FASTA header
  -o,  --output FILE   Output file (default: stdout)
  -h,  --help          Show this help message
  -v,  --version       Show version information

Output:
  A '#'-prefixed header line followed by each input record extended with:

    %AT  %GC  #A  #C  #G  #T  #N  #oth  seq_len  [seq]  [pattern_count]

  Columns enclosed in [] appear only when -seq / -pattern are set.

Notes:
  - Coordinates are 0-based, half-open.
  - Features whose chromosome is missing, or whose end exceeds the contig
    length, are skipped with a warning on stderr (parity with upstream).
`

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, "Error:", err)
		os.Exit(1)
	}
}

func run(argv []string, stdout, stderr *os.File) error {
	fs := flag.NewFlagSet("bednuc", flag.ContinueOnError)
	fs.SetOutput(stderr)

	var (
		fastaPath, bedPath, output string
		pattern                    string
		patternSet                 bool
		forceStrand, printSeq      bool
		ignoreCase, fullHeader     bool
		help, showVer              bool
	)

	cliflag.StringVar(fs, &fastaPath, "", "fi", "", "FASTA reference (required)")
	cliflag.StringVar(fs, &bedPath, "", "bed", "", "BED file (required, '-' for stdin)")
	cliflag.BoolVar(fs, &forceStrand, "s", "strand", false, "Reverse-complement '-' strand intervals")
	cliflag.BoolVar(fs, &printSeq, "", "seq", false, "Emit extracted sequence")
	// -seq is a single-dash multi-letter flag in upstream; register it
	// directly so `-seq` works (Go's flag package accepts `--seq` too).
	fs.BoolVar(&printSeq, "seq", false, "")
	cliflag.StringVar(fs, &pattern, "", "pattern", "", "Substring to count per interval")
	fs.StringVar(&pattern, "pattern", "", "")
	cliflag.BoolVar(fs, &ignoreCase, "C", "ignorecase", false, "Ignore case for -pattern")
	cliflag.BoolVar(fs, &fullHeader, "", "fullHeader", false, "Match contigs by full FASTA header")
	fs.BoolVar(&fullHeader, "fullHeader", false, "")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout)")
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVer, "v", "version", false, "Show version")

	fs.Usage = func() { fmt.Fprintf(stderr, "%s", usage) }

	// Detect whether -pattern was provided so we know to emit the column
	// even when the user passes an empty pattern explicitly. We treat
	// absence-of-flag as no-pattern and any explicit `-pattern X` as
	// pattern mode on.
	for _, a := range argv {
		if a == "-pattern" || a == "--pattern" ||
			strings.HasPrefix(a, "-pattern=") ||
			strings.HasPrefix(a, "--pattern=") {
			patternSet = true
			break
		}
	}

	if err := fs.Parse(argv); err != nil {
		return err
	}
	if help {
		fmt.Fprintf(stderr, "%s", usage)
		return nil
	}
	if showVer {
		fmt.Fprintf(stdout, "bednuc version %s\n", version)
		return nil
	}
	if fastaPath == "" || bedPath == "" {
		return fmt.Errorf("both -fi and -bed are required (use -h for help)")
	}
	if pattern != "" {
		patternSet = true
	}

	bedR, err := iohelper.OpenReader(bedPath)
	if err != nil {
		return fmt.Errorf("opening BED: %w", err)
	}
	defer bedR.Close()
	w, err := iohelper.OpenWriter(output)
	if err != nil {
		return fmt.Errorf("opening output: %w", err)
	}
	defer w.Close()

	opts := bednuc.Options{
		PrintSeq:    printSeq,
		Pattern:     pattern,
		HasPattern:  patternSet,
		IgnoreCase:  ignoreCase,
		ForceStrand: forceStrand,
		FullHeader:  fullHeader,
	}
	if _, err := bednuc.Run(bedR, fastaPath, w, stderr, opts); err != nil {
		return err
	}
	return nil
}
