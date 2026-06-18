// htsfile — print one-line summaries of bioinformatics file formats.
// Mirrors htslib's htsfile binary at a behavioural level. Detects
// SAM/BAM/CRAM/VCF/BCF/FASTA/FASTQ/BED/GFF wrapped in plain text,
// gzip, or BGZF.
//
// Usage:
//
//	htsfile FILE [FILE...]
//
// Each FILE is sniffed and a one-line summary is printed in the form
// "FILE: <format-name> [version <X.Y>] <compression> <kind>", e.g.:
//
//	$ htsfile alignments.bam variants.vcf.gz
//	alignments.bam: BAM BGZF-compressed sequence data
//	variants.vcf.gz: VCF version 4.2 BGZF-compressed variant calling data
//
// The "-c"/"--view" mode writes the canonical textual form of each file
// (the format-aware htslib round-trip): a VCF/BCF is re-serialised as VCF
// text. Other formats report that view is unsupported rather than emitting
// non-canonical bytes.
//
// Differences from upstream htsfile (intentional):
//   - We don't link against libhts; sniffing is a pure-Go peek that
//     never decompresses more than the first BGZF block.
//   - The "--copy" mode (`-C`) raw byte copy is not implemented.
//   - -h/-v stay bound to help/version (project CLI convention), so
//     upstream's -h/--header-only / -H/--no-header view modifiers are
//     not offered.

package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/htsfile/pkg/htsfile"
)

const usage = `htsfile - print bioinformatics file format information

Usage:
  htsfile FILE [FILE...]

Each FILE is identified and a one-line summary is printed. With "-"
as a filename the input is read from stdin (still sniffed, but the
output prefix becomes "<stdin>:").

Flags:
  -h, --help     show this help
  -v, --version  print version and exit
`

const version = "htsfile 0.1.0 (pure-Go, htsgo)"

func main() {
	os.Exit(run(os.Args[1:]))
}

func run(args []string) int {
	fs := flag.NewFlagSet("htsfile", flag.ContinueOnError)
	fs.SetOutput(io.Discard) // we print usage ourselves
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	var (
		helpFlag     bool
		versionFlag  bool
		copyMode     bool
		copyNoDecode bool
	)
	cliflag.BoolVar(fs, &helpFlag, "h", "help", false, "show help")
	cliflag.BoolVar(fs, &versionFlag, "v", "version", false, "show version")
	// -c ("--view") is the format-aware view: each file is routed through the
	// matching format reader and re-serialised (a VCF/BCF becomes canonical VCF
	// text), NOT raw-decompressed. -C ("--copy") is the raw byte copy. Our -h/-v
	// stay bound to help/version per the project CLI convention, so upstream's
	// -h/--header-only / -H/--no-header view modifiers are not offered.
	cliflag.BoolVar(fs, &copyMode, "c", "view", false, "Write the textual (format-aware) form of FILEs to stdout")
	cliflag.BoolVar(fs, &copyNoDecode, "C", "copy", false, "Ignored: raw copy mode not implemented (legacy)")

	// Route through cliflag.Parse so POSIX getopt-style short-flag bundling
	// works the way upstream htsfile's getopt parser accepts it.
	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, usage)
		return 2
	}

	if helpFlag {
		fmt.Print(usage)
		return 0
	}
	if versionFlag {
		fmt.Println(version)
		return 0
	}
	_ = copyNoDecode

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "htsfile: no input file specified")
		fmt.Fprint(os.Stderr, usage)
		return 1
	}

	failed := 0
	for _, path := range rest {
		var err error
		if copyMode {
			err = viewOne(path)
		} else {
			err = runOne(path)
		}
		if err != nil {
			fmt.Fprintf(os.Stderr, "htsfile: %s: %v\n", path, err)
			failed++
		}
	}
	if failed > 0 {
		return 1
	}
	return 0
}

// viewOne identifies path and writes its canonical textual form to stdout,
// implementing `htsfile -c`.
func viewOne(path string) error {
	var (
		f   *htsfile.Format
		err error
	)
	if path == "-" {
		// Identification consumes the sniff bytes from stdin; viewing stdin
		// would need to re-read them, which this tool does not buffer.
		return fmt.Errorf("htsfile -c does not support stdin")
	}
	f, err = htsfile.Identify(path)
	if err != nil {
		return err
	}
	return htsfile.View(path, f, os.Stdout)
}

func runOne(path string) error {
	var (
		f   *htsfile.Format
		err error
	)
	label := path
	if path == "-" {
		label = "<stdin>"
		f, err = htsfile.IdentifyReader(os.Stdin)
	} else {
		f, err = htsfile.Identify(path)
	}
	if err != nil {
		return err
	}
	// htslib's htsfile separates the path and description with a TAB (":\t"),
	// not a space.
	_, werr := fmt.Fprintf(os.Stdout, "%s:\t%s\n", label, f.Describe())
	if werr != nil && werr != io.EOF {
		return werr
	}
	return nil
}
