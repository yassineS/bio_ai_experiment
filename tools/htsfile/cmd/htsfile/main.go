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
// Differences from upstream htsfile (intentional):
//   - We don't link against libhts; sniffing is a pure-Go peek that
//     never decompresses more than the first BGZF block.
//   - The "--copy" output mode (`-c`) is not implemented; the v1
//     scope is identification only. Use shell redirection if you
//     need the raw bytes.

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
	// Upstream htsfile (htsfile.c getopt "cChHv") accepts -c/-C copy modes.
	// The v1 scope is identification only, so these are accepted no-ops kept
	// for backward compatibility — and so bundled clusters that include them
	// still parse. The -H spelling is upstream's "print header" copy variant;
	// we map it to the same accepted no-op (our -h remains help, registered
	// above as the long --help short form).
	cliflag.BoolVar(fs, &copyMode, "c", "", false, "Ignored: copy mode not implemented (legacy)")
	cliflag.BoolVar(fs, &copyNoDecode, "C", "", false, "Ignored: copy mode not implemented (legacy)")

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
	_ = copyMode
	_ = copyNoDecode

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "htsfile: no input file specified")
		fmt.Fprint(os.Stderr, usage)
		return 1
	}

	failed := 0
	for _, path := range rest {
		if err := runOne(path); err != nil {
			fmt.Fprintf(os.Stderr, "htsfile: %s: %v\n", path, err)
			failed++
		}
	}
	if failed > 0 {
		return 1
	}
	return 0
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
	_, werr := fmt.Fprintf(os.Stdout, "%s: %s\n", label, f.Describe())
	if werr != nil && werr != io.EOF {
		return werr
	}
	return nil
}
