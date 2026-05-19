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
	flag.Usage = func() { fmt.Fprint(os.Stderr, usage) }
	helpFlag := flag.Bool("help", false, "")
	helpShort := flag.Bool("h", false, "")
	versionFlag := flag.Bool("version", false, "")
	versionShort := flag.Bool("v", false, "")
	flag.Parse()

	if *helpFlag || *helpShort {
		fmt.Print(usage)
		return
	}
	if *versionFlag || *versionShort {
		fmt.Println(version)
		return
	}

	args := flag.Args()
	if len(args) == 0 {
		fmt.Fprintln(os.Stderr, "htsfile: no input file specified")
		fmt.Fprint(os.Stderr, usage)
		os.Exit(1)
	}

	failed := 0
	for _, path := range args {
		if err := runOne(path); err != nil {
			fmt.Fprintf(os.Stderr, "htsfile: %s: %v\n", path, err)
			failed++
		}
	}
	if failed > 0 {
		os.Exit(1)
	}
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
