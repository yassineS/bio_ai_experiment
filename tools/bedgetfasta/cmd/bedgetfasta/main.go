// Command bedgetfasta is a pure-Go reimplementation of `bedtools getfasta`.
// For each BED interval, it pulls the corresponding FASTA subsequence from
// a FAI-indexed reference and writes the result as FASTA (or TSV with
// `-tab`). See pkg/bedgetfasta for behaviour.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedgetfasta/pkg/bedgetfasta"
)

const version = "0.1.0"

const usage = `bedgetfasta - extract FASTA subsequences for BED intervals.

Usage:
  bedgetfasta -fi FASTA -bed BED [-fo OUT] [-name|-name+|-nameOnly] [-tab|-bedOut] [-s] [-split] [-rna]

Required:
  -fi, --fasta FILE       FASTA reference. Must have a sibling .fai index
                          (one will be built on the fly if missing).
  -bed, --bed FILE        BED file ('-' = stdin). Transparent gzip.

Optional:
  -fo, --output FILE      Output (default stdout). '-' = stdout.
  -name, --name           Header is '<name>::<chrom>:<start>-<end>'.
  -name+                  Deprecated alias of -name (identical header).
  -nameOnly, --nameOnly   Header is just '<name>'.
  -tab, --tab             Emit TSV ('header<TAB>seq') instead of FASTA.
  -bedOut, --bedOut       Re-emit the BED record with a trailing sequence
                          column (tab-delimited) instead of FASTA.
  -s, --strand            Reverse-complement '-' strand intervals
                          (case-preserving, IUPAC-aware).
  -split, --split         Concatenate the blocks of BED12 records before
                          emission (per-block stranded with -s).
  -rna, --rna             Emit U/u in place of T/t (after -s).
  -fullHeader, --full-header
                          Index FASTA contigs by the full header line
                          (whitespace included). Lets a BED row keyed by
                          a multi-word contig name match the corresponding
                          FASTA sequence. The sibling .fai (which only
                          stores first-token names) is rebuilt in-memory.

Standard:
  -h, --help              Show this help.
  -v, --version           Show version.
`

func main() {
	fs := flag.NewFlagSet("bedgetfasta", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		fastaPath   string
		bedPath     string
		outPath     string
		nameFlag    bool
		namePlus    bool
		nameOnly    bool
		tab         bool
		bedOut      bool
		strand      bool
		split       bool
		rna         bool
		fullHeader  bool
		showHelp    bool
		showVersion bool
	)
	// NOTE: cliflag registers a short name and a long name separately on the
	// FlagSet, so a name must not be supplied as BOTH short and long (that
	// would register the same flag twice and panic). For upstream-only flag
	// names (-bed, -name, -tab, ...) we register the single canonical name;
	// Go's flag package accepts it under either one or two leading dashes.
	cliflag.StringVar(fs, &fastaPath, "fi", "fasta", "", "FASTA reference")
	cliflag.StringVar(fs, &bedPath, "bed", "", "", "BED file")
	cliflag.StringVar(fs, &outPath, "fo", "output", "", "Output file")
	cliflag.BoolVar(fs, &nameFlag, "name", "", false, "Use BED name as header")
	// -name+ is upstream's deprecated alias of -name. The literal '+' is a
	// valid Go flag name, so register it directly on the FlagSet.
	fs.BoolVar(&namePlus, "name+", false, "(deprecated) same as -name")
	cliflag.BoolVar(fs, &nameOnly, "nameOnly", "", false, "Header is just <name>")
	cliflag.BoolVar(fs, &tab, "tab", "", false, "TSV output")
	cliflag.BoolVar(fs, &bedOut, "bedOut", "", false, "Emit BED record + trailing seq column")
	cliflag.BoolVar(fs, &strand, "s", "strand", false, "Reverse-complement '-' strand")
	cliflag.BoolVar(fs, &split, "split", "", false, "Concatenate BED12 blocks")
	cliflag.BoolVar(fs, &rna, "rna", "", false, "Emit U in place of T")
	cliflag.BoolVar(fs, &fullHeader, "fullHeader", "full-header", false, "Index by full FASTA header line")
	cliflag.BoolVar(fs, &showHelp, "h", "help", false, "Help")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Version")

	if err := fs.Parse(os.Args[1:]); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, usage)
		os.Exit(2)
	}
	if showHelp {
		fmt.Print(usage)
		return
	}
	if showVersion {
		fmt.Println(version)
		return
	}
	if fastaPath == "" {
		fmt.Fprintln(os.Stderr, "bedgetfasta: -fi FASTA is required")
		fmt.Fprint(os.Stderr, usage)
		os.Exit(2)
	}
	if bedPath == "" {
		fmt.Fprintln(os.Stderr, "bedgetfasta: -bed BED is required")
		fmt.Fprint(os.Stderr, usage)
		os.Exit(2)
	}

	bed, err := iohelper.OpenReader(bedPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bedgetfasta: open bed: %v\n", err)
		os.Exit(1)
	}
	defer bed.Close()

	out, err := iohelper.OpenWriter(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bedgetfasta: open output: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if _, err := bedgetfasta.Run(bed, fastaPath, out, os.Stderr, bedgetfasta.Options{
		Name:       nameFlag,
		NamePlus:   namePlus,
		NameOnly:   nameOnly,
		Tab:        tab,
		BedOut:     bedOut,
		Strand:     strand,
		Split:      split,
		RNA:        rna,
		FullHeader: fullHeader,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "bedgetfasta: %v\n", err)
		os.Exit(1)
	}
}
