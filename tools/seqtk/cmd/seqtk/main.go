// seqtk: A fast FASTA/Q processor implemented in Go
//
// This is a Go reimplementation of the popular seqtk tool with enhanced performance
// and better error handling. It provides common operations on FASTA and FASTQ files.
//
// Usage:
//
//	seqtk <command> [options] <input>
//
// Commands:
//
//	seq        Transform sequences (reverse complement, etc.)
//	subseq     Extract subsequences
//	sample     Subsample sequences
//	trimfq     Trim FASTQ sequences based on quality
//	fq2fa      Convert FASTQ to FASTA
//	comp       Get sequence composition statistics
//	mergepe    Interleave two paired-end FASTA/FASTQ files
//	cutN       Cut sequences at runs of N
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/seqtk/pkg/seqtk"
)

const version = "1.0.0"

func main() {
	if len(os.Args) < 2 {
		printUsage()
		os.Exit(1)
	}

	command := os.Args[1]

	switch command {
	case "seq":
		seqCommand()
	case "subseq":
		subseqCommand()
	case "sample":
		sampleCommand()
	case "trimfq":
		trimfqCommand()
	case "fq2fa":
		fq2faCommand()
	case "comp":
		compCommand()
	case "mergepe":
		mergePECommand()
	case "cutN":
		cutNCommand()
	case "version", "-v", "--version":
		fmt.Printf("seqtk version %s\n", version)
	case "help", "-h", "--help":
		printUsage()
	default:
		fmt.Fprintf(os.Stderr, "Unknown command: %s\n", command)
		printUsage()
		os.Exit(1)
	}
}

func printUsage() {
	fmt.Fprintf(os.Stderr, `seqtk - A fast FASTA/Q processor (Go implementation)

Usage: seqtk <command> [options] <input>

Commands:
  seq        Transform sequences (reverse complement, filtering)
  subseq     Extract subsequences
  sample     Subsample sequences
  trimfq     Trim FASTQ sequences based on quality
  fq2fa      Convert FASTQ to FASTA
  comp       Get sequence composition statistics
  mergepe    Interleave two paired-end FASTA/FASTQ files
  cutN       Cut sequences at runs of N (>= -n bases)
  version    Show version information
  help       Show this help message

Use 'seqtk <command> -h' for help on a specific command.

Examples:
  seqtk comp sequences.fasta
  seqtk fq2fa reads.fastq > reads.fasta
  seqtk seq -r sequences.fasta > rev_comp.fasta
  seqtk seq -l 100 -L 500 reads.fastq > filtered.fastq
  seqtk subseq genome.fa regions.bed > regions.fa
  seqtk subseq genome.fa names.txt > selected.fa
  seqtk sample reads.fastq 0.1 > sample.fastq
  seqtk trimfq -q 20 reads.fastq > trimmed.fastq
  seqtk mergepe r1.fq r2.fq > interleaved.fq
  seqtk cutN -n 10 genome.fa > fragments.fa

`)
}

func seqCommand() {
	fs := flag.NewFlagSet("seq", flag.ExitOnError)
	var revComp, phred64 bool
	var output string
	var minLen, maxLen int
	var pattern string

	cliflag.BoolVar(fs, &revComp, "r", "reverse", false, "Reverse complement sequences")
	cliflag.BoolVar(fs, &phred64, "6", "phred64", false, "Use Phred+64 quality encoding (default: Phred+33)")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout)")
	cliflag.IntVar(fs, &minLen, "l", "min-len", 0, "Minimum sequence length (0 = no filter)")
	cliflag.IntVar(fs, &maxLen, "L", "max-len", 0, "Maximum sequence length (0 = no filter)")
	cliflag.StringVar(fs, &pattern, "n", "name", "", "Filter by sequence name pattern")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk seq [options] <input>

Transform sequences.

Arguments:
  <input>    Input file (use '-' for stdin, supports .gz and .bz2)

Options:
  -r, --reverse          Reverse complement sequences
  -l, --min-len INT      Minimum sequence length (0 = no filter)
  -L, --max-len INT      Maximum sequence length (0 = no filter)
  -n, --name PATTERN     Filter by sequence name pattern
  -6, --phred64          Use Phred+64 quality encoding (default: Phred+33)
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk seq -r reads.fastq > rev_comp.fastq
  seqtk seq -l 100 -L 500 reads.fastq > filtered.fastq
  seqtk seq -n "chr1" sequences.fasta > chr1.fasta

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 1 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)

	// Detect file type (before opening to avoid consuming stdin)
	var isFastq bool
	var err error
	if inputFile != "-" {
		isFastq, err = seqtk.GetFileType(inputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error detecting file type: %v\n", err)
			os.Exit(1)
		}
	} else {
		// For stdin, assume FASTQ (could be improved with buffered detection)
		isFastq = true
	}

	// Open input with compression support
	input, err := seqtk.OpenInput(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	// Open output with compression support
	out, err := seqtk.OpenOutput(output)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	encoding := fastq.Phred33
	if phred64 {
		encoding = fastq.Phred64
	}

	// Check if any transformation or filter is specified
	hasFilter := minLen > 0 || maxLen > 0 || pattern != ""

	if revComp {
		if err := seqtk.ReverseComplement(input, out, isFastq, encoding); err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v\n", err)
			os.Exit(1)
		}
	} else if hasFilter {
		opts := seqtk.FilterOptions{
			MinLength: minLen,
			MaxLength: maxLen,
			Pattern:   pattern,
		}
		if err := seqtk.Filter(input, out, opts, isFastq, encoding); err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v\n", err)
			os.Exit(1)
		}
	} else {
		// Just copy through (could add more transformations here)
		fmt.Fprintf(os.Stderr, "No transformation specified. Use -r for reverse complement or -l/-L/-n for filtering.\n")
		os.Exit(1)
	}
}

func subseqCommand() {
	fs := flag.NewFlagSet("subseq", flag.ExitOnError)
	var output string
	var lineLen int

	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout)")
	cliflag.IntVar(fs, &lineLen, "l", "line-length", 0, "Wrap output sequence lines at INT characters (0 = no wrap)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk subseq [options] <in.fa> <name.list | reg.bed>

Extract subsequences from a FASTA/FASTQ file given a list of sequence names or a
BED file of regions. The second argument's format is auto-detected: if its first
non-comment line splits into at least three whitespace/tab fields whose second
and third fields are integers it is treated as BED, otherwise as a name list.
Output is always FASTA.

Arguments:
  <in.fa>    Input FASTA/FASTQ file (use '-' for stdin, supports .gz and .bz2)
  <name.list | reg.bed>
             Either one sequence name per line, or a BED file
             (chrom<TAB>start<TAB>end; 0-based half-open; extra columns ignored).
             For each BED region a record named "chrom:start+1-end" is emitted.

Options:
  -l, --line-length INT  Wrap output sequence lines at INT characters (0 = no wrap)
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk subseq genome.fa names.txt > selected.fa
  seqtk subseq genome.fa regions.bed > regions.fa
  cat reads.fq.gz | seqtk subseq - names.txt > selected.fa

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 2 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)
	specFile := fs.Arg(1)

	// Open input with compression support.
	input, err := seqtk.OpenInput(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	// Open the name list / BED file with compression support.
	spec, err := seqtk.OpenInput(specFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening region/name file: %v\n", err)
		os.Exit(1)
	}
	defer spec.Close()

	// Open output with compression support.
	out, err := seqtk.OpenOutput(output)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if err := seqtk.Subseq(input, spec, out, lineLen); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func sampleCommand() {
	fs := flag.NewFlagSet("sample", flag.ExitOnError)
	var phred64 bool
	var output string

	cliflag.BoolVar(fs, &phred64, "6", "phred64", false, "Use Phred+64 quality encoding (default: Phred+33)")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk sample [options] <input> <fraction>

Subsample sequences randomly.

Arguments:
  <input>      Input FASTA/FASTQ file (use '-' for stdin, supports .gz and .bz2)
  <fraction>   Fraction of sequences to sample (0.0-1.0)

Options:
  -6, --phred64          Use Phred+64 quality encoding (default: Phred+33)
  -o, --output FILE      Output file (default: stdout, supports .gz)

Example:
  seqtk sample reads.fastq 0.1 > sample.fastq
  cat reads.fastq.gz | seqtk sample - 0.1 > sample.fastq

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 2 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)
	var fraction float64
	if _, err := fmt.Sscanf(fs.Arg(1), "%f", &fraction); err != nil {
		fmt.Fprintf(os.Stderr, "Invalid fraction: %v\n", err)
		os.Exit(1)
	}

	// Detect file type
	var isFastq bool
	var err error
	if inputFile != "-" {
		isFastq, err = seqtk.GetFileType(inputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error detecting file type: %v\n", err)
			os.Exit(1)
		}
	} else {
		isFastq = true // Default to FASTQ for stdin
	}

	// Open input with compression support
	input, err := seqtk.OpenInput(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	// Open output with compression support
	out, err := seqtk.OpenOutput(output)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	encoding := fastq.Phred33
	if phred64 {
		encoding = fastq.Phred64
	}

	if err := seqtk.Sample(input, out, fraction, isFastq, encoding); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func trimfqCommand() {
	fs := flag.NewFlagSet("trimfq", flag.ExitOnError)
	var quality int
	var phred64 bool
	var output string

	cliflag.IntVar(fs, &quality, "q", "quality", 20, "Minimum quality threshold")
	cliflag.BoolVar(fs, &phred64, "6", "phred64", false, "Use Phred+64 quality encoding (default: Phred+33)")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk trimfq [options] <input>

Trim FASTQ sequences based on quality.

Arguments:
  <input>    Input FASTQ file (use '-' for stdin, supports .gz and .bz2)

Options:
  -q, --quality INT      Minimum quality threshold (default: 20)
  -6, --phred64          Use Phred+64 quality encoding (default: Phred+33)
  -o, --output FILE      Output file (default: stdout, supports .gz)

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 1 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)

	// Open input with compression support
	input, err := seqtk.OpenInput(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	// Open output with compression support
	out, err := seqtk.OpenOutput(output)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	encoding := fastq.Phred33
	if phred64 {
		encoding = fastq.Phred64
	}

	if err := seqtk.TrimQuality(input, out, quality, encoding); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func fq2faCommand() {
	fs := flag.NewFlagSet("fq2fa", flag.ExitOnError)
	var phred64 bool
	var output string

	cliflag.BoolVar(fs, &phred64, "6", "phred64", false, "Use Phred+64 quality encoding (default: Phred+33)")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk fq2fa [options] <input>

Convert FASTQ to FASTA.

Arguments:
  <input>    Input FASTQ file (use '-' for stdin, supports .gz and .bz2)

Options:
  -6, --phred64          Use Phred+64 quality encoding (default: Phred+33)
  -o, --output FILE      Output file (default: stdout, supports .gz)

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 1 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)

	// Open input with compression support
	input, err := seqtk.OpenInput(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	// Open output with compression support
	out, err := seqtk.OpenOutput(output)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	encoding := fastq.Phred33
	if phred64 {
		encoding = fastq.Phred64
	}

	if err := seqtk.ConvertFastqToFasta(input, out, encoding); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func mergePECommand() {
	fs := flag.NewFlagSet("mergepe", flag.ExitOnError)
	var output string

	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout, supports .gz)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk mergepe [options] <in1> <in2>

Interleave two paired-end FASTA/FASTQ files, writing
read1[0], read2[0], read1[1], read2[1], ... to the output stream.

Both inputs must have the same format and the same number of records;
if the counts differ, an error identifying the shorter input and the
pair index where the mismatch was detected is returned.

Arguments:
  <in1>      First mate file (use '-' for stdin, supports .gz)
  <in2>      Second mate file (use '-' for stdin, supports .gz)

Options:
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk mergepe r1.fq r2.fq > interleaved.fq
  seqtk mergepe r1.fa.gz r2.fa.gz -o pairs.fa.gz
  zcat r1.fq.gz | seqtk mergepe - r2.fq > interleaved.fq

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 2 {
		fs.Usage()
		os.Exit(1)
	}

	in1Name := fs.Arg(0)
	in2Name := fs.Arg(1)

	if in1Name == "-" && in2Name == "-" {
		fmt.Fprintf(os.Stderr, "Error: both inputs cannot be stdin ('-')\n")
		os.Exit(1)
	}

	in1, err := seqtk.OpenInput(in1Name)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening %s: %v\n", in1Name, err)
		os.Exit(1)
	}
	defer in1.Close()

	in2, err := seqtk.OpenInput(in2Name)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening %s: %v\n", in2Name, err)
		os.Exit(1)
	}
	defer in2.Close()

	out, err := seqtk.OpenOutput(output)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if err := seqtk.MergePE(in1, in2, out); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func cutNCommand() {
	fs := flag.NewFlagSet("cutN", flag.ExitOnError)
	var minN int
	var gaps bool
	var output string

	// -1 sentinel so we can detect "not provided" — required flag with no default.
	cliflag.IntVar(fs, &minN, "n", "min-n", -1, "Minimum N-run length to cut at (required)")
	cliflag.BoolVar(fs, &gaps, "g", "gaps", false, "Print BED-format records of cut N-runs to stderr")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout, supports .gz)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk cutN -n INT [options] <input>

Cut sequences at runs of 'N' or 'n' of length >= -n, writing the
resulting fragments as new FASTA records named "<orig-name>:<start>-<end>"
where coordinates are 1-based inclusive (start = position of first
retained base, end = position of last retained base).

Records with no qualifying N-run are emitted unchanged with their
original name. All-N sequences (or those with only leading/trailing
N-runs) produce no output for that record.

Output is always FASTA. Input may be FASTA or FASTQ (auto-detected).

Arguments:
  <input>    Input FASTA/FASTQ file (use '-' for stdin, supports .gz)

Options:
  -n, --min-n INT        Minimum N-run length to cut at (required)
  -g, --gaps             Print cut N-runs to stderr in BED format
                         (chrom<TAB>start0<TAB>end<TAB>N; 0-based half-open)
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk cutN -n 10 genome.fa > fragments.fa
  seqtk cutN -n 5 -g genome.fa > fragments.fa 2> gaps.bed

`)
	}

	fs.Parse(os.Args[2:])

	if minN < 1 {
		fmt.Fprintf(os.Stderr, "Error: -n/--min-n is required and must be >= 1\n")
		fs.Usage()
		os.Exit(1)
	}

	if fs.NArg() < 1 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)

	input, err := seqtk.OpenInput(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	out, err := seqtk.OpenOutput(output)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	opts := seqtk.CutNOptions{MinN: minN, EmitGaps: gaps}
	if gaps {
		opts.GapsW = os.Stderr
	}

	if err := seqtk.CutN(input, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func compCommand() {
	fs := flag.NewFlagSet("comp", flag.ExitOnError)
	var phred64 bool

	cliflag.BoolVar(fs, &phred64, "6", "phred64", false, "Use Phred+64 quality encoding for FASTQ (default: Phred+33)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk comp [options] <input>

Get sequence composition statistics.

Arguments:
  <input>    Input file (use '-' for stdin, supports .gz and .bz2)

Options:
  -6, --phred64          Use Phred+64 quality encoding for FASTQ (default: Phred+33)

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 1 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)

	// Detect file type
	var isFastq bool
	var err error
	if inputFile != "-" {
		isFastq, err = seqtk.GetFileType(inputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error detecting file type: %v\n", err)
			os.Exit(1)
		}
	} else {
		isFastq = true // Default to FASTQ for stdin
	}

	// Open input with compression support
	input, err := seqtk.OpenInput(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	var stats *seqtk.Stats
	if isFastq {
		encoding := fastq.Phred33
		if phred64 {
			encoding = fastq.Phred64
		}
		stats, err = seqtk.CalculateFastqStats(input, encoding)
	} else {
		stats, err = seqtk.CalculateFastaStats(input)
	}

	if err != nil {
		fmt.Fprintf(os.Stderr, "Error calculating statistics: %v\n", err)
		os.Exit(1)
	}

	fmt.Printf("Number of sequences: %d\n", stats.NumSequences)
	fmt.Printf("Total bases: %d\n", stats.TotalBases)
	fmt.Printf("Min length: %d\n", stats.MinLength)
	fmt.Printf("Max length: %d\n", stats.MaxLength)
	fmt.Printf("Average length: %.2f\n", stats.AvgLength)
	fmt.Printf("GC content: %.2f%%\n", stats.GCContent)
	if isFastq {
		fmt.Printf("Average quality: %.2f\n", stats.AvgQuality)
	}
}
