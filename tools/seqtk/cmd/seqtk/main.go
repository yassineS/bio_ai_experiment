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
//	mutfa      Apply point mutations to a FASTA file
//	randbase   Randomly resolve IUPAC ambiguity codes
//	hpc        Homopolymer-compress sequences
//	gap        Find gap (non-ACGT) regions in FASTA
//	gc         Find GC-rich (or AT-rich) regions in FASTA
//	dropse     Drop unpaired reads from an interleaved FASTA/Q
//	rename     Rename records (renumber, optional prefix)
//	split      Split a FASTA/Q into N round-robin output files
//	size       Print total record count and total sequence length
//	famask     Apply a FASTA mask to a source FASTA (X = keep, x =
//	           soft-mask, other = overwrite)
//	mergefa    Merge two FASTA/FASTQ files base-by-base via IUPAC
package main

import (
	"flag"
	"fmt"
	"os"
	"time"

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
	case "mutfa":
		mutfaCommand()
	case "randbase":
		randbaseCommand()
	case "hpc":
		hpcCommand()
	case "gap":
		gapCommand()
	case "gc":
		gcCommand()
	case "dropse":
		dropseCommand()
	case "rename":
		renameCommand()
	case "split":
		splitCommand()
	case "size":
		sizeCommand()
	case "famask":
		famaskCommand()
	case "mergefa":
		mergefaCommand()
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
  mutfa      Apply point mutations to a FASTA file from a TSV mutation list
  randbase   Replace IUPAC ambiguity bases with a random expansion
  hpc        Homopolymer-compress sequences (collapse runs of identical bases)
  gap        Find gap (non-ACGT) regions in FASTA, emit BED3
  gc         Find GC-rich (or AT-rich) regions in FASTA, emit BED4
  dropse     Drop unpaired reads from an interleaved FASTA/Q stream
  rename     Rename records as <prefix><N>; pairs share N
  split      Split input into N round-robin output files <prefix>.NNNNN.fa
  size       Print '<num_records>\t<total_bases>' (upstream summary form)
  famask     Apply a FASTA mask file to a source FASTA
  mergefa    Merge two FASTA/FASTQ inputs base-by-base via IUPAC codes
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
  seqtk mutfa ref.fa muts.tsv > mutated.fa
  seqtk randbase -s 42 ambig.fa > resolved.fa
  seqtk hpc reads.fa > collapsed.fa
  seqtk gap -l 10 genome.fa > gaps.bed
  seqtk gc -f 0.7 -l 50 genome.fa > gc_rich.bed
  seqtk dropse interleaved.fq > paired.fq
  seqtk rename reads.fq SAMPLE_ > renamed.fq
  seqtk split  -n 4 part reads.fq      # writes part.00001.fa .. part.00004.fa
  seqtk size   genome.fa               # "<num_records>\t<total_bases>"
  seqtk famask src.fa mask.fa > out.fa # X=keep, x=soft-mask, else=overwrite
  seqtk mergefa a.fa b.fa > merged.fa  # IUPAC merge; -i / -m / -h / -r / -q
                                       # select the merge mode

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
	var phred64, summary bool

	cliflag.BoolVar(fs, &phred64, "6", "phred64", false, "Use Phred+64 quality encoding for FASTQ (default: Phred+33)")
	cliflag.BoolVar(fs, &summary, "", "summary", false, "Emit summary statistics (legacy) instead of upstream per-record output")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk comp [options] <input>

Get per-record nucleotide composition (upstream "seqtk comp" output):
one tab-separated row per input record with the columns:

    name\tlen\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts

Arguments:
  <input>    Input file (use '-' for stdin, supports .gz and .bz2)

Options:
  -6, --phred64          Use Phred+64 quality encoding for FASTQ (default: Phred+33)
      --summary          Emit summary statistics instead of upstream-format
                         per-record rows (legacy behaviour from before the
                         2026-05-14 parity audit).

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

	if !summary {
		// Upstream-compatible per-record output (the default).
		if err := seqtk.Comp(input, os.Stdout); err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v\n", err)
			os.Exit(1)
		}
		return
	}

	// --summary: legacy aggregate stats; still useful for quick eyeballing.
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

func mutfaCommand() {
	fs := flag.NewFlagSet("mutfa", flag.ExitOnError)
	var output string

	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout, supports .gz)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk mutfa [options] <in.fa> <mutfile>

Apply point mutations to a FASTA file. <mutfile> is a TSV with at least
three whitespace/tab-separated columns per line:

  chrom    pos(1-based)    base

For compatibility with upstream seqtk's four-column format
  chrom    pos    ref    alt
the new base is taken from column 4 when there are four or more columns.
Lines starting with '#' and blank lines are ignored.

Output is FASTA written to stdout, preserving the line-width layout of the
input. Mutations are applied on the forward strand. Names from the mutation
file that are not present in the input, and positions past the end of their
sequence, produce a warning on stderr and are skipped.

Arguments:
  <in.fa>    Input FASTA file (use '-' for stdin, supports .gz)
  <mutfile>  TSV mutation list (use '-' for stdin, supports .gz)

Options:
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk mutfa ref.fa muts.tsv > mutated.fa
  seqtk mutfa ref.fa.gz muts.tsv -o mutated.fa.gz

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 2 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)
	mutFile := fs.Arg(1)

	if inputFile == "-" && mutFile == "-" {
		fmt.Fprintf(os.Stderr, "Error: both inputs cannot be stdin ('-')\n")
		os.Exit(1)
	}

	input, err := seqtk.OpenInput(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	mut, err := seqtk.OpenInput(mutFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening mutation file: %v\n", err)
		os.Exit(1)
	}
	defer mut.Close()

	out, err := seqtk.OpenOutput(output)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if err := seqtk.Mutfa(input, mut, out); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func randbaseCommand() {
	fs := flag.NewFlagSet("randbase", flag.ExitOnError)
	var output string
	var seed int64

	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout, supports .gz)")
	// Negative sentinel means "no -s provided"; we'll time-seed in that case.
	// Use Int64Var directly since cliflag doesn't expose an int64 helper.
	fs.Int64Var(&seed, "s", -1, "Random seed for reproducibility (default: time-seeded)")
	fs.Int64Var(&seed, "seed", -1, "Random seed for reproducibility (default: time-seeded)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk randbase [options] <in.fa>

Replace every IUPAC ambiguity base in the input FASTA with one of the
unambiguous bases it represents, chosen uniformly at random. Case is
preserved (e.g. 'r' becomes 'a' or 'g'). The IUPAC table used is:

  R -> A,G    Y -> C,T    S -> G,C    W -> A,T
  K -> G,T    M -> A,C    B -> C,G,T  D -> A,G,T
  H -> A,C,T  V -> A,C,G  N -> A,C,G,T

Output is FASTA written to stdout, preserving the line-width layout of
the input.

Arguments:
  <in.fa>    Input FASTA file (use '-' for stdin, supports .gz)

Options:
  -s, --seed INT         Random seed for reproducibility (default: time-seeded)
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk randbase -s 42 ambig.fa > resolved.fa
  seqtk randbase ambig.fa.gz -o resolved.fa.gz

`)
	}

	fs.Parse(os.Args[2:])

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

	effectiveSeed := seed
	if effectiveSeed < 0 {
		effectiveSeed = time.Now().UnixNano()
	}

	if err := seqtk.Randbase(input, out, effectiveSeed); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func hpcCommand() {
	fs := flag.NewFlagSet("hpc", flag.ExitOnError)
	var output string

	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout, supports .gz)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk hpc [options] <in.fa>

Homopolymer-compress sequences: collapse every maximal run of identical
bases to a single base. Sequence names are preserved; the compressed
sequence is written on a single line with no wrapping, matching upstream
seqtk hpc.

Input may be FASTA or FASTQ (auto-detected via the first non-whitespace
byte: '>' => FASTA, '@' => FASTQ). Output is always FASTA.

Arguments:
  <in.fa>    Input FASTA/FASTQ file (use '-' for stdin, supports .gz)

Options:
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk hpc reads.fa > collapsed.fa
  seqtk hpc reads.fq.gz -o collapsed.fa.gz

`)
	}

	fs.Parse(os.Args[2:])

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

	if err := seqtk.HPC(input, out); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func gapCommand() {
	fs := flag.NewFlagSet("gap", flag.ExitOnError)
	var minSize int
	var output string

	cliflag.IntVar(fs, &minSize, "l", "min-size", seqtk.DefaultGapMinSize, "Minimum gap-run length to report")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout, supports .gz)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk gap [options] <in.fa>

Find gap regions in a FASTA file. A "gap" is a maximal run of bytes that are
not A, C, G or T (case-insensitive) — i.e. N's, IUPAC ambiguity codes and any
other non-ACGT byte all count. For every gap of length >= -l a BED3 record is
written to stdout: chrom\tstart\tend (0-based half-open).

Arguments:
  <in.fa>    Input FASTA file (use '-' for stdin, supports .gz)

Options:
  -l, --min-size INT     Minimum gap-run length to report (default: %d)
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk gap genome.fa > gaps.bed
  seqtk gap -l 10 genome.fa > short_gaps.bed
  zcat genome.fa.gz | seqtk gap - > gaps.bed

`, seqtk.DefaultGapMinSize)
	}

	fs.Parse(os.Args[2:])

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

	if err := seqtk.Gap(input, out, seqtk.GapOptions{MinSize: minSize}); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func gcCommand() {
	fs := flag.NewFlagSet("gc", flag.ExitOnError)
	var minLen int
	var minFrac, xDropoff float64
	var isAT bool
	var output string

	cliflag.IntVar(fs, &minLen, "l", "min-length", seqtk.DefaultGCMinLength, "Minimum region length to report")
	cliflag.Float64Var(fs, &minFrac, "f", "min-frac", seqtk.DefaultGCMinFrac, "Min GC fraction (or AT fraction with -w)")
	cliflag.Float64Var(fs, &xDropoff, "x", "x-dropoff", seqtk.DefaultGCXDropoff, "X-dropoff threshold")
	cliflag.BoolVar(fs, &isAT, "w", "at", false, "Identify high-AT regions instead of high-GC")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout, supports .gz)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk gc [options] <in.fa>

Find GC-rich (or, with -w, AT-rich) regions in a FASTA file using upstream
seqtk's X-dropoff scoring algorithm. Every "hit" base contributes
+(1-f)/f to a running score, every non-hit contributes -1, and a region is
emitted whenever the score drops below zero or X below its running maximum,
provided the region is at least -l bases long.

Output is BED4 (0-based half-open): chrom\tstart\tend\thits, where hits is
the number of GC (or AT) positions in [start, end).

Arguments:
  <in.fa>    Input FASTA file (use '-' for stdin, supports .gz)

Options:
  -w, --at               Identify high-AT regions instead of high-GC
  -f, --min-frac FLOAT   Min GC fraction (or AT fraction with -w) [%.2f]
  -l, --min-length INT   Min region length to output [%d]
  -x, --x-dropoff FLOAT  X-dropoff threshold [%.1f]
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk gc genome.fa > gc_rich.bed
  seqtk gc -f 0.75 -l 100 genome.fa > strong_gc.bed
  seqtk gc -w -f 0.7 genome.fa > at_rich.bed

`, seqtk.DefaultGCMinFrac, seqtk.DefaultGCMinLength, seqtk.DefaultGCXDropoff)
	}

	fs.Parse(os.Args[2:])

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

	opts := seqtk.GCOptions{
		MinLength: minLen,
		MinFrac:   minFrac,
		XDropoff:  xDropoff,
		IsAT:      isAT,
	}
	if err := seqtk.GC(input, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func dropseCommand() {
	fs := flag.NewFlagSet("dropse", flag.ExitOnError)
	var output string

	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout, supports .gz)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk dropse [options] <in.fq>

Drop unpaired (singleton) reads from an interleaved FASTA/FASTQ stream.
Two adjacent records are considered mates if their names are identical
after stripping a trailing "/<digit>" (e.g. "/1" or "/2"); any record
whose immediate neighbour does not match this rule is silently dropped,
matching upstream "seqtk dropse" byte-for-byte.

Upstream surface (verified against reference_code/seqtk v1.5-r133):
the subcommand takes no flags, only a single positional input. The
"-o/--output FILE" option here is a Go-port convenience covered by
the project's CLI conventions (it does not alter parity).

Arguments:
  <in.fq>    Interleaved input (use '-' for stdin, supports .gz)

Options:
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk dropse interleaved.fq > paired.fq
  cat reads.fq.gz | seqtk dropse - > paired.fq

`)
	}

	fs.Parse(os.Args[2:])

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

	if err := seqtk.Dropse(input, out); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func renameCommand() {
	fs := flag.NewFlagSet("rename", flag.ExitOnError)
	var output string

	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout, supports .gz)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk rename <in.fq> [prefix]

Rename FASTA/FASTQ records to "<prefix><N>" where N is a 1-based counter
that advances per "fragment". Two adjacent records sharing a name modulo
a trailing "/<digit>" suffix are treated as a pair and share the same N
(mirroring upstream "seqtk rename" byte-for-byte). The prefix is
optional; when omitted, names are bare integers ("1", "2", ...).

Comments in the original header (anything past the first whitespace
following the name) are preserved verbatim after the new name.

Output format mirrors the input (FASTA -> FASTA, FASTQ -> FASTQ);
sequence and quality are emitted on a single un-wrapped line, matching
upstream stk_printseq_renamed(..., line_len=0).

Upstream surface (verified against reference_code/seqtk v1.5-r133): the
subcommand takes no flags and accepts only the positional arguments
<in.fa> and an optional <prefix>. The "-o/--output FILE" option here is
the project-wide Go-port convenience; it does not affect parity.

Arguments:
  <in.fq>    Input FASTA/FASTQ (use '-' for stdin, supports .gz)
  [prefix]   Optional name prefix (default: empty, so names are
             bare integers)

Options:
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk rename reads.fq SAMPLE_ > renamed.fq
  seqtk rename contigs.fa > numbered.fa

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 1 {
		fs.Usage()
		os.Exit(1)
	}

	inputFile := fs.Arg(0)
	prefix := ""
	if fs.NArg() > 1 {
		prefix = fs.Arg(1)
	}

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

	if err := seqtk.Rename(input, out, prefix); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func splitCommand() {
	fs := flag.NewFlagSet("split", flag.ExitOnError)
	var n, lineLen int

	cliflag.IntVar(fs, &n, "n", "num", seqtk.DefaultSplitN, "Number of output files")
	cliflag.IntVar(fs, &lineLen, "l", "line-length", 0, "Sequence/quality line length (0 = no wrap)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk split [options] <prefix> <in.fa>

Round-robin every record from <in.fa> across N output files named
"<prefix>.<5-digit 1-based>.fa" (note the literal ".fa" suffix --
upstream uses it even for FASTQ input). The first record goes to
<prefix>.00001.fa, the second to <prefix>.00002.fa, ..., the (N+1)-th
wraps to <prefix>.00001.fa, and so on.

Within each output file the input format is preserved (FASTA stays
FASTA, FASTQ stays FASTQ). With "-l INT" sequence (and FASTQ quality)
lines are wrapped at INT characters; the upstream default of 0 keeps
each sequence/quality on a single line.

Output files are written uncompressed even when their name ends in
".fa", matching upstream byte-for-byte.

Upstream surface (verified against reference_code/seqtk v1.5-r133):
flags -n INT (default %d) and -l INT (default 0). Positional arguments
are <prefix> followed by <in.fa>.

Arguments:
  <prefix>   Output file-name prefix (e.g. "part")
  <in.fa>    Input FASTA/FASTQ file (use '-' for stdin, supports .gz)

Options:
  -n, --num INT          Number of output files [%d]
  -l, --line-length INT  Sequence/quality line length (0 = no wrap)

Examples:
  seqtk split -n 4 part reads.fq      # writes part.00001.fa .. part.00004.fa
  seqtk split -n 8 -l 60 chunk genome.fa

`, seqtk.DefaultSplitN, seqtk.DefaultSplitN)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 2 {
		fs.Usage()
		os.Exit(1)
	}

	prefix := fs.Arg(0)
	inputFile := fs.Arg(1)

	input, err := seqtk.OpenInput(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	opts := seqtk.SplitOptions{N: n, LineLen: lineLen, Prefix: prefix}
	if err := seqtk.Split(input, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func sizeCommand() {
	fs := flag.NewFlagSet("size", flag.ExitOnError)
	var output string

	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout, supports .gz)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk size <in.fq>

Print one tab-separated line on stdout with the total number of records
and the total number of bases across the input:

    <num_records>\t<total_bases>\n

Matches upstream "seqtk size" byte-for-byte (verified against
reference_code/seqtk v1.5-r133).

Upstream surface: the subcommand takes no flags. The "-o/--output FILE"
option here is the project-wide Go-port convenience.

Arguments:
  <in.fq>    Input FASTA/FASTQ file (use '-' for stdin, supports .gz)

Options:
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk size genome.fa
  seqtk size reads.fq.gz

`)
	}

	fs.Parse(os.Args[2:])

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

	if err := seqtk.Size(input, out); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func famaskCommand() {
	fs := flag.NewFlagSet("famask", flag.ExitOnError)
	var output string

	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout, supports .gz)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk famask <src.fa> <mask.fa>

Apply a FASTA-format mask to a source FASTA, byte-for-byte:

  mask byte == 'X'   keep the source base unchanged
  mask byte == 'x'   lowercase the source base (soft-mask)
  any other byte     overwrite the source base with the mask byte

Records are paired by stream order; name mismatches and length
mismatches print a warning to stderr (matching upstream) and the
shorter length is used. Output is FASTA wrapped at 60 bases per line,
matching upstream "seqtk famask" byte-for-byte.

Upstream surface (verified against reference_code/seqtk/seqtk.c v1.5-r133
line 872, getopt("")) — the subcommand takes NO flags whatsoever, only
the two positional inputs. The "-o/--output FILE" option here is the
project-wide Go-port convenience and does not affect parity.

Arguments:
  <src.fa>   Source FASTA (use '-' for stdin, supports .gz)
  <mask.fa>  Mask FASTA  (use '-' for stdin, supports .gz)

Options:
  -o, --output FILE      Output file (default: stdout, supports .gz)

Examples:
  seqtk famask genome.fa repeats.fa > masked.fa
  seqtk famask src.fa.gz mask.fa.gz -o out.fa.gz

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 2 {
		fs.Usage()
		os.Exit(1)
	}

	srcFile := fs.Arg(0)
	maskFile := fs.Arg(1)

	if srcFile == "-" && maskFile == "-" {
		fmt.Fprintf(os.Stderr, "Error: both inputs cannot be stdin ('-')\n")
		os.Exit(1)
	}

	src, err := seqtk.OpenInput(srcFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening source FASTA: %v\n", err)
		os.Exit(1)
	}
	defer src.Close()

	mask, err := seqtk.OpenInput(maskFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening mask FASTA: %v\n", err)
		os.Exit(1)
	}
	defer mask.Close()

	out, err := seqtk.OpenOutput(output)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer out.Close()

	if err := seqtk.Famask(src, mask, out); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func mergefaCommand() {
	fs := flag.NewFlagSet("mergefa", flag.ExitOnError)
	var quality int
	var intersect, haploid, mask, randhet bool
	var output string

	cliflag.IntVar(fs, &quality, "q", "quality", 0, "Quality threshold (PHRED+33; lowercase below this)")
	cliflag.BoolVar(fs, &intersect, "i", "intersect", false, "Take intersection (c0 & c1)")
	cliflag.BoolVar(fs, &haploid, "h", "haploid", false, "Suppress hets in the input")
	cliflag.BoolVar(fs, &mask, "m", "mask", false, "Lowercase when one input is N")
	cliflag.BoolVar(fs, &randhet, "r", "rand-het", false, "Pick a random allele from het")
	cliflag.StringVar(fs, &output, "o", "output", "", "Output file (default: stdout, supports .gz)")

	fs.Usage = func() {
		fmt.Fprintf(os.Stderr, `Usage: seqtk mergefa [options] <in1.fa> <in2.fa>

Merge two FASTA (or FASTQ) inputs base-by-base. For every paired
position the two bases are looked up in upstream's seq_nt16 table and
combined into a single IUPAC code:

  default     OR the two codes (c0 | c1). Het pairs collapse to the
              matching IUPAC ambiguity code (e.g. A+G -> R, C+T -> Y).
  -i          intersect: c0 & c1 (no overlap -> 'x').
  -m          mask: like -i but lowercase whenever either side is N.
  -r          pick a random allele on hets (uses Go math/rand).
  -h          haploid: heterozygous merges are lowercased.

Output case encodes confidence: a base is uppercase only when both
inputs are uppercase (or in the OR-modes when either is). With FASTQ
input, bases whose PHRED+33 quality is below -q are lowercased before
merging.

A "[mergefa] (same,diff,hom-het,het-hom,het-het)=(...)" summary is
written to stderr after the last record, matching upstream
byte-for-byte.

Upstream surface (verified against reference_code/seqtk/seqtk.c v1.5-r133
line 767, getopt("himrq:")):
  -q INT   quality threshold [0]
  -i       take intersection
  -m       lowercase when one of the inputs is N
  -r       pick a random allele from het
  -h       suppress hets in the input

The "-o/--output FILE" option here is the project-wide Go-port
convenience and does not affect parity.

Arguments:
  <in1.fa>   First FASTA/FASTQ  (use '-' for stdin, supports .gz)
  <in2.fa>   Second FASTA/FASTQ (use '-' for stdin, supports .gz)

Examples:
  seqtk mergefa a.fa b.fa > merged.fa
  seqtk mergefa -i a.fa b.fa > intersect.fa
  seqtk mergefa -q 20 a.fq b.fq > merged.fa

`)
	}

	fs.Parse(os.Args[2:])

	if fs.NArg() < 2 {
		fs.Usage()
		os.Exit(1)
	}

	if intersect && mask {
		fmt.Fprintf(os.Stderr, "Error: -i and -m cannot be applied at the same time\n")
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

	opts := seqtk.MergefaOptions{
		Quality:   quality,
		Intersect: intersect,
		Haploid:   haploid,
		Mask:      mask,
		RandHet:   randhet,
	}
	if err := seqtk.Mergefa(in1, in2, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}
