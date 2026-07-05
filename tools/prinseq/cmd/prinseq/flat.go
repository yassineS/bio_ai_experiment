package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/tools/prinseq/pkg/prinseq"
)

// isFlatInvocation reports whether the program was invoked with upstream
// prinseq-lite.pl's flat-flag CLI rather than one of our subcommands. The flat
// CLI starts with an option (e.g. "-fastq", "-out_good"); a subcommand is a
// bare word ("filter", "stats", ...). Global help/version switches are handled
// by main() before this is consulted, so any leading dash here means flat mode.
func isFlatInvocation(args []string) bool {
	if len(args) == 0 {
		return false
	}
	return strings.HasPrefix(args[0], "-")
}

// runFlat implements upstream prinseq-lite.pl's flat-flag CLI as a drop-in
// front end over our existing filter machinery. It accepts the common
// filter/trim flags in their upstream single-dash spelling (-fastq, -fasta,
// -out_good, -out_bad, -min_len, -trim_qual_right, ...) and writes output files
// using upstream's "<prefix>.fastq"/"<prefix>.fasta" naming convention
// (prinseq-lite.pl:1227-1301). The special prefixes "null" (suppress output)
// and "stdout" (write to STDOUT) are honoured for both -out_good and -out_bad.
//
// The numeric filter/trim semantics, record formatting and byte-exact output
// are entirely the existing prinseq.Filter / prinseq.FilterPaired code paths;
// runFlat only re-maps the CLI surface so that command lines written for
// upstream prinseq run unchanged against this port.
func runFlat(args []string) {
	fs := flag.NewFlagSet("prinseq", flag.ContinueOnError)
	fs.SetOutput(os.Stderr)

	// Input selection (upstream -fastq/-fasta and their paired -fastq2/-fasta2).
	var fastq, fasta, fastq2, fasta2 string
	fs.StringVar(&fastq, "fastq", "", "Input FASTQ file (use '-' for stdin)")
	fs.StringVar(&fasta, "fasta", "", "Input FASTA file (use '-' for stdin)")
	fs.StringVar(&fastq2, "fastq2", "", "Paired-end input FASTQ file 2")
	fs.StringVar(&fasta2, "fasta2", "", "Paired-end input FASTA file 2")

	// Output prefixes (upstream -out_good/-out_bad). Empty means "unset"; the
	// special values "null" and "stdout" are resolved after parsing.
	var outGood, outBad string
	fs.StringVar(&outGood, "out_good", "", "Output filename prefix for good reads (or 'null'/'stdout')")
	fs.StringVar(&outBad, "out_bad", "", "Output filename prefix for bad reads (or 'null'/'stdout')")

	// Length / GC / N filters.
	var minLen, maxLen, minGC, maxGC, maxNsN int
	fs.IntVar(&minLen, "min_len", 0, "Minimum sequence length")
	fs.IntVar(&maxLen, "max_len", 0, "Maximum sequence length")
	fs.IntVar(&minGC, "min_gc", 0, "Minimum GC content percentage")
	fs.IntVar(&maxGC, "max_gc", 0, "Maximum GC content percentage")
	fs.IntVar(&maxNsN, "ns_max_n", 0, "Maximum number of Ns allowed")
	var maxNsP int
	fs.IntVar(&maxNsP, "ns_max_p", 0, "Maximum percentage of Ns allowed")
	var rangeLen, rangeGC string
	fs.StringVar(&rangeLen, "range_len", "", "Keep lengths within ranges (e.g. 50-100,250-300)")
	fs.StringVar(&rangeGC, "range_gc", "", "Keep GC%% within ranges (e.g. 50-60,75-90)")

	// Mean-quality filters.
	var minQualMean, maxQualMean int
	fs.IntVar(&minQualMean, "min_qual_mean", 0, "Minimum mean quality score")
	fs.IntVar(&maxQualMean, "max_qual_mean", 0, "Maximum mean quality score")
	// Per-base quality-score filters are accepted for CLI compatibility but are
	// a known parity gap in this port (the package filters on mean quality, not
	// the per-base min/max). They are parsed so command lines do not error.
	var minQualScore, maxQualScore int
	fs.IntVar(&minQualScore, "min_qual_score", 0, "Minimum per-base quality score (accepted; see notes)")
	fs.IntVar(&maxQualScore, "max_qual_score", 0, "Maximum per-base quality score (accepted; see notes)")

	// Trimming.
	var trimLeft, trimRight, trimLeftP, trimRightP int
	var trimQualL, trimQualR int
	var trimNsLeft, trimNsRight, trimTailLeft, trimTailRight int
	var trimToLen int
	fs.IntVar(&trimLeft, "trim_left", 0, "Trim bases from 5' end")
	fs.IntVar(&trimRight, "trim_right", 0, "Trim bases from 3' end")
	fs.IntVar(&trimLeftP, "trim_left_p", 0, "Trim percentage from 5' end")
	fs.IntVar(&trimRightP, "trim_right_p", 0, "Trim percentage from 3' end")
	fs.IntVar(&trimQualL, "trim_qual_left", 0, "Quality threshold for 5' trimming")
	fs.IntVar(&trimQualR, "trim_qual_right", 0, "Quality threshold for 3' trimming")
	fs.IntVar(&trimNsLeft, "trim_ns_left", 0, "Trim poly-N from 5' end (min length)")
	fs.IntVar(&trimNsRight, "trim_ns_right", 0, "Trim poly-N from 3' end (min length)")
	fs.IntVar(&trimTailLeft, "trim_tail_left", 0, "Trim poly-A/T from 5' end (min length)")
	fs.IntVar(&trimTailRight, "trim_tail_right", 0, "Trim poly-A/T from 3' end (min length)")
	fs.IntVar(&trimToLen, "trim_to_len", 0, "Hard-trim reads to at most this length")

	// Sliding-window quality trimming.
	var trimQualWindow, trimQualStep int
	var trimQualType, trimQualRule string
	fs.IntVar(&trimQualWindow, "trim_qual_window", 0, "Sliding-window size for quality trimming (default 1)")
	fs.IntVar(&trimQualStep, "trim_qual_step", 0, "Step size for the quality-trim window (default 1)")
	fs.StringVar(&trimQualType, "trim_qual_type", "", "Window score: min (default), max, mean, sum")
	fs.StringVar(&trimQualRule, "trim_qual_rule", "", "Window rule: lt (default), gt, et")

	// Duplicate removal.
	var derep, derepMin int
	fs.IntVar(&derep, "derep", 0, "Remove duplicates: 1=exact, 4=revcomp, 5=both")
	fs.IntVar(&derepMin, "derep_min", 2, "Minimum occurrences to keep")
	var exactOnly bool
	fs.BoolVar(&exactOnly, "exact_only", false, "Restrict duplicate detection to exact dups")

	// Complexity filtering.
	var lcMethod string
	var lcThreshold int
	fs.StringVar(&lcMethod, "lc_method", "", "Low complexity method: dust or entropy")
	fs.IntVar(&lcThreshold, "lc_threshold", 0, "Low complexity threshold (default: 7 dust, 70 entropy)")

	// Sequence/header transforms and misc knobs.
	var noniupac, phred64, rmHeader, noQualHeader bool
	fs.BoolVar(&noniupac, "noniupac", false, "Filter sequences with bases outside ACGTN")
	fs.BoolVar(&phred64, "phred64", false, "Input FASTQ uses Phred+64 encoding")
	fs.BoolVar(&rmHeader, "rm_header", false, "Drop the original header comment")
	fs.BoolVar(&noQualHeader, "no_qual_header", false, "Emit a bare '+' line in FASTQ output")
	var seqCase, dnaRna, seqID, seqIDMappings, customParams, paramsFile string
	var lineWidth, seqNum, outFormat int
	fs.StringVar(&seqCase, "seq_case", "", "Force sequence case: upper|lower")
	fs.StringVar(&dnaRna, "dna_rna", "", "Convert T<->U: dna|rna")
	fs.StringVar(&seqID, "seq_id", "", "Rename passing records to '<PREFIX><N>'")
	fs.StringVar(&seqIDMappings, "seq_id_mappings", "", "Write '<orig>\\t<new>' TSV (requires -seq_id)")
	fs.StringVar(&customParams, "custom_params", "", "Dinucleotide-odds complexity rules")
	fs.StringVar(&paramsFile, "params", "", "Read parameters from a file")
	fs.IntVar(&lineWidth, "line_width", 0, "Wrap FASTA/QUAL output at N chars (0 = no wrap)")
	fs.IntVar(&seqNum, "seq_num", 0, "Keep only the first N passing records")
	fs.IntVar(&outFormat, "out_format", 0, "1=FASTA, 2=FASTA+QUAL, 3=FASTQ, 4=FASTQ+FASTA, 5=FASTQ+FASTA+QUAL")

	// Flat stats-reporting flags (upstream -stats_*). When ANY of these is
	// set, prinseq emits the summary-statistics TSV to STDOUT and skips the
	// filter/output path entirely (prinseq-lite.pl:676-752, 1944-2048).
	var statsInfo, statsLen, statsDupl, statsTag, statsDinuc, statsNs, statsAssembly, statsAll bool
	fs.BoolVar(&statsInfo, "stats_info", false, "Report read/base counts")
	fs.BoolVar(&statsLen, "stats_len", false, "Report length-distribution statistics")
	fs.BoolVar(&statsDupl, "stats_dupl", false, "Report duplicate-sequence statistics")
	fs.BoolVar(&statsTag, "stats_tag", false, "Report tag-sequence probability statistics")
	fs.BoolVar(&statsDinuc, "stats_dinuc", false, "Report dinucleotide-odds statistics")
	fs.BoolVar(&statsNs, "stats_ns", false, "Report ambiguous-base (N) statistics")
	fs.BoolVar(&statsAssembly, "stats_assembly", false, "Report assembly (Nx) statistics")
	fs.BoolVar(&statsAll, "stats_all", false, "Report every -stats_* group")

	// Switches we accept but treat as no-ops in the flat front end (they affect
	// only the stats/graph reporting paths, not the filtered read stream).
	var verbose bool
	var logFile string
	fs.BoolVar(&verbose, "verbose", false, "Verbose progress (accepted, no-op)")
	fs.StringVar(&logFile, "log", "", "Log file (accepted, no-op)")

	fs.Usage = func() {
		fmt.Fprint(os.Stderr, `Usage: prinseq -fastq FILE -out_good PREFIX [options]

Drop-in front end for prinseq-lite.pl's flat-flag CLI. Reads a FASTA/FASTQ
file, applies the requested filter/trim options, and writes the passing reads
to "<out_good>.fastq" (or .fasta) and the rejected reads to "<out_bad>.fastq".

Input:
  -fastq FILE / -fasta FILE        Input file (use '-' for stdin)
  -fastq2 FILE / -fasta2 FILE      Paired-end mate file

Output:
  -out_good PREFIX                 Good-read prefix (or 'null'/'stdout')
  -out_bad PREFIX                  Bad-read prefix (or 'null'/'stdout')
  -out_format INT                  1..5 (default: mirror input format)

Common filters:  -min_len -max_len -min_gc -max_gc -ns_max_n -ns_max_p
                 -min_qual_mean -max_qual_mean -range_len -range_gc -noniupac
                 -lc_method -lc_threshold -derep -derep_min
Common trimming: -trim_left -trim_right -trim_left_p -trim_right_p
                 -trim_qual_left -trim_qual_right -trim_qual_type -trim_qual_rule
                 -trim_qual_window -trim_qual_step -trim_to_len
                 -trim_ns_left -trim_ns_right -trim_tail_left -trim_tail_right
Misc:            -phred64 -seq_case -dna_rna -rm_header -no_qual_header
                 -line_width -seq_num -seq_id -seq_id_mappings -custom_params

(The subcommand form "prinseq filter ..." remains available.)
`)
	}

	if err := fs.Parse(args); err != nil {
		// flag already printed the error and usage under ContinueOnError; mirror
		// upstream's pod2usage(2) exit code for a bad option.
		os.Exit(2)
	}

	// Resolve input file and format. Upstream accepts exactly one of
	// -fastq / -fasta for the primary file (prinseq-lite.pl:777-789).
	var input1, input2 string
	var isFastq bool
	switch {
	case fastq != "":
		input1, isFastq = fastq, true
		input2 = fastq2
	case fasta != "":
		input1, isFastq = fasta, false
		input2 = fasta2
	default:
		fmt.Fprintln(os.Stderr, "Error: need an input file via -fastq or -fasta")
		os.Exit(1)
	}
	isPaired := input2 != ""

	// Flat stats-reporting path. When any -stats_* flag is set, emit the
	// summary-statistics TSV to STDOUT and skip filtering entirely, matching
	// upstream prinseq-lite.pl (lines 676-752, 1944-2048). -stats_all enables
	// every group.
	groups := prinseq.StatsGroups{
		Info:     statsInfo,
		Len:      statsLen,
		Dupl:     statsDupl,
		Tag:      statsTag,
		Dinuc:    statsDinuc,
		Ns:       statsNs,
		Assembly: statsAssembly,
	}
	if statsAll {
		groups = prinseq.StatsGroupsAll()
	}
	if groups.Any() {
		reader, err := openInput(input1)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening file: %v\n", err)
			os.Exit(1)
		}
		defer reader.Close()
		lines, err := prinseq.CollectFlatStats(reader, isFastq, groups)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error collecting statistics: %v\n", err)
			os.Exit(1)
		}
		w := bufio.NewWriter(os.Stdout)
		for _, l := range lines {
			fmt.Fprintln(w, l)
		}
		if err := w.Flush(); err != nil {
			fmt.Fprintf(os.Stderr, "Error writing statistics: %v\n", err)
			os.Exit(1)
		}
		return
	}

	// Default quality-trim window/step to 1 when quality trimming is active,
	// mirroring upstream's lazy default-fill and our subcommand path.
	if trimQualL > 0 || trimQualR > 0 {
		if trimQualWindow <= 0 {
			trimQualWindow = 1
		}
		if trimQualStep <= 0 {
			trimQualStep = 1
		}
	}

	// Validate trim_qual_rule / trim_qual_type, matching upstream's
	// "invalid value for ..." errors (prinseq-lite.pl:880-894).
	if trimQualRule != "" && trimQualRule != "lt" && trimQualRule != "gt" && trimQualRule != "et" {
		fmt.Fprintln(os.Stderr, "Error: invalid value for -trim_qual_rule (expected lt, gt or et)")
		os.Exit(1)
	}
	if trimQualType != "" && trimQualType != "min" && trimQualType != "max" && trimQualType != "mean" && trimQualType != "sum" {
		fmt.Fprintln(os.Stderr, "Error: invalid value for -trim_qual_type (expected min, max, mean or sum)")
		os.Exit(1)
	}

	// Quality encoding: -phred64 selects Illumina (Phred+64) encoding.
	qualType := "sanger"
	if phred64 {
		qualType = "illumina"
	}

	// -lc_method default thresholds (prinseq-lite.pl: 7 for dust, 70 entropy).
	lcThresholdF := float64(lcThreshold)
	if lcMethod != "" && lcThreshold == 0 {
		if lcMethod == "dust" {
			lcThresholdF = 7
		} else if lcMethod == "entropy" {
			lcThresholdF = 70
		}
	}

	// -seq_id_mappings requires -seq_id (prinseq-lite.pl:945-946).
	if seqIDMappings != "" && seqID == "" {
		fmt.Fprintln(os.Stderr, "Error: option -seq_id_mappings requires option -seq_id")
		os.Exit(1)
	}

	// Validate -out_format (prinseq-lite.pl:771).
	if outFormat < 0 || outFormat > 5 {
		fmt.Fprintln(os.Stderr, "Error: -out_format must be an integer between 1 and 5")
		os.Exit(1)
	}

	// Resolve the effective output format for filename selection. When unset
	// (0), upstream mirrors the input: FASTQ input -> format 3, FASTA -> 1
	// (prinseq-lite.pl:785-789). The file extension is .fastq for formats
	// 3/4/5 and .fasta otherwise (prinseq-lite.pl:1245-1246).
	effFormat := outFormat
	if effFormat == 0 {
		if isFastq {
			effFormat = 3
		} else {
			effFormat = 1
		}
	}
	goodExt := ".fasta"
	if effFormat == 3 || effFormat == 4 || effFormat == 5 {
		goodExt = ".fastq"
	}

	// -out_good / -out_bad equal-prefix check (prinseq-lite.pl:797-798).
	if outGood != "" && outBad != "" && outGood == outBad &&
		outGood != "null" && outGood != "stdout" {
		fmt.Fprintln(os.Stderr, "Error: the output names for -out_good and -out_bad have to be different")
		os.Exit(1)
	}

	// Resolve the effective line width for FASTA/QUAL wrapping, mirroring the
	// subcommand path: pure-FASTQ never wraps; an explicit -line_width wins;
	// otherwise the default is 60 (prinseq-lite.pl:931-939).
	lineWidthSet := flagSet(fs, "line_width")
	if effFormat == 3 {
		lineWidth = 0
	} else if !lineWidthSet {
		lineWidth = 60
	}

	// Validate -seq_case / -dna_rna domains (prinseq-lite.pl:912-924).
	if seqCase != "" && seqCase != "upper" && seqCase != "lower" {
		fmt.Fprintln(os.Stderr, "Error: -seq_case must be 'upper' or 'lower'")
		os.Exit(1)
	}
	if dnaRna != "" && dnaRna != "dna" && dnaRna != "rna" {
		fmt.Fprintln(os.Stderr, "Error: -dna_rna must be 'dna' or 'rna'")
		os.Exit(1)
	}

	opts := prinseq.FilterOptions{
		MinLen:        minLen,
		MaxLen:        maxLen,
		MinGC:         float64(minGC),
		MaxGC:         float64(maxGC),
		MinQualMean:   float64(minQualMean),
		MaxQualMean:   float64(maxQualMean),
		MaxNsP:        float64(maxNsP),
		MaxNsN:        maxNsN,
		TrimLeft:      trimLeft,
		TrimRight:     trimRight,
		TrimLeftP:     trimLeftP,
		TrimRightP:    trimRightP,
		TrimQualL:     trimQualL,
		TrimQualR:     trimQualR,
		TrimNsLeft:    trimNsLeft,
		TrimNsRight:   trimNsRight,
		TrimTailLeft:  trimTailLeft,
		TrimTailRight: trimTailRight,
		Derep:         derep,
		DerepMin:      derepMin,
		QualType:      qualType,
		LcMethod:      lcMethod,
		LcThreshold:   lcThresholdF,
		NonIUPAC:      noniupac,
		SeqID:         seqID,
		OutFormat:     outFormat,

		TrimQualWindow: trimQualWindow,
		TrimQualStep:   trimQualStep,
		TrimQualType:   trimQualType,
		TrimQualRule:   trimQualRule,
		TrimToLen:      trimToLen,
		RangeLen:       rangeLen,
		RangeGC:        rangeGC,

		SeqCase:       seqCase,
		DNARNA:        dnaRna,
		RmHeader:      rmHeader,
		NoQualHeader:  noQualHeader,
		SeqNum:        seqNum,
		ExactOnly:     exactOnly,
		QualLineWidth: lineWidth,
		CustomParams:  prinseq.ParseCustomParams(customParams),
	}

	// seq_id_mappings TSV writer (truncate on each run, like upstream).
	if seqIDMappings != "" {
		mapW, err := os.Create(seqIDMappings)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating seq_id_mappings file: %v\n", err)
			os.Exit(1)
		}
		defer mapW.Close()
		opts.SeqIDMap = mapW
	}

	// Resolve the bad-read writer. Upstream defaults out_bad to a temp file
	// when the flag is omitted; for the flat front end we treat an omitted
	// -out_bad like "null" (suppress) so a bare command does not litter the
	// working directory, while honouring an explicit prefix / null / stdout.
	badWriter, badCloser, badIsStdout, err := openFlatOutput(outBad, goodExt, true)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening bad output file: %v\n", err)
		os.Exit(1)
	}
	if badCloser != nil {
		defer badCloser.Close()
	}
	opts.OutBad = badWriter

	// Resolve the good-read writer (the primary output stream).
	goodWriter, goodCloser, goodIsStdout, err := openFlatOutput(outGood, goodExt, false)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening good output file: %v\n", err)
		os.Exit(1)
	}
	if goodCloser != nil {
		defer goodCloser.Close()
	}
	if badIsStdout && goodIsStdout {
		fmt.Fprintln(os.Stderr, "Error: -out_good and -out_bad cannot both write to stdout")
		os.Exit(1)
	}

	// Multi-stream out_format (2/4/5) needs the .fasta/.qual siblings, derived
	// from the out_good prefix (prinseq-lite.pl:1302-1348). These cannot be
	// streamed to stdout.
	if effFormat == 2 || effFormat == 4 || effFormat == 5 {
		if outGood == "" || outGood == "null" || outGood == "stdout" {
			fmt.Fprintln(os.Stderr, "Error: -out_format 2/4/5 require -out_good to provide a filename prefix")
			os.Exit(1)
		}
		if effFormat == 4 || effFormat == 5 {
			f, err := os.Create(outGood + ".fasta")
			if err != nil {
				fmt.Fprintf(os.Stderr, "Error creating .fasta output: %v\n", err)
				os.Exit(1)
			}
			defer f.Close()
			opts.FastaOut = f
		}
		if effFormat == 2 || effFormat == 5 {
			f, err := os.Create(outGood + ".qual")
			if err != nil {
				fmt.Fprintf(os.Stderr, "Error creating .qual output: %v\n", err)
				os.Exit(1)
			}
			defer f.Close()
			opts.QualOut = f
		}
	}

	if isPaired {
		reader1, err := openInput(input1)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening file 1: %v\n", err)
			os.Exit(1)
		}
		defer reader1.Close()
		reader2, err := openInput(input2)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening file 2: %v\n", err)
			os.Exit(1)
		}
		defer reader2.Close()

		// Upstream writes paired good reads to "<prefix>_1.fastq" and
		// "<prefix>_2.fastq" (prinseq-lite.pl:1235-1243). The second mate
		// stream is derived from the out_good prefix here.
		if outGood == "" || outGood == "null" || outGood == "stdout" {
			fmt.Fprintln(os.Stderr, "Error: paired-end input requires -out_good to provide a filename prefix")
			os.Exit(1)
		}
		w2, err := os.Create(outGood + "_2" + goodExt)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file 2: %v\n", err)
			os.Exit(1)
		}
		defer w2.Close()

		if err := prinseq.FilterPaired(reader1, reader2, goodWriter, w2, isFastq, opts); err != nil {
			fmt.Fprintf(os.Stderr, "Error filtering paired sequences: %v\n", err)
			os.Exit(1)
		}
		return
	}

	reader, err := openInput(input1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening file: %v\n", err)
		os.Exit(1)
	}
	defer reader.Close()

	if err := prinseq.Filter(reader, goodWriter, isFastq, opts); err != nil {
		fmt.Fprintf(os.Stderr, "Error filtering sequences: %v\n", err)
		os.Exit(1)
	}
}

// openFlatOutput resolves an upstream -out_good / -out_bad prefix into a
// destination writer. It honours the special values "null" (suppress the
// stream — returns a nil writer) and "stdout" (write to STDOUT). A plain prefix
// is turned into "<prefix><ext>" (e.g. "<prefix>.fastq"), matching upstream's
// single-end naming (prinseq-lite.pl:1245-1246, 1286-1287).
//
// When isBad is true an empty prefix is treated as "null" (the flat front end
// suppresses bad output by default so a bare invocation does not litter the
// working directory). When isBad is false an empty good prefix defaults to
// "stdout" so "prinseq -fastq in.fq -min_len 40" prints to the terminal.
//
// The returned closer is non-nil only when a real file was opened (callers must
// Close it); writer is nil for the suppressed ("null") case.
func openFlatOutput(prefix, ext string, isBad bool) (writer io.Writer, closer io.Closer, isStdout bool, err error) {
	switch prefix {
	case "null":
		return nil, nil, false, nil
	case "stdout":
		return os.Stdout, nil, true, nil
	case "":
		if isBad {
			return nil, nil, false, nil
		}
		return os.Stdout, nil, true, nil
	default:
		f, e := os.Create(prefix + ext)
		if e != nil {
			return nil, nil, false, e
		}
		return f, f, false, nil
	}
}
