// bedintersect finds intersecting intervals between feature files, a drop-in
// re-implementation of `bedtools intersect`.
package main

import (
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedintersect/pkg/bedintersect"
)

const version = "1.0.0"

const usage = `bedintersect - Report overlaps between two feature files (BED/GFF/VCF/BAM)

Usage:
  bedintersect [options] -a <fileA> -b <fileB> [fileB2 ...]

Description:
  Report intervals in A that overlap intervals in B. -b may be followed by
  multiple database files. By default, reports the overlapping portion of each
  pair. The options below control what is reported and the overlap required.

Output options:
  -wa            Write the original A entry for each overlap.
  -wb            Write the original B entry for each overlap.
  -loj           Left outer join: report every A; a NULL B when no overlap.
  -wo            Write A and B plus the number of overlapping bases (overlaps only).
  -wao           Like -wo, but also report A with a NULL B and 0 when no overlap.
  -u             Write each original A entry once if any overlap is found.
  -c             For each A, report the total number of overlaps with B.
  -C             For each A, report the number of overlaps with each B file separately.
  -v             Only report A entries that have NO overlap with B.

Overlap options:
  -f FLOAT       Minimum overlap as a fraction of A (default 1bp).
  -F FLOAT       Minimum overlap as a fraction of B (default 1bp).
  -r             Require the fraction overlap be reciprocal for A AND B.
  -e             Require the minimum fraction be satisfied for A OR B.
  -s             Require the same strand.
  -S             Require opposite strands.
  -split         Treat split BAM/BED12 entries as distinct intervals.

Input / sorting options:
  -a FILE        Input file A (BED/GFF/VCF/BAM; '-' for stdin). Required.
  -b FILE...     Input file(s) B. Required.
  -abam FILE     Alias for -a with a BAM file.
  -ibam FILE     Alias for -a with a BAM file.
  -bed           With BAM/CRAM input, write output as BED instead of the default
                 binary alignments (a BAM/CRAM query writes BAM, or CRAM when a
                 reference is given, by default).
  -ubam          Write uncompressed (level-0) BAM output. The format choice
                 still follows the CRAM reference, matching upstream; -ubam
                 affects only the BAM compression mode.
  --cram-ref FA  CRAM reference FASTA. A CRAM query writes CRAM output (rather
                 than BAM) only when this (or CRAM_REFERENCE) is set.
  -names ...     Aliases for each B file (printed instead of a numeric file id).
  -filenames     Print each B file's name instead of a numeric file id.
  -sortout       Sort the per-A DB hits by position across all B files.
  -header        Print the header from the A file prior to results.
  -sorted        Validate that the inputs are coordinate-sorted (error if not).
  -g FILE        Genome file fixing the chromosome order for -sorted validation.
  -nonamecheck   Suppress the chromosome naming-convention warning.
  -o FILE        Output file (default: stdout).
  -m INT         Minimum overlap in bp (bedintersect extension, default 1).
  -d             Report distance to the nearest B feature (bedintersect extension).
  -k             Report the closest B feature for each A (bedintersect extension).
  -t             Use an interval tree for large B files (bedintersect extension).
  --stats        Print summary statistics to stderr (bedintersect extension).
  -h, --help     Show this help message.
  --version      Show version information and exit.

Notes:
  - Coordinates are 0-based, half-open [start, end).
  - Multiple B hits generate multiple output lines (unless -c/-C/-u/-v).
  - With multiple -b files, the DB-id column (-names/-filenames/numeric) is added.
`

// options collects the parsed command-line state.
type options struct {
	inputA    string
	inputB    []string
	output    string
	names     []string
	filenames bool

	genome      string
	sorted      bool
	noNameCheck bool

	minOverlap int
	fractionA  float64
	fractionB  float64
	reciprocal bool
	either     bool
	strand     bool
	opposite   bool
	split      bool

	writeA          bool
	writeB          bool
	leftJoin        bool
	writeOverlap    bool
	writeAllOverlap bool
	unique          bool
	count           bool
	countEach       bool
	invert          bool
	sortOut         bool
	header          bool

	distance bool
	closest  bool
	useTree  bool
	stats    bool

	// bedOutput records the -bed flag. With BAM/CRAM query input it forces BED
	// text output (the upstream default for BAM-A would be BAM); with text input
	// it has no effect. It is the switch that distinguishes "emit BAM" from
	// "emit BED" for a BAM/CRAM query file.
	bedOutput bool

	// uncompressedBAM records the -ubam flag. It selects an uncompressed
	// (level-0) BAM writer for the binary-output path: the BAM is still
	// BGZF-framed, just with stored DEFLATE blocks. As upstream, -ubam does NOT
	// turn a CRAM query into BAM — the BAM/CRAM format choice is driven solely
	// by the CRAM reference, so a CRAM query with a reference still writes CRAM
	// under -ubam, and -ubam then only affects the (BAM) compression mode.
	uncompressedBAM bool

	// cramRef names the CRAM reference FASTA. Upstream takes it from the global
	// `--cram-ref <fa>` flag or the CRAM_REFERENCE environment variable, and
	// emits CRAM (rather than BAM) output for a CRAM query only when it is set.
	cramRef string

	help        bool
	showVersion bool
}

func main() {
	opts, err := parseArgs(os.Args[1:])
	if err != nil {
		var unrec *unrecognizedParamError
		var rng *rangeError
		if errors.As(err, &unrec) || errors.As(err, &rng) {
			// Upstream prints just this banner (as the final stderr line) and exits.
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		fmt.Fprintln(os.Stderr, "Use -h for help")
		os.Exit(1)
	}

	if opts.help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if opts.showVersion {
		fmt.Printf("bedintersect version %s\n", version)
		os.Exit(0)
	}

	if opts.inputA == "" || len(opts.inputB) == 0 {
		fmt.Fprintln(os.Stderr, "Error: Both -a and -b are required")
		fmt.Fprint(os.Stderr, "\nUsage: bedintersect -a fileA -b fileB [fileB2 ...] [options]\n")
		fmt.Fprintln(os.Stderr, "Use -h for help")
		os.Exit(1)
	}

	// Upstream accepts the literal "stdin" as a synonym for "-".
	opts.inputA = normalizeStdin(opts.inputA)
	for i := range opts.inputB {
		opts.inputB[i] = normalizeStdin(opts.inputB[i])
	}

	readerACloser, err := iohelper.OpenReader(opts.inputA)
	if err != nil {
		// Match upstream's exact "Unable to open file" message for file errors.
		fmt.Fprintf(os.Stderr, "Error: Unable to open file %s. Exiting.\n", opts.inputA)
		os.Exit(1)
	}
	defer readerACloser.Close()

	// Classify the query file. Upstream determines the output type from the
	// query (-a) file: a BAM/CRAM query writes binary alignments by default
	// (unless -bed), and some flags are gated accordingly. The sniffed reader
	// replaces readerA so the bytes consumed while probing are not lost.
	queryFormat, readerA, err := bedintersect.ClassifyQueryInput(readerACloser)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: Unable to open file %s. Exiting.\n", opts.inputA)
		os.Exit(1)
	}
	// CRAM_REFERENCE supplements --cram-ref (upstream reads both; the flag wins).
	if opts.cramRef == "" {
		opts.cramRef = os.Getenv("CRAM_REFERENCE")
	}
	// binaryOutput is true when the surviving alignments must be written back out
	// as binary (BAM or CRAM): the query is BAM/CRAM and the user did not force
	// BED text with -bed. outputFormat selects CRAM vs BAM, mirroring upstream
	// (RecordOutputMgr + BamWriter::Open): a CRAM query writes CRAM only when a
	// CRAM reference is available, otherwise BAM. Upstream gates the format purely
	// on the reference — -ubam controls only the BAM compression mode (a no-op in
	// upstream's writer) and does NOT force BAM when a reference is set, so we do
	// not let it override the format either.
	binaryOutput := queryFormat != bedintersect.QueryText && !opts.bedOutput
	outputFormat := bedintersect.OutputBAM
	if queryFormat == bedintersect.QueryCRAM && opts.cramRef != "" {
		outputFormat = bedintersect.OutputCRAM
	}
	if binaryOutput {
		if err := validateBAMOutputFlags(opts); err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		emitBAMOutputWarnings(opts)
	}

	readersB := make([]io.Reader, 0, len(opts.inputB))
	for _, path := range opts.inputB {
		rb, err := iohelper.OpenReader(path)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error: Unable to open file %s. Exiting.\n", path)
			os.Exit(1)
		}
		defer rb.Close()
		readersB = append(readersB, rb)
	}

	writer, err := iohelper.OpenWriter(opts.output)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		os.Exit(1)
	}
	defer writer.Close()

	var genomeOrder []string
	if opts.genome != "" {
		genomeOrder, err = readGenomeOrder(opts.genome)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading genome file: %v\n", err)
			os.Exit(1)
		}
	}

	iopts := bedintersect.IntersectOptions{
		MinOverlap:      opts.minOverlap,
		FractionA:       opts.fractionA,
		FractionB:       opts.fractionB,
		StrandSpec:      opts.strand,
		ForceOpposite:   opts.opposite,
		NoOverlap:       opts.invert,
		WriteA:          opts.writeA,
		WriteB:          opts.writeB,
		Count:           opts.count,
		CountEach:       opts.countEach,
		Unique:          opts.unique,
		Reciprocal:      opts.reciprocal,
		EitherFraction:  opts.either,
		Distance:        opts.distance,
		Closest:         opts.closest,
		UseTree:         opts.useTree,
		LeftJoin:        opts.leftJoin,
		WriteOverlap:    opts.writeOverlap,
		WriteAllOverlap: opts.writeAllOverlap,
		Split:           opts.split,
		Names:           opts.names,
		FileNames:       opts.filenames,
		FilePaths:       opts.inputB,
		SortOut:         opts.sortOut,
		Header:          opts.header,
		Sorted:          opts.sorted,
		GenomeOrder:     genomeOrder,
		GenomeFile:      opts.genome,
		NoNameCheck:     opts.noNameCheck,
		NameA:           opts.inputA,
		NameB:           firstName(opts.inputB),
		Warnings:        os.Stderr,
	}

	if opts.stats {
		// --stats only supports a single B file in this port.
		stats, err := bedintersect.IntersectWithStats(readerA, readersB[0], writer, iopts)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error finding intersections: %v\n", err)
			os.Exit(1)
		}
		fmt.Fprintf(os.Stderr, "Intervals in A: %d\n", stats.IntervalsA)
		fmt.Fprintf(os.Stderr, "Intervals in B: %d\n", stats.IntervalsB)
		fmt.Fprintf(os.Stderr, "A intervals with hits: %d\n", stats.IntervalsAHit)
		fmt.Fprintf(os.Stderr, "A intervals with no hits: %d\n", stats.IntervalsAMiss)
		fmt.Fprintf(os.Stderr, "Total overlaps: %d\n", stats.Overlaps)
		return
	}

	if binaryOutput {
		// BAM/CRAM query without -bed: write the surviving alignments back out as
		// binary (BAM, or CRAM when a reference is available), matching upstream's
		// default behaviour. -header is ignored here (the alignment file carries
		// its own header); the emitBAMOutputWarnings call above has already issued
		// the upstream warning when -header/-wb/-loj were given.
		alnOut := bedintersect.AlnOutputOptions{Format: outputFormat, ReferenceFASTA: opts.cramRef, Uncompressed: opts.uncompressedBAM}
		if _, err := bedintersect.IntersectBinaryOutput(readerA, readersB, writer, iopts, alnOut); err != nil {
			fmt.Fprintf(os.Stderr, "Error finding intersections: %v\n", err)
			os.Exit(1)
		}
		return
	}

	if _, err := bedintersect.IntersectMulti(readerA, readersB, writer, iopts); err != nil {
		if bedintersect.IsVerbatimError(err) {
			// Upstream prints these verbatim (sort/field-count messages) and exits 1.
			fmt.Fprintln(os.Stderr, bedintersect.VerbatimMessage(err))
			os.Exit(1)
		}
		fmt.Fprintf(os.Stderr, "Error finding intersections: %v\n", err)
		os.Exit(1)
	}
}

// validateBAMOutputFlags reproduces upstream ContextIntersect::isValidState's
// gating for a BAM/CRAM query without -bed: the flags whose output is only
// meaningful as BED text are rejected, with the same ERROR banner upstream
// prints (and exit 1). Only -c (writeCount) and -wo/-wao (writeOverlap /
// writeAllOverlap) fall into this group. Notably -C (per-database counts) is
// NOT gated by upstream: with a BAM query it stays "printable" and so falls
// through to the default BAM-output rule (each A alignment with >=1 overlap),
// not a count. Everything else (-u/-v/-wa, default, and the warn-and-ignore
// -wb/-loj/-header) is allowed and produces BAM output.
func validateBAMOutputFlags(o *options) error {
	if o.count {
		return fmt.Errorf("***** ERROR: writeCount option is not valid with BAM query input, unless bed output is specified with -bed option. *****")
	}
	if o.writeOverlap || o.writeAllOverlap {
		return fmt.Errorf("***** ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option. *****")
	}
	return nil
}

// emitBAMOutputWarnings reproduces upstream's warn-and-ignore diagnostics for a
// BAM/CRAM query without -bed: -wb/-loj and -header are ignored (output is BAM),
// but upstream prints a stderr warning for each. The warnings are emitted
// verbatim so a `2>&1` capture matches upstream byte-for-byte.
func emitBAMOutputWarnings(o *options) {
	if o.writeB || o.leftJoin {
		fmt.Fprintln(os.Stderr, "\n*****\n*****WARNING: -wb and -loj are ignored with bam input, unless bed output is specified with -bed option.\n*****")
	}
	if o.header {
		fmt.Fprintln(os.Stderr, "\n*****\n*****WARNING: -header option is not valid for BAM input.\n*****")
	}
}

// boolFlags maps every boolean flag name (without the leading dash) to a pointer
// setter on the parsed options. Both single-dash and double-dash spellings are
// accepted, matching the upstream single-dash multi-character convention and the
// project's GNU long-flag convention.
func boolFlags(o *options) map[string]*bool {
	return map[string]*bool{
		"wa": &o.writeA, "write-a": &o.writeA,
		"wb": &o.writeB, "write-b": &o.writeB,
		"loj": &o.leftJoin, "left-outer-join": &o.leftJoin,
		"wo": &o.writeOverlap, "write-overlap": &o.writeOverlap,
		"wao": &o.writeAllOverlap, "write-all-overlap": &o.writeAllOverlap,
		"u": &o.unique, "unique": &o.unique,
		"c": &o.count, "count": &o.count,
		"C": &o.countEach,
		"v": &o.invert, "invert": &o.invert,
		"s": &o.strand, "strand": &o.strand,
		"S": &o.opposite,
		"r": &o.reciprocal, "reciprocal": &o.reciprocal,
		"e": &o.either, "either": &o.either,
		"split":       &o.split,
		"filenames":   &o.filenames,
		"sortout":     &o.sortOut,
		"header":      &o.header,
		"sorted":      &o.sorted,
		"nonamecheck": &o.noNameCheck,
		"bed":         &o.bedOutput,       // with BAM/CRAM query input, force BED text output
		"nobuf":       new(bool),          // accepted; output is buffered regardless
		"ubam":        &o.uncompressedBAM, // accepted; format follows the CRAM reference, BAM stays BGZF-compressed
		"d":           &o.distance, "distance": &o.distance,
		"k": &o.closest, "closest": &o.closest,
		"t": &o.useTree, "tree": &o.useTree,
		"stats": &o.stats,
		"h":     &o.help, "help": &o.help,
		"version": &o.showVersion,
	}
}

// parseArgs hand-parses the bedtools-style command line: single-dash
// multi-character flag names, "-b f1 f2 ..." multi-value flags, and "-names a
// b ..." alias lists. It returns the populated options or an error.
func parseArgs(args []string) (*options, error) {
	o := &options{minOverlap: 1}
	bools := boolFlags(o)

	i := 0
	for i < len(args) {
		arg := args[i]
		if arg == "" || arg == "-" || !strings.HasPrefix(arg, "-") {
			return nil, fmt.Errorf("unexpected argument %q", arg)
		}
		name := strings.TrimLeft(arg, "-")
		// Support --flag=value and -m=value spellings.
		var inlineVal string
		hasInline := false
		if eq := strings.IndexByte(name, '='); eq >= 0 {
			inlineVal = name[eq+1:]
			name = name[:eq]
			hasInline = true
		}

		// Multi-value collectors come first (-b and -names consume following
		// non-flag tokens).
		switch name {
		case "b":
			vals, next, err := collectValues(args, i, hasInline, inlineVal)
			if err != nil {
				return nil, fmt.Errorf("-b: %w", err)
			}
			o.inputB = append(o.inputB, vals...)
			i = next
			continue
		case "names":
			vals, next, err := collectValues(args, i, hasInline, inlineVal)
			if err != nil {
				return nil, fmt.Errorf("-names: %w", err)
			}
			o.names = append(o.names, vals...)
			i = next
			continue
		}

		// Single-value string/number flags.
		if ptr, ok := stringTargets(o)[name]; ok {
			val, next, err := singleValue(args, i, hasInline, inlineVal)
			if err != nil {
				return nil, fmt.Errorf("-%s: %w", name, err)
			}
			*ptr = val
			i = next
			continue
		}
		switch name {
		case "m", "min-overlap":
			val, next, err := singleValue(args, i, hasInline, inlineVal)
			if err != nil {
				return nil, fmt.Errorf("-%s: %w", name, err)
			}
			n, perr := strconv.Atoi(val)
			if perr != nil {
				return nil, fmt.Errorf("-%s: invalid integer %q", name, val)
			}
			o.minOverlap = n
			i = next
			continue
		case "f", "fraction-a", "F", "fraction-b":
			val, next, err := singleValue(args, i, hasInline, inlineVal)
			if err != nil {
				return nil, fmt.Errorf("-%s: %w", name, err)
			}
			fv, perr := strconv.ParseFloat(val, 64)
			if perr != nil {
				return nil, fmt.Errorf("-%s: invalid float %q", name, val)
			}
			if name == "f" || name == "fraction-a" {
				if fv <= 0.0 || fv > 1.0 {
					return nil, &rangeError{"-f"}
				}
				o.fractionA = fv
			} else {
				if fv <= 0.0 || fv > 1.0 {
					return nil, &rangeError{"-F"}
				}
				o.fractionB = fv
			}
			i = next
			continue
		case "iobuf":
			// Accepted for compatibility; consumes a value but has no effect.
			_, next, err := singleValue(args, i, hasInline, inlineVal)
			if err != nil {
				return nil, fmt.Errorf("-iobuf: %w", err)
			}
			i = next
			continue
		}

		// Boolean flags.
		if ptr, ok := bools[name]; ok {
			if hasInline {
				return nil, fmt.Errorf("-%s takes no value", name)
			}
			*ptr = true
			i++
			continue
		}

		return nil, &unrecognizedParamError{arg}
	}
	return o, nil
}

// unrecognizedParamError mirrors upstream's exact "Unrecognized parameter"
// banner so an unknown flag produces the same stderr line.
type unrecognizedParamError struct{ param string }

func (e *unrecognizedParamError) Error() string {
	return fmt.Sprintf("***** ERROR: Unrecognized parameter: %s *****", e.param)
}

// rangeError mirrors upstream's exact banner for an out-of-range -f/-F value,
// which must lie in (0.0, 1.0].
type rangeError struct{ flag string }

func (e *rangeError) Error() string {
	return fmt.Sprintf("***** ERROR: %s must be in the range (0.0, 1.0]. *****", e.flag)
}

// stringTargets maps the value-taking string flags to their option fields.
func stringTargets(o *options) map[string]*string {
	return map[string]*string{
		"a": &o.inputA, "input-a": &o.inputA,
		"abam": &o.inputA, "ibam": &o.inputA,
		"o": &o.output, "output": &o.output,
		"g": &o.genome, "genome": &o.genome,
		"cram-ref": &o.cramRef, // CRAM reference FASTA (upstream's global --cram-ref)
	}
}

// singleValue returns the value for a flag at args[i], either inline (-f=0.5) or
// as the next token (-f 0.5), and the index of the next unprocessed argument.
func singleValue(args []string, i int, hasInline bool, inlineVal string) (string, int, error) {
	if hasInline {
		return inlineVal, i + 1, nil
	}
	if i+1 >= len(args) {
		return "", 0, fmt.Errorf("missing value")
	}
	return args[i+1], i + 2, nil
}

// normalizeStdin maps upstream's "stdin" alias (and "/dev/stdin") to iohelper's
// "-" stdin marker, leaving every other path unchanged.
func normalizeStdin(path string) string {
	if path == "stdin" || path == "/dev/stdin" {
		return "-"
	}
	return path
}

// firstName returns the first element of names, or "" when empty. Used to label
// the (first) B file in -sorted error messages.
func firstName(names []string) string {
	if len(names) == 0 {
		return ""
	}
	return names[0]
}

// readGenomeOrder reads a genome file (chrom<TAB>size per line) and returns the
// chromosome names in file order, which fixes the required chromosome order for
// -sorted -g validation. Blank and comment lines are skipped.
func readGenomeOrder(path string) ([]string, error) {
	r, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, err
	}
	defer r.Close()
	data, err := io.ReadAll(r)
	if err != nil {
		return nil, err
	}
	var order []string
	for _, line := range strings.Split(string(data), "\n") {
		line = strings.TrimRight(line, "\r")
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.SplitN(line, "\t", 2)
		order = append(order, fields[0])
	}
	return order, nil
}

// collectValues gathers one or more values for a multi-value flag (-b, -names):
// every following token until the next flag (a token starting with '-', other
// than the lone "-" stdin marker). An inline value (=v) contributes a single
// value and stops there.
func collectValues(args []string, i int, hasInline bool, inlineVal string) ([]string, int, error) {
	if hasInline {
		return []string{inlineVal}, i + 1, nil
	}
	var vals []string
	j := i + 1
	for j < len(args) {
		tok := args[j]
		if tok == "-" || !strings.HasPrefix(tok, "-") {
			vals = append(vals, tok)
			j++
			continue
		}
		break
	}
	if len(vals) == 0 {
		return nil, 0, fmt.Errorf("missing value")
	}
	return vals, j, nil
}
