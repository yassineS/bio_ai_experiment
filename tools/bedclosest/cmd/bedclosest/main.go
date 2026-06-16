// bedclosest finds, for each A interval, the closest B interval (mirrors
// `bedtools closest`).
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/bedclosest/pkg/bedclosest"
)

const version = "1.0.0"

const usage = `bedclosest - Find the closest B interval for each A interval

Usage:
  bedclosest -a <fileA.bed> -b <fileB.bed> [options]

Description:
  For each interval in A (sorted), find the closest interval in B (also
  sorted) on the same chromosome and report A's columns + B's columns.
  With -d the unsigned distance is appended; with -D the signed distance
  is appended. Distance is 0 when A and B overlap. For tied distances the
  default is to emit one row per tied B (-t all).

  Both inputs MUST be sorted on (chrom, start). bedclosest errors out
  clearly when they are not. Output matches bedtools closest byte-for-byte.

Options:
  -a, --a FILE             Input BED file A (sorted; use '-' for stdin)
  -b, --b FILE...          One or more sorted BED database files (use '-' for
                           stdin). With multiple files a database-label column
                           (the 1-based file index, or -names/-filenames) is
                           inserted between A's and B's columns.
  -names NAME...           Labels for the -b databases (one per file, in order);
                           replaces the numeric file-index column. Mutually
                           exclusive with -filenames.
  -filenames               Use each -b file's name as its database-label column.
  -mdb each|all            Multi-database mode: 'each' (default) reports the
                           closest feature from every database on its own row;
                           'all' reports the single overall closest across all
                           databases.
  -o, --output FILE        Output BED file ('-' for stdout, default: stdout)
  -d                       Append the (unsigned) distance to the closest B.
  -D MODE                  Append the signed distance, with the sign convention:
                           ref (downstream of A positive), a (relative to A's
                           strand), or b (relative to B's strand).
  -k N                     Report the N closest hits per A interval (default 1).
  -N                       Require the closest B to have a different name
                           (BED column 4) than A.
  -io                      Ignore B features that overlap A.
  -iu                      Ignore B features upstream of A (requires -D).
  -id                      Ignore B features downstream of A (requires -D).
  -fu                      Force the closest upstream B feature (requires -D).
  -fd                      Force the closest downstream B feature (requires -D).
  -s                       Require the closest B to be on the SAME strand as A.
  -S                       Require the closest B to be on the OPPOSITE strand.
                           -s and -S are mutually exclusive.
  -t MODE                  Tie-break mode for equally close B's:
                             all   - emit one row per tied B (default)
                             first - emit only the first (in B's input order)
                             last  - emit only the last
  -header                  Print A's header lines before the results.
  -h, --help               Show this help message and exit
  -v, --version            Show version information and exit

Examples:
  # Closest peak for each gene (both sorted)
  bedclosest -a genes.sorted.bed -b peaks.sorted.bed > out.bed

  # Append the distance column
  bedclosest -a a.bed -b b.bed -d > out.bed

  # Report the 5 closest features
  bedclosest -a a.bed -b b.bed -k 5 > out.bed

  # Single hit per A (first in B input order on ties)
  bedclosest -a a.bed -b b.bed -t first > out.bed

Format:
  Input: BED format (tab-delimited, minimum 3 columns), sorted on
         (chrom, start).
  Output: A's columns, then B's columns, then distance (if -d or -D).
`

func main() {
	fs := flag.CommandLine
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	// -b and -names are variadic (`-b f1 f2 f3`, `-names a b c`), mirroring
	// upstream's space-separated multi-value syntax. They are pulled out of the
	// argument list before the remaining flags are parsed: everything up to the
	// next dash-prefixed token belongs to the flag.
	argv := os.Args[1:]
	bFiles, argv := extractVarArg(argv, []string{"-b", "--b"})
	nameLabels, argv := extractVarArg(argv, []string{"-names", "--names"})

	var aFile, outputFile string
	// "-a" and "--a" are treated identically by Go's flag package; register once.
	fs.StringVar(&aFile, "a", "", "Input BED file A (sorted)")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output BED file")

	var useFilenames bool
	fs.BoolVar(&useFilenames, "filenames", false, "Use each -b file's name as its database-label column")

	var mdbMode string
	fs.StringVar(&mdbMode, "mdb", "each", "Multi-database mode: each|all")

	// -d appends the unsigned distance; -D appends the signed distance with a
	// chosen sign convention. They are independent of each other.
	var printAbsDist bool
	fs.BoolVar(&printAbsDist, "d", false, "Append the unsigned distance to the closest B")

	var distMode string
	fs.StringVar(&distMode, "D", "", "Append the signed distance: ref|a|b")

	var kClosest int
	fs.IntVar(&kClosest, "k", 1, "Report the N closest hits per A interval")

	var differentNames bool
	fs.BoolVar(&differentNames, "N", false, "Require the closest B to have a different name (column 4) than A")

	var ignoreOverlaps, ignoreUp, ignoreDown, forceUp, forceDown bool
	fs.BoolVar(&ignoreOverlaps, "io", false, "Ignore features in B that overlap A")
	fs.BoolVar(&ignoreUp, "iu", false, "Ignore features in B that are upstream of A (requires -D)")
	fs.BoolVar(&ignoreDown, "id", false, "Ignore features in B that are downstream of A (requires -D)")
	fs.BoolVar(&forceUp, "fu", false, "Force the closest upstream feature in B (requires -D)")
	fs.BoolVar(&forceDown, "fd", false, "Force the closest downstream feature in B (requires -D)")

	var tieMode string
	fs.StringVar(&tieMode, "t", "all", "Tie-break: all|first|last")

	// Strand filters mirror upstream's bare -s/-S (no long form upstream).
	var sameStrand, oppositeStrand bool
	fs.BoolVar(&sameStrand, "s", false, "Require the closest B to be on the same strand as A")
	fs.BoolVar(&oppositeStrand, "S", false, "Require the closest B to be on the opposite strand to A")

	var printHeader bool
	fs.BoolVar(&printHeader, "header", false, "Print A's header lines before the results")

	var help, showVersion bool
	cliflag.BoolVar(fs, &help, "h", "help", false, "Show help message")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version information")

	if err := fs.Parse(argv); err != nil {
		os.Exit(2)
	}

	if help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}
	if showVersion {
		fmt.Printf("bedclosest version %s\n", version)
		os.Exit(0)
	}

	if aFile == "" || len(bFiles) == 0 {
		fmt.Fprintln(os.Stderr, "Error: -a and -b are required")
		os.Exit(2)
	}
	stdinCount := 0
	if aFile == "-" {
		stdinCount++
	}
	for _, b := range bFiles {
		if b == "-" {
			stdinCount++
		}
	}
	if stdinCount > 1 {
		fmt.Fprintln(os.Stderr, "Error: at most one input may be '-' (stdin)")
		os.Exit(2)
	}

	// -d and -D are mutually exclusive (matching upstream). The distance mode
	// is DistanceRef (unsigned) for -d, or the parsed signed mode for -D.
	dGiven := false
	dDGiven := false
	fs.Visit(func(f *flag.Flag) {
		switch f.Name {
		case "d":
			dGiven = true
		case "D":
			dDGiven = true
		}
	})
	if dGiven && dDGiven {
		fmt.Fprintln(os.Stderr, "Error: Request either -d OR -D, not both.")
		os.Exit(2)
	}
	reportDistance := dGiven || dDGiven
	dm := bedclosest.DistanceRef
	if dDGiven {
		parsed, err := parseDistanceMode(distMode)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v\n", err)
			os.Exit(2)
		}
		dm = parsed
	}

	tb, err := parseTieBreak(tieMode)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(2)
	}
	mm, err := parseMultiDBMode(mdbMode)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(2)
	}

	// -names and -filenames are mutually exclusive and, when set, must label
	// each -b database, mirroring upstream's checks.
	if len(nameLabels) > 0 && useFilenames {
		fmt.Fprintln(os.Stderr, "Error: -names and -filenames are mutually exclusive.")
		os.Exit(2)
	}
	var dbLabels []string
	switch {
	case len(nameLabels) > 0:
		if len(nameLabels) != len(bFiles) {
			fmt.Fprintf(os.Stderr, "Error: the number of -names (%d) does not match the number of -b files (%d)\n",
				len(nameLabels), len(bFiles))
			os.Exit(2)
		}
		dbLabels = nameLabels
	case useFilenames:
		dbLabels = append([]string(nil), bFiles...)
	}

	// Validate the directional flags against upstream's ContextClosest rules.
	if ignoreUp && ignoreDown {
		fmt.Fprintln(os.Stderr, "Error: Request either -iu OR -id, not both.")
		os.Exit(2)
	}
	if forceUp && forceDown {
		fmt.Fprintln(os.Stderr, "Error: Request either -fu OR -fd, not both.")
		os.Exit(2)
	}
	if ignoreUp && forceUp {
		fmt.Fprintln(os.Stderr, "Error: Request either -iu OR -fu, not both.")
		os.Exit(2)
	}
	if ignoreDown && forceDown {
		fmt.Fprintln(os.Stderr, "Error: Request either -id OR -fd, not both.")
		os.Exit(2)
	}
	// -iu/-id/-fu/-fd require an explicit stranded distance mode (-D).
	if (ignoreUp || ignoreDown || forceUp || forceDown) && !dDGiven {
		fmt.Fprintln(os.Stderr, "Error: When requesting -iu, -id, -fu, or -fd, you also need to specify -D.")
		os.Exit(2)
	}
	if sameStrand && oppositeStrand {
		fmt.Fprintln(os.Stderr, "Error: Request either -s OR -S, not both.")
		os.Exit(2)
	}

	readerA, err := iohelper.OpenReader(aFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening A: %v\n", err)
		os.Exit(1)
	}
	defer readerA.Close()

	readersB := make([]io.Reader, 0, len(bFiles))
	var bClosers []io.Closer
	defer func() {
		for _, c := range bClosers {
			c.Close()
		}
	}()
	for _, b := range bFiles {
		rb, err := iohelper.OpenReader(b)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening B %q: %v\n", b, err)
			os.Exit(1)
		}
		bClosers = append(bClosers, rb)
		readersB = append(readersB, rb)
	}

	writer, err := iohelper.OpenWriter(outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening output: %v\n", err)
		os.Exit(1)
	}
	defer writer.Close()

	opts := bedclosest.Options{
		ReportDistance:   reportDistance,
		DistanceMode:     dm,
		KClosest:         kClosest,
		TieBreak:         tb,
		IgnoreOverlaps:   ignoreOverlaps,
		IgnoreUpstream:   ignoreUp,
		IgnoreDownstream: ignoreDown,
		ForceUpstream:    forceUp,
		ForceDownstream:  forceDown,
		SameStrand:       sameStrand,
		OppositeStrand:   oppositeStrand,
		DifferentNames:   differentNames,
		MultiDBMode:      mm,
		DBLabels:         dbLabels,
		PrintHeader:      printHeader,
		// Enable upstream's cross-file chromosome sort-order / naming-convention
		// validation, writing its WARNING/ERROR text to stderr. The file-name
		// labels are the names exactly as given on the command line, matching
		// upstream's getInputFileName.
		WarnWriter: os.Stderr,
		QueryName:  aFile,
		DBNames:    bFiles,
	}

	if _, err := bedclosest.ClosestMulti(readerA, readersB, writer, opts); err != nil {
		// A validation ERROR has already written its exact upstream text to
		// stderr; just propagate the non-zero exit. Other errors get the
		// generic prefix.
		if !bedclosest.IsValidationError(err) {
			fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		}
		os.Exit(1)
	}
}

// extractVarArg pulls the values following any of the given trigger flags until
// the next dash-prefixed token (or end of input), returning the collected
// values and the remaining arguments. It implements upstream's space-separated
// multi-value flag syntax (e.g. `-b a.bed b.bed`, `-names a b c`).
func extractVarArg(argv, triggers []string) (values, rest []string) {
	rest = make([]string, 0, len(argv))
	for i := 0; i < len(argv); i++ {
		a := argv[i]
		matched := false
		for _, t := range triggers {
			if a == t {
				matched = true
				break
			}
		}
		if !matched {
			rest = append(rest, a)
			continue
		}
		i++
		for ; i < len(argv); i++ {
			v := argv[i]
			if len(v) > 0 && v[0] == '-' && v != "-" {
				i--
				break
			}
			values = append(values, v)
		}
	}
	return values, rest
}

// parseMultiDBMode converts the -mdb flag string into a bedclosest.MultiDBMode.
func parseMultiDBMode(s string) (bedclosest.MultiDBMode, error) {
	switch s {
	case "each":
		return bedclosest.MultiDBEach, nil
	case "all":
		return bedclosest.MultiDBAll, nil
	default:
		return 0, fmt.Errorf("invalid -mdb value %q (expected each|all)", s)
	}
}

// parseDistanceMode converts the -D flag string into a signed bedclosest
// DistanceMode.
func parseDistanceMode(s string) (bedclosest.DistanceMode, error) {
	switch s {
	case "ref":
		return bedclosest.DistanceSignedRef, nil
	case "a":
		return bedclosest.DistanceA, nil
	case "b":
		return bedclosest.DistanceB, nil
	default:
		return 0, fmt.Errorf("invalid -D value %q (expected ref|a|b)", s)
	}
}

// parseTieBreak converts the -t flag string into a bedclosest.TieBreak.
func parseTieBreak(s string) (bedclosest.TieBreak, error) {
	switch s {
	case "all":
		return bedclosest.TieAll, nil
	case "first":
		return bedclosest.TieFirst, nil
	case "last":
		return bedclosest.TieLast, nil
	default:
		return 0, fmt.Errorf("invalid -t value %q (expected all|first|last)", s)
	}
}
