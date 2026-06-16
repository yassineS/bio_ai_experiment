// Package bedmerge provides functionality to merge overlapping or adjacent BED
// intervals, mirroring `bedtools merge`. It accepts BED, GFF, VCF, and BAM
// input (auto-detected), supports the -c/-o column-operation family, strand
// modes (-s/-S), the merge distance (-d), a custom list delimiter (-delim) and
// numeric output precision (-prec).
package bedmerge

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
)

// DefaultPrecision is the number of significant digits used when formatting
// floating-point column-operation results, matching upstream bedtools'
// KeyListOps DEFAULT_PRECISION (10). The -prec/--precision flag overrides it.
const DefaultPrecision = 10

// OutputFields specifies which fields to include in the output.
type OutputFields struct {
	Name     bool // Include name field
	Score    bool // Include score field
	Strand   bool // Include strand field
	Count    bool // Include count of merged intervals as name field
	BedGraph bool // Output in bedGraph format (chrom, start, end, score)
}

// MergeOptions contains options for merging BED intervals.
type MergeOptions struct {
	MaxDistance  int          // Maximum distance between intervals to merge (default: 0)
	StrandSpec   bool         // Merge only intervals on the same strand (-s)
	OutputFields OutputFields // Fields to include in output (legacy convenience flags)
	// StrandFilter ("-S <strand>"), when non-empty, must be "+" or "-": only
	// records on that strand are kept, and the survivors are merged positionally
	// (NOT per-strand), matching upstream bedtools merge -S. It is mutually
	// exclusive with StrandSpec ("-s").
	StrandFilter string
	Streaming    bool // Use streaming mode for large files (reserved; see note below)
	// ColumnOps, when non-nil, requests bedtools-merge-style aggregation of input
	// columns (the -c/-o options). When set, output columns are chrom, start, end
	// followed by one aggregated value per requested column.
	ColumnOps *ColumnOps
	// Precision is the number of significant digits for floating-point
	// column-op output (-prec). Zero means DefaultPrecision.
	Precision int
	// Warn, when set, receives upstream-style stderr warnings (e.g. the
	// "Non numeric value" message). When nil, warnings are discarded.
	Warn io.Writer
}

// precision returns the effective output precision (DefaultPrecision when 0).
func (o MergeOptions) precision() int {
	if o.Precision <= 0 {
		return DefaultPrecision
	}
	return o.Precision
}

// Merge reads intervals (BED/GFF/VCF/BAM, auto-detected), sorts them, merges
// overlapping/adjacent intervals, and writes the result. Returns the number of
// merged intervals written.
func Merge(reader io.Reader, writer io.Writer, opts MergeOptions) (int, error) {
	if err := validateStrandOptions(opts); err != nil {
		return 0, err
	}
	recs, format, err := readInput(reader, opts)
	if err != nil {
		return 0, err
	}
	if err := validateBAMColumns(format, opts.ColumnOps); err != nil {
		return 0, err
	}
	return mergeRecords(recs, writer, opts)
}

// bamNumFields is the number of SAM fields a BAM record exposes to -c column
// operations (QNAME..QUAL), matching upstream BamFileReader::getNumFields.
const bamNumFields = 11

// BAMColumnError reports a -c column request that is invalid for BAM input: a
// column outside 1..11, or column 2 (the unsupported Flags field). The CLI
// formats the matching upstream message, which includes the input file name.
type BAMColumnError struct {
	Column int
	Flags  bool // true when the column is 2 (the Flags field)
}

func (e *BAMColumnError) Error() string {
	if e.Flags {
		return "requested column 2 of a BAM file, which is the Flags field"
	}
	return fmt.Sprintf("requested column %d, but BAM input only has fields 1 - %d", e.Column, bamNumFields)
}

// validateBAMColumns enforces the BAM-specific -c column constraints (upstream
// KeyListOps::init): a column must be within 1..11 and column 2 (Flags) is
// unsupported. It is a no-op for non-BAM input or when no column ops are set.
func validateBAMColumns(format inputFormat, co *ColumnOps) error {
	if format != fmtBAM || co == nil {
		return nil
	}
	for _, col := range co.Columns {
		if col < 1 || col > bamNumFields {
			return &BAMColumnError{Column: col}
		}
	}
	for _, col := range co.Columns {
		if col == 2 {
			return &BAMColumnError{Column: 2, Flags: true}
		}
	}
	return nil
}

// mergeRecords runs the strand bucketing, position merge, and output for an
// already-parsed record slice. It is shared by Merge and MergeWithStats.
func mergeRecords(recs []record, writer io.Writer, opts MergeOptions) (int, error) {
	groups := buildGroups(recs, opts)
	bw := bufio.NewWriter(writer)
	count := 0
	for _, g := range groups {
		if err := writeGroup(bw, g, opts); err != nil {
			return 0, err
		}
		count++
	}
	if err := bw.Flush(); err != nil {
		return 0, fmt.Errorf("error flushing output: %w", err)
	}
	return count, nil
}

// buildGroups filters by strand, sorts, and merges into groups of records that
// collapse into one output interval. Under -s, '+' and '-' are merged
// independently and the two streams re-merged by (chrom, start, end);
// unknown-strand ('.'/”) records are dropped. Under -S <strand>, only that
// strand survives and is merged positionally.
func buildGroups(recs []record, opts MergeOptions) [][]record {
	if opts.StrandFilter != "" {
		kept := recs[:0:0]
		for _, r := range recs {
			if r.strand == opts.StrandFilter {
				kept = append(kept, r)
			}
		}
		recs = kept
	}

	if opts.StrandSpec {
		var plus, minus []record
		for _, r := range recs {
			switch r.strand {
			case "+":
				plus = append(plus, r)
			case "-":
				minus = append(minus, r)
				// "." and "" are dropped (upstream FileRecordMergeMgr behaviour).
			}
		}
		sortRecords(plus)
		sortRecords(minus)
		pg := positionGroups(plus, opts)
		mg := positionGroups(minus, opts)
		return mergeSortedGroups(pg, mg)
	}

	sortRecords(recs)
	return positionGroups(recs, opts)
}

// positionGroups runs a single-pass position-only merge over sorted records
// (strand is ignored; the caller handles any strand bucketing). Each output
// group is the slice of input records that merge together.
func positionGroups(recs []record, opts MergeOptions) [][]record {
	if len(recs) == 0 {
		return nil
	}
	var out [][]record
	group := []record{recs[0]}
	curChrom := recs[0].chrom
	curEnd := recs[0].end
	for i := 1; i < len(recs); i++ {
		r := recs[i]
		if r.chrom == curChrom && r.start <= curEnd+opts.MaxDistance {
			if r.end > curEnd {
				curEnd = r.end
			}
			group = append(group, r)
			continue
		}
		out = append(out, group)
		group = []record{r}
		curChrom = r.chrom
		curEnd = r.end
	}
	out = append(out, group)
	return out
}

// mergeSortedGroups two-way-merges two slices of already-sorted groups into a
// single stream ordered by (chrom, start, end). Used to recombine the '+' and
// '-' outputs of a stranded merge.
func mergeSortedGroups(a, b [][]record) [][]record {
	out := make([][]record, 0, len(a)+len(b))
	i, j := 0, 0
	groupEnd := func(g []record) int {
		e := g[0].end
		for _, r := range g[1:] {
			if r.end > e {
				e = r.end
			}
		}
		return e
	}
	less := func(x, y []record) bool {
		if x[0].chrom != y[0].chrom {
			return x[0].chrom < y[0].chrom
		}
		if x[0].start != y[0].start {
			return x[0].start < y[0].start
		}
		return groupEnd(x) < groupEnd(y)
	}
	for i < len(a) && j < len(b) {
		if less(a[i], b[j]) {
			out = append(out, a[i])
			i++
		} else {
			out = append(out, b[j])
			j++
		}
	}
	out = append(out, a[i:]...)
	out = append(out, b[j:]...)
	return out
}

// writeGroup writes one merged output line for a group of records. With
// ColumnOps it emits chrom/start/end plus one aggregated value per column;
// otherwise it honours the legacy OutputFields convenience flags (count,
// bedGraph) and defaults to BED3.
func writeGroup(w *bufio.Writer, group []record, opts MergeOptions) error {
	chrom := group[0].chrom
	start := group[0].start
	end := group[0].end
	for _, r := range group {
		if r.end > end {
			end = r.end
		}
	}

	if opts.ColumnOps != nil {
		out := []string{chrom, strconv.Itoa(start), strconv.Itoa(end)}
		vals, warn := applyColumnOps(opts.ColumnOps, group, opts.precision())
		out = append(out, vals...)
		if warn != "" && opts.Warn != nil {
			fmt.Fprintln(opts.Warn, warn)
		}
		_, err := fmt.Fprintln(w, joinTab(out))
		return err
	}

	first := group[0].fields

	// Legacy bedmerge convenience modes (not part of upstream bedtools merge,
	// but retained as bedmerge extras). bedGraph emits the 4-col score (column
	// 4); the Name/Score/Strand flags echo those columns from the first record.
	if opts.OutputFields.BedGraph {
		score := "0"
		if len(first) > 3 {
			score = first[3]
		}
		_, err := fmt.Fprintf(w, "%s\t%d\t%d\t%s\n", chrom, start, end, score)
		return err
	}

	if opts.OutputFields.Name || opts.OutputFields.Score || opts.OutputFields.Strand {
		out := []string{chrom, strconv.Itoa(start), strconv.Itoa(end)}
		if opts.OutputFields.Name {
			out = append(out, fieldOr(first, 3, ""))
		}
		if opts.OutputFields.Score {
			out = append(out, fieldOr(first, 4, "0"))
		}
		if opts.OutputFields.Strand {
			out = append(out, fieldOr(first, 5, "."))
		}
		_, err := fmt.Fprintln(w, joinTab(out))
		return err
	}

	if opts.OutputFields.Count {
		_, err := fmt.Fprintf(w, "%s\t%d\t%d\t%d\n", chrom, start, end, len(group))
		return err
	}

	_, err := fmt.Fprintf(w, "%s\t%d\t%d\n", chrom, start, end)
	return err
}

// fieldOr returns the value at 0-based index i of fields, or def when absent.
func fieldOr(fields []string, i int, def string) string {
	if i < len(fields) {
		return fields[i]
	}
	return def
}

// joinTab joins fields with a tab. Inline helper for the column-op hot path.
func joinTab(fields []string) string {
	if len(fields) == 0 {
		return ""
	}
	n := len(fields) - 1
	for _, f := range fields {
		n += len(f)
	}
	b := make([]byte, 0, n)
	for i, f := range fields {
		if i > 0 {
			b = append(b, '\t')
		}
		b = append(b, f...)
	}
	return string(b)
}

// validateStrandOptions rejects an invalid -S argument and the illegal
// combination of -s and -S, mirroring upstream bedtools merge.
func validateStrandOptions(opts MergeOptions) error {
	if opts.StrandFilter == "" {
		return nil
	}
	if opts.StrandFilter != "+" && opts.StrandFilter != "-" {
		return errBadStrandArg
	}
	if opts.StrandSpec {
		return fmt.Errorf("-s and -S are mutually exclusive")
	}
	return nil
}

// ErrBadStrandArg is the sentinel for an invalid -S argument; the CLI prints the
// upstream-formatted message ("-S option must be followed by + or -").
var ErrBadStrandArg = errBadStrandArg

var errBadStrandArg = fmt.Errorf("-S option must be followed by + or -")

// ErrStrandedVCF is the sentinel for the unsupported "-s with VCF" combination;
// the CLI prints the upstream-formatted message including the file name.
var ErrStrandedVCF = errStrandedVCF

// Stats contains statistics about the merge operation.
type Stats struct {
	InputIntervals  int
	OutputIntervals int
	MergedCount     int
}

// MergeWithStats performs a merge and returns detailed statistics.
func MergeWithStats(reader io.Reader, writer io.Writer, opts MergeOptions) (*Stats, error) {
	if err := validateStrandOptions(opts); err != nil {
		return nil, err
	}
	recs, format, err := readInput(reader, opts)
	if err != nil {
		return nil, err
	}
	if err := validateBAMColumns(format, opts.ColumnOps); err != nil {
		return nil, err
	}
	stats := &Stats{InputIntervals: len(recs)}
	out, err := mergeRecords(recs, writer, opts)
	if err != nil {
		return nil, err
	}
	stats.OutputIntervals = out
	stats.MergedCount = stats.InputIntervals - stats.OutputIntervals
	return stats, nil
}
