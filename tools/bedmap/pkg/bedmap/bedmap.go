// Package bedmap applies column aggregation ops to A's intervals using B as
// the value source, mirroring upstream `bedtools map`.
//
// For each interval in A the package finds the overlapping B records, takes
// the values from one or more B columns, and runs the requested aggregation
// (e.g. sum, mean, count, distinct, collapse, ...). The aggregated values are
// appended to A's columns. When no B records overlap, the configured Null
// string is emitted instead (default ".").
//
// B may be BED, VCF, GFF, or BAM; the format is auto-detected exactly as
// upstream's BedFile::parseLine (see input.go). Column-op semantics for the
// list/count family come from bedmerge.ApplyOp; the numeric family (sum, mean,
// min, max, absmin, absmax, median, stdev, sstdev) is computed locally so that
// non-numeric values produce upstream's WARNING + null behaviour and results
// are formatted with upstream's precision (default 10 significant digits).
package bedmap

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	"github.com/yassineS/bio_ai_experiment/tools/bedmerge/pkg/bedmerge"
)

// Options controls Map.
type Options struct {
	// Columns is the 1-based list of B columns to aggregate. Default {5}
	// (the score column) matches upstream `bedtools map`.
	Columns []int
	// Ops is one op per column, OR a single op applied to all columns. Default
	// {"sum"} matches upstream.
	Ops []string

	// Null is the placeholder emitted when A has no overlapping B record (and
	// when a numeric op produces NaN). Default ".".
	Null string
	// Delim is the separator used by collapse / distinct (and other
	// "concatenate" ops). Default ",".
	Delim string
	// Precision is the number of significant digits used to format numeric
	// op results, matching upstream's -prec (default 10).
	Precision int

	// SameStrand: only consider B records on the same strand as A.
	SameStrand bool
	// OppositeStrand: only consider B records on the opposite strand.
	OppositeStrand bool

	// FractionA: minimum fraction of A that must overlap a single B record
	// before that B contributes. 0 disables the check.
	FractionA float64
	// FractionB: minimum fraction of B that must overlap A.
	FractionB float64
	// Reciprocal: when true with -f and -F, both thresholds must hold.
	Reciprocal bool

	// Split: when true, treat BED12/BAM A and B records as their constituent
	// blocks for overlap detection (only blocks overlap blocks), matching
	// upstream `bedtools map -split`.
	Split bool

	// Header: when true, echo A's leading comment/track/browser header lines
	// to the output verbatim, matching upstream `bedtools map -header`.
	Header bool

	// BFileName is the name of the B file as the user supplied it. It is used
	// only to reproduce upstream's column-range error text
	// ("... database file <name> only has fields 1 - N."). When empty, "B"
	// is used.
	BFileName string

	// WarnWriter receives the per-row "WARNING: Non numeric value ..." lines
	// upstream prints to stderr. When nil, warnings are discarded.
	WarnWriter io.Writer
}

// Validate fills in defaults and rejects obviously bad combinations. It
// mutates opts in place.
func (opts *Options) Validate() error {
	if len(opts.Columns) == 0 {
		opts.Columns = []int{5}
	}
	// Out-of-range columns (including <= 0) are reported by Map against the
	// database's field count, matching upstream's column-range error text, so
	// they are intentionally not rejected here.
	if len(opts.Ops) == 0 {
		opts.Ops = []string{"sum"}
	}
	switch {
	case len(opts.Ops) == len(opts.Columns):
		// ok
	case len(opts.Ops) == 1:
		one := opts.Ops[0]
		opts.Ops = make([]string, len(opts.Columns))
		for i := range opts.Ops {
			opts.Ops[i] = one
		}
	case len(opts.Columns) == 1:
		// One column, many ops: replicate the column for each op.
		col := opts.Columns[0]
		opts.Columns = make([]int, len(opts.Ops))
		for i := range opts.Columns {
			opts.Columns[i] = col
		}
	default:
		// Match upstream's exact stderr text for the columns/operations
		// mismatch (KeyListOps context validation), including its three-line
		// guidance suffix.
		return fmt.Errorf("\n*****\n***** ERROR: There are %d columns given, but there are %d operations.\n"+
			"\tPlease provide either a single operation that will be applied to all listed columns, \n"+
			"\ta single column to which all operations will be applied,\n"+
			"\tor an operation for each column.", len(opts.Columns), len(opts.Ops))
	}
	if opts.Null == "" {
		opts.Null = "."
	}
	if opts.Delim == "" {
		opts.Delim = ","
	}
	if opts.Precision == 0 {
		opts.Precision = 10
	}
	if opts.SameStrand && opts.OppositeStrand {
		return fmt.Errorf("-s and -S are mutually exclusive")
	}
	return nil
}

// rawRecord pairs a parsed BED record with its raw whitespace-split columns,
// so column-value extraction by 1-based index works on the original text. order
// is the record's 0-based position in B-file load order, used to break ties on
// equal (chrom, start) so the order-sensitive ops (collapse, distinct) emit
// values in upstream's input order rather than by chromEnd.
type rawRecord struct {
	rec    *bed.Record
	fields []string
	order  int
}

// numericOps lists the operations that convert their column values to numbers,
// warn on non-numeric input, and print the null value on a NaN result, matching
// upstream KeyListOps. mode / antimode operate on string frequency maps and are
// therefore excluded here (delegated to bedmerge.ApplyOp).
var numericOps = map[string]bool{
	"sum":    true,
	"mean":   true,
	"min":    true,
	"max":    true,
	"absmin": true,
	"absmax": true,
	"median": true,
	"stdev":  true,
	"sstdev": true,
}

// Map runs the per-A column aggregation against B and writes the results to
// writer. Returns the number of A records processed.
func Map(readerA, readerB io.Reader, writer io.Writer, opts Options) (int, error) {
	if err := opts.Validate(); err != nil {
		return 0, err
	}

	// Read B with format auto-detection (BED/VCF/GFF/BAM).
	bRaw, maxFields, err := readBRecords(readerB)
	if err != nil {
		return 0, fmt.Errorf("error reading B: %w", err)
	}

	// Validate the requested columns against the database's field count,
	// reproducing upstream's exact error text. Upstream reports the maximum
	// number of fields it saw across the file.
	bName := opts.BFileName
	if bName == "" {
		bName = "B"
	}
	for _, c := range opts.Columns {
		if c < 1 || c > maxFields {
			return 0, fmt.Errorf("\n*****\n***** ERROR: Requested column %d, but database file %s only has fields 1 - %d.", c, bName, maxFields)
		}
	}

	bByChrom := map[string][]rawRecord{}
	for i := range bRaw {
		// Stamp B-file load order so the per-chrom and per-A match sorts can
		// preserve input order on equal (chrom, start) keys.
		bRaw[i].order = i
		rr := bRaw[i]
		bByChrom[rr.rec.Chrom] = append(bByChrom[rr.rec.Chrom], rr)
	}
	for chrom := range bByChrom {
		recs := bByChrom[chrom]
		// Sort by chromStart only, preserving input order on ties. Upstream
		// `bedtools map` consumes a (chrom, start)-sorted stream and never uses
		// chromEnd as a tie-break, so equal-start records keep input order.
		sort.SliceStable(recs, func(i, j int) bool {
			return recs[i].rec.ChromStart < recs[j].rec.ChromStart
		})
		bByChrom[chrom] = recs
	}
	// Build trees + an index from *bed.Record back to its rawRecord so we can
	// recover the fields after a tree query.
	trees := map[string]*bed.IntervalTree{}
	rawIndex := map[*bed.Record]rawRecord{}
	for chrom, recs := range bByChrom {
		recPtrs := make([]*bed.Record, len(recs))
		for i := range recs {
			recPtrs[i] = recs[i].rec
			rawIndex[recs[i].rec] = recs[i]
		}
		trees[chrom] = bed.NewIntervalTree(recPtrs)
	}

	// Stream A line-by-line: we want to preserve A's original text columns
	// verbatim (`bedtools map` echoes A's full record then appends).
	bw := bufio.NewWriter(writer)
	defer bw.Flush()
	scanner := bufio.NewScanner(readerA)
	scanner.Buffer(make([]byte, 0, 64*1024), 4*1024*1024)
	count := 0
	for scanner.Scan() {
		line := strings.TrimRight(scanner.Text(), "\r")
		trimmed := strings.TrimSpace(line)
		if trimmed == "" {
			continue
		}
		if strings.HasPrefix(trimmed, "#") || strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			// Under -header, echo A's leading header lines verbatim.
			if opts.Header {
				if _, err := fmt.Fprintln(bw, line); err != nil {
					return count, err
				}
			}
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 3 {
			return count, fmt.Errorf("A record must have at least 3 columns, got %d", len(fields))
		}
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return count, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return count, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
		}
		recA := &bed.Record{Chrom: fields[0], ChromStart: start, ChromEnd: end}
		if len(fields) > 5 {
			recA.Strand = fields[5]
		}
		var aBlocks []block
		if opts.Split {
			aBlocks = bed12Blocks(fields, start)
		}

		// Query B.
		var matches []rawRecord
		if tree, ok := trees[recA.Chrom]; ok {
			candidates := tree.Query(recA)
			for _, c := range candidates {
				if !strandPass(recA, c, opts) {
					continue
				}
				rr, ok := rawIndex[c]
				if !ok {
					continue
				}
				if opts.Split && !splitOverlap(aBlocks, rr) {
					continue
				}
				if !fractionPass(recA, c, opts) {
					continue
				}
				matches = append(matches, rr)
			}
			// Restore upstream's stream order; the tree query may return
			// candidates out of order. The key is (chromStart, B-file order) —
			// NOT chromEnd — so equal-start records keep input order, matching
			// upstream's collapse/distinct value order.
			sort.SliceStable(matches, func(i, j int) bool {
				if matches[i].rec.ChromStart != matches[j].rec.ChromStart {
					return matches[i].rec.ChromStart < matches[j].rec.ChromStart
				}
				return matches[i].order < matches[j].order
			})
		}

		// Compute aggregated columns.
		extras := make([]string, len(opts.Columns))
		var nonNumVal string // last non-numeric value seen, for the per-row warning
		var nonNumCol int
		nonNum := false
		for i, col := range opts.Columns {
			op := opts.Ops[i]
			vals := make([]string, len(matches))
			for j := range matches {
				vals[j] = matches[j].fields[col-1]
			}
			if numericOps[op] {
				res, lastBad, badCol, sawBad := applyNumericOp(op, col, vals, opts)
				if sawBad {
					nonNum = true
					nonNumVal = lastBad
					nonNumCol = badCol
				}
				extras[i] = res
				continue
			}
			if len(matches) == 0 {
				// count / count_distinct produce 0; everything else is null.
				if op == "count" || op == "count_distinct" {
					extras[i] = "0"
				} else {
					extras[i] = opts.Null
				}
				continue
			}
			res, err := bedmerge.ApplyOp(op, col, vals)
			if err != nil {
				return count, err
			}
			if (op == "collapse" || op == "distinct") && opts.Delim != "," {
				res = strings.ReplaceAll(res, ",", opts.Delim)
			}
			extras[i] = res
		}

		// Write A's original columns + extras.
		out := append([]string(nil), fields...)
		out = append(out, extras...)
		if _, err := fmt.Fprintln(bw, strings.Join(out, "\t")); err != nil {
			return count, err
		}
		// After all columns of this row, emit the single non-numeric WARNING
		// upstream prints when a numeric op encountered a non-numeric value.
		if nonNum && opts.WarnWriter != nil {
			fmt.Fprintf(opts.WarnWriter, " ***** WARNING: Non numeric value %s in %d.\n", nonNumVal, nonNumCol)
		}
		count++
	}
	if err := scanner.Err(); err != nil {
		return count, fmt.Errorf("error reading A: %w", err)
	}
	return count, nil
}

// applyNumericOp computes a numeric aggregation over vals, mirroring upstream
// KeyListOps: non-numeric values become NaN and are still included in the
// running computation (so e.g. sum over any non-numeric value is NaN), the
// last non-numeric value + column are returned so the caller can warn, and a
// NaN result (including the empty-group case) is rendered as the null value.
func applyNumericOp(op string, col int, vals []string, opts Options) (result, lastBad string, badCol int, sawBad bool) {
	nums := make([]float64, len(vals))
	for i, v := range vals {
		f, err := strconv.ParseFloat(strings.TrimSpace(v), 64)
		if err != nil {
			nums[i] = math.NaN()
			lastBad = v
			badCol = col
			sawBad = true
			continue
		}
		nums[i] = f
	}
	val := computeNumeric(op, nums)
	if math.IsNaN(val) {
		return opts.Null, lastBad, badCol, sawBad
	}
	return formatPrec(val, opts.Precision), lastBad, badCol, sawBad
}

// computeNumeric returns the numeric op's value, or NaN for the empty group
// (upstream returns NaN when the value list is empty). NaN values in nums
// propagate through sum/mean/stdev (any NaN ⇒ NaN), matching upstream's
// behaviour of summing atof()'d-to-NaN non-numeric values.
func computeNumeric(op string, nums []float64) float64 {
	if len(nums) == 0 {
		return math.NaN()
	}
	switch op {
	case "sum":
		s := 0.0
		for _, n := range nums {
			s += n
		}
		return s
	case "mean":
		s := 0.0
		for _, n := range nums {
			s += n
		}
		return s / float64(len(nums))
	case "min":
		m := nums[0]
		for _, n := range nums[1:] {
			if n < m {
				m = n
			}
		}
		return m
	case "max":
		m := nums[0]
		for _, n := range nums[1:] {
			if n > m {
				m = n
			}
		}
		return m
	case "absmin":
		m := math.Abs(nums[0])
		for _, n := range nums[1:] {
			if a := math.Abs(n); a < m {
				m = a
			}
		}
		return m
	case "absmax":
		m := math.Abs(nums[0])
		for _, n := range nums[1:] {
			if a := math.Abs(n); a > m {
				m = a
			}
		}
		return m
	case "median":
		sorted := append([]float64(nil), nums...)
		sort.Float64s(sorted)
		n := len(sorted)
		if n%2 == 1 {
			return sorted[n/2]
		}
		return (sorted[n/2-1] + sorted[n/2]) / 2
	case "stdev":
		mean := 0.0
		for _, n := range nums {
			mean += n
		}
		mean /= float64(len(nums))
		sq := 0.0
		for _, n := range nums {
			d := n - mean
			sq += d * d
		}
		return math.Sqrt(sq / float64(len(nums)))
	case "sstdev":
		if len(nums) == 1 {
			return math.NaN()
		}
		mean := 0.0
		for _, n := range nums {
			mean += n
		}
		mean /= float64(len(nums))
		sq := 0.0
		for _, n := range nums {
			d := n - mean
			sq += d * d
		}
		return math.Sqrt(sq / float64(len(nums)-1))
	}
	return math.NaN()
}

// formatPrec formats v with prec significant digits, matching upstream's C++
// `std::setprecision(prec) << val` (the default float format, equivalent to
// Go's 'g' verb): integer-valued results print with no decimal point and other
// values carry up to prec significant digits with no trailing-zero noise.
func formatPrec(v float64, prec int) string {
	return strconv.FormatFloat(v, 'g', prec, 64)
}

// splitOverlap reports whether any of A's blocks overlaps any of B's blocks,
// the overlap rule upstream `bedtools map -split` uses. B's block list is
// derived on demand from its BED12 columns (or its whole span when not BED12).
func splitOverlap(aBlocks []block, b rawRecord) bool {
	bBlocks := bed12Blocks(b.fields, b.rec.ChromStart)
	for _, ab := range aBlocks {
		for _, bb := range bBlocks {
			s := ab.start
			if bb.start > s {
				s = bb.start
			}
			e := ab.end
			if bb.end < e {
				e = bb.end
			}
			if e > s {
				return true
			}
		}
	}
	return false
}

// bed12Blocks derives the sub-block intervals (0-based, absolute) of a BED12
// line whose chromStart is start. Columns 10/11/12 hold blockCount, the
// comma-separated blockSizes, and the comma-separated blockStarts (relative to
// chromStart). For non-BED12 lines it falls back to the whole [start,end) span.
func bed12Blocks(fields []string, start int) []block {
	if len(fields) < 12 {
		return []block{{start: start, end: spanEnd(fields, start)}}
	}
	sizes := strings.Split(strings.TrimRight(fields[10], ","), ",")
	starts := strings.Split(strings.TrimRight(fields[11], ","), ",")
	n := len(sizes)
	if len(starts) < n {
		n = len(starts)
	}
	var out []block
	for i := 0; i < n; i++ {
		sz, e1 := strconv.Atoi(strings.TrimSpace(sizes[i]))
		bs, e2 := strconv.Atoi(strings.TrimSpace(starts[i]))
		if e1 != nil || e2 != nil {
			continue
		}
		out = append(out, block{start: start + bs, end: start + bs + sz})
	}
	if len(out) == 0 {
		out = append(out, block{start: start, end: spanEnd(fields, start)})
	}
	return out
}

// spanEnd returns the chromEnd from column 3, or start when it cannot be parsed.
func spanEnd(fields []string, start int) int {
	if len(fields) >= 3 {
		if e, err := strconv.Atoi(fields[2]); err == nil {
			return e
		}
	}
	return start
}

// parseIntervalLine parses one B record, auto-detecting BED vs GFF (and VCF).
// Retained for unit tests; the streaming reader uses readBRecords (input.go)
// which also handles BAM. A BED record has a numeric chromStart in column 2; a
// GFF record (8 or 9 columns) has a textual source in column 2 and numeric
// start/end in columns 4/5.
func parseIntervalLine(fields []string) (rawRecord, error) {
	f, ok := detectTextFormat(fields)
	if !ok {
		second := ""
		if len(fields) > 1 {
			second = fields[1]
		}
		return rawRecord{}, fmt.Errorf("invalid record: column 2 %q is not a numeric BED start and the line is not a GFF feature", second)
	}
	return parseTextRecord(fields, f)
}

// strandPass enforces -s / -S.
func strandPass(a, b *bed.Record, opts Options) bool {
	if opts.SameStrand {
		if a.Strand == "" || b.Strand == "" {
			return false
		}
		if a.Strand != b.Strand {
			return false
		}
	}
	if opts.OppositeStrand {
		if a.Strand == "" || b.Strand == "" {
			return false
		}
		if a.Strand == b.Strand {
			return false
		}
	}
	return true
}

// fractionPass enforces -f / -F (and -r). Returns true if B should be
// considered against A.
func fractionPass(a, b *bed.Record, opts Options) bool {
	if opts.FractionA == 0 && opts.FractionB == 0 {
		return true
	}
	overlapStart := a.ChromStart
	if b.ChromStart > overlapStart {
		overlapStart = b.ChromStart
	}
	overlapEnd := a.ChromEnd
	if b.ChromEnd < overlapEnd {
		overlapEnd = b.ChromEnd
	}
	ov := overlapEnd - overlapStart
	if ov <= 0 {
		return false
	}
	lenA := a.ChromEnd - a.ChromStart
	lenB := b.ChromEnd - b.ChromStart
	passA := opts.FractionA == 0 || (lenA > 0 && float64(ov)/float64(lenA) >= opts.FractionA)
	passB := opts.FractionB == 0 || (lenB > 0 && float64(ov)/float64(lenB) >= opts.FractionB)
	return passA && passB
}
