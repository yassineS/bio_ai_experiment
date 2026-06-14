// Package bedmap applies column aggregation ops to A's intervals using B as
// the value source, mirroring upstream `bedtools map`.
//
// For each interval in A the package finds the overlapping B records, takes
// the values from one or more B columns, and runs the requested aggregation
// (e.g. sum, mean, count, distinct, collapse, ...). The aggregated values are
// appended to A's columns. When no B records overlap, the configured Null
// string is emitted instead (default ".").
//
// Column-op semantics come from bedmerge.ApplyOp, so the same vocabulary
// (sum, min, max, mean, median, count, count_distinct, distinct, collapse,
// first, last, mode, antimode) is supported.
package bedmap

import (
	"bufio"
	"fmt"
	"io"
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

	// Null is the placeholder emitted when A has no overlapping B record.
	// Default ".".
	Null string
	// Delim is the separator used by collapse / distinct (and other
	// "concatenate" ops). Default ",".
	Delim string

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
}

// Validate fills in defaults and rejects obviously bad combinations. It
// mutates opts in place.
func (opts *Options) Validate() error {
	if len(opts.Columns) == 0 {
		opts.Columns = []int{5}
	}
	for _, c := range opts.Columns {
		if c < 1 {
			return fmt.Errorf("column numbers must be >= 1, got %d", c)
		}
	}
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
		return fmt.Errorf("number of ops (%d) must be 1 or equal to number of columns (%d)", len(opts.Ops), len(opts.Columns))
	}
	if opts.Null == "" {
		opts.Null = "."
	}
	if opts.Delim == "" {
		opts.Delim = ","
	}
	if opts.SameStrand && opts.OppositeStrand {
		return fmt.Errorf("-s and -S are mutually exclusive")
	}
	return nil
}

// rawRecord pairs a parsed BED record with its raw whitespace-split columns,
// so column-value extraction by 1-based index works on the original text.
type rawRecord struct {
	rec    *bed.Record
	fields []string
}

// Map runs the per-A column aggregation against B and writes the results to
// writer. Returns the number of A records processed.
func Map(readerA, readerB io.Reader, writer io.Writer, opts Options) (int, error) {
	if err := opts.Validate(); err != nil {
		return 0, err
	}

	// Read B into per-chromosome sorted slices, then build per-chromosome
	// interval trees. We need both the parsed record (for the tree) and the
	// raw fields (for arbitrary column extraction), so we read B as plain
	// text and parse a Record per line ourselves.
	bRaw, err := readRawRecords(readerB)
	if err != nil {
		return 0, fmt.Errorf("error reading B: %w", err)
	}
	bByChrom := map[string][]rawRecord{}
	for _, rr := range bRaw {
		// Validate the requested columns exist on every B record.
		for _, c := range opts.Columns {
			if c > len(rr.fields) {
				return 0, fmt.Errorf("B record at %s:%d-%d has %d columns but column %d requested",
					rr.rec.Chrom, rr.rec.ChromStart, rr.rec.ChromEnd, len(rr.fields), c)
			}
		}
		bByChrom[rr.rec.Chrom] = append(bByChrom[rr.rec.Chrom], rr)
	}
	for chrom := range bByChrom {
		recs := bByChrom[chrom]
		sort.SliceStable(recs, func(i, j int) bool {
			if recs[i].rec.ChromStart != recs[j].rec.ChromStart {
				return recs[i].rec.ChromStart < recs[j].rec.ChromStart
			}
			return recs[i].rec.ChromEnd < recs[j].rec.ChromEnd
		})
		bByChrom[chrom] = recs
	}
	// Build trees + an index from *bed.Record back to its rawRecord so we can
	// recover the fields after a tree query.
	trees := map[string]*bed.IntervalTree{}
	rawIndex := map[*bed.Record]rawRecord{}
	for chrom, recs := range bByChrom {
		recPtrs := make([]*bed.Record, len(recs))
		for i, rr := range recs {
			recPtrs[i] = rr.rec
			rawIndex[rr.rec] = rr
		}
		trees[chrom] = bed.NewIntervalTree(recPtrs)
	}

	// Stream A line-by-line: we want to preserve A's original text columns
	// verbatim too (`bedtools map` echoes A's full record then appends).
	bw := bufio.NewWriter(writer)
	defer bw.Flush()
	scanner := bufio.NewScanner(readerA)
	scanner.Buffer(make([]byte, 0, 64*1024), 4*1024*1024)
	count := 0
	for scanner.Scan() {
		line := scanner.Text()
		trimmed := strings.TrimSpace(line)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") || strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
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

		// Query B.
		var matches []rawRecord
		if tree, ok := trees[recA.Chrom]; ok {
			candidates := tree.Query(recA)
			for _, c := range candidates {
				if !strandPass(recA, c, opts) {
					continue
				}
				if !fractionPass(recA, c, opts) {
					continue
				}
				if rr, ok := rawIndex[c]; ok {
					matches = append(matches, rr)
				}
			}
		}

		// Compute aggregated columns.
		extras := make([]string, len(opts.Columns))
		for i, col := range opts.Columns {
			op := opts.Ops[i]
			if len(matches) == 0 {
				// `count` and `count_distinct` always produce a number
				// (zero) when there are no matches; everything else gets
				// the null placeholder. Matches upstream `bedtools map`.
				if op == "count" || op == "count_distinct" {
					extras[i] = "0"
				} else {
					extras[i] = opts.Null
				}
				continue
			}
			vals := make([]string, len(matches))
			for j, m := range matches {
				vals[j] = m.fields[col-1]
			}
			res, err := bedmerge.ApplyOp(op, col, vals)
			if err != nil {
				return count, err
			}
			// Apply the configured delimiter to collapse/distinct outputs
			// (ApplyOp uses "," by default). Swap only if user asked.
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
		count++
	}
	if err := scanner.Err(); err != nil {
		return count, fmt.Errorf("error reading A: %w", err)
	}
	return count, nil
}

// readRawRecords reads a BED-like stream as (parsed Record, raw fields)
// pairs. The parsed Record is needed for interval-tree queries; the raw
// fields are needed so any column can be extracted by index.
func readRawRecords(r io.Reader) ([]rawRecord, error) {
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 4*1024*1024)
	var out []rawRecord
	for scanner.Scan() {
		line := scanner.Text()
		trimmed := strings.TrimSpace(line)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") || strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 3 {
			return nil, fmt.Errorf("BED record must have at least 3 fields, got %d", len(fields))
		}
		rr, err := parseIntervalLine(fields)
		if err != nil {
			return nil, err
		}
		out = append(out, rr)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// parseIntervalLine parses one B record, auto-detecting BED vs GFF. A BED
// record has a numeric chromStart in column 2; a GFF record (9 columns, 1-based)
// has a textual source in column 2 and numeric start/end in columns 4/5. The
// original fields are preserved verbatim so the -c column extraction reads the
// literal source columns (for GFF these are the 1-based GFF columns, matching
// upstream bedtools map -b <gff>).
func parseIntervalLine(fields []string) (rawRecord, error) {
	// BED: columns 2 and 3 are the 0-based half-open coordinates.
	if start, err := strconv.Atoi(fields[1]); err == nil {
		end, err2 := strconv.Atoi(fields[2])
		if err2 != nil {
			return rawRecord{}, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err2)
		}
		rec := &bed.Record{Chrom: fields[0], ChromStart: start, ChromEnd: end}
		if len(fields) > 5 {
			rec.Strand = fields[5]
		}
		return rawRecord{rec: rec, fields: fields}, nil
	}
	// GFF: column 2 (source) is non-numeric; columns 4/5 are 1-based start/end.
	if len(fields) >= 8 {
		gstart, err := strconv.Atoi(fields[3])
		if err == nil {
			gend, err2 := strconv.Atoi(fields[4])
			if err2 == nil {
				rec := &bed.Record{Chrom: fields[0], ChromStart: gstart - 1, ChromEnd: gend}
				if len(fields) > 6 {
					rec.Strand = fields[6]
				}
				return rawRecord{rec: rec, fields: fields}, nil
			}
		}
	}
	return rawRecord{}, fmt.Errorf("invalid record: column 2 %q is not a numeric BED start and the line is not a GFF feature", fields[1])
}

// strandPass: same as bedcoverage; duplicated here to avoid a cross-tool
// import (the two ports stay independent at the package layer).
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
	if opts.Reciprocal {
		return passA && passB
	}
	return passA && passB
}
