package bedmerge

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"
)

// ColumnOps describes a set of bedtools-merge-style column aggregations
// requested via the -c (columns) and -o (operations) flags.
type ColumnOps struct {
	// Columns holds the 1-based input column numbers to aggregate over.
	Columns []int
	// Ops holds the operation name for each column in Columns.
	Ops []string
}

// validOps is the set of supported aggregation operations. Matches the
// vocabulary of bedtools' `KeyListOps` (see
// reference_code/bedtools/src/utils/KeyListOps/KeyListOps.cpp).
var validOps = map[string]bool{
	"sum":            true,
	"min":            true,
	"max":            true,
	"absmin":         true,
	"absmax":         true,
	"mean":           true,
	"median":         true,
	"stdev":          true,
	"sstdev":         true,
	"count":          true,
	"count_distinct": true,
	"distinct":       true,
	"collapse":       true,
	"cat":            true,
	"cat_uniq":       true,
	"first":          true,
	"last":           true,
	"mode":           true,
	"antimode":       true,
}

// numericOps is the set of operations that require their column values to
// parse as numbers.
var numericOps = map[string]bool{
	"sum":      true,
	"min":      true,
	"max":      true,
	"absmin":   true,
	"absmax":   true,
	"mean":     true,
	"median":   true,
	"stdev":    true,
	"sstdev":   true,
	"mode":     true,
	"antimode": true,
}

// ParseColumnOps parses the comma-separated -c columns string and -o operations
// string into a ColumnOps. Either both must be empty (returns nil, nil) or both
// non-empty. If one op is given it is applied to every column; otherwise the
// number of ops must equal the number of columns.
func ParseColumnOps(colsStr, opsStr string) (*ColumnOps, error) {
	colsStr = strings.TrimSpace(colsStr)
	opsStr = strings.TrimSpace(opsStr)
	if colsStr == "" && opsStr == "" {
		return nil, nil
	}
	if colsStr == "" {
		return nil, fmt.Errorf("-o/--operations requires -c/--columns")
	}
	if opsStr == "" {
		return nil, fmt.Errorf("-c/--columns requires -o/--operations")
	}

	colParts := strings.Split(colsStr, ",")
	cols := make([]int, 0, len(colParts))
	for _, p := range colParts {
		p = strings.TrimSpace(p)
		n, err := strconv.Atoi(p)
		if err != nil {
			return nil, fmt.Errorf("invalid column number %q in -c/--columns", p)
		}
		if n < 1 {
			return nil, fmt.Errorf("column numbers in -c/--columns must be >= 1, got %d", n)
		}
		cols = append(cols, n)
	}

	opParts := strings.Split(opsStr, ",")
	for i, op := range opParts {
		opParts[i] = strings.TrimSpace(op)
		if !validOps[opParts[i]] {
			return nil, fmt.Errorf("unsupported operation %q in -o/--operations", opParts[i])
		}
	}

	var ops []string
	switch {
	case len(opParts) == len(cols):
		ops = opParts
	case len(opParts) == 1:
		ops = make([]string, len(cols))
		for i := range ops {
			ops[i] = opParts[0]
		}
	case len(cols) == 1:
		// Single column, many ops: apply each op to the same column, producing
		// one output column per op. This matches upstream bedtools merge.
		ops = opParts
		expanded := make([]int, len(opParts))
		for i := range expanded {
			expanded[i] = cols[0]
		}
		cols = expanded
	default:
		return nil, fmt.Errorf("number of operations (%d) must be 1 or equal to number of columns (%d)", len(opParts), len(cols))
	}

	return &ColumnOps{Columns: cols, Ops: ops}, nil
}

// colInterval is an interval kept with all of its raw input columns so that
// arbitrary columns can be aggregated after merging.
type colInterval struct {
	chrom  string
	start  int
	end    int
	strand string
	fields []string // all raw columns of the original line
}

// mergeWithColumnOps performs a merge while aggregating the requested input
// columns over each merged group, mirroring `bedtools merge -c ... -o ...`.
func mergeWithColumnOps(reader io.Reader, writer io.Writer, opts MergeOptions) (int, error) {
	co := opts.ColumnOps

	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	var intervals []colInterval
	for scanner.Scan() {
		line := strings.TrimRight(scanner.Text(), "\r\n")
		trimmed := strings.TrimSpace(line)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") || strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 3 {
			return 0, fmt.Errorf("BED record must have at least 3 fields, got %d", len(fields))
		}
		start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
		if err != nil {
			return 0, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
		}
		end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
		if err != nil {
			return 0, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
		}
		strand := ""
		if len(fields) > 5 {
			strand = fields[5]
		}
		// -S strand filter: drop records whose strand doesn't match.
		if opts.StrandFilter != "" && strand != opts.StrandFilter {
			continue
		}
		// Verify the requested columns exist.
		for _, c := range co.Columns {
			if c > len(fields) {
				return 0, fmt.Errorf("requested column %d but input line has only %d columns", c, len(fields))
			}
		}
		intervals = append(intervals, colInterval{
			chrom:  fields[0],
			start:  start,
			end:    end,
			strand: strand,
			fields: fields,
		})
	}
	if err := scanner.Err(); err != nil {
		return 0, fmt.Errorf("error reading BED input: %w", err)
	}

	if len(intervals) == 0 {
		return 0, nil
	}

	sortColIntervals(intervals)

	w := bufio.NewWriter(writer)

	// Under -s, upstream `bedtools merge -s` drops UNKNOWN-strand ("."
	// or "") records entirely and merges `+` / `-` independently, then
	// outputs both streams in (chrom, start, end) order. See
	// reference_code/bedtools/src/utils/FileRecordTools/FileRecordMergeMgr.cpp
	// lines 47-58 + 96-129. Replicate that here so that `merge.t15`
	// matches upstream byte-for-byte.
	if opts.StrandSpec {
		var plus, minus []colInterval
		for _, iv := range intervals {
			switch iv.strand {
			case "+":
				plus = append(plus, iv)
			case "-":
				minus = append(minus, iv)
				// "" and "." dropped.
			}
		}
		// Collect merged groups (don't write yet) so we can merge-sort
		// the two strand outputs by (chrom, start, end) before writing.
		pg := groupMergedColIntervals(plus, opts)
		mg := groupMergedColIntervals(minus, opts)
		all := mergeSortedColGroupsByPos(pg, mg)
		outCount := 0
		for _, g := range all {
			wrote, err := flushColumnGroup(w, co, g, opts.Delim)
			if err != nil {
				return 0, err
			}
			outCount += wrote
		}
		if err := w.Flush(); err != nil {
			return 0, fmt.Errorf("error flushing output: %w", err)
		}
		return outCount, nil
	}

	groups := groupMergedColIntervals(intervals, opts)
	outCount := 0
	for _, g := range groups {
		wrote, err := flushColumnGroup(w, co, g, opts.Delim)
		if err != nil {
			return 0, err
		}
		outCount += wrote
	}

	if err := w.Flush(); err != nil {
		return 0, fmt.Errorf("error flushing output: %w", err)
	}

	return outCount, nil
}

// sortColIntervals sorts colIntervals by (chrom, start, end) — used as a
// shared prerequisite for both the strand-aware and strand-agnostic
// merge paths.
func sortColIntervals(intervals []colInterval) {
	sort.SliceStable(intervals, func(i, j int) bool {
		if intervals[i].chrom != intervals[j].chrom {
			return intervals[i].chrom < intervals[j].chrom
		}
		if intervals[i].start != intervals[j].start {
			return intervals[i].start < intervals[j].start
		}
		return intervals[i].end < intervals[j].end
	})
}

// groupMergedColIntervals runs a single-pass position-only merge over
// sorted intervals (strand is ignored — the caller is responsible for
// any strand bucketing). The output is a list of groups, each of which
// is a slice of all the input intervals that merged together.
func groupMergedColIntervals(intervals []colInterval, opts MergeOptions) [][]colInterval {
	if len(intervals) == 0 {
		return nil
	}
	var out [][]colInterval
	group := []colInterval{intervals[0]}
	curChrom := intervals[0].chrom
	curEnd := intervals[0].end
	for i := 1; i < len(intervals); i++ {
		iv := intervals[i]
		if iv.chrom == curChrom && iv.start <= curEnd+opts.MaxDistance {
			if iv.end > curEnd {
				curEnd = iv.end
			}
			group = append(group, iv)
			continue
		}
		out = append(out, group)
		group = []colInterval{iv}
		curChrom = iv.chrom
		curEnd = iv.end
	}
	out = append(out, group)
	return out
}

// mergeSortedColGroupsByPos two-way-merges two slices of column-merged
// groups (each one already sorted by chrom+start) into a single stream,
// preserving (chrom, start, end) order.
func mergeSortedColGroupsByPos(a, b [][]colInterval) [][]colInterval {
	out := make([][]colInterval, 0, len(a)+len(b))
	i, j := 0, 0
	less := func(x, y []colInterval) bool {
		if x[0].chrom != y[0].chrom {
			return x[0].chrom < y[0].chrom
		}
		if x[0].start != y[0].start {
			return x[0].start < y[0].start
		}
		// Compute end of each group to break ties.
		xe := x[0].end
		for _, iv := range x[1:] {
			if iv.end > xe {
				xe = iv.end
			}
		}
		ye := y[0].end
		for _, iv := range y[1:] {
			if iv.end > ye {
				ye = iv.end
			}
		}
		return xe < ye
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

// flushColumnGroup writes one merged output line aggregating the requested
// columns over the given group of intervals. An empty group is a no-op and
// returns 0 written rows; otherwise it returns 1 on success.
func flushColumnGroup(w io.Writer, co *ColumnOps, group []colInterval, delim string) (int, error) {
	if len(group) == 0 {
		return 0, nil
	}
	chrom := group[0].chrom
	start := group[0].start
	end := group[0].end
	for _, iv := range group {
		if iv.end > end {
			end = iv.end
		}
	}
	out := []string{chrom, strconv.Itoa(start), strconv.Itoa(end)}
	for i, col := range co.Columns {
		op := co.Ops[i]
		vals := make([]string, len(group))
		for j, iv := range group {
			vals[j] = iv.fields[col-1]
		}
		res, err := applyOpDelim(op, col, vals, delim)
		if err != nil {
			return 0, err
		}
		out = append(out, res)
	}
	if _, err := fmt.Fprintln(w, strings.Join(out, "\t")); err != nil {
		return 0, err
	}
	return 1, nil
}

// applyOp applies a single aggregation operation to the slice of column values
// taken from a merged group (in input/sorted order). col is the 1-based column
// number, used only for error messages.
// ApplyOp applies a column op (see validOps) to a slice of string values
// originating from column col (1-based, used only for error messages).
// Exported so other tools (bedgroupby, bedmap, bedcoverage) can reuse the
// same op vocabulary as bedmerge. Uses the default "," join separator for
// collapse / distinct / freqdesc / freqasc; see ApplyOpDelim for an override.
func ApplyOp(op string, col int, vals []string) (string, error) {
	return applyOpDelim(op, col, vals, ",")
}

// ApplyOpDelim is the variant of ApplyOp that lets the caller override the
// join separator used by collapse / distinct / freqdesc / freqasc. Pass ""
// to fall back to ",". Mirrors `bedtools merge -delim`.
func ApplyOpDelim(op string, col int, vals []string, delim string) (string, error) {
	return applyOpDelim(op, col, vals, delim)
}

func applyOp(op string, col int, vals []string) (string, error) {
	return applyOpDelim(op, col, vals, ",")
}

func applyOpDelim(op string, col int, vals []string, delim string) (string, error) {
	if delim == "" {
		delim = ","
	}
	if numericOps[op] {
		nums := make([]float64, len(vals))
		for i, v := range vals {
			f, err := strconv.ParseFloat(strings.TrimSpace(v), 64)
			if err != nil {
				return "", fmt.Errorf("operation %q requires numeric values, but column %d contains non-numeric value %q", op, col, v)
			}
			nums[i] = f
		}
		switch op {
		case "sum":
			s := 0.0
			for _, n := range nums {
				s += n
			}
			return formatNum(s), nil
		case "min":
			m := nums[0]
			for _, n := range nums[1:] {
				if n < m {
					m = n
				}
			}
			return formatNum(m), nil
		case "max":
			m := nums[0]
			for _, n := range nums[1:] {
				if n > m {
					m = n
				}
			}
			return formatNum(m), nil
		case "absmin":
			// Upstream's `getAbsMin` (KeyListOpsMethods.cpp) returns the
			// minimum of `|x|` over the group — the sign of the original
			// value is intentionally dropped.
			m := math.Abs(nums[0])
			for _, n := range nums[1:] {
				if a := math.Abs(n); a < m {
					m = a
				}
			}
			return formatNum(m), nil
		case "absmax":
			// Upstream's `getAbsMax`: max of `|x|`; sign is dropped.
			m := math.Abs(nums[0])
			for _, n := range nums[1:] {
				if a := math.Abs(n); a > m {
					m = a
				}
			}
			return formatNum(m), nil
		case "mean":
			s := 0.0
			for _, n := range nums {
				s += n
			}
			return formatNum(s / float64(len(nums))), nil
		case "median":
			sorted := append([]float64(nil), nums...)
			sort.Float64s(sorted)
			n := len(sorted)
			var med float64
			if n%2 == 1 {
				med = sorted[n/2]
			} else {
				med = (sorted[n/2-1] + sorted[n/2]) / 2
			}
			return formatNum(med), nil
		case "stdev":
			// Population standard deviation: sqrt(Σ(x-μ)² / n).
			// Matches upstream `getStddev` in
			// reference_code/bedtools/src/utils/KeyListOps/KeyListOpsMethods.cpp.
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
			v := math.Sqrt(sq / float64(len(nums)))
			if math.IsNaN(v) {
				return ".", nil
			}
			return formatNum(v), nil
		case "sstdev":
			// Sample standard deviation: sqrt(Σ(x-μ)² / (n-1)).
			// Upstream returns NaN (printed as "." via getNullValue) when
			// n == 1; we replicate that.
			if len(nums) == 1 {
				return ".", nil
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
			v := math.Sqrt(sq / float64(len(nums)-1))
			if math.IsNaN(v) {
				return ".", nil
			}
			return formatNum(v), nil
		case "mode":
			return modeOrAntimode(vals, true), nil
		case "antimode":
			return modeOrAntimode(vals, false), nil
		}
	}

	switch op {
	case "count":
		return strconv.Itoa(len(vals)), nil
	case "count_distinct":
		seen := map[string]bool{}
		n := 0
		for _, v := range vals {
			if !seen[v] {
				seen[v] = true
				n++
			}
		}
		return strconv.Itoa(n), nil
	case "distinct":
		seen := map[string]bool{}
		var out []string
		for _, v := range vals {
			if !seen[v] {
				seen[v] = true
				out = append(out, v)
			}
		}
		return strings.Join(out, delim), nil
	case "collapse":
		return strings.Join(vals, delim), nil
	case "cat":
		// Concatenate all values with no separator. Mirrors upstream's
		// CONCAT op (`getConcat` in KeyListOpsMethods.cpp), which is the
		// `concat` operation in bedtools merge/groupby; the `cat` /
		// `cat_uniq` names are the friendlier aliases used by bedmap.
		return strings.Join(vals, ""), nil
	case "cat_uniq":
		// Concatenate unique values (first-appearance order) with no
		// separator.
		seen := map[string]bool{}
		var out []string
		for _, v := range vals {
			if !seen[v] {
				seen[v] = true
				out = append(out, v)
			}
		}
		return strings.Join(out, ""), nil
	case "first":
		return vals[0], nil
	case "last":
		return vals[len(vals)-1], nil
	}
	return "", fmt.Errorf("unsupported operation %q", op)
}

// modeOrAntimode returns the most (mode=true) or least (mode=false) frequent
// value; ties are broken by first-seen order.
func modeOrAntimode(vals []string, mode bool) string {
	counts := map[string]int{}
	var order []string
	for _, v := range vals {
		if _, ok := counts[v]; !ok {
			order = append(order, v)
		}
		counts[v]++
	}
	best := order[0]
	bestCount := counts[best]
	for _, v := range order[1:] {
		c := counts[v]
		if mode {
			if c > bestCount {
				best, bestCount = v, c
			}
		} else {
			if c < bestCount {
				best, bestCount = v, c
			}
		}
	}
	return best
}

// formatNum formats a float so that integer-valued results print without a
// decimal point and other values print with up to ~10 significant digits and no
// trailing-zero noise, matching bedtools' %g-ish output.
func formatNum(v float64) string {
	if v == float64(int64(v)) {
		return strconv.FormatInt(int64(v), 10)
	}
	return strconv.FormatFloat(v, 'g', -1, 64)
}
