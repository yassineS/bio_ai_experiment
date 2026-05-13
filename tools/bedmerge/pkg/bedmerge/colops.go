package bedmerge

import (
	"bufio"
	"fmt"
	"io"
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

// validOps is the set of supported aggregation operations.
var validOps = map[string]bool{
	"sum":            true,
	"min":            true,
	"max":            true,
	"mean":           true,
	"median":         true,
	"count":          true,
	"count_distinct": true,
	"distinct":       true,
	"collapse":       true,
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
	"mean":     true,
	"median":   true,
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

	sort.SliceStable(intervals, func(i, j int) bool {
		if intervals[i].chrom != intervals[j].chrom {
			return intervals[i].chrom < intervals[j].chrom
		}
		if intervals[i].start != intervals[j].start {
			return intervals[i].start < intervals[j].start
		}
		return intervals[i].end < intervals[j].end
	})

	w := bufio.NewWriter(writer)
	outCount := 0

	flushGroup := func(group []colInterval) error {
		if len(group) == 0 {
			return nil
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
			res, err := applyOp(op, col, vals)
			if err != nil {
				return err
			}
			out = append(out, res)
		}
		if _, err := fmt.Fprintln(w, strings.Join(out, "\t")); err != nil {
			return err
		}
		outCount++
		return nil
	}

	var group []colInterval
	group = append(group, intervals[0])
	curChrom := intervals[0].chrom
	curStrand := intervals[0].strand
	curEnd := intervals[0].end

	for i := 1; i < len(intervals); i++ {
		iv := intervals[i]
		canMerge := false
		if iv.chrom == curChrom {
			if !opts.StrandSpec || iv.strand == curStrand || curStrand == "." || iv.strand == "." {
				if iv.start <= curEnd+opts.MaxDistance {
					canMerge = true
				}
			}
		}
		if canMerge {
			if iv.end > curEnd {
				curEnd = iv.end
			}
			group = append(group, iv)
		} else {
			if err := flushGroup(group); err != nil {
				return 0, err
			}
			group = []colInterval{iv}
			curChrom = iv.chrom
			curStrand = iv.strand
			curEnd = iv.end
		}
	}
	if err := flushGroup(group); err != nil {
		return 0, err
	}

	if err := w.Flush(); err != nil {
		return 0, fmt.Errorf("error flushing output: %w", err)
	}

	return outCount, nil
}

// applyOp applies a single aggregation operation to the slice of column values
// taken from a merged group (in input/sorted order). col is the 1-based column
// number, used only for error messages.
func applyOp(op string, col int, vals []string) (string, error) {
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
		return strings.Join(out, ","), nil
	case "collapse":
		return strings.Join(vals, ","), nil
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
