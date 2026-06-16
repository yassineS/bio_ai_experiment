package bedmerge

import (
	"fmt"
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
	// Delim is the delimiter joining the values produced by the list-style
	// operations (collapse, distinct, distinct_only, distinct_sort_num[_desc],
	// freqasc, freqdesc) — upstream `bedtools merge -delim`. An empty value is
	// treated as the default ",". The concat/cat family always joins with no
	// delimiter regardless of this setting, matching upstream getConcat.
	Delim string
}

// validOps is the set of supported aggregation operations. Matches the
// vocabulary of bedtools' `KeyListOps` (see
// reference_code/bedtools/src/utils/KeyListOps/KeyListOps.cpp).
var validOps = map[string]bool{
	"sum":                    true,
	"min":                    true,
	"max":                    true,
	"absmin":                 true,
	"absmax":                 true,
	"mean":                   true,
	"median":                 true,
	"stdev":                  true,
	"sstdev":                 true,
	"count":                  true,
	"count_distinct":         true,
	"distinct":               true,
	"collapse":               true,
	"cat":                    true,
	"cat_uniq":               true,
	"first":                  true,
	"last":                   true,
	"mode":                   true,
	"antimode":               true,
	"concat":                 true,
	"distinct_only":          true,
	"distinct_sort_num":      true,
	"distinct_sort_num_desc": true,
	"freqasc":                true,
	"freqdesc":               true,
}

// numericOps is the set of operations that consume their column values as
// numbers. A non-numeric value contributes NaN (and triggers a warning) rather
// than aborting, matching upstream KeyListOpsMethods::getColValNum.
var numericOps = map[string]bool{
	"sum":                    true,
	"min":                    true,
	"max":                    true,
	"absmin":                 true,
	"absmax":                 true,
	"mean":                   true,
	"median":                 true,
	"stdev":                  true,
	"sstdev":                 true,
	"mode":                   true,
	"antimode":               true,
	"distinct_sort_num":      true,
	"distinct_sort_num_desc": true,
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

// delim returns the effective list-join delimiter, defaulting to "," when
// unset.
func (co *ColumnOps) delim() string {
	if co == nil || co.Delim == "" {
		return ","
	}
	return co.Delim
}

// applyColumnOps aggregates each requested column over a merged group of
// records, returning one output string per column op. It mirrors upstream's
// KeyListOps::getOpVals: a single warning string is returned holding the LAST
// "Non numeric value" message seen across all columns/ops in this group (the
// upstream _errMsg is overwritten per value and printed once per record).
func applyColumnOps(co *ColumnOps, group []record, precision int) ([]string, string) {
	out := make([]string, len(co.Columns))
	warn := ""
	for i, col := range co.Columns {
		op := co.Ops[i]
		vals := make([]string, len(group))
		for j, r := range group {
			vals[j] = bamSafeField(r, col)
		}
		res, w := applyOpDelim(op, col, vals, co.delim(), precision)
		if w != "" {
			warn = w
		}
		out[i] = res
	}
	return out, warn
}

// bamSafeField returns the value of 1-based column col from a record's fields.
// For BAM records an empty or absent field becomes the null value "." (upstream
// KeyListOpsMethods::getColVal returns _nullVal for an empty BAM field), so
// list ops like collapse render missing mate info as ".". For text records the
// raw value is returned (a column past the width yields "").
func bamSafeField(r record, col int) string {
	v := ""
	if col-1 < len(r.fields) {
		v = r.fields[col-1]
	}
	if r.isBAM && v == "" {
		return "."
	}
	return v
}

// ApplyOp applies a column op (see validOps) to a slice of string values
// originating from column col (1-based, used only for error messages).
// Exported so other tools (bedgroupby, bedmap, bedcoverage) can reuse the same
// op vocabulary as bedmerge. It uses the "," list delimiter and the default
// output precision, and returns an error when a numeric op meets a non-numeric
// value (the standalone-helper contract, distinct from the merge path which
// warns and emits the null value).
func ApplyOp(op string, col int, vals []string) (string, error) {
	res, warn := applyOpDelim(op, col, vals, ",", DefaultPrecision)
	if warn != "" {
		return "", fmt.Errorf("operation %q requires numeric values, but column %d contains a non-numeric value", op, col)
	}
	return res, nil
}

// applyOpDelim applies op to vals, joining list results with delim and
// formatting floats to precision significant digits. For numeric ops a
// non-numeric value contributes NaN and sets the returned warning to the
// upstream-formatted message for that value; a NaN aggregate result prints as
// the null value ".".
func applyOpDelim(op string, col int, vals []string, delim string, precision int) (result, warn string) {
	if numericOps[op] {
		nums := make([]float64, len(vals))
		for i, v := range vals {
			if isNumericUpstream(v) {
				f, _ := strconv.ParseFloat(strings.TrimSpace(v), 64)
				nums[i] = f
			} else {
				nums[i] = math.NaN()
				warn = fmt.Sprintf(" ***** WARNING: Non numeric value %s in %d.", v, col)
			}
		}
		return numericResult(op, nums, delim, precision), warn
	}

	switch op {
	case "count":
		return strconv.Itoa(len(vals)), ""
	case "count_distinct":
		seen := map[string]bool{}
		n := 0
		for _, v := range vals {
			if !seen[v] {
				seen[v] = true
				n++
			}
		}
		return strconv.Itoa(n), ""
	case "distinct":
		// Upstream getDistinct iterates its std::map<string,int> freqMap, so the
		// unique values come out in ascending value-string order.
		return strings.Join(sortedKeys(freqCounts(vals)), delim), ""
	case "collapse":
		return strings.Join(vals, delim), ""
	case "cat", "concat":
		// Concatenate all values with no separator (upstream getConcat).
		return strings.Join(vals, ""), ""
	case "cat_uniq":
		seen := map[string]bool{}
		var out []string
		for _, v := range vals {
			if !seen[v] {
				seen[v] = true
				out = append(out, v)
			}
		}
		return strings.Join(out, ""), ""
	case "distinct_only":
		counts := freqCounts(vals)
		var out []string
		for _, k := range sortedKeys(counts) {
			if counts[k] == 1 {
				out = append(out, k)
			}
		}
		return strings.Join(out, delim), ""
	case "freqasc", "freqdesc":
		counts := freqCounts(vals)
		keys := sortedKeys(counts)
		sort.SliceStable(keys, func(i, j int) bool {
			if op == "freqasc" {
				return counts[keys[i]] < counts[keys[j]]
			}
			return counts[keys[i]] > counts[keys[j]]
		})
		var out []string
		for _, k := range keys {
			out = append(out, fmt.Sprintf("%s:%d", k, counts[k]))
		}
		return strings.Join(out, delim), ""
	case "first":
		if len(vals) == 0 {
			return ".", ""
		}
		return vals[0], ""
	case "last":
		if len(vals) == 0 {
			return ".", ""
		}
		return vals[len(vals)-1], ""
	}
	return "", ""
}

// numericResult computes the value of a numeric op over nums, returning the
// null value "." when the aggregate is NaN (e.g. a non-numeric value was
// present, or sstdev with a single element).
func numericResult(op string, nums []float64, delim string, precision int) string {
	switch op {
	case "sum":
		s := 0.0
		for _, n := range nums {
			s += n
		}
		return formatNum(s, precision)
	case "min":
		m := nums[0]
		for _, n := range nums[1:] {
			if n < m {
				m = n
			}
		}
		return formatNum(m, precision)
	case "max":
		m := nums[0]
		for _, n := range nums[1:] {
			if n > m {
				m = n
			}
		}
		return formatNum(m, precision)
	case "absmin":
		m := math.Abs(nums[0])
		for _, n := range nums[1:] {
			if a := math.Abs(n); a < m {
				m = a
			}
		}
		return formatNum(m, precision)
	case "absmax":
		m := math.Abs(nums[0])
		for _, n := range nums[1:] {
			if a := math.Abs(n); a > m {
				m = a
			}
		}
		return formatNum(m, precision)
	case "mean":
		s := 0.0
		for _, n := range nums {
			s += n
		}
		return formatNum(s/float64(len(nums)), precision)
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
		return formatNum(med, precision)
	case "stdev":
		// Population standard deviation: sqrt(Σ(x-μ)² / n).
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
		return formatNum(math.Sqrt(sq/float64(len(nums))), precision)
	case "sstdev":
		// Sample standard deviation: sqrt(Σ(x-μ)² / (n-1)); NaN -> "." when n==1.
		if len(nums) == 1 {
			return "."
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
		return formatNum(math.Sqrt(sq/float64(len(nums)-1)), precision)
	case "distinct_sort_num", "distinct_sort_num_desc":
		sorted := make([]float64, len(nums))
		copy(sorted, nums)
		sort.Float64s(sorted)
		if op == "distinct_sort_num_desc" {
			for i, j := 0, len(sorted)-1; i < j; i, j = i+1, j-1 {
				sorted[i], sorted[j] = sorted[j], sorted[i]
			}
		}
		var out []string
		for i, n := range sorted {
			if i > 0 && n == sorted[i-1] {
				continue
			}
			out = append(out, formatNum(n, precision))
		}
		return strings.Join(out, delim)
	case "mode":
		return modeOrAntimodeNum(nums, true, precision)
	case "antimode":
		return modeOrAntimodeNum(nums, false, precision)
	}
	return "."
}

// freqCounts tallies how many times each value appears.
func freqCounts(vals []string) map[string]int {
	counts := make(map[string]int, len(vals))
	for _, v := range vals {
		counts[v]++
	}
	return counts
}

// sortedKeys returns the keys of counts in ascending string order, matching the
// iteration order of upstream bedtools' std::map<string,int> freqMap.
func sortedKeys(counts map[string]int) []string {
	keys := make([]string, 0, len(counts))
	for k := range counts {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

// modeOrAntimodeNum returns the most (mode) or least (antimode) frequent numeric
// value, formatting it to the requested precision. Ties are broken by first-seen
// order, matching upstream.
func modeOrAntimodeNum(nums []float64, mode bool, precision int) string {
	counts := map[float64]int{}
	var order []float64
	for _, v := range nums {
		if _, ok := counts[v]; !ok {
			order = append(order, v)
		}
		counts[v]++
	}
	if len(order) == 0 {
		return "."
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
	return formatNum(best, precision)
}

// isNumericUpstream replicates upstream ParseTools::isNumeric: a string is
// numeric if every character is a digit, '+', '-', '.', 'e', or 'E', and it
// contains at least one digit.
func isNumericUpstream(s string) bool {
	hasDigit := false
	for i := 0; i < len(s); i++ {
		c := s[i]
		switch {
		case c >= '0' && c <= '9':
			hasDigit = true
		case c == '+' || c == '-' || c == '.' || c == 'e' || c == 'E':
		default:
			return false
		}
	}
	return hasDigit
}

// formatNum formats a float to `precision` significant digits using Go's 'g'
// verb, which matches C++ std::ostream << setprecision(precision) << val (the
// default float format upstream uses). A NaN result prints as the null value
// ".".
func formatNum(v float64, precision int) string {
	if math.IsNaN(v) {
		return "."
	}
	if precision <= 0 {
		precision = DefaultPrecision
	}
	return strconv.FormatFloat(v, 'g', precision, 64)
}
