// Package bedexpand implements `bedtools expand`: it takes a tab-delimited
// file (typically BED) and expands one or more columns whose values are
// comma-separated lists, producing one output row per element. When multiple
// columns are expanded together they are zipped in lock-step: the i-th
// element of each expanded column appears together on the i-th output row.
//
// Mirrors upstream's `bedtools expand -c COLS` (1-based comma-separated
// column list). Per the upstream algorithm (see
// reference_code/bedtools/src/expand/expand.cpp), output rows iterate
// through the row's original columns in their original positions; when a
// column is in -c, the value emitted there is taken from the n-th list
// element of the column referenced by the matching slot of -c (in -c
// order). So `-c 5,4` substitutes column-5 elements at position 4 and
// column-4 elements at position 5 — i.e. it swaps the two expanded
// columns. This matches expand.t3 in the upstream test suite.
package bedexpand

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Options configures Expand.
type Options struct {
	// Columns is the 1-based set of columns to expand, in the user-specified
	// order. The k-th column listed supplies the value emitted at the k-th
	// expanded position when walking the row left to right.
	Columns []int
}

// ParseColumns parses a comma-separated list of 1-based column indices, e.g.
// "4" or "5,4". Whitespace around individual numbers is tolerated.
func ParseColumns(s string) ([]int, error) {
	if strings.TrimSpace(s) == "" {
		return nil, fmt.Errorf("-c must list at least one column")
	}
	parts := strings.Split(s, ",")
	out := make([]int, 0, len(parts))
	for _, p := range parts {
		t := strings.TrimSpace(p)
		if t == "" {
			return nil, fmt.Errorf("empty column in -c %q", s)
		}
		n, err := strconv.Atoi(t)
		if err != nil {
			return nil, fmt.Errorf("invalid column %q in -c: %v", t, err)
		}
		if n <= 0 {
			return nil, fmt.Errorf("column must be 1-based (got %d)", n)
		}
		out = append(out, n)
	}
	return out, nil
}

// Expand reads tab-delimited records from r and writes the expanded rows to
// w. Returns the number of output records (data rows) written. Each input
// row must contain enough columns to cover every requested column index,
// and within a row every requested column must have the same number of
// comma-separated list elements (matching upstream).
//
// Empty lines are passed through unchanged; lines that start with '#',
// 'track', or 'browser' are emitted verbatim as headers.
func Expand(r io.Reader, w io.Writer, opts Options) (int, error) {
	if len(opts.Columns) == 0 {
		return 0, fmt.Errorf("Expand: at least one column required")
	}
	for _, c := range opts.Columns {
		if c <= 0 {
			return 0, fmt.Errorf("Expand: column must be 1-based (got %d)", c)
		}
	}

	// Build a quick membership set: position-in-row -> position-in-opts.Columns.
	// Walking the row left to right, the k-th expanded column we encounter
	// (k counted by occurrence) consumes opts.Columns[k] as its source list.
	expSet := make(map[int]bool, len(opts.Columns))
	for _, c := range opts.Columns {
		expSet[c] = true
	}

	bw := bufio.NewWriter(w)
	defer bw.Flush()

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	written := 0
	lineNo := 0
	for sc.Scan() {
		lineNo++
		raw := sc.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" {
			if _, err := bw.WriteString(raw); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return written, err
			}
			continue
		}
		if strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") ||
			strings.HasPrefix(trimmed, "browser") {
			if _, err := bw.WriteString(raw); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return written, err
			}
			continue
		}

		fields := strings.Split(raw, "\t")
		// Bounds check first: every requested column must exist.
		for _, c := range opts.Columns {
			if c > len(fields) {
				return written, fmt.Errorf("line %d: requested column %d but row has %d field(s)", lineNo, c, len(fields))
			}
		}

		// Split each requested column on commas (in opts.Columns order) and
		// verify the lists agree in length within the row.
		lists := make([][]string, len(opts.Columns))
		n := -1
		for i, c := range opts.Columns {
			lists[i] = strings.Split(fields[c-1], ",")
			if n == -1 {
				n = len(lists[i])
			} else if len(lists[i]) != n {
				return written, fmt.Errorf("line %d: each expanded column must have the same number of elements (%d vs %d)", lineNo, n, len(lists[i]))
			}
		}

		// Emit n rows. Upstream's algorithm: walk every column in the row's
		// natural order; for non-expanded columns print verbatim; for the
		// k-th expanded column encountered (k starts at 0), print lists[k][j].
		for j := 0; j < n; j++ {
			numExpSeen := 0
			row := make([]string, 0, len(fields))
			for ci, f := range fields {
				if expSet[ci+1] {
					row = append(row, lists[numExpSeen][j])
					numExpSeen++
				} else {
					row = append(row, f)
				}
			}
			if _, err := bw.WriteString(strings.Join(row, "\t")); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return written, err
			}
			written++
		}
	}
	if err := sc.Err(); err != nil {
		return written, err
	}
	return written, nil
}
