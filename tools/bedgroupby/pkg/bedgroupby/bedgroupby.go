// Package bedgroupby groups records from a tabular file (BED, BEDPLUS,
// VCF-like, or any TSV) by one or more grouping columns and applies aggregation
// operations to each group, mirroring upstream `bedtools groupby`.
package bedgroupby

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/tools/bedmerge/pkg/bedmerge"
)

// Options controls bedgroupby behaviour.
type Options struct {
	// GroupCols is the 1-based list of columns that define a group. Records
	// with identical values in these columns (consecutive in input order)
	// form one group. Defaults to {1, 2, 3} when nil/empty.
	GroupCols []int
	// AggCols is the 1-based list of columns to aggregate. Required.
	AggCols []int
	// Ops is the per-column operation name. Must be either len(AggCols) long
	// or a single op (applied to every AggCol).
	Ops []string
	// Full, when true, emits every original column of the first record in
	// each group before appending the aggregation outputs, matching
	// `bedtools groupby -full`.
	Full bool
	// IgnoreCase, when true, performs case-insensitive comparison of the
	// grouping-column values (the original case from the first record in
	// each group is preserved on output, matching upstream `-ignorecase`).
	IgnoreCase bool
	// InHeader, when true, drops the first non-empty data line as a header
	// even if it does not start with one of the recognised marker prefixes.
	InHeader bool
	// OutHeader, when true, prints either the input header (if one was
	// detected/forced) or a synthetic `col_1\tcol_2\t...` line.
	OutHeader bool
	// Header, when true, prints either the input header (if any) or a
	// synthetic positional header, matching `-header` upstream which is
	// stronger than `-outheader` for unmarked inputs.
	Header bool
}

// ParseGroupSpec parses a comma-separated -g column list (with ranges like
// "2-4") into a slice of 1-based column numbers. An empty string yields nil.
func ParseGroupSpec(spec string) ([]int, error) {
	spec = strings.TrimSpace(spec)
	if spec == "" {
		return nil, nil
	}
	var cols []int
	for _, part := range strings.Split(spec, ",") {
		part = strings.TrimSpace(part)
		if part == "" {
			continue
		}
		if dash := strings.IndexByte(part, '-'); dash > 0 {
			loStr := strings.TrimSpace(part[:dash])
			hiStr := strings.TrimSpace(part[dash+1:])
			lo, err := strconv.Atoi(loStr)
			if err != nil {
				return nil, fmt.Errorf("invalid column range %q in -g: %v", part, err)
			}
			hi, err := strconv.Atoi(hiStr)
			if err != nil {
				return nil, fmt.Errorf("invalid column range %q in -g: %v", part, err)
			}
			if lo < 1 || hi < lo {
				return nil, fmt.Errorf("invalid column range %q: must be >= 1 and increasing", part)
			}
			for c := lo; c <= hi; c++ {
				cols = append(cols, c)
			}
			continue
		}
		n, err := strconv.Atoi(part)
		if err != nil {
			return nil, fmt.Errorf("invalid column %q in -g: %v", part, err)
		}
		if n < 1 {
			return nil, fmt.Errorf("invalid column %q in -g: must be >= 1", part)
		}
		cols = append(cols, n)
	}
	return cols, nil
}

// Validate fills in defaults and sanity-checks the options. Mutates opts.
func (opts *Options) Validate() error {
	if len(opts.GroupCols) == 0 {
		opts.GroupCols = []int{1, 2, 3}
	}
	for _, c := range opts.GroupCols {
		if c < 1 {
			return fmt.Errorf("group column numbers must be >= 1, got %d", c)
		}
	}
	if len(opts.AggCols) == 0 {
		return fmt.Errorf("-c/--columns is required")
	}
	for _, c := range opts.AggCols {
		if c < 1 {
			return fmt.Errorf("aggregation column numbers must be >= 1, got %d", c)
		}
	}
	if len(opts.Ops) == 0 {
		// bedtools groupby's default op is `sum`.
		opts.Ops = []string{"sum"}
	}
	switch {
	case len(opts.Ops) == 1:
		// Broadcast the single op to every column.
		single := opts.Ops[0]
		opts.Ops = make([]string, len(opts.AggCols))
		for i := range opts.Ops {
			opts.Ops[i] = single
		}
	case len(opts.Ops) != len(opts.AggCols):
		return fmt.Errorf("there are %d columns given, but there are %d operations", len(opts.AggCols), len(opts.Ops))
	}
	return nil
}

// headerPrefixes are the leading-character prefixes that mark a line as a
// header/comment in bedtools' inputs.
var headerPrefixes = []string{"#", "track", "browser", "@"}

func isHeaderLine(line string) bool {
	for _, p := range headerPrefixes {
		if strings.HasPrefix(line, p) {
			return true
		}
	}
	return false
}

// Group reads records from reader, groups consecutive records sharing the
// configured grouping-column values, and writes one TSV row per group to
// writer. It returns the number of output rows.
func Group(reader io.Reader, writer io.Writer, opts Options) (int, error) {
	if err := opts.Validate(); err != nil {
		return 0, err
	}
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 64*1024), 8*1024*1024)
	w := bufio.NewWriter(writer)

	var (
		header         string
		headerEmitted  bool
		sawFirstRecord bool
		curKey         string
		curFirst       []string
		curGroupVals   []string   // per GroupCol: first-seen original-case value
		curVals        [][]string // per AggCol: values in input order
		outCount       int
	)

	flush := func() error {
		if !sawFirstRecord {
			return nil
		}
		if !headerEmitted && (opts.OutHeader || opts.Header) {
			switch {
			case header != "":
				if _, err := fmt.Fprintln(w, header); err != nil {
					return err
				}
			case opts.OutHeader:
				// Synthesise a `col_1\tcol_2\t...` header when there is no
				// marked/forced header. Matches upstream `groupby.t9`.
				cols := make([]string, len(curFirst))
				for i := range cols {
					cols[i] = "col_" + strconv.Itoa(i+1)
				}
				if _, err := fmt.Fprintln(w, strings.Join(cols, "\t")); err != nil {
					return err
				}
			}
			headerEmitted = true
		}

		var out []string
		if opts.Full {
			out = append(out, curFirst...)
		} else {
			// curGroupVals captures the first-seen value (with its
			// original case) for each grouping column. Under
			// -ignorecase that means the group's first row sets the
			// emitted spelling, which matches upstream.
			out = append(out, curGroupVals...)
		}
		for i, ac := range opts.AggCols {
			res, err := bedmerge.ApplyOp(opts.Ops[i], ac, curVals[i])
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

	firstDataLine := true
	for scanner.Scan() {
		raw := strings.TrimRight(scanner.Text(), "\r\n")
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" {
			continue
		}
		if isHeaderLine(raw) {
			if header == "" {
				header = raw
			}
			continue
		}
		// Treat the first unmarked data line as a header when -inheader or
		// -header is set.
		if firstDataLine {
			firstDataLine = false
			if (opts.InHeader || opts.Header) && header == "" {
				header = raw
				continue
			}
		}
		fields := strings.Split(raw, "\t")
		// Verify all referenced columns exist.
		maxCol := 0
		for _, c := range opts.GroupCols {
			if c > maxCol {
				maxCol = c
			}
		}
		for _, c := range opts.AggCols {
			if c > maxCol {
				maxCol = c
			}
		}
		if maxCol > len(fields) {
			return outCount, fmt.Errorf("requested column %d but input line has only %d columns", maxCol, len(fields))
		}

		key := groupKey(fields, opts.GroupCols, opts.IgnoreCase)
		if !sawFirstRecord {
			curKey = key
			curFirst = fields
			curGroupVals = pickGroupVals(fields, opts.GroupCols)
			curVals = make([][]string, len(opts.AggCols))
			for i, ac := range opts.AggCols {
				curVals[i] = append(curVals[i], fields[ac-1])
			}
			sawFirstRecord = true
			continue
		}
		if key == curKey {
			for i, ac := range opts.AggCols {
				curVals[i] = append(curVals[i], fields[ac-1])
			}
			continue
		}
		// New group: flush the previous one.
		if err := flush(); err != nil {
			return outCount, err
		}
		curKey = key
		curFirst = fields
		curGroupVals = pickGroupVals(fields, opts.GroupCols)
		curVals = make([][]string, len(opts.AggCols))
		for i, ac := range opts.AggCols {
			curVals[i] = append(curVals[i], fields[ac-1])
		}
	}
	if err := scanner.Err(); err != nil {
		return outCount, fmt.Errorf("error reading input: %w", err)
	}
	if err := flush(); err != nil {
		return outCount, err
	}
	if err := w.Flush(); err != nil {
		return outCount, fmt.Errorf("error flushing output: %w", err)
	}
	return outCount, nil
}

// pickGroupVals copies the grouping-column values out of a record's
// fields, preserving each value's original case. It is used to seed
// curGroupVals at the start of each group so the output row spells the
// grouping key the way the first row spelled it — matching upstream's
// `-ignorecase` first-seen-case behaviour.
func pickGroupVals(fields []string, cols []int) []string {
	out := make([]string, len(cols))
	for i, c := range cols {
		out[i] = fieldOr(fields, c)
	}
	return out
}

func fieldOr(fields []string, col int) string {
	if col < 1 || col > len(fields) {
		return ""
	}
	return fields[col-1]
}

func groupKey(fields []string, cols []int, ignoreCase bool) string {
	parts := make([]string, len(cols))
	for i, c := range cols {
		v := fieldOr(fields, c)
		if ignoreCase {
			v = strings.ToLower(v)
		}
		parts[i] = v
	}
	return strings.Join(parts, "\x00")
}
