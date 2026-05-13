// Package bedsubtract implements bedtools-subtract-style interval subtraction
// for BED files.
//
// For each interval in A, any overlap with intervals in B is subtracted,
// emitting the remaining segments of A (in input order). When an interval
// in B punches a hole in the middle of an A interval, the A interval is
// split into multiple output rows. Original A columns (beyond chrom/start/end)
// are preserved through the split.
package bedsubtract

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// row is a tab-delimited BED-ish record that preserves the original raw fields
// of the input line, plus parsed chrom/start/end (and strand if present).
type row struct {
	fields []string // raw tab-separated columns from the input line
	chrom  string   // fields[0]
	start  int      // parsed from fields[1]
	end    int      // parsed from fields[2]
	strand string   // fields[5] if present, else ""
}

// length returns the size of the interval in base pairs.
func (r *row) length() int { return r.end - r.start }

// clone returns a copy of r with an independent fields slice.
func (r *row) clone() *row {
	cp := *r
	cp.fields = append([]string(nil), r.fields...)
	return &cp
}

// withSpan returns a copy of r with the given [start, end) interval, updating
// fields[1] and fields[2] to match.
func (r *row) withSpan(start, end int) *row {
	c := r.clone()
	c.start = start
	c.end = end
	c.fields[1] = strconv.Itoa(start)
	c.fields[2] = strconv.Itoa(end)
	return c
}

// readRows scans BED-style lines from r, skipping blank/comment/track/browser
// header lines. Each returned row preserves the raw column count of its input.
func readRows(r io.Reader) ([]*row, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var rows []*row
	lineNum := 0
	for sc.Scan() {
		lineNum++
		line := sc.Text()
		trimmed := strings.TrimSpace(line)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 3 {
			return nil, fmt.Errorf("line %d: BED record must have at least 3 fields, got %d", lineNum, len(fields))
		}
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("line %d: invalid start %q: %v", lineNum, fields[1], err)
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return nil, fmt.Errorf("line %d: invalid end %q: %v", lineNum, fields[2], err)
		}
		if end < start {
			return nil, fmt.Errorf("line %d: end < start (%d < %d)", lineNum, end, start)
		}
		rr := &row{fields: fields, chrom: fields[0], start: start, end: end}
		if len(fields) >= 6 {
			rr.strand = fields[5]
		}
		rows = append(rows, rr)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return rows, nil
}

// writeRow writes r as a tab-separated line followed by a newline.
func writeRow(w *bufio.Writer, r *row) error {
	for i, f := range r.fields {
		if i > 0 {
			if err := w.WriteByte('\t'); err != nil {
				return err
			}
		}
		if _, err := w.WriteString(f); err != nil {
			return err
		}
	}
	return w.WriteByte('\n')
}
