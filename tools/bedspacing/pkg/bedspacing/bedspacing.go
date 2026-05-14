// Package bedspacing implements `bedtools spacing`: it walks a BED file (or
// any tab-delimited file with chrom/start/end in columns 1-3) and appends a
// column that reports the spacing between adjacent intervals on the same
// chromosome.
//
// Following upstream, the spacing token is one of:
//
//   - "."  — first interval on its chromosome.
//   - "-1" — overlaps the previous interval on the same chromosome.
//   - "0"  — exactly abuts the previous interval (prev.end == this.start).
//   - N    — otherwise: this.start - prev.end (positive gap in bases).
//
// The "previous" interval is tracked per-chromosome; ordering is whatever
// the input provides — bedspacing does not sort. To get the conventional
// genome-sorted spacing report, pipe a sorted BED in (`bedsort` or
// `sort -k1,1 -k2,2n`).
//
// The input columns are preserved verbatim; the spacing token is appended as
// a new trailing tab-separated column. Header lines (`#`, `track`, `browser`)
// and blank lines are passed through unchanged. Mirrors upstream's
// `bedtools spacing -i a.bed`.
package bedspacing

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Spacing reads BED-like records from r, computes the spacing to the
// previous interval on the same chromosome, and writes each record back to
// w with the spacing appended as a new column. Returns the number of data
// records written.
func Spacing(r io.Reader, w io.Writer) (int, error) {
	bw := bufio.NewWriter(w)
	defer bw.Flush()

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	// Tracks the most recent end seen per chromosome.
	prevEnd := make(map[string]int)
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
		if len(fields) < 3 {
			return written, fmt.Errorf("line %d: BED record must have at least 3 columns: %q", lineNo, raw)
		}
		chrom := fields[0]
		start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
		if err != nil {
			return written, fmt.Errorf("line %d: invalid chromStart %q: %v", lineNo, fields[1], err)
		}
		end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
		if err != nil {
			return written, fmt.Errorf("line %d: invalid chromEnd %q: %v", lineNo, fields[2], err)
		}

		var spacing string
		if pe, ok := prevEnd[chrom]; !ok {
			spacing = "."
		} else if start < pe {
			spacing = "-1"
		} else if start == pe {
			spacing = "0"
		} else {
			spacing = strconv.Itoa(start - pe)
		}

		if _, err := bw.WriteString(raw); err != nil {
			return written, err
		}
		if err := bw.WriteByte('\t'); err != nil {
			return written, err
		}
		if _, err := bw.WriteString(spacing); err != nil {
			return written, err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return written, err
		}
		// Upstream stores only the *immediately preceding* record per chrom
		// (see reference_code/bedtools/src/spacingFile/spacingFile.cpp). We
		// match that exactly: each record's spacing is computed against the
		// previous record's end on the same chromosome, then the previous
		// pointer advances. Row 5 of spacing.t01 (75-100 after 60-80)
		// therefore reports -1 (overlap), and row 6 (105-110) reports
		// 105 - 100 = 5.
		prevEnd[chrom] = end
		written++
	}
	if err := sc.Err(); err != nil {
		return written, err
	}
	return written, nil
}
