// Package bedoverlap computes the amount of overlap (positive) or distance
// (negative) between two intervals described by four user-specified columns of
// each input line, appending the result as a new trailing column. It mirrors
// the behaviour of `bedtools overlap` (aka getOverlap).
//
// For each input line the four 1-based columns named by Cols (start1, end1,
// start2, end2) are parsed as integers and the overlap is computed as
//
//	min(end1, end2) - max(start1, start2)
//
// A positive result is the number of overlapping bases; a non-positive result
// is the (negative) gap between the two intervals. The original line is emitted
// verbatim with the overlap appended after a tab, exactly as upstream does.
package bedoverlap

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Cols identifies the four 1-based column indices that hold the start and end
// coordinates of the two intervals on each input line.
type Cols struct {
	// S1 is the 1-based column of the first interval's start.
	S1 int
	// E1 is the 1-based column of the first interval's end.
	E1 int
	// S2 is the 1-based column of the second interval's start.
	S2 int
	// E2 is the 1-based column of the second interval's end.
	E2 int
}

// ParseCols parses a comma-separated "s1,e1,s2,e2" specification (as supplied
// to the -cols flag) into a Cols value. It returns an error unless exactly four
// integer fields are present, matching upstream's requirement.
func ParseCols(spec string) (Cols, error) {
	parts := strings.Split(spec, ",")
	if len(parts) != 4 {
		return Cols{}, fmt.Errorf("Please specify 4, comma-separated position columns.")
	}
	vals := make([]int, 4)
	for i, p := range parts {
		n, err := strconv.Atoi(strings.TrimSpace(p))
		if err != nil {
			return Cols{}, fmt.Errorf("invalid column %q: %v", p, err)
		}
		vals[i] = n
	}
	return Cols{S1: vals[0], E1: vals[1], S2: vals[2], E2: vals[3]}, nil
}

// overlap returns min(e1, e2) - max(s1, s2): the number of overlapping bases
// when positive, or the (negative) distance between the intervals otherwise.
func overlap(s1, e1, s2, e2 int64) int64 {
	end := e1
	if e2 < end {
		end = e2
	}
	start := s1
	if s2 > start {
		start = s2
	}
	return end - start
}

// parseInt mimics C strtol(base 10): it parses the leading integer prefix of s
// and reports whether any digits were consumed. A leading sign and surrounding
// whitespace are tolerated. Upstream rejects a column only when no integer
// prefix is present (its strtol end pointer would equal the start pointer).
func parseInt(s string) (int64, bool) {
	i := 0
	n := len(s)
	for i < n && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == '\r' || s[i] == '\v' || s[i] == '\f') {
		i++
	}
	start := i
	if i < n && (s[i] == '+' || s[i] == '-') {
		i++
	}
	digitsStart := i
	for i < n && s[i] >= '0' && s[i] <= '9' {
		i++
	}
	if i == digitsStart {
		// No digits consumed: strtol would leave the end pointer at the
		// original start, which upstream treats as non-numeric.
		return 0, false
	}
	v, err := strconv.ParseInt(s[start:i], 10, 64)
	if err != nil {
		return 0, false
	}
	return v, true
}

// Overlap reads input lines from in, computes the overlap/distance for each
// line using the columns in cols, and writes each original line followed by a
// tab and the overlap value to out.
//
// Lines that tokenize (on tabs) to a single field or fewer are skipped without
// output, matching upstream. If any of the four target columns on a processed
// line lacks a leading integer, Overlap returns an error describing the 1-based
// line number, mirroring upstream's "non-numeric" diagnostic and exit.
func Overlap(in io.Reader, out io.Writer, cols Cols) error {
	scanner := bufio.NewScanner(in)
	scanner.Buffer(make([]byte, 64*1024), 64*1024*1024)
	bw := bufio.NewWriter(out)
	defer bw.Flush()

	lineNum := 0
	for scanner.Scan() {
		line := scanner.Text()
		lineNum++
		fields := strings.Split(line, "\t")
		if len(fields) <= 1 {
			continue
		}
		s1, ok1 := getCol(fields, cols.S1)
		e1, ok2 := getCol(fields, cols.E1)
		s2, ok3 := getCol(fields, cols.S2)
		e2, ok4 := getCol(fields, cols.E2)
		if !ok1 || !ok2 || !ok3 || !ok4 {
			// Upstream emits the message followed by a blank line ("endl <<
			// endl"); the trailing "\n" here combines with main's "%v\n" to
			// reproduce that exactly.
			return fmt.Errorf("One of your columns appears to be non-numeric at line %d. Exiting...\n", lineNum)
		}
		ov := overlap(s1, e1, s2, e2)
		if _, err := fmt.Fprintf(bw, "%s\t%d\n", line, ov); err != nil {
			return err
		}
	}
	return scanner.Err()
}

// getCol returns the integer value of the 1-based column col from fields and
// whether it parsed. An out-of-range column is reported as non-numeric, which
// is how upstream's strtol-on-empty-string behaves.
func getCol(fields []string, col int) (int64, bool) {
	if col < 1 || col > len(fields) {
		return 0, false
	}
	return parseInt(fields[col-1])
}
