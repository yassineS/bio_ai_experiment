// Package bedsort sorts BED intervals using various sort modes.
//
// The default sort is by chromosome (lexicographic), then by chromStart
// ascending, then by chromEnd ascending. Alternative modes mirror the
// upstream bedtools sort tool: by interval size (asc/desc), by chromosome
// then by size or score, and by an external faidx/genome file that fixes the
// chromosome order.
//
// Records are kept as raw text lines so the full set of input columns (BED3,
// BED6, BED12, or any number of extra fields) round-trips through the sort
// unchanged.
package bedsort

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"
)

// SortMode selects the ordering applied to records.
type SortMode int

const (
	// ModeChrom is the default: by chromosome (lexicographic), then start asc,
	// then end asc.
	ModeChrom SortMode = iota
	// ModeSizeA sorts by interval size ascending only.
	ModeSizeA
	// ModeSizeD sorts by interval size descending only.
	ModeSizeD
	// ModeChrThenSizeA sorts by chromosome then by interval size ascending.
	ModeChrThenSizeA
	// ModeChrThenSizeD sorts by chromosome then by interval size descending.
	ModeChrThenSizeD
	// ModeChrThenScoreA sorts by chromosome then by score (column 5) ascending.
	ModeChrThenScoreA
	// ModeChrThenScoreD sorts by chromosome then by score (column 5) descending.
	ModeChrThenScoreD
)

// Options bundles the sort knobs for Sort.
type Options struct {
	// Mode selects the ordering. Default ModeChrom.
	Mode SortMode
	// ChromOrder, if non-empty, gives an explicit chromosome ordering. Records
	// on chromosomes not present in ChromOrder are sorted after the listed
	// ones in lexicographic order. Used by ModeChrom and the ChrThen* modes.
	ChromOrder []string
}

// record is the parsed view of one BED line used by the sort.
type record struct {
	line  string // original line (without trailing newline)
	chrom string
	start int
	end   int
	score int
	// hasScore is false when the line has fewer than 5 columns (or the score
	// column did not parse as an integer). Records without a score sort as if
	// they had score 0 in the ChrThenScore* modes but we still try to keep the
	// ordering stable.
	hasScore bool
}

// ReadAll reads every non-empty, non-comment BED line from r.
//
// Comment and header lines (those starting with "#", "track", or "browser")
// are filtered out (matching the behaviour of upstream bedtools sort, which
// discards them). The returned slice contains the records in input order; the
// caller is expected to feed it to Sort.
func ReadAll(r io.Reader) ([]record, error) {
	scanner := bufio.NewScanner(r)
	// Allow long BED lines.
	scanner.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var out []record
	for scanner.Scan() {
		raw := scanner.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" {
			continue
		}
		if strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") ||
			strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(raw, "\t")
		if len(fields) < 3 {
			return nil, fmt.Errorf("BED record must have at least 3 fields, got %d: %q", len(fields), raw)
		}
		start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
		if err != nil {
			return nil, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
		}
		end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
		if err != nil {
			return nil, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
		}
		rec := record{
			line:  raw,
			chrom: fields[0],
			start: start,
			end:   end,
		}
		if len(fields) >= 5 {
			if s, err := strconv.Atoi(strings.TrimSpace(fields[4])); err == nil {
				rec.score = s
				rec.hasScore = true
			}
		}
		out = append(out, rec)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// Sort sorts records in place according to opts.
//
// Sort is stable, and size/score modes additionally break ties on the default
// (chrom, start, end) ordering so the output matches upstream bedtools sort
// deterministically rather than depending on input order.
func Sort(records []record, opts Options) {
	chromRank := buildChromRank(opts.ChromOrder)
	cmpChrom := func(a, b record) int {
		if chromRank != nil {
			ra, oka := chromRank[a.chrom]
			rb, okb := chromRank[b.chrom]
			switch {
			case oka && okb:
				if ra != rb {
					if ra < rb {
						return -1
					}
					return 1
				}
				return 0
			case oka && !okb:
				return -1
			case !oka && okb:
				return 1
			}
		}
		if a.chrom == b.chrom {
			return 0
		}
		if a.chrom < b.chrom {
			return -1
		}
		return 1
	}
	// tieBreak compares two records on the default (chrom asc, start asc, end
	// asc) ordering. It is used to break ties for size/score sorts so the
	// output matches upstream bedtools sort, which secondary-sorts ties
	// deterministically rather than relying on input order.
	tieBreak := func(a, b record) int {
		if c := cmpChrom(a, b); c != 0 {
			return c
		}
		if a.start != b.start {
			if a.start < b.start {
				return -1
			}
			return 1
		}
		if a.end != b.end {
			if a.end < b.end {
				return -1
			}
			return 1
		}
		return 0
	}
	less := func(i, j int) bool {
		a, b := records[i], records[j]
		switch opts.Mode {
		case ModeSizeA:
			la, lb := a.end-a.start, b.end-b.start
			if la != lb {
				return la < lb
			}
			return tieBreak(a, b) < 0
		case ModeSizeD:
			la, lb := a.end-a.start, b.end-b.start
			if la != lb {
				return la > lb
			}
			return tieBreak(a, b) < 0
		case ModeChrThenSizeA:
			if c := cmpChrom(a, b); c != 0 {
				return c < 0
			}
			la, lb := a.end-a.start, b.end-b.start
			if la != lb {
				return la < lb
			}
			return tieBreak(a, b) < 0
		case ModeChrThenSizeD:
			if c := cmpChrom(a, b); c != 0 {
				return c < 0
			}
			la, lb := a.end-a.start, b.end-b.start
			if la != lb {
				return la > lb
			}
			return tieBreak(a, b) < 0
		case ModeChrThenScoreA:
			if c := cmpChrom(a, b); c != 0 {
				return c < 0
			}
			if a.score != b.score {
				return a.score < b.score
			}
			return tieBreak(a, b) < 0
		case ModeChrThenScoreD:
			if c := cmpChrom(a, b); c != 0 {
				return c < 0
			}
			if a.score != b.score {
				return a.score > b.score
			}
			return tieBreak(a, b) < 0
		default: // ModeChrom
			if c := cmpChrom(a, b); c != 0 {
				return c < 0
			}
			if a.start != b.start {
				return a.start < b.start
			}
			return a.end < b.end
		}
	}
	sort.SliceStable(records, less)
}

// Write writes records to w as their original input lines, one per line.
func Write(w io.Writer, records []record) error {
	bw := bufio.NewWriter(w)
	for _, r := range records {
		if _, err := bw.WriteString(r.line); err != nil {
			return err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
	}
	return bw.Flush()
}

// Run reads BED records from r, sorts them according to opts, and writes the
// sorted output to w. It is the convenience entry point used by the CLI.
func Run(r io.Reader, w io.Writer, opts Options) error {
	records, err := ReadAll(r)
	if err != nil {
		return err
	}
	Sort(records, opts)
	return Write(w, records)
}

// ReadFaidx parses a .fai (samtools faidx) or genome (chrom-sizes) file and
// returns the chromosome order it implies. Only the first whitespace-separated
// column on each non-empty, non-comment line is used, so the same routine
// works for `chrom\tlen` chrom-sizes files and full .fai files (which add more
// columns after the length).
func ReadFaidx(r io.Reader) ([]string, error) {
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var order []string
	seen := make(map[string]struct{})
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		// Split on any whitespace so we accept both tab- and space-separated
		// columns.
		fields := strings.Fields(line)
		if len(fields) == 0 {
			continue
		}
		chrom := fields[0]
		if _, ok := seen[chrom]; ok {
			continue
		}
		seen[chrom] = struct{}{}
		order = append(order, chrom)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return order, nil
}

// buildChromRank returns a map from chromosome name to its 0-based rank in
// order. It returns nil when order is empty so callers can use lexicographic
// ordering as the default.
func buildChromRank(order []string) map[string]int {
	if len(order) == 0 {
		return nil
	}
	m := make(map[string]int, len(order))
	for i, c := range order {
		if _, ok := m[c]; !ok {
			m[c] = i
		}
	}
	return m
}
