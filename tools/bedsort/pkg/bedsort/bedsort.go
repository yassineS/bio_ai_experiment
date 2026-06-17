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

	"github.com/yassineS/bio_ai_experiment/pkg/cppsort"
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
	// Header, when true, preserves leading header lines (`#`-prefixed comments,
	// `track ` directives, and `browser ` directives) at the top of the
	// output, before the sorted body. Matches upstream `bedtools sort -header`.
	// Header lines are emitted verbatim in the order they appeared in the
	// input. Header lines that appear mid-file (after at least one record)
	// are dropped, matching upstream's "leading header only" semantics.
	Header bool
}

// record is the parsed view of one BED line used by the sort.
type record struct {
	line  string // original line (without trailing newline)
	chrom string
	start int
	end   int
	score int
	// scoreStr is the raw column-5 text. Upstream bedtools stores BED.score as a
	// std::string and sortByScore{Asc,Desc} compares those strings directly, so
	// the score sort is LEXICOGRAPHIC ("10" < "100" < "101" < "2"), not numeric.
	// We keep the raw string to reproduce that ordering byte-for-byte.
	scoreStr string
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
	recs, _, err := readAll(r, false)
	return recs, err
}

// readAll reads BED records from r. When keepHeader is true, the leading
// `#`-prefix / `track ` / `browser ` lines (before the first data record) are
// returned verbatim in headers, in input order. Header-style lines that
// appear after the first data record are dropped, matching upstream
// `bedtools sort -header`.
func readAll(r io.Reader, keepHeader bool) (recs []record, headers []string, err error) {
	scanner := bufio.NewScanner(r)
	// Allow long BED lines.
	scanner.Buffer(make([]byte, 64*1024), 16*1024*1024)
	seenData := false
	for scanner.Scan() {
		raw := scanner.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" {
			continue
		}
		if strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") ||
			strings.HasPrefix(trimmed, "browser") {
			if keepHeader && !seenData {
				headers = append(headers, raw)
			}
			continue
		}
		fields := strings.Split(raw, "\t")
		if len(fields) < 3 {
			return nil, nil, fmt.Errorf("BED record must have at least 3 fields, got %d: %q", len(fields), raw)
		}
		start, errStart := strconv.Atoi(strings.TrimSpace(fields[1]))
		if errStart != nil {
			return nil, nil, fmt.Errorf("invalid chromStart %q: %v", fields[1], errStart)
		}
		end, errEnd := strconv.Atoi(strings.TrimSpace(fields[2]))
		if errEnd != nil {
			return nil, nil, fmt.Errorf("invalid chromEnd %q: %v", fields[2], errEnd)
		}
		rec := record{
			line:  raw,
			chrom: fields[0],
			start: start,
			end:   end,
		}
		if len(fields) >= 5 {
			rec.scoreStr = fields[4]
			if s, errScore := strconv.Atoi(strings.TrimSpace(fields[4])); errScore == nil {
				rec.score = s
				rec.hasScore = true
			}
		}
		recs = append(recs, rec)
		seenData = true
	}
	if err := scanner.Err(); err != nil {
		return nil, nil, err
	}
	return recs, headers, nil
}

// Sort sorts records in place according to opts.
//
// Sort faithfully reproduces upstream bedtools sort (the legacy sortBed tool),
// which works in two stages:
//
//  1. Records are bucketed by chromosome (a std::map, so chromosomes come out
//     in lexicographic order, or in opts.ChromOrder order when a faidx/genome
//     is supplied) and each bucket is sorted by chromStart ALONE. That sort
//     (upstream's sortByStart) compares only the start coordinate, so records
//     with an equal (chrom, start) key retain their original input order. It
//     does NOT use chromEnd as a tie-break.
//
//  2. For the size/score modes a second sort is applied on top of that
//     start-ordered arrangement. Upstream's size-descending and score
//     comparators compare only the size (or score) with no secondary key, so
//     equal-keyed records keep the order established in stage 1 (i.e. input
//     order on equal (chrom, start)). The ascending size comparator adds a
//     (chrom, start) tie-break, which is already the stage-1 order, so the two
//     are equivalent here.
//
// Both stages use a stable sort so the input-order tie-break that upstream
// gets from its start-only ordering is reproduced byte-for-byte. The previous
// implementation incorrectly used chromEnd as a tie-break, which diverged from
// upstream for the default, -sizeD, -chrThenSizeD, -chrThenScoreA and
// -chrThenScoreD modes whenever the input was not already ordered by end.
func Sort(records []record, opts Options) {
	chromRank := buildChromRank(opts.ChromOrder)
	// cmpChrom compares two records by chromosome, honouring an explicit
	// ChromOrder (faidx/genome) when present and falling back to lexicographic
	// order. Chromosomes absent from ChromOrder sort after the listed ones, in
	// lexicographic order amongst themselves.
	cmpChrom := func(a, b record) int {
		if chromRank != nil {
			ra, oka := chromRank[a.chrom]
			rb, okb := chromRank[b.chrom]
			switch {
			case oka && okb:
				switch {
				case ra < rb:
					return -1
				case ra > rb:
					return 1
				default:
					return 0
				}
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

	// Stage 1: bucket records by chromosome (preserving input order within each
	// bucket via a stable group), then run libstdc++'s std::sort on each bucket
	// with a start-ONLY comparator. This reproduces upstream sortBed exactly:
	// loadBedFileIntoMapNoBin pushes records into a std::map<chrom, vector<BED>>
	// (so the per-chrom vector is in input order) and then calls
	// `std::sort(vec, sortByStart)`. std::sort is NOT stable, so records that
	// tie on start come out in introsort's order — neither input order nor any
	// secondary key — which Go's sort.Slice does not reproduce. cppsort.Sort is
	// a byte-faithful port of that introsort, and startCorrected mirrors
	// sortByStart's zeroLength (start==end) +1 adjustment.
	startCorrected := func(r record) int {
		if r.start == r.end {
			return r.start + 1
		}
		return r.start
	}
	startLess := func(a, b record) bool { return startCorrected(a) < startCorrected(b) }

	// Stable group by chromosome to recover each chrom's input-order vector,
	// then emit chromosomes in cmpChrom order with each bucket introsorted.
	sortByChromBucketsStart := func() {
		buckets := map[string][]record{}
		var order []record // one representative per chrom, in first-seen order
		for _, r := range records {
			if _, ok := buckets[r.chrom]; !ok {
				order = append(order, r)
			}
			buckets[r.chrom] = append(buckets[r.chrom], r)
		}
		// Order the chromosome keys exactly as upstream's std::map iteration
		// (or the genome/faidx order) does, via cmpChrom.
		sort.SliceStable(order, func(i, j int) bool { return cmpChrom(order[i], order[j]) < 0 })
		out := records[:0]
		for _, rep := range order {
			b := buckets[rep.chrom]
			cppsort.Sort(b, startLess)
			out = append(out, b...)
		}
	}
	sortByChromBucketsStart()

	if opts.Mode == ModeChrom {
		return
	}

	// Stage 2: apply the mode key with libstdc++ std::sort, using upstream's
	// exact comparators (utils/bedFile/bedFile.cpp). Each is a partial order:
	//   - sortBySizeAsc: by length asc, ties broken by byChromThenStart (chrom
	//     then start, lexicographic — NOT genome order).
	//   - sortBySizeDesc / sortByScore{Asc,Desc}: length / score only, no
	//     secondary key (so equal-key records land in introsort's order).
	// -sizeA/-sizeD sort the GLOBAL masterList (stage-1 output concatenated
	// across chromosomes); -chrThenSize*/-chrThenScore* sort EACH chromosome's
	// contiguous run in place. cppsort reproduces std::sort byte-for-byte.
	size := func(r record) int { return r.end - r.start }
	byChromThenStart := func(a, b record) bool {
		if a.chrom != b.chrom {
			return a.chrom < b.chrom
		}
		return a.start < b.start
	}
	sizeAsc := func(a, b record) bool {
		la, lb := size(a), size(b)
		if la < lb {
			return true
		}
		if la > lb {
			return false
		}
		return byChromThenStart(a, b)
	}
	sizeDesc := func(a, b record) bool { return size(a) > size(b) }
	scoreAsc := func(a, b record) bool { return a.scoreStr < b.scoreStr }
	scoreDesc := func(a, b record) bool { return a.scoreStr > b.scoreStr }

	// perChrom runs less over each maximal contiguous same-chromosome run of
	// records (stage 1 already left them grouped in chromosome order).
	perChrom := func(less func(a, b record) bool) {
		i := 0
		for i < len(records) {
			j := i + 1
			for j < len(records) && records[j].chrom == records[i].chrom {
				j++
			}
			cppsort.Sort(records[i:j], less)
			i = j
		}
	}

	switch opts.Mode {
	case ModeSizeA:
		cppsort.Sort(records, sizeAsc)
	case ModeSizeD:
		cppsort.Sort(records, sizeDesc)
	case ModeChrThenSizeA:
		perChrom(sizeAsc)
	case ModeChrThenSizeD:
		perChrom(sizeDesc)
	case ModeChrThenScoreA:
		perChrom(scoreAsc)
	case ModeChrThenScoreD:
		perChrom(scoreDesc)
	}
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
//
// When opts.Header is true, leading `#` / `track ` / `browser ` lines are
// emitted verbatim before the sorted body, mirroring upstream
// `bedtools sort -header`.
func Run(r io.Reader, w io.Writer, opts Options) error {
	records, headers, err := readAll(r, opts.Header)
	if err != nil {
		return err
	}
	Sort(records, opts)
	if opts.Header && len(headers) > 0 {
		bw := bufio.NewWriter(w)
		for _, h := range headers {
			if _, err := bw.WriteString(h); err != nil {
				return err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return err
			}
		}
		if err := bw.Flush(); err != nil {
			return err
		}
	}
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
