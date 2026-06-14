// Package bedcomplement emits the genomic regions NOT covered by a sorted BED
// file. It mirrors the behaviour of `bedtools complement`.
//
// For each chromosome in the chrom-sizes file the complement is the set of
// half-open intervals between consecutive sorted input intervals, plus the
// leading gap [0, first.start) and trailing gap [last.end, chromSize).
// Chromosomes that appear in the chrom-sizes file but have no intervals in
// the input emit a single record `chrom\t0\tchromSize`.
//
// The input is required to be sorted on (chrom, start); the package detects
// out-of-order input on the fly and aborts with an error.
package bedcomplement

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"
)

// ChromSizes maps chromosome name to its total length in bases.
type ChromSizes map[string]int

// ReadChromSizes parses a chrom-sizes file (one `chrom<TAB>size` per line).
// It also accepts samtools-style .fai files (uses the first two whitespace-
// separated columns). It returns the size map plus the chromosome order in
// which the entries were observed (used to make the output deterministic).
// Blank lines and comments (`#`) are skipped.
func ReadChromSizes(r io.Reader) (ChromSizes, []string, error) {
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 64*1024), 16*1024*1024)
	sizes := make(ChromSizes)
	var order []string
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, nil, fmt.Errorf("chrom-sizes line %q must have at least 2 fields", line)
		}
		size, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, nil, fmt.Errorf("invalid size %q for chromosome %q: %v", fields[1], fields[0], err)
		}
		if size < 0 {
			return nil, nil, fmt.Errorf("negative size %d for chromosome %q", size, fields[0])
		}
		if _, dup := sizes[fields[0]]; !dup {
			order = append(order, fields[0])
		}
		sizes[fields[0]] = size
	}
	if err := scanner.Err(); err != nil {
		return nil, nil, err
	}
	return sizes, order, nil
}

// Complement reads sorted BED records from in and writes the complementary
// intervals (BED3) to out using the chromosome sizes in sizes. The
// `chromOrder` argument controls the order in which chromosomes are emitted;
// any chromosomes in sizes that are not listed in chromOrder are appended in
// lexicographic order.
//
// Intervals on chromosomes that are not in sizes are skipped (a single
// warning per chromosome is emitted to warn). The function returns an error
// if the input is not sorted by (chrom, start) or if a record is malformed.
//
// When limitToInput is true (upstream `bedtools complement -L`), only
// chromosomes that had at least one input record are emitted; otherwise every
// chromosome in sizes is emitted (a chromosome with no input records yields a
// single full-length gap).
//
// Complement returns the number of complementary intervals written.
func Complement(in io.Reader, out io.Writer, warn io.Writer, sizes ChromSizes, chromOrder []string, limitToInput bool) (int, error) {
	scanner := bufio.NewScanner(in)
	scanner.Buffer(make([]byte, 64*1024), 16*1024*1024)

	// Group intervals by chromosome, in the order each chromosome is first
	// seen in the (sorted) input. Within a chromosome, intervals are required
	// to be ordered by start.
	grouped := make(map[string][]interval)
	groupOrder := []string{}
	missing := make(map[string]bool)
	var prevChrom string
	prevHasIvl := false
	var prevEnd, prevStart int
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
			return 0, fmt.Errorf("BED record must have at least 3 fields, got %d: %q", len(fields), raw)
		}
		start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
		if err != nil {
			return 0, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
		}
		end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
		if err != nil {
			return 0, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
		}
		if end < start {
			return 0, fmt.Errorf("invalid interval %s\t%d\t%d: end < start", fields[0], start, end)
		}
		chrom := fields[0]
		if _, ok := sizes[chrom]; !ok {
			if !missing[chrom] {
				missing[chrom] = true
				if warn != nil {
					fmt.Fprintf(warn, "warning: chromosome %q not in genome file; skipping intervals on it\n", chrom)
				}
			}
			continue
		}
		if chrom != prevChrom {
			// Starting a new chromosome. Verify this chromosome wasn't seen
			// before (which would mean the input is not chrom-grouped, i.e.
			// unsorted).
			if _, seen := grouped[chrom]; seen {
				return 0, fmt.Errorf("input not sorted: chromosome %q reappears after %q", chrom, prevChrom)
			}
			groupOrder = append(groupOrder, chrom)
			prevHasIvl = false
		} else if prevHasIvl {
			// Same chromosome: require non-decreasing start (ties on start are
			// fine; we re-sort within chrom defensively below to handle them).
			if start < prevStart {
				return 0, fmt.Errorf("input not sorted: %s\t%d\t%d follows %s\t%d\t%d",
					chrom, start, end, chrom, prevStart, prevEnd)
			}
		}
		grouped[chrom] = append(grouped[chrom], interval{start: start, end: end})
		prevChrom = chrom
		prevStart = start
		prevEnd = end
		prevHasIvl = true
	}
	if err := scanner.Err(); err != nil {
		return 0, err
	}

	// Build the output chromosome order: chromOrder first, then any other
	// chromosomes from sizes in lexicographic order. This makes the output
	// deterministic regardless of how chromOrder was supplied.
	emitOrder := buildEmitOrder(sizes, chromOrder)
	// Under -L only emit chromosomes that actually had input records.
	if limitToInput {
		filtered := emitOrder[:0:0]
		for _, chrom := range emitOrder {
			if _, ok := grouped[chrom]; ok {
				filtered = append(filtered, chrom)
			}
		}
		emitOrder = filtered
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()
	written := 0
	for _, chrom := range emitOrder {
		chromSize := sizes[chrom]
		ivls := grouped[chrom]
		// Within a chromosome, sort intervals by start (defensive: input was
		// already validated to be non-decreasing, but ties on start need
		// secondary ordering by end for correct merging below).
		sort.SliceStable(ivls, func(i, j int) bool {
			if ivls[i].start != ivls[j].start {
				return ivls[i].start < ivls[j].start
			}
			return ivls[i].end < ivls[j].end
		})
		// Merge any overlapping or touching intervals so the complement gaps
		// don't include negative or empty regions.
		merged := mergeIntervals(ivls)
		// Emit complementary gaps clipped to [0, chromSize].
		prev := 0
		for _, iv := range merged {
			s := iv.start
			if s < 0 {
				s = 0
			}
			if s > chromSize {
				s = chromSize
			}
			e := iv.end
			if e > chromSize {
				e = chromSize
			}
			if s > prev {
				if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\n", chrom, prev, s); err != nil {
					return written, err
				}
				written++
			}
			if e > prev {
				prev = e
			}
		}
		if prev < chromSize {
			if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\n", chrom, prev, chromSize); err != nil {
				return written, err
			}
			written++
		}
	}
	return written, nil
}

// interval is a half-open BED interval used internally.
type interval struct{ start, end int }

// mergeIntervals collapses overlapping or touching intervals (assumes the
// slice is already sorted by start). It returns a new slice.
func mergeIntervals(in []interval) []interval {
	if len(in) == 0 {
		return nil
	}
	out := []interval{in[0]}
	for _, iv := range in[1:] {
		top := &out[len(out)-1]
		if iv.start <= top.end {
			if iv.end > top.end {
				top.end = iv.end
			}
		} else {
			out = append(out, iv)
		}
	}
	return out
}

// buildEmitOrder returns the chromosome emission order: chromOrder first,
// then any other chromosomes in sizes in lexicographic order. Duplicates in
// chromOrder are de-duplicated.
func buildEmitOrder(sizes ChromSizes, chromOrder []string) []string {
	seen := make(map[string]bool, len(sizes))
	out := make([]string, 0, len(sizes))
	for _, c := range chromOrder {
		if _, ok := sizes[c]; !ok {
			continue
		}
		if seen[c] {
			continue
		}
		seen[c] = true
		out = append(out, c)
	}
	var rest []string
	for c := range sizes {
		if !seen[c] {
			rest = append(rest, c)
		}
	}
	sort.Strings(rest)
	out = append(out, rest...)
	return out
}
