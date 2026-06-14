// Package bedmerge provides functionality to merge overlapping or adjacent BED intervals.
package bedmerge

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// OutputFields specifies which fields to include in the output.
type OutputFields struct {
	Name     bool // Include name field
	Score    bool // Include score field
	Strand   bool // Include strand field
	Count    bool // Include count of merged intervals as name field
	BedGraph bool // Output in bedGraph format (chrom, start, end, score)
}

// MergeOptions contains options for merging BED intervals.
type MergeOptions struct {
	MaxDistance  int          // Maximum distance between intervals to merge (default: 0)
	StrandSpec   bool         // Merge only intervals on the same strand
	OutputFields OutputFields // Fields to include in output
	// StrandFilter ("-S <strand>"), when non-empty, must be "+" or "-": only
	// records on that strand are kept, and the survivors are merged
	// positionally (NOT per-strand), matching upstream bedtools merge -S.
	// It is mutually exclusive with StrandSpec ("-s").
	StrandFilter string
	Streaming    bool // Use streaming mode for large files
	// ColumnOps, when non-nil, requests bedtools-merge-style aggregation of
	// input columns (the -c/-o options). When set, output columns are
	// chrom, start, end followed by one aggregated value per requested column.
	ColumnOps *ColumnOps
}

// Merge reads BED intervals, sorts them, and merges overlapping/adjacent intervals.
// Returns the number of merged intervals.
func Merge(reader io.Reader, writer io.Writer, opts MergeOptions) (int, error) {
	if err := validateStrandOptions(opts); err != nil {
		return 0, err
	}
	reader, err := gffAwareReader(reader)
	if err != nil {
		return 0, err
	}
	// Column-aggregation mode (bedtools merge -c/-o style).
	if opts.ColumnOps != nil {
		return mergeWithColumnOps(reader, writer, opts)
	}

	// Use streaming mode if requested
	if opts.Streaming {
		return streamingMerge(reader, writer, opts)
	}

	// Read all intervals
	bedReader := bed.NewReader(reader)
	var intervals []*bed.Record

	for {
		record, err := bedReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, fmt.Errorf("error reading BED record: %w", err)
		}

		// -S single-strand filter: drop records on the other strand.
		if strandFiltered(record, opts) {
			continue
		}

		// If bedGraph input mode, parse the name field as score
		if opts.OutputFields.BedGraph && record.Name != "" && record.Score == 0 {
			// Try to parse name as score
			var score int
			if _, err := fmt.Sscanf(record.Name, "%d", &score); err == nil {
				record.Score = score
				record.Name = ""
			}
		}

		intervals = append(intervals, record)
	}

	if len(intervals) == 0 {
		return 0, nil
	}

	// Sort intervals by chromosome, then start position
	sort.Slice(intervals, func(i, j int) bool {
		if intervals[i].Chrom != intervals[j].Chrom {
			return intervals[i].Chrom < intervals[j].Chrom
		}
		if intervals[i].ChromStart != intervals[j].ChromStart {
			return intervals[i].ChromStart < intervals[j].ChromStart
		}
		return intervals[i].ChromEnd < intervals[j].ChromEnd
	})

	// Merge overlapping/adjacent intervals
	merged := mergeIntervals(intervals, opts)

	// Write merged intervals
	if err := writeIntervals(writer, merged, opts); err != nil {
		return 0, err
	}

	return len(merged), nil
}

// gffAwareReader auto-detects GFF or VCF input and, when found, returns a
// reader whose lines have been transformed into BED — so the rest of the merge
// pipeline (which is BED-only) merges those features just like upstream
// bedtools merge -i <gff|vcf>. A BED input is returned unchanged (the peeked
// header/first-data bytes are preserved). GFF is detected by the per-record
// heuristic (non-numeric source in column 2, numeric 1-based start/end in
// columns 4/5); VCF is detected by a `##fileformat=VCF` or `#CHROM` header.
func gffAwareReader(r io.Reader) (io.Reader, error) {
	br := bufio.NewReader(r)
	var header strings.Builder // consumed comment/header lines, replayed for BED
	isVCF := false
	for {
		line, err := br.ReadString('\n')
		trimmed := strings.TrimSpace(line)
		if strings.HasPrefix(trimmed, "##fileformat=VCF") ||
			strings.HasPrefix(trimmed, "#CHROM\tPOS") {
			isVCF = true
		}
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			header.WriteString(line)
			if err != nil {
				// EOF before any data line: replay what we consumed.
				return strings.NewReader(header.String()), nil
			}
			continue
		}
		// First data line. Replay the consumed header + this line + the rest.
		rest := io.MultiReader(strings.NewReader(header.String()), strings.NewReader(line), br)
		if isVCF {
			return transformVCF(rest)
		}
		if looksLikeGFF(strings.TrimRight(line, "\r\n")) {
			return transformGFF(rest)
		}
		return rest, nil
	}
}

// transformVCF rewrites every VCF data line into a BED3 line: the interval is
// [POS-1, POS-1+len(REF)), matching upstream bedtools merge -i <vcf>. Header
// lines ('#') are skipped.
func transformVCF(r io.Reader) (io.Reader, error) {
	var out strings.Builder
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 1<<20)
	for sc.Scan() {
		line := sc.Text()
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		f := strings.Split(line, "\t")
		if len(f) < 4 {
			return nil, fmt.Errorf("invalid VCF line: %q", line)
		}
		pos, err := strconv.Atoi(f[1])
		if err != nil {
			return nil, fmt.Errorf("invalid VCF POS %q: %w", f[1], err)
		}
		start := pos - 1
		end := start + len(f[3]) // len(REF)
		fmt.Fprintf(&out, "%s\t%d\t%d\n", f[0], start, end)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return strings.NewReader(out.String()), nil
}

// looksLikeGFF reports whether a data line is a GFF feature: at least 8
// tab-separated fields, a non-numeric column 2 (source), and numeric 1-based
// start/end in columns 4 and 5.
func looksLikeGFF(line string) bool {
	f := strings.Split(line, "\t")
	if len(f) < 8 {
		return false
	}
	if _, err := strconv.Atoi(f[1]); err == nil {
		return false // column 2 numeric -> BED start, not GFF
	}
	if _, err := strconv.Atoi(f[3]); err != nil {
		return false
	}
	if _, err := strconv.Atoi(f[4]); err != nil {
		return false
	}
	return true
}

// transformGFF rewrites every GFF data line into a BED6 line (comments and
// blank lines pass through). GFF is 1-based inclusive; BED is 0-based
// half-open, so start becomes col4-1 and end stays col5.
func transformGFF(r io.Reader) (io.Reader, error) {
	var out strings.Builder
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 1<<20)
	for sc.Scan() {
		line := sc.Text()
		trimmed := strings.TrimSpace(line)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		f := strings.Split(line, "\t")
		if len(f) < 8 {
			return nil, fmt.Errorf("invalid GFF line: %q", line)
		}
		start, err := strconv.Atoi(f[3])
		if err != nil {
			return nil, fmt.Errorf("invalid GFF start %q: %w", f[3], err)
		}
		end, err := strconv.Atoi(f[4])
		if err != nil {
			return nil, fmt.Errorf("invalid GFF end %q: %w", f[4], err)
		}
		name := f[2]   // feature type
		strand := f[6] // GFF strand column
		// The BED reader requires a numeric score; GFF scores are often "."
		// (and merge ignores the score anyway), so normalise to 0.
		score := "0"
		if _, err := strconv.Atoi(f[5]); err == nil {
			score = f[5]
		}
		fmt.Fprintf(&out, "%s\t%d\t%d\t%s\t%s\t%s\n", f[0], start-1, end, name, score, strand)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return strings.NewReader(out.String()), nil
}

// validateStrandOptions rejects an invalid -S argument and the illegal
// combination of -s and -S, mirroring upstream bedtools merge.
func validateStrandOptions(opts MergeOptions) error {
	if opts.StrandFilter == "" {
		return nil
	}
	if opts.StrandFilter != "+" && opts.StrandFilter != "-" {
		return fmt.Errorf("invalid strand for -S: %q (must be + or -)", opts.StrandFilter)
	}
	if opts.StrandSpec {
		return fmt.Errorf("-s and -S are mutually exclusive")
	}
	return nil
}

// strandFiltered reports whether record should be dropped under the -S
// single-strand filter. With no filter active it always returns false.
func strandFiltered(record *bed.Record, opts MergeOptions) bool {
	return opts.StrandFilter != "" && record.Strand != opts.StrandFilter
}

// mergedInterval represents a merged interval with metadata.
type mergedInterval struct {
	*bed.Record
	count int // Number of original intervals merged into this one
}

// mergeIntervals performs the actual merging of sorted intervals.
//
// When opts.StrandSpec is true ("bedtools merge -s"), the merge runs
// strictly per-strand: records with `.` or empty strand are DROPPED
// (matching upstream's FileRecordMergeMgr behaviour, which deletes
// UNKNOWN-strand records in stranded mode — see
// reference_code/bedtools/src/utils/FileRecordTools/FileRecordMergeMgr.cpp
// lines 47-58 and 96-129). The `+` and `-` groups are merged
// independently and then re-merged into a single (chrom, start, end)
// sorted output stream — which is what upstream `bedtools merge -s`
// emits and what `merge.t15` asserts.
func mergeIntervals(intervals []*bed.Record, opts MergeOptions) []mergedInterval {
	if len(intervals) == 0 {
		return nil
	}

	if opts.StrandSpec {
		// Split into per-strand buckets, merge each, then merge-sort the
		// two outputs back into one sorted stream.
		var plus, minus []*bed.Record
		for _, iv := range intervals {
			switch iv.Strand {
			case "+":
				plus = append(plus, iv)
			case "-":
				minus = append(minus, iv)
				// "." and "" are intentionally dropped (see doc comment).
			}
		}
		// Disable StrandSpec on the recursive calls — within each bucket
		// all records share the same strand and the simple single-pass
		// merge below is the right thing.
		inner := opts
		inner.StrandSpec = false
		mp := mergeIntervals(plus, inner)
		mm := mergeIntervals(minus, inner)
		return mergeSortedMergedByPos(mp, mm)
	}

	merged := []mergedInterval{}
	current := mergedInterval{
		Record: &bed.Record{
			Chrom:      intervals[0].Chrom,
			ChromStart: intervals[0].ChromStart,
			ChromEnd:   intervals[0].ChromEnd,
			Strand:     intervals[0].Strand,
			Name:       intervals[0].Name,
			Score:      intervals[0].Score,
		},
		count: 1,
	}

	for i := 1; i < len(intervals); i++ {
		interval := intervals[i]

		// Check if we can merge with current
		canMerge := false

		// Same chromosome?
		if interval.Chrom == current.Chrom {
			// Check if overlapping or within max distance
			if interval.ChromStart <= current.ChromEnd+opts.MaxDistance {
				canMerge = true
			}
		}

		if canMerge {
			// Extend current interval
			if interval.ChromEnd > current.ChromEnd {
				current.ChromEnd = interval.ChromEnd
			}
			current.count++

			// For bedGraph, we might want to average or sum scores
			// For now, we keep the first score
		} else {
			// Save current and start new interval
			merged = append(merged, current)
			current = mergedInterval{
				Record: &bed.Record{
					Chrom:      interval.Chrom,
					ChromStart: interval.ChromStart,
					ChromEnd:   interval.ChromEnd,
					Strand:     interval.Strand,
					Name:       interval.Name,
					Score:      interval.Score,
				},
				count: 1,
			}
		}
	}

	// Add the last interval
	merged = append(merged, current)

	return merged
}

// mergeSortedMergedByPos two-way-merges two already-sorted slices of
// merged intervals by (chrom, start, end). Used to recombine the `+` and
// `-` outputs of strand-aware merge into a single sorted stream.
func mergeSortedMergedByPos(a, b []mergedInterval) []mergedInterval {
	out := make([]mergedInterval, 0, len(a)+len(b))
	i, j := 0, 0
	less := func(x, y mergedInterval) bool {
		if x.Chrom != y.Chrom {
			return x.Chrom < y.Chrom
		}
		if x.ChromStart != y.ChromStart {
			return x.ChromStart < y.ChromStart
		}
		return x.ChromEnd < y.ChromEnd
	}
	for i < len(a) && j < len(b) {
		if less(a[i], b[j]) {
			out = append(out, a[i])
			i++
		} else {
			out = append(out, b[j])
			j++
		}
	}
	out = append(out, a[i:]...)
	out = append(out, b[j:]...)
	return out
}

// writeIntervals writes merged intervals according to output options.
func writeIntervals(writer io.Writer, merged []mergedInterval, opts MergeOptions) error {
	// Handle bedGraph format separately (requires score without name)
	if opts.OutputFields.BedGraph {
		return writeBedGraph(writer, merged)
	}

	bedWriter := bed.NewWriter(writer)

	for _, m := range merged {
		record := m.Record

		// Handle count field
		if opts.OutputFields.Count {
			record.Name = fmt.Sprintf("%d", m.count)
		} else if !opts.OutputFields.Name {
			record.Name = ""
		}

		if !opts.OutputFields.Score {
			record.Score = 0
		}

		if !opts.OutputFields.Strand {
			record.Strand = ""
		}

		if err := bedWriter.Write(record); err != nil {
			return fmt.Errorf("error writing BED record: %w", err)
		}
	}

	if err := bedWriter.Flush(); err != nil {
		return fmt.Errorf("error flushing output: %w", err)
	}

	return nil
}

// writeBedGraph writes intervals in bedGraph format (chrom, start, end, score).
func writeBedGraph(writer io.Writer, merged []mergedInterval) error {
	for _, m := range merged {
		line := fmt.Sprintf("%s\t%d\t%d\t%d\n",
			m.Chrom, m.ChromStart, m.ChromEnd, m.Score)
		if _, err := writer.Write([]byte(line)); err != nil {
			return fmt.Errorf("error writing bedGraph record: %w", err)
		}
	}
	return nil
}

// streamingMerge performs streaming merge for very large files.
// This uses a sliding window approach to avoid loading all intervals into memory.
func streamingMerge(reader io.Reader, writer io.Writer, opts MergeOptions) (int, error) {
	bedReader := bed.NewReader(reader)

	// For streaming mode, we need to read intervals in chunks per chromosome
	// We'll collect all intervals for one chromosome, merge them, write output,
	// then move to the next chromosome. This assumes input is roughly sorted by chromosome.

	var currentChrom string
	var chromIntervals []*bed.Record
	outputCount := 0

	flushChrom := func() error {
		if len(chromIntervals) == 0 {
			return nil
		}

		// Sort intervals for this chromosome
		sort.Slice(chromIntervals, func(i, j int) bool {
			if chromIntervals[i].ChromStart != chromIntervals[j].ChromStart {
				return chromIntervals[i].ChromStart < chromIntervals[j].ChromStart
			}
			return chromIntervals[i].ChromEnd < chromIntervals[j].ChromEnd
		})

		// Merge and write
		merged := mergeIntervals(chromIntervals, opts)
		if err := writeIntervals(writer, merged, opts); err != nil {
			return err
		}

		outputCount += len(merged)
		chromIntervals = nil
		return nil
	}

	for {
		record, err := bedReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, fmt.Errorf("error reading BED record: %w", err)
		}

		// If bedGraph input mode, parse the name field as score
		if opts.OutputFields.BedGraph && record.Name != "" && record.Score == 0 {
			// Try to parse name as score
			var score int
			if _, err := fmt.Sscanf(record.Name, "%d", &score); err == nil {
				record.Score = score
				record.Name = ""
			}
		}

		// -S single-strand filter: drop records on the other strand.
		if strandFiltered(record, opts) {
			continue
		}

		// If chromosome changed, flush previous chromosome
		if record.Chrom != currentChrom && currentChrom != "" {
			if err := flushChrom(); err != nil {
				return 0, err
			}
		}

		currentChrom = record.Chrom
		chromIntervals = append(chromIntervals, record)
	}

	// Flush last chromosome
	if err := flushChrom(); err != nil {
		return 0, err
	}

	return outputCount, nil
}

// Stats contains statistics about the merge operation.
type Stats struct {
	InputIntervals  int
	OutputIntervals int
	MergedCount     int
}

// MergeWithStats performs merge and returns detailed statistics.
func MergeWithStats(reader io.Reader, writer io.Writer, opts MergeOptions) (*Stats, error) {
	if err := validateStrandOptions(opts); err != nil {
		return nil, err
	}
	reader, err := gffAwareReader(reader)
	if err != nil {
		return nil, err
	}
	// Column-aggregation mode: report only the output count.
	if opts.ColumnOps != nil {
		count, err := mergeWithColumnOps(reader, writer, opts)
		if err != nil {
			return nil, err
		}
		return &Stats{OutputIntervals: count}, nil
	}

	// Streaming mode doesn't support detailed stats tracking currently
	if opts.Streaming {
		count, err := streamingMerge(reader, writer, opts)
		if err != nil {
			return nil, err
		}
		return &Stats{
			OutputIntervals: count,
		}, nil
	}

	// Read all intervals
	bedReader := bed.NewReader(reader)
	var intervals []*bed.Record

	for {
		record, err := bedReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("error reading BED record: %w", err)
		}

		// -S single-strand filter: drop records on the other strand.
		if strandFiltered(record, opts) {
			continue
		}

		// If bedGraph input mode, parse the name field as score
		if opts.OutputFields.BedGraph && record.Name != "" && record.Score == 0 {
			// Try to parse name as score
			var score int
			if _, err := fmt.Sscanf(record.Name, "%d", &score); err == nil {
				record.Score = score
				record.Name = ""
			}
		}

		intervals = append(intervals, record)
	}

	stats := &Stats{
		InputIntervals: len(intervals),
	}

	if len(intervals) == 0 {
		return stats, nil
	}

	// Sort intervals
	sort.Slice(intervals, func(i, j int) bool {
		if intervals[i].Chrom != intervals[j].Chrom {
			return intervals[i].Chrom < intervals[j].Chrom
		}
		if intervals[i].ChromStart != intervals[j].ChromStart {
			return intervals[i].ChromStart < intervals[j].ChromStart
		}
		return intervals[i].ChromEnd < intervals[j].ChromEnd
	})

	// Merge overlapping/adjacent intervals
	merged := mergeIntervals(intervals, opts)
	stats.OutputIntervals = len(merged)
	stats.MergedCount = stats.InputIntervals - stats.OutputIntervals

	// Write merged intervals
	if err := writeIntervals(writer, merged, opts); err != nil {
		return nil, err
	}

	return stats, nil
}
