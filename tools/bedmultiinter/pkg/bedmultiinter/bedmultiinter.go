// Package bedmultiinter implements `bedtools multiinter`: a multi-way
// intersection across N BED files. For each event-segment along each
// chromosome, the package emits a row giving the segment's bounds, the
// number of input files contributing, a comma-list of contributing file
// names, and N 0/1 indicator columns (one per input).
//
// Each input file may be BED, VCF, or GFF: the format is autodetected from
// the first non-header data line (matching upstream BedFile), and
// coordinates are converted to BED-style 0-based half-open spans.
//
// Algorithm: a per-chromosome event sweep. For each chromosome we read
// every input file's intervals (merging adjacent / overlapping records
// within a single file, as upstream does), then sort the START and END
// events across all files. Walking left-to-right, the active set of
// files is the set whose START has been seen but whose END has not. A
// segment is the half-open interval between two consecutive event
// coordinates; we emit one row per segment.
//
// With `Empty=true` and a genome file, the gaps at the chromosome head
// and tail are also emitted as `num=0 / list=none` rows.
//
// With `Cluster=true`, adjacent same-active-set segments are collapsed
// into one row (mirroring upstream `-cluster`).
package bedmultiinter

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"
)

// Options configures Run.
type Options struct {
	// Names is the column label per input file. When empty, list entries
	// are 1-based integer indices and the header row (if requested) uses
	// filenames as supplied.
	Names []string
	// Filenames are the path strings used in the default header. They
	// are needed independently from Names because, in upstream, when no
	// -names are given, the header still shows the filenames but the
	// `list` field uses indices.
	Filenames []string
	// Empty, when true, emits regions with num=0 (no file contributing)
	// at the leading and trailing edges of each chrom. Requires
	// ChromSizes to be populated.
	Empty bool
	// ChromSizes maps chromosome -> length, sourced from a genome file.
	// Required when Empty=true; otherwise unused.
	ChromSizes map[string]int
	// Cluster, when true, collapses adjacent segments that share the
	// same active set into a single row.
	Cluster bool
	// Header, when true, emits a header line before the data rows.
	Header bool
	// Filler is the indicator for "file not contributing" cells. Upstream
	// default is "0"; some workflows use "N/A". Use Filler="" to take the
	// upstream default.
	Filler string
}

// Run reads each B file from bRs and writes the multi-intersection to w.
// Returns the number of data rows emitted.
func Run(bRs []io.Reader, w io.Writer, opts Options) (int, error) {
	if len(bRs) < 2 {
		return 0, fmt.Errorf("multiinter requires at least 2 input files, got %d", len(bRs))
	}
	if len(opts.Filenames) != len(bRs) {
		return 0, fmt.Errorf("Filenames count %d does not match input count %d",
			len(opts.Filenames), len(bRs))
	}
	if len(opts.Names) > 0 && len(opts.Names) != len(bRs) {
		return 0, fmt.Errorf("Names count %d does not match input count %d",
			len(opts.Names), len(bRs))
	}
	if opts.Empty && opts.ChromSizes == nil {
		return 0, fmt.Errorf("-empty requires a genome (-g) for chromosome sizes")
	}
	filler := opts.Filler
	if filler == "" {
		filler = "0"
	}

	// Load each file's intervals, merging within-file overlaps.
	perFile := make([]map[string][][2]int, len(bRs))
	chromSet := map[string]struct{}{}
	for i, r := range bRs {
		recs, err := readAndMerge(r)
		if err != nil {
			return 0, fmt.Errorf("reading input %d (%s): %w", i+1, opts.Filenames[i], err)
		}
		perFile[i] = recs
		for chrom := range recs {
			chromSet[chrom] = struct{}{}
		}
	}
	chroms := make([]string, 0, len(chromSet))
	for chrom := range chromSet {
		chroms = append(chroms, chrom)
	}
	sort.Strings(chroms)

	bw := bufio.NewWriter(w)
	defer bw.Flush()

	if opts.Header {
		writeHeader(bw, opts)
	}

	rows := 0
	for _, chrom := range chroms {
		n, err := emitChrom(bw, chrom, perFile, opts, filler)
		if err != nil {
			return rows, err
		}
		rows += n
	}
	return rows, nil
}

// segment is one emitted region. Active is the indicator vector (len N).
type segment struct {
	start, end int
	active     []bool
	num        int
}

// emitChrom walks one chromosome's events, accumulates segments, then
// writes them according to opts.Cluster / opts.Empty.
func emitChrom(bw *bufio.Writer, chrom string, perFile []map[string][][2]int,
	opts Options, filler string) (int, error) {
	n := len(perFile)
	type event struct {
		coord, src int
		kind       int // +1 start, -1 end
	}
	events := []event{}
	for i, m := range perFile {
		for _, iv := range m[chrom] {
			events = append(events, event{coord: iv[0], src: i, kind: +1})
			events = append(events, event{coord: iv[1], src: i, kind: -1})
		}
	}
	if len(events) == 0 {
		// No data on this chrom; if -empty + size known, emit one gap row.
		if opts.Empty {
			if size, ok := opts.ChromSizes[chrom]; ok && size > 0 {
				active := make([]bool, n)
				return emitSegments(bw, chrom, []segment{{start: 0, end: size, active: active}},
					opts, filler)
			}
		}
		return 0, nil
	}
	// Sort by coord; END before START at the same coord so that a record
	// ending exactly where another starts does not generate a phantom
	// overlap.
	sort.SliceStable(events, func(i, j int) bool {
		if events[i].coord != events[j].coord {
			return events[i].coord < events[j].coord
		}
		return events[i].kind < events[j].kind // -1 (end) before +1 (start)
	})

	active := make([]bool, n)
	depth := 0
	var segs []segment
	prev := events[0].coord
	// Emit leading gap if -empty.
	if opts.Empty {
		if prev > 0 {
			emptyActive := make([]bool, n)
			segs = append(segs, segment{start: 0, end: prev, active: emptyActive})
		}
	}
	for _, e := range events {
		if e.coord > prev && depth > 0 {
			// Snapshot active set as a segment.
			snap := make([]bool, n)
			copy(snap, active)
			segs = append(segs, segment{start: prev, end: e.coord, active: snap, num: depth})
		} else if e.coord > prev && depth == 0 && opts.Empty {
			snap := make([]bool, n)
			segs = append(segs, segment{start: prev, end: e.coord, active: snap})
		}
		if e.kind == +1 {
			active[e.src] = true
			depth++
		} else {
			active[e.src] = false
			depth--
		}
		prev = e.coord
	}
	// Trailing gap if -empty.
	if opts.Empty {
		if size, ok := opts.ChromSizes[chrom]; ok && size > prev {
			emptyActive := make([]bool, n)
			segs = append(segs, segment{start: prev, end: size, active: emptyActive})
		}
	}

	if opts.Cluster {
		segs = clusterSegments(segs)
	}
	return emitSegments(bw, chrom, segs, opts, filler)
}

// clusterSegments collapses adjacent segments that have the same active
// set into a single segment (mirrors `multiinter -cluster`).
func clusterSegments(in []segment) []segment {
	if len(in) == 0 {
		return in
	}
	out := []segment{in[0]}
	for i := 1; i < len(in); i++ {
		cur := in[i]
		prev := &out[len(out)-1]
		if prev.end == cur.start && sameActive(prev.active, cur.active) {
			prev.end = cur.end
			continue
		}
		out = append(out, cur)
	}
	return out
}

// sameActive reports whether two active-set vectors are identical.
func sameActive(a, b []bool) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// emitSegments writes one row per segment using the configured filler.
func emitSegments(bw *bufio.Writer, chrom string, segs []segment, opts Options,
	filler string) (int, error) {
	n := 0
	for _, s := range segs {
		if s.num == 0 && !opts.Empty {
			continue
		}
		// Build the list field.
		var list string
		if s.num == 0 {
			list = "none"
		} else {
			parts := make([]string, 0, s.num)
			for i, on := range s.active {
				if !on {
					continue
				}
				if len(opts.Names) > 0 {
					parts = append(parts, opts.Names[i])
				} else {
					parts = append(parts, strconv.Itoa(i+1))
				}
			}
			list = strings.Join(parts, ",")
		}
		// chrom start end num list <per-file>
		if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%d\t%s",
			chrom, s.start, s.end, s.num, list); err != nil {
			return n, err
		}
		for _, on := range s.active {
			if on {
				if _, err := bw.WriteString("\t1"); err != nil {
					return n, err
				}
			} else {
				if _, err := bw.WriteString("\t" + filler); err != nil {
					return n, err
				}
			}
		}
		if err := bw.WriteByte('\n'); err != nil {
			return n, err
		}
		n++
	}
	return n, nil
}

// writeHeader emits the upstream-style header line. With Names empty,
// upstream uses the raw filenames as supplied; we mirror that.
func writeHeader(bw *bufio.Writer, opts Options) {
	bw.WriteString("chrom\tstart\tend\tnum\tlist")
	if len(opts.Names) > 0 {
		for _, n := range opts.Names {
			bw.WriteByte('\t')
			bw.WriteString(n)
		}
	} else {
		for _, fn := range opts.Filenames {
			bw.WriteByte('\t')
			bw.WriteString(fn)
		}
	}
	bw.WriteByte('\n')
}

// interval is one parsed [start,end) span on a chromosome, format-agnostic.
type interval struct {
	chrom      string
	start, end int
}

// inputFormat tags the autodetected format of an input file.
type inputFormat int

const (
	formatUnknown inputFormat = iota
	formatBED
	formatVCF
	formatGFF
)

// isInteger reports whether s is a base-10 integer, matching upstream's
// isInteger (which accepts an optional leading sign).
func isInteger(s string) bool {
	if s == "" {
		return false
	}
	_, err := strconv.Atoi(s)
	return err == nil
}

// detectFormat classifies a tokenized data line as BED, VCF, or GFF using
// upstream BedFile::parseLine's precedence: BED first (cols 2,3 integer),
// then VCF (col 2 integer and ≥8 cols), then GFF (8 or 9 cols with cols 4,5
// integer). Returns formatUnknown when none match.
func detectFormat(fields []string) inputFormat {
	n := len(fields)
	if n < 3 {
		return formatUnknown
	}
	if isInteger(fields[1]) && isInteger(fields[2]) {
		return formatBED
	}
	if isInteger(fields[1]) && n >= 8 {
		return formatVCF
	}
	if (n == 8 || n == 9) && isInteger(fields[3]) && isInteger(fields[4]) {
		return formatGFF
	}
	return formatUnknown
}

// parseTyped converts a tokenized data line into a 0-based half-open
// interval according to the resolved format. The coordinate conventions
// mirror upstream's parseVcfLine / parseGffLine / parseBedLine exactly:
//
//   - BED: start = col2, end = col3 (already 0-based half-open).
//   - VCF: start = POS-1, end = start + len(REF) (col4, the affected
//     reference allele); a 1-based, REF-length span.
//   - GFF: start = col4-1, end = col5 (1-based inclusive → 0-based
//     half-open).
func parseTyped(fields []string, format inputFormat) (interval, error) {
	switch format {
	case formatBED:
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return interval{}, fmt.Errorf("invalid BED start %q: %w", fields[1], err)
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return interval{}, fmt.Errorf("invalid BED end %q: %w", fields[2], err)
		}
		return interval{chrom: fields[0], start: start, end: end}, nil
	case formatVCF:
		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			return interval{}, fmt.Errorf("invalid VCF POS %q: %w", fields[1], err)
		}
		start := pos - 1
		end := start + len(fields[3])
		if start < 0 || end < start {
			return interval{}, fmt.Errorf("malformed VCF entry: start=%d end=%d", start, end)
		}
		return interval{chrom: fields[0], start: start, end: end}, nil
	case formatGFF:
		gstart, err := strconv.Atoi(fields[3])
		if err != nil {
			return interval{}, fmt.Errorf("invalid GFF start %q: %w", fields[3], err)
		}
		gend, err := strconv.Atoi(fields[4])
		if err != nil {
			return interval{}, fmt.Errorf("invalid GFF end %q: %w", fields[4], err)
		}
		start := gstart - 1
		end := gend
		if start < 0 || end < start {
			return interval{}, fmt.Errorf("malformed GFF entry: start=%d end=%d", start, end)
		}
		return interval{chrom: fields[0], start: start, end: end}, nil
	default:
		return interval{}, fmt.Errorf("unknown input format")
	}
}

// readIntervals reads every data line from r, autodetecting BED/VCF/GFF
// from the first non-header line and applying that format to the remainder
// of the file (matching upstream, which locks the file type on the first
// data record). Header lines (`#`, `track`, `browser`) and blank lines are
// skipped. `##fileformat=VCF` in the header forces VCF detection up front,
// as upstream's GetHeader does.
func readIntervals(r io.Reader) ([]interval, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var out []interval
	format := formatUnknown
	headerForcedVCF := false
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r")
		trimmed := strings.TrimSpace(line)
		if trimmed == "" {
			continue
		}
		if strings.HasPrefix(trimmed, "#") || strings.HasPrefix(trimmed, "track") ||
			strings.HasPrefix(trimmed, "browser") {
			if strings.HasPrefix(line, "##fileformat=VCF") {
				headerForcedVCF = true
			}
			continue
		}
		fields := strings.Split(line, "\t")
		if format == formatUnknown {
			if headerForcedVCF {
				format = formatVCF
			} else {
				format = detectFormat(fields)
			}
			if format == formatUnknown {
				return nil, fmt.Errorf("unexpected file format: please use tab-delimited BED, GFF, or VCF (line: %q)", line)
			}
		}
		iv, err := parseTyped(fields, format)
		if err != nil {
			return nil, err
		}
		out = append(out, iv)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// readAndMerge reads every record from r, groups by chrom, and merges
// records that overlap or abut within a single file. Returns a map
// chrom -> sorted, merged [start,end) intervals.
//
// The input format (BED, VCF, or GFF) is autodetected from the first
// non-header data line, mirroring upstream's BedFile::parseLine: a line
// whose 2nd and 3rd columns are integers is BED; a line whose 2nd column
// is an integer with ≥8 columns is VCF; a line with exactly 8 or 9 columns
// whose 4th and 5th columns are integers is GFF. Coordinates are converted
// to BED-style 0-based half-open spans exactly as upstream does.
func readAndMerge(r io.Reader) (map[string][][2]int, error) {
	intervals, err := readIntervals(r)
	if err != nil {
		return nil, err
	}
	byChrom := map[string][][2]int{}
	for _, iv := range intervals {
		byChrom[iv.chrom] = append(byChrom[iv.chrom], [2]int{iv.start, iv.end})
	}
	for chrom, ivs := range byChrom {
		sort.Slice(ivs, func(i, j int) bool {
			if ivs[i][0] != ivs[j][0] {
				return ivs[i][0] < ivs[j][0]
			}
			return ivs[i][1] < ivs[j][1]
		})
		merged := ivs[:0]
		for _, iv := range ivs {
			if len(merged) > 0 && iv[0] <= merged[len(merged)-1][1] {
				if iv[1] > merged[len(merged)-1][1] {
					merged[len(merged)-1][1] = iv[1]
				}
				continue
			}
			merged = append(merged, iv)
		}
		byChrom[chrom] = merged
	}
	return byChrom, nil
}

// ReadGenomeSizes parses a chrom-sizes file (tab-separated `chrom\tsize`
// lines, comments starting with '#'). Returned map matches Options.ChromSizes.
func ReadGenomeSizes(r io.Reader) (map[string]int, error) {
	out := map[string]int{}
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 1024*1024)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, fmt.Errorf("genome file line lacks 2 columns: %q", line)
		}
		n, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("invalid size %q for %q: %v", fields[1], fields[0], err)
		}
		out[fields[0]] = n
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// DefaultNames returns the upstream default label per path: the file's
// basename with the last extension stripped (matches `stl_basename`).
func DefaultNames(paths []string) []string {
	out := make([]string, len(paths))
	for i, p := range paths {
		base := p
		if idx := strings.LastIndex(p, "/"); idx >= 0 {
			base = p[idx+1:]
		}
		if dot := strings.LastIndex(base, "."); dot > 0 {
			base = base[:dot]
		}
		out[i] = base
	}
	return out
}
