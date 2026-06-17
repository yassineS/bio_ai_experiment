// Package bedannotate implements `bedtools annotate`: for each interval in
// a primary BED file (`-i`), it annotates with overlap statistics drawn
// from N additional BED files (`-files`). Three output modes mirror the
// upstream flags:
//
//   - default        — emit the fraction of A covered by each B (as %f)
//   - -counts        — emit the count of overlapping records per B
//   - -both          — emit `<count>\t<fraction>` per B
//
// Strand filters (-s same-strand, -S opposite-strand) restrict which
// overlaps count. The output preserves A's original columns and appends
// the per-B columns; an optional column-header line ("#…") is emitted
// when names are supplied (either explicitly via -names or implicitly
// from -files basenames).
//
// Internally we read each B file into a per-chromosome interval tree
// (`pkg/htsgo/bed.IntervalTree`) and stream A line by line, so
// the working set is O(sum(|B_i|)) rather than O(|A|·|B_i|).
package bedannotate

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// Mode controls which columns are appended per B file.
type Mode int

const (
	// ModeFraction emits a single float column per B (default).
	ModeFraction Mode = iota
	// ModeCounts emits a single int column per B (-counts).
	ModeCounts
	// ModeBoth emits `<count>\t<fraction>` per B (-both).
	ModeBoth
)

// Options configures Run.
type Options struct {
	// Mode selects the output shape. Default = ModeFraction.
	Mode Mode
	// Names holds the per-B header labels supplied via -names. A header line
	// is emitted ONLY when this is non-empty, matching upstream's
	// `if (_annoTitles.size() > 0) PrintHeader()`. The CLI must leave this nil
	// when -names is absent (file basenames do NOT trigger a header).
	Names []string
	// BedType is the column count of the main (-i) file, used to pad the
	// header's leading hash (upstream prints bedType-1 tabs after '#'). When
	// zero it is inferred from the first A record's column count.
	BedType int
	// SameStrand: -s, require A.Strand == B.Strand.
	SameStrand bool
	// OppositeStrand: -S, require A.Strand != B.Strand.
	OppositeStrand bool
}

// ucscBin replicates bedtools' getBin (src/utils/bedFile/bedFile.h): it maps a
// [start, end) interval to a UCSC binning-scheme bin. The annotate output is
// grouped per chromosome then ordered by ascending bin (and insertion order
// within a bin), so reproducing this function is required for byte-for-byte
// record ordering.
func ucscBin(start, end int) int {
	const (
		binFirstShift = 14
		binNextShift  = 3
		binLevels     = 8
	)
	binOffsetsExtended := [binLevels]int{
		262144 + 32678 + 4096 + 512 + 64 + 8 + 1,
		32678 + 4096 + 512 + 64 + 8 + 1,
		4096 + 512 + 64 + 8 + 1,
		512 + 64 + 8 + 1,
		64 + 8 + 1,
		8 + 1,
		1,
		0,
	}
	end--
	s := start >> binFirstShift
	e := end >> binFirstShift
	for i := 0; i < binLevels; i++ {
		if s == e {
			return binOffsetsExtended[i] + s
		}
		s >>= binNextShift
		e >>= binNextShift
	}
	return 0
}

// Run reads A from aR and the N B files from bRs in order, indexes each B
// into per-chromosome interval trees, then streams A and emits one row per
// A record with the per-B columns appended. Returns the number of A
// records processed.
func Run(aR io.Reader, bRs []io.Reader, out io.Writer, opts Options) (int, error) {
	if opts.SameStrand && opts.OppositeStrand {
		return 0, fmt.Errorf("cannot combine -s and -S")
	}
	// Build a tree per B file.
	trees := make([]map[string]*bed.IntervalTree, len(bRs))
	for i, br := range bRs {
		t, err := indexB(br)
		if err != nil {
			return 0, fmt.Errorf("file %d: %w", i+1, err)
		}
		trees[i] = t
	}

	// Read every A record, preserving its raw columns and recording its UCSC
	// bin. Upstream loads the main file into a map<chrom, map<bin, vector>> and
	// reports in that order: chromosome (lexicographic), then ascending bin,
	// then insertion order within a bin. We reproduce that exact ordering.
	type aRecord struct {
		rec    *bed.Record
		fields []string
		bin    int
		seq    int // input order, used as the within-bin tiebreaker
	}
	byChromBin := map[string][]*aRecord{}
	var chroms []string
	bedType := 0

	sc := bufio.NewScanner(aR)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	lineNo := 0
	seq := 0
	for sc.Scan() {
		lineNo++
		raw := sc.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(raw, "\t")
		if len(fields) < 3 {
			return 0, fmt.Errorf("line %d: BED record needs >=3 columns: %q", lineNo, raw)
		}
		rec, err := parseRecord(fields)
		if err != nil {
			return 0, fmt.Errorf("line %d: %w", lineNo, err)
		}
		if bedType == 0 {
			bedType = len(fields)
		}
		if _, seen := byChromBin[rec.Chrom]; !seen {
			chroms = append(chroms, rec.Chrom)
		}
		byChromBin[rec.Chrom] = append(byChromBin[rec.Chrom], &aRecord{
			rec: rec, fields: fields, bin: ucscBin(rec.ChromStart, rec.ChromEnd), seq: seq,
		})
		seq++
	}
	if err := sc.Err(); err != nil {
		return 0, err
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()

	// Header line: emitted ONLY when explicit -names labels were given.
	if len(opts.Names) > 0 {
		ht := opts.BedType
		if ht == 0 {
			ht = bedType
		}
		if ht == 0 {
			ht = 1
		}
		writeHeader(bw, opts.Names, opts.Mode, ht)
	}

	// Iterate chromosomes lexicographically.
	sort.Strings(chroms)
	count := 0
	for _, chrom := range chroms {
		recs := byChromBin[chrom]
		// Stable order by (bin asc, input order). Upstream's std::map orders by
		// bin; within a bin the vector keeps insertion (input) order.
		sort.SliceStable(recs, func(i, j int) bool {
			if recs[i].bin != recs[j].bin {
				return recs[i].bin < recs[j].bin
			}
			return recs[i].seq < recs[j].seq
		})
		for _, ar := range recs {
			if _, err := bw.WriteString(strings.Join(ar.fields, "\t")); err != nil {
				return count, err
			}
			for _, t := range trees {
				matches := selectOverlapping(ar.rec, t[ar.rec.Chrom], opts)
				cnt := len(matches)
				frac := coveredFraction(ar.rec, matches)
				switch opts.Mode {
				case ModeCounts:
					fmt.Fprintf(bw, "\t%d", cnt)
				case ModeBoth:
					fmt.Fprintf(bw, "\t%d\t%f", cnt, frac)
				default:
					fmt.Fprintf(bw, "\t%f", frac)
				}
			}
			if err := bw.WriteByte('\n'); err != nil {
				return count, err
			}
			count++
		}
	}
	return count, nil
}

// writeHeader emits the "#<TAB>...<TAB>name1<TAB>name2..." header line,
// matching upstream's PrintHeader: a '#' followed by bedType-1 empty tabs (so
// the first label aligns under the first appended column), then one tab plus
// each label. With -both, each label is split into `_cnt`/`_pct`.
func writeHeader(w *bufio.Writer, names []string, mode Mode, bedType int) {
	w.WriteByte('#')
	for i := 1; i < bedType; i++ {
		w.WriteByte('\t')
	}
	if mode == ModeBoth {
		for _, n := range names {
			fmt.Fprintf(w, "\t%s_cnt\t%s_pct", n, n)
		}
	} else {
		for _, n := range names {
			fmt.Fprintf(w, "\t%s", n)
		}
	}
	w.WriteByte('\n')
}

// parseRecord parses the minimum subset of a BED line we need for overlap
// + strand filtering. Extra columns are preserved by the caller as raw
// fields.
func parseRecord(fields []string) (*bed.Record, error) {
	start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
	if err != nil {
		return nil, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
	}
	end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
	if err != nil {
		return nil, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
	}
	r := &bed.Record{
		Chrom:      fields[0],
		ChromStart: start,
		ChromEnd:   end,
	}
	if len(fields) >= 6 {
		r.Strand = fields[5]
	}
	return r, nil
}

// indexB reads a B file fully into memory and returns a per-chrom tree.
func indexB(r io.Reader) (map[string]*bed.IntervalTree, error) {
	rd := bed.NewReader(r)
	all, err := rd.ReadAll()
	if err != nil {
		return nil, err
	}
	byChrom := map[string][]*bed.Record{}
	for _, x := range all {
		byChrom[x.Chrom] = append(byChrom[x.Chrom], x)
	}
	out := make(map[string]*bed.IntervalTree, len(byChrom))
	for chrom, recs := range byChrom {
		sort.SliceStable(recs, func(i, j int) bool {
			if recs[i].ChromStart != recs[j].ChromStart {
				return recs[i].ChromStart < recs[j].ChromStart
			}
			return recs[i].ChromEnd < recs[j].ChromEnd
		})
		out[chrom] = bed.NewIntervalTree(recs)
	}
	return out, nil
}

// selectOverlapping returns the B records overlapping a (after applying
// the strand filters in opts).
func selectOverlapping(a *bed.Record, t *bed.IntervalTree, opts Options) []*bed.Record {
	if t == nil {
		return nil
	}
	cand := t.Query(a)
	if len(cand) == 0 {
		return nil
	}
	out := cand[:0:0]
	for _, b := range cand {
		if !strandOK(a, b, opts) {
			continue
		}
		out = append(out, b)
	}
	return out
}

// strandOK applies the -s / -S filters. Upstream annotate compares the strand
// strings directly (strands_are_same = a.strand == b.strand) with no special
// handling for missing strands, so two BED3 records (both with the default
// strand) compare as same-strand. We replicate that raw comparison: normalise
// the empty default to "." (the value upstream's BED parser uses) so BED3 vs
// BED3 counts as same-strand under -s.
func strandOK(a, b *bed.Record, opts Options) bool {
	if !opts.SameStrand && !opts.OppositeStrand {
		return true
	}
	sa, sb := a.Strand, b.Strand
	if sa == "" {
		sa = "."
	}
	if sb == "" {
		sb = "."
	}
	same := sa == sb
	if opts.SameStrand {
		return same
	}
	return !same // OppositeStrand
}

// coveredFraction returns the fraction of A's length covered by at least
// one of the matches (depth >= 1). Used for default and -both modes.
func coveredFraction(a *bed.Record, matches []*bed.Record) float64 {
	lenA := a.ChromEnd - a.ChromStart
	if lenA <= 0 {
		return 0
	}
	if len(matches) == 0 {
		return 0
	}
	covered := make([]bool, lenA)
	for _, b := range matches {
		start := b.ChromStart - a.ChromStart
		end := b.ChromEnd - a.ChromStart
		if start < 0 {
			start = 0
		}
		if end > lenA {
			end = lenA
		}
		for i := start; i < end; i++ {
			covered[i] = true
		}
	}
	n := 0
	for _, c := range covered {
		if c {
			n++
		}
	}
	return float64(n) / float64(lenA)
}
