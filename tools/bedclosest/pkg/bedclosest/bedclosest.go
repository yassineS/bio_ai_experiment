// Package bedclosest finds, for each interval in A, the closest interval in B
// (mirrors `bedtools closest`).
//
// Both inputs MUST be sorted on (chrom, start). For each A interval the
// closest B on the same chromosome (lowest signed distance; 0 if they overlap)
// is reported. The output line is A's columns + B's columns + the signed
// distance. For tied distances, one row per tied B is emitted by default
// (in B's input order); this is configurable via Options.TieBreak. The sign
// of the distance is controlled by Options.DistanceMode.
package bedclosest

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"
)

// DistanceMode selects how the distance between A and B is computed and
// reported.
type DistanceMode int

const (
	// DistanceRef computes the sign on the reference: downstream is positive
	// (B.start > A.end -> positive distance, B.end < A.start -> negative).
	DistanceRef DistanceMode = iota
	// DistanceA computes the sign relative to A's strand (BED6 col 6). On a
	// '-' strand A, what was "downstream on the reference" becomes "upstream"
	// and the sign is flipped.
	DistanceA
	// DistanceB computes the sign relative to B's strand.
	DistanceB
	// DistanceAbsolute reports the unsigned distance only. This matches
	// upstream `bedtools closest -d`, which omits the sign; `-D <mode>` is
	// the upstream flag for signed distance.
	DistanceAbsolute
)

// TieBreak controls how ties among multiple equally-close B intervals are
// handled.
type TieBreak int

const (
	// TieAll emits one output row per tied B.
	TieAll TieBreak = iota
	// TieFirst emits only the first tied B (lowest position in B's input).
	TieFirst
	// TieLast emits only the last tied B.
	TieLast
)

// Options configures Closest.
type Options struct {
	// PrintDistance controls whether the trailing signed-distance column is
	// emitted. Default true (bedclosest deviates from upstream bedtools here,
	// where the distance column is opt-in).
	PrintDistance bool
	// DistanceMode controls the sign convention for the distance column.
	DistanceMode DistanceMode
	// RequireOverlap, when set, only emits rows where A and B overlap. All
	// non-overlapping B intervals are treated as infinity (skipped).
	// NOTE: this is distinct from upstream `-N`; see RequireDifferentNames.
	RequireOverlap bool
	// RequireDifferentNames (`-N`) forces the closest B to have a different
	// name field (BED column 4) than A. B's with the same name are ignored.
	RequireDifferentNames bool
	// SameStrand (`-s`) restricts the search to B's whose strand matches A's
	// strand (BED column 6). If A has no strand, all B's are eligible.
	SameStrand bool
	// OppositeStrand (`-S`) restricts the search to B's whose strand is the
	// opposite of A's strand. Mutually exclusive with SameStrand at upstream;
	// here the caller is expected not to set both.
	OppositeStrand bool
	// TieBreak controls how ties (multiple equally-close B intervals) are
	// resolved. Default TieAll.
	TieBreak TieBreak
}

// Row is the parsed representation of a BED line. The original fields are
// preserved so the full input column count round-trips to the output.
type Row struct {
	Fields []string
	Chrom  string
	Start  int
	End    int
	Strand string // fields[5] if len(fields) >= 6, else "+"
}

// MissingRow is a sentinel "no closest B" row used when A's chromosome doesn't
// exist in B. It contains "." for chrom and -1 for the coordinates, matching
// bedtools convention. The Fields slice is padded by Closest to match B's
// observed column width using upstream's Bed{3,4,5,6,12}::printNull format.
var MissingRow = &Row{Fields: []string{".", "-1", "-1"}, Chrom: ".", Start: -1, End: -1, Strand: "+"}

// missingRowForCols returns a MissingRow padded to bCols using upstream's
// Bed{3,4,5,6,12}::printNull format.
func missingRowForCols(bCols int) *Row {
	var fields []string
	switch {
	case bCols >= 12:
		fields = []string{".", "-1", "-1", ".", "-1", ".", ".", ".", ".", ".", ".", "."}
	case bCols >= 6:
		fields = []string{".", "-1", "-1", ".", "-1", "."}
	case bCols >= 5:
		fields = []string{".", "-1", "-1", ".", "-1"}
	case bCols >= 4:
		fields = []string{".", "-1", "-1", "."}
	default:
		fields = []string{".", "-1", "-1"}
	}
	return &Row{Fields: fields, Chrom: ".", Start: -1, End: -1, Strand: "+"}
}

// ReadAll parses every BED record from r into Row values. Lines that begin
// with '#', "track", or "browser" are skipped, as are blank lines.
func ReadAll(r io.Reader) ([]*Row, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var rows []*Row
	lineNum := 0
	for sc.Scan() {
		lineNum++
		raw := sc.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(raw, "\t")
		if len(fields) < 3 {
			return nil, fmt.Errorf("line %d: BED record must have at least 3 fields", lineNum)
		}
		start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
		if err != nil {
			return nil, fmt.Errorf("line %d: invalid start %q: %v", lineNum, fields[1], err)
		}
		end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
		if err != nil {
			return nil, fmt.Errorf("line %d: invalid end %q: %v", lineNum, fields[2], err)
		}
		if end < start {
			return nil, fmt.Errorf("line %d: end < start (%d < %d)", lineNum, end, start)
		}
		row := &Row{Fields: fields, Chrom: fields[0], Start: start, End: end, Strand: "+"}
		if len(fields) >= 6 {
			row.Strand = fields[5]
		}
		rows = append(rows, row)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return rows, nil
}

// CheckSorted returns an error naming the offending line if rows are not
// sorted on (chrom, start). Equal chromosomes must have non-decreasing starts.
func CheckSorted(rows []*Row, label string) error {
	for i := 1; i < len(rows); i++ {
		prev, cur := rows[i-1], rows[i]
		if cur.Chrom == prev.Chrom {
			if cur.Start < prev.Start {
				return fmt.Errorf("%s is not sorted on (chrom, start) at record %d: %s\t%d came after %s\t%d",
					label, i+1, cur.Chrom, cur.Start, prev.Chrom, prev.Start)
			}
		}
	}
	return nil
}

// Closest reads BED records from readerA and readerB, finds the closest B for
// each A, and writes the result to writer. Both inputs MUST be sorted on
// (chrom, start); otherwise it returns an error. Returns the number of output
// rows.
func Closest(readerA, readerB io.Reader, writer io.Writer, opts Options) (int, error) {
	aRows, err := ReadAll(readerA)
	if err != nil {
		return 0, fmt.Errorf("error reading A: %w", err)
	}
	bRows, err := ReadAll(readerB)
	if err != nil {
		return 0, fmt.Errorf("error reading B: %w", err)
	}
	if err := CheckSorted(aRows, "A"); err != nil {
		return 0, err
	}
	if err := CheckSorted(bRows, "B"); err != nil {
		return 0, err
	}

	// Bucket B by chromosome; within each chromosome B is already sorted by
	// start (CheckSorted guarantees that, and we keep original ordering). We
	// also precompute the running max .End for each chrom slice so closestFor
	// can prune its left walk safely even when a long interval upstream might
	// still reach back to A.
	type chromB struct {
		rows       []*Row
		maxEndPref []int // maxEndPref[i] = max(rows[0..i].End)
	}
	bByChrom := make(map[string]*chromB, 16)
	bMaxCols := 3
	for _, b := range bRows {
		c := bByChrom[b.Chrom]
		if c == nil {
			c = &chromB{}
			bByChrom[b.Chrom] = c
		}
		c.rows = append(c.rows, b)
		if n := len(b.Fields); n > bMaxCols {
			bMaxCols = n
		}
	}
	for _, c := range bByChrom {
		c.maxEndPref = make([]int, len(c.rows))
		max := 0
		for i, r := range c.rows {
			if r.End > max {
				max = r.End
			}
			c.maxEndPref[i] = max
		}
	}
	missing := missingRowForCols(bMaxCols)

	bw := bufio.NewWriter(writer)
	defer bw.Flush()

	count := 0
	for _, a := range aRows {
		var bs []*Row
		var maxEnd []int
		if c := bByChrom[a.Chrom]; c != nil {
			bs = c.rows
			maxEnd = c.maxEndPref
		}
		hits := closestFor(a, bs, maxEnd, opts, missing)
		for _, h := range hits {
			if err := writeRow(bw, a, h.b, h.dist, opts); err != nil {
				return count, fmt.Errorf("error writing output: %w", err)
			}
			count++
		}
	}
	return count, nil
}

// hit holds a candidate B together with its signed distance to A.
type hit struct {
	b    *Row
	dist int64
}

// ClosestMulti is the multiple-database equivalent of Closest. For each A row
// it emits one output line per database, in the order the databases were
// supplied. The line is:
//
//	<A's columns>  <label>  <closest-B-from-that-db's columns>  [<distance>]
//
// where <label> is dbLabels[i] (use the empty string to omit the label column;
// pass strconv.Itoa(i+1) for upstream's default "1/2/3" labelling, the basename
// for -filenames, or the user-supplied -names tokens). When a database has no
// B on A's chromosome it is reported with the per-db null placeholder
// (Bed3/4/5/6/12 width detected per database).
func ClosestMulti(readerA io.Reader, dbReaders []io.Reader, writer io.Writer, opts Options, dbLabels []string) (int, error) {
	if len(dbReaders) == 0 {
		return 0, fmt.Errorf("ClosestMulti: no databases supplied")
	}
	if len(dbLabels) != 0 && len(dbLabels) != len(dbReaders) {
		return 0, fmt.Errorf("ClosestMulti: dbLabels (%d) must match dbReaders (%d) or be empty",
			len(dbLabels), len(dbReaders))
	}

	aRows, err := ReadAll(readerA)
	if err != nil {
		return 0, fmt.Errorf("error reading A: %w", err)
	}
	if err := CheckSorted(aRows, "A"); err != nil {
		return 0, err
	}

	type dbIndex struct {
		rows       map[string][]*Row
		maxEndPref map[string][]int
		missing    *Row
	}
	dbs := make([]*dbIndex, len(dbReaders))
	for di, r := range dbReaders {
		bRows, err := ReadAll(r)
		if err != nil {
			return 0, fmt.Errorf("error reading db %d: %w", di+1, err)
		}
		if err := CheckSorted(bRows, fmt.Sprintf("db %d", di+1)); err != nil {
			return 0, err
		}
		bByChrom := make(map[string][]*Row, 16)
		maxCols := 3
		for _, b := range bRows {
			bByChrom[b.Chrom] = append(bByChrom[b.Chrom], b)
			if n := len(b.Fields); n > maxCols {
				maxCols = n
			}
		}
		maxEndByChrom := make(map[string][]int, len(bByChrom))
		for chrom, rs := range bByChrom {
			pref := make([]int, len(rs))
			m := 0
			for i, r := range rs {
				if r.End > m {
					m = r.End
				}
				pref[i] = m
			}
			maxEndByChrom[chrom] = pref
		}
		dbs[di] = &dbIndex{rows: bByChrom, maxEndPref: maxEndByChrom, missing: missingRowForCols(maxCols)}
	}

	bw := bufio.NewWriter(writer)
	defer bw.Flush()

	count := 0
	for _, a := range aRows {
		for di, db := range dbs {
			bs := db.rows[a.Chrom]
			maxEnd := db.maxEndPref[a.Chrom]
			hits := closestFor(a, bs, maxEnd, opts, db.missing)
			label := ""
			if len(dbLabels) > 0 {
				label = dbLabels[di]
			}
			for _, h := range hits {
				if err := writeRowLabeled(bw, a, label, h.b, h.dist, opts); err != nil {
					return count, fmt.Errorf("error writing output: %w", err)
				}
				count++
			}
		}
	}
	return count, nil
}

// closestFor returns the closest hits for A on its chromosome, taking
// Options.TieBreak into account.
//
// The implementation does a binary search to locate the first B with
// Start >= a.Start, then walks outward (left then right) collecting
// candidates. Once the gap on the corresponding side exceeds the best
// absolute distance seen so far, that side's walk can stop because
// monotonically increasing |B.Start - a.Start| (left) / B.Start - a.End
// (right) is then a strict lower bound on the gap for all remaining B's on
// that side. To handle the case where a B further left has a long End that
// could overlap A, we widen the left walk: any B that overlaps a.Start (i.e.
// B.End > a.Start) is also considered.
func closestFor(a *Row, bs []*Row, maxEndPref []int, opts Options, missing *Row) []hit {
	if len(bs) == 0 {
		if opts.RequireOverlap {
			return nil
		}
		return []hit{{b: missing, dist: -1}}
	}

	idx := sort.Search(len(bs), func(i int) bool { return bs[i].Start >= a.Start })

	bestAbs := int64(math.MaxInt64)
	type cand struct {
		idx    int
		signed int64
	}
	var cands []cand
	// nameOf returns the BED column-4 name, or "" if the record has < 4 fields.
	nameOf := func(r *Row) string {
		if len(r.Fields) >= 4 {
			return r.Fields[3]
		}
		return ""
	}
	consider := func(i int) {
		b := bs[i]
		// -N: skip B's whose name equals A's.
		if opts.RequireDifferentNames && nameOf(a) == nameOf(b) {
			return
		}
		// -s / -S: strand filters. Only applied when B has a strand column
		// (matches upstream which treats a missing strand as a mismatch).
		if opts.SameStrand {
			if a.Strand == "" || b.Strand == "" || a.Strand != b.Strand {
				return
			}
		}
		if opts.OppositeStrand {
			if a.Strand == "" || b.Strand == "" || a.Strand == b.Strand {
				return
			}
			if a.Strand != "+" && a.Strand != "-" {
				return
			}
		}
		signed := signedDistance(a, b, opts)
		if opts.RequireOverlap && signed != 0 {
			return
		}
		abs := signed
		if abs < 0 {
			abs = -abs
		}
		if abs < bestAbs {
			bestAbs = abs
			cands = cands[:0]
			cands = append(cands, cand{idx: i, signed: signed})
		} else if abs == bestAbs {
			cands = append(cands, cand{idx: i, signed: signed})
		}
	}

	// Walk left from idx-1. maxEndPref[i] is the largest End of bs[0..i]; the
	// minimum possible reference gap from any B in bs[0..i] is
	//   max(0, a.Start - maxEndPref[i])
	// so once that strictly exceeds bestAbs we can stop walking left. We use
	// '>' (not '>=') so that we still walk through ties at the same distance.
	for i := idx - 1; i >= 0; i-- {
		consider(i)
		if i > 0 {
			minPossibleGap := int64(a.Start - maxEndPref[i-1])
			if minPossibleGap < 0 {
				minPossibleGap = 0
			}
			if minPossibleGap > bestAbs {
				break
			}
		}
	}
	// Walk right from idx. bs[i].Start is monotonically non-decreasing; once
	// bs[i].Start - a.End > bestAbs we can stop. Again use '>' to keep ties.
	for i := idx; i < len(bs); i++ {
		consider(i)
		gap := int64(bs[i].Start - a.End)
		if gap < 0 {
			gap = 0
		}
		if gap > bestAbs {
			break
		}
	}

	if len(cands) == 0 {
		if opts.RequireOverlap {
			return nil
		}
		return []hit{{b: missing, dist: -1}}
	}

	// Sort candidates by their index in B (input order) for deterministic output.
	sort.Slice(cands, func(i, j int) bool { return cands[i].idx < cands[j].idx })

	switch opts.TieBreak {
	case TieFirst:
		first := cands[0]
		return []hit{{b: bs[first.idx], dist: first.signed}}
	case TieLast:
		last := cands[len(cands)-1]
		return []hit{{b: bs[last.idx], dist: last.signed}}
	default: // TieAll
		out := make([]hit, 0, len(cands))
		for _, c := range cands {
			out = append(out, hit{b: bs[c.idx], dist: c.signed})
		}
		return out
	}
}

// signedDistance returns the signed distance between A and B according to
// opts.DistanceMode. 0 means they overlap (truly share at least one base on a
// 0-based half-open interval); touching records (b.start == a.end or
// a.start == b.end) report distance 1, matching `bedtools closest -d`'s
// (b.start - a.end) + 1 formula. The sign indicates whether B is downstream
// (positive) or upstream (negative) of A under the chosen DistanceMode.
func signedDistance(a, b *Row, opts Options) int64 {
	// Overlap on a 0-based half-open interval requires a.Start < b.End AND
	// b.Start < a.End.
	if a.Start < b.End && b.Start < a.End {
		return 0
	}
	var refSigned int64
	if b.Start >= a.End {
		// B is downstream of A. Upstream uses (B.start - A.end) + 1, so a 0-bp
		// gap (touching intervals) is reported as 1.
		refSigned = int64(b.Start-a.End) + 1
	} else {
		// B is upstream of A.
		refSigned = -(int64(a.Start-b.End) + 1)
	}
	switch opts.DistanceMode {
	case DistanceA:
		if a.Strand == "-" {
			return -refSigned
		}
		return refSigned
	case DistanceB:
		if b.Strand == "-" {
			return -refSigned
		}
		return refSigned
	case DistanceAbsolute:
		if refSigned < 0 {
			return -refSigned
		}
		return refSigned
	default: // DistanceRef
		return refSigned
	}
}

// writeRowLabeled writes A's columns, an optional label column (between A
// and B), B's columns, and optionally the trailing distance column. The
// label is omitted entirely when it is empty.
func writeRowLabeled(bw *bufio.Writer, a *Row, label string, b *Row, dist int64, opts Options) error {
	if _, err := bw.WriteString(strings.Join(a.Fields, "\t")); err != nil {
		return err
	}
	if label != "" {
		if err := bw.WriteByte('\t'); err != nil {
			return err
		}
		if _, err := bw.WriteString(label); err != nil {
			return err
		}
	}
	if err := bw.WriteByte('\t'); err != nil {
		return err
	}
	if _, err := bw.WriteString(strings.Join(b.Fields, "\t")); err != nil {
		return err
	}
	if opts.PrintDistance {
		if err := bw.WriteByte('\t'); err != nil {
			return err
		}
		if _, err := bw.WriteString(strconv.FormatInt(dist, 10)); err != nil {
			return err
		}
	}
	return bw.WriteByte('\n')
}

// writeRow writes one output row: A's columns, then B's columns, then the
// signed distance.
func writeRow(bw *bufio.Writer, a, b *Row, dist int64, opts Options) error {
	if _, err := bw.WriteString(strings.Join(a.Fields, "\t")); err != nil {
		return err
	}
	if err := bw.WriteByte('\t'); err != nil {
		return err
	}
	if _, err := bw.WriteString(strings.Join(b.Fields, "\t")); err != nil {
		return err
	}
	if opts.PrintDistance {
		if err := bw.WriteByte('\t'); err != nil {
			return err
		}
		if _, err := bw.WriteString(strconv.FormatInt(dist, 10)); err != nil {
			return err
		}
	}
	return bw.WriteByte('\n')
}
