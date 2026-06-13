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
	// RequireOverlap (`-N`) only emits rows where A and B overlap. All
	// non-overlapping B intervals are treated as infinity (skipped).
	RequireOverlap bool
	// TieBreak controls how ties (multiple equally-close B intervals) are
	// resolved. Default TieAll.
	TieBreak TieBreak

	// IgnoreUpstream (`-iu`) drops B features that are upstream of A; the
	// closest B is selected only among downstream and overlapping features.
	// IgnoreDownstream (`-id`) is the symmetric flag. Upstream/downstream is
	// determined relative to the active stranded distance mode (`-D`), so both
	// require DistanceMode != DistanceRef-by-default — the CLI enforces that
	// they are only set together with an explicit `-D`.
	IgnoreUpstream   bool
	IgnoreDownstream bool

	// ForceUpstream (`-fu`) reports the closest upstream B in preference to any
	// equally-or-closer downstream/overlapping B; ForceDownstream (`-fd`) is the
	// symmetric flag. Like the ignore flags, these are meaningful only with an
	// explicit `-D`.
	ForceUpstream   bool
	ForceDownstream bool

	// SameStrand (`-s`) restricts the candidate B intervals to those on the
	// same strand as A; SameStrand and OppositeStrand are mutually exclusive.
	SameStrand bool
	// OppositeStrand (`-S`) restricts the candidate B intervals to those on
	// the opposite strand to A.
	OppositeStrand bool

	// DifferentNames (`-N`) requires the reported closest B to have a
	// different name (BED column 4) than A; a B sharing A's name is excluded
	// from candidate consideration. Matches upstream bedtools closest -N.
	DifferentNames bool
}

// Validate reports configuration errors that cannot be captured by the type
// system, mirroring upstream's mutually-exclusive flag checks.
func (o Options) Validate() error {
	if o.SameStrand && o.OppositeStrand {
		return fmt.Errorf("-s and -S are mutually exclusive")
	}
	return nil
}

// strandMatch reports whether B is an eligible candidate for A under the strand
// filters. With neither -s nor -S set every B is eligible. With -s, only B's on
// the same strand as A qualify; with -S, only B's on the opposite strand. A
// missing or unknown strand (empty or ".") on either side cannot be classified
// as same or opposite, so such a B is excluded, matching upstream bedtools.
func strandMatch(a, b *Row, opts Options) bool {
	if !opts.SameStrand && !opts.OppositeStrand {
		return true
	}
	if a.Strand == "" || a.Strand == "." || b.Strand == "" || b.Strand == "." {
		return false
	}
	if opts.SameStrand {
		return a.Strand == b.Strand
	}
	return a.Strand != b.Strand
}

// nameOf returns a row's BED name (column 4) or "" when absent.
func nameOf(r *Row) string {
	if len(r.Fields) >= 4 {
		return r.Fields[3]
	}
	return ""
}

// nameEligible reports whether B is an eligible candidate for A under the
// -N (different-names) filter. With -N off every B qualifies; with -N on a B
// sharing A's name (column 4) is excluded, matching upstream bedtools closest.
func nameEligible(a, b *Row, opts Options) bool {
	if !opts.DifferentNames {
		return true
	}
	return nameOf(a) != nameOf(b)
}

// streamDir classifies a non-overlapping B hit as upstream or downstream of A
// under the active distance mode, mirroring upstream CloseSweep::considerRecord.
type streamDir int

const (
	streamOverlap streamDir = iota
	streamUpstream
	streamDownstream
)

// classifyStream returns whether B is upstream, downstream, or overlapping A,
// mirroring the UPSTREAM/DOWNSTREAM assignment in CloseSweep::considerRecord.
// A B to the right of A (b.Start >= a.End) is normally DOWNSTREAM, but becomes
// UPSTREAM when (A_dist && A is reverse-strand) or (B_dist && B is
// forward-strand); a B to the left is the mirror image. Overlaps are neither.
func classifyStream(a, b *Row, opts Options) streamDir {
	if a.Start < b.End && b.Start < a.End {
		return streamOverlap
	}
	flip := (opts.DistanceMode == DistanceA && a.Strand == "-") ||
		(opts.DistanceMode == DistanceB && b.Strand == "+")
	if b.Start >= a.End {
		// B is to the right of A.
		if flip {
			return streamUpstream
		}
		return streamDownstream
	}
	// B is to the left of A.
	if flip {
		return streamDownstream
	}
	return streamUpstream
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
// bedtools convention.
var MissingRow = &Row{Fields: []string{".", "-1", "-1"}, Chrom: ".", Start: -1, End: -1, Strand: "+"}

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
	if err := opts.Validate(); err != nil {
		return 0, err
	}
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
	for _, b := range bRows {
		c := bByChrom[b.Chrom]
		if c == nil {
			c = &chromB{}
			bByChrom[b.Chrom] = c
		}
		c.rows = append(c.rows, b)
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
		hits := closestFor(a, bs, maxEnd, opts)
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
func closestFor(a *Row, bs []*Row, maxEndPref []int, opts Options) []hit {
	if len(bs) == 0 {
		if opts.RequireOverlap {
			return nil
		}
		return []hit{{b: MissingRow, dist: -1}}
	}

	// The directional flags (-iu/-id/-fu/-fd) change which candidates are
	// eligible and the selection priority, so they take a dedicated stream-aware
	// path rather than the plain closest-by-absolute-distance scan below.
	if opts.IgnoreUpstream || opts.IgnoreDownstream || opts.ForceUpstream || opts.ForceDownstream {
		return closestForDirectional(a, bs, opts)
	}

	idx := sort.Search(len(bs), func(i int) bool { return bs[i].Start >= a.Start })

	bestAbs := int64(math.MaxInt64)
	type cand struct {
		idx    int
		signed int64
	}
	var cands []cand
	consider := func(i int) {
		// Strand- and name-ineligible B's are skipped before they can influence
		// bestAbs or the candidate set, so the closest is chosen purely among the
		// eligible subset (upstream bug281 cache-purge semantics).
		if !strandMatch(a, bs[i], opts) || !nameEligible(a, bs[i], opts) {
			return
		}
		signed := signedDistance(a, bs[i], opts)
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
		return []hit{{b: MissingRow, dist: -1}}
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

// refGap returns the absolute reference gap between non-overlapping A and B,
// using upstream's (gap + 1) convention so touching intervals report 1. It
// assumes A and B do not overlap (the caller classifies overlaps separately).
func refGap(a, b *Row) int64 {
	if b.Start >= a.End {
		return int64(b.Start-a.End) + 1
	}
	return int64(a.Start-b.End) + 1
}

// dcand is a directional-path candidate: its index in B, signed distance,
// absolute distance, and stream classification relative to A.
type dcand struct {
	idx    int
	signed int64
	abs    int64
	stream streamDir
}

// closestForDirectional selects the closest B for A honouring the directional
// flags, mirroring CloseSweep::finalizeSelections for the single-closest
// (k=1) case this tool reports:
//
//   - Each B on A's chromosome is classified as overlap/upstream/downstream
//     (classifyStream), and -iu/-id drop the corresponding stream entirely.
//   - -fu (forceUpstream) consumes the closest upstream feature(s) first, even
//     ahead of overlaps; -fd does the same for downstream. Only if the forced
//     stream is empty does selection fall through to overlaps, then to the
//     closest of the remaining upstream/downstream features.
//   - Without a force flag, overlaps win (distance 0), then the closest of the
//     surviving upstream/downstream features, with ties resolved per TieBreak.
func closestForDirectional(a *Row, bs []*Row, opts Options) []hit {
	var overlaps, ups, downs []dcand
	for i, b := range bs {
		if a.Chrom != b.Chrom {
			continue
		}
		// Skip strand- and name-ineligible B's entirely so they never enter the
		// candidate pools (overlaps/ups/downs), matching the non-directional path.
		if !strandMatch(a, b, opts) || !nameEligible(a, b, opts) {
			continue
		}
		stream := classifyStream(a, b, opts)
		abs := refGap(a, b)
		// Mirror CloseSweep::finalizeSelections, where the reported distance sign
		// is derived from the stream classification: upstream features are
		// negative, downstream positive, overlaps zero. signedDistance reuses the
		// same classifyStream decision, so the directional and non-directional
		// paths agree on the sign (including the -D b mode).
		switch stream {
		case streamOverlap:
			overlaps = append(overlaps, dcand{i, 0, 0, streamOverlap})
		case streamUpstream:
			if !opts.IgnoreUpstream && !opts.RequireOverlap {
				ups = append(ups, dcand{i, -abs, abs, streamUpstream})
			}
		case streamDownstream:
			if !opts.IgnoreDownstream && !opts.RequireOverlap {
				downs = append(downs, dcand{i, abs, abs, streamDownstream})
			}
		}
	}

	// closestOf returns the minimum-abs-distance candidates from a stream group
	// (all ties at that distance), in B input order.
	closestOf := func(group []dcand) []dcand {
		if len(group) == 0 {
			return nil
		}
		best := int64(math.MaxInt64)
		for _, c := range group {
			if c.abs < best {
				best = c.abs
			}
		}
		var out []dcand
		for _, c := range group {
			if c.abs == best {
				out = append(out, c)
			}
		}
		return out
	}

	var chosen []dcand
	switch {
	case opts.ForceUpstream:
		if chosen = closestOf(ups); chosen == nil {
			if chosen = overlaps; len(chosen) == 0 {
				chosen = closestOf(downs)
			} else {
				chosen = closestOf(overlaps)
			}
		}
	case opts.ForceDownstream:
		if chosen = closestOf(downs); chosen == nil {
			if chosen = overlaps; len(chosen) == 0 {
				chosen = closestOf(ups)
			} else {
				chosen = closestOf(overlaps)
			}
		}
	default:
		// No force flag (-iu / -id only): overlaps win, then the closest of the
		// surviving upstream/downstream features.
		if len(overlaps) > 0 {
			chosen = closestOf(overlaps)
		} else {
			rest := append(append([]dcand(nil), ups...), downs...)
			chosen = closestOf(rest)
		}
	}

	if len(chosen) == 0 {
		if opts.RequireOverlap {
			return nil
		}
		return []hit{{b: MissingRow, dist: -1}}
	}

	sort.Slice(chosen, func(i, j int) bool { return chosen[i].idx < chosen[j].idx })
	switch opts.TieBreak {
	case TieFirst:
		return []hit{{b: bs[chosen[0].idx], dist: chosen[0].signed}}
	case TieLast:
		last := chosen[len(chosen)-1]
		return []hit{{b: bs[last.idx], dist: last.signed}}
	default:
		out := make([]hit, 0, len(chosen))
		for _, c := range chosen {
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
	// Magnitude of the gap: (b.start - a.end) + 1 to the right, or
	// (a.start - b.end) + 1 to the left, matching upstream's currDist.
	var magnitude int64
	if b.Start >= a.End {
		magnitude = int64(b.Start-a.End) + 1
	} else {
		magnitude = int64(a.Start-b.End) + 1
	}
	if opts.DistanceMode == DistanceAbsolute {
		return magnitude
	}
	// The sign is negative when the hit is classified UPSTREAM of A and
	// positive when DOWNSTREAM, exactly as upstream CloseSweep applies
	// `0 - dist` to upstream records and `+dist` to downstream ones. Reuse
	// classifyStream so the non-directional and directional paths agree on
	// the convention (-D a flips on A reverse-strand, -D b on B forward-strand).
	if classifyStream(a, b, opts) == streamUpstream {
		return -magnitude
	}
	return magnitude
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
