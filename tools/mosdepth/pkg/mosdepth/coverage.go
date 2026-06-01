package mosdepth

import (
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// covEvent is a single signed depth event at a 0-based reference position:
// +1 when an alignment run starts, -1 when it ends. The accumulator sweeps
// these in order to recover per-position depth without ever materialising a
// depth-per-base array.
type covEvent struct {
	pos   int
	delta int32
}

// covAccum accumulates per-position depth events for one reference
// sequence. Events may be inserted out of order; Sort sorts them in place
// before emission. The struct is intentionally simple — the higher-level
// processor owns one of these per reference being scanned.
type covAccum struct {
	events []covEvent
	// refLen caps the highest position we will emit. Events beyond refLen
	// are clamped to refLen so the closing event still subtracts the
	// open count at the boundary.
	refLen int
}

// newCovAccum constructs an accumulator for a reference of refLen bases.
// refLen <= 0 means "unknown" and the accumulator will emit up to the
// largest event position observed.
func newCovAccum(refLen int) *covAccum {
	return &covAccum{refLen: refLen}
}

// add inserts a single (+/-1) event pair at [start, end) on the reference,
// half-open 0-based. start < end is required; out-of-range coordinates are
// clamped to [0, refLen]. If the resulting interval is empty after
// clamping, no events are added.
func (a *covAccum) add(start, end int) {
	a.addSigned(start, end, 1)
}

// addSigned inserts a (+sign, -sign) event pair at [start, end). It is used
// by the overlap-pair detector to cancel one copy of depth across the
// region shared by two mates (sign = -1). The semantics mirror add — same
// boundary clamping, same drop-empty rule.
func (a *covAccum) addSigned(start, end int, sign int32) {
	if a.refLen > 0 {
		if start < 0 {
			start = 0
		}
		if end > a.refLen {
			end = a.refLen
		}
	} else {
		if start < 0 {
			start = 0
		}
	}
	if end <= start {
		return
	}
	a.events = append(a.events, covEvent{pos: start, delta: sign})
	a.events = append(a.events, covEvent{pos: end, delta: -sign})
}

// addRecord walks rec's CIGAR and inserts one event pair per contiguous
// reference-consuming run (M/=/X). Deletions and reference-skips break a
// run because they do not increment depth. In fast mode the whole read
// span from POS to POS+ReferenceLength is added as a single run, skipping
// the CIGAR walk entirely.
func (a *covAccum) addRecord(rec *sam.Record, fast bool) {
	for _, iv := range recordRefIntervals(rec, fast) {
		a.add(iv[0], iv[1])
	}
}

// recordRefIntervals returns the list of half-open [start, end) reference
// intervals contributed by rec under the current mode. Empty result means
// the read does not contribute any depth (unmapped, empty CIGAR, etc.).
//
// In fast mode the result is the single span POS..POS+ReferenceLength; in
// the default mode it is one interval per contiguous M/=/X run (D/N break
// the run, I/S/H/P do not advance the reference). It is exported as a
// package-private helper so the overlap-pair detector can reuse it without
// duplicating the CIGAR walk.
func recordRefIntervals(rec *sam.Record, fast bool) [][2]int {
	if rec == nil || rec.Pos <= 0 {
		return nil
	}
	if fast {
		start := int(rec.Pos) - 1
		refLen := rec.Cigar.ReferenceLength()
		if refLen == 0 {
			refLen = len(rec.Seq)
		}
		if refLen <= 0 {
			return nil
		}
		return [][2]int{{start, start + refLen}}
	}
	refPos := int(rec.Pos) - 1
	var out [][2]int
	for _, op := range rec.Cigar {
		l := int(op.Length())
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			if l > 0 {
				out = append(out, [2]int{refPos, refPos + l})
			}
			refPos += l
		case sam.CigarDeletion, sam.CigarSkipped:
			refPos += l
		case sam.CigarInsertion, sam.CigarSoftClip, sam.CigarHardClip, sam.CigarPadding:
			// No reference advance.
		}
	}
	return out
}

// fragmentIntervals returns the reference intervals contributed by rec
// when --fragment-mode is in effect. The convention matches upstream
// mosdepth's fragment view:
//
//   - Paired reads with TLEN > 0 own the whole [POS-1, POS-1+TLEN) span;
//     their mate (TLEN < 0) returns no intervals so the fragment is
//     counted exactly once.
//   - Singletons, mate-unmapped reads, and any record with TLEN == 0
//     fall back to the CIGAR-walk view (single contiguous span if mapped).
//
// This produces depth that counts each *fragment* once at every base it
// physically covers — overlapping mates collapse for free.
func fragmentIntervals(rec *sam.Record) [][2]int {
	if rec == nil || rec.Pos <= 0 {
		return nil
	}
	if rec.Flag&sam.FlagPaired != 0 && rec.Flag&sam.FlagMateUnmapped == 0 && rec.TLen != 0 {
		// Same-chrom mates carry RNext == "=" or RNext == RName.
		if rec.RNext == "=" || rec.RNext == "" || rec.RNext == rec.RName {
			if rec.TLen < 0 {
				// Right mate — left mate already covered the fragment.
				return nil
			}
			start := int(rec.Pos) - 1
			end := start + int(rec.TLen)
			if end <= start {
				return nil
			}
			return [][2]int{{start, end}}
		}
	}
	// Singleton / cross-chrom / mate-unmapped / TLEN==0 fallback.
	start := int(rec.Pos) - 1
	refLen := rec.Cigar.ReferenceLength()
	if refLen == 0 {
		refLen = len(rec.Seq)
	}
	if refLen <= 0 {
		return nil
	}
	return [][2]int{{start, start + refLen}}
}

// overlapIntervals returns the intersection (in half-open [start, end)
// coordinates) of two interval lists. Both inputs are expected to be
// position-ascending (the natural order from a CIGAR walk). The result is
// also ascending; empty when there is no overlap.
//
// Used by the overlap-pair detector to compute the bases shared by two
// mates of the same fragment so the second contribution can be cancelled.
func overlapIntervals(a, b [][2]int) [][2]int {
	if len(a) == 0 || len(b) == 0 {
		return nil
	}
	var out [][2]int
	i, j := 0, 0
	for i < len(a) && j < len(b) {
		s := a[i][0]
		if b[j][0] > s {
			s = b[j][0]
		}
		e := a[i][1]
		if b[j][1] < e {
			e = b[j][1]
		}
		if s < e {
			out = append(out, [2]int{s, e})
		}
		if a[i][1] < b[j][1] {
			i++
		} else {
			j++
		}
	}
	return out
}

// sortEvents sorts events by position ascending. Equal positions keep their
// relative order so emit() applies all deltas at a position atomically; a
// stable sort isn't required because emit() collapses ties.
func (a *covAccum) sortEvents() {
	// In-place insertion sort would be cheap for the often-mostly-sorted
	// input we receive from coordinate-sorted BAMs, but stdlib's sort.Sort
	// is simple and good enough.
	sortEventSlice(a.events)
}

// emit walks the (sorted) event list and invokes fn for every 0-based
// position in [0, refLen) with the depth value at that position. Positions
// where the depth equals the previous emitted depth ARE still emitted —
// callers (the per-base writer) collapse runs themselves; emit just
// guarantees ordered, complete position coverage.
//
// When the depth across a contiguous run is constant, the per-base writer
// fuses them into a single BED interval matching upstream mosdepth's
// per-base.bed.gz format.
func (a *covAccum) emit(fn func(pos int, depth int32)) {
	a.sortEvents()
	upper := a.refLen
	if upper <= 0 {
		// Find max event position.
		for _, ev := range a.events {
			if ev.pos > upper {
				upper = ev.pos
			}
		}
	}
	if upper <= 0 {
		return
	}
	var depth int32
	idx := 0
	for pos := 0; pos < upper; pos++ {
		// Apply all events whose pos == this position.
		for idx < len(a.events) && a.events[idx].pos == pos {
			depth += a.events[idx].delta
			idx++
		}
		// Apply events with pos < this position too (defensive: handles
		// clipped events).
		for idx < len(a.events) && a.events[idx].pos < pos {
			depth += a.events[idx].delta
			idx++
		}
		fn(pos, depth)
	}
}

// emitRuns is a convenience wrapper around emit that collapses contiguous
// equal-depth positions into [start, end, depth) tuples. This matches the
// per-base.bed.gz format upstream mosdepth produces: BED records with the
// 4th column being the integer depth and zero-depth runs included.
func (a *covAccum) emitRuns(fn func(start, end int, depth int32)) {
	var runStart int = -1
	var runDepth int32
	emitter := func(pos int, depth int32) {
		if runStart < 0 {
			runStart = pos
			runDepth = depth
			return
		}
		if depth != runDepth {
			fn(runStart, pos, runDepth)
			runStart = pos
			runDepth = depth
		}
	}
	a.emit(emitter)
	if runStart >= 0 {
		// Flush trailing run.
		end := a.refLen
		if end <= runStart {
			end = runStart + 1
		}
		fn(runStart, end, runDepth)
	}
}

// regionStats computes the per-base depth statistics across the half-open
// interval [beg0, end0) on this accumulator, returning the sum of depths
// (so caller can derive a mean) and per-threshold counts (number of bases
// whose depth >= threshold[i]). The minimum and maximum observed depth are
// also returned. emitFn, when non-nil, is invoked once per [runStart,
// runEnd, depth] collapsed run inside the region.
//
// The implementation is a separate sweep so it can be applied to regions
// that don't span the whole reference without re-emitting the per-base
// output.
func (a *covAccum) regionStats(beg0, end0 int, thresholds []int, emitFn func(start, end int, depth int32)) (sum int64, perThreshold []int64, minD, maxD int32) {
	if beg0 < 0 {
		beg0 = 0
	}
	if a.refLen > 0 && end0 > a.refLen {
		end0 = a.refLen
	}
	if end0 <= beg0 {
		perThreshold = make([]int64, len(thresholds))
		return
	}
	perThreshold = make([]int64, len(thresholds))
	// Sort events once.
	a.sortEvents()
	var depth int32
	idx := 0
	// Advance to depth at beg0 by applying every event with pos <= beg0.
	for idx < len(a.events) && a.events[idx].pos <= beg0 {
		depth += a.events[idx].delta
		idx++
	}
	pos := beg0
	first := true
	for pos < end0 {
		// Next event position inside the region defines a constant run.
		nextPos := end0
		if idx < len(a.events) && a.events[idx].pos < end0 {
			nextPos = a.events[idx].pos
		}
		if nextPos <= pos {
			nextPos = pos + 1
		}
		runLen := nextPos - pos
		sum += int64(depth) * int64(runLen)
		for ti, th := range thresholds {
			if int(depth) >= th {
				perThreshold[ti] += int64(runLen)
			}
		}
		if first {
			minD = depth
			maxD = depth
			first = false
		} else {
			if depth < minD {
				minD = depth
			}
			if depth > maxD {
				maxD = depth
			}
		}
		if emitFn != nil {
			emitFn(pos, nextPos, depth)
		}
		pos = nextPos
		// Apply every event at pos.
		for idx < len(a.events) && a.events[idx].pos == pos {
			depth += a.events[idx].delta
			idx++
		}
	}
	return sum, perThreshold, minD, maxD
}

// medianHistSize is the number of histogram bins mosdepth allocates for a
// CountStat when --use-median is requested (initCountStat[uint32](size =
// 65536)). Depth values >= medianHistSize-1 are clamped into the final bin,
// exactly mirroring upstream's `c.counts[min(c.counts.high, value)].inc`.
const medianHistSize = 65536

// regionMedian computes the median per-base depth across the half-open
// interval [beg0, end0) using the identical algorithm to upstream mosdepth's
// CountStat.median:
//
//	stop_n = int(0.5 + n*0.5)        // round-half-up of n/2
//	walk depth bins ascending, accumulating counts;
//	return the first bin whose cumulative count >= stop_n.
//
// For an even number of bases this yields the LOWER of the two middle values
// (mosdepth does not average the two central order statistics). Depths are
// bucketed into a histogram capped at medianHistSize-1 bins so that pathological
// ultra-high depths clamp the same way upstream's fixed-size CountStat does.
//
// The returned value is the integer median depth as a float64 so callers can
// format it with the same precision as the mean (e.g. "1.00"). When the region
// is empty the result is 0, matching upstream's `if start > len: return 0`.
func (a *covAccum) regionMedian(beg0, end0 int) float64 {
	if beg0 < 0 {
		beg0 = 0
	}
	if a.refLen > 0 && end0 > a.refLen {
		end0 = a.refLen
	}
	if end0 <= beg0 {
		return 0
	}
	a.sortEvents()
	// Build the per-base depth histogram for the region. The histogram is
	// grown lazily so the common shallow-coverage case stays cheap; bins are
	// still capped at medianHistSize-1 to match upstream's CountStat sizing.
	hist := make([]int64, 0, 16)
	bump := func(depth int32, count int) {
		d := int(depth)
		if d < 0 {
			d = 0
		}
		if d >= medianHistSize-1 {
			d = medianHistSize - 1
		}
		if d >= len(hist) {
			grown := make([]int64, d+1)
			copy(grown, hist)
			hist = grown
		}
		hist[d] += int64(count)
	}
	var depth int32
	idx := 0
	for idx < len(a.events) && a.events[idx].pos <= beg0 {
		depth += a.events[idx].delta
		idx++
	}
	pos := beg0
	for pos < end0 {
		nextPos := end0
		if idx < len(a.events) && a.events[idx].pos < end0 {
			nextPos = a.events[idx].pos
		}
		if nextPos <= pos {
			nextPos = pos + 1
		}
		bump(depth, nextPos-pos)
		pos = nextPos
		for idx < len(a.events) && a.events[idx].pos == pos {
			depth += a.events[idx].delta
			idx++
		}
	}
	n := end0 - beg0
	stopN := int(0.5 + float64(n)*0.5)
	cum := 0
	for d, c := range hist {
		cum += int(c)
		if cum >= stopN {
			return float64(d)
		}
	}
	return 0
}

// sortEventSlice sorts a slice of covEvent by ascending position.
//
// This is a small, allocation-free quicksort tuned for the mostly-sorted
// data that comes out of coordinate-sorted BAMs. Falling back to insertion
// sort below a small threshold keeps it cache-friendly.
func sortEventSlice(s []covEvent) {
	if len(s) < 2 {
		return
	}
	// Use the stdlib sort.Slice; the project tolerates the small extra
	// overhead in exchange for simplicity.
	// (Mosdepth produces O(reads) events, ~1e7 for a 30x WGS BAM — well
	// within stdlib sort's comfort range.)
	stdSortEvents(s)
}
