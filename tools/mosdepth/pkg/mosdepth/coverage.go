package mosdepth

import (
	"sort"

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
	a.events = append(a.events, covEvent{pos: start, delta: 1})
	a.events = append(a.events, covEvent{pos: end, delta: -1})
}

// sub inserts a single (-/+1) event pair at [start, end) on the reference,
// half-open 0-based, subtracting one copy of depth over the interval. It uses
// the same clamping as add, so out-of-range coordinates are clamped to
// [0, refLen] and an empty interval after clamping adds no events. The
// prefix-sum in emit / emitRuns / regionStats accumulates the resulting
// negative deltas naturally; depth never actually goes below 0 because sub is
// only ever applied where two copies were previously added (overlapping mate
// pairs).
func (a *covAccum) sub(start, end int) {
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
	a.events = append(a.events, covEvent{pos: start, delta: -1})
	a.events = append(a.events, covEvent{pos: end, delta: 1})
}

// refBlocks returns the maximal runs of reference-consuming, depth-contributing
// CIGAR ops (M/=/X) as 0-based half-open [start, end) intervals, splitting on
// both deletions (D) and reference-skips (N). This is mosdepth's depth model:
// deleted reference bases get depth 0 (a run breaks on D and N alike), which
// deliberately differs from bedtools' block model (where only N breaks). It is
// the single source of truth for both addRecord and applyOverlap so the two
// stay consistent. Returns nil for an empty / unmapped record.
func refBlocks(rec *sam.Record) [][2]int {
	if rec.Pos <= 0 {
		return nil
	}
	var blocks [][2]int
	refPos := int(rec.Pos) - 1
	for _, op := range rec.Cigar {
		l := int(op.Length())
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			blocks = append(blocks, [2]int{refPos, refPos + l})
			refPos += l
		case sam.CigarDeletion, sam.CigarSkipped:
			refPos += l
		case sam.CigarInsertion, sam.CigarSoftClip, sam.CigarHardClip, sam.CigarPadding:
			// No reference advance.
		}
	}
	return blocks
}

// addRecord walks rec's CIGAR and inserts one event pair per contiguous
// reference-consuming run (M/=/X). Deletions and reference-skips break a
// run because they do not increment depth. In fast mode the whole read
// span from POS to POS+ReferenceLength is added as a single run, skipping
// the CIGAR walk entirely.
func (a *covAccum) addRecord(rec *sam.Record, fast bool) {
	if rec.Pos <= 0 {
		return
	}
	if fast {
		start := int(rec.Pos) - 1
		refLen := rec.Cigar.ReferenceLength()
		if refLen == 0 {
			// Fall back to len(SEQ) when CIGAR is "*" but the read has
			// a sequence — better than nothing in fast mode.
			refLen = len(rec.Seq)
		}
		if refLen <= 0 {
			return
		}
		a.add(start, start+refLen)
		return
	}
	for _, b := range refBlocks(rec) {
		a.add(b[0], b[1])
	}
}

// applyOverlap implements mosdepth's default-mode overlapping mate-pair
// coverage correction: where the two mates of a properly-paired, same-chromosome
// fragment overlap on the reference, the overlapped bases must contribute depth
// 1 rather than 2, so one copy is subtracted. It is a faithful port of upstream
// mosdepth's seen-cache algorithm (gen_start_ends in mosdepth.nim).
//
// seen caches the left mate of each fragment by QNAME; it must be reset per
// chromosome by the caller. applyOverlap is gated to default mode only — the
// caller must not invoke it in fast mode or fragment mode. Records must arrive
// in ascending start coordinate (coordinate-sorted), so the left mate is always
// seen before the right.
func (a *covAccum) applyOverlap(rec *sam.Record, seen map[string]*sam.Record) {
	if !rec.IsProperPair() || rec.IsSupplementary() {
		return
	}
	// The per-chrom loop already restricts rec to one chromosome, so the mate
	// is on the same chromosome iff RNEXT is "=" (the BAM reader sets RNEXT to
	// "=" when nextRefID == refID).
	if rec.RNext != "=" {
		return
	}
	start := int(rec.Pos) - 1                   // 0-based read start.
	stop := start + rec.Cigar.ReferenceLength() // 0-based exclusive read stop.
	matePos := int(rec.PNext) - 1               // 0-based mate start.
	// Left/cache condition: this is the left (cached) mate when it starts at or
	// before its mate and the pair overlaps (stop > matePos). The equal-start
	// tie is broken by which mate was seen first.
	if stop > matePos &&
		(start < matePos || (start == matePos && seen[rec.QName] == nil)) {
		cp := *rec
		seen[rec.QName] = &cp
		return
	}
	mate, ok := seen[rec.QName]
	if !ok {
		return
	}
	delete(seen, rec.QName)
	mateStart := int(mate.Pos) - 1
	mateStop := mateStart + mate.Cigar.ReferenceLength()
	if len(rec.Cigar) == 1 && len(mate.Cigar) == 1 {
		// Single-op fast path: the overlap is exactly [rec.start, mate.stop).
		a.sub(start, mateStop)
		return
	}
	// General path: merge both mates' M/=/X blocks, sweep them, and subtract one
	// copy wherever the pair depth is 2.
	type ovlEvent struct {
		pos   int
		delta int
	}
	blocks := append(refBlocks(rec), refBlocks(mate)...)
	ses := make([]ovlEvent, 0, len(blocks)*2)
	for _, b := range blocks {
		ses = append(ses, ovlEvent{pos: b[0], delta: 1})
		ses = append(ses, ovlEvent{pos: b[1], delta: -1})
	}
	sort.SliceStable(ses, func(i, j int) bool {
		if ses[i].pos != ses[j].pos {
			return ses[i].pos < ses[j].pos
		}
		return ses[i].delta < ses[j].delta
	})
	pairDepth := 0
	lastPos := 0
	for _, ev := range ses {
		if ev.delta == -1 && pairDepth == 2 {
			a.sub(lastPos, ev.pos)
		}
		pairDepth += ev.delta
		lastPos = ev.pos
	}
}

// addFragment adds full-fragment coverage for a single properly-paired,
// non-supplementary read1 record, mirroring upstream mosdepth's
// --fragment-mode. It covers the whole template between the mates: the span
// starts at min(read start, mate start) — both 0-based — and extends for the
// absolute insert size (|TLEN|). The caller is responsible for gating on the
// read1 / proper-pair / supplementary flags before calling this.
func (a *covAccum) addFragment(rec *sam.Record) {
	if rec.Pos <= 0 {
		return
	}
	start := int(rec.Pos) - 1     // 0-based read start.
	matePos := int(rec.PNext) - 1 // 0-based mate start (PNext is 1-based; 0 when unset).
	if matePos >= 0 && matePos < start {
		start = matePos
	}
	isize := int(rec.TLen)
	if isize < 0 {
		isize = -isize
	}
	if isize <= 0 {
		return
	}
	a.add(start, start+isize)
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

// medianHistCap mirrors upstream mosdepth's CountStat histogram size of
// 65536 (initCountStat[uint32](size = 65536) when --use-median is set). Depth
// values are clamped to the highest index, so any per-base depth at or above
// medianHistCap-1 contributes to the top bucket — identical to upstream's
// `c.counts[min(c.counts.high, value)].inc`.
const medianHistCap = 65536

// medianHist is a depth histogram that mirrors upstream mosdepth's
// depthstat.CountStat: counts[d] is the number of bases observed at depth d,
// with depths >= medianHistCap-1 folded into the top bucket. It is fed one
// collapsed [start, end, depth] run at a time (matching regionStats' emitFn
// signature) so a region's median can be derived from the same single sweep
// that computes the mean and per-threshold columns.
type medianHist struct {
	counts []int64
	n      int64
}

// addRun records a constant-depth run of length end-start at the given depth,
// folding out-of-range depths into the histogram's endpoints exactly as
// upstream's `c.counts[min(c.counts.high, value)].inc` does. Its signature
// matches regionStats' emitFn so it can be passed directly as the callback.
func (h *medianHist) addRun(start, end int, depth int32) {
	dv := int(depth)
	if dv < 0 {
		dv = 0
	}
	if dv > medianHistCap-1 {
		dv = medianHistCap - 1
	}
	if dv >= len(h.counts) {
		grown := make([]int64, dv+1)
		copy(grown, h.counts)
		h.counts = grown
	}
	runLen := int64(end - start)
	h.counts[dv] += runLen
	h.n += runLen
}

// median returns the histogram median as a float64, matching upstream
// depthstat.CountStat.median: it walks depths ascending and returns the first
// depth at which the cumulative count reaches stop_n = int(0.5 + n*0.5)
// (round-half-up of n/2). An empty histogram yields 0, matching upstream's
// behaviour when no bases contribute.
func (h *medianHist) median() float64 {
	if h.n == 0 {
		return 0
	}
	stopN := int64(0.5 + float64(h.n)*0.5)
	var cum int64
	for d, cnt := range h.counts {
		cum += cnt
		if cum >= stopN {
			return float64(d)
		}
	}
	return 0
}

// regionMedian computes the per-base depth median across the half-open
// interval [beg0, end0) on this accumulator, returning the integer median
// depth as a float64 (upstream prints it through the same float formatter as
// the mean).
//
// The per-base depth profile is obtained from the same regionStats sweep that
// produces the mean/threshold columns: its emitFn fires once per collapsed
// [start, end, depth] run, so the histogram is built in that single pass
// rather than re-walking the event list. This keeps median and mean computed
// against an identical depth profile by construction.
func (a *covAccum) regionMedian(beg0, end0 int) float64 {
	var h medianHist
	a.regionStats(beg0, end0, nil, h.addRun)
	return h.median()
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
