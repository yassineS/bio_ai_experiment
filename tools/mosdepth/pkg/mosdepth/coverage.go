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
	refPos := int(rec.Pos) - 1
	for _, op := range rec.Cigar {
		l := int(op.Length())
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			a.add(refPos, refPos+l)
			refPos += l
		case sam.CigarDeletion, sam.CigarSkipped:
			refPos += l
		case sam.CigarInsertion, sam.CigarSoftClip, sam.CigarHardClip, sam.CigarPadding:
			// No reference advance.
		}
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

// addRecords feeds every record for one reference into the accumulator,
// applying mosdepth's mode-specific coverage rule:
//
//   - fragmentMode: only read1 of a proper, non-supplementary pair contributes,
//     covering the whole template (see addFragment). Callers have already
//     filtered to those reads, so addFragment is called unconditionally.
//   - fastMode: each read covers POS..POS+ReferenceLength with no CIGAR walk and
//     no overlap correction (upstream's --fast-mode skips both).
//   - default mode: each read contributes its CIGAR-aware coverage AND
//     overlapping mate pairs are de-duplicated so a base covered by both mates
//     of a template counts once, exactly as upstream mosdepth does.
//
// recs must be in the BAM's stream (coordinate) order for the overlap detector
// to pair mates correctly, which matches how upstream consumes the file.
func (a *covAccum) addRecords(recs []*sam.Record, fastMode, fragmentMode bool) {
	if fragmentMode {
		for _, rec := range recs {
			a.addFragment(rec)
		}
		return
	}
	if fastMode {
		for _, rec := range recs {
			a.addRecord(rec, true)
		}
		return
	}
	// Default mode: CIGAR-aware coverage plus overlap-pair de-duplication.
	// seen holds the lower-coordinate mate of a still-open overlapping pair,
	// keyed by QName, mirroring upstream's `seen` table.
	seen := map[string]*sam.Record{}
	for _, rec := range recs {
		a.addRecord(rec, false)
		// Overlap handling applies only to proper, non-supplementary pairs whose
		// mate is on the same reference. We approximate upstream's
		// `rec.b.core.tid == rec.b.core.mtid` test with RNext in {"=", RName}.
		if rec.Flag&sam.FlagProperPair == 0 || rec.Flag&sam.FlagSupplementary != 0 {
			continue
		}
		if !mateOnSameRef(rec) {
			continue
		}
		recStart := int(rec.Pos) - 1                      // 0-based start.
		recStop := recStart + rec.Cigar.ReferenceLength() // 0-based exclusive end.
		matePos := int(rec.PNext) - 1                     // 0-based mate start.
		// Mirror upstream's single if/else exactly. The "store the lower mate"
		// branch requires that rec extends past the mate start AND that rec is
		// the lower (or first-seen equal-start) read of the pair.
		_, alreadySeen := seen[rec.QName]
		store := recStop > matePos &&
			(recStart < matePos || (recStart == matePos && !alreadySeen))
		if store {
			cp := *rec
			seen[rec.QName] = &cp
			continue
		}
		// else: this is the higher-coordinate mate (or a non-overlapping read).
		// If its partner was stored, apply the overlap correction and discard it
		// — exactly upstream's `seen.take(rec.qname, mate)`.
		if mate, ok := seen[rec.QName]; ok {
			delete(seen, rec.QName)
			mateStart := int(mate.Pos) - 1
			a.addOverlapCorrection(rec, mate, recStart, mateStart)
		}
	}
}

// mateOnSameRef reports whether rec's mate maps to the same reference, matching
// upstream mosdepth's `rec.b.core.tid == rec.b.core.mtid` test. RNext == "="
// means "same as RName" in SAM/BAM; an explicit RNext equal to RName counts too.
func mateOnSameRef(rec *sam.Record) bool {
	return rec.RNext == "=" || (rec.RNext != "" && rec.RNext == rec.RName)
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

// startEnd is a signed depth event used by the overlap-pair detector: pos is
// the 0-based reference coordinate and value is +1 (a covered run begins) or
// -1 (a covered run ends). It mirrors upstream mosdepth's `pair` tuple emitted
// by gen_start_ends.
type startEnd struct {
	pos   int
	value int32
}

// genStartEnds returns the +1/-1 depth events for one alignment's CIGAR,
// anchored at the 0-based reference start ipos. It is a faithful port of
// upstream mosdepth's gen_start_ends iterator: contiguous reference- AND
// query-consuming ops (M/=/X) are fused into a single covered run, while ops
// that consume reference but not query (D/N) break the run so the gap is not
// counted. Soft/hard clips, insertions, and padding advance neither.
//
// The returned slice always pairs each +1 with a later -1, so the cumulative
// value over the slice is zero. It is the exact event geometry the overlap
// subtractor needs to find the doubly-counted region between two mates.
func genStartEnds(c sam.Cigar, ipos int) []startEnd {
	// Single-M fast path, matching upstream's `c.len == 1 and c[0].op == match`.
	if len(c) == 1 && c[0].Op() == sam.CigarMatch {
		l := int(c[0].Length())
		return []startEnd{{ipos, 1}, {ipos + l, -1}}
	}
	var out []startEnd
	pos := ipos
	lastStop := -1
	for _, op := range c {
		o := op.Op()
		consumesRef := o == sam.CigarMatch || o == sam.CigarEqual || o == sam.CigarMismatch ||
			o == sam.CigarDeletion || o == sam.CigarSkipped
		if !consumesRef {
			continue
		}
		olen := int(op.Length())
		consumesQuery := o == sam.CigarMatch || o == sam.CigarEqual || o == sam.CigarMismatch
		if consumesQuery {
			if pos != lastStop {
				out = append(out, startEnd{pos, 1})
				if lastStop != -1 {
					out = append(out, startEnd{lastStop, -1})
				}
			}
			lastStop = pos + olen
		}
		pos += olen
	}
	if lastStop != -1 {
		out = append(out, startEnd{lastStop, -1})
	}
	return out
}

// addOverlapCorrection removes the depth double-counted where the two mates of
// a read pair overlap, mirroring upstream mosdepth's default-mode overlap
// handling. mate is the lower-coordinate read (already seen) and rec is the
// higher-coordinate read of the same template; both have already contributed
// their full per-base coverage via addRecord. This inserts the compensating
// -1/+1 events so the overlapping bases are counted once, not twice.
//
// recStart and mateStart are the 0-based reference starts of rec and mate.
// When both reads have a single CIGAR op (the common gapless case) upstream
// takes a shortcut: subtract one copy across [recStart, mateStop). Otherwise it
// merges both reads' gen_start_ends events, sorts them, and subtracts one copy
// of every span where the combined pair depth reaches 2.
func (a *covAccum) addOverlapCorrection(rec, mate *sam.Record, recStart, mateStart int) {
	if len(rec.Cigar) == 1 && len(mate.Cigar) == 1 {
		// mate:   --------------
		// rec:             ------------
		// decrement:       -----  (from rec.start to mate.stop). Upstream does
		// dec(arr[rec.start]); inc(arr[mate.stop]); i.e. a -1 run over the
		// overlap span.
		mateStop := mateStart + mate.Cigar.ReferenceLength()
		a.events = append(a.events,
			covEvent{pos: recStart, delta: -1},
			covEvent{pos: mateStop, delta: 1})
		return
	}
	ses := genStartEnds(rec.Cigar, recStart)
	ses = append(ses, genStartEnds(mate.Cigar, mateStart)...)
	sortStartEnds(ses)
	var pairDepth int32
	lastPos := 0
	for _, p := range ses {
		// When pair depth is 2 and a run closes, [lastPos, p.pos) is the
		// doubly-covered span: subtract exactly one copy of it.
		if p.value == -1 && pairDepth == 2 {
			a.events = append(a.events,
				covEvent{pos: lastPos, delta: -1},
				covEvent{pos: p.pos, delta: 1})
		}
		pairDepth += p.value
		lastPos = p.pos
	}
}

// sortStartEnds sorts a slice of startEnd by ascending position, matching
// upstream's pair_sort (which orders solely on pos). Ties keep input order,
// which is sufficient because the overlap walk only reacts to -1 events while
// pair depth is already 2.
func sortStartEnds(s []startEnd) {
	if len(s) < 2 {
		return
	}
	// Simple insertion sort: these slices have at most a handful of events
	// (one per CIGAR block across two reads).
	for i := 1; i < len(s); i++ {
		for j := i; j > 0 && s[j-1].pos > s[j].pos; j-- {
			s[j-1], s[j] = s[j], s[j-1]
		}
	}
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
