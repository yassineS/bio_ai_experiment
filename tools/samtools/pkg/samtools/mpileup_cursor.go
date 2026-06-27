package samtools

import (
	"bufio"
	"fmt"
	"io"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/baq"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// runMpileupStreamingCursor piles up a single coordinate-sorted input with
// mate-overlap removal ON, driving emission through a true htslib bam_plp
// push-then-emit cursor rather than the eager tile walk.
//
// The distinction matters only for deletion/refskip placeholders. Their '*'
// (or '>'/'<') borrows the quality of the read's first base AFTER the gap; that
// base may sit in a mate-overlap region and be zeroed/summed by
// tweak_overlap_quality. htslib applies a pair's tweak when the LATER mate is
// pushed (overlap_push, sam.c:6127) and reads the placeholder's borrowed
// quality LIVE at the moment the placeholder's column is emitted
// (bam_get_qual(p->b)[p->qpos], bam_plcmd.c:643). A column at reference position
// P is finalised once a read whose start is > P has been pushed
// (bam_plp64_next's `max_pos > iter->pos`, sam.c:6018). For these placeholders
// the later mate is usually pushed AFTER that trigger, so the borrowed quality
// is still the raw (pre-tweak) value — which is exactly why upstream keeps a '*'
// the eager whole-window tweak would have dropped.
//
// This cursor reproduces that ordering byte-for-byte: it pushes records one at a
// time in coordinate (file) order, applies BAQ and the overlap tweak at the push
// point, accumulates each read's events with a LIVE quality reference
// (accumulateRecordEventsLazy), and finalises each reference column the instant
// a read starting past it is pushed — reading every borrowed quality at that
// emit time. For everything other than the placeholder borrow timing the output
// is identical to the eager walk (aligned-base tweaks always precede their
// column's emit), so the only behavioural change is the residual it closes.
func runMpileupStreamingCursor(rd sam.Reader, out io.Writer, opts MpileupOptions, refFA *fasta.RandomAccess, regionByChrom map[string][][2]int, posFilter *positionFilter) error {
	hdr := rd.Header()
	bw := bufio.NewWriter(out)
	defer bw.Flush()

	doBAQ := refFA != nil && !opts.NoBAQ
	baqFlag := textMpileupBAQFlag(opts)
	src := newMpileupSource(rd, opts, hdr)

	dropChrom := func(chrom string) {
		for p := src.peek(); p != nil && p.RName == chrom; p = src.peek() {
			src.pop()
		}
	}

	for {
		head := src.peek()
		if head == nil {
			break
		}
		chrom := head.RName
		refLen := int(refLengthForName(hdr, chrom))
		if refLen <= 0 {
			dropChrom(chrom)
			continue
		}

		var windows [][2]int
		if regionByChrom != nil {
			ws, ok := regionByChrom[chrom]
			if !ok {
				dropChrom(chrom)
				continue
			}
			windows = ws
		} else {
			windows = [][2]int{{0, refLen}}
		}

		// Fetch the whole contig once for BAQ and the per-row reference base.
		var contig []byte
		if refFA != nil {
			seq, err := refFA.Fetch(chrom, 0, refFA.Length(chrom))
			if err != nil {
				return fmt.Errorf("samtools mpileup: fetch %s: %w", chrom, err)
			}
			contig = seq
		}

		cur := &mpileupCursor{
			bw:        bw,
			chrom:     chrom,
			refLen:    refLen,
			contig:    contig,
			refFA:     refFA,
			posFilter: posFilter,
			opts:      opts,
			overlaps:  make(map[string]*sam.Record),
		}

		for _, w := range windows {
			wBeg, wEnd := w[0], w[1]
			if wBeg < 0 {
				wBeg = 0
			}
			if wEnd > refLen {
				wEnd = refLen
			}
			if wBeg >= wEnd {
				continue
			}
			if err := cur.runWindow(src, chrom, wBeg, wEnd, doBAQ, baqFlag); err != nil {
				return err
			}
		}

		// Advance to the next chromosome.
		dropChrom(chrom)
		if src.err != nil {
			return src.err
		}
	}
	return src.err
}

// mpileupCursor holds the per-chromosome state of the push-then-emit walk: the
// finalisable column buffer (a sliding ring keyed by absolute 0-based reference
// position), the pending-overlap-mate map, and the output context. One cursor
// serves all of a chromosome's windows so a read spanning a window boundary is
// accumulated once.
type mpileupCursor struct {
	bw        *bufio.Writer
	chrom     string
	refLen    int
	contig    []byte
	refFA     *fasta.RandomAccess
	posFilter *positionFilter
	opts      MpileupOptions

	// ring is a sliding window of per-position event buffers, indexed by
	// (abs-ringBase) modulo cap(ring). Columns are accumulated as reads are
	// pushed and freed (truncated for reuse) the instant they are emitted, so
	// purge/advance is O(1) instead of the O(map) churn the previous map[int]
	// keyed buffer incurred on the hot path. htslib's bam_plp streams the same
	// way (a ring indexed by pos-base). readIdx is assigned in push order (==
	// coordinate order), so within any single column events are appended in
	// ascending readIdx — i.e. already in upstream's render order (sortEvents is
	// then a no-op fast-path; see sortEvents).
	ring     [][]pileupEvent // physical backing slots; ring[i] reused after emit
	ringHead int             // physical index of the logical slot 0 (== ringBase)
	ringBase int             // absolute 0-based position mapped to logical slot 0
	ringLen  int             // number of logical slots currently live

	// scratch is the per-read column matrix handed to accumulateRecordEventsLazy.
	// Carrying it on the cursor (reset, never re-allocated per read) removes the
	// per-read make([][]pileupEvent, width) that dominated the push hot path.
	scratch  [][]pileupEvent
	overlaps map[string]*sam.Record
	readIdx  int32
}

// ringSlot returns a pointer to the event buffer for absolute position abs,
// growing the ring as needed so abs is addressable. It assumes abs >= ringBase
// (the cursor never writes to an already-emitted column: reads arrive in
// coordinate order and every column strictly below the latest start is drained
// before the next read is accumulated). When the ring is empty it (re)anchors
// ringBase at abs.
func (c *mpileupCursor) ringSlot(abs int) *[]pileupEvent {
	if c.ringLen == 0 {
		// Empty ring: anchor the head at abs so the first live column is slot 0.
		c.ringBase = abs
		c.ringHead = 0
		if cap(c.ring) == 0 {
			c.ring = make([][]pileupEvent, 16)
		}
		c.ringLen = 1
		idx := c.ringHead
		c.ring[idx] = c.ring[idx][:0]
		return &c.ring[idx]
	}
	off := abs - c.ringBase
	if off >= c.ringLen {
		c.growRing(off + 1)
	}
	idx := c.ringHead + off
	if n := cap(c.ring); idx >= n {
		idx -= n
	}
	return &c.ring[idx]
}

// growRing ensures the ring can address at least need logical slots, extending
// the live span (and growing/re-linearising the backing array when the capacity
// is exhausted). Newly exposed slots are truncated for reuse.
func (c *mpileupCursor) growRing(need int) {
	if need <= cap(c.ring) {
		// Capacity suffices; just extend the live span and clear the new slots.
		for c.ringLen < need {
			idx := c.ringHead + c.ringLen
			if n := cap(c.ring); idx >= n {
				idx -= n
			}
			c.ring[idx] = c.ring[idx][:0]
			c.ringLen++
		}
		return
	}
	// Grow: copy the live slots into a larger, head-aligned backing array.
	newCap := cap(c.ring) * 2
	if newCap < need {
		newCap = need
	}
	nr := make([][]pileupEvent, newCap)
	for i := 0; i < c.ringLen; i++ {
		idx := c.ringHead + i
		if n := cap(c.ring); idx >= n {
			idx -= n
		}
		nr[i] = c.ring[idx]
	}
	c.ring = nr
	c.ringHead = 0
	for c.ringLen < need {
		c.ring[c.ringLen] = c.ring[c.ringLen][:0]
		c.ringLen++
	}
}

// runWindow emits every requested column of [wBeg, wEnd) on the cursor's
// chromosome, pushing records from src as needed. It mirrors bam_plp64_auto:
// push a record (applying BAQ and the overlap tweak), then drain every column
// the push made finalisable (those strictly below the just-pushed start). At the
// window's end any records starting beyond it remain buffered in src for the
// next window or chromosome.
func (c *mpileupCursor) runWindow(src *mpileupSource, chrom string, wBeg, wEnd int, doBAQ bool, baqFlag int) error {
	// Drop any already-buffered events left of the window (a prior window or the
	// leading gap before a -r region). Coordinates only advance, so anything
	// below wBeg can never be emitted now.
	c.purgeBefore(wBeg)

	emitTo := wBeg // next column to emit (0-based)
	for {
		p := src.peek()
		// A record on another chromosome (or EOF) ends this chromosome's pushes;
		// flush every remaining column of the window first.
		if p == nil || p.RName != chrom {
			break
		}
		beg := int(p.Pos) - 1
		if beg >= wEnd {
			// The record starts past this window. It still finalises every column
			// up to wEnd (its start is > any of them); emit them, but leave the
			// record in src for the next window/chromosome.
			break
		}

		// Push this record: BAQ, then the overlap tweak against an already-seen
		// mate (mutating qualities in place), then accumulate its lazy events.
		// Records ending at/below wBeg cannot contribute to this or any later
		// window and are popped without being held.
		if int(p.EndPosition()) > wBeg {
			if doBAQ {
				if r := baq.SamProbRealn(p, c.contig, baqFlag); r < -3 {
					return fmt.Errorf("samtools mpileup: BAQ alignment failed for read %q", p.QName)
				}
			}
			overlapPush(p, c.overlaps)

			// Drain every column strictly below this read's start BEFORE
			// accumulating it. The just-pushed start is the new max_pos, so those
			// columns are finalised now — and the overlap tweak above has already
			// run, so their borrowed qualities reflect it exactly as htslib's
			// bam_plp would at this emit point. Draining first keeps the sliding
			// ring no wider than a single read's span: this read only contributes
			// columns at/above its start, none of which are being drained, so the
			// emitted bytes are unchanged while the ring never has to span the
			// (possibly large) coverage gap below the read.
			drainTo := beg
			if drainTo > wEnd {
				drainTo = wEnd
			}
			for ; emitTo < drainTo; emitTo++ {
				if err := c.emitColumn(emitTo); err != nil {
					return err
				}
			}

			c.accumulate(p, wBeg, wEnd)
			src.pop()
			continue
		}
		src.pop()

		// A record that does not reach the window still advances max_pos: every
		// column strictly below its start is finalised. Emit them up to wEnd.
		drainTo := beg
		if drainTo > wEnd {
			drainTo = wEnd
		}
		for ; emitTo < drainTo; emitTo++ {
			if err := c.emitColumn(emitTo); err != nil {
				return err
			}
		}
	}

	// EOF / end-of-chromosome / record-past-window: every remaining column of
	// the window is finalised.
	for ; emitTo < wEnd; emitTo++ {
		if err := c.emitColumn(emitTo); err != nil {
			return err
		}
	}
	return nil
}

// accumulate appends rec's lazy pileup events to the per-position column buffer.
// The accumulation window is clamped to the read's OWN reference span (capped to
// the emit window [wBeg, wEnd)) rather than the whole window, so a single read
// allocates a scratch proportional to its length — never to the contig. The ^/$
// markers are unaffected: accumulateRecordEventsLazy decides readStart/readEnd
// from the read's first/last in-window event versus its true start/end, which is
// identical whether the window is the read's span or the whole contig (positions
// outside the read contribute no events either way). The merged result is
// byte-identical to the buffered walk.
func (c *mpileupCursor) accumulate(rec *sam.Record, wBeg, wEnd int) {
	beg0 := int(rec.Pos) - 1
	if beg0 < wBeg {
		beg0 = wBeg
	}
	end0 := int(rec.EndPosition()) // 1-based inclusive end == 0-based exclusive end
	if end0 > wEnd {
		end0 = wEnd
	}
	c.readIdx++
	if beg0 >= end0 {
		return
	}
	width := end0 - beg0
	// Reuse a pooled per-read scratch matrix rather than allocating one per read:
	// grow it to the read's span and reset the used columns to empty. The columns
	// hold their own backing arrays across reads, so the per-read append churn is
	// amortised after the first few reads of a contig.
	if cap(c.scratch) < width {
		c.scratch = make([][]pileupEvent, width)
	}
	scratch := c.scratch[:width]
	for col := 0; col < width; col++ {
		scratch[col] = scratch[col][:0]
	}
	accumulateRecordEventsLazy(rec, int(c.readIdx-1), beg0, end0, scratch, c.opts, c.refFA, c.chrom, true, c.contig)
	for col := 0; col < width; col++ {
		if len(scratch[col]) == 0 {
			continue
		}
		abs := beg0 + col
		slot := c.ringSlot(abs)
		*slot = append(*slot, scratch[col]...)
	}
}

// emitColumn writes (or, when uncovered and unrequested, skips) the row for the
// absolute 0-based reference position pos0, then frees its event buffer. The
// emit logic mirrors emitMpileupWindow for a single column so output is
// byte-identical to the buffered renderer; the only difference is that every
// borrowed quality is read live here, at the true emit time.
func (c *mpileupCursor) emitColumn(pos0 int) error {
	// Read this column's events from the ring (nil when no read reached it), then
	// retire every slot up to and including pos0. Emit is strictly monotonic, so
	// once pos0 is rendered every slot at or below it is dead and reusable.
	evs := c.colEvents(pos0)
	c.retireRingThrough(pos0)

	d := liveDepth(evs, c.opts.MinBaseQ)
	if c.opts.MaxDepth > 0 && d > c.opts.MaxDepth {
		d = c.opts.MaxDepth
	}
	spanned := hasSpanEvent(evs)

	pos1 := pos0 + 1
	if c.posFilter != nil && !c.posFilter.contains(c.chrom, pos1) {
		return nil
	}
	if d == 0 && !spanned {
		// The streaming path never serves -a/-aa (runMpileup excludes those), so
		// the only reason to print an empty position is an explicit positions
		// file; otherwise skip it exactly as the buffered walk does.
		if c.posFilter == nil {
			return nil
		}
	}

	ref := byte('N')
	if c.contig != nil && pos0 < len(c.contig) {
		ref = c.contig[pos0]
	} else if c.refFA != nil {
		if b, err := c.refFA.Fetch(c.chrom, int64(pos0), int64(pos0+1)); err == nil && len(b) > 0 {
			ref = b[0]
		}
	}

	bw := c.bw
	if _, err := bw.WriteString(c.chrom); err != nil {
		return err
	}
	bw.WriteByte('\t')
	bw.WriteString(strconv.Itoa(pos1))
	bw.WriteByte('\t')
	bw.WriteByte(ref)
	bw.WriteByte('\t')
	bw.WriteString(strconv.Itoa(d))
	bw.WriteByte('\t')
	if d == 0 {
		bw.WriteByte('*')
	} else {
		writeBasesColumn(bw, evs, ref, c.opts)
	}
	bw.WriteByte('\t')
	if d == 0 {
		bw.WriteByte('*')
	} else {
		writeQualsColumn(bw, evs, c.opts)
	}
	if c.opts.OutputMapQ {
		bw.WriteByte('\t')
		if d == 0 {
			bw.WriteByte('*')
		} else {
			writeMapqColumn(bw, evs, c.opts)
		}
	}
	if c.opts.OutputBP {
		bw.WriteByte('\t')
		if d == 0 {
			bw.WriteByte('*')
		} else {
			writeReadBPColumn(bw, evs, c.opts)
		}
	}
	bw.WriteByte('\n')
	return nil
}

// colEvents returns the accumulated events for absolute position pos0, or nil
// when pos0 lies outside the ring's currently live span (no read has reached it,
// or it has already been retired).
func (c *mpileupCursor) colEvents(pos0 int) []pileupEvent {
	if c.ringLen == 0 {
		return nil
	}
	off := pos0 - c.ringBase
	if off < 0 || off >= c.ringLen {
		return nil
	}
	idx := c.ringHead + off
	if n := cap(c.ring); idx >= n {
		idx -= n
	}
	return c.ring[idx]
}

// retireRingThrough frees (truncates for reuse) and advances past every live
// ring slot whose absolute position is at or below pos0. After it returns the
// ring's head sits at pos0+1 (or the ring is empty if it had no slots beyond
// pos0). This is the O(1)-amortised analogue of the old O(map) purge.
func (c *mpileupCursor) retireRingThrough(pos0 int) {
	for c.ringLen > 0 && c.ringBase <= pos0 {
		c.ring[c.ringHead] = c.ring[c.ringHead][:0]
		c.ringHead++
		if c.ringHead >= cap(c.ring) {
			c.ringHead = 0
		}
		c.ringBase++
		c.ringLen--
	}
	if c.ringLen == 0 {
		// Re-anchor so the next accumulate starts a fresh span at its own start.
		c.ringHead = 0
	}
}

// purgeBefore discards any buffered column events strictly below pos0; they can
// never be emitted once the cursor has moved past them. With the sliding ring
// this just retires the live slots below pos0 (O(retired) work).
func (c *mpileupCursor) purgeBefore(pos0 int) {
	if c.ringLen == 0 {
		return
	}
	if pos0 > c.ringBase {
		c.retireRingThrough(pos0 - 1)
	}
}
