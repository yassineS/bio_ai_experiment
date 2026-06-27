package samtools

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// tmpEvent is the per-read CIGAR-walk scratch record built by
// accumulateRecordEventsLazy: one (pos0, kind, ...) tuple per reference-consuming
// op before first/last (^/$) markers are resolved. It is package scope (rather
// than a function-local type) so its backing slice can be pooled across reads.
type tmpEvent struct {
	pos0            int
	col             int // column index into evs (-1 when out of window)
	kind            pileupEventKind
	base            byte
	qual            byte
	qpos            int // query index whose quality this event borrows (lazy path)
	readBP          int
	insAfter        string
	delAfter        int
	delBases        string
	refSkipBoundary bool
}

// tmpEventPool recycles the per-read tmpEvent scratch slice. accumulateRecord-
// EventsLazy runs once per read (across the whole contig in the streaming
// cursor), so allocating its scratch fresh each call was a leading mallocgc
// source; the pool keeps that allocation off the hot path while staying safe for
// the single in-flight call per goroutine.
var tmpEventPool = sync.Pool{New: func() any { s := make([]tmpEvent, 0, 16); return &s }}

// pileupEventKind tags what a pileup event represents at a given position.
type pileupEventKind uint8

const (
	// pileupEventBase is a base read aligned to the reference at this
	// position via an M/=/X CIGAR op. Base/Qual/MapQ are valid.
	pileupEventBase pileupEventKind = iota
	// pileupEventDel is a placeholder base for a position skipped by a D
	// or N CIGAR op (the read covers the position but contributes no
	// base). Reported as '*' in the bases column. Qual is `!` (0+33).
	pileupEventDel
	// pileupEventRefSkip is a position inside an N CIGAR run. Reported
	// the same as a deletion ('*') in upstream pileup output but tracked
	// separately because it does not generate a "-<len>" indel
	// annotation at the prior position.
	pileupEventRefSkip
)

// pileupEvent is one read's contribution to one reference position.
//
// Field order is deliberate: the two (usually empty) string fields and the
// 32-bit integer fields are grouped first so the struct packs to 56 bytes with
// no interior padding. A whole-chromosome scan keeps a persistent column matrix
// of these events, so every byte trimmed here multiplies across the resident
// set; the integer fields are 32-bit because a read index, a within-read base
// position, and a deletion length all fit comfortably in int32 (a contig holds
// far fewer than 2^31 reads, and read/deletion lengths are tiny). This is purely
// a memory-layout choice and never affects the emitted bytes.
type pileupEvent struct {
	// insAfter, when non-empty, is the inserted sequence to annotate as
	// "+<len><seq>" immediately after this event's base. Only valid on
	// pileupEventBase.
	insAfter string
	// delBases holds the reference bases of a deletion that begins at the
	// position immediately AFTER this event's base (see delAfter), used for
	// accurate "-<len><seq>" rendering when a FASTA reference is loaded;
	// otherwise the renderer falls back to 'N' x delAfter, matching upstream.
	delBases string
	// readIdx is the 0-based index of the originating record in the
	// per-chrom record slice; used to order events stably within a
	// position so multi-read columns match upstream's "in-order"
	// rendering.
	readIdx int32
	// readBP is the 1-based position of this base within the read's
	// SEQ (i.e. the original query base index + 1, regardless of strand).
	// 0 for non-base events.
	readBP int32
	// delAfter, when > 0, is the length of a deletion that begins at the
	// position immediately AFTER this event's base.
	delAfter int32
	// qrec, when non-nil, makes the event's effective quality a LAZY read of
	// qrec.Qual[qpos] at emit time instead of the frozen qual field. The
	// streaming overlap-active cursor (runMpileupStreamingCursor) sets this so a
	// deletion/refskip placeholder borrows its post-gap base quality AS OF the
	// moment its column is emitted — i.e. reflecting exactly the mate-overlap
	// tweaks htslib's bam_plp had applied by that column (tweaks fire at the
	// later mate's push, which for these placeholders happens after the column
	// emits). The buffered walk and the consensus caller leave qrec nil, so they
	// read the frozen qual and stay byte-identical. Only the cursor populates it.
	qrec *sam.Record
	// qpos is the query-sequence index whose quality this event borrows when
	// qrec is non-nil (the base used by htslib's bam_get_qual(p->b)[p->qpos]).
	qpos int32
	kind pileupEventKind
	// base is the read base (uppercase if forward, lowercase if reverse)
	// at this reference position. '*' for pileupEventDel/RefSkip.
	base byte
	// qual is the Phred quality of the read base; 0 for non-base events. When
	// qrec is non-nil this is ignored in favour of the lazy read (see qrec).
	qual byte
	// mapq is the read's MAPQ.
	mapq uint8
	// isReverse reports the read's reverse-strand flag bit.
	isReverse bool
	// readStart is true when this is the very first event emitted for
	// this read; the output gets a leading "^<mapq>" marker.
	readStart bool
	// readEnd is true when this is the very last event emitted for this
	// read; the output gets a trailing "$".
	readEnd bool
	// dropped, when true, indicates this event was deselected by post-
	// hoc filtering (overlap-pair removal or max-depth cap). Stays in
	// the slice so readStart/readEnd markers can still be skipped/moved.
	dropped bool
	// noSeq is true when the originating record had SEQ=* (no sequence
	// info); the bases column is always rendered as explicit N/n in
	// that case rather than the dot/comma match form.
	noSeq bool
	// refSkipBoundary marks an aligned base (pileupEventBase) that sits
	// immediately before or after a ref-skip (CIGAR N) run for this read.
	// Upstream's consensus pileup engine flags such a base with p->ref_skip
	// (consensus_pileup.c:239-244, 251-260): the Gap5/bayesian caller then
	// excludes it from its depth (bam_consensus.c:1333), so an exon base at
	// the very edge of an intron does not count toward the bayesian
	// consensus. The base is still displayed verbatim in the pileup seq
	// column (basic_pileup emits p->base, the real base), and the SIMPLE
	// caller ignores the flag entirely — only the bayesian depth is
	// affected. This is read solely by the consensus bayesian path; mpileup
	// does not consult it.
	refSkipBoundary bool
}

// emitMpileupWindow walks every position in [beg0, end0) on chrom, gathers
// per-input pileup events, and emits one text mpileup line per position
// that has at least one read (or every position when opts.AllPositions /
// AllPositionsAllChroms is set, or the position is in posFilter).
// mpileupScratch holds the per-window buffers emitMpileupWindow reuses across
// tiles: the per-input, per-position event matrix and the per-position depth
// slice. Reusing them keeps tiling from re-allocating (and re-zeroing) a fresh
// event matrix for every tile, which would otherwise dominate the tiled path's
// CPU.
type mpileupScratch struct {
	events [][][]pileupEvent // [input][column][event]
	depths []int
}

// mpileupColumnCapBound is the per-column event-array capacity above which
// reset reallocates the column small instead of recycling it. The scratch is
// persistent for a whole-chromosome scan, and reset only re-slices each
// column's backing array to [:0]; without this clamp every column's backing
// array would pin the deepest coverage ever seen at its tile offset, and since
// the same width-wide scratch is reused across every tile of the contig, that
// high-water mark ratchets toward the GLOBAL maximum depth on the whole contig.
// A region of normal depth ~tens then carries the resident cost of the single
// deepest pileup spike anywhere on the chromosome, for every column, for the
// rest of the run.
//
// 64 sits comfortably above the typical per-position read depth of WGS data
// (single-to-low-double-digit coverage), so the overwhelmingly common shallow
// column still recycles its small backing array cheaply via [:0]; only columns
// that transiently held a deep pileup are dropped and re-grown on demand. The
// re-grow is a rare, amortised geometric reallocation confined to genuine
// coverage spikes, so wall time is unaffected while the resident event matrix
// tracks the local — not the global — depth. This changes only how backing
// arrays are recycled, never the emitted bytes.
const mpileupColumnCapBound = 64

// reset prepares the scratch for a window of nIn inputs and the given width,
// growing the buffers when needed and truncating each per-position event slice
// to length zero while keeping its backing array for reuse. A column whose
// backing array has ratcheted far past mpileupColumnCapBound is dropped (set to
// nil) instead of recycled, so a single deep column cannot pin an oversized
// array for the rest of a long whole-chromosome scan.
func (s *mpileupScratch) reset(nIn, width int) {
	if cap(s.depths) < nIn {
		s.depths = make([]int, nIn)
	}
	s.depths = s.depths[:nIn]
	if cap(s.events) < nIn {
		s.events = make([][][]pileupEvent, nIn)
	}
	s.events = s.events[:nIn]
	for i := 0; i < nIn; i++ {
		if cap(s.events[i]) < width {
			s.events[i] = make([][]pileupEvent, width)
		} else {
			s.events[i] = s.events[i][:width]
		}
		for c := 0; c < width; c++ {
			// Common case: small backing array, recycle it cheaply.
			// Memory-ratchet guard: if this column once held a
			// pathologically deep pileup, drop the oversized array so a
			// single deep column does not pin memory for the whole run.
			if cap(s.events[i][c]) > mpileupColumnCapBound {
				s.events[i][c] = nil
			} else {
				s.events[i][c] = s.events[i][c][:0]
			}
		}
	}
}

// emitMpileupWindow writes the pileup rows for [beg0, end0). contig, when
// non-nil, is the whole-chromosome reference sequence (uppercased, newline
// stripped) and is sliced for the per-row reference base instead of re-fetching
// from refFA each tile; sc is reusable scratch shared across the window's tiles.
func emitMpileupWindow(bw *bufio.Writer, chrom string, beg0, end0, refLen int,
	perInputChromRecs [][]*sam.Record, contig []byte, refFA *fasta.RandomAccess,
	posFilter *positionFilter, opts MpileupOptions, sc *mpileupScratch) error {

	nIn := len(perInputChromRecs)

	// Build a per-input slice of per-position event slices for the window.
	// We accumulate by walking each read's CIGAR once and appending each
	// event to events[pos-beg0]. width can be very large for chrom-wide
	// scans; for typical region queries it is small.
	width := end0 - beg0
	sc.reset(nIn, width)
	events := sc.events
	for i := 0; i < nIn; i++ {
		for ridx, rec := range perInputChromRecs[i] {
			accumulateRecordEvents(rec, ridx, beg0, end0, events[i], opts, refFA, chrom, contig)
		}
	}

	// Mate-pair overlap removal is applied to the records' qualities before this
	// per-tile accumulation (see runMpileup / runMpileupStreaming), matching
	// htslib's bam_plp, so nothing is dropped here.

	// Optionally fetch a reference slab so each line can carry the right
	// refbase column. Empty when no FASTA was supplied (upstream emits 'N').
	var refSlab []byte
	switch {
	case contig != nil && end0 <= len(contig):
		// Slice the pre-fetched whole-chromosome reference (no per-tile I/O).
		refSlab = contig[beg0:end0]
	case refFA != nil:
		// Coerce [beg0, end0) onto the contig length already done by the
		// caller, so this Fetch is always in-range.
		b, err := refFA.Fetch(chrom, int64(beg0), int64(end0))
		if err != nil {
			return fmt.Errorf("samtools mpileup: fetch %s:%d-%d: %w", chrom, beg0, end0, err)
		}
		refSlab = b
	}

	depths := sc.depths
	for pos0 := beg0; pos0 < end0; pos0++ {
		col := pos0 - beg0
		// Gather depths per input (depths is reused; every entry is overwritten
		// below before use).
		any := false
		spanned := false
		for i := 0; i < nIn; i++ {
			d := liveDepth(events[i][col], opts.MinBaseQ)
			if opts.MaxDepth > 0 && d > opts.MaxDepth {
				d = opts.MaxDepth
			}
			depths[i] = d
			if d > 0 {
				any = true
			}
			// A position is "spanned" when at least one read physically
			// overlaps it (any non-dropped pileup event), regardless of the
			// base-quality filter. Upstream's pileup iterator yields exactly
			// these positions: bam_mplp64_auto returns every position that
			// some read covers (n_plp > 0), and the row is always printed —
			// the per-base min-baseQ filter (bam_plcmd.c:646) only lowers the
			// printed depth (cnt), it never suppresses the row. A position
			// whose only covering base(s) fail -Q therefore prints depth 0
			// with '*' columns, not nothing.
			if hasSpanEvent(events[i][col]) {
				spanned = true
			}
		}

		// Decide whether to emit this position.
		pos1 := pos0 + 1
		if posFilter != nil && !posFilter.contains(chrom, pos1) {
			continue
		}
		if !any && !spanned {
			switch {
			case opts.AllPositionsAllChroms:
				// emit zero-depth row
			case opts.AllPositions:
				// Upstream `-a` emits EVERY position (1..LN) of every chrom
				// that carries at least one read — not merely the covered
				// extent (bam_plcmd.c: with conf->all==1 the missing-portion
				// loops at lines 603-631 fill the leading/interior gaps and
				// the post-loop flush at 845-857 fills the trailing gap up to
				// sam_hdr_tid2len, i.e. the full contig). Chroms with no reads
				// are excluded from chromsToWalk, so reaching this branch
				// already implies the contig is read-bearing; emit the row.
				// emit zero-depth row
			case posFilter != nil:
				// In positions-file mode, emit zero-depth rows so the
				// caller can see "no coverage here", matching upstream.
			default:
				continue
			}
		}

		// Refbase column.
		ref := byte('N')
		if refSlab != nil && col < len(refSlab) {
			ref = refSlab[col]
		}

		if _, err := bw.WriteString(chrom); err != nil {
			return err
		}
		bw.WriteByte('\t')
		bw.WriteString(strconv.Itoa(pos1))
		bw.WriteByte('\t')
		bw.WriteByte(ref)
		for i := 0; i < nIn; i++ {
			bw.WriteByte('\t')
			bw.WriteString(strconv.Itoa(depths[i]))
			bw.WriteByte('\t')
			if depths[i] == 0 {
				bw.WriteByte('*')
			} else {
				writeBasesColumn(bw, events[i][col], ref, opts)
			}
			bw.WriteByte('\t')
			if depths[i] == 0 {
				bw.WriteByte('*')
			} else {
				writeQualsColumn(bw, events[i][col], opts)
			}
			if opts.OutputMapQ {
				bw.WriteByte('\t')
				if depths[i] == 0 {
					bw.WriteByte('*')
				} else {
					writeMapqColumn(bw, events[i][col], opts)
				}
			}
			if opts.OutputBP {
				bw.WriteByte('\t')
				if depths[i] == 0 {
					bw.WriteByte('*')
				} else {
					writeReadBPColumn(bw, events[i][col], opts)
				}
			}
		}
		bw.WriteByte('\n')
	}
	return nil
}

// effectiveQual returns the Phred quality that drives this event's -Q filter
// and quals-column rendering. When qrec is non-nil the value is read live from
// qrec.Qual[qpos] (the streaming overlap cursor's lazy borrow, matching
// htslib's emit-time bam_get_qual(p->b)[p->qpos]); otherwise the frozen qual is
// used, so the buffered and consensus callers are byte-unchanged.
func (e *pileupEvent) effectiveQual() byte {
	if e.qrec != nil {
		if int(e.qpos) >= 0 && int(e.qpos) < len(e.qrec.Qual) {
			return e.qrec.Qual[e.qpos]
		}
		return 0
	}
	return e.qual
}

// hasSpanEvent reports whether any non-dropped pileup event covers this
// position, i.e. at least one read physically overlaps it (an aligned base,
// a deletion '*', or a reference-skip '<'/'>' placeholder). This is the
// "position exists in the pileup" test, independent of the -Q base-quality
// filter: upstream emits a row for every such position, printing depth 0
// (and '*' columns) when every covering base is filtered out. Events
// dropped by overlap removal (-x) do not keep the position alive on their
// own — they were never part of the displayed pileup — but the surviving
// half of the pair still does.
func hasSpanEvent(evs []pileupEvent) bool {
	for i := range evs {
		if !evs[i].dropped {
			return true
		}
	}
	return false
}

// liveDepth returns the count of "live" (not filtered, base-quality OK)
// events at a position. pileupEventDel/RefSkip occupy a depth slot too (they
// are written as '*'/'>'/'<'), but — like aligned bases — only when their
// quality clears the -Q/--min-BQ threshold: upstream applies min_baseQ to
// bam_get_qual[p->qpos] for every entry, and our del/refskip events carry that
// post-gap base's quality, so a placeholder flanking a low-quality deletion is
// excluded from depth exactly as upstream excludes it.
func liveDepth(evs []pileupEvent, minBQ uint8) int {
	n := 0
	for i := range evs {
		if evs[i].dropped {
			continue
		}
		if minBQ > 0 && evs[i].effectiveQual() < minBQ {
			continue
		}
		n++
	}
	return n
}

// writeBasesColumn writes the 5th column ("bases" / "read-bases string") of
// the mpileup record. The fiddly encoding rules are mirrored verbatim from
// upstream `bam_plcmd.c::pileup_seq`.
func writeBasesColumn(bw *bufio.Writer, evs []pileupEvent, ref byte, opts MpileupOptions) {
	sortEvents(evs)
	upRef := upper(ref)
	for _, e := range evs {
		if e.dropped {
			continue
		}
		// The -Q/--min-BQ filter applies to EVERY pileup entry, not just
		// aligned bases: upstream (bam_plcmd.c) tests bam_get_qual[p->qpos]
		// against min_baseQ for del/refskip placeholders too, where qpos is
		// the post-gap base. Our del/refskip events carry exactly that base's
		// quality (accumulateRecordEvents' gapQual = rec.Qual[queryPos]), so a
		// '*'/'>'/'<' flanking a low-quality deletion drops out just as it does
		// upstream. Gating this on pileupEventBase over-retained those
		// placeholders, inflating depth near indels.
		if opts.MinBaseQ > 0 && e.effectiveQual() < opts.MinBaseQ {
			continue
		}
		if e.readStart {
			bw.WriteByte('^')
			// `^` is followed by a single char encoding MAPQ as
			// `mapq + 33`. Clamp to printable range to stay within
			// ASCII 33..126.
			c := int(e.mapq) + 33
			if c > 126 {
				c = 126
			}
			bw.WriteByte(byte(c))
		}
		switch e.kind {
		case pileupEventBase:
			b := upper(e.base)
			// When the originating read had SEQ=* (no sequence), upstream
			// renders an explicit 'N'/'n' regardless of whether the ref
			// is also N — the placeholder is not a match. We carry this
			// through with the noSeq flag.
			if b == upRef && !e.noSeq {
				if e.isReverse {
					bw.WriteByte(',')
				} else {
					bw.WriteByte('.')
				}
			} else {
				if e.isReverse {
					bw.WriteByte(lower(b))
				} else {
					bw.WriteByte(b)
				}
			}
		case pileupEventDel:
			// D ops render as '*' regardless of strand.
			bw.WriteByte('*')
		case pileupEventRefSkip:
			// N ops render as '>' (forward) / '<' (reverse), matching
			// upstream `bam_plcmd.c::pileup_seq` for refskip placeholders.
			if e.isReverse {
				bw.WriteByte('<')
			} else {
				bw.WriteByte('>')
			}
		}
		// Indel annotations attached AFTER this event.
		if e.insAfter != "" {
			bw.WriteByte('+')
			bw.WriteString(strconv.Itoa(len(e.insAfter)))
			if e.isReverse {
				bw.WriteString(strings.ToLower(e.insAfter))
			} else {
				bw.WriteString(strings.ToUpper(e.insAfter))
			}
		}
		if e.delAfter > 0 {
			bw.WriteByte('-')
			bw.WriteString(strconv.Itoa(int(e.delAfter)))
			delSeq := e.delBases
			if delSeq == "" {
				delSeq = strings.Repeat("N", int(e.delAfter))
			}
			if e.isReverse {
				bw.WriteString(strings.ToLower(delSeq))
			} else {
				bw.WriteString(strings.ToUpper(delSeq))
			}
		}
		if e.readEnd {
			bw.WriteByte('$')
		}
	}
}

// writeQualsColumn writes the 6th column (read-base qualities, Phred+33).
// Deletion / refskip placeholders carry the quality of the next base on
// the originating read (or 0 if none), matching upstream's
// `bam_plcmd.c::pileup_seq` rendering.
func writeQualsColumn(bw *bufio.Writer, evs []pileupEvent, opts MpileupOptions) {
	sortEvents(evs)
	for _, e := range evs {
		if e.dropped {
			continue
		}
		// The -Q/--min-BQ filter applies to EVERY pileup entry, not just
		// aligned bases: upstream (bam_plcmd.c) tests bam_get_qual[p->qpos]
		// against min_baseQ for del/refskip placeholders too, where qpos is
		// the post-gap base. Our del/refskip events carry exactly that base's
		// quality (accumulateRecordEvents' gapQual = rec.Qual[queryPos]), so a
		// '*'/'>'/'<' flanking a low-quality deletion drops out just as it does
		// upstream. Gating this on pileupEventBase over-retained those
		// placeholders, inflating depth near indels.
		if opts.MinBaseQ > 0 && e.effectiveQual() < opts.MinBaseQ {
			continue
		}
		bw.WriteByte(e.effectiveQual() + 33)
	}
}

// writeMapqColumn writes the optional MAPQ column (-s) using Phred+33
// encoding (matching upstream).
func writeMapqColumn(bw *bufio.Writer, evs []pileupEvent, opts MpileupOptions) {
	sortEvents(evs)
	for _, e := range evs {
		if e.dropped {
			continue
		}
		// The -Q/--min-BQ filter applies to EVERY pileup entry, not just
		// aligned bases: upstream (bam_plcmd.c) tests bam_get_qual[p->qpos]
		// against min_baseQ for del/refskip placeholders too, where qpos is
		// the post-gap base. Our del/refskip events carry exactly that base's
		// quality (accumulateRecordEvents' gapQual = rec.Qual[queryPos]), so a
		// '*'/'>'/'<' flanking a low-quality deletion drops out just as it does
		// upstream. Gating this on pileupEventBase over-retained those
		// placeholders, inflating depth near indels.
		if opts.MinBaseQ > 0 && e.effectiveQual() < opts.MinBaseQ {
			continue
		}
		c := int(e.mapq) + 33
		if c > 126 {
			c = 126
		}
		bw.WriteByte(byte(c))
	}
}

// writeReadBPColumn writes the optional per-read base-position column (-O):
// comma-separated 1-based read positions. Deletion / ref-skip placeholders
// report the post-gap base's 1-based read position (qpos+1), matching
// upstream's bam_plcmd.c:716.
func writeReadBPColumn(bw *bufio.Writer, evs []pileupEvent, opts MpileupOptions) {
	sortEvents(evs)
	first := true
	for _, e := range evs {
		if e.dropped {
			continue
		}
		// The -Q/--min-BQ filter applies to EVERY pileup entry, not just
		// aligned bases: upstream (bam_plcmd.c) tests bam_get_qual[p->qpos]
		// against min_baseQ for del/refskip placeholders too, where qpos is
		// the post-gap base. Our del/refskip events carry exactly that base's
		// quality (accumulateRecordEvents' gapQual = rec.Qual[queryPos]), so a
		// '*'/'>'/'<' flanking a low-quality deletion drops out just as it does
		// upstream. Gating this on pileupEventBase over-retained those
		// placeholders, inflating depth near indels.
		if opts.MinBaseQ > 0 && e.effectiveQual() < opts.MinBaseQ {
			continue
		}
		if !first {
			bw.WriteByte(',')
		}
		first = false
		bw.WriteString(strconv.Itoa(int(e.readBP)))
	}
}

// sortEvents stably orders events by their originating read index so that
// the multi-read column is rendered in a deterministic order (matches
// upstream's "iterate in pileup-walker order").
//
// A single O(n) ordered check guards the sort: the streaming cursor appends a
// column's events in ascending readIdx already (reads are pushed in coordinate
// == push order and each touches a column at most once), so its columns are
// pre-sorted and skip the comparison-sort entirely. The buffered walk's columns
// are likewise built in read order, so the fast-path applies there too — and
// when they are NOT ordered (a future caller, or a non-monotone build) the
// stable sort still runs, producing the identical ordering it always did. The
// fast-path is therefore output-neutral: it only avoids redundant work on input
// that is already in the exact order the stable sort would leave it.
func sortEvents(evs []pileupEvent) {
	if sortedByReadIdx(evs) {
		return
	}
	sort.SliceStable(evs, func(i, j int) bool {
		return evs[i].readIdx < evs[j].readIdx
	})
}

// sortedByReadIdx reports whether evs is already in non-decreasing readIdx
// order — exactly the order a stable sort by readIdx would leave it, so when it
// holds the sort is a guaranteed no-op and can be skipped.
func sortedByReadIdx(evs []pileupEvent) bool {
	for i := 1; i < len(evs); i++ {
		if evs[i].readIdx < evs[i-1].readIdx {
			return false
		}
	}
	return true
}

// accumulateRecordEvents walks rec's CIGAR and appends per-position events
// to evs[pos-beg0] for every reference position the record covers within
// [beg0, end0). It also threads through readStart/readEnd markers and
// attaches +/- indel annotations.
func accumulateRecordEvents(rec *sam.Record, readIdx, beg0, end0 int, evs [][]pileupEvent,
	opts MpileupOptions, refFA *fasta.RandomAccess, chrom string, contig []byte) {
	accumulateRecordEventsLazy(rec, readIdx, beg0, end0, evs, opts, refFA, chrom, false, contig)
}

// accumulateRecordEventsLazy is accumulateRecordEvents with an explicit lazyQual
// switch. When lazyQual is true every emitted event carries a live reference
// (qrec/qpos) to the originating record's quality array instead of a frozen
// byte, so the -Q filter and quals column reflect the record's quality AS OF
// emit time — the streaming overlap cursor relies on this so a deletion/refskip
// placeholder borrows its post-gap base quality before its later mate's overlap
// tweak has run (matching htslib's emit-time bam_get_qual(p->b)[p->qpos]). The
// geometry (^/$ markers, +/- indel annotations, ref-skip boundaries, which
// base/qpos each event uses) is identical to the eager path; only whether the
// quality is read live differs, so lazyQual=false reproduces the original bytes.
//
// contig, when non-nil, is the whole-chromosome reference (uppercased, the exact
// bytes refFA.Fetch returns) and is sliced for a deletion's "-<len><seq>" ref
// bases instead of issuing a fresh refFA.Fetch per D op — that per-deletion Fetch
// re-decompresses the FASTA bgzf block and dominated the streaming walk. Slicing
// the already-loaded contig is byte-identical (same uppercased source bytes); the
// refFA.Fetch fallback is kept for callers without a loaded contig.
func accumulateRecordEventsLazy(rec *sam.Record, readIdx, beg0, end0 int, evs [][]pileupEvent,
	opts MpileupOptions, refFA *fasta.RandomAccess, chrom string, lazyQual bool, contig []byte) {
	if rec.Pos <= 0 {
		return
	}
	refPos := int(rec.Pos) - 1 // 0-based
	queryPos := 0              // 0-based index into rec.Seq
	isReverse := rec.Flag&sam.FlagReverse != 0
	mapq := rec.MapQ
	noSeq := rec.Seq == ""

	// First, identify the very first and very last "outputtable" event for
	// this record (M/=/X bases, D/N gap placeholders). We make a single
	// CIGAR pass and collect (pos0, kind) tuples for the in-window events;
	// then we know which one is first and which is last so we can attach
	// ^/$ markers. The scratch slice is pooled (tmpEventPool) so this CIGAR
	// pass — run once per read across the whole contig in the streaming cursor —
	// no longer heap-allocates a fresh slice per read (it was a top mallocgc
	// source on the streaming hot path).
	tmpp := tmpEventPool.Get().(*[]tmpEvent)
	tmp := (*tmpp)[:0]
	defer func() { *tmpp = tmp[:0]; tmpEventPool.Put(tmpp) }()

	pendingInsTarget := -1 // index in tmp to attach a pending insertion to
	pendingDelTarget := -1 // index in tmp to attach a pending deletion to

	// Walk CIGAR. We collect events for every reference-consuming op so
	// we can later compute first/last markers. Out-of-window events are
	// still recorded (with col == -1) so the first/last markers fall in
	// the right place even when the record spills outside the window.
	for _, op := range rec.Cigar {
		l := int(op.Length())
		o := op.Op()
		switch o {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			for k := 0; k < l; k++ {
				p := refPos + k
				q := queryPos + k
				col := -1
				if p >= beg0 && p < end0 {
					col = p - beg0
				}
				var base byte = 'N'
				if q < len(rec.Seq) {
					base = rec.Seq[q]
				}
				var qual byte
				if q < len(rec.Qual) {
					qual = rec.Qual[q]
				}
				tmp = append(tmp, tmpEvent{
					pos0:   p,
					col:    col,
					kind:   pileupEventBase,
					base:   base,
					qual:   qual,
					qpos:   q,
					readBP: q + 1,
				})
				// If we had a pending insertion or deletion attached to
				// the prior event, fire it now (we've already attached at
				// queue time, so nothing to do here).
				_ = pendingInsTarget
				_ = pendingDelTarget
			}
			refPos += l
			queryPos += l
		case sam.CigarInsertion:
			// Attach a "+<len><seq>" annotation to the most recent
			// outputtable event. If there is none (insertion at the very
			// start of the read), drop the annotation — matches upstream's
			// "starting with I isn't handled" comment in mp_I.sam.
			ins := ""
			if queryPos+l <= len(rec.Seq) {
				ins = rec.Seq[queryPos : queryPos+l]
			}
			if len(tmp) > 0 {
				tmp[len(tmp)-1].insAfter = ins
			}
			queryPos += l
		case sam.CigarDeletion:
			// Each deleted reference base becomes a pileupEventDel.
			// We also attach a "-<len><seq>" annotation to the previous
			// outputtable event, matching upstream.
			//
			// Upstream attaches the quality of the *next* M/=/X base
			// after the deletion to each '*' placeholder; if no next
			// base exists the placeholder gets Phred 0 ('!').
			var delBases string
			if contig != nil && refPos >= 0 && refPos+l <= len(contig) {
				// Slice the already-loaded, uppercased contig — identical bytes to
				// refFA.Fetch but with no per-deletion FASTA re-decompression.
				delBases = string(contig[refPos : refPos+l])
			} else if refFA != nil {
				if b, err := refFA.Fetch(chrom, int64(refPos), int64(refPos+l)); err == nil {
					delBases = string(b)
				}
			}
			if len(tmp) > 0 {
				tmp[len(tmp)-1].delAfter = l
				tmp[len(tmp)-1].delBases = delBases
			}
			// Look ahead in the CIGAR to find the qual of the next M-base
			// on this read; default to 0 if none.
			var gapQual byte
			if queryPos < len(rec.Qual) {
				gapQual = rec.Qual[queryPos]
			}
			for k := 0; k < l; k++ {
				p := refPos + k
				col := -1
				if p >= beg0 && p < end0 {
					col = p - beg0
				}
				tmp = append(tmp, tmpEvent{
					pos0: p,
					col:  col,
					kind: pileupEventDel,
					qual: gapQual,
					qpos: queryPos,
					// Upstream's per-read-position column (-O) prints p->qpos+1
					// for a deletion placeholder (bam_plcmd.c:716), where qpos is
					// the post-gap base index — the same base whose quality the
					// '*' borrows. Match it rather than emitting 0.
					readBP: queryPos + 1,
				})
			}
			refPos += l
		case sam.CigarSkipped:
			// N — reference skip; render '>'/'<' at each position. Like
			// deletions, the qual carried alongside the placeholder is
			// the qual of the next M-base, or 0 if none.
			var gapQual byte
			if queryPos < len(rec.Qual) {
				gapQual = rec.Qual[queryPos]
			}
			for k := 0; k < l; k++ {
				p := refPos + k
				col := -1
				if p >= beg0 && p < end0 {
					col = p - beg0
				}
				tmp = append(tmp, tmpEvent{
					pos0: p,
					col:  col,
					kind: pileupEventRefSkip,
					qual: gapQual,
					qpos: queryPos,
					// As for deletions, upstream's -O column prints p->qpos+1 for
					// a ref-skip placeholder (bam_plcmd.c:716).
					readBP: queryPos + 1,
				})
			}
			refPos += l
		case sam.CigarSoftClip:
			queryPos += l
		case sam.CigarHardClip, sam.CigarPadding:
			// Neither consumes ref nor query (in our accounting).
		}
	}
	if len(tmp) == 0 {
		return
	}

	// Mark the positions bordering each ref-skip (CIGAR N) run as ref-skip
	// boundaries. Upstream's consensus pileup engine sets p->ref_skip on the
	// position immediately before an N (consensus_pileup.c:251-260, when the
	// next op is CREF_SKIP) and on the first position after an N
	// (lines 239-244, `if (p->eof && p->base != '.')`); its Gap5/bayesian
	// caller then excludes those positions from the consensus depth
	// (bam_consensus.c:1333). Crucially the upstream test is `p->base != '.'`
	// — ANY non-ref-skip position, so a deletion ('*') or pad directly
	// abutting the N is flagged too, not only an aligned base. tmp is in
	// reference order, so any non-RefSkip entry adjacent to a RefSkip entry
	// is exactly such a boundary. (The simple caller and the displayed
	// pileup column ignore the flag; only the bayesian depth uses it.)
	for i := range tmp {
		if tmp[i].kind != pileupEventRefSkip {
			continue
		}
		if i > 0 && tmp[i-1].kind != pileupEventRefSkip {
			tmp[i-1].refSkipBoundary = true
		}
		if i+1 < len(tmp) && tmp[i+1].kind != pileupEventRefSkip {
			tmp[i+1].refSkipBoundary = true
		}
	}

	// Find first and last in-window indices so ^/$ markers go on a visible
	// event (we follow upstream: ^ goes on the very first emitted event for
	// the read; if the leading events are out of window we still mark the
	// first visible one to communicate "read starts here-ish" — upstream
	// does the same when seeking with a region).
	firstVisible := -1
	lastVisible := -1
	for i := range tmp {
		if tmp[i].col >= 0 {
			if firstVisible == -1 {
				firstVisible = i
			}
			lastVisible = i
		}
	}
	if firstVisible == -1 {
		return
	}

	// Append the in-window tmp events to the per-position event slices.
	for i := firstVisible; i <= lastVisible; i++ {
		t := tmp[i]
		if t.col < 0 {
			continue
		}
		ev := pileupEvent{
			kind:            t.kind,
			readIdx:         int32(readIdx),
			base:            t.base,
			qual:            t.qual,
			mapq:            mapq,
			readBP:          int32(t.readBP),
			isReverse:       isReverse,
			noSeq:           noSeq,
			readStart:       i == firstVisible && tmp[firstVisible].pos0 == int(rec.Pos)-1,
			readEnd:         i == lastVisible && lastVisible == len(tmp)-1,
			insAfter:        t.insAfter,
			delAfter:        int32(t.delAfter),
			delBases:        t.delBases,
			refSkipBoundary: t.refSkipBoundary,
		}
		if lazyQual {
			// Read this event's quality live from the record at emit time, so
			// mate-overlap tweaks applied after the event was built (but before
			// its column emits) are reflected — and, crucially, tweaks applied
			// AFTER its column emits are not. See pileupEvent.qrec.
			ev.qrec = rec
			ev.qpos = int32(t.qpos)
		}
		evs[t.col] = append(evs[t.col], ev)
	}
}

// upper / lower are tiny byte ASCII case helpers (we want zero-alloc).
func upper(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 32
	}
	return b
}
func lower(b byte) byte {
	if b >= 'A' && b <= 'Z' {
		return b + 32
	}
	return b
}

// mpileupNopReader returned for safety in zero-input edge cases — exported
// for tests where we want to ensure a nil reader can't crash the pipeline.
var mpileupNopReader = strings.NewReader("")

// Ensure the io import is exercised in this file (the linter would
// otherwise complain when the only uses are conditional). The empty stub
// below stays compiled.
var _ io.Reader = mpileupNopReader
