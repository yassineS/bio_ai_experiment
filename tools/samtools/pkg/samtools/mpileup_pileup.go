package samtools

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

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
type pileupEvent struct {
	kind pileupEventKind
	// readIdx is the 0-based index of the originating record in the
	// per-chrom record slice; used to order events stably within a
	// position so multi-read columns match upstream's "in-order"
	// rendering.
	readIdx int
	// base is the read base (uppercase if forward, lowercase if reverse)
	// at this reference position. '*' for pileupEventDel/RefSkip.
	base byte
	// qual is the Phred quality of the read base; 0 for non-base events.
	qual byte
	// mapq is the read's MAPQ.
	mapq uint8
	// readBP is the 1-based position of this base within the read's
	// SEQ (i.e. the original query base index + 1, regardless of strand).
	// 0 for non-base events.
	readBP int
	// isReverse reports the read's reverse-strand flag bit.
	isReverse bool
	// readStart is true when this is the very first event emitted for
	// this read; the output gets a leading "^<mapq>" marker.
	readStart bool
	// readEnd is true when this is the very last event emitted for this
	// read; the output gets a trailing "$".
	readEnd bool
	// insAfter, when non-nil, is the inserted sequence to annotate as
	// "+<len><seq>" immediately after this event's base. Only valid on
	// pileupEventBase.
	insAfter string
	// delAfter, when > 0, is the length of a deletion that begins at the
	// position immediately AFTER this event's base. The reference bases
	// being deleted are recorded in delBases for accurate "-<len><seq>"
	// rendering when a FASTA reference is loaded; otherwise they fall
	// back to 'N' x delAfter, matching upstream.
	delAfter int
	delBases string
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
func emitMpileupWindow(bw *bufio.Writer, chrom string, beg0, end0, refLen int,
	perInputChromRecs [][]*sam.Record, refFA *fasta.RandomAccess, posFilter *positionFilter,
	opts MpileupOptions) error {

	nIn := len(perInputChromRecs)

	// Build a per-input slice of per-position event slices for the window.
	// We accumulate by walking each read's CIGAR once and appending each
	// event to events[pos-beg0]. width can be very large for chrom-wide
	// scans; for typical region queries it is small.
	width := end0 - beg0
	events := make([][][]pileupEvent, nIn)
	for i := 0; i < nIn; i++ {
		events[i] = make([][]pileupEvent, width)
		for ridx, rec := range perInputChromRecs[i] {
			accumulateRecordEvents(rec, ridx, beg0, end0, events[i], opts, refFA, chrom)
		}
	}

	// Optionally remove overlapping mate-pair contributions (-x).
	if opts.IgnoreOverlaps {
		for i := 0; i < nIn; i++ {
			dropOverlapEvents(events[i], perInputChromRecs[i])
		}
	}

	// Optionally fetch a reference slab so each line can carry the right
	// refbase column. Empty when no FASTA was supplied (upstream emits 'N').
	var refSlab []byte
	if refFA != nil {
		// Coerce [beg0, end0) onto the contig length already done by the
		// caller, so this Fetch is always in-range.
		b, err := refFA.Fetch(chrom, int64(beg0), int64(end0))
		if err != nil {
			return fmt.Errorf("samtools mpileup: fetch %s:%d-%d: %w", chrom, beg0, end0, err)
		}
		refSlab = b
	}

	for pos0 := beg0; pos0 < end0; pos0++ {
		col := pos0 - beg0
		// Gather depths per input.
		depths := make([]int, nIn)
		any := false
		for i := 0; i < nIn; i++ {
			d := liveDepth(events[i][col], opts.MinBaseQ)
			if opts.MaxDepth > 0 && d > opts.MaxDepth {
				d = opts.MaxDepth
			}
			depths[i] = d
			if d > 0 {
				any = true
			}
		}

		// Decide whether to emit this position.
		pos1 := pos0 + 1
		if posFilter != nil && !posFilter.contains(chrom, pos1) {
			continue
		}
		if !any {
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

// liveDepth returns the count of "live" (not filtered, base-quality OK)
// events at a position. We count pileupEventDel/RefSkip too because they
// occupy a depth slot in upstream output (they're written as '*').
func liveDepth(evs []pileupEvent, minBQ uint8) int {
	n := 0
	for i := range evs {
		if evs[i].dropped {
			continue
		}
		if evs[i].kind == pileupEventBase {
			if minBQ > 0 && evs[i].qual < minBQ {
				continue
			}
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
		if e.kind == pileupEventBase {
			if opts.MinBaseQ > 0 && e.qual < opts.MinBaseQ {
				continue
			}
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
			bw.WriteString(strconv.Itoa(e.delAfter))
			delSeq := e.delBases
			if delSeq == "" {
				delSeq = strings.Repeat("N", e.delAfter)
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
		if e.kind == pileupEventBase {
			if opts.MinBaseQ > 0 && e.qual < opts.MinBaseQ {
				continue
			}
		}
		bw.WriteByte(e.qual + 33)
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
		if e.kind == pileupEventBase {
			if opts.MinBaseQ > 0 && e.qual < opts.MinBaseQ {
				continue
			}
		}
		c := int(e.mapq) + 33
		if c > 126 {
			c = 126
		}
		bw.WriteByte(byte(c))
	}
}

// writeReadBPColumn writes the optional per-read base-position column (-O):
// comma-separated 1-based read positions. Deletion events get a "0".
func writeReadBPColumn(bw *bufio.Writer, evs []pileupEvent, opts MpileupOptions) {
	sortEvents(evs)
	first := true
	for _, e := range evs {
		if e.dropped {
			continue
		}
		if e.kind == pileupEventBase {
			if opts.MinBaseQ > 0 && e.qual < opts.MinBaseQ {
				continue
			}
		}
		if !first {
			bw.WriteByte(',')
		}
		first = false
		bw.WriteString(strconv.Itoa(e.readBP))
	}
}

// sortEvents stably orders events by their originating read index so that
// the multi-read column is rendered in a deterministic order (matches
// upstream's "iterate in pileup-walker order").
func sortEvents(evs []pileupEvent) {
	sort.SliceStable(evs, func(i, j int) bool {
		return evs[i].readIdx < evs[j].readIdx
	})
}

// accumulateRecordEvents walks rec's CIGAR and appends per-position events
// to evs[pos-beg0] for every reference position the record covers within
// [beg0, end0). It also threads through readStart/readEnd markers and
// attaches +/- indel annotations.
func accumulateRecordEvents(rec *sam.Record, readIdx, beg0, end0 int, evs [][]pileupEvent,
	opts MpileupOptions, refFA *fasta.RandomAccess, chrom string) {
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
	// ^/$ markers.
	type tmpEvent struct {
		pos0            int
		col             int // column index into evs (-1 when out of window)
		kind            pileupEventKind
		base            byte
		qual            byte
		readBP          int
		insAfter        string
		delAfter        int
		delBases        string
		refSkipBoundary bool
	}
	tmp := make([]tmpEvent, 0, 8)

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
			if refFA != nil {
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
			readIdx:         readIdx,
			base:            t.base,
			qual:            t.qual,
			mapq:            mapq,
			readBP:          t.readBP,
			isReverse:       isReverse,
			noSeq:           noSeq,
			readStart:       i == firstVisible && tmp[firstVisible].pos0 == int(rec.Pos)-1,
			readEnd:         i == lastVisible && lastVisible == len(tmp)-1,
			insAfter:        t.insAfter,
			delAfter:        t.delAfter,
			delBases:        t.delBases,
			refSkipBoundary: t.refSkipBoundary,
		}
		evs[t.col] = append(evs[t.col], ev)
	}
}

// dropOverlapEvents implements `-x / --ignore-overlaps`: when two records
// with the same QName cover the same reference position (i.e. mate-pair
// overlap), one half's events are masked. The convention matches upstream
// `bam_plp_overlap`: keep the higher-quality base, drop the lower; on
// ties keep the first-seen.
func dropOverlapEvents(events [][]pileupEvent, recs []*sam.Record) {
	if len(recs) < 2 {
		return
	}
	// Build a quick QName -> []recordIdx map so we know which read pairs
	// could overlap.
	byQName := map[string][]int{}
	for i, r := range recs {
		if r.Flag&sam.FlagPaired == 0 {
			continue
		}
		byQName[r.QName] = append(byQName[r.QName], i)
	}
	if len(byQName) == 0 {
		return
	}
	// Resolve the overlap one position at a time.
	for col := range events {
		evs := events[col]
		if len(evs) < 2 {
			continue
		}
		// Group events by QName via their readIdx.
		seenByQName := map[string]int{}
		for i := range evs {
			if evs[i].dropped {
				continue
			}
			qname := recs[evs[i].readIdx].QName
			if recs[evs[i].readIdx].Flag&sam.FlagPaired == 0 {
				continue
			}
			prev, ok := seenByQName[qname]
			if !ok {
				seenByQName[qname] = i
				continue
			}
			// We have two events for the same QName at this position
			// — keep the higher-quality base, drop the other.
			if eventQual(evs[i]) > eventQual(evs[prev]) {
				evs[prev].dropped = true
				seenByQName[qname] = i
			} else {
				evs[i].dropped = true
			}
		}
	}
}

// eventQual returns the qual for ordering purposes; non-base events get
// 0 so a real base wins.
func eventQual(e pileupEvent) byte {
	if e.kind == pileupEventBase {
		return e.qual
	}
	return 0
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
