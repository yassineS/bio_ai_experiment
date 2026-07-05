package samtools

import (
	"container/heap"
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// CoverageOptions configures Coverage. Defaults match `samtools coverage`'s
// tabular mode: one row per reference (or per region when Regions is set).
type CoverageOptions struct {
	// Regions restricts output to the given chr[:start-end] specs. When
	// empty, every @SQ entry is reported.
	Regions []string
	// MinMAPQ skips records with MAPQ < MinMAPQ.
	MinMAPQ uint8
	// MinBaseQ skips bases with quality < MinBaseQ from depth/baseq stats.
	MinBaseQ uint8
	// IncludeFlags requires every flag bit in the mask.
	IncludeFlags uint16
	// ExcludeFlags drops records with any of these flag bits. Defaults to
	// FlagUnmapped|FlagSecondary|FlagSupplementary|FlagDuplicate when zero.
	ExcludeFlags uint16
	// MinReadLen mirrors upstream's `-l`: records whose query length is
	// shorter than MinReadLen are skipped (and do not count toward
	// n_selected_reads). Zero disables the filter.
	MinReadLen int
	// MinDepth mirrors upstream's `--min-depth` (default 1): a position with
	// depth below MinDepth is not counted as covered. A zero value is treated
	// as 1, matching upstream which only overrides the default for atoi(opt)>0.
	MinDepth int
	// Histogram requests the ASCII/UTF histogram output mode. Upstream enables
	// this for any of -m / -A / -D / -w.
	Histogram bool
	// PlotDepth selects upstream's `-D` "plot depth" mode (histogram of summed
	// depth per bin instead of breadth of coverage). Only meaningful with
	// Histogram.
	PlotDepth bool
	// ASCII selects upstream's `-A` mode: render the histogram with the
	// two-character ASCII ramp (".", ":") rather than the eight UTF-8 block
	// glyphs. Only meaningful with Histogram.
	ASCII bool
	// NBins is the histogram column count (upstream `-w`). Zero means "use the
	// default 40 columns" (upstream derives this from the terminal width minus
	// 40, falling back to 40 when there is no TTY; the Go port has no TTY when
	// writing to a pipe/file so 40 is the byte-parity default).
	NBins int
	// HeaderOff (-H) suppresses the column-header line in tabular mode.
	HeaderOff bool
}

// CoverageRow is one line of `samtools coverage`'s tabular output.
type CoverageRow struct {
	RName     string
	StartPos  int32
	EndPos    int32
	NumReads  uint64
	CovBases  uint64
	Coverage  float64
	MeanDepth float64
	MeanBaseQ float64
	MeanMapQ  float64
}

// coverageRefState is the per-reference accumulator used during the
// streaming pass.
type coverageRefState struct {
	length int32
	// nReads counts every record with this tid (before MAPQ/length/flag
	// filters), mirroring upstream's stats_aux_t.n_reads.
	nReads uint64
	// nSelectedReads counts records that passed every filter, mirroring
	// upstream's stats_aux_t.n_selected_reads.
	nSelectedReads uint64
	// posDelta maps 0-based reference position -> delta of in-flight depth.
	posDelta map[int]int
	// baseQSumAtPos / baseQCntAtPos accumulate, per 0-based reference
	// position, the sum and count of qualifying base qualities. Upstream only
	// folds these into the mean baseQ at positions counted as covered (depth
	// >= mindepth), so they are kept per-position and summed in
	// computeCoverageStats. They are only allocated when perPosBaseQ is set
	// (mindepth > 1 or a region restriction) — otherwise the whole-contig,
	// mindepth<=1 case uses the global sums below, since every base then sits at
	// a covered position. Keeping per-position maps only when required avoids an
	// O(reference length) blow-up (a whole-chromosome OOM on deep BAMs).
	baseQSumAtPos map[int]uint64
	baseQCntAtPos map[int]uint64
	// baseQSumGlobal / baseQCntGlobal accumulate the same quantities over ALL
	// qualifying bases; used directly when perPosBaseQ is false.
	baseQSumGlobal uint64
	baseQCntGlobal uint64
	// perPosBaseQ records whether the per-position baseQ maps are populated.
	perPosBaseQ bool
	mapQSum     uint64
	// regionRestricted, regBeg and regEnd mirror upstream's BAI-iterator
	// behaviour under `-r`: only reads OVERLAPPING the requested window
	// [regBeg, regEnd) (0-based, half-open) are returned by sam_itr_next, so
	// only they count toward n_reads / n_selected_reads / summed_mapQ. When
	// regionRestricted is false (whole-reference tabular/histogram output) every
	// record with this tid is counted, matching upstream's plain sam_read1 loop.
	regionRestricted bool
	regBeg, regEnd   int
	// streamed reports that the covered-base / summed-coverage aggregates were
	// produced by the streaming interval sweep (the common tabular path) rather
	// than by resolving posDelta. When set, posDelta is never populated and
	// streamCovBases / streamSummedCov below carry the totals.
	streamed        bool
	streamCovBases  uint64
	streamSummedCov uint64
}

// defaultCoverageNBins is the histogram column count used when no terminal
// width is available (the Go port writes to a pipe/file, so upstream's
// terminal-width probe yields its 40-column fallback). Upstream: coverage.c
// `opt_n_bins = 40`.
const defaultCoverageNBins = 40

// Coverage walks records from in (a single BAM/SAM stream — multi-file
// support layered by the CLI) and emits either the tabular summary (one row
// per @SQ entry / region) or, when opts.Histogram is set, the per-reference
// ASCII/UTF histogram.
func Coverage(in io.Reader, w io.Writer, opts CoverageOptions) error {
	r, err := sam.NewReader(in)
	if err != nil {
		return err
	}
	hdr := r.Header()

	excl := opts.ExcludeFlags
	if excl == 0 {
		excl = sam.FlagUnmapped | sam.FlagSecondary | sam.FlagSupplementary | sam.FlagDuplicate
	}
	minDepth := opts.MinDepth
	if minDepth <= 0 {
		minDepth = 1
	}

	// Per-position baseQ bookkeeping is only needed when a base can sit at an
	// UNcovered position (mindepth > 1) or the output is region-restricted; the
	// common whole-contig, mindepth<=1 case uses the global sums instead.
	perPosBaseQ := minDepth > 1 || len(opts.Regions) > 0
	// The common tabular path (no histogram, mindepth<=1, whole references)
	// streams a sorted interval sweep instead of materialising a per-position
	// delta map for the entire contig — the whole-contig posDelta map was the
	// last O(reference length) structure and dominated resident memory on deep,
	// coordinate-sorted BAMs. The other output modes keep the delta map because
	// they need per-position depth (histogram) or per-position baseQ (mindepth>1
	// / region restriction) that the streaming totals do not retain.
	streaming := !opts.Histogram && minDepth <= 1 && len(opts.Regions) == 0
	states := make([]*coverageRefState, len(hdr.Refs))
	for i, ref := range hdr.Refs {
		st := &coverageRefState{
			length:      ref.Length,
			perPosBaseQ: perPosBaseQ,
			streamed:    streaming,
		}
		if !streaming {
			st.posDelta = map[int]int{}
		}
		if perPosBaseQ {
			st.baseQSumAtPos = map[int]uint64{}
			st.baseQCntAtPos = map[int]uint64{}
		}
		states[i] = st
	}

	// Resolve regions list — when empty, every @SQ counts as the full
	// reference. Resolved up-front (before the read loop) so region windows can
	// be stamped onto the per-reference states: under `-r` upstream fetches
	// reads through a BAI iterator, so only reads overlapping the window are
	// counted, and we reproduce that by gating the read counters below.
	type spec struct {
		name     string
		idx      int
		startPos int32
		endPos   int32
	}
	var specs []spec
	if len(opts.Regions) == 0 {
		for i, ref := range hdr.Refs {
			specs = append(specs, spec{name: ref.Name, idx: i, startPos: 1, endPos: ref.Length})
		}
	} else {
		for _, rr := range opts.Regions {
			rg, perr := region.ParseRegion(rr)
			if perr != nil {
				return perr
			}
			i := hdr.RefIndex(rg.Chrom)
			if i < 0 {
				return fmt.Errorf("samtools coverage: unknown ref %q", rg.Chrom)
			}
			startPos := int32(1)
			endPos := hdr.Refs[i].Length
			if rg.Beg > 0 {
				startPos = int32(rg.Beg)
			}
			if rg.End > 0 {
				endPos = int32(rg.End)
			}
			specs = append(specs, spec{name: rg.Chrom, idx: i, startPos: startPos, endPos: endPos})
			// Stamp the window onto the state so the read loop only counts reads
			// overlapping it. Upstream supports a single -r region; when several
			// specs share a reference the last window wins (an over-fetch our
			// multi-region extension accepts — upstream never hits this case).
			states[i].regionRestricted = true
			states[i].regBeg = int(startPos) - 1
			states[i].regEnd = int(endPos)
		}
	}

	// Streaming interval sweep state (common path only). sw is the sweep for
	// the currently active reference; activeIdx is its tid, and finalised marks
	// references whose sweep has already been drained so a coordinate-sort
	// violation (a reference re-appearing after we moved past it, or a record
	// preceding the sweep cursor) is caught rather than silently miscounted.
	var sw *coverageSweep
	activeIdx := -1
	var finalised []bool
	if streaming {
		finalised = make([]bool, len(states))
	}

	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		idx := -1
		if rec.RName != "" && rec.RName != "*" {
			idx = hdr.RefIndex(rec.RName)
		}
		if idx < 0 {
			continue
		}
		st := states[idx]
		// Under `-r`, upstream's BAI iterator only yields reads overlapping the
		// requested window, so reads falling entirely outside it are counted
		// toward NOTHING (neither n_reads nor n_selected_reads nor the depth
		// stats). Depth/baseQ are additionally clipped to the window at summary
		// time, so gating them here is not strictly required, but the read
		// counters (n_reads / n_selected_reads / summed_mapQ) are NOT window-
		// clipped later and must be restricted here to match upstream.
		if st.regionRestricted {
			recStart := int(rec.Pos) - 1
			recEnd := recStart + cigarRefLen(rec.Cigar)
			if recStart >= st.regEnd || recEnd <= st.regBeg {
				continue // no overlap with the region window
			}
		}
		// n_reads counts every record with a valid tid, before any filter
		// (upstream read_bam increments n_reads first).
		st.nReads++
		if rec.Flag&excl != 0 {
			continue
		}
		if opts.IncludeFlags != 0 && rec.Flag&opts.IncludeFlags != opts.IncludeFlags {
			continue
		}
		if rec.MapQ < opts.MinMAPQ {
			continue
		}
		if opts.MinReadLen > 0 && cigarQueryLen(rec.Cigar) < opts.MinReadLen {
			continue
		}
		st.nSelectedReads++
		st.mapQSum += uint64(rec.MapQ)

		if streaming {
			recStart := int(rec.Pos) - 1
			if idx != activeIdx {
				if sw != nil {
					sw.advanceTo(math.MaxInt)
					states[activeIdx].streamCovBases = sw.covBases
					states[activeIdx].streamSummedCov = sw.summedCov
					finalised[activeIdx] = true
				}
				if finalised[idx] {
					return fmt.Errorf("samtools coverage: input is not coordinate-sorted (reference %q re-appears); sort the input first", rec.RName)
				}
				sw = newCoverageSweep(0, int(st.length), minDepth)
				activeIdx = idx
			}
			if recStart < sw.cursor {
				return fmt.Errorf("samtools coverage: input is not coordinate-sorted (record at %s:%d precedes cursor); sort the input first", rec.RName, recStart+1)
			}
			// Finalise every segment strictly before this record's leftmost
			// reference position: coordinate order guarantees no later record
			// (or run) contributes there.
			sw.advanceTo(recStart)
		}

		// Walk CIGAR, accumulating per-base depth and baseQ stats. Bases
		// below MinBaseQ contribute neither to depth nor baseq sums
		// (matches upstream).
		refPos := int(rec.Pos) - 1
		queryPos := 0
		for _, op := range rec.Cigar {
			length := int(op.Length())
			switch op.Op() {
			case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
				if opts.MinBaseQ == 0 {
					// No base-quality filter: every base of the run counts toward
					// depth, so record the depth delta ONCE for the whole run
					// instead of a posDelta entry per reference position (the
					// telescoping intermediate entries are net-zero — same sweep,
					// O(1) instead of O(length)).
					if streaming {
						sw.pushInterval(refPos, refPos+length)
					} else {
						st.posDelta[refPos]++
						st.posDelta[refPos+length]--
					}
					for i := 0; i < length; i++ {
						var q uint8 = 0xff
						if queryPos+i < len(rec.Qual) {
							q = rec.Qual[queryPos+i]
						}
						if q != 0xff {
							st.baseQSumGlobal += uint64(q)
							st.baseQCntGlobal++
							if st.perPosBaseQ {
								st.baseQSumAtPos[refPos+i] += uint64(q)
								st.baseQCntAtPos[refPos+i]++
							}
						}
					}
				} else {
					for i := 0; i < length; i++ {
						var q uint8 = 0xff
						if queryPos+i < len(rec.Qual) {
							q = rec.Qual[queryPos+i]
						}
						if q < opts.MinBaseQ {
							continue
						}
						pos := refPos + i
						if q != 0xff {
							st.baseQSumGlobal += uint64(q)
							st.baseQCntGlobal++
							if st.perPosBaseQ {
								st.baseQSumAtPos[pos] += uint64(q)
								st.baseQCntAtPos[pos]++
							}
						}
						if streaming {
							sw.pushInterval(pos, pos+1)
						} else {
							st.posDelta[pos]++
							st.posDelta[pos+1]--
						}
					}
				}
				refPos += length
				queryPos += length
			case sam.CigarInsertion, sam.CigarSoftClip:
				queryPos += length
			case sam.CigarDeletion, sam.CigarSkipped:
				refPos += length
			case sam.CigarHardClip, sam.CigarPadding:
				// no movement on query or ref
			}
		}
	}

	// Drain the last active reference's sweep into its totals.
	if streaming && sw != nil {
		sw.advanceTo(math.MaxInt)
		states[activeIdx].streamCovBases = sw.covBases
		states[activeIdx].streamSummedCov = sw.summedCov
	}

	if opts.Histogram {
		// Upstream only prints a histogram for references the pileup actually
		// visits — i.e. those with at least one selected read producing a
		// pileup column. References with no selected reads are silently
		// skipped (unlike the tabular form, which lists every @SQ). A
		// region-restricted run (-r) always renders its single requested
		// reference, matching upstream which seeds the iterator's tid.
		regionRestricted := len(opts.Regions) > 0
		first := true
		for _, s := range specs {
			if !regionRestricted && states[s.idx].nSelectedReads == 0 {
				continue
			}
			if !first {
				fmt.Fprint(w, "\n")
			}
			first = false
			printCoverageHist(w, s.name, s.startPos, s.endPos, states[s.idx], opts, minDepth)
		}
		return nil
	}

	if !opts.HeaderOff {
		fmt.Fprintln(w, "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq")
	}
	// Emission order matches upstream coverage.c exactly. In whole-reference
	// (non-region) tabular mode upstream prints every reference the pileup
	// visits FIRST — in the order the coordinate-sorted pileup encounters them,
	// which is header/@SQ order for sorted input (print_tabular_line is called
	// on each tid transition, coverage.c:597/676) — and only then walks the
	// remaining, never-covered references in header order (the trailing
	// `if (!opt_reg)` loop, coverage.c:681-688). A reference is "covered" iff at
	// least one selected read reached the pileup (nSelectedReads > 0), the same
	// criterion the histogram path uses. Region-restricted output (`-r`) keeps
	// the requested order: upstream never runs its reorder loop under opt_reg,
	// so specs are emitted as given. Only the O(#references) summary rows are
	// buffered/reordered here; the per-position streaming accumulators are
	// untouched, preserving the bounded-memory RSS fix.
	order := make([]int, 0, len(specs))
	if len(opts.Regions) == 0 {
		for i, s := range specs {
			if states[s.idx].nSelectedReads > 0 {
				order = append(order, i)
			}
		}
		for i, s := range specs {
			if states[s.idx].nSelectedReads == 0 {
				order = append(order, i)
			}
		}
	} else {
		for i := range specs {
			order = append(order, i)
		}
	}
	for _, oi := range order {
		s := specs[oi]
		row := summariseCoverage(s.name, s.startPos, s.endPos, states[s.idx], minDepth)
		// Match upstream coverage.c print_tabular_line exactly: %g for the
		// coverage percentage and mean depth, %.3g for mean baseQ / mapQ.
		fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
			row.RName, row.StartPos, row.EndPos, row.NumReads, row.CovBases,
			formatGShortest(row.Coverage), formatGShortest(row.MeanDepth),
			formatG(row.MeanBaseQ, 3), formatG(row.MeanMapQ, 3))
	}
	return nil
}

// coverageSweep computes the covered-base and summed-coverage aggregates for a
// single reference by sweeping half-open depth intervals in coordinate order,
// without ever materialising a per-position map. It reproduces exactly the
// segment arithmetic of computeCoverageStats (same [lo,hi) clip and depth >=
// minDepth gate), fed incrementally: pushInterval adds a read's covered run(s),
// and advanceTo resolves every segment strictly before a frontier once the
// coordinate-sorted stream guarantees nothing earlier remains.
type coverageSweep struct {
	lo, hi    int // half-open clip window [lo, hi) (0-based)
	minDepth  int
	cursor    int // next reference position not yet accounted
	curDepth  int
	starts    covStartHeap // pending intervals whose start >= cursor
	ends      covEndHeap   // end positions of currently-open intervals
	covBases  uint64
	summedCov uint64
}

// newCoverageSweep returns a sweep clipped to the half-open window [lo, hi).
func newCoverageSweep(lo, hi, minDepth int) *coverageSweep {
	return &coverageSweep{lo: lo, hi: hi, minDepth: minDepth, cursor: lo}
}

// pushInterval queues the half-open covered interval [start, end). The caller
// guarantees start >= cursor (coordinate order), so the interval waits in the
// start heap until the sweep reaches it.
func (s *coverageSweep) pushInterval(start, end int) {
	if end <= start {
		return
	}
	heap.Push(&s.starts, covInterval{start: start, end: end})
}

// advanceTo resolves all segments and events strictly before frontier, folding
// each covered segment (depth >= minDepth) into covBases / summedCov with the
// same [lo,hi) clipping computeCoverageStats applies. Passing math.MaxInt drains
// every remaining interval.
func (s *coverageSweep) advanceTo(frontier int) {
	for {
		next := frontier
		if s.starts.Len() > 0 && s.starts[0].start < next {
			next = s.starts[0].start
		}
		if s.ends.Len() > 0 && s.ends[0] < next {
			next = s.ends[0]
		}
		if s.curDepth >= s.minDepth {
			from, to := s.cursor, next
			if from < s.lo {
				from = s.lo
			}
			if to > s.hi {
				to = s.hi
			}
			if to > from {
				n := uint64(to - from)
				s.covBases += n
				s.summedCov += uint64(s.curDepth) * n
			}
		}
		s.cursor = next
		if next == frontier {
			return
		}
		// Apply every event coincident with `next`. Start/end ordering at a
		// shared position does not affect any segment or the final depth.
		for s.starts.Len() > 0 && s.starts[0].start == next {
			iv := heap.Pop(&s.starts).(covInterval)
			heap.Push(&s.ends, iv.end)
			s.curDepth++
		}
		for s.ends.Len() > 0 && s.ends[0] == next {
			heap.Pop(&s.ends)
			s.curDepth--
		}
	}
}

// covInterval is a half-open covered interval awaiting the sweep cursor.
type covInterval struct{ start, end int }

// covStartHeap is a min-heap of pending intervals ordered by start position.
type covStartHeap []covInterval

func (h covStartHeap) Len() int           { return len(h) }
func (h covStartHeap) Less(i, j int) bool { return h[i].start < h[j].start }
func (h covStartHeap) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }
func (h *covStartHeap) Push(x any)        { *h = append(*h, x.(covInterval)) }
func (h *covStartHeap) Pop() any {
	old := *h
	n := len(old)
	v := old[n-1]
	*h = old[:n-1]
	return v
}

// covEndHeap is a min-heap of open-interval end positions.
type covEndHeap []int

func (h covEndHeap) Len() int           { return len(h) }
func (h covEndHeap) Less(i, j int) bool { return h[i] < h[j] }
func (h covEndHeap) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }
func (h *covEndHeap) Push(x any)        { *h = append(*h, x.(int)) }
func (h *covEndHeap) Pop() any {
	old := *h
	n := len(old)
	v := old[n-1]
	*h = old[:n-1]
	return v
}

// cigarQueryLen returns the number of query bases consumed by cig (M/I/S/=/X),
// matching htslib's bam_cigar2qlen used by coverage's -l filter.
func cigarQueryLen(cig []sam.CigarOp) int {
	n := 0
	for _, op := range cig {
		switch op.Op() {
		case sam.CigarMatch, sam.CigarInsertion, sam.CigarSoftClip, sam.CigarEqual, sam.CigarMismatch:
			n += int(op.Length())
		}
	}
	return n
}

// cigarRefLen returns the number of reference bases spanned by cig
// (M/D/N/=/X), matching htslib's bam_cigar2rlen. It gives the read's reference
// end offset used to test overlap with a `-r` region window.
func cigarRefLen(cig []sam.CigarOp) int {
	n := 0
	for _, op := range cig {
		switch op.Op() {
		case sam.CigarMatch, sam.CigarDeletion, sam.CigarSkipped, sam.CigarEqual, sam.CigarMismatch:
			n += int(op.Length())
		}
	}
	return n
}

// coverageStats holds the aggregate per-region statistics shared by the
// tabular and histogram renderers. covDepth maps each covered 0-based
// position within [start-1, end) to its (filtered) depth. baseQSum / baseQCnt
// accumulate the per-base quality sums at covered positions only, matching
// upstream which folds baseQ into the mean only where depth >= mindepth.
type coverageStats struct {
	covBases       uint64
	summedCoverage uint64
	baseQSum       uint64
	baseQCnt       uint64
	covDepth       map[int]int
}

// walkDepths resolves the running delta map into a per-position depth map for
// the half-open 0-based window [start-1, end), keeping only positions whose
// depth is >= minDepth.
func walkDepths(start, end int32, st *coverageRefState, minDepth int) map[int]int {
	out := map[int]int{}
	if len(st.posDelta) == 0 {
		return out
	}
	keys := make([]int, 0, len(st.posDelta))
	for k := range st.posDelta {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	lo := int(start) - 1
	hi := int(end)
	depth := 0
	prev := keys[0]
	for _, k := range keys {
		if depth >= minDepth && k > prev {
			from := prev
			to := k
			if from < lo {
				from = lo
			}
			if to > hi {
				to = hi
			}
			for p := from; p < to; p++ {
				out[p] = depth
			}
		}
		depth += st.posDelta[k]
		prev = k
	}
	return out
}

// computeCoverageStats resolves the covered-base / summed-coverage / baseQ
// totals for the half-open 0-based window [start-1, end) by SWEEPING the depth
// delta intervals — it never materialises a per-position depth map for the
// tabular row (the O(reference length) blow-up). Positions with depth < minDepth
// are not counted as covered, mirroring upstream's `depth >= mindepth` gate.
// baseQ is taken from the global sums when the per-position maps were not kept
// (the whole-contig, mindepth<=1 case — every base then sits at a covered
// position), else summed per covered position. When wantCovDepth is set (the
// histogram renderer) the per-position depth map is materialised as before.
func computeCoverageStats(start, end int32, st *coverageRefState, minDepth int, wantCovDepth bool) coverageStats {
	var cs coverageStats
	if wantCovDepth {
		cs.covDepth = map[int]int{}
	}
	lo := int(start) - 1
	hi := int(end)
	if len(st.posDelta) > 0 {
		keys := make([]int, 0, len(st.posDelta))
		for k := range st.posDelta {
			keys = append(keys, k)
		}
		sort.Ints(keys)
		depth := 0
		prev := keys[0]
		for _, k := range keys {
			if depth >= minDepth && k > prev {
				from, to := prev, k
				if from < lo {
					from = lo
				}
				if to > hi {
					to = hi
				}
				if to > from {
					n := uint64(to - from)
					cs.covBases += n
					cs.summedCoverage += uint64(depth) * n
					if st.perPosBaseQ {
						for p := from; p < to; p++ {
							cs.baseQSum += st.baseQSumAtPos[p]
							cs.baseQCnt += st.baseQCntAtPos[p]
						}
					}
					if wantCovDepth {
						for p := from; p < to; p++ {
							cs.covDepth[p] = depth
						}
					}
				}
			}
			depth += st.posDelta[k]
			prev = k
		}
	}
	if !st.perPosBaseQ {
		cs.baseQSum = st.baseQSumGlobal
		cs.baseQCnt = st.baseQCntGlobal
	}
	return cs
}

// summariseCoverage computes the per-region tabular row from a refState.
func summariseCoverage(name string, start, end int32, st *coverageRefState, minDepth int) CoverageRow {
	row := CoverageRow{
		RName:    name,
		StartPos: start,
		EndPos:   end,
		NumReads: st.nSelectedReads,
	}
	if st.nSelectedReads > 0 {
		row.MeanMapQ = float64(st.mapQSum) / float64(st.nSelectedReads)
	}
	var cs coverageStats
	if st.streamed {
		// Common path: covered-base / summed-coverage totals were produced by
		// the streaming interval sweep, and (perPosBaseQ being false there)
		// baseQ comes from the global running sums.
		cs = coverageStats{
			covBases:       st.streamCovBases,
			summedCoverage: st.streamSummedCov,
			baseQSum:       st.baseQSumGlobal,
			baseQCnt:       st.baseQCntGlobal,
		}
	} else {
		cs = computeCoverageStats(start, end, st, minDepth, false)
	}
	if cs.baseQCnt > 0 {
		row.MeanBaseQ = float64(cs.baseQSum) / float64(cs.baseQCnt)
	}
	row.CovBases = cs.covBases
	regionLen := float64(end - (start - 1))
	if regionLen > 0 {
		row.Coverage = 100.0 * float64(cs.covBases) / regionLen
		row.MeanDepth = float64(cs.summedCoverage) / regionLen
	}
	return row
}

// blockChars8 are the eight UTF-8 block glyphs (LOWER ONE EIGHTH BLOCK …
// FULL BLOCK) used by upstream's default histogram. blockChars2 is the ASCII
// fallback selected by -A.
var (
	blockChars8 = []string{"▁", "▂", "▃", "▄", "▅", "▆", "▇", "█"}
	blockChars2 = []string{".", ":"}
)

const verticalLine = "│" // BOX DRAWINGS LIGHT VERTICAL

// printCoverageHist renders the per-reference histogram, mirroring upstream
// coverage.c print_hist. start/end are 1-based inclusive region bounds.
func printCoverageHist(w io.Writer, name string, start, end int32, st *coverageRefState, opts CoverageOptions, minDepth int) {
	beg := int64(start) - 1 // 0-based, half-open with end
	endPos := int64(end)
	regionLen := float64(endPos - beg)

	nBins := opts.NBins
	if nBins <= 0 {
		nBins = defaultCoverageNBins
	}
	if int64(nBins) > endPos-beg {
		nBins = int(endPos - beg)
	}
	if nBins < 1 {
		nBins = 1
	}
	binWidth := (endPos - beg) / int64(nBins)
	if binWidth < 1 {
		binWidth = 1
	}

	cs := computeCoverageStats(start, end, st, minDepth, true)

	// Fill the per-bin histogram counts. Under -D (plot depth) upstream sums
	// the depth at every position with depth >= 1 — BEFORE the mindepth gate
	// (coverage.c:648) — so use a depth-1 walk there. The default breadth-of-
	// coverage histogram counts only positions counted as covered (depth >=
	// mindepth), so it reuses the covered set.
	hist := make([]uint64, nBins)
	if opts.PlotDepth {
		for pos, depth := range walkDepths(start, end, st, 1) {
			bin := (int64(pos) - beg) / binWidth
			if bin < 0 || bin >= int64(nBins) {
				continue
			}
			hist[bin] += uint64(depth)
		}
	} else {
		for pos := range cs.covDepth {
			bin := (int64(pos) - beg) / binWidth
			if bin < 0 || bin >= int64(nBins) {
				continue
			}
			hist[bin]++
		}
	}

	full := !opts.ASCII
	blockChars := blockChars8
	blockLen := 8
	if !full {
		blockChars = blockChars2
		blockLen = 2
	}
	const nRows = 10

	histData := make([]float64, nBins)
	maxVal := 0.0
	scale := 100.0
	if opts.PlotDepth {
		scale = 1.0
	}
	for i := 0; i < nBins; i++ {
		histData[i] = scale * float64(hist[i]) / float64(binWidth)
		if histData[i] > maxVal {
			maxVal = histData[i]
		}
	}

	fmt.Fprintf(w, "%s (%sbp)\n", name, readableBPs(float64(st.length)))

	rowBinSize := maxVal / float64(nRows)
	for i := nRows - 1; i >= 0; i-- {
		currentBin := rowBinSize * float64(i)
		if opts.PlotDepth {
			fmt.Fprintf(w, ">%8.1f ", float64(i)*rowBinSize)
		} else {
			fmt.Fprintf(w, ">%7.2f%% ", currentBin)
		}
		if full {
			fmt.Fprint(w, verticalLine)
		} else {
			fmt.Fprint(w, "|")
		}
		for col := 0; col < nBins; col++ {
			// Upstream: int cur_val_diff = round(blockchar_len * (hist_data[col]
			// - current_bin) / row_bin_size) - 1; (coverage.c:256). When the
			// histogram is empty (max_val == 0, so row_bin_size == 0) this is
			// (int)(round(NaN) - 1) == (int)NaN — platform-dependent undefined
			// behaviour: 0 on ARM64 (FCVTZS), so the bar shows blockChars[0]
			// (▁); INT_MIN on x86-64 (CVTTSD2SI), which is < 0, so it shows a
			// space. Replicate upstream's cast with a float64->int32 conversion
			// (same hardware instruction as C) so the bar matches the upstream
			// binary byte-for-byte on either platform rather than pinning one.
			curValDiff := int(int32(math.Round(float64(blockLen)*(histData[col]-currentBin)/rowBinSize) - 1))
			if curValDiff < 0 {
				fmt.Fprint(w, " ")
			} else {
				if curValDiff >= blockLen {
					curValDiff = blockLen - 1
				}
				fmt.Fprint(w, blockChars[curValDiff])
			}
		}
		if full {
			fmt.Fprint(w, verticalLine)
		} else {
			fmt.Fprint(w, "|")
		}
		fmt.Fprint(w, " ")
		switch i {
		case 9:
			fmt.Fprintf(w, "Number of reads: %d", st.nSelectedReads)
		case 8:
			if st.nReads-st.nSelectedReads > 0 {
				fmt.Fprintf(w, "    (%d filtered)", st.nReads-st.nSelectedReads)
			}
		case 7:
			fmt.Fprintf(w, "Covered bases:   %sbp", readableBPs(float64(cs.covBases)))
		case 6:
			fmt.Fprintf(w, "Percent covered: %s%%", formatG(100.0*float64(cs.covBases)/regionLen, 4))
		case 5:
			fmt.Fprintf(w, "Mean coverage:   %sx", formatG(float64(cs.summedCoverage)/regionLen, 3))
		case 4:
			mbq := 0.0
			if cs.baseQCnt > 0 {
				mbq = float64(cs.baseQSum) / float64(cs.baseQCnt)
			}
			fmt.Fprintf(w, "Mean baseQ:      %s", formatG(mbq, 3))
		case 3:
			mmq := 0.0
			if st.nSelectedReads > 0 {
				mmq = float64(st.mapQSum) / float64(st.nSelectedReads)
			}
			fmt.Fprintf(w, "Mean mapQ:       %s", formatG(mmq, 3))
		case 1:
			fmt.Fprintf(w, "Histo bin width: %sbp", readableBPs(float64(binWidth)))
		case 0:
			if opts.PlotDepth {
				fmt.Fprintf(w, "Histo max cov:   %s", formatG(maxVal, 5))
			} else {
				fmt.Fprintf(w, "Histo max bin:   %s%%", formatG(maxVal, 5))
			}
		}
		fmt.Fprint(w, "\n")
	}

	// x-axis labels, centered in 10-char fields. Upstream prints a label at
	// beg+1, then every tenth bin, then the end coordinate.
	fmt.Fprintf(w, "     %s", centerText(readableBPs(float64(beg+1)), 10))
	for rest := 10; rest < 10*(nBins/10); rest += 10 {
		fmt.Fprint(w, centerText(readableBPs(float64(beg+binWidth*int64(rest))), 10))
	}
	lastPadding := nBins % 10
	fmt.Fprintf(w, "%*s%s", lastPadding, " ", centerText(readableBPs(float64(endPos)), 10))
	fmt.Fprint(w, "\n")
}

// readableBPs renders a base-pair count with K/M/G/T units, mirroring
// upstream coverage.c readable_bps (the decimal precision equals the chosen
// unit index, so "1000" -> "1K", "1500000" -> "1.50M").
func readableBPs(bp float64) string {
	units := []string{"", "K", "M", "G", "T"}
	i := 0
	for bp >= 1000 && i < len(units)-1 {
		bp /= 1000
		i++
	}
	return strconv.FormatFloat(bp, 'f', i, 64) + units[i]
}

// centerText centers text in a field of the given width, prefixing a single
// leading space, mirroring upstream coverage.c center_text. When the text is
// at least as wide as the field it is returned unpadded.
func centerText(text string, width int) string {
	l := len(text)
	if l > width {
		// Upstream asserts len <= width; defensively return the text as-is.
		return text
	}
	padding := (width - l) / 2
	paddingEx := (width - l) % 2
	if padding >= 1 {
		// " %*s%*s": a leading space, then text right-justified in l+padding,
		// then (padding-1+paddingEx) trailing spaces.
		return " " + fmt.Sprintf("%*s%*s", l+padding, text, padding-1+paddingEx, " ")
	}
	return text
}

// formatG renders v using C's "%.<prec>g" semantics (shortest of %e/%f with
// the given number of significant digits, trailing zeros stripped), matching
// the histogram side-panel statistics printed by upstream coverage.c.
func formatG(v float64, prec int) string {
	return strconv.FormatFloat(v, 'g', prec, 64)
}

// formatGShortest renders v using C's bare "%g" semantics, which defaults to
// 6 significant digits. Used by the tabular coverage / mean-depth columns.
func formatGShortest(v float64) string {
	return strconv.FormatFloat(v, 'g', 6, 64)
}
