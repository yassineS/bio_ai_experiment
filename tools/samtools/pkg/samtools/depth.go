package samtools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"runtime"
	"sort"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// ReadDecodeThreads resolves the BGZF inflate worker count for a streaming
// scan. A caller-supplied positive thread count (CLI -@) is honoured verbatim;
// 0 — the default when no -@ is given — opts into parallel inflate across the
// machine's cores (capped at 8 to avoid oversubscribing on many-core hosts).
// Only BGZF block inflation is parallelised, so the decoded record stream — and
// therefore every tool's output — is byte-identical for any thread count.
func ReadDecodeThreads(n int) int {
	if n > 0 {
		return n
	}
	c := runtime.NumCPU()
	if c > 8 {
		c = 8
	}
	if c < 1 {
		c = 1
	}
	return c
}

// DefaultDepthExcludeFlags matches upstream samtools depth's default
// filter-out flag list, UNMAP,SECONDARY,QCFAIL,DUP (see
// reference_code/samtools/bam2depth.c: .flag = BAM_FUNMAP | BAM_FSECONDARY
// | BAM_FDUP | BAM_FQCFAIL).
const DefaultDepthExcludeFlags uint16 = sam.FlagUnmapped | sam.FlagSecondary | sam.FlagQCFail | sam.FlagDuplicate

// DepthOptions configures the behaviour of Depth.
type DepthOptions struct {
	// AllPositions, when true, emits positions with zero depth that fall
	// inside the regions actually covered by at least one read (matches
	// upstream `-a/--all`).
	AllPositions bool
	// AllTransPositions, when true, emits every position of every reference
	// in the header (matches upstream `-A/--all-trans`).
	AllTransPositions bool
	// Regions are samtools-style "chr:start-end" specifiers (CLI `-r`).
	Regions []string
	// BedPath, when non-empty, restricts emitted positions to the union of
	// these BED intervals (CLI `-b`).
	BedPath string
	// MinMAPQ skips reads with MAPQ < this value. Upstream exposes this as
	// `-Q`/`--min-MQ` (bam2depth.c opt.min_mqual).
	MinMAPQ uint8
	// MinBaseQ skips bases with quality < this value. Upstream exposes this
	// as `-q`/`--min-BQ` (bam2depth.c opt.min_qual).
	MinBaseQ uint8
	// MinReadLen skips reads shorter than this (after CIGAR query-length)
	// (CLI `-l`).
	MinReadLen int
	// IncludeFlags requires ALL these flag bits to be set (CLI `-f`).
	IncludeFlags uint16
	// ExcludeFlags drops reads with ANY of these flag bits set (CLI `-F`,
	// default 0x4 unmapped).
	ExcludeFlags uint16
	// MaxDepth, when > 0, caps reported depth (CLI `-d`).
	MaxDepth int
	// Threads is accepted for compatibility; v1 is single-threaded.
	Threads int
}

// depthRegion is the resolved set of intervals (per reference, 0-based half-
// open) we will report on. nil means "every position of every reference"
// (i.e. `-A`).
type depthRegion struct {
	// byRef maps a reference name to the sorted, non-overlapping list of
	// half-open 0-based [beg, end) intervals to include. A nil entry means
	// "all positions on this reference".
	byRef map[string][][2]int
	// names is the ordered set of reference names to consider, drawn from
	// the input header order. When AllTrans is true we report every
	// reference; otherwise we report only references that appear in either
	// the region list or the BED.
	names []string
}

// Depth runs the depth computation across one or more BAM/SAM inputs. It
// emits one line per emitted position to out:
//
//	chrom\tpos\tdepth1[\tdepth2 ...]\n
//
// Where pos is 1-based and depth_k is the integer depth in the k-th input
// (parallel positional iteration across inputs).
//
// All inputs must share an identical @SQ ordering for the output to be
// well-defined; this is the same constraint as upstream samtools.
//
// The implementation streams: it mirrors upstream bam2depth.c by holding only
// a sliding ring buffer of per-input depth/span counters spanning the active
// read window (O(active reads), not O(positions)). Each output line is written
// to the buffered writer as soon as its position is flushed, so memory stays
// bounded even for chromosome-wide, coordinate-sorted inputs — independent of
// how many positions are emitted.
func Depth(inputs []io.Reader, out io.Writer, opts DepthOptions) error {
	if len(inputs) == 0 {
		return fmt.Errorf("samtools depth: no inputs")
	}
	readers := make([]sam.Reader, len(inputs))
	for i, r := range inputs {
		rd, err := alnio.NewReaderThreaded(r, "", ReadDecodeThreads(opts.Threads))
		if err != nil {
			return fmt.Errorf("samtools depth: input %d: %w", i, err)
		}
		if rc, ok := rd.(io.Closer); ok {
			defer rc.Close()
		}
		readers[i] = rd
	}
	return depthFromReaders(readers, out, opts)
}

// depthFromReaders runs the streaming depth engine over already-constructed
// per-input sam.Readers. It is shared by Depth (linear, whole-file readers) and
// DepthFile (the indexed seek-and-scan readers), so the aggregation, region
// clamping and output are byte-identical regardless of how the records were
// sourced.
func depthFromReaders(readers []sam.Reader, out io.Writer, opts DepthOptions) error {
	if len(readers) == 0 {
		return fmt.Errorf("samtools depth: no inputs")
	}
	hdr := readers[0].Header()
	for i := 1; i < len(readers); i++ {
		if !sameRefOrder(hdr, readers[i].Header()) {
			return fmt.Errorf("samtools depth: input %d has a different @SQ ordering than input 0", i)
		}
	}

	region, err := resolveDepthRegion(opts, hdr)
	if err != nil {
		return err
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()

	// Wrap each input in a one-record-ahead peeking source that applies the
	// depth-level filters and tags every record with its header reference
	// index, so the streaming engine can merge inputs in coordinate order
	// without buffering the whole stream.
	srcs := make([]*depthSource, len(readers))
	for i, rd := range readers {
		srcs[i] = newDepthSource(rd, hdr, opts)
	}

	st := newDepthStream(bw, hdr, region, opts, len(readers))
	return st.run(srcs)
}

// depthSource is a one-record-lookahead view over a single input. It applies
// the per-read filters (flag include/exclude, MAPQ, min read length) as it
// reads, so dropped reads never reach the streaming engine, and resolves each
// kept record's reference index once for the coordinate-order merge.
//
// It decodes each record in place into a single reused buffer (rec) via the
// reader's depth-tailored decode, which skips the read name, mate reference,
// SEQ expansion and aux stream — none of which depth reads — and the QUAL block
// unless a base-quality filter is active. This removes the dominant per-record
// decode and allocation cost (aux parsing alone is ~a third of the BAM decode
// on real data) while leaving the emitted depth byte-identical.
type depthSource struct {
	rd      sam.Reader
	hdr     *sam.Header
	opts    DepthOptions
	cur     *sam.Record // points at rec when a kept record is buffered, else nil
	rec     sam.Record  // reused decode target for the next record
	tid     int         // header reference index of cur
	needQ   bool        // decode QUAL (only when a base-quality floor is set)
	lastRef string      // cache key for lastTid
	lastTid int         // header index resolved for lastRef (avoids re-scanning)
	err     error
	done    bool

	// readDepth is the reader's depth-tailored decode when it exposes one,
	// else nil (SAM / CRAM fall back to readInto / Read).
	readDepth func(*sam.Record, bool) error
	readInto  func(*sam.Record) error
}

// newDepthSource builds a depthSource and primes its first kept record.
func newDepthSource(rd sam.Reader, hdr *sam.Header, opts DepthOptions) *depthSource {
	s := &depthSource{rd: rd, hdr: hdr, opts: opts, lastTid: -1}
	s.needQ = opts.MinBaseQ > 0
	if rd, ok := rd.(interface {
		ReadDepthInto(*sam.Record, bool) error
	}); ok {
		s.readDepth = rd.ReadDepthInto
	}
	if ri, ok := rd.(interface{ ReadInto(*sam.Record) error }); ok {
		s.readInto = ri.ReadInto
	}
	s.advance()
	return s
}

// resolveTid maps a record's reference name to its header index, caching the
// last (name → index) pair. Input is coordinate-sorted, so a reference's reads
// are contiguous and consecutive records overwhelmingly share a name, turning
// the per-record header scan into a single string compare in the common case.
func (s *depthSource) resolveTid(name string) int {
	if name == s.lastRef && s.lastTid >= 0 {
		return s.lastTid
	}
	tid := s.hdr.RefIndex(name)
	s.lastRef, s.lastTid = name, tid
	return tid
}

// advance loads the next kept record into s.cur (or marks the source done /
// failed). Filtered-out records are skipped here so the engine only ever sees
// records that contribute to depth, mirroring the filter loop in
// bam2depth.c's fastdepth_core.
func (s *depthSource) advance() {
	for {
		var err error
		switch {
		case s.readDepth != nil:
			err = s.readDepth(&s.rec, s.needQ)
		case s.readInto != nil:
			err = s.readInto(&s.rec)
		default:
			var rec *sam.Record
			rec, err = s.rd.Read()
			if err == nil {
				s.rec = *rec
			}
		}
		if err == io.EOF {
			s.cur, s.done = nil, true
			return
		}
		if err != nil {
			s.cur, s.err, s.done = nil, err, true
			return
		}
		if !keepDepthRecord(&s.rec, s.opts) {
			continue
		}
		tid := s.resolveTid(s.rec.RName)
		if tid < 0 {
			// Reference not in the header ordering; skip (matches upstream's
			// silent handling of records on unknown references).
			continue
		}
		s.cur, s.tid = &s.rec, tid
		return
	}
}

// depthStream is the streaming aggregator. It holds a per-input ring buffer of
// depth and span counters spanning the active read window for the reference
// currently being processed, exactly mirroring upstream bam2depth.c's
// depth_hist: positions are flushed (and lines written) as reads advance, so
// only O(active reads) state is retained, never O(positions).
type depthStream struct {
	bw     *bufio.Writer
	hdr    *sam.Header
	region depthRegion
	opts   DepthOptions
	n      int // number of inputs

	// Ring buffers, one per input, indexed by (pos & mask). depth holds the
	// base-quality-filtered count; span holds the number of reads physically
	// spanning the position (M/=/X/D/N), used to emit interior depth-0 rows.
	size  int
	mask  int
	depth [][]int32
	span  [][]int32
	// endPos[i] is one past the furthest reference position any read from
	// input i has reached on the current reference (absolute coordinate).
	endPos []int

	lastOutput int // next absolute position not yet emitted on the current ref
	curRef     int // header index of the reference being processed, or -1
	refName    string

	// Per-reference interval cursor over the merged include intervals. When
	// nil the whole reference is included (sentinel for `-A` / default).
	intervals [][2]int
	ivIdx     int

	line []byte // reused output line buffer
}

// newDepthStream constructs the streaming aggregator.
func newDepthStream(bw *bufio.Writer, hdr *sam.Header, region depthRegion, opts DepthOptions, n int) *depthStream {
	return &depthStream{
		bw:         bw,
		hdr:        hdr,
		region:     region,
		opts:       opts,
		n:          n,
		curRef:     -1,
		lastOutput: 0,
		endPos:     make([]int, n),
	}
}

// run merges the per-input sources in coordinate order and feeds each record
// to the aggregator, then flushes the trailing reference.
func (st *depthStream) run(srcs []*depthSource) error {
	for {
		// Pick the next record across all inputs by (tid, pos), matching the
		// merge in bam2depth.c's main loop.
		best := -1
		bestTid, bestPos := 0, 0
		for i, s := range srcs {
			if s.err != nil {
				return s.err
			}
			if s.cur == nil {
				continue
			}
			pos := int(s.cur.Pos) - 1
			if best < 0 || s.tid < bestTid || (s.tid == bestTid && pos < bestPos) {
				best, bestTid, bestPos = i, s.tid, pos
			}
		}
		if best < 0 {
			break // all inputs exhausted
		}
		if err := st.addRecord(srcs[best].cur, srcs[best].tid, best); err != nil {
			return err
		}
		srcs[best].advance()
	}
	// Flush the final reference (and, under -aa, any references after it).
	return st.finish()
}

// orderIndex returns the position of ref tid within region.names, or -1 if
// the reference is not selected for output.
func (st *depthStream) selected(tid int) bool {
	name := st.hdr.Refs[tid].Name
	_, ok := st.region.byRef[name]
	return ok
}

// addRecord feeds a single coordinate-ordered record (from input `file`, on
// header reference `tid`) into the ring buffer, flushing any positions that
// precede it. Mirrors bam2depth.c add_depth.
//
// Records on references that are not in the output set (e.g. a `-r`/BED
// restriction to other contigs) are dropped without disturbing any ring or
// output state: input is coordinate-sorted, so all of a reference's reads are
// contiguous and an unselected reference never interleaves a selected one.
func (st *depthStream) addRecord(rec *sam.Record, tid, file int) error {
	if !st.selected(tid) {
		return nil
	}
	if tid != st.curRef {
		// Close out the previous reference, then any skipped references that
		// -aa must zero-fill, then open the new one.
		if err := st.closeRef(tid); err != nil {
			return err
		}
		if err := st.openRef(tid, int(rec.Pos)-1); err != nil {
			return err
		}
	} else {
		pos0 := int(rec.Pos) - 1
		if st.lastOutput < pos0 {
			if err := st.flushTo(pos0); err != nil {
				return err
			}
		}
	}
	return st.accumulate(rec, file)
}

// openRef initialises the ring buffer state for a freshly seen reference whose
// first record starts at firstPos0 (0-based), zero-filling the head of the ref
// under -a/-aa.
func (st *depthStream) openRef(tid, firstPos0 int) error {
	st.curRef = tid
	st.refName = st.hdr.Refs[tid].Name
	for i := range st.endPos {
		st.endPos[i] = 0
	}
	// Reset ring contents lazily: we re-zero positions on demand as reads
	// extend the window, so a fresh ref simply restarts lastOutput.
	st.resetRing()

	st.loadIntervals(tid)
	begClamp := st.refBeg()
	st.lastOutput = firstPos0
	if begClamp >= 0 && st.lastOutput < begClamp {
		st.lastOutput = begClamp
	}

	if st.includeZeros() {
		// Zero-fill the start of the reference up to the first read.
		if err := st.zeroRegion(0, firstPos0); err != nil {
			return err
		}
	}
	return nil
}

// closeRef flushes the tail of the current reference (positions still covered
// by a read's span past lastOutput) and, under -a/-aa, the zero-depth tail to
// the reference end. nextTid is the reference about to be opened (or
// len(Refs) when finishing); under -aa any wholly skipped references between
// the two are zero-filled.
func (st *depthStream) closeRef(nextTid int) error {
	if st.curRef >= 0 {
		// Flush positions that remain inside some input's read span.
		i := st.lastOutput
		for {
			covered := false
			for f := 0; f < st.n; f++ {
				if i < st.endPos[f] {
					covered = true
					break
				}
			}
			if !covered {
				break
			}
			if err := st.emitPos(i); err != nil {
				return err
			}
			i++
		}
		st.lastOutput = i
		if st.includeZeros() {
			refLen := int(st.hdr.Refs[st.curRef].Length)
			if err := st.zeroRegion(i, refLen); err != nil {
				return err
			}
		}
	}

	// Under -aa (AllTrans) without a region restriction, zero-fill any
	// references wholly skipped between curRef and nextTid.
	if st.opts.AllTransPositions && len(st.opts.Regions) == 0 && st.opts.BedPath == "" {
		from := 0
		if st.curRef >= 0 {
			from = st.curRef + 1
		}
		for r := from; r < nextTid && r < len(st.hdr.Refs); r++ {
			st.curRef = r
			st.refName = st.hdr.Refs[r].Name
			st.loadIntervals(r)
			st.lastOutput = 0
			refLen := int(st.hdr.Refs[r].Length)
			if err := st.zeroRegion(0, refLen); err != nil {
				return err
			}
		}
	}
	return nil
}

// finish flushes the final reference and, under -aa, any trailing references.
func (st *depthStream) finish() error {
	if err := st.closeRef(len(st.hdr.Refs)); err != nil {
		return err
	}
	// -a/-aa with a region but no reads at all: zero-fill the region. Matches
	// bam2depth.c's "-a or -aa without a single read being output yet" branch.
	if st.curRef < 0 && st.includeZeros() && (len(st.opts.Regions) > 0 || st.opts.BedPath != "" || st.opts.AllTransPositions) {
		for _, name := range st.region.names {
			tid := st.hdr.RefIndex(name)
			if tid < 0 {
				continue
			}
			st.curRef = tid
			st.refName = name
			st.loadIntervals(tid)
			st.lastOutput = 0
			refLen := int(st.hdr.Refs[tid].Length)
			if err := st.zeroRegion(0, refLen); err != nil {
				return err
			}
		}
	}
	return st.bw.Flush()
}

// accumulate records one read's depth/span contribution into the ring buffer,
// growing the ring if the read's span exceeds the current capacity.
func (st *depthStream) accumulate(rec *sam.Record, file int) error {
	if rec.Pos <= 0 {
		return nil
	}
	pos0 := int(rec.Pos) - 1
	endPos := pos0
	if n := rec.Cigar.ReferenceLength(); n > 0 {
		endPos = pos0 + n // 0-based, one past end
	}
	// Clip the read end to the region end if one is set.
	if hi := st.refEnd(); hi >= 0 && endPos > hi {
		endPos = hi
	}

	st.ensureCapacity(endPos - pos0)

	// Zero any newly seen ring slots between the old endPos[file] and this
	// read's end so accumulation starts from a clean count (upstream zeroes
	// the same window before incrementing).
	from := st.endPos[file]
	if pos0 > from {
		from = pos0
	}
	for i := from; i < endPos; i++ {
		st.depth[file][i&st.mask] = 0
		st.span[file][i&st.mask] = 0
	}

	st.addReadDepthRing(rec, file, pos0, endPos)

	if st.endPos[file] < endPos {
		st.endPos[file] = endPos
	}
	return nil
}

// addReadDepthRing walks the CIGAR and increments the ring buffer in place,
// applying the per-base quality filter. It mirrors addReadDepth's original
// CIGAR accounting (M/=/X add filtered depth and span; D/N add span only;
// I/S consume query only; H/P ignored) but writes directly into the sliding
// ring instead of a per-interval difference array.
func (st *depthStream) addReadDepthRing(rec *sam.Record, file, beg0, regionEnd int) {
	depth := st.depth[file]
	span := st.span[file]
	mask := st.mask
	hasQual := st.opts.MinBaseQ > 0 && len(rec.Qual) > 0
	minQ := st.opts.MinBaseQ

	refPos := int(rec.Pos) - 1
	queryPos := 0
	for _, op := range rec.Cigar {
		l := int(op.Length())
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			runBeg := refPos
			runEnd := refPos + l
			lo, hi := clip(runBeg, runEnd, regionEnd)
			for p := lo; p < hi; p++ {
				span[p&mask]++
				if hasQual {
					qIdx := queryPos + (p - runBeg)
					if qIdx >= 0 && qIdx < len(rec.Qual) && rec.Qual[qIdx] < minQ {
						continue
					}
				}
				depth[p&mask]++
			}
			refPos += l
			queryPos += l
		case sam.CigarInsertion, sam.CigarSoftClip:
			queryPos += l
		case sam.CigarDeletion, sam.CigarSkipped:
			// Deletions and reference skips consume reference and extend the
			// read's covered span (they print as depth 0) but add no depth.
			runBeg := refPos
			runEnd := refPos + l
			lo, hi := clip(runBeg, runEnd, regionEnd)
			for p := lo; p < hi; p++ {
				span[p&mask]++
			}
			refPos += l
		case sam.CigarHardClip, sam.CigarPadding:
			// Neither consumes ref nor query in our accounting.
		}
		if refPos >= regionEnd {
			return
		}
	}
}

// clip intersects [runBeg, runEnd) with [0, regionEnd) (the region end is
// already clamped into regionEnd by the caller; the low edge is implicitly 0
// because positions before a read's start cannot be touched here).
func clip(runBeg, runEnd, regionEnd int) (int, int) {
	lo := runBeg
	hi := runEnd
	if hi > regionEnd {
		hi = regionEnd
	}
	return lo, hi
}

// flushTo emits every position in [lastOutput, upto) that is still covered by
// some input's read span, zero-filling holes under -a/-aa, then advances
// lastOutput to upto. Mirrors the in-ref flush loop of bam2depth.c add_depth.
func (st *depthStream) flushTo(upto int) error {
	i := st.lastOutput
	for ; i < upto; i++ {
		covered := false
		for f := 0; f < st.n; f++ {
			if i < st.endPos[f] {
				covered = true
				break
			}
		}
		if !covered {
			break
		}
		if err := st.emitPos(i); err != nil {
			return err
		}
	}
	if st.includeZeros() && i < upto {
		if err := st.zeroRegion(i, upto); err != nil {
			return err
		}
	}
	st.lastOutput = upto
	return nil
}

// emitPos writes the depth line for a single covered position (0-based) if it
// falls inside an include interval, reading each input's count from the ring.
func (st *depthStream) emitPos(pos0 int) error {
	if !st.inInterval(pos0) {
		return nil
	}
	st.line = append(st.line[:0], st.refName...)
	st.line = append(st.line, '\t')
	st.line = strconv.AppendInt(st.line, int64(pos0+1), 10)
	mask := st.mask
	for f := 0; f < st.n; f++ {
		var d int32
		if pos0 < st.endPos[f] {
			d = st.depth[f][pos0&mask]
		}
		if st.opts.MaxDepth > 0 && d > int32(st.opts.MaxDepth) {
			d = int32(st.opts.MaxDepth)
		}
		st.line = append(st.line, '\t')
		st.line = strconv.AppendInt(st.line, int64(d), 10)
	}
	st.line = append(st.line, '\n')
	_, err := st.bw.Write(st.line)
	return err
}

// zeroRegion emits depth-0 rows for [start, end) (0-based, half-open) clamped
// to the active region, for positions inside an include interval. Mirrors
// bam2depth.c zero_region.
func (st *depthStream) zeroRegion(start, end int) error {
	if begClamp := st.refBeg(); begClamp >= 0 && start < begClamp {
		start = begClamp
	}
	if endClamp := st.refEnd(); endClamp >= 0 && end > endClamp {
		end = endClamp
	}
	for i := start; i < end; i++ {
		if !st.inInterval(i) {
			continue
		}
		st.line = append(st.line[:0], st.refName...)
		st.line = append(st.line, '\t')
		st.line = strconv.AppendInt(st.line, int64(i+1), 10)
		for f := 0; f < st.n; f++ {
			st.line = append(st.line, '\t', '0')
		}
		st.line = append(st.line, '\n')
		if _, err := st.bw.Write(st.line); err != nil {
			return err
		}
	}
	return nil
}

// includeZeros reports whether zero-depth positions inside covered regions are
// emitted (-a or -aa).
func (st *depthStream) includeZeros() bool {
	return st.opts.AllPositions || st.opts.AllTransPositions
}

// ensureCapacity grows the per-input ring buffers so they can hold a window of
// `need` positions, mirroring the geometric growth in bam2depth.c. Growing
// preserves the live window [lastOutput, lastOutput+oldSize) by re-keying.
func (st *depthStream) ensureCapacity(need int) {
	if need+1 < st.size && st.size > 0 {
		return
	}
	old := st.size
	newSize := st.size
	if newSize == 0 {
		newSize = 2048
	}
	for need+1 >= newSize {
		newSize *= 2
	}
	newMask := newSize - 1
	for f := 0; f < st.n; f++ {
		nd := make([]int32, newSize)
		ns := make([]int32, newSize)
		if old > 0 {
			for i := st.lastOutput; i < st.lastOutput+old; i++ {
				nd[i&newMask] = st.depth[f][i&st.mask]
				ns[i&newMask] = st.span[f][i&st.mask]
			}
		}
		st.depth[f] = nd
		st.span[f] = ns
	}
	st.size = newSize
	st.mask = newMask
}

// resetRing allocates the per-input ring slices on first use (they are reused
// verbatim across references; positions are re-zeroed lazily as reads extend
// the window, so no per-reference clearing is needed).
func (st *depthStream) resetRing() {
	if st.depth == nil {
		st.depth = make([][]int32, st.n)
		st.span = make([][]int32, st.n)
	}
}

// loadIntervals selects the merged include-interval set for reference tid and
// resets the interval cursor. A nil interval set means "whole reference".
func (st *depthStream) loadIntervals(tid int) {
	name := st.hdr.Refs[tid].Name
	ivs := st.region.byRef[name]
	if ivs == nil {
		st.intervals = nil
	} else {
		st.intervals = mergeIntervals(ivs)
	}
	st.ivIdx = 0
}

// refBeg returns the lower clamp (0-based) for the active reference when a
// single contiguous region restriction is in force, or -1 for no clamp. The
// per-position interval test (inInterval) enforces the actual membership; this
// clamp only mirrors upstream's dh->beg fast path for the common single-region
// case so zero-fill loops start at the right place.
func (st *depthStream) refBeg() int {
	if len(st.intervals) == 1 {
		return st.intervals[0][0]
	}
	return -1
}

// refEnd returns the upper clamp (0-based, exclusive) for the active reference
// under a single contiguous region restriction, or -1 for no clamp.
func (st *depthStream) refEnd() int {
	if len(st.intervals) == 1 {
		return st.intervals[0][1]
	}
	return -1
}

// inInterval reports whether 0-based position pos0 falls inside an include
// interval for the active reference. A nil interval set includes everything.
// The cursor st.ivIdx advances monotonically because positions are emitted in
// ascending order.
func (st *depthStream) inInterval(pos0 int) bool {
	if st.intervals == nil {
		return true
	}
	for st.ivIdx < len(st.intervals) && pos0 >= st.intervals[st.ivIdx][1] {
		st.ivIdx++
	}
	if st.ivIdx >= len(st.intervals) {
		return false
	}
	return pos0 >= st.intervals[st.ivIdx][0]
}

// keepDepthRecord applies the depth-level filters: flag include/exclude,
// MAPQ, minimum read-length on query bases.
func keepDepthRecord(rec *sam.Record, opts DepthOptions) bool {
	if rec.Pos <= 0 || rec.RName == "" {
		return false
	}
	if opts.IncludeFlags != 0 && rec.Flag&opts.IncludeFlags != opts.IncludeFlags {
		return false
	}
	if opts.ExcludeFlags != 0 && rec.Flag&opts.ExcludeFlags != 0 {
		return false
	}
	if opts.MinMAPQ > 0 && rec.MapQ < opts.MinMAPQ {
		return false
	}
	if opts.MinReadLen > 0 && rec.Cigar.QueryLength() < opts.MinReadLen {
		return false
	}
	return true
}

// resolveDepthRegion produces the set of [chrom, beg0, end0) intervals we
// will emit depth for, based on opts.Regions, opts.BedPath, opts.AllTrans,
// and the input header.
func resolveDepthRegion(opts DepthOptions, hdr *sam.Header) (depthRegion, error) {
	out := depthRegion{byRef: map[string][][2]int{}}
	// Build a header order index for stable output.
	orderIdx := map[string]int{}
	for i, r := range hdr.Refs {
		orderIdx[r.Name] = i
	}

	add := func(chrom string, beg0, end0 int) {
		if _, ok := orderIdx[chrom]; !ok {
			// Skip unknown chromosomes silently — matches upstream.
			return
		}
		out.byRef[chrom] = append(out.byRef[chrom], [2]int{beg0, end0})
	}

	switch {
	case opts.AllTransPositions:
		for _, r := range hdr.Refs {
			out.byRef[r.Name] = nil // sentinel: "all positions"
		}
	case opts.BedPath != "" || len(opts.Regions) > 0:
		if opts.BedPath != "" {
			f, err := os.Open(opts.BedPath)
			if err != nil {
				return out, fmt.Errorf("samtools depth: open BED: %w", err)
			}
			defer f.Close()
			rd := bed.NewReader(f)
			for {
				rec, err := rd.Read()
				if err == io.EOF {
					break
				}
				if err != nil {
					return out, fmt.Errorf("samtools depth: read BED: %w", err)
				}
				add(rec.Chrom, rec.ChromStart, rec.ChromEnd)
			}
		}
		if len(opts.Regions) > 0 {
			resolved, _, rerr := region.ResolveRegions(opts.Regions, func(name string) int {
				return hdr.RefIndex(name)
			})
			if rerr != nil {
				return out, rerr
			}
			for _, r := range resolved {
				end0 := r.End0
				if end0 > 1<<29 {
					// open-ended — clamp to the reference length.
					end0 = int(refLength(hdr, r.Region.Chrom))
				}
				add(r.Region.Chrom, r.Beg0, end0)
			}
		}
	default:
		// "Whatever the reads cover" — emit every reference; the
		// streaming path will skip the zero-depth positions unless `-a`
		// is set.
		for _, r := range hdr.Refs {
			out.byRef[r.Name] = nil
		}
	}

	// Compose ordered names from header order.
	for _, r := range hdr.Refs {
		if _, ok := out.byRef[r.Name]; ok {
			out.names = append(out.names, r.Name)
		}
	}
	return out, nil
}

// refLength returns the @SQ LN length for the named reference, or 0 if it
// is not present.
func refLength(hdr *sam.Header, name string) int {
	for _, r := range hdr.Refs {
		if r.Name == name {
			return int(r.Length)
		}
	}
	return 0
}

// sameRefOrder reports whether two headers list references in the same
// order (matching name and order).
func sameRefOrder(a, b *sam.Header) bool {
	if len(a.Refs) != len(b.Refs) {
		return false
	}
	for i := range a.Refs {
		if a.Refs[i].Name != b.Refs[i].Name {
			return false
		}
	}
	return true
}

// mergeIntervals returns the union of overlapping/adjacent 0-based half-
// open intervals, sorted by start.
func mergeIntervals(in [][2]int) [][2]int {
	if len(in) <= 1 {
		return in
	}
	cp := make([][2]int, len(in))
	copy(cp, in)
	sort.Slice(cp, func(i, j int) bool { return cp[i][0] < cp[j][0] })
	out := cp[:0]
	cur := cp[0]
	for i := 1; i < len(cp); i++ {
		if cp[i][0] <= cur[1] {
			if cp[i][1] > cur[1] {
				cur[1] = cp[i][1]
			}
		} else {
			out = append(out, cur)
			cur = cp[i]
		}
	}
	out = append(out, cur)
	return out
}
