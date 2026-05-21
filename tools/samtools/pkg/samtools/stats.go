// Package samtools — stats implementation.
//
// Mirrors `samtools stats`. The Summary Numbers (SN) block, the read-length
// sections (RL/FRL/LRL), the MAPQ and IS histograms, the per-cycle quality
// histograms (FFQ/LFQ), the GC-content sections (GCF/GCL), the ACGT-content
// sections (GCC/GCT), the indel sections (IC/ID), the leading CHK checksum
// block and the COV coverage-distribution histogram are all byte-faithful to
// upstream. PARITY_VALIDATION.md tracks which sections are byte-faithful vs.
// deferred.
//
// Skipped intentionally for v1 (documented):
//   - GCD GC-depth distribution (requires reference bases).
//   - OXC oxidation-context counts (requires reference bases).
//   - --target-regions BED restriction.
//   - --remove-overlaps (single-record stats are unaffected by overlap
//     removal for the counters we emit).
package samtools

import (
	"bufio"
	"fmt"
	"hash/crc32"
	"io"
	"math"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// StatsOptions configures the Stats run.
type StatsOptions struct {
	// RefSeq path; accepted but only the optional sections that need a
	// reference are gated on it (currently none in v1).
	RefSeq string
	// Coverage is the raw "MIN,MAX,STEP" string passed via -c; parsed to
	// configure the COV coverage-distribution bins (default "1,1000,1").
	Coverage string
	// RequiredFlag requires ALL bits set on each record (-l/--required-flag).
	RequiredFlag uint16
	// FilteringFlag drops records with ANY bit set (-F/--filtering-flag).
	FilteringFlag uint16
	// MaxDepth is accepted for CLI compatibility; the COV histogram bins
	// every position's depth without an upstream-style depth cap.
	MaxDepth int
	// MinMAPQ skips records below this MAPQ.
	MinMAPQ uint8
	// RemoveDups skips records with the duplicate flag set.
	RemoveDups bool
	// RemoveOverlaps is accepted but unused (no overlap-walk in v1).
	RemoveOverlaps bool
	// MaxInsertSize is the upper bound for the IS section (default 8000).
	MaxInsertSize int
	// Sparse omits empty placeholder section bodies.
	Sparse bool
	// TargetBED restricts stats to regions; accepted but unused.
	TargetBED string
	// Threads is accepted; v1 is single-threaded.
	Threads int
}

// StatsCounters holds the accumulator state for one stats run.
type StatsCounters struct {
	// Identity / sort
	IsSorted int // 1 if observed records are in coordinate order or @HD SO:coordinate
	// internal sort-detect state
	sortBroken bool
	lastRef    string
	lastPos    int32

	// Counters per upstream's bam_stats.c. All counted on the QC-passed
	// branch unless the field name says otherwise.
	RawTotal             int64 // primary records only
	FilteredOut          int64 // dropped by RequiredFlag/FilteringFlag
	Sequences            int64 // RawTotal - FilteredOut
	FirstFrag            int64
	LastFrag             int64
	ReadsMapped          int64
	ReadsMappedAndPaired int64
	ReadsUnmapped        int64
	ReadsProperlyPaired  int64
	ReadsPaired          int64
	ReadsDuplicated      int64
	ReadsMQ0             int64
	ReadsQCFailed        int64
	NonPrimary           int64
	Supplementary        int64
	TotalLength          int64
	TotalFirstFragLength int64
	TotalLastFragLength  int64
	BasesMapped          int64 // ignores clipping (just seq.length for mapped reads)
	BasesMappedCigar     int64 // sums M/=/X bases from cigar
	BasesTrimmed         int64 // soft-clipped bases
	BasesDuplicated      int64
	Mismatches           int64 // sum of NM aux fields
	MaxLength            int64
	MaxFirstFragLength   int64
	MaxLastFragLength    int64
	TotalQual            int64 // sum of (qual byte) over all bases of QC-passed reads (255 for "*")
	TotalQualBases       int64
	InwardPairs          int64
	OutwardPairs         int64
	OtherOrientPairs     int64
	DiffChromosomePairs  int64
	InsertSumAbs         int64 // |TLEN| over pairs accepted into IS
	InsertSumSqAbs       int64 // |TLEN|^2
	InsertPairsCounted   int64

	// histograms keyed on small integer ranges
	RL      map[int64]int64 // read length distribution
	FRL     map[int64]int64
	LRL     map[int64]int64
	MAPQ    [256]int64
	ISInw   map[int64]int64 // insert size: inward
	ISOutw  map[int64]int64 // outward
	ISOther map[int64]int64 // other

	// Per-cycle and base-content accumulators (upstream stats.c). Cycle is
	// 0-based internally and printed 1-based.
	//
	// qualsFirst/qualsLast hold, per cycle, the count of each quality value
	// 0..statsNQuals-1. Indexed [cycle][qual].
	qualsFirst []cycleQuals
	qualsLast  []cycleQuals
	maxQual    int // highest observed quality value
	maxLen1st  int // longest first-fragment unclipped read
	maxLen2nd  int // longest last-fragment unclipped read
	maxLen     int // longest unclipped read overall

	// gcFirst/gcLast are GC-content histograms over statsNGC bins.
	gcFirst [statsNGC]int64
	gcLast  [statsNGC]int64

	// acgtCycles1st/2nd hold as-sequenced ACGT/N/other counts per cycle;
	// acgtRevcomp holds read-oriented counts (reverse reads complemented).
	acgtCycles1st []acgtNoCount
	acgtCycles2nd []acgtNoCount
	acgtRevcomp   []acgtNoCount

	// Indel distributions. insertions/deletions are keyed by indel length-1;
	// ins/del cycle buffers are keyed by cycle index.
	insertions  map[int]int64
	deletions   map[int]int64
	insCycle1st map[int]int64
	insCycle2nd map[int]int64
	delCycle1st map[int]int64
	delCycle2nd map[int]int64

	// Pair tracking by qname (mate's flag/pos for orientation classification).
	// Cleared as pairs are observed (memory-bounded by max pending mates).
	mates map[string]mateInfo

	// CHK checksum block: per-record crc32 values summed with 32-bit
	// overflow, mirroring upstream's update_checksum (stats.c:755).
	ChkNames uint32
	ChkReads uint32
	ChkQuals uint32

	// COV coverage distribution. Rather than retaining every covered
	// reference position (which would be O(genome) and OOM on whole-genome
	// BAMs), depth is accumulated in a bounded sliding window for the
	// *current contig only* and finalized positions are binned into cov
	// incrementally — mirroring upstream's fixed-size cov_rbuf ring buffer
	// and round_buffer_flush (stats.c:326). covWindow holds depth for
	// positions >= covFlushed on contig covContig; covContig/covFlushed
	// track the streaming flush frontier. The cov bin array is sized once
	// up front from parseCoverageBins so binning can happen during the
	// flush. COV is emitted only for coordinate-sorted input, matching
	// upstream's is_sorted gating (stats.c:1848).
	covWindow  map[int32]int32 // per-position depth for the current contig
	covContig  string          // RName the window currently tracks ("" = none yet)
	covFlushed int32           // positions < this have been flushed for covContig
	cov        []int64         // COV bin array, sized from covMin/covMax/covStep
	covMin     int
	covMax     int
	covStep    int
	ncov       int
	covBinsSet bool // true once the cov bins have been initialized
}

// Buffer dimensions mirroring upstream stats.c initial allocation.
const (
	statsNQuals = 256 // quality values 0..255
	statsNGC    = 200 // GC-content histogram bins
)

// cycleQuals is the per-cycle quality histogram (one count per quality value).
type cycleQuals [statsNQuals]int64

// acgtNoCount mirrors upstream's acgtno_count_t: per-cycle base tallies.
type acgtNoCount struct {
	a, c, g, t, n, other int64
}

type mateInfo struct {
	pos    int32
	end    int32
	rname  string
	flag   uint16
	tlen   int32
	insert int32
}

// newStatsCounters returns a fresh counter set with maps pre-allocated.
func newStatsCounters() *StatsCounters {
	return &StatsCounters{
		RL:          make(map[int64]int64),
		FRL:         make(map[int64]int64),
		LRL:         make(map[int64]int64),
		ISInw:       make(map[int64]int64),
		ISOutw:      make(map[int64]int64),
		ISOther:     make(map[int64]int64),
		insertions:  make(map[int]int64),
		deletions:   make(map[int]int64),
		insCycle1st: make(map[int]int64),
		insCycle2nd: make(map[int]int64),
		delCycle1st: make(map[int]int64),
		delCycle2nd: make(map[int]int64),
		mates:       make(map[string]mateInfo),
		covWindow:   make(map[int32]int32),
	}
}

// initCovBins computes the COV bin geometry once from the -c option and
// allocates the cov bin array. Calling it more than once is a no-op so the
// streaming flush can rely on the bins being ready.
func (c *StatsCounters) initCovBins(opts StatsOptions) {
	if c.covBinsSet {
		return
	}
	c.covMin, c.covMax, c.covStep, c.ncov = parseCoverageBins(opts.Coverage)
	c.cov = make([]int64, c.ncov)
	c.covBinsSet = true
}

// Stats reads a SAM/BAM stream, accumulates statistics, and writes the
// upstream-style text report to out.
func Stats(in io.Reader, out io.Writer, opts StatsOptions) error {
	if opts.MaxInsertSize <= 0 {
		opts.MaxInsertSize = 8000
	}
	br, err := sam.NewReader(in)
	if err != nil {
		return err
	}
	_ = br.Header() // header is read for side-effects (consumed from stream)
	c := newStatsCounters()
	// Compute the COV bin geometry once up front so depth can be binned
	// incrementally during the streaming flush.
	c.initCovBins(opts)
	// Match upstream stats.c:2327 — is_sorted starts at 1 and is demoted to
	// 0 the first time an out-of-order record appears. The @HD SO:coordinate
	// line is NOT consulted; an empty BAM is reported as sorted.
	c.IsSorted = 1

	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		c.observe(rec, opts)
	}
	if c.sortBroken {
		c.IsSorted = 0
	}
	// Flush whatever coverage depth remains in the window for the last
	// contig before the report is written.
	c.flushCoverageWindow(math.MaxInt32)
	return c.Write(out, opts)
}

// observe folds one record into the accumulator.
func (c *StatsCounters) observe(rec *sam.Record, opts StatsOptions) {
	// Sort-order detector — track whether records arrive coordinate-sorted.
	if !c.sortBroken && rec.IsMapped() {
		if c.lastRef == "" {
			c.lastRef = rec.RName
			c.lastPos = rec.Pos
		} else if rec.RName == c.lastRef {
			if rec.Pos < c.lastPos {
				c.sortBroken = true
			} else {
				c.lastPos = rec.Pos
			}
		} else {
			// New reference; positions only need to be sorted within a ref.
			c.lastRef = rec.RName
			c.lastPos = rec.Pos
		}
	}
	// CHK checksum — upstream's update_checksum (stats.c:755) runs for every
	// record that passed the user flag filters, BEFORE the secondary/
	// supplementary classification, so it must be folded in on all three
	// record branches. Note: upstream applies the -F/-f/-l flag filters
	// before update_checksum, whereas this code applies those filters only
	// on the primary path below; computing CHK for every record is a known,
	// pre-existing flag-ordering divergence that is correct for the default
	// invocation (no flag filters) and is not addressed here.
	c.updateChecksum(rec)

	// Branch 1: classify non-primary reads first (they don't count toward
	// "raw total sequences" upstream, they only bump their own counters).
	if rec.IsSecondary() {
		c.NonPrimary++
		return
	}
	if rec.IsSupplementary() {
		c.Supplementary++
		// Upstream's bases_mapped_cigar counter folds in supplementary
		// alignments — stats.c walks the cigar of EVERY mapped record
		// (primary + supplementary) into nbases_mapped_cigar.
		for _, op := range rec.Cigar {
			ch := op.Char()
			ln := int64(op.Length())
			switch ch {
			case 'M', '=', 'X', 'I':
				c.BasesMappedCigar += ln
			}
		}
		// Upstream count_indels runs for every mapped record, including
		// supplementary alignments (only secondary reads are excluded).
		if rec.IsMapped() {
			c.countIndels(rec)
			// COV folds in supplementary alignments too; upstream walks
			// the CIGAR of every record reaching collect_stats (only
			// secondary reads return early).
			c.accumulateCoverage(rec)
		}
		return
	}
	// QC-fail records ARE counted toward RawTotal/Sequences/totals per
	// upstream stats.c (collect_orig_read_stats runs for every primary
	// record, QCFAIL or not). Only the MAPQ histogram excludes QCFAIL+DUP
	// (see gating below).
	if rec.IsQCFail() {
		c.ReadsQCFailed++
	}

	// Apply user-specified flag filters BEFORE counting RawTotal so the
	// "filtered sequences" delta is meaningful (matches upstream).
	c.RawTotal++
	if opts.RequiredFlag != 0 && rec.Flag&opts.RequiredFlag != opts.RequiredFlag {
		c.FilteredOut++
		return
	}
	if opts.FilteringFlag != 0 && rec.Flag&opts.FilteringFlag != 0 {
		c.FilteredOut++
		return
	}
	if opts.RemoveDups && rec.IsDuplicate() {
		c.FilteredOut++
		return
	}
	if opts.MinMAPQ > 0 && rec.MapQ < opts.MinMAPQ {
		c.FilteredOut++
		return
	}
	c.Sequences++

	// Length (raw, ignores clipping)
	rlen := int64(len(rec.Seq))
	c.TotalLength += rlen
	if rlen > c.MaxLength {
		c.MaxLength = rlen
	}
	c.RL[rlen]++
	if rec.IsRead2() {
		c.LastFrag++
		c.TotalLastFragLength += rlen
		if rlen > c.MaxLastFragLength {
			c.MaxLastFragLength = rlen
		}
		c.LRL[rlen]++
	} else {
		// upstream treats not-Read2 (incl. unpaired) as "1st fragment".
		c.FirstFrag++
		c.TotalFirstFragLength += rlen
		if rlen > c.MaxFirstFragLength {
			c.MaxFirstFragLength = rlen
		}
		c.FRL[rlen]++
	}

	// Per-cycle quality, base-content and GC accumulation. Mirrors upstream
	// collect_orig_read_stats — runs for every QC-passed primary record.
	c.collectCycleStats(rec)

	if rec.IsMapped() {
		c.ReadsMapped++
		if rec.IsPaired() && !rec.IsMateUnmapped() {
			c.ReadsMappedAndPaired++
		}
		c.BasesMapped += rlen
		// MAPQ histogram excludes UNMAP|SEC|SUPP|QCFAIL|DUP
		// (upstream stats.c:1239). UNMAP/SEC/SUPP are already excluded by
		// outer branches; gate QCFAIL+DUP here.
		if !rec.IsQCFail() && !rec.IsDuplicate() {
			c.MAPQ[rec.MapQ]++
		}
		// ReadsMQ0 (the SN counter) INCLUDES QCFAIL+DUP — upstream's
		// nreads_mq0 lives inside collect_orig_read_stats which is called
		// for every primary record regardless of QCFAIL/DUP.
		if rec.MapQ == 0 {
			c.ReadsMQ0++
		}
	} else {
		c.ReadsUnmapped++
	}
	if rec.IsPaired() {
		c.ReadsPaired++
		if rec.IsProperPair() {
			c.ReadsProperlyPaired++
		}
	}
	if rec.IsDuplicate() {
		c.ReadsDuplicated++
		c.BasesDuplicated += rlen
	}

	// Cigar-aware "bases mapped (cigar)" — matches upstream `nbases_mapped_cigar`
	// which counts M, =, X *and* I bases (see stats.c line 1581 comment).
	// Note: `bases trimmed` in upstream tracks BWA-style quality trimming
	// (see bwa_trim_read). We don't implement trim-quality scanning, so
	// BasesTrimmed stays 0 — matching upstream's default behaviour when
	// `--trim-quality` (`-q`) is not passed. Documented in PARITY_ROADMAP.md.
	for _, op := range rec.Cigar {
		ch := op.Char()
		ln := int64(op.Length())
		switch ch {
		case 'M', '=', 'X', 'I':
			c.BasesMappedCigar += ln
		}
	}

	// Indel distribution and per-cycle indel counts (mapped reads only).
	c.countIndels(rec)

	// COV coverage-distribution depth (mapped reads only).
	if rec.IsMapped() {
		c.accumulateCoverage(rec)
	}

	// Mismatches via NM aux tag.
	if a, ok := rec.GetAux("NM"); ok {
		if v, ok := a.Int(); ok {
			c.Mismatches += v
		}
	}

	// Quality sum (only for QC-passed sequences).
	// Upstream uses 255 ("missing") for "*" qualities, which gives the
	// large "average quality: 255.0" observed in the fixture.
	if len(rec.Qual) == 0 {
		// "*" qual — encode 255 per base, matching upstream's missing-qual handling.
		c.TotalQual += 255 * rlen
	} else {
		for _, q := range rec.Qual {
			c.TotalQual += int64(q)
		}
	}
	c.TotalQualBases += rlen

	// Pair orientation + insert size classification.
	if rec.IsPaired() && rec.IsMapped() && !rec.IsMateUnmapped() {
		c.classifyPair(rec, opts)
	}
}

// classifyPair counts each pair once on the SECOND-observed mate. The
// insert size source is upstream's `bam_line->core.isize` (TLEN, abs value
// capped at MaxInsertSize) per stats.c:1265-1268, NOT a position-derived
// computation. Orientation is classified from the leftmost mate's strand
// (inward = leftward forward + rightward reverse).
func (c *StatsCounters) classifyPair(rec *sam.Record, opts StatsOptions) {
	if prev, ok := c.mates[rec.QName]; ok {
		delete(c.mates, rec.QName)
		// Different chromosome?
		if rec.RName != "" && prev.rname != "" && rec.RName != prev.rname {
			c.DiffChromosomePairs++
			return
		}
		// Same chromosome — classify orientation from the leftmost.
		left, right := prev, mateInfo{pos: rec.Pos, end: rec.EndPosition(), rname: rec.RName, flag: rec.Flag, tlen: rec.TLen}
		if rec.Pos < prev.pos {
			left, right = right, prev
		}
		leftRev := left.flag&sam.FlagReverse != 0
		rightRev := right.flag&sam.FlagReverse != 0
		var bucket map[int64]int64
		switch {
		case !leftRev && rightRev:
			c.InwardPairs++
			bucket = c.ISInw
		case leftRev && !rightRev:
			c.OutwardPairs++
			bucket = c.ISOutw
		default:
			c.OtherOrientPairs++
			bucket = c.ISOther
		}
		// Insert size: |TLEN| of the second-observed record. Upstream
		// caps at -i/--max-insert (default 8000) by clamping, not skipping.
		ins := int64(rec.TLen)
		if ins < 0 {
			ins = -ins
		}
		if ins > int64(opts.MaxInsertSize) {
			ins = int64(opts.MaxInsertSize)
		}
		if ins > 0 {
			bucket[ins]++
			c.InsertSumAbs += ins
			c.InsertSumSqAbs += ins * ins
			c.InsertPairsCounted++
		}
		return
	}
	c.mates[rec.QName] = mateInfo{
		pos:   rec.Pos,
		end:   rec.EndPosition(),
		rname: rec.RName,
		flag:  rec.Flag,
		tlen:  rec.TLen,
	}
}

// Read-order classifications mirroring upstream READ_ORDER_FIRST/LAST.
const (
	orderOther = 0
	orderFirst = 1
	orderLast  = 2
)

// readOrder classifies a record as first/last fragment or "other", matching
// upstream stats.c: unpaired reads count as first fragments; a paired read is
// first if READ1 is set, last if READ2 is set, otherwise "other".
func readOrder(rec *sam.Record) int {
	if !rec.IsPaired() {
		return orderFirst
	}
	switch {
	case rec.IsRead1() && !rec.IsRead2():
		return orderFirst
	case rec.IsRead2() && !rec.IsRead1():
		return orderLast
	default:
		return orderOther
	}
}

// unclippedLength returns the read length including hard-clipped bases, the
// length upstream uses to size per-cycle buffers (unclipped_length in stats.c).
func unclippedLength(rec *sam.Record) int {
	n := len(rec.Seq)
	for _, op := range rec.Cigar {
		if op.Char() == 'H' {
			n += int(op.Length())
		}
	}
	return n
}

// growCycles ensures sl has at least n entries.
func growCycles(sl []acgtNoCount, n int) []acgtNoCount {
	for len(sl) < n {
		sl = append(sl, acgtNoCount{})
	}
	return sl
}

// growQuals ensures sl has at least n entries.
func growQuals(sl []cycleQuals, n int) []cycleQuals {
	for len(sl) < n {
		sl = append(sl, cycleQuals{})
	}
	return sl
}

// collectCycleStats folds one QC-passed primary record into the per-cycle
// quality, base-content and GC-content accumulators (upstream
// collect_orig_read_stats).
func (c *StatsCounters) collectCycleStats(rec *sam.Record) {
	seqLen := len(rec.Seq)
	order := readOrder(rec)
	reverse := rec.Flag&sam.FlagReverse != 0

	ulen := unclippedLength(rec)
	if ulen > c.maxLen {
		c.maxLen = ulen
	}
	if order == orderFirst && ulen > c.maxLen1st {
		c.maxLen1st = ulen
	}
	if order == orderLast && ulen > c.maxLen2nd {
		c.maxLen2nd = ulen
	}
	if seqLen == 0 {
		return
	}

	c.acgtRevcomp = growCycles(c.acgtRevcomp, seqLen)
	var cycles []acgtNoCount
	switch order {
	case orderFirst:
		c.acgtCycles1st = growCycles(c.acgtCycles1st, seqLen)
		cycles = c.acgtCycles1st
	case orderLast:
		c.acgtCycles2nd = growCycles(c.acgtCycles2nd, seqLen)
		cycles = c.acgtCycles2nd
	}

	gcCount := 0
	for i := 0; i < seqLen; i++ {
		readCycle := i
		if reverse {
			readCycle = seqLen - i - 1
		}
		switch rec.Seq[i] {
		case 'A', 'a':
			if cycles != nil {
				cycles[readCycle].a++
			}
			if reverse {
				c.acgtRevcomp[readCycle].t++
			} else {
				c.acgtRevcomp[readCycle].a++
			}
		case 'C', 'c':
			if cycles != nil {
				cycles[readCycle].c++
			}
			if reverse {
				c.acgtRevcomp[readCycle].g++
			} else {
				c.acgtRevcomp[readCycle].c++
			}
			gcCount++
		case 'G', 'g':
			if cycles != nil {
				cycles[readCycle].g++
			}
			if reverse {
				c.acgtRevcomp[readCycle].c++
			} else {
				c.acgtRevcomp[readCycle].g++
			}
			gcCount++
		case 'T', 't':
			if cycles != nil {
				cycles[readCycle].t++
			}
			if reverse {
				c.acgtRevcomp[readCycle].a++
			} else {
				c.acgtRevcomp[readCycle].t++
			}
		case 'N', 'n':
			if cycles != nil {
				cycles[readCycle].n++
			}
		default:
			if cycles != nil {
				cycles[readCycle].other++
			}
		}
	}

	// GC-content histogram: spread the read's GC fraction over a [min,max) bin
	// range, matching upstream's gc_idx_min/gc_idx_max integer arithmetic.
	gcIdxMin := gcCount * (statsNGC - 1) / seqLen
	gcIdxMax := (gcCount + 1) * (statsNGC - 1) / seqLen
	if gcIdxMax >= statsNGC {
		gcIdxMax = statsNGC - 1
	}
	switch order {
	case orderFirst:
		for i := gcIdxMin; i < gcIdxMax; i++ {
			c.gcFirst[i]++
		}
	case orderLast:
		for i := gcIdxMin; i < gcIdxMax; i++ {
			c.gcLast[i]++
		}
	}

	// Quality histogram per cycle. A "*" quality string (rec.Qual empty) is
	// treated as quality 255 per base, matching upstream's missing-qual value.
	var quals []cycleQuals
	switch order {
	case orderFirst:
		c.qualsFirst = growQuals(c.qualsFirst, seqLen)
		quals = c.qualsFirst
	case orderLast:
		c.qualsLast = growQuals(c.qualsLast, seqLen)
		quals = c.qualsLast
	}
	if quals != nil {
		for i := 0; i < seqLen; i++ {
			q := 255
			if len(rec.Qual) == seqLen {
				idx := i
				if reverse {
					idx = seqLen - i - 1
				}
				q = int(rec.Qual[idx])
			}
			if q >= statsNQuals {
				q = statsNQuals - 1
			}
			if q > c.maxQual {
				c.maxQual = q
			}
			quals[i][q]++
		}
	}
}

// countIndels folds insertion/deletion cigar operations of one mapped record
// into the indel-length and per-cycle indel accumulators (upstream
// count_indels).
func (c *StatsCounters) countIndels(rec *sam.Record) {
	isFwd := rec.Flag&sam.FlagReverse == 0
	order := readOrder(rec)
	readLen := len(rec.Seq)
	icycle := 0
	for _, op := range rec.Cigar {
		ch := op.Char()
		ncig := int(op.Length())
		if ncig == 0 {
			continue
		}
		switch ch {
		case 'I':
			idx := icycle
			if !isFwd {
				idx = readLen - icycle - ncig
			}
			if idx >= 0 {
				switch order {
				case orderFirst:
					c.insCycle1st[idx]++
				case orderLast:
					c.insCycle2nd[idx]++
				}
			}
			icycle += ncig
			c.insertions[ncig-1]++
		case 'D':
			idx := icycle - 1
			if !isFwd {
				idx = readLen - icycle - 1
			}
			if idx >= 0 {
				switch order {
				case orderFirst:
					c.delCycle1st[idx]++
				case orderLast:
					c.delCycle2nd[idx]++
				}
			}
			c.deletions[ncig-1]++
		case 'N', 'H', 'P':
			// Reference skips, hard clips and padding do not advance the cycle.
		default:
			icycle += ncig
		}
	}
}

// updateChecksum folds one record into the CHK checksum block. It mirrors
// upstream's update_checksum (stats.c:755): the qname, the BAM 4-bit-packed
// sequence and the quality bytes are each crc32'd with the IEEE polynomial
// (identical to zlib's crc32) and summed into the running totals with 32-bit
// overflow. A record with no sequence (SEQ "*") contributes only its name.
func (c *StatsCounters) updateChecksum(rec *sam.Record) {
	c.ChkNames += crc32.ChecksumIEEE([]byte(rec.QName))

	seqLen := len(rec.Seq)
	if seqLen == 0 {
		return
	}
	// PackedSeq yields the BAM nibble encoding byte-identical to htslib's
	// bam_get_seq buffer, so crc32 over it matches upstream exactly.
	c.ChkReads += crc32.ChecksumIEEE(rec.PackedSeq())

	if len(rec.Qual) == seqLen {
		c.chkQualsAdd(rec.Qual)
		return
	}
	// QUAL "*" — upstream's bam_get_qual buffer is all 0xff for seq_len bytes.
	missing := make([]byte, seqLen)
	for i := range missing {
		missing[i] = 0xff
	}
	c.chkQualsAdd(missing)
}

// chkQualsAdd folds one quality buffer into the CHK qualities checksum.
func (c *StatsCounters) chkQualsAdd(qual []byte) {
	c.ChkQuals += crc32.ChecksumIEEE(qual)
}

// accumulateCoverage folds one mapped record's per-position depth into the
// bounded coverage window. Only M, = and X CIGAR operations contribute depth;
// soft-clips, insertions, deletions and reference skips do not (stats.c
// comment at line 29 and the cov_rbuf insertion loop at stats.c:1457).
//
// Because COV is emitted only for coordinate-sorted input, records arrive in
// non-decreasing position order: once a record starts at rec.Pos every
// reference position < rec.Pos is final and can be binned and dropped. This
// mirrors upstream's fixed-size cov_rbuf ring buffer + round_buffer_flush
// (stats.c:326), keeping COV memory O(longest read span) instead of
// O(genome). On a contig change the previous contig's window is flushed in
// full first. Defensive against out-of-order input: a backwards jump never
// rewinds covFlushed, so the window and flush logic cannot be corrupted (the
// COV section is suppressed anyway when the input is not sorted).
func (c *StatsCounters) accumulateCoverage(rec *sam.Record) {
	if rec.RName != c.covContig {
		// New contig: every windowed position of the previous contig is
		// final. Flush it all, then switch contigs.
		c.flushCoverageWindow(math.MaxInt32)
		c.covContig = rec.RName
		c.covFlushed = 0
	}
	// Every position strictly below rec.Pos can receive no further depth.
	c.flushCoverageWindow(rec.Pos)

	pos := rec.Pos
	for _, op := range rec.Cigar {
		ch := op.Char()
		ln := int32(op.Length())
		switch ch {
		case 'M', '=', 'X':
			for i := int32(0); i < ln; i++ {
				p := pos + i
				// Defensive: never re-add depth to an already-flushed
				// position (only possible on out-of-order input, where
				// COV is suppressed regardless).
				if p >= c.covFlushed {
					c.covWindow[p]++
				}
			}
			pos += ln
		case 'D', 'N':
			// Consume reference without contributing depth.
			pos += ln
		default:
			// I, S, H, P do not consume reference.
		}
	}
}

// flushCoverageWindow bins and removes every windowed position strictly below
// upTo, advancing the flush frontier. Passing math.MaxInt32 flushes the whole
// window (used on a contig change and at end of input). It is the streaming
// equivalent of upstream's round_buffer_flush (stats.c:326).
func (c *StatsCounters) flushCoverageWindow(upTo int32) {
	if len(c.covWindow) == 0 {
		if upTo > c.covFlushed {
			c.covFlushed = upTo
		}
		return
	}
	for p, d := range c.covWindow {
		if p < upTo {
			c.cov[coverageIdx(c.covMin, c.covMax, c.ncov, c.covStep, int(d))]++
			delete(c.covWindow, p)
		}
	}
	if upTo > c.covFlushed {
		c.covFlushed = upTo
	}
}

// coverageIdx maps a read depth to its COV bin index, mirroring upstream's
// coverage_idx (stats.c:310): bin 0 collects depths below covMin, bin n-1
// collects depths above covMax, and the middle bins are stepped by covStep.
func coverageIdx(min, max, n, step, depth int) int {
	if depth < min {
		return 0
	}
	if depth > max {
		return n - 1
	}
	return 1 + (depth-min)/step
}

// parseCoverageBins parses the -c "MIN,MAX,STEP" string and returns the
// effective covMin, covMax, covStep and bin count, applying the same step and
// max adjustments upstream performs at stats.c:2359-2366. An empty string
// yields the upstream defaults 1,1000,1.
func parseCoverageBins(spec string) (covMin, covMax, covStep, ncov int) {
	covMin, covMax, covStep = 1, 1000, 1
	if spec != "" {
		var mn, mx, st int
		if n, _ := fmt.Sscanf(spec, "%d,%d,%d", &mn, &mx, &st); n == 3 {
			covMin, covMax, covStep = mn, mx, st
		}
	}
	if covStep > covMax-covMin+1 {
		covStep = covMax - covMin
		if covStep <= 0 {
			covStep = 1
		}
	}
	ncov = 3 + (covMax-covMin)/covStep
	covMax = covMin + ((covMax-covMin)/covStep+1)*covStep - 1
	return covMin, covMax, covStep, ncov
}

// Write emits the upstream-compatible text report to w.
func (c *StatsCounters) Write(w io.Writer, opts StatsOptions) error {
	bw := bufio.NewWriter(w)
	// CHK block — emitted first, ahead of SN, matching upstream stats.c:1557.
	c.writeCHK(bw)
	// SN block — full upstream parity.
	c.writeSN(bw)
	// Histogram and per-cycle section bodies are emitted only in the
	// non-sparse path. FFQ/LFQ/GCF/GCL/GCC/GCT/IC/ID carry real data and are
	// byte-faithful to upstream; section order matches stats.c.
	if !opts.Sparse {
		c.writeFFQLFQ(bw)
		c.writeGCFGCL(bw)
		c.writeGCC(bw)
		c.writeGCT(bw)
		c.writeRL(bw)
		c.writeMAPQ(bw)
		c.writeIS(bw, opts)
		c.writeID(bw)
		c.writeIC(bw)
		c.writeCOV(bw, opts)
	}
	return bw.Flush()
}

// writeCHK emits the leading CRC32 checksum block (read names, sequences,
// qualities), matching upstream stats.c:1557-1559.
func (c *StatsCounters) writeCHK(bw *bufio.Writer) {
	fmt.Fprintln(bw, "# CHK, Checksum\t[2]Read Names\t[3]Sequences\t[4]Qualities")
	fmt.Fprintln(bw, "# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)")
	fmt.Fprintf(bw, "CHK\t%08x\t%08x\t%08x\n", c.ChkNames, c.ChkReads, c.ChkQuals)
}

// writeCOV emits the coverage-distribution histogram. It is only emitted for
// coordinate-sorted input, matching upstream's is_sorted gating
// (stats.c:1848). The cov bin array is populated incrementally during the
// streaming flush (see accumulateCoverage/flushCoverageWindow); writeCOV just
// formats it, printing only non-empty bins (stats.c:1850-1856).
func (c *StatsCounters) writeCOV(bw *bufio.Writer, opts StatsOptions) {
	if c.IsSorted != 1 {
		return
	}
	c.initCovBins(opts)
	covMin, covMax, covStep, ncov := c.covMin, c.covMax, c.covStep, c.ncov
	_ = covMax
	cov := c.cov
	fmt.Fprintln(bw, "# Coverage distribution. Use `grep ^COV | cut -f 2-` to extract this part.")
	if cov[0] > 0 {
		fmt.Fprintf(bw, "COV\t[<%d]\t%d\t%d\n", covMin, covMin-1, cov[0])
	}
	for icov := 1; icov < ncov-1; icov++ {
		if cov[icov] == 0 {
			continue
		}
		lo := covMin + (icov-1)*covStep
		hi := covMin + icov*covStep - 1
		fmt.Fprintf(bw, "COV\t[%d-%d]\t%d\t%d\n", lo, hi, hi, cov[icov])
	}
	if cov[ncov-1] > 0 {
		edge := covMin + (ncov-2)*covStep - 1
		fmt.Fprintf(bw, "COV\t[%d<]\t%d\t%d\n", edge, edge, cov[ncov-1])
	}
}

// writeSN emits the Summary Numbers section.
func (c *StatsCounters) writeSN(bw *bufio.Writer) {
	fmt.Fprintln(bw, "# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.")
	emit := func(key string, val interface{}, comment string) {
		if comment == "" {
			fmt.Fprintf(bw, "SN\t%s:\t%v\n", key, val)
			return
		}
		fmt.Fprintf(bw, "SN\t%s:\t%v\t%s\n", key, val, comment)
	}
	emit("raw total sequences", c.RawTotal, "# excluding supplementary and secondary reads")
	emit("filtered sequences", c.FilteredOut, "")
	emit("sequences", c.Sequences, "")
	if c.IsSorted == 1 {
		emit("is sorted", 1, "# sorted by coordinate")
	} else {
		emit("is sorted", 0, "# not sorted by coordinate")
	}
	emit("1st fragments", c.FirstFrag, "")
	emit("last fragments", c.LastFrag, "")
	emit("reads mapped", c.ReadsMapped, "")
	emit("reads mapped and paired", c.ReadsMappedAndPaired, "# paired-end technology bit set + both mates mapped")
	emit("reads unmapped", c.ReadsUnmapped, "")
	emit("reads properly paired", c.ReadsProperlyPaired, "# proper-pair bit set")
	emit("reads paired", c.ReadsPaired, "# paired-end technology bit set")
	emit("reads duplicated", c.ReadsDuplicated, "# PCR or optical duplicate bit set")
	emit("reads MQ0", c.ReadsMQ0, "# mapped and MQ=0")
	emit("reads QC failed", c.ReadsQCFailed, "")
	emit("non-primary alignments", c.NonPrimary, "")
	emit("supplementary alignments", c.Supplementary, "")
	emit("total length", c.TotalLength, "# ignores clipping")
	emit("total first fragment length", c.TotalFirstFragLength, "# ignores clipping")
	emit("total last fragment length", c.TotalLastFragLength, "# ignores clipping")
	emit("bases mapped", c.BasesMapped, "# ignores clipping")
	emit("bases mapped (cigar)", c.BasesMappedCigar, "# more accurate")
	emit("bases trimmed", c.BasesTrimmed, "")
	emit("bases duplicated", c.BasesDuplicated, "")
	emit("mismatches", c.Mismatches, "# from NM fields")
	// error rate uses scientific notation, 6-digit precision matching
	// upstream's `%.6e` printout.
	errRate := 0.0
	if c.BasesMappedCigar > 0 {
		errRate = float64(c.Mismatches) / float64(c.BasesMappedCigar)
	}
	emit("error rate", fmt.Sprintf("%.6e", errRate), "# mismatches / bases mapped (cigar)")
	avgLen := int64(0)
	if c.Sequences > 0 {
		avgLen = c.TotalLength / c.Sequences
	}
	emit("average length", avgLen, "")
	avgFirst := int64(0)
	if c.FirstFrag > 0 {
		avgFirst = c.TotalFirstFragLength / c.FirstFrag
	}
	emit("average first fragment length", avgFirst, "")
	avgLast := int64(0)
	if c.LastFrag > 0 {
		avgLast = c.TotalLastFragLength / c.LastFrag
	}
	emit("average last fragment length", avgLast, "")
	emit("maximum length", c.MaxLength, "")
	emit("maximum first fragment length", c.MaxFirstFragLength, "")
	emit("maximum last fragment length", c.MaxLastFragLength, "")
	// Average quality: 1-decimal precision per upstream.
	avgQual := 0.0
	if c.TotalQualBases > 0 {
		avgQual = float64(c.TotalQual) / float64(c.TotalQualBases)
	}
	emit("average quality", fmt.Sprintf("%.1f", avgQual), "")
	// Insert-size mean and stddev: |TLEN| over inward+outward+other pairs.
	insAvg := 0.0
	insSD := 0.0
	if c.InsertPairsCounted > 0 {
		insAvg = float64(c.InsertSumAbs) / float64(c.InsertPairsCounted)
		mean2 := float64(c.InsertSumSqAbs) / float64(c.InsertPairsCounted)
		variance := mean2 - insAvg*insAvg
		if variance > 0 {
			insSD = math.Sqrt(variance)
		}
	}
	emit("insert size average", fmt.Sprintf("%.1f", insAvg), "")
	emit("insert size standard deviation", fmt.Sprintf("%.1f", insSD), "")
	emit("inward oriented pairs", c.InwardPairs, "")
	emit("outward oriented pairs", c.OutwardPairs, "")
	emit("pairs with other orientation", c.OtherOrientPairs, "")
	emit("pairs on different chromosomes", c.DiffChromosomePairs, "")
	// Proper-pair % uses Sequences (= 1st + 2nd + other) as the denominator,
	// matching upstream stats.c:1606. NOT ReadsPaired — they differ when
	// records have the paired bit but no first/last fragment classification.
	pp := 0.0
	if c.Sequences > 0 {
		pp = 100.0 * float64(c.ReadsProperlyPaired) / float64(c.Sequences)
	}
	emit("percentage of properly paired reads (%)", fmt.Sprintf("%.1f", pp), "")
}

// effectiveMaxQual returns the highest quality column index to print. Upstream
// bumps max_qual by one (so a trailing all-zero column is shown) provided that
// stays within the quality buffer.
func (c *StatsCounters) effectiveMaxQual() int {
	mq := c.maxQual
	if mq+1 < statsNQuals {
		mq++
	}
	return mq
}

// writeFFQLFQ emits the per-cycle quality histograms for first and last
// fragments.
func (c *StatsCounters) writeFFQLFQ(bw *bufio.Writer) {
	maxQual := c.effectiveMaxQual()
	fmt.Fprintln(bw, "# First Fragment Qualities. Use `grep ^FFQ | cut -f 2-` to extract this part.")
	fmt.Fprintln(bw, "# Columns correspond to qualities and rows to cycles. First column is the cycle number.")
	c.writeQualCycles(bw, "FFQ", c.qualsFirst, c.maxLen1st, maxQual)
	fmt.Fprintln(bw, "# Last Fragment Qualities. Use `grep ^LFQ | cut -f 2-` to extract this part.")
	fmt.Fprintln(bw, "# Columns correspond to qualities and rows to cycles. First column is the cycle number.")
	c.writeQualCycles(bw, "LFQ", c.qualsLast, c.maxLen2nd, maxQual)
}

// writeQualCycles writes one quality-histogram block (tag is FFQ or LFQ): one
// row per cycle 1..maxLen, each row listing counts for qualities 0..maxQual.
func (c *StatsCounters) writeQualCycles(bw *bufio.Writer, tag string, quals []cycleQuals, maxLen, maxQual int) {
	for ibase := 0; ibase < maxLen; ibase++ {
		fmt.Fprintf(bw, "%s\t%d", tag, ibase+1)
		for iqual := 0; iqual <= maxQual; iqual++ {
			var v int64
			if ibase < len(quals) {
				v = quals[ibase][iqual]
			}
			fmt.Fprintf(bw, "\t%d", v)
		}
		fmt.Fprintln(bw)
	}
}

// writeGCFGCL emits the GC-content sections for first and last fragments.
// Consecutive bins with an equal count are collapsed exactly as upstream does.
func (c *StatsCounters) writeGCFGCL(bw *bufio.Writer) {
	fmt.Fprintln(bw, "# GC Content of first fragments. Use `grep ^GCF | cut -f 2-` to extract this part.")
	writeGCHist(bw, "GCF", c.gcFirst[:])
	fmt.Fprintln(bw, "# GC Content of last fragments. Use `grep ^GCL | cut -f 2-` to extract this part.")
	writeGCHist(bw, "GCL", c.gcLast[:])
}

// writeGCHist writes one GC-content histogram. The GC% of each emitted row is
// the midpoint of the run of equal-count bins, per upstream's formula.
func writeGCHist(bw *bufio.Writer, tag string, gc []int64) {
	prev := 0
	for ibase := 0; ibase < len(gc); ibase++ {
		if gc[ibase] == gc[prev] {
			continue
		}
		pct := float64(ibase+prev) * 0.5 * 100.0 / float64(statsNGC-1)
		fmt.Fprintf(bw, "%s\t%.2f\t%d\n", tag, pct, gc[prev])
		prev = ibase
	}
}

// writeGCC emits the as-sequenced ACGT-content-per-cycle section, summing the
// first- and last-fragment cycle buffers.
func (c *StatsCounters) writeGCC(bw *bufio.Writer) {
	fmt.Fprintln(bw, "# ACGT content per cycle. Use `grep ^GCC | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]")
	for ibase := 0; ibase < c.maxLen; ibase++ {
		var first, last acgtNoCount
		if ibase < len(c.acgtCycles1st) {
			first = c.acgtCycles1st[ibase]
		}
		if ibase < len(c.acgtCycles2nd) {
			last = c.acgtCycles2nd[ibase]
		}
		sum := first.a + first.c + first.g + first.t + last.a + last.c + last.g + last.t
		if sum == 0 {
			continue
		}
		fs := float64(sum)
		fmt.Fprintf(bw, "GCC\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", ibase+1,
			100.0*float64(first.a+last.a)/fs,
			100.0*float64(first.c+last.c)/fs,
			100.0*float64(first.g+last.g)/fs,
			100.0*float64(first.t+last.t)/fs,
			100.0*float64(first.n+last.n)/fs,
			100.0*float64(first.other+last.other)/fs)
	}
}

// writeGCT emits the read-oriented ACGT-content-per-cycle section, where
// reverse-strand reads have contributed reverse-complemented bases.
func (c *StatsCounters) writeGCT(bw *bufio.Writer) {
	fmt.Fprintln(bw, "# ACGT content per cycle, read oriented. Use `grep ^GCT | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]")
	for ibase := 0; ibase < c.maxLen; ibase++ {
		var rc acgtNoCount
		if ibase < len(c.acgtRevcomp) {
			rc = c.acgtRevcomp[ibase]
		}
		sum := rc.a + rc.c + rc.g + rc.t
		if sum == 0 {
			continue
		}
		fs := float64(sum)
		fmt.Fprintf(bw, "GCT\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n", ibase+1,
			100.0*float64(rc.a)/fs,
			100.0*float64(rc.c)/fs,
			100.0*float64(rc.g)/fs,
			100.0*float64(rc.t)/fs)
	}
}

// writeID emits the indel-distribution section: insertion and deletion counts
// keyed by indel length.
func (c *StatsCounters) writeID(bw *bufio.Writer) {
	fmt.Fprintln(bw, "# Indel distribution. Use `grep ^ID | cut -f 2-` to extract this part. The columns are: length, number of insertions, number of deletions")
	maxLen := 0
	for k := range c.insertions {
		if k+1 > maxLen {
			maxLen = k + 1
		}
	}
	for k := range c.deletions {
		if k+1 > maxLen {
			maxLen = k + 1
		}
	}
	for ilen := 0; ilen < maxLen; ilen++ {
		ins, del := c.insertions[ilen], c.deletions[ilen]
		if ins > 0 || del > 0 {
			fmt.Fprintf(bw, "ID\t%d\t%d\t%d\n", ilen+1, ins, del)
		}
	}
}

// writeIC emits the indels-per-cycle section: insertion and deletion counts
// (forward and reverse) keyed by cycle.
func (c *StatsCounters) writeIC(bw *bufio.Writer) {
	fmt.Fprintln(bw, "# Indels per cycle. Use `grep ^IC | cut -f 2-` to extract this part. The columns are: cycle, number of insertions (fwd), .. (rev) , number of deletions (fwd), .. (rev)")
	maxLen := 0
	for _, m := range []map[int]int64{c.insCycle1st, c.insCycle2nd, c.delCycle1st, c.delCycle2nd} {
		for k := range m {
			if k+1 > maxLen {
				maxLen = k + 1
			}
		}
	}
	for ilen := 0; ilen < maxLen; ilen++ {
		insF, insR := c.insCycle1st[ilen], c.insCycle2nd[ilen]
		delF, delR := c.delCycle1st[ilen], c.delCycle2nd[ilen]
		if insF > 0 || insR > 0 || delF > 0 || delR > 0 {
			fmt.Fprintf(bw, "IC\t%d\t%d\t%d\t%d\t%d\n", ilen+1, insF, insR, delF, delR)
		}
	}
}

// writeRL emits the Read Length sections.
func (c *StatsCounters) writeRL(bw *bufio.Writer) {
	fmt.Fprintln(bw, "# Read lengths. Use `grep ^RL | cut -f 2-` to extract this part. The columns are: read length, count")
	for _, k := range sortedKeys(c.RL) {
		fmt.Fprintf(bw, "RL\t%d\t%d\n", k, c.RL[k])
	}
	fmt.Fprintln(bw, "# Read lengths - first fragments. Use `grep ^FRL | cut -f 2-` to extract this part. The columns are: read length, count")
	for _, k := range sortedKeys(c.FRL) {
		fmt.Fprintf(bw, "FRL\t%d\t%d\n", k, c.FRL[k])
	}
	fmt.Fprintln(bw, "# Read lengths - last fragments. Use `grep ^LRL | cut -f 2-` to extract this part. The columns are: read length, count")
	for _, k := range sortedKeys(c.LRL) {
		fmt.Fprintf(bw, "LRL\t%d\t%d\n", k, c.LRL[k])
	}
}

// writeMAPQ emits the MAPQ histogram section.
func (c *StatsCounters) writeMAPQ(bw *bufio.Writer) {
	fmt.Fprintln(bw, "# Mapping qualities for reads !(UNMAP|SECOND|SUPPL|QCFAIL|DUP). Use `grep ^MAPQ | cut -f 2-` to extract this part. The columns are: mapq, count")
	for i := 0; i < len(c.MAPQ); i++ {
		if c.MAPQ[i] == 0 {
			continue
		}
		fmt.Fprintf(bw, "MAPQ\t%d\t%d\n", i, c.MAPQ[i])
	}
}

// writeIS emits the Insert Size section.
func (c *StatsCounters) writeIS(bw *bufio.Writer, opts StatsOptions) {
	fmt.Fprintln(bw, "# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part. The columns are: insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs")
	// Union of all insert-size keys, sorted.
	seen := make(map[int64]struct{})
	for k := range c.ISInw {
		seen[k] = struct{}{}
	}
	for k := range c.ISOutw {
		seen[k] = struct{}{}
	}
	for k := range c.ISOther {
		seen[k] = struct{}{}
	}
	keys := make([]int64, 0, len(seen))
	for k := range seen {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })
	for _, k := range keys {
		inw := c.ISInw[k]
		outw := c.ISOutw[k]
		oth := c.ISOther[k]
		total := inw + outw + oth
		fmt.Fprintf(bw, "IS\t%d\t%d\t%d\t%d\t%d\n", k, total, inw, outw, oth)
	}
}

// sortedKeys returns the keys of m in ascending order.
func sortedKeys(m map[int64]int64) []int64 {
	ks := make([]int64, 0, len(m))
	for k := range m {
		ks = append(ks, k)
	}
	sort.Slice(ks, func(i, j int) bool { return ks[i] < ks[j] })
	return ks
}
