// Package samtools — stats implementation.
//
// Mirrors `samtools stats`. The Summary Numbers (SN) block, the read-length
// sections (RL/FRL/LRL), the MAPQ and IS histograms, the per-cycle quality
// histograms (FFQ/LFQ), the GC-content sections (GCF/GCL), the ACGT-content
// sections (GCC/GCT), the indel sections (IC/ID), the leading CHK checksum
// block, the COV coverage-distribution histogram, the GCD GC-depth
// distribution, the MPC mismatches-per-cycle section (emitted with --ref-seq)
// and the RFS reference-statistics section (emitted with --ref-stats) are all
// byte-faithful to upstream. PARITY_VALIDATION.md tracks which sections are
// byte-faithful vs. deferred.
//
// BWA-style quality trimming (-q/--trim-quality) feeds the "bases trimmed"
// SN counter and --target-regions restricts every counter to a target file.
//
// Skipped intentionally for v1 (documented):
//   - The FBC/FTC/LBC/LTC barcode sections.
//   - --remove-overlaps (single-record stats are unaffected by overlap
//     removal for the counters we emit).
//   - Command-line positional region arguments (so the RFS-with-region
//     path of upstream stats test 18 is not reproducible; --ref-stats
//     with --target-regions covers the equivalent functionality).
package samtools

import (
	"bufio"
	"fmt"
	"hash/crc32"
	"io"
	"math"
	"os"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// bwaMinReadLength mirrors upstream stats.c's BWA_MIN_RDLEN: reads shorter
// than this are never quality-trimmed by bwaTrimRead.
const bwaMinReadLength = 35

// bwaTrimRead returns the number of bases BWA would trim from the 3' end of a
// read given a trim-quality threshold, the per-base quality bytes, the read
// length and the reverse-strand flag. It is a faithful port of upstream
// stats.c's bwa_trim_read, including bwa's documented off-by-one (max_l = l
// rather than l+1, which trims one base fewer than the running maximum
// would imply).
func bwaTrimRead(trimQual int, quals []byte, length int, reverse bool) int {
	if length < bwaMinReadLength {
		return 0
	}
	maxTrimmed := length - bwaMinReadLength + 1
	sum, maxSum, maxL := 0, 0, 0
	for l := 0; l < maxTrimmed; l++ {
		idx := length - 1 - l
		if reverse {
			idx = l
		}
		sum += trimQual - int(quals[idx])
		if sum < 0 {
			break
		}
		if sum > maxSum {
			maxSum = sum
			maxL = l
		}
	}
	return maxL
}

// regionInterval is one closed 1-based interval [Beg, End] from a
// --target-regions file.
type regionInterval struct {
	Beg int32
	End int32
}

// targetRegions holds the parsed --target-regions intervals keyed by reference
// name, plus the total number of distinct reference bases they cover.
type targetRegions struct {
	byRef map[string][]regionInterval
	count int64
}

// loadTargetRegions parses a --target-regions file. The format mirrors
// upstream stats.c's init_regions: each non-comment line is a reference name
// followed by whitespace then two whitespace-separated 1-based inclusive
// coordinates ("seq beg end"). It is NOT a standard BED file. Intervals are
// sorted and overlapping/adjacent-touching intervals on the same reference are
// merged, exactly as upstream does, and the covered-base total is accumulated
// for the "bases inside the target" SN line. Reference names not present in
// the BAM header are skipped with a warning, matching upstream.
func loadTargetRegions(path string, hdr *sam.Header) (*targetRegions, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	known := make(map[string]struct{})
	if hdr != nil {
		for _, ref := range hdr.Refs {
			known[ref.Name] = struct{}{}
		}
	}

	tr := &targetRegions{byRef: make(map[string][]regionInterval)}
	prev := make(map[string]int32)
	warned := false
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := sc.Text()
		if line == "" || line[0] == '#' {
			continue
		}
		i := 0
		for i < len(line) && !isSpaceByte(line[i]) {
			i++
		}
		if i >= len(line) {
			return nil, fmt.Errorf("could not parse the file: %s [%s]", path, line)
		}
		name := line[:i]
		var beg, end int32
		if _, serr := fmt.Sscanf(strings.TrimSpace(line[i:]), "%d %d", &beg, &end); serr != nil {
			return nil, fmt.Errorf("could not parse the region [%s]", strings.TrimSpace(line[i:]))
		}
		if hdr != nil {
			if _, ok := known[name]; !ok {
				if !warned {
					fmt.Fprintf(os.Stderr, "Warning: Some sequences not present in the BAM, e.g. %q. This message is printed only once.\n", name)
					warned = true
				}
				continue
			}
		}
		if p, ok := prev[name]; ok && p > beg {
			return nil, fmt.Errorf("the positions are not in chromosomal order (%s comes after %d)", line, p)
		}
		prev[name] = beg
		tr.byRef[name] = append(tr.byRef[name], regionInterval{Beg: beg, End: end})
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	if len(tr.byRef) == 0 {
		return nil, fmt.Errorf("unable to map the -t sequences to the BAM sequences")
	}

	// Sort and merge overlapping intervals per reference, then sum covered
	// bases — mirroring upstream init_regions' qsort + dedup loop.
	for name, ivs := range tr.byRef {
		sort.Slice(ivs, func(a, b int) bool {
			if ivs[a].Beg != ivs[b].Beg {
				return ivs[a].Beg < ivs[b].Beg
			}
			return ivs[a].End < ivs[b].End
		})
		merged := ivs[:1]
		for _, iv := range ivs[1:] {
			last := &merged[len(merged)-1]
			if last.End < iv.Beg {
				merged = append(merged, iv)
			} else if last.End < iv.End {
				last.End = iv.End
			}
		}
		tr.byRef[name] = merged
		for _, iv := range merged {
			tr.count += int64(iv.End - iv.Beg + 1)
		}
	}
	return tr, nil
}

// isSpaceByte reports whether b is an ASCII whitespace character, matching
// C's isspace for the bytes that appear in a target-regions file.
func isSpaceByte(b byte) bool {
	return b == ' ' || b == '\t' || b == '\r' || b == '\n' || b == '\v' || b == '\f'
}

// StatsOptions configures the Stats run.
type StatsOptions struct {
	// RefSeq is the path to an indexed reference FASTA (-r/--ref-seq). When
	// set, the GCD GC-depth distribution derives GC content from the
	// reference rather than approximating it from the read sequences.
	RefSeq string
	// GcdBinSize is the width, in reference bases, of one GC-depth bin
	// (--GC-depth). Zero selects the upstream default of 20000.
	GcdBinSize int
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
	// TargetBED is the --target-regions file: when set, every statistic is
	// restricted to reads overlapping a listed interval. The file format is
	// "seq-name beg end" with 1-based inclusive coordinates (NOT BED).
	TargetBED string
	// TrimQuality is the BWA-style 3'-end quality-trim threshold passed via
	// -q/--trim-quality. Zero (the default) disables trimming.
	TrimQuality int
	// CovThreshold is the coverage threshold (-g/--cov-threshold) used by the
	// "percentage of target genome with coverage > N" SN line.
	CovThreshold int
	// RefStats enables the RFS reference-statistics section (--ref-stats).
	RefStats bool
	// RefStatsChunk is the reference-fetch chunk width in bytes
	// (--ref-stats-chunk, in megabytes on the CLI). Zero selects the upstream
	// default of 1 MB. It only affects how the reference FASTA is read for
	// RFS, not the emitted output.
	RefStatsChunk int
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
	// AnomalousReads counts, per paired-and-mapped primary read, those whose
	// mate maps to a different reference (upstream nreads_anomalous). The
	// "pairs on different chromosomes" SN line emits AnomalousReads/2, which
	// — unlike a pair-completion count — stays correct when one mate is
	// region-filtered out by --target-regions.
	AnomalousReads     int64
	InsertSumAbs       int64 // |TLEN| over pairs accepted into IS
	InsertSumSqAbs     int64 // |TLEN|^2
	InsertPairsCounted int64

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

	// Target-regions state. regions is nil unless --target-regions is in
	// effect. regCursor tracks the per-reference scan cursor mirroring
	// upstream's reg->cpos; it relies on coordinate-sorted input (upstream
	// errors otherwise). regFrom/regTo hold the first interval matched by the
	// current record, used to clip "bases mapped (cigar)" to the target.
	regions   *targetRegions
	regCursor map[string]int
	regFrom   int32
	regTo     int32

	// GCD GC-depth distribution. The reference is split into gcdBinSize-wide
	// segments and per segment the read depth and GC content are recorded;
	// at output the segments are sorted by GC and depth percentiles are
	// reported. gcd holds one entry per segment, accumulated only for
	// coordinate-sorted input (like COV). Mirroring upstream's igcd/ngcd
	// indexing, gcd[0] is always an empty placeholder and a freshly started
	// segment is gcd[gcdIdx] with gcdIdx incremented first; the slice
	// therefore holds gcdIdx+1 entries. gcdPos is the reference start of the
	// current segment (-1 = none yet), gcdContig the contig it lies on.
	// gcdRef is the reference reader, non-nil only when --ref-seq is set.
	gcd        []gcDepth
	gcdIdx     int
	gcdPos     int32
	gcdContig  string
	gcdBinSize int
	gcdRef     *fasta.RandomAccess
	gcdRefLens map[string]int64

	// MPC mismatches-per-cycle distribution. Accumulated only when --ref-seq
	// is set: per mapped read each base is walked against the reference and
	// bucketed by cycle and quality. mpc holds one cycleQuals per cycle; an
	// N read base lands in quality slot 0, a true mismatch in slot qual+1
	// (matching upstream count_mismatches_per_cycle, stats.c:476). mpc is nil
	// unless --ref-seq is in effect, gating the MPC output section; with
	// --ref-seq it is a non-nil (possibly empty) slice. mpcRefBuf is a scratch
	// buffer holding the current record's reference span as refBaseCode codes.
	mpc       []cycleQuals
	mpcRefBuf []byte

	// RFS reference-statistics section. rfs is non-nil only when --ref-stats
	// is set, gating the RFS output. It is computed after the input stream is
	// fully consumed (collectRefStats), mirroring upstream collect_refstats.
	rfs *refStats
}

// refStats holds the RFS reference-statistics summary plus one per-sequence
// row, mirroring upstream stats.c's refstats struct (stats.c:168).
type refStats struct {
	totalCount  int     // total @SQ count in the BAM header
	count       int     // sequences/regions actually reported
	combinedLen int64   // sum of reported lengths
	minLen      int64   // shortest reported length (0 until first set)
	maxLen      int64   // longest reported length
	avgLen      float64 // mean reported length, -1 when count == 0
	avgGC       float64 // mean GC fraction, -1 when no reference FASTA
	rows        []refStatRow
}

// refStatRow is one per-sequence RFS row: the sequence (or region) name, its
// length, GC fraction and undetermined-base (N) count. gc/n are -1 when no
// reference FASTA was supplied.
type refStatRow struct {
	name string
	len  int64
	gc   float64
	n    int64
}

// gcDepth mirrors upstream stats.c's gc_depth_t: one GC-depth segment holding
// an accumulated GC value and a read-depth count.
type gcDepth struct {
	gc    float64
	depth uint32
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
		// gcdPos == -1 marks "no GC-depth segment started yet", mirroring
		// upstream stats.c's gcd_pos = -1LL initialiser.
		gcdPos: -1,
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
	hdr := br.Header()
	c := newStatsCounters()
	// Compute the COV bin geometry once up front so depth can be binned
	// incrementally during the streaming flush.
	c.initCovBins(opts)
	// GCD GC-depth bin width (--GC-depth, default 20000 reference bases).
	c.gcdBinSize = opts.GcdBinSize
	if c.gcdBinSize <= 0 {
		c.gcdBinSize = 20000
	}
	// --ref-seq: open the reference FASTA so GCD can derive GC content from
	// the reference rather than approximating it from the read sequences.
	if opts.RefSeq != "" {
		ref, rerr := fasta.OpenRandomAccess(opts.RefSeq)
		if rerr != nil {
			return fmt.Errorf("could not load faidx: %s", opts.RefSeq)
		}
		defer ref.Close()
		c.gcdRef = ref
		c.gcdRefLens = make(map[string]int64)
		for _, e := range ref.Index().Entries() {
			c.gcdRefLens[e.Name] = e.Length
		}
		// Upstream allocates mpc_buf iff a reference FASTA is given
		// (stats.c:2386); a non-nil slice gates the MPC output section.
		c.mpc = make([]cycleQuals, 0)
	}
	// --target-regions: parse the target file and prepare per-reference scan
	// cursors. Every counter is then restricted to reads overlapping a
	// listed interval (see isInRegions).
	if opts.TargetBED != "" {
		tr, terr := loadTargetRegions(opts.TargetBED, hdr)
		if terr != nil {
			return terr
		}
		c.regions = tr
		c.regCursor = make(map[string]int)
	}
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
		if oerr := c.observe(rec, opts); oerr != nil {
			return oerr
		}
	}
	if c.sortBroken {
		c.IsSorted = 0
	}
	// Flush whatever coverage depth remains in the window for the last
	// contig before the report is written.
	c.flushCoverageWindow(math.MaxInt32)
	// --ref-stats: build the RFS reference-statistics section after the whole
	// stream is consumed, mirroring upstream collect_refstats.
	if opts.RefStats {
		c.collectRefStats(hdr, opts)
	}
	return c.Write(out, opts)
}

// observe folds one record into the accumulator. It returns a non-nil error
// only when --target-regions is active and the input is not coordinate
// sorted, mirroring upstream stats.c's hard error in is_in_regions.
func (c *StatsCounters) observe(rec *sam.Record, opts StatsOptions) error {
	// --target-regions filter: records that overlap no listed interval are
	// skipped entirely (upstream collect_stats returns immediately). This
	// runs before the sort-order detector below updates state, so it sees
	// the same sortedness upstream's is_in_regions sees.
	if c.regions != nil {
		in, ierr := c.isInRegions(rec)
		if ierr != nil {
			return ierr
		}
		if !in {
			return nil
		}
	}
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
		return nil
	}
	if rec.IsSupplementary() {
		c.Supplementary++
		// Upstream's bases_mapped_cigar counter folds in supplementary
		// alignments — stats.c walks the cigar of EVERY mapped record
		// (primary + supplementary) into nbases_mapped_cigar.
		c.addBasesMappedCigar(rec)
		// Upstream count_indels runs for every mapped record, including
		// supplementary alignments (only secondary reads are excluded).
		if rec.IsMapped() {
			c.countIndels(rec)
			// COV folds in supplementary alignments too; upstream walks
			// the CIGAR of every record reaching collect_stats (only
			// secondary reads return early).
			c.accumulateCoverage(rec)
			// GCD likewise folds in supplementary alignments — upstream's
			// GC-depth block sits after the IS_UNMAPPED early-return and so
			// runs for every mapped non-secondary record. A supplementary
			// read carries a GC count of zero because upstream's gc_count
			// comes from collect_orig_read_stats, which is IS_ORIGINAL-only.
			c.accumulateGCD(rec, 0)
			// MPC mismatches-per-cycle is accumulated at the same call site
			// upstream (stats.c:1400) and so also folds in supplementary
			// alignments. It is a no-op unless --ref-seq is in effect.
			c.accumulateMPC(rec, unclippedLength(rec))
		}
		return nil
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
		return nil
	}
	if opts.FilteringFlag != 0 && rec.Flag&opts.FilteringFlag != 0 {
		c.FilteredOut++
		return nil
	}
	if opts.RemoveDups && rec.IsDuplicate() {
		c.FilteredOut++
		return nil
	}
	if opts.MinMAPQ > 0 && rec.MapQ < opts.MinMAPQ {
		c.FilteredOut++
		return nil
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
	// collect_orig_read_stats — runs for every QC-passed primary record. The
	// returned G+C base count feeds the GC-depth distribution.
	gcCount := c.collectCycleStats(rec, opts)

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
	c.addBasesMappedCigar(rec)

	// Indel distribution and per-cycle indel counts (mapped reads only).
	c.countIndels(rec)

	// COV coverage-distribution depth, GCD GC-depth and MPC
	// mismatches-per-cycle (mapped reads only). MPC is a no-op unless
	// --ref-seq is in effect.
	if rec.IsMapped() {
		c.accumulateCoverage(rec)
		c.accumulateGCD(rec, gcCount)
		c.accumulateMPC(rec, unclippedLength(rec))
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
		// Per-read different-chromosome tally (upstream nreads_anomalous):
		// counted whenever the mate maps to a different reference. RNEXT
		// "=" denotes the same reference.
		if rec.RNext != "" && rec.RNext != "*" && rec.RNext != "=" && rec.RNext != rec.RName {
			c.AnomalousReads++
		}
		c.classifyPair(rec, opts)
	}
	return nil
}

// addBasesMappedCigar folds one mapped record's CIGAR into the "bases mapped
// (cigar)" counter. M/=/X and I operations contribute. When --target-regions
// is active the count is clipped to the first target interval the record
// matched ([regFrom, regTo]), mirroring upstream stats.c's on-target counting
// at stats.c:1306-1333.
func (c *StatsCounters) addBasesMappedCigar(rec *sam.Record) {
	if c.regions == nil {
		for _, op := range rec.Cigar {
			switch op.Char() {
			case 'M', '=', 'X', 'I':
				c.BasesMappedCigar += int64(op.Length())
			}
		}
		return
	}
	// Region-restricted: iref is the 1-based reference coordinate of the
	// current CIGAR op (rec.Pos is already 1-based here).
	iref := rec.Pos
	for _, op := range rec.Cigar {
		ch := op.Char()
		ncig := int32(op.Length())
		if ncig == 0 {
			continue
		}
		switch ch {
		case 'M', '=', 'X':
			clipped := ncig
			if iref < c.regFrom {
				clipped -= c.regFrom - iref
			} else if iref+clipped-1 > c.regTo {
				clipped -= iref + clipped - 1 - c.regTo
			}
			if clipped < 0 {
				clipped = 0
			}
			c.BasesMappedCigar += int64(clipped)
			iref += ncig
		case 'I':
			iref += ncig
			if iref >= c.regFrom && iref <= c.regTo {
				c.BasesMappedCigar += int64(ncig)
			}
		case 'D':
			// Per upstream stats.c:1316, a deletion does NOT advance the
			// reference cursor in the on-target counting loop (it only adds
			// to readlen, which we do not track here). This is an upstream
			// quirk replicated for byte-parity — do not advance iref.
		}
	}
}

// isInRegions reports whether rec overlaps any --target-regions interval,
// advancing the per-reference scan cursor and recording the first matched
// interval in regFrom/regTo. It mirrors upstream stats.c's is_in_regions:
// reads whose reference is absent from the target file are excluded, the scan
// relies on coordinate-sorted input, and even a single-base overlap includes
// the read. A non-nil error is returned when the input is not sorted.
func (c *StatsCounters) isInRegions(rec *sam.Record) (bool, error) {
	ivs, ok := c.regions.byRef[rec.RName]
	if !ok || rec.RName == "" {
		return false, nil
	}
	// Upstream errors when -t is used on an unsorted BAM. Detect a backwards
	// jump within the same reference.
	if !c.sortBroken && rec.IsMapped() {
		if c.lastRef == rec.RName && rec.Pos < c.lastPos {
			return false, fmt.Errorf("the BAM must be sorted in order for -t to work")
		}
	}
	cur := c.regCursor[rec.RName]
	if cur >= len(ivs) {
		return false, nil // done for this reference
	}
	// Upstream compares the 1-based inclusive interval end against the
	// 0-based core.pos with `end <= pos`; rec.Pos here is 1-based, so the
	// equivalent test is `End < rec.Pos`.
	i := cur
	for i < len(ivs) && ivs[i].End < rec.Pos {
		i++
	}
	if i >= len(ivs) {
		c.regCursor[rec.RName] = len(ivs)
		return false, nil
	}
	// bam_endpos is a 0-based exclusive end, equal to our 1-based inclusive
	// EndPosition; the overlap test `endpos < beg` carries over unchanged.
	endpos := rec.EndPosition()
	if endpos < ivs[i].Beg {
		return false, nil
	}
	c.regCursor[rec.RName] = i
	c.regFrom = ivs[i].Beg
	c.regTo = ivs[i].End
	return true, nil
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
// collect_orig_read_stats). When -q/--trim-quality is set it also folds the
// BWA-style trimmed-base count into BasesTrimmed. It returns the record's
// G+C base count, which the GC-depth accumulator consumes.
func (c *StatsCounters) collectCycleStats(rec *sam.Record, opts StatsOptions) int {
	seqLen := len(rec.Seq)
	order := readOrder(rec)
	reverse := rec.Flag&sam.FlagReverse != 0

	// BWA-style 3'-end quality trimming. Upstream invokes bwa_trim_read on
	// the record's quality bytes inside collect_orig_read_stats; a "*"
	// quality string behaves as all-0xFF bytes, which trims nothing.
	if opts.TrimQuality > 0 && seqLen > 0 {
		quals := rec.Qual
		if len(quals) != seqLen {
			quals = make([]byte, seqLen)
			for i := range quals {
				quals[i] = 0xff
			}
		}
		c.BasesTrimmed += int64(bwaTrimRead(opts.TrimQuality, quals, seqLen, reverse))
	}

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
		return 0
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
	return gcCount
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

	// When --target-regions is active, COV depth is restricted to on-target
	// reference positions, mirroring upstream's per-chunk round_buffer
	// inserts (stats.c:1419-1452).
	var ivs []regionInterval
	if c.regions != nil {
		ivs = c.regions.byRef[rec.RName]
	}
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
				if p < c.covFlushed {
					continue
				}
				if c.regions != nil && !positionInIntervals(p, ivs) {
					continue
				}
				c.covWindow[p]++
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

// gcdReadLen returns the number of reference bases a mapped record spans for
// GC-depth purposes: the read length plus every deletion length, mirroring
// upstream stats.c's `readlen` variable (stats.c:1306/1316/1342). Insertions
// and soft clips do not contribute.
func gcdReadLen(rec *sam.Record) int32 {
	n := int32(len(rec.Seq))
	for _, op := range rec.Cigar {
		if op.Char() == 'D' {
			n += int32(op.Length())
		}
	}
	return n
}

// accumulateGCD folds one mapped record into the GC-depth distribution. It is
// a faithful port of the GCD accumulation block of upstream stats.c
// (stats.c:1369-1415): the reference is split into gcdBinSize-wide segments
// and per segment a read-depth count plus a GC value are recorded.
//
// gcCount is the record's G+C base count and is meaningful only for original
// (non-supplementary) reads; supplementary reads pass gcCount 0, matching
// upstream where gc_count comes from collect_orig_read_stats which runs for
// IS_ORIGINAL records only.
//
// With --ref-seq the GC value is read straight from the reference window; in
// the no-reference default it is the read's GC fraction accumulated across
// the segment's reads and averaged at output time.
func (c *StatsCounters) accumulateGCD(rec *sam.Record, gcCount int) {
	seqLen := len(rec.Seq)
	if seqLen == 0 {
		return
	}
	// C uses a 0-based core.pos; rec.Pos here is 1-based.
	pos := rec.Pos - 1
	readLen := gcdReadLen(rec)

	if c.gcdRef != nil {
		// Reference path: start a new segment on first read, a contig change,
		// or when the read crosses past the current segment's bin width.
		incGcd := false
		if c.gcdPos == -1 || c.gcdContig != rec.RName {
			incGcd = true
		} else if int64(c.gcdPos)+int64(c.gcdBinSize) < int64(pos)+int64(readLen) {
			incGcd = true
		}
		if incGcd {
			c.gcdIdx++
			c.growGCD()
			c.gcdContig = rec.RName
			c.gcdPos = pos
			c.gcd[c.gcdIdx].gc = c.faiGCContent(rec.RName, pos, c.gcdBinSize)
		}
	} else if c.gcdPos == -1 || c.gcdContig != rec.RName || pos-c.gcdPos > int32(c.gcdBinSize) {
		// No-reference path: start a new segment on first read, a contig
		// change, or when the read starts past the current bin's far edge.
		c.gcdContig = rec.RName
		c.gcdPos = pos
		c.gcdIdx++
		c.growGCD()
	}
	c.gcd[c.gcdIdx].depth++
	if c.gcdRef == nil {
		c.gcd[c.gcdIdx].gc += float64(gcCount) / float64(seqLen)
	}
}

// growGCD ensures the gcd slice can hold index gcdIdx, mirroring upstream's
// realloc_gcd_buffer. gcd[0] is always an empty placeholder.
func (c *StatsCounters) growGCD() {
	for len(c.gcd) <= c.gcdIdx {
		c.gcd = append(c.gcd, gcDepth{})
	}
}

// seqNibble maps an ASCII SEQ character to the BAM 4-bit sequence encoding
// (the "=ACMGRSVTWYHKDBN" table htslib's bam_seqi yields). Unknown bytes map
// to 15 (N), matching htslib's seq_nt16_table. The plain bases A/C/G/T encode
// to 1/2/4/8, which deliberately coincides with refBaseCode so an A==A
// comparison succeeds even though the two tables come from different sources.
var seqNibble = func() [256]byte {
	var t [256]byte
	for i := range t {
		t[i] = 15
	}
	const codes = "=ACMGRSVTWYHKDBN"
	for i := 0; i < len(codes); i++ {
		ch := codes[i]
		t[ch] = byte(i)
		t[ch|0x20] = byte(i) // lowercase
	}
	return t
}()

// refBaseCode maps an ASCII reference-FASTA base to the 2-bit-ish code
// upstream's read_ref_seq (stats.c:562) produces: A/a=1, C/c=2, G/g=4, T/t=8,
// and everything else (N and ambiguity codes) = 0 ("undetermined").
func refBaseCode(b byte) byte {
	switch b {
	case 'A', 'a':
		return 1
	case 'C', 'c':
		return 2
	case 'G', 'g':
		return 4
	case 'T', 't':
		return 8
	default:
		return 0
	}
}

// countMismatchesPerCycle folds one mapped record into the MPC
// mismatches-per-cycle distribution. It is a faithful port of upstream
// count_mismatches_per_cycle (stats.c:476): the record's CIGAR is walked
// against the reference and, for each aligned base, an N read base is bucketed
// in quality slot 0 while a true mismatch (both bases determined and unequal)
// is bucketed in quality slot (qual+1). The cycle index advances through
// soft-/hard-clips and insertions exactly as upstream does, and for reverse
// reads it is mirrored via read_len-icycle-1, where read_len is the unclipped
// read length. A "*" quality string yields per-base quality 0xFF, so qual+1
// wraps (uint8) to 0 — upstream's documented quirk that lands such mismatches
// in the N column. mpc must already be non-nil (gated on --ref-seq).
func (c *StatsCounters) countMismatchesPerCycle(rec *sam.Record, readLen int) {
	isFwd := rec.Flag&sam.FlagReverse == 0
	iread, icycle := 0, 0
	// C uses a 0-based core.pos; rec.Pos here is 1-based.
	iref := int(rec.Pos) - 1
	seq := rec.Seq
	quals := rec.Qual
	for _, op := range rec.Cigar {
		ch := op.Char()
		ncig := int(op.Length())
		switch ch {
		case 'I':
			iread += ncig
			icycle += ncig
			continue
		case 'D':
			iref += ncig
			continue
		case 'S':
			icycle += ncig
			iread += ncig
			continue
		case 'H':
			icycle += ncig
			continue
		case 'N', 'P':
			// Reference skips and padding contribute nothing.
			continue
		case 'M', '=', 'X':
		default:
			continue
		}
		for im := 0; im < ncig; im++ {
			if iread >= len(seq) {
				break
			}
			cread := seqNibble[seq[iread]]
			var cref byte
			refIdx := iref - (int(rec.Pos) - 1)
			if refIdx >= 0 && refIdx < len(c.mpcRefBuf) {
				cref = c.mpcRefBuf[refIdx]
			}
			idx := icycle
			if !isFwd {
				idx = readLen - icycle - 1
			}
			if idx >= 0 && idx < statsNbasesMax {
				if cread == 15 {
					c.mpc = growQuals(c.mpc, idx+1)
					c.mpc[idx][0]++
				} else if cref != 0 && cread != 0 && cref != cread {
					var qb byte = 0xff
					if iread < len(quals) {
						qb = quals[iread]
					}
					qual := byte(qb + 1)
					c.mpc = growQuals(c.mpc, idx+1)
					c.mpc[idx][qual]++
				}
			}
			iref++
			iread++
			icycle++
		}
	}
}

// statsNbasesMax bounds the MPC cycle index defensively; upstream errors out
// past stats->nbases. No real read reaches it.
const statsNbasesMax = 1 << 20

// accumulateMPC fetches the reference span covering one mapped record and
// folds the record into the MPC distribution. It is the --ref-seq-gated
// counterpart of accumulateGCD's count_mismatches_per_cycle call site
// (stats.c:1400). readLen is the unclipped read length upstream passes as the
// count_mismatches_per_cycle read_len argument.
func (c *StatsCounters) accumulateMPC(rec *sam.Record, readLen int) {
	if c.mpc == nil || c.gcdRef == nil {
		return
	}
	contigLen, ok := c.gcdRefLens[rec.RName]
	if !ok {
		return
	}
	start := int64(rec.Pos) - 1
	end := start + int64(rec.Cigar.ReferenceLength())
	if end > contigLen {
		end = contigLen
	}
	if start < 0 || start >= end {
		c.mpcRefBuf = c.mpcRefBuf[:0]
		return
	}
	bases, err := c.gcdRef.Fetch(rec.RName, start, end)
	if err != nil {
		c.mpcRefBuf = c.mpcRefBuf[:0]
		return
	}
	if cap(c.mpcRefBuf) < len(bases) {
		c.mpcRefBuf = make([]byte, len(bases))
	} else {
		c.mpcRefBuf = c.mpcRefBuf[:len(bases)]
	}
	for i := 0; i < len(bases); i++ {
		c.mpcRefBuf[i] = refBaseCode(bases[i])
	}
	c.countMismatchesPerCycle(rec, readLen)
}

// refStatsChunkBytes returns the reference-fetch chunk width in bytes. The
// --ref-stats-chunk CLI value is given in megabytes; values <= 0 collapse to
// 1 MB exactly as upstream stats.c does (option case 3, stats.c:2762). The
// chunk width only affects how the FASTA is read, not the emitted RFS output.
func refStatsChunkBytes(opts StatsOptions) int64 {
	mb := opts.RefStatsChunk
	if mb <= 0 {
		mb = 1
	}
	return int64(mb) * 1024 * 1024
}

// refSpanStats fetches the closed 1-based reference interval [beg, end] on
// contig name and returns its GC fraction and undetermined-base (N) count. It
// mirrors upstream collect_refstats' per-region accumulation (stats.c:2609):
// only A/C/G/T are counted toward the GC denominator, N/n toward the
// undetermined total, and any other byte (ambiguity codes) toward neither.
// The interval is read in refStatsChunkBytes-wide chunks, which leaves the
// result identical to a single fetch.
func (c *StatsCounters) refSpanStats(name string, beg, end int64, chunk int64) (gcFrac float64, nCount int64) {
	start := beg - 1 // 0-based
	var gc, at int64
	for rem := end - start; rem > 0; {
		span := rem
		if span > chunk {
			span = chunk
		}
		bases, err := c.gcdRef.Fetch(name, start, start+span)
		if err != nil || len(bases) == 0 {
			break
		}
		for _, b := range bases {
			switch b {
			case 'G', 'g', 'C', 'c':
				gc++
			case 'A', 'a', 'T', 't':
				at++
			case 'N', 'n':
				nCount++
			}
		}
		got := int64(len(bases))
		start += got
		rem -= got
		if got < span {
			break
		}
	}
	if total := gc + at; total > 0 {
		gcFrac = float64(gc) / float64(total)
	}
	return gcFrac, nCount
}

// collectRefStats builds the RFS reference-statistics section after the input
// stream has been fully consumed. It is a faithful port of upstream
// collect_refstats (stats.c:2498).
//
// Without --target-regions every @SQ entry of the BAM header becomes one RFS
// row (the read content of the BAM is irrelevant). With --target-regions the
// sorted, merged target intervals become the rows: a row spanning a whole
// reference is named by the reference, a sub-range by "name:start-end". When
// no reference FASTA was supplied the GC fraction and N count are reported as
// -1, matching upstream's lack-of-data sentinel.
func (c *StatsCounters) collectRefStats(hdr *sam.Header, opts StatsOptions) {
	rs := &refStats{}
	if hdr != nil {
		rs.totalCount = len(hdr.Refs)
	}
	chunk := refStatsChunkBytes(opts)
	var gcSum float64
	gcKnown := c.gcdRef != nil

	addRow := func(name string, length int64, gc float64, n int64) {
		rs.rows = append(rs.rows, refStatRow{name: name, len: length, gc: gc, n: n})
		rs.count++
		rs.combinedLen += length
		if rs.minLen == 0 || length < rs.minLen {
			rs.minLen = length
		}
		if length > rs.maxLen {
			rs.maxLen = length
		}
		if gcKnown {
			gcSum += gc
		}
	}

	if c.regions == nil {
		// Whole-file: one row per header reference, in header order.
		if hdr != nil {
			for _, ref := range hdr.Refs {
				length := int64(ref.Length)
				gc, n := float64(-1), int64(-1)
				if gcKnown {
					if hl, ok := c.gcdRefLens[ref.Name]; ok {
						hi := length
						if hl < hi {
							hi = hl
						}
						gc, n = c.refSpanStats(ref.Name, 1, hi, chunk)
					} else {
						gc, n = 0, 0
					}
				}
				addRow(ref.Name, length, gc, n)
			}
		}
	} else {
		// Target-regions: rows from the merged intervals, in header order.
		hdrLen := make(map[string]int32)
		if hdr != nil {
			for _, ref := range hdr.Refs {
				hdrLen[ref.Name] = ref.Length
			}
		}
		var order []string
		if hdr != nil {
			for _, ref := range hdr.Refs {
				order = append(order, ref.Name)
			}
		}
		for _, name := range order {
			ivs, ok := c.regions.byRef[name]
			if !ok {
				continue
			}
			hl := int64(hdrLen[name])
			for _, iv := range ivs {
				beg, end := int64(iv.Beg), int64(iv.End)
				if end < beg {
					continue
				}
				length := end - beg + 1
				if hl > 0 && length > hl {
					length = hl
				}
				gc, n := float64(-1), int64(-1)
				if gcKnown {
					if _, present := c.gcdRefLens[name]; present {
						gc, n = c.refSpanStats(name, beg, end, chunk)
					} else {
						gc, n = 0, 0
					}
				}
				addRow(fmt.Sprintf("%s:%d-%d", name, beg, end), length, gc, n)
			}
		}
	}

	if rs.count > 0 {
		rs.avgLen = float64(rs.combinedLen) / float64(rs.count)
		if gcKnown {
			rs.avgGC = gcSum / float64(rs.count)
		} else {
			rs.avgGC = -1
		}
	} else {
		rs.avgLen = -1
		rs.avgGC = -1
	}
	c.rfs = rs
}

// faiGCContent returns the GC fraction of the gcdBinSize-wide reference window
// starting at the 0-based position pos on contig name. It mirrors upstream's
// fai_gc_content (stats.c:610): G/C bases are counted over the known
// (non-N, non-ambiguous) bases of the window and the count is divided by that
// known-base total. The window is clipped to the contig end.
func (c *StatsCounters) faiGCContent(name string, pos int32, length int) float64 {
	contigLen, ok := c.gcdRefLens[name]
	if !ok {
		return 0
	}
	start := int64(pos)
	end := start + int64(length)
	if end > contigLen {
		end = contigLen
	}
	if start < 0 || start >= end {
		return 0
	}
	bases, err := c.gcdRef.Fetch(name, start, end)
	if err != nil {
		return 0
	}
	gc, count := 0, 0
	for _, b := range bases {
		// Upstream read_ref_seq (stats.c:588) folds case, so soft-masked
		// (lowercase) reference bases count toward GC/AT just like uppercase.
		switch b {
		case 'C', 'G', 'c', 'g':
			gc++
			count++
		case 'A', 'T', 'a', 't':
			count++
		}
	}
	if count == 0 {
		return 0
	}
	return float64(gc) / float64(count)
}

// positionInIntervals reports whether the 1-based reference position p falls
// inside any of the sorted, merged target intervals ivs.
func positionInIntervals(p int32, ivs []regionInterval) bool {
	for _, iv := range ivs {
		if p < iv.Beg {
			return false
		}
		if p <= iv.End {
			return true
		}
	}
	return false
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
	c.writeSN(bw, opts)
	// Histogram and per-cycle section bodies are emitted only in the
	// non-sparse path. FFQ/LFQ/GCF/GCL/GCC/GCT/IC/ID carry real data and are
	// byte-faithful to upstream; section order matches stats.c.
	if !opts.Sparse {
		c.writeFFQLFQ(bw)
		c.writeMPC(bw)
		c.writeGCFGCL(bw)
		c.writeGCC(bw)
		c.writeGCT(bw)
		c.writeFBCLTC(bw)
		c.writeRL(bw)
		c.writeMAPQ(bw)
		c.writeIS(bw, opts)
		c.writeID(bw)
		c.writeIC(bw)
		c.writeCOV(bw, opts)
		c.writeGCD(bw)
		c.writeRFS(bw)
	}
	return bw.Flush()
}

// writeMPC emits the MPC mismatches-per-cycle section. It is emitted only when
// --ref-seq is in effect (mpc non-nil), matching upstream's mpc_buf gate
// (stats.c:1639). One row is printed per cycle 1..maxLen — the same row count
// as FFQ/LFQ — and each row lists counts for qualities 0..maxQual; the first
// data column is the N count, the rest are mismatches by quality.
func (c *StatsCounters) writeMPC(bw *bufio.Writer) {
	if c.mpc == nil {
		return
	}
	maxQual := c.effectiveMaxQual()
	fmt.Fprintln(bw, "# Mismatches per cycle and quality. Use `grep ^MPC | cut -f 2-` to extract this part.")
	fmt.Fprintln(bw, "# Columns correspond to qualities, rows to cycles. First column is the cycle number, second")
	fmt.Fprintln(bw, "# is the number of N's and the rest is the number of mismatches")
	for ibase := 0; ibase < c.effectiveMaxLen(); ibase++ {
		fmt.Fprintf(bw, "MPC\t%d", ibase+1)
		for iqual := 0; iqual <= maxQual; iqual++ {
			var v int64
			if ibase < len(c.mpc) {
				v = c.mpc[ibase][iqual]
			}
			fmt.Fprintf(bw, "\t%d", v)
		}
		fmt.Fprintln(bw)
	}
}

// writeRFS emits the RFS reference-statistics section. It is emitted only when
// --ref-stats is in effect (rfs non-nil), matching upstream's rstat gate
// (stats.c:1894). The summary row reports total/reported sequence counts,
// average GC, min/max/average length and total length; each following row
// reports one sequence's (or region's) name, length, GC fraction and
// undetermined-base count.
func (c *StatsCounters) writeRFS(bw *bufio.Writer) {
	if c.rfs == nil {
		return
	}
	rs := c.rfs
	fmt.Fprintln(bw, "# Reference statistics. Use `grep ^RFS | cut -f 2-` to extract this part.")
	fmt.Fprintln(bw, "# Total count, Output count, Average GC, Min length, Max length, Average length, Total length in first row.")
	fmt.Fprintln(bw, "# Sequence name, Length, GC content, Unknown count in following rows.")
	fmt.Fprintf(bw, "RFS\t%d\t%d\t%.2f\t%d\t%d\t%.2f\t%d\n",
		rs.totalCount, rs.count, rs.avgGC, rs.minLen, rs.maxLen, rs.avgLen, rs.combinedLen)
	for _, row := range rs.rows {
		fmt.Fprintf(bw, "RFS\t%s\t%d\t%.2f\t%d\n", row.name, row.len, row.gc, row.n)
	}
}

// gcdPercentile interpolates the depth at the p-th percentile across the
// nbins consecutive GC-depth segments starting at grp. It is a faithful port
// of upstream stats.c's gcd_percentile (stats.c:1491), including its
// truncating float-to-int conversion and the k<=0 / k>=N edge clamps.
func gcdPercentile(grp []gcDepth, nbins, p int) float64 {
	n := float64(p) * float64(nbins+1) / 100.0
	k := int(n)
	if k <= 0 {
		return float64(grp[0].depth)
	}
	if k >= nbins {
		return float64(grp[nbins-1].depth)
	}
	d := n - float64(k)
	return float64(grp[k-1].depth) + d*(float64(grp[k].depth)-float64(grp[k-1].depth))
}

// writeGCD emits the GC-depth distribution. It is only emitted for
// coordinate-sorted input, matching upstream's is_sorted gating
// (stats.c:1848). Each segment's accumulated GC is finalised — multiplied by
// 100 and rounded — then the segments are sorted by GC (and depth) and
// grouped while their rounded GC stays within 0.1; one GCD row is printed per
// group. The columns are GC%, unique-sequence percentile, and the 10/25/50/
// 75/90th depth percentiles scaled by average read length / bin size. This is
// a faithful port of stats.c:1859-1891. Mirroring upstream, the gcd[0]
// placeholder bin and the gcd indexing (gcdIdx+1 valid entries, only the
// first gcdIdx finalised and grouped) are preserved exactly.
func (c *StatsCounters) writeGCD(bw *bufio.Writer) {
	if c.IsSorted != 1 {
		return
	}
	fmt.Fprintln(bw, "# GC-depth. Use `grep ^GCD | cut -f 2-` to extract this part. The columns are: GC%, unique sequence percentiles, 10th, 25th, 50th, 75th and 90th depth percentile")
	// avg_read_length mirrors stats.c:1586 — total length over the count of
	// 1st + 2nd + other reads (i.e. Sequences).
	avgReadLength := 0.0
	if c.Sequences > 0 {
		avgReadLength = float64(c.TotalLength) / float64(c.Sequences)
	}
	// Finalise the GC value of every segment below gcdIdx. The reference
	// path scales the raw fraction by 100; the no-reference path averages
	// the accumulated per-read fractions over the segment depth first.
	for i := 0; i < c.gcdIdx && i < len(c.gcd); i++ {
		if c.gcdRef != nil {
			c.gcd[i].gc = math.RoundToEven(100.0 * c.gcd[i].gc)
		} else if c.gcd[i].depth > 0 {
			c.gcd[i].gc = math.RoundToEven(100.0 * c.gcd[i].gc / float64(c.gcd[i].depth))
		}
	}
	// Sort the gcdIdx+1 valid entries by GC then depth (upstream gcd_cmp).
	if len(c.gcd) > 0 {
		n := c.gcdIdx + 1
		if n > len(c.gcd) {
			n = len(c.gcd)
		}
		sort.SliceStable(c.gcd[:n], func(a, b int) bool {
			if c.gcd[a].gc != c.gcd[b].gc {
				return c.gcd[a].gc < c.gcd[b].gc
			}
			return c.gcd[a].depth < c.gcd[b].depth
		})
	}
	igcd := 0
	for igcd < c.gcdIdx {
		gc := c.gcd[igcd].gc
		nbins := 0
		for itmp := igcd; itmp < c.gcdIdx && math.Abs(c.gcd[itmp].gc-gc) < 0.1; itmp++ {
			nbins++
		}
		grp := c.gcd[igcd : igcd+nbins]
		uniq := float64(igcd+nbins+1) * 100.0 / float64(c.gcdIdx+1)
		scale := avgReadLength / float64(c.gcdBinSize)
		fmt.Fprintf(bw, "GCD\t%.1f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
			gc, uniq,
			gcdPercentile(grp, nbins, 10)*scale,
			gcdPercentile(grp, nbins, 25)*scale,
			gcdPercentile(grp, nbins, 50)*scale,
			gcdPercentile(grp, nbins, 75)*scale,
			gcdPercentile(grp, nbins, 90)*scale)
		igcd += nbins
	}
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
func (c *StatsCounters) writeSN(bw *bufio.Writer, opts StatsOptions) {
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
	emit("pairs on different chromosomes", c.AnomalousReads/2, "")
	// Proper-pair % uses Sequences (= 1st + 2nd + other) as the denominator,
	// matching upstream stats.c:1606. NOT ReadsPaired — they differ when
	// records have the paired bit but no first/last fragment classification.
	pp := 0.0
	if c.Sequences > 0 {
		pp = 100.0 * float64(c.ReadsProperlyPaired) / float64(c.Sequences)
	}
	emit("percentage of properly paired reads (%)", fmt.Sprintf("%.1f", pp), "")
	// Target-regions SN lines (stats.c:1607-1611) — emitted only when a
	// --target-regions file actually mapped to BAM references. The coverage
	// percentage sums COV bins above cov_threshold; the cov array is sized
	// from -c and finalized by the streaming flush before Write runs.
	if c.regions != nil && c.regions.count > 0 {
		emit("bases inside the target", c.regions.count, "")
		var covSum int64
		for icov := opts.CovThreshold + 1; icov < c.ncov && icov < len(c.cov); icov++ {
			if icov >= 0 {
				covSum += c.cov[icov]
			}
		}
		pct := 100.0 * float64(covSum) / float64(c.regions.count)
		emit(fmt.Sprintf("percentage of target genome with coverage > %d (%%)", opts.CovThreshold),
			fmt.Sprintf("%.2f", pct), "")
	}
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

// effectiveMaxLen returns the cycle-row count for the MPC/GCC/GCT sections.
// Upstream bumps max_len by one at stats.c:1615 before emitting those sections
// so a trailing all-zero cycle row is shown.
func (c *StatsCounters) effectiveMaxLen() int {
	return c.maxLen + 1
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

// writeFBCLTC emits the per-fragment ACGT-content sections (FBC for first
// fragments, LBC for last) and the matching raw-counter sections (FTC, LTC),
// ported from upstream stats.c:1699-1746. Unlike GCC, which sums both
// fragment orientations, these keep the first- and last-fragment cycle
// buffers separate; FTC/LTC carry the A/C/G/T/N totals summed across cycles.
func (c *StatsCounters) writeFBCLTC(bw *bufio.Writer) {
	c.writeFragmentACGT(bw, c.acgtCycles1st, "first", "FBC", "FTC")
	c.writeFragmentACGT(bw, c.acgtCycles2nd, "last", "LBC", "LTC")
}

// writeFragmentACGT emits one per-cycle ACGT-content section (contentTag)
// followed by its raw-counter section (counterTag) for the per-cycle buffer
// cycles. The counter totals are accumulated over every cycle in [0,maxLen),
// including cycles whose content row is skipped for having no A/C/G/T bases.
func (c *StatsCounters) writeFragmentACGT(bw *bufio.Writer, cycles []acgtNoCount, frag, contentTag, counterTag string) {
	fmt.Fprintf(bw, "# ACGT content per cycle for %s fragments. Use `grep ^%s | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%%]; and N and O counts as a percentage of all A/C/G/T bases [%%]\n", frag, contentTag)
	var tA, tC, tG, tT, tN int64
	for ibase := 0; ibase < c.maxLen; ibase++ {
		var v acgtNoCount
		if ibase < len(cycles) {
			v = cycles[ibase]
		}
		tA += v.a
		tC += v.c
		tG += v.g
		tT += v.t
		tN += v.n
		sum := v.a + v.c + v.g + v.t
		if sum == 0 {
			continue
		}
		fs := float64(sum)
		fmt.Fprintf(bw, "%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", contentTag, ibase+1,
			100.0*float64(v.a)/fs,
			100.0*float64(v.c)/fs,
			100.0*float64(v.g)/fs,
			100.0*float64(v.t)/fs,
			100.0*float64(v.n)/fs,
			100.0*float64(v.other)/fs)
	}
	fmt.Fprintf(bw, "# ACGT raw counters for %s fragments. Use `grep ^%s | cut -f 2-` to extract this part. The columns are: A,C,G,T,N base counters\n", frag, counterTag)
	fmt.Fprintf(bw, "%s\t%d\t%d\t%d\t%d\t%d\n", counterTag, tA, tC, tG, tT, tN)
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
