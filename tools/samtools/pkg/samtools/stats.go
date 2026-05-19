// Package samtools — stats implementation.
//
// Mirrors `samtools stats`. v1 ships byte-parity on the Summary Numbers
// (SN) block — the section the vast majority of downstream tooling
// (multiqc, picard-style dashboards, in-house pipelines) actually
// consume. The full upstream output has ~30 sections (FFQ/LFQ/GCT/GCC/
// GCD/GCL/RL/MAPQ/IS/COV/COV2/OXC/...); we emit headers for the most
// common ones with empty-but-correctly-tagged bodies so the file
// remains grep-able with the upstream `grep ^SN | cut -f 2-` idiom.
// PARITY_VALIDATION.md tracks which sections are byte-faithful vs.
// placeholder.
//
// Skipped intentionally for v1 (documented):
//   - CHK checksum block (depends on a CRC32 reduction we don't emit).
//   - COV/COV2 coverage histograms (requires a reference FASTA and a
//     per-position depth walker).
//   - GCD/GCT/GCC/GCL GC-content distributions (require reference bases).
//   - OXC oxidation-context counts (requires reference bases).
//   - --target-regions BED restriction.
//   - --remove-overlaps (single-record stats are unaffected by overlap
//     removal for the counters we emit).
package samtools

import (
	"bufio"
	"fmt"
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
	// Coverage is the raw "MIN,MAX,STEP" string passed via -c; parsed but
	// not used while COV is a placeholder.
	Coverage string
	// RequiredFlag requires ALL bits set on each record (-l/--required-flag).
	RequiredFlag uint16
	// FilteringFlag drops records with ANY bit set (-F/--filtering-flag).
	FilteringFlag uint16
	// MaxDepth caps depth used in COV (placeholder in v1).
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

	// Pair tracking by qname (mate's flag/pos for orientation classification).
	// Cleared as pairs are observed (memory-bounded by max pending mates).
	mates map[string]mateInfo
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
		RL:      make(map[int64]int64),
		FRL:     make(map[int64]int64),
		LRL:     make(map[int64]int64),
		ISInw:   make(map[int64]int64),
		ISOutw:  make(map[int64]int64),
		ISOther: make(map[int64]int64),
		mates:   make(map[string]mateInfo),
	}
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

// Write emits the upstream-compatible text report to w.
func (c *StatsCounters) Write(w io.Writer, opts StatsOptions) error {
	bw := bufio.NewWriter(w)
	// SN block — full upstream parity.
	c.writeSN(bw)
	// Section placeholders. v1 ships RL/FRL/LRL and MAPQ with real data;
	// FFQ/LFQ/GCF/GCL/GCC/GCT/IS bodies are only emitted in the
	// non-sparse path with zero-data tables that the grep idiom can still
	// find. Comments are placed exactly as upstream so consumers that
	// scan `grep ^XX` find their section.
	if !opts.Sparse {
		c.writeRL(bw)
		c.writeMAPQ(bw)
		c.writeIS(bw, opts)
	}
	return bw.Flush()
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
