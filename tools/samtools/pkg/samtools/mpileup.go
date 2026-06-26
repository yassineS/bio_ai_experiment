package samtools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/baq"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// DefaultMpileupMaxDepth is the upstream samtools mpileup -d default. Reads
// past this per-position threshold are silently dropped.
const DefaultMpileupMaxDepth = 8000

// DefaultMpileupMinBaseQ is the upstream samtools mpileup -Q default. Bases
// with quality below this are silently dropped.
const DefaultMpileupMinBaseQ = 13

// MpileupOptions configures the behaviour of Mpileup.
//
// The zero value is a reasonable default ("min-mapq 0, min-baseq 0, max-depth
// 0 meaning no cap") but for upstream parity callers should set MinBaseQ to
// DefaultMpileupMinBaseQ and MaxDepth to DefaultMpileupMaxDepth — the CLI
// does so automatically.
type MpileupOptions struct {
	// Inputs is the list of BAM/SAM paths to pile up. Multi-BAM input emits
	// parallel `depth/bases/quals` columns per BAM.
	Inputs []string
	// FastaRef is the reference FASTA used to fill the per-position
	// reference-base column (the 3rd output column). When empty every
	// position's refbase is `N`, matching upstream.
	FastaRef string
	// Regions are samtools-style "chr[:beg-end]" specifiers (CLI -r).
	Regions []string
	// PositionsFile is a BED (3+ columns) or a 2-column "chrom\tpos" listing
	// that restricts emitted positions (CLI -l).
	PositionsFile string
	// MinMAPQ skips reads with MAPQ below this value (CLI -q).
	MinMAPQ uint8
	// MinBaseQ skips bases with quality below this value (CLI -Q).
	MinBaseQ uint8
	// MaxDepth caps reads per position. 0 means no cap.
	MaxDepth int
	// CountOrphans includes reads whose mate is unmapped (CLI -A). By
	// default unmapped-mate reads in a proper-pair flag context are
	// dropped, matching upstream.
	CountOrphans bool
	// IgnoreOverlaps discards the overlapping half of an overlapping mate
	// pair, matching upstream `-x` and `samtools depth -s`.
	IgnoreOverlaps bool
	// AllPositions emits zero-depth positions inside the covered range of
	// each chromosome touched by at least one read (CLI -a).
	AllPositions bool
	// AllPositionsAllChroms emits every reference position of every chrom
	// in the header (CLI -aa).
	AllPositionsAllChroms bool
	// OutputMapQ appends a per-position MAPQs column (CLI -s).
	OutputMapQ bool
	// OutputBP appends a per-read read-position column (CLI -O).
	OutputBP bool
	// NoBAQ disables BAQ (per-Base Alignment Quality) realignment — CLI
	// -B. BAQ is applied by default whenever a reference (FastaRef) is
	// supplied, matching upstream `bam_plcmd.c:442` (the MPLP_REALN
	// default). With -B the raw input base qualities are used unchanged.
	NoBAQ bool
	// RedoBAQ recomputes BAQ from scratch, discarding any pre-existing BQ
	// aux tag — CLI -E. It maps to the extended-BAQ flag with the redo bit
	// set (upstream passes 7 instead of 3 to sam_prob_realn).
	RedoBAQ bool
	// Output is the path to write to (default stdout). CLI -o.
	Output string
	// Threads is accepted; v1 is single-threaded. CLI -@.
	Threads int
	// BamList, when non-empty, is the path of a file listing one BAM path
	// per line (CLI -b). Lines starting with '#' and blank lines are
	// ignored. The resolved paths are appended to Inputs.
	BamList string
}

// MpileupFile is the file-path entry point for the CLI. It opens each
// BAM/SAM input, the reference FASTA (when given), the positions BED (when
// given), and emits text mpileup records to out. The BAM list (-b) is
// resolved up front and merged with Inputs in the order
// inputs-then-list-file.
//
// When the caller supplies regions (-r) AND every input is a BAM file
// with a sibling .bai index, MpileupFile uses the BAI to seek to the
// relevant chunks; otherwise it falls back to a linear scan. The
// per-position emission logic is identical either way.
func MpileupFile(opts MpileupOptions, out io.Writer) error {
	// Resolve the BAM list (one path per line).
	if opts.BamList != "" {
		extra, err := readBamList(opts.BamList)
		if err != nil {
			return err
		}
		opts.Inputs = append(opts.Inputs, extra...)
	}
	if len(opts.Inputs) == 0 {
		return fmt.Errorf("samtools mpileup: no input files")
	}

	// Open the reference if given.
	var refFA *fasta.RandomAccess
	if opts.FastaRef != "" {
		ra, err := fasta.OpenRandomAccess(opts.FastaRef)
		if err != nil {
			return fmt.Errorf("samtools mpileup: open reference: %w", err)
		}
		defer ra.Close()
		refFA = ra
	}

	// Load positions file if given.
	var posFilter *positionFilter
	if opts.PositionsFile != "" {
		pf, err := loadPositionsFile(opts.PositionsFile)
		if err != nil {
			return fmt.Errorf("samtools mpileup: %w", err)
		}
		posFilter = pf
	}

	// Indexed region fast path: a single coordinate-sorted BAM restricted to
	// -r regions and carrying a .csi/.bai index is read by seeking to the
	// region's index chunks, so only the region's reads are decoded — matching
	// upstream's indexed seek instead of linearly scanning a whole-genome BAM
	// (which churned through every record). openBAMRegionReader returns nil for
	// anything it cannot seek (no index, CRAM/SAM, unsorted), so the linear path
	// below is the unchanged fallback. Output is identical: the streaming pileup
	// still tiles the same region windows.
	if len(opts.Inputs) == 1 && len(opts.Regions) > 0 && !opts.AllPositions && !opts.AllPositionsAllChroms {
		rr, rerr := openBAMRegionReader(opts.Inputs[0], opts.Regions)
		if rerr != nil {
			return fmt.Errorf("samtools mpileup: %w", rerr)
		}
		if rr != nil {
			defer rr.Close()
			return runMpileup([]sam.Reader{rr}, out, opts, refFA, posFilter)
		}
	}

	// Linear path: open every input through alnio (BAM/CRAM/SAM) and scan it.
	readers := make([]sam.Reader, len(opts.Inputs))
	closers := make([]io.Closer, len(opts.Inputs))
	for i, path := range opts.Inputs {
		f, err := openSeekable(path)
		if err != nil {
			closeAll(closers)
			return fmt.Errorf("samtools mpileup: %w", err)
		}
		closers[i] = f
		rd, err := alnio.NewReader(f)
		if err != nil {
			closeAll(closers)
			return fmt.Errorf("samtools mpileup: %s: %w", path, err)
		}
		readers[i] = rd
	}
	defer closeAll(closers)

	return runMpileup(readers, out, opts, refFA, posFilter)
}

// Mpileup is the streaming entry point used by tests. Inputs are already-open
// io.Readers; the output writer receives text mpileup records. The reference
// (refFA) and positions filter (posFilter) may be nil.
func Mpileup(inputs []io.Reader, out io.Writer, opts MpileupOptions, refFA *fasta.RandomAccess, posFilter *positionFilter) error {
	if len(inputs) == 0 {
		return fmt.Errorf("samtools mpileup: no inputs")
	}
	readers := make([]sam.Reader, len(inputs))
	for i, r := range inputs {
		rd, err := alnio.NewReader(r)
		if err != nil {
			return fmt.Errorf("samtools mpileup: input %d: %w", i, err)
		}
		readers[i] = rd
	}
	return runMpileup(readers, out, opts, refFA, posFilter)
}

// runMpileup is the shared pipeline used by both MpileupFile and Mpileup.
// It groups records by reference, then walks each reference one position
// at a time emitting pileup rows.
func runMpileup(readers []sam.Reader, out io.Writer, opts MpileupOptions, refFA *fasta.RandomAccess, posFilter *positionFilter) error {
	hdr := readers[0].Header()
	for i := 1; i < len(readers); i++ {
		// We do NOT require identical @SQ ordering across inputs — upstream
		// samtools mpileup just uses the first input's header to drive the
		// loop over references. Records on a chrom missing from input 0's
		// header are silently skipped, matching upstream.
	}

	if opts.MaxDepth <= 0 {
		opts.MaxDepth = DefaultMpileupMaxDepth
	}

	// Resolve region restrictions up front so both the streaming and buffered
	// walks share them.
	regions, _, err := region.ResolveRegions(opts.Regions, func(name string) int { return hdr.RefIndex(name) })
	if err != nil {
		return err
	}
	// Pre-bucket regions by chrom (using the header's name, not the index)
	// so the per-chrom loop can quickly filter.
	regionByChrom := map[string][][2]int{}
	for _, r := range regions {
		end0 := r.End0
		if end0 > 1<<29 {
			end0 = int(refLengthForName(hdr, r.Region.Chrom))
		}
		regionByChrom[r.Region.Chrom] = append(regionByChrom[r.Region.Chrom], [2]int{r.Beg0, end0})
	}

	// Fast path: a single coordinate-sorted input — WITH OR WITHOUT a region /
	// positions restriction (but not all-positions mode) — is piled up by
	// streaming the records tile by tile, so peak memory is O(tile width x
	// depth) rather than the whole file. Previously any -r/-l restriction fell
	// through to the buffered walk below, which loads every read of the file
	// into memory first (an 11 GB OOM on a whole-genome BAM restricted to one
	// chromosome). The streaming walk feeds emitMpileupWindow the same per-tile
	// record sets and position filter, so the output is byte-identical.
	if len(readers) == 1 && !opts.AllPositions && !opts.AllPositionsAllChroms && headerIsCoordinateSorted(hdr) {
		var rbc map[string][][2]int
		if len(regions) > 0 {
			rbc = mergeRegionByChrom(regionByChrom)
		}
		return runMpileupStreaming(readers[0], out, opts, refFA, rbc, posFilter)
	}

	// Pull records from every input, bucketed by chrom in the order of
	// input 0's header.
	perInputRecs := make([]map[string][]*sam.Record, len(readers))
	for i, rd := range readers {
		recs, rerr := bucketByChrom(rd, opts, hdr)
		if rerr != nil {
			return rerr
		}
		perInputRecs[i] = recs
	}

	// Apply BAQ (per-Base Alignment Quality) realignment. Upstream
	// (bam_plcmd.c:442) realigns every read with sam_prob_realn whenever a
	// reference is loaded and -B was not given; -E adds the redo bit. This
	// lowers the per-base qualities in place, which then feed the quality
	// column and the min-base-quality (-Q) depth filter.
	if refFA != nil && !opts.NoBAQ {
		if err := applyTextMpileupBAQ(perInputRecs, refFA, opts); err != nil {
			return err
		}
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()

	// Decide which chroms to walk.
	var chromsToWalk []string
	switch {
	case len(opts.Regions) > 0:
		// Just the regions' chroms, in the order they were requested,
		// deduplicated to keep emission stable.
		seen := map[string]struct{}{}
		for _, r := range regions {
			if _, ok := seen[r.Region.Chrom]; ok {
				continue
			}
			seen[r.Region.Chrom] = struct{}{}
			chromsToWalk = append(chromsToWalk, r.Region.Chrom)
		}
	case opts.AllPositionsAllChroms:
		for _, ref := range hdr.Refs {
			chromsToWalk = append(chromsToWalk, ref.Name)
		}
	case posFilter != nil:
		// Just the chroms named in the positions file, in header order.
		seen := posFilter.chroms()
		for _, ref := range hdr.Refs {
			if _, ok := seen[ref.Name]; ok {
				chromsToWalk = append(chromsToWalk, ref.Name)
			}
		}
	default:
		// All chroms that any input actually has records on, in header order.
		hit := map[string]struct{}{}
		for _, recs := range perInputRecs {
			for chrom := range recs {
				hit[chrom] = struct{}{}
			}
		}
		for _, ref := range hdr.Refs {
			if _, ok := hit[ref.Name]; ok {
				chromsToWalk = append(chromsToWalk, ref.Name)
			}
		}
	}

	var sc mpileupScratch
	for _, chrom := range chromsToWalk {
		refLen := int(refLengthForName(hdr, chrom))
		if refLen <= 0 {
			continue
		}
		// Walk windows.
		var windows [][2]int
		if iv, ok := regionByChrom[chrom]; ok {
			windows = mergeIntervals(iv)
		} else {
			windows = [][2]int{{0, refLen}}
		}

		// Compute the union of per-chrom records across inputs once.
		var perInputChromRecs [][]*sam.Record
		for _, recs := range perInputRecs {
			perInputChromRecs = append(perInputChromRecs, recs[chrom])
		}

		// Mate-pair overlap removal (on by default, disabled by -x): adjust each
		// pair's qualities once, in coordinate order, before accumulation — the
		// same point and order as bam_plp, so the result is byte-identical.
		if !opts.IgnoreOverlaps {
			for i := range perInputChromRecs {
				applyOverlapRemoval(perInputChromRecs[i])
			}
		}

		// Fetch the contig once for the per-row reference base, so the tiled
		// emit slices it instead of issuing a Fetch per tile.
		var contig []byte
		if refFA != nil {
			b, err := refFA.Fetch(chrom, 0, refFA.Length(chrom))
			if err != nil {
				return fmt.Errorf("samtools mpileup: fetch %s: %w", chrom, err)
			}
			contig = b
		}

		for _, w := range windows {
			beg0 := w[0]
			end0 := w[1]
			if end0 > refLen {
				end0 = refLen
			}
			if beg0 < 0 {
				beg0 = 0
			}
			if beg0 >= end0 {
				continue
			}
			// Tile the window into fixed-width sub-windows so the per-position
			// event matrix emitMpileupWindow builds is bounded by the tile width
			// rather than the whole contig (which for a chromosome-wide scan was
			// the dominant memory cost). Records are coordinate-sorted, so a
			// per-input sliding active set yields exactly the reads overlapping
			// each tile; reads spanning a tile boundary are carried into the next
			// tile, where accumulateRecordEvents places only their in-tile
			// columns — making the tiled output byte-identical to a single
			// whole-window pass. BAQ is already applied to every record above,
			// so tiling does not perturb base qualities.
			nextIdx := make([]int, len(perInputChromRecs))
			active := make([][]*sam.Record, len(perInputChromRecs))
			for tBeg := beg0; tBeg < end0; tBeg += mpileupTileWidth {
				tEnd := tBeg + mpileupTileWidth
				if tEnd > end0 {
					tEnd = end0
				}
				for i, recs := range perInputChromRecs {
					active[i] = pruneEndedBefore(active[i], tBeg)
					for nextIdx[i] < len(recs) && int(recs[nextIdx[i]].Pos)-1 < tEnd {
						active[i] = append(active[i], recs[nextIdx[i]])
						nextIdx[i]++
					}
				}
				if err := emitMpileupWindow(bw, chrom, tBeg, tEnd, refLen,
					active, contig, refFA, posFilter, opts, &sc); err != nil {
					return err
				}
			}
		}
	}
	return nil
}

// mpileupTileWidth is the reference span of one pileup tile. It bounds the
// per-position event matrix (and the live read set) emitMpileupWindow holds at
// once, so peak memory is O(tile width x depth) rather than O(contig length).
const mpileupTileWidth = 16 * 1024

// runMpileupStreaming piles up a single coordinate-sorted input without ever
// buffering the whole file: records are pulled from a peekable, pre-filtered
// source tile by tile, BAQ is applied to each read as it enters the active set,
// and reads are dropped once the pileup column passes their end. Peak memory is
// therefore the reads overlapping one tile, matching upstream's streaming
// pileup. The per-tile record sets handed to emitMpileupWindow are identical to
// those the buffered path produces for a coordinate-sorted input, so output is
// byte-for-byte the same.
// mergeRegionByChrom returns the per-chrom windows sorted ascending and with
// overlapping/adjacent intervals merged, so the forward-only streaming walk can
// emit them in a single left-to-right pass (matching upstream's coordinate
// order). The input map is left unmodified.
func mergeRegionByChrom(in map[string][][2]int) map[string][][2]int {
	out := make(map[string][][2]int, len(in))
	for chrom, ws := range in {
		cp := append([][2]int(nil), ws...)
		sort.Slice(cp, func(i, j int) bool { return cp[i][0] < cp[j][0] })
		merged := cp[:0]
		for _, w := range cp {
			if n := len(merged); n > 0 && w[0] <= merged[n-1][1] {
				if w[1] > merged[n-1][1] {
					merged[n-1][1] = w[1]
				}
				continue
			}
			merged = append(merged, w)
		}
		out[chrom] = merged
	}
	return out
}

func runMpileupStreaming(rd sam.Reader, out io.Writer, opts MpileupOptions, refFA *fasta.RandomAccess, regionByChrom map[string][][2]int, posFilter *positionFilter) error {
	hdr := rd.Header()
	bw := bufio.NewWriter(out)
	defer bw.Flush()

	doBAQ := refFA != nil && !opts.NoBAQ
	doOverlap := !opts.IgnoreOverlaps // mate-pair overlap removal is on by default
	baqFlag := textMpileupBAQFlag(opts)
	src := newMpileupSource(rd, opts, hdr)
	var sc mpileupScratch

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
			// Chromosome absent from the header (or zero length): its records
			// are never emitted, exactly as in the buffered walk. Drop them.
			dropChrom(chrom)
			continue
		}

		// Determine the coordinate windows to emit on this chrom: the requested
		// regions clamped to the contig (sorted, non-overlapping) when region-
		// restricted, else the whole contig. A chrom outside every region is
		// streamed past without emitting, so peak memory stays O(tile x depth)
		// even with -r — the buffered walk would have loaded the whole file.
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

		// Fetch the whole contig once: it serves both BAQ and the per-row
		// reference base (emitMpileupWindow slices it instead of re-fetching).
		var contig []byte
		if refFA != nil {
			seq, err := refFA.Fetch(chrom, 0, refFA.Length(chrom))
			if err != nil {
				return fmt.Errorf("samtools mpileup: fetch %s: %w", chrom, err)
			}
			contig = seq
		}

		var active []*sam.Record
		// overlaps holds the earlier-seen mate of each pair until the later one
		// arrives (reset per chrom), so overlap removal is applied as reads enter
		// the pileup, in coordinate order — the same order htslib's bam_plp uses.
		overlaps := make(map[string]*sam.Record)
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
			active = pruneEndedBefore(active, wBeg)
			for tBeg := wBeg; tBeg < wEnd; tBeg += mpileupTileWidth {
				tEnd := tBeg + mpileupTileWidth
				if tEnd > wEnd {
					tEnd = wEnd
				}
				active = pruneEndedBefore(active, tBeg)
				for {
					p := src.peek()
					if p == nil || p.RName != chrom || int(p.Pos)-1 >= tEnd {
						break
					}
					// Consume the record regardless, but only retain (and BAQ-
					// realign) it when it actually overlaps this tile. Reads that
					// ended before the tile — e.g. everything between the contig
					// start and a deep -r region — are popped without being held,
					// which is what keeps `active` small when streaming to a region
					// that starts far into the contig.
					if int(p.EndPosition()) > tBeg {
						if doBAQ {
							if r := baq.SamProbRealn(p, contig, baqFlag); r < -3 {
								return fmt.Errorf("samtools mpileup: BAQ alignment failed for read %q", p.QName)
							}
						}
						// Overlap removal runs after BAQ (matching bam_plp), tweaking
						// this read's qualities against an already-seen mate.
						if doOverlap {
							overlapPush(p, overlaps)
						}
						active = append(active, p)
					}
					src.pop()
				}
				if err := emitMpileupWindow(bw, chrom, tBeg, tEnd, refLen,
					[][]*sam.Record{active}, contig, refFA, posFilter, opts, &sc); err != nil {
					return err
				}
			}
		}
		// Discard any remaining records of this chrom (before the first window,
		// between windows, or beyond the last) so the loop advances to the next
		// chromosome rather than spinning.
		dropChrom(chrom)
		if src.err != nil {
			return src.err
		}
	}
	return src.err
}

// mpileupSource is a peekable stream of mpileup-eligible records from one input.
// advance applies the same record-level filters as the buffered bucketByChrom,
// so the pileup sees an identical record set.
type mpileupSource struct {
	rd   sam.Reader
	opts MpileupOptions
	hdr0 *sam.Header
	head *sam.Record // next record to yield; nil at EOF or after an error
	err  error
}

func newMpileupSource(rd sam.Reader, opts MpileupOptions, hdr0 *sam.Header) *mpileupSource {
	s := &mpileupSource{rd: rd, opts: opts, hdr0: hdr0}
	s.advance()
	return s
}

func (s *mpileupSource) advance() {
	for {
		rec, err := s.rd.Read()
		if err == io.EOF {
			s.head = nil
			return
		}
		if err != nil {
			s.err = err
			s.head = nil
			return
		}
		if !keepMpileupRecord(rec, s.opts, s.hdr0) {
			continue
		}
		s.head = rec
		return
	}
}

func (s *mpileupSource) peek() *sam.Record { return s.head }

func (s *mpileupSource) pop() *sam.Record {
	r := s.head
	s.advance()
	return r
}

// headerIsCoordinateSorted reports whether the @HD line declares SO:coordinate,
// the precondition for the streaming pileup (records arrive in reference order).
func headerIsCoordinateSorted(hdr *sam.Header) bool {
	for _, f := range hdr.HDFields {
		if f.Tag == "SO" {
			return f.Value == "coordinate"
		}
	}
	return false
}

// pruneEndedBefore drops records whose alignment ends at or before tBeg (0-based
// half-open end), i.e. reads that no longer overlap the upcoming tile. It
// compacts in place, preserving the coordinate order of the survivors.
func pruneEndedBefore(active []*sam.Record, tBeg int) []*sam.Record {
	kept := active[:0]
	for _, r := range active {
		if int(r.EndPosition()) > tBeg {
			kept = append(kept, r)
		}
	}
	return kept
}

// textMpileupBAQFlag derives the realn flag passed to baq.SamProbRealn.
// Upstream samtools text mpileup always realigns in apply+extend mode
// (bam_plcmd.c passes 3 = BAQ_APPLY|BAQ_EXTEND), and -E adds the redo bit
// (passing 7 = BAQ_APPLY|BAQ_EXTEND|BAQ_REDO).
func textMpileupBAQFlag(opts MpileupOptions) int {
	flag := baq.FlagApply | baq.FlagExtend
	if opts.RedoBAQ {
		flag |= baq.FlagRedo
	}
	return flag
}

// applyTextMpileupBAQ realigns every bucketed read's base qualities in
// place using BAQ, contig by contig. Each contig's reference sequence is
// fetched once (records are bucketed per chrom) and reused across inputs,
// mirroring samtools calmd's single-slot reference cache. A read with no
// BAQ work (unmapped, no M/=/X, etc.) is left untouched: SamProbRealn
// returns a benign -1/-3 in those cases, and only a hard alignment failure
// (-4 and below) is surfaced as an error.
func applyTextMpileupBAQ(perInputRecs []map[string][]*sam.Record, refFA *fasta.RandomAccess, opts MpileupOptions) error {
	flag := textMpileupBAQFlag(opts)
	// Collect the set of contigs that carry records, so we fetch each
	// reference span exactly once.
	contigSeq := map[string][]byte{}
	for _, recs := range perInputRecs {
		for chrom := range recs {
			if _, ok := contigSeq[chrom]; ok {
				continue
			}
			seq, err := refFA.Fetch(chrom, 0, refFA.Length(chrom))
			if err != nil {
				return fmt.Errorf("samtools mpileup: BAQ fetch %s: %w", chrom, err)
			}
			contigSeq[chrom] = seq
		}
	}
	for _, recs := range perInputRecs {
		for chrom, list := range recs {
			ref := contigSeq[chrom]
			for _, rec := range list {
				if r := baq.SamProbRealn(rec, ref, flag); r < -3 {
					return fmt.Errorf("samtools mpileup: BAQ alignment failed for read %q", rec.QName)
				}
			}
		}
	}
	return nil
}

// bucketByChrom drains rd into a per-chrom record map, applying the
// chromosome-agnostic read filters (flag bits, MAPQ, etc).
func bucketByChrom(rd sam.Reader, opts MpileupOptions, hdr0 *sam.Header) (map[string][]*sam.Record, error) {
	out := map[string][]*sam.Record{}
	for {
		rec, err := rd.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if !keepMpileupRecord(rec, opts, hdr0) {
			continue
		}
		out[rec.RName] = append(out[rec.RName], rec)
	}
	// Pre-sort by Pos to make the sliding window deterministic; SAM/BAM
	// files are usually coordinate-sorted already, but tests sometimes pass
	// hand-written SAM where the records are in a different order.
	for k := range out {
		sort.SliceStable(out[k], func(i, j int) bool {
			return out[k][i].Pos < out[k][j].Pos
		})
	}
	return out, nil
}

// keepMpileupRecord applies the record-level filters (flag bits, MAPQ).
func keepMpileupRecord(rec *sam.Record, opts MpileupOptions, hdr0 *sam.Header) bool {
	if rec.Pos <= 0 || rec.RName == "" {
		return false
	}
	// Match upstream defaults: drop unmapped, secondary, supplementary,
	// QCfail, duplicate.
	if rec.Flag&(sam.FlagUnmapped|sam.FlagSecondary|sam.FlagQCFail|sam.FlagDuplicate|sam.FlagSupplementary) != 0 {
		return false
	}
	// `-A` accepts reads whose mate is unmapped; without it, anomalous
	// pairs (paired-but-mate-unmapped, or paired-but-not-proper) are
	// dropped — same heuristic as upstream `bam_plcmd.c::mplp_get_read`.
	if !opts.CountOrphans && rec.Flag&sam.FlagPaired != 0 {
		if rec.Flag&sam.FlagMateUnmapped != 0 {
			return false
		}
		// Drop paired-but-not-proper (anomalous) reads, matching the
		// upstream `-A` default.
		if rec.Flag&sam.FlagProperPair == 0 {
			return false
		}
	}
	if opts.MinMAPQ > 0 && rec.MapQ < opts.MinMAPQ {
		return false
	}
	if hdr0 != nil && hdr0.RefIndex(rec.RName) < 0 {
		// Chrom missing from the leading input header — silently drop.
		return false
	}
	return true
}

// refLengthForName looks up a contig length on hdr; we kept the existing
// refLength helper but redirect through it for clarity.
func refLengthForName(hdr *sam.Header, name string) int32 {
	for _, r := range hdr.Refs {
		if r.Name == name {
			return r.Length
		}
	}
	return 0
}

// closeAll closes every non-nil io.Closer in slc.
func closeAll(slc []io.Closer) {
	for _, c := range slc {
		if c != nil {
			_ = c.Close()
		}
	}
}

// readBamList reads a one-path-per-line BAM list file.
func readBamList(path string) ([]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open bam list: %w", err)
	}
	defer f.Close()
	var out []string
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		ln := strings.TrimSpace(sc.Text())
		if ln == "" || strings.HasPrefix(ln, "#") {
			continue
		}
		out = append(out, ln)
	}
	return out, sc.Err()
}

// positionFilter restricts emission to a set of (chrom, 1-based position)
// pairs. It accepts BED records (chrom\tbegin0\tend0) and 2-column
// chrom\tpos lines transparently.
type positionFilter struct {
	// byChrom[chrom] is the sorted, merged list of 1-based [beg,end]
	// inclusive intervals that are allowed.
	byChrom map[string][][2]int
}

// chroms returns the set of chroms the filter touches.
func (p *positionFilter) chroms() map[string]struct{} {
	out := make(map[string]struct{}, len(p.byChrom))
	for k := range p.byChrom {
		out[k] = struct{}{}
	}
	return out
}

// contains reports whether pos1 (1-based) on chrom is inside the filter.
func (p *positionFilter) contains(chrom string, pos1 int) bool {
	if p == nil {
		return true
	}
	iv, ok := p.byChrom[chrom]
	if !ok {
		return false
	}
	// Binary search would be nicer; linear is fine for typical positions
	// files which are dozens-to-thousands of intervals.
	for _, r := range iv {
		if pos1 >= r[0] && pos1 <= r[1] {
			return true
		}
	}
	return false
}

// loadPositionsFile reads a BED or 2-column position list and returns a
// positionFilter. Each line is split on tab and dispatched by column count:
//
//   - 2 columns: "chrom\tpos1" — a single 1-based position
//   - >=3 columns: "chrom\tbegin0\tend0[\t...]" — a half-open 0-based BED
//     interval
//
// Comment lines (starting with '#' or 'track' / 'browser') are skipped.
func loadPositionsFile(path string) (*positionFilter, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("open positions file: %w", err)
	}
	defer f.Close()
	pf := &positionFilter{byChrom: map[string][][2]int{}}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 1<<20)
	for sc.Scan() {
		ln := strings.TrimRight(sc.Text(), "\r\n")
		if ln == "" || strings.HasPrefix(ln, "#") || strings.HasPrefix(ln, "track") || strings.HasPrefix(ln, "browser") {
			continue
		}
		fields := strings.Split(ln, "\t")
		switch {
		case len(fields) >= 3:
			chrom := fields[0]
			b, perr := strconv.Atoi(fields[1])
			if perr != nil {
				return nil, fmt.Errorf("positions file: parse begin %q: %w", fields[1], perr)
			}
			e, perr := strconv.Atoi(fields[2])
			if perr != nil {
				return nil, fmt.Errorf("positions file: parse end %q: %w", fields[2], perr)
			}
			// BED is half-open 0-based; mpileup positions are 1-based
			// inclusive. Translate [b, e) to [b+1, e].
			pf.byChrom[chrom] = append(pf.byChrom[chrom], [2]int{b + 1, e})
		case len(fields) == 2:
			chrom := fields[0]
			p, perr := strconv.Atoi(fields[1])
			if perr != nil {
				return nil, fmt.Errorf("positions file: parse pos %q: %w", fields[1], perr)
			}
			pf.byChrom[chrom] = append(pf.byChrom[chrom], [2]int{p, p})
		default:
			return nil, fmt.Errorf("positions file: bad line %q (need 2 or >=3 columns)", ln)
		}
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	// Sort + merge intervals per chrom for stable membership checks.
	for k, iv := range pf.byChrom {
		sort.Slice(iv, func(i, j int) bool { return iv[i][0] < iv[j][0] })
		merged := iv[:0]
		cur := iv[0]
		for i := 1; i < len(iv); i++ {
			if iv[i][0] <= cur[1]+1 {
				if iv[i][1] > cur[1] {
					cur[1] = iv[i][1]
				}
			} else {
				merged = append(merged, cur)
				cur = iv[i]
			}
		}
		merged = append(merged, cur)
		pf.byChrom[k] = merged
	}
	return pf, nil
}
