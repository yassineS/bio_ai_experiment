package samtools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

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
	// MinMAPQ skips reads with MAPQ < this value (CLI `-q`).
	MinMAPQ uint8
	// MinBaseQ skips bases with quality < this value (CLI `-Q`).
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
func Depth(inputs []io.Reader, out io.Writer, opts DepthOptions) error {
	if len(inputs) == 0 {
		return fmt.Errorf("samtools depth: no inputs")
	}
	readers := make([]sam.Reader, len(inputs))
	for i, r := range inputs {
		rd, err := alnio.NewReaderThreaded(r, "", opts.Threads)
		if err != nil {
			return fmt.Errorf("samtools depth: input %d: %w", i, err)
		}
		if rc, ok := rd.(io.Closer); ok {
			defer rc.Close()
		}
		readers[i] = rd
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

	// We compute depth one reference at a time. Pre-buffer records per
	// reference so we can do a streaming sliding-window across positions.
	perRefRecords := make([][]map[string][]*sam.Record, len(readers))
	for i, rd := range readers {
		recs, rerr := groupRecordsByRef(rd, opts)
		if rerr != nil {
			return fmt.Errorf("samtools depth: input %d: %w", i, rerr)
		}
		perRefRecords[i] = []map[string][]*sam.Record{recs}
	}

	// Iterate references in header order (subset by region.names).
	for _, ref := range region.names {
		// Pull intervals for this ref.
		intervals, ok := region.byRef[ref]
		if !ok {
			continue
		}
		// Build per-input counters as a slice of position→depth maps.
		// To keep memory bounded we sort intervals and walk them in turn.
		if intervals == nil {
			// "All positions of this reference" — use the @SQ length.
			refLen := refLength(hdr, ref)
			if refLen <= 0 {
				continue
			}
			intervals = [][2]int{{0, refLen}}
		}
		// Sort + dedupe overlaps.
		intervals = mergeIntervals(intervals)

		// Gather per-input records for this reference.
		perInputRecs := make([][]*sam.Record, len(readers))
		for i := range readers {
			perInputRecs[i] = perRefRecords[i][0][ref]
		}

		for _, iv := range intervals {
			if err := emitDepthInterval(bw, ref, iv[0], iv[1], perInputRecs, opts); err != nil {
				return err
			}
		}
	}
	return nil
}

// emitDepthInterval emits per-position depth lines for one [beg0, end0)
// half-open interval on the named reference, across the parallel per-input
// record slices.
func emitDepthInterval(bw *bufio.Writer, ref string, beg0, end0 int, perInputRecs [][]*sam.Record, opts DepthOptions) error {
	if end0 <= beg0 {
		return nil
	}
	n := len(perInputRecs)
	// Pre-compute per-input cumulative depth across the interval using a
	// difference array, then walk it. This is O(reads + width) per input.
	width := end0 - beg0
	diffs := make([][]int32, n)
	for i := range perInputRecs {
		diffs[i] = make([]int32, width+1)
		for _, rec := range perInputRecs[i] {
			addReadDepth(rec, beg0, end0, opts, diffs[i])
		}
	}

	// Stream the depth values position by position.
	cur := make([]int32, n)
	for pos0 := beg0; pos0 < end0; pos0++ {
		// Apply diffs[i][pos0 - beg0] to cur[i].
		idx := pos0 - beg0
		any := false
		for i := range cur {
			cur[i] += diffs[i][idx]
			if cur[i] > 0 {
				any = true
			}
		}
		if !any && !opts.AllPositions && !opts.AllTransPositions {
			continue
		}
		// Emit pos+1 (SAM is 1-based).
		if _, err := bw.WriteString(ref); err != nil {
			return err
		}
		if err := bw.WriteByte('\t'); err != nil {
			return err
		}
		bw.WriteString(strconv.Itoa(pos0 + 1))
		for i := 0; i < n; i++ {
			d := cur[i]
			if opts.MaxDepth > 0 && d > int32(opts.MaxDepth) {
				d = int32(opts.MaxDepth)
			}
			bw.WriteByte('\t')
			bw.WriteString(strconv.FormatInt(int64(d), 10))
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
	}
	return nil
}

// addReadDepth walks a record's CIGAR and increments diff[start-beg0] /
// decrements diff[end-beg0] for each contiguous reference-consuming run
// (M/=/X). Soft/hard clips and insertions consume query bases but not
// reference; deletions and reference skips consume reference but do not add
// depth.
//
// Bases below opts.MinBaseQ are skipped on a per-base basis when SEQ/QUAL
// information is available.
func addReadDepth(rec *sam.Record, beg0, end0 int, opts DepthOptions, diff []int32) {
	if rec.Pos <= 0 {
		return
	}
	refPos := int(rec.Pos) - 1
	queryPos := 0
	for _, op := range rec.Cigar {
		l := int(op.Length())
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			runBeg := refPos
			runEnd := refPos + l
			// Intersect with the requested interval before recording.
			clipBeg := runBeg
			if clipBeg < beg0 {
				clipBeg = beg0
			}
			clipEnd := runEnd
			if clipEnd > end0 {
				clipEnd = end0
			}
			if clipEnd > clipBeg {
				if opts.MinBaseQ > 0 && len(rec.Qual) > 0 {
					// Per-base baseQ filter: bump diff in single-position
					// steps so each base can be accepted or rejected.
					for p := clipBeg; p < clipEnd; p++ {
						qIdx := queryPos + (p - runBeg)
						if qIdx >= 0 && qIdx < len(rec.Qual) && rec.Qual[qIdx] < opts.MinBaseQ {
							continue
						}
						diff[p-beg0]++
						diff[p+1-beg0]--
					}
				} else {
					diff[clipBeg-beg0]++
					diff[clipEnd-beg0]--
				}
			}
			refPos += l
			queryPos += l
		case sam.CigarInsertion, sam.CigarSoftClip:
			queryPos += l
		case sam.CigarDeletion, sam.CigarSkipped:
			refPos += l
		case sam.CigarHardClip, sam.CigarPadding:
			// Neither consumes ref nor query (in our accounting).
		}
		if refPos >= end0 {
			return
		}
	}
}

// groupRecordsByRef drains a single sam.Reader into per-reference record
// slices, applying flag / mapq / readlen filters as it goes so we never hold
// records we are going to drop.
func groupRecordsByRef(rd sam.Reader, opts DepthOptions) (map[string][]*sam.Record, error) {
	out := map[string][]*sam.Record{}
	for {
		rec, err := rd.Read()
		if err == io.EOF {
			return out, nil
		}
		if err != nil {
			return nil, err
		}
		if !keepDepthRecord(rec, opts) {
			continue
		}
		out[rec.RName] = append(out[rec.RName], rec)
	}
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
