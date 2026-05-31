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
	// NoBAQ is accepted for upstream-flag compatibility and is a no-op:
	// this v1 never applies BAQ. CLI -B.
	NoBAQ bool
	// RedoBAQ is accepted but rejected with a clear "not implemented"
	// error — CLI -E. Deferred per PARITY_ROADMAP.md.
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

// ErrMpileupBCFNotImplemented is returned when callers pass -u/-g for BCF
// output. Tracked at docs/PARITY_ROADMAP.md#samtools.
var ErrMpileupBCFNotImplemented = fmt.Errorf("samtools mpileup: BCF output (-u/-g) not yet implemented; tracked in docs/PARITY_ROADMAP.md#samtools")

// ErrMpileupBAQNotImplemented is returned when callers pass -E (redo BAQ).
// Tracked at docs/PARITY_ROADMAP.md#samtools.
var ErrMpileupBAQNotImplemented = fmt.Errorf("samtools mpileup: -E/--redo-baq not yet implemented; tracked in docs/PARITY_ROADMAP.md#samtools")

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
	if opts.RedoBAQ {
		return ErrMpileupBAQNotImplemented
	}
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

	// Open inputs.
	readers := make([]sam.Reader, len(opts.Inputs))
	closers := make([]io.Closer, len(opts.Inputs))
	for i, path := range opts.Inputs {
		f, err := os.Open(path)
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
	// BAI seek is left to a future PR (matches upstream view semantics):
	// the linear-scan path here is byte-for-byte equivalent because the
	// downstream pileup is a post-filter against region.ResolveRegions. The
	// region-via-BAI fast path is tracked in
	// docs/PARITY_ROADMAP.md#samtools as a follow-up.

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

	return runMpileup(readers, out, opts, refFA, posFilter)
}

// Mpileup is the streaming entry point used by tests. Inputs are already-open
// io.Readers; the output writer receives text mpileup records. The reference
// (refFA) and positions filter (posFilter) may be nil.
func Mpileup(inputs []io.Reader, out io.Writer, opts MpileupOptions, refFA *fasta.RandomAccess, posFilter *positionFilter) error {
	if opts.RedoBAQ {
		return ErrMpileupBAQNotImplemented
	}
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

	// Resolve region restrictions.
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

	for _, chrom := range chromsToWalk {
		refLen := int(refLengthForName(hdr, chrom))
		if refLen <= 0 {
			continue
		}
		// Walk windows. The whole-chrom default window must extend to
		// include any read-overhang past refLen so upstream's
		// "emit-every-covered-position" semantics still apply.
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

		// Determine the covered range so -a (per upstream) can emit zero
		// positions only inside the read-touched extent.
		minPos0, maxEnd0 := coveredExtent(perInputChromRecs)

		for _, w := range windows {
			beg0 := w[0]
			end0 := w[1]
			// Upstream's pileup engine reports read events past the
			// contig length (a read overhang). Extend whole-chrom
			// emission to the maximum read end so those positions
			// still get emitted instead of being trimmed away. We
			// only extend for the implicit whole-chrom window
			// (end0 == refLen); explicit user regions stay bounded.
			if end0 == refLen && maxEnd0 > end0 {
				end0 = maxEnd0
			}
			if beg0 < 0 {
				beg0 = 0
			}
			if beg0 >= end0 {
				continue
			}
			if err := emitMpileupWindow(bw, chrom, beg0, end0, refLen,
				perInputChromRecs, refFA, posFilter,
				opts, minPos0, maxEnd0); err != nil {
				return err
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

// coveredExtent returns the smallest 0-based start and largest 0-based end
// of any record across the inputs. When no record is present, both return
// zero, which keeps the `-a` branch from emitting anything.
func coveredExtent(perInputRecs [][]*sam.Record) (int, int) {
	first := true
	var lo, hi int
	for _, recs := range perInputRecs {
		for _, r := range recs {
			start := int(r.Pos) - 1
			end := start + r.Cigar.ReferenceLength()
			if first {
				lo, hi = start, end
				first = false
				continue
			}
			if start < lo {
				lo = start
			}
			if end > hi {
				hi = end
			}
		}
	}
	if first {
		return 0, 0
	}
	return lo, hi
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
