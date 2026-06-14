package mosdepth

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	bgzf "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// DefaultExcludeFlag is mosdepth's default `-F` value: drop reads that are
// unmapped (0x4), secondary (0x100), QC-fail (0x200), or duplicate (0x400).
// OR'd together this is 0x704 = 1796 — the exact value upstream mosdepth
// prints when invoked with --help. Note it does NOT include supplementary
// (0x800): upstream's default keeps supplementary alignments, so excluding
// them here would diverge from upstream per-base/region depth.
const DefaultExcludeFlag uint16 = sam.FlagUnmapped | sam.FlagSecondary |
	sam.FlagQCFail | sam.FlagDuplicate

// Options configures a single mosdepth run.
type Options struct {
	// Prefix is the output-file prefix. Files written are:
	//   <prefix>.mosdepth.global.dist.txt
	//   <prefix>.mosdepth.summary.txt
	//   <prefix>.per-base.bed.gz                (and .csi)         — unless NoPerBase
	//   <prefix>.regions.bed.gz                 (and .csi)         — only when ByBED/ByWindow is set
	//   <prefix>.thresholds.bed.gz              (and .csi)         — only when Thresholds is non-empty
	Prefix string

	// ByBED is a path to a BED file of regions. Mutually exclusive with
	// ByWindow.
	ByBED string
	// ByWindow is a fixed integer window size in bases. When >0, mosdepth
	// emits a regions.bed.gz tiled across each reference at this stride.
	// Mutually exclusive with ByBED.
	ByWindow int

	// MinMAPQ excludes reads with MAPQ < MinMAPQ (CLI `-Q`).
	MinMAPQ uint8
	// ExcludeFlag drops a read with ANY of these bits set (CLI `-F`).
	ExcludeFlag uint16
	// IncludeFlag keeps only reads with ALL these bits set (CLI `-i`).
	IncludeFlag uint16

	// FastMode skips CIGAR walking and treats each read as covering
	// POS..POS+ReferenceLength. ~3x faster, slightly inaccurate around
	// indels.
	FastMode bool

	// NoPerBase suppresses the per-base.bed.gz output entirely.
	NoPerBase bool

	// Thresholds is the parsed list of integer depths (sorted ascending);
	// when non-empty mosdepth emits a thresholds.bed.gz file with the
	// count of bases at or above each threshold per region.
	Thresholds []int

	// Chrom restricts processing to the named chromosome. Other chromosomes
	// are skipped entirely. Empty means "all chromosomes".
	Chrom string

	// D4Output, when true, writes the per-base depth track to
	// <prefix>.per-base.d4 in the dense D4 binary format instead of the
	// bgzipped per-base BED (matching upstream mosdepth's --d4 behavior).
	// It has no effect when NoPerBase is set or when --by selects regions.
	D4Output bool

	// ReadGroups, when non-empty, keeps reads whose RG aux tag is in the
	// set. The special prefix "OPS:" filters on the OPS aux tag instead
	// (matching upstream mosdepth's --read-groups OPS:X,Y).
	ReadGroups []string

	// MinFragLen / MaxFragLen filter by absolute TLEN. <=0 means no
	// minimum / maximum, respectively.
	MinFragLen int
	MaxFragLen int

	// FragmentMode, when true, counts coverage across the whole template
	// (fragment) of each properly-paired read pair rather than across the
	// aligned reads only — matching upstream mosdepth's --fragment-mode.
	// Only read1 of a proper, non-supplementary pair contributes; it covers
	// [min(read,mate) start, +|TLEN|). Mutually exclusive with FastMode in
	// upstream's hot loop (fragment-mode takes precedence when both are set
	// upstream; this port rejects the combination at the CLI).
	FragmentMode bool

	// Quantize, when non-empty, holds the parsed quantize bin boundaries
	// (see ParseQuantize). When set, mosdepth writes a
	// <prefix>.quantized.bed.gz binning each base's depth into the
	// corresponding labelled segment, matching upstream's --quantize output.
	Quantize []int

	// Threads sets the number of BAM/BGZF decompression worker threads. The
	// decoded output is identical regardless of this value; it only affects
	// throughput. 0 or 1 means single-threaded. See OpenAndRun.
	Threads int

	// Fasta names a FASTA reference (with a sibling .fai) used to decode
	// reference-backed CRAM input (CLI -f/--fasta / --reference). It is
	// ignored for BAM and SAM input, which carry their sequence inline. When
	// empty, a CRAM input still decodes — reference-derived bases fall back to
	// 'N' — and the REF_CACHE environment variable is honoured as an
	// additional CRAM reference source regardless of this value. Depth only
	// depends on alignment coordinates, so CRAM depth output is identical to
	// the equivalent BAM whether or not a reference is supplied.
	Fasta string

	// UseMedian, when true, makes the per-region (--by) output report the
	// median per-base depth instead of the mean, matching upstream mosdepth's
	// -m/--use-median. It affects ONLY the regions.bed.gz depth column; the
	// summary, distribution, thresholds, quantized, and per-base outputs are
	// unchanged — upstream routes only the regions mean through its
	// median-capable imean(), leaving the sum-based summary mean intact.
	UseMedian bool
}

// ErrByConflict is returned when both ByBED and ByWindow are set.
var ErrByConflict = errors.New("mosdepth: -b/--by cannot specify both a BED file and an integer window")

// ErrBadFragLenBounds is returned when MaxFragLen is a positive cap lower than
// MinFragLen, which can never keep any read. Upstream mosdepth rejects this
// combination with `[mosdepth] error --max-frag-len was lower than
// --min-frag-len.` and exit code 2; this error carries the same message so the
// CLI can print it verbatim.
var ErrBadFragLenBounds = errors.New("[mosdepth] error --max-frag-len was lower than --min-frag-len.")

// ChromNotFoundError reports that a --chrom restriction named a reference that
// is absent from the input header. Upstream mosdepth prints `[mosdepth]
// chromosome <name> not found` and exits 1; the CLI formats this error the same
// way. The Chrom field holds the offending name.
type ChromNotFoundError struct{ Chrom string }

// Error formats the message exactly as upstream mosdepth's check_chrom does.
func (e *ChromNotFoundError) Error() string {
	return "[mosdepth] chromosome " + e.Chrom + " not found"
}

// Run executes a full mosdepth pipeline against the SAM/BAM bytes streaming
// in from in. The header is read first, then records are streamed and depth
// is accumulated per reference. When the cursor moves to a new reference
// (or input EOF) the previous reference's outputs are emitted.
//
// The function intentionally takes an io.Reader rather than a file path so
// that tests can drive it from an in-memory BAM buffer; the CLI front-end
// opens the file and passes its handle. CRAM input, which needs its own
// reference-aware decoder, is routed through OpenAndRun / runWithReader.
func Run(in io.Reader, opts Options) error {
	if err := validateOptions(opts); err != nil {
		return err
	}
	rd, err := sam.NewReader(in)
	if err != nil {
		return fmt.Errorf("mosdepth: open BAM: %w", err)
	}
	return runWithReader(rd, opts)
}

// validateOptions checks the option-level invariants that do not depend on the
// input stream, so they can be reported before any reader is opened. It is
// called by every entry point (Run, OpenAndRun) to keep error ordering stable
// regardless of input format.
func validateOptions(opts Options) error {
	if opts.ByBED != "" && opts.ByWindow > 0 {
		return ErrByConflict
	}
	if opts.Prefix == "" {
		return fmt.Errorf("mosdepth: empty output prefix")
	}
	// Reject a max-fragment-length cap that sits below the minimum, mirroring
	// upstream mosdepth's `max_len < min_len` guard. The port uses 0 to mean
	// "unset" for both bounds, so the check only fires when both are positive
	// (an actual cap below an actual floor).
	if opts.MaxFragLen > 0 && opts.MinFragLen > 0 && opts.MaxFragLen < opts.MinFragLen {
		return ErrBadFragLenBounds
	}
	return nil
}

// runWithReader executes the mosdepth pipeline against an already-opened
// alignment reader. Depth is computed identically whether rd decodes BAM,
// SAM, or CRAM, so CRAM input yields byte-for-byte the same outputs as the
// equivalent BAM. Callers must validate opts (see validateOptions) before
// opening rd.
func runWithReader(rd sam.Reader, opts Options) error {
	if err := validateOptions(opts); err != nil {
		return err
	}
	hdr := rd.Header()
	if hdr == nil {
		return fmt.Errorf("mosdepth: BAM has no header")
	}

	// A --chrom restriction that names a reference absent from the header is a
	// hard error upstream (check_chrom → exit 1), not a silent empty run.
	if opts.Chrom != "" {
		found := false
		for _, r := range hdr.Refs {
			if r.Name == opts.Chrom {
				found = true
				break
			}
		}
		if !found {
			return &ChromNotFoundError{Chrom: opts.Chrom}
		}
	}

	// Resolve regions per-chrom up front. perChromRegions[chrom] is the
	// sorted list of [beg0, end0) intervals to summarise; nil for "no
	// regions". An empty map means "no regions configured" — i.e. only
	// per-base output is requested.
	perChromRegions, regionNames, err := resolveRegions(hdr, opts)
	if err != nil {
		return err
	}

	// Per-base output is emitted only when --by isn't set (matches
	// upstream). When --by IS set, upstream skips per-base unless
	// explicitly requested; we follow the same rule. When --d4 is set the
	// per-base track is written to <prefix>.per-base.d4 in the dense D4
	// binary format instead of the bgzipped BED.
	wantPerBase := !opts.NoPerBase && len(perChromRegions) == 0

	var perBaseW *bedGzWriter
	var d4W *d4Writer
	if wantPerBase && opts.D4Output {
		chroms := make([]d4Chrom, 0, len(hdr.Refs))
		for _, r := range hdr.Refs {
			if opts.Chrom != "" && r.Name != opts.Chrom {
				continue
			}
			chroms = append(chroms, d4Chrom{Name: r.Name, Length: int64(r.Length)})
		}
		w, perr := newD4Writer(opts.Prefix+".per-base.d4", chroms)
		if perr != nil {
			return perr
		}
		d4W = w
	} else if wantPerBase {
		p, perr := newBedGzWriter(opts.Prefix + ".per-base.bed.gz")
		if perr != nil {
			return perr
		}
		perBaseW = p
	}

	var quantW *bedGzWriter
	var quantLabels []string
	if len(opts.Quantize) > 0 {
		q, qerr := newBedGzWriter(opts.Prefix + ".quantized.bed.gz")
		if qerr != nil {
			if perBaseW != nil {
				_ = perBaseW.Close()
			}
			if d4W != nil {
				_ = d4W.Close()
			}
			return qerr
		}
		quantW = q
		quantLabels = quantizeLabels(opts.Quantize)
	}

	var regionsW *bedGzWriter
	if len(perChromRegions) > 0 {
		p, perr := newBedGzWriter(opts.Prefix + ".regions.bed.gz")
		if perr != nil {
			if perBaseW != nil {
				_ = perBaseW.Close()
			}
			return perr
		}
		regionsW = p
	}

	var thresholdsW *bedGzWriter
	if len(opts.Thresholds) > 0 && len(perChromRegions) > 0 {
		p, perr := newBedGzWriter(opts.Prefix + ".thresholds.bed.gz")
		if perr != nil {
			if perBaseW != nil {
				_ = perBaseW.Close()
			}
			if regionsW != nil {
				_ = regionsW.Close()
			}
			return perr
		}
		thresholdsW = p
		if err := writeThresholdHeader(thresholdsW, opts.Thresholds); err != nil {
			return err
		}
	}

	// Bucket records by reference so we can scan in @SQ order. For each
	// chrom we build a covAccum, feed every record, then emit outputs and
	// drop the accum.
	perChromHist := map[string][]int64{}
	summaryRows := make([]summaryRow, 0, len(hdr.Refs))

	// Pull records once, grouping by chrom. We don't trust that the BAM
	// is sorted — buffering per-chrom slices keeps the algorithm correct
	// even on a name-sorted input. Memory is O(records that match
	// filters) which is acceptable for the typical mosdepth workload.
	byChrom, err := groupRecords(rd, opts)
	if err != nil {
		return err
	}

	chromOrder := make([]string, 0, len(hdr.Refs))
	for _, r := range hdr.Refs {
		if opts.Chrom != "" && r.Name != opts.Chrom {
			continue
		}
		chromOrder = append(chromOrder, r.Name)
	}

	for _, r := range hdr.Refs {
		if opts.Chrom != "" && r.Name != opts.Chrom {
			continue
		}
		recs := byChrom[r.Name]
		accum := newCovAccum(int(r.Length))
		accum.addRecords(recs, opts.FastMode, opts.FragmentMode)
		// Per-base emission: collapse runs of equal depth.
		hist := accumHistogram(accum)
		perChromHist[r.Name] = hist

		// Compute mean/min/max across the whole chromosome — also used in
		// the summary file.
		sum, _, minD, maxD := accum.regionStats(0, int(r.Length), nil, nil)
		row := summaryRow{
			chrom:  r.Name,
			length: int64(r.Length),
			bases:  sum,
			minD:   minD,
			maxD:   maxD,
		}
		if r.Length > 0 {
			row.mean = float64(sum) / float64(r.Length)
		}
		summaryRows = append(summaryRows, row)

		if perBaseW != nil {
			if err := emitPerBase(perBaseW, r.Name, accum); err != nil {
				return err
			}
		}
		if d4W != nil {
			if err := d4W.writeChrom(r.Name, d4DenseDepths(accum)); err != nil {
				return err
			}
		}
		if quantW != nil {
			if err := emitQuantized(quantW, r.Name, accum, opts.Quantize, quantLabels); err != nil {
				return err
			}
		}

		if regionsW != nil {
			ivs := perChromRegions[r.Name]
			for _, iv := range ivs {
				// When --use-median is set, build the depth histogram in the
				// same regionStats sweep that computes the thresholds, so the
				// median and mean share a single pass over the event list and
				// an identical per-base depth profile (see medianHist).
				var mh medianHist
				var emit func(start, end int, depth int32)
				if opts.UseMedian {
					emit = mh.addRun
				}
				bsum, perTh, mn, mx := accum.regionStats(iv.beg, iv.end, opts.Thresholds, emit)
				_ = mn
				_ = mx
				width := iv.end - iv.beg
				var stat float64
				if opts.UseMedian {
					// Upstream routes the regions depth column through
					// imean(), which returns the histogram median when
					// --use-median is set. Only this column changes.
					stat = mh.median()
				} else if width > 0 {
					stat = float64(bsum) / float64(width)
				}
				extras := []string{}
				if iv.name != "" {
					extras = append(extras, iv.name)
				}
				extras = append(extras, formatMean(stat))
				if err := regionsW.writeBED(r.Name, iv.beg, iv.end, extras...); err != nil {
					return err
				}
				if thresholdsW != nil {
					th := []string{}
					if iv.name != "" {
						th = append(th, iv.name)
					} else {
						th = append(th, fmt.Sprintf("%s:%d-%d", r.Name, iv.beg, iv.end))
					}
					for i := range opts.Thresholds {
						th = append(th, strconv.FormatInt(perTh[i], 10))
					}
					if err := thresholdsW.writeBED(r.Name, iv.beg, iv.end, th...); err != nil {
						return err
					}
				}
			}
		}
		// Free the accumulator's events explicitly to keep memory bounded.
		accum.events = nil
	}

	// Close output writers + build CSI indexes (matching upstream mosdepth,
	// which emits .csi alongside each bgzipped BED output).
	if perBaseW != nil {
		if err := perBaseW.Close(); err != nil {
			return err
		}
		if err := buildBedCsi(perBaseW.path); err != nil {
			return err
		}
	}
	if d4W != nil {
		if err := d4W.Close(); err != nil {
			return err
		}
	}
	if quantW != nil {
		if err := quantW.Close(); err != nil {
			return err
		}
		if err := buildBedCsi(quantW.path); err != nil {
			return err
		}
	}
	if regionsW != nil {
		if err := regionsW.Close(); err != nil {
			return err
		}
		if err := buildBedCsi(regionsW.path); err != nil {
			return err
		}
	}
	if thresholdsW != nil {
		if err := thresholdsW.Close(); err != nil {
			return err
		}
		if err := buildBedCsi(thresholdsW.path); err != nil {
			return err
		}
	}

	if err := writeDistribution(opts.Prefix+".mosdepth.global.dist.txt", perChromHist, chromOrder); err != nil {
		return err
	}
	if err := writeSummary(opts.Prefix+".mosdepth.summary.txt", summaryRows); err != nil {
		return err
	}
	_ = regionNames // silence linter; used implicitly via perChromRegions ordering above.
	return nil
}

// region is a single (chrom, beg0, end0[, name]) interval.
type region struct {
	beg, end int
	name     string
}

// resolveRegions returns the per-chromosome region list derived from
// opts.ByBED / opts.ByWindow.
//
// regionNames is the union of unique region names (for the regions file
// 4th column); when the BED has no name column, regionNames is nil.
func resolveRegions(hdr *sam.Header, opts Options) (map[string][]region, []string, error) {
	out := map[string][]region{}
	var names []string
	switch {
	case opts.ByBED != "":
		f, err := os.Open(opts.ByBED)
		if err != nil {
			return nil, nil, fmt.Errorf("mosdepth: open BED: %w", err)
		}
		defer f.Close()
		rd := bed.NewReader(f)
		for {
			rec, err := rd.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return nil, nil, fmt.Errorf("mosdepth: read BED: %w", err)
			}
			if opts.Chrom != "" && rec.Chrom != opts.Chrom {
				continue
			}
			out[rec.Chrom] = append(out[rec.Chrom], region{beg: rec.ChromStart, end: rec.ChromEnd, name: rec.Name})
			if rec.Name != "" {
				names = append(names, rec.Name)
			}
		}
	case opts.ByWindow > 0:
		size := opts.ByWindow
		for _, r := range hdr.Refs {
			if opts.Chrom != "" && r.Name != opts.Chrom {
				continue
			}
			n := int(r.Length)
			for s := 0; s < n; s += size {
				e := s + size
				if e > n {
					e = n
				}
				out[r.Name] = append(out[r.Name], region{beg: s, end: e})
			}
		}
	}
	return out, stringSliceUnique(names), nil
}

// disableMapqFastPath, when set, forces groupRecords onto the general
// keep-predicate path even when the fast path would otherwise apply. It
// exists solely so tests can prove the fast and general paths produce
// byte-identical output; production code never touches it.
var disableMapqFastPath bool

// mapqFastPath reports whether the MAPQ-free fast path may be taken for opts.
// Upstream mosdepth special-cases the common `--mapq 0` invocation (no MAPQ
// filter) by dropping the per-read MAPQ comparison from the hot loop. The
// fast path is purely a performance optimisation: its output is identical to
// the general path because, when MinMAPQ == 0, the MAPQ predicate
// `rec.MapQ < 0` is unsatisfiable and never rejects a read.
func mapqFastPath(opts Options) bool { return opts.MinMAPQ == 0 && !disableMapqFastPath }

// groupRecords drains rd into a map keyed by reference name, applying all
// configured read-level filters. Reads with unknown / "*" RNAME are
// dropped silently — they cannot contribute to depth on any reference.
//
// When no MAPQ filter is in effect (opts.MinMAPQ == 0, see mapqFastPath) the
// per-read keep predicate is bound once to keepRecordNoMapq, which omits the
// MAPQ comparison from the hot loop. This mirrors upstream mosdepth's
// `--mapq 0` fast path and is byte-for-byte equivalent to the general path.
func groupRecords(rd sam.Reader, opts Options) (map[string][]*sam.Record, error) {
	keep := keepRecord
	if mapqFastPath(opts) {
		keep = keepRecordNoMapq
	}
	out := map[string][]*sam.Record{}
	for {
		rec, err := rd.Read()
		if err == io.EOF {
			return out, nil
		}
		if err != nil {
			return nil, fmt.Errorf("mosdepth: read record: %w", err)
		}
		if !keep(rec, opts) {
			continue
		}
		out[rec.RName] = append(out[rec.RName], rec)
	}
}

// keepRecord applies every per-read filter configured on opts — including the
// MAPQ floor — and returns true if the read should contribute to depth.
func keepRecord(rec *sam.Record, opts Options) bool {
	if !keepRecordCommon(rec, opts) {
		return false
	}
	if opts.MinMAPQ > 0 && rec.MapQ < opts.MinMAPQ {
		return false
	}
	return keepRecordTail(rec, opts)
}

// keepRecordNoMapq is the `--mapq 0` fast-path keep predicate: it applies
// every filter except the MAPQ floor, which is a no-op when MinMAPQ == 0.
// Callers must only use it when mapqFastPath(opts) is true.
func keepRecordNoMapq(rec *sam.Record, opts Options) bool {
	if !keepRecordCommon(rec, opts) {
		return false
	}
	return keepRecordTail(rec, opts)
}

// keepRecordCommon applies the position, chromosome, and SAM-flag filters
// shared by both the general and fast-path keep predicates.
func keepRecordCommon(rec *sam.Record, opts Options) bool {
	if rec.Pos <= 0 || rec.RName == "" || rec.RName == "*" {
		return false
	}
	if opts.Chrom != "" && rec.RName != opts.Chrom {
		return false
	}
	if opts.ExcludeFlag != 0 && rec.Flag&opts.ExcludeFlag != 0 {
		return false
	}
	if opts.IncludeFlag != 0 && rec.Flag&opts.IncludeFlag != opts.IncludeFlag {
		return false
	}
	return true
}

// keepRecordTail applies the fragment-length and read-group filters that run
// after the MAPQ check in both keep predicates.
func keepRecordTail(rec *sam.Record, opts Options) bool {
	if opts.MinFragLen > 0 || opts.MaxFragLen > 0 {
		t := int(rec.TLen)
		if t < 0 {
			t = -t
		}
		if opts.MinFragLen > 0 && t < opts.MinFragLen {
			return false
		}
		if opts.MaxFragLen > 0 && t > opts.MaxFragLen {
			return false
		}
	}
	if len(opts.ReadGroups) > 0 {
		if !readGroupMatches(rec, opts.ReadGroups) {
			return false
		}
	}
	if opts.FragmentMode {
		// Upstream --fragment-mode counts only read1 of a proper,
		// non-supplementary pair; everything else is skipped before it can
		// contribute fragment coverage.
		if rec.Flag&sam.FlagRead2 != 0 ||
			rec.Flag&sam.FlagProperPair == 0 ||
			rec.Flag&sam.FlagSupplementary != 0 {
			return false
		}
	}
	return true
}

// readGroupMatches returns true when rec's RG (or OPS) aux tag is in the
// configured allow-list. The first element may take the form "OPS:X,Y" to
// indicate the OPS tag should be filtered instead of RG; in that case the
// remaining elements (after stripping the OPS: prefix on the first) are
// the allowed values.
func readGroupMatches(rec *sam.Record, allow []string) bool {
	tag := "RG"
	values := allow
	if len(allow) > 0 && strings.HasPrefix(allow[0], "OPS:") {
		tag = "OPS"
		// Reconstruct the list with the prefix stripped on element 0.
		stripped := append([]string{}, allow...)
		stripped[0] = strings.TrimPrefix(stripped[0], "OPS:")
		values = stripped
	}
	a, ok := rec.GetAux(tag)
	if !ok {
		return false
	}
	got, ok := a.String()
	if !ok {
		return false
	}
	for _, v := range values {
		if v == got {
			return true
		}
	}
	return false
}

// accumHistogram sweeps an accumulator and returns a depth histogram where
// hist[d] is the number of bases observed at exactly depth d.
//
// The histogram is the engine of writeDistribution; building it during the
// same pass we already do for per-base emission saves a second traversal.
func accumHistogram(a *covAccum) []int64 {
	hist := []int64{}
	a.emit(func(_ int, depth int32) {
		if int(depth) >= len(hist) {
			grown := make([]int64, int(depth)+1)
			copy(grown, hist)
			hist = grown
		}
		hist[depth]++
	})
	return hist
}

// emitPerBase emits per-base BED runs for a single chromosome to writer.
func emitPerBase(w *bedGzWriter, chrom string, a *covAccum) error {
	var emitErr error
	a.emitRuns(func(start, end int, depth int32) {
		if emitErr != nil {
			return
		}
		emitErr = w.writeBED(chrom, start, end, strconv.FormatInt(int64(depth), 10))
	})
	return emitErr
}

// ParseByArg parses the user-supplied `-b/--by` argument. If it parses as a
// positive integer, the result is treated as a window size; otherwise it is
// treated as a path to a BED file. Empty input means "no regions".
//
// Returns (windowSize, bedPath, err). At most one of windowSize/bedPath is
// non-zero/non-empty.
func ParseByArg(arg string) (int, string, error) {
	arg = strings.TrimSpace(arg)
	if arg == "" {
		return 0, "", nil
	}
	if v, err := strconv.Atoi(arg); err == nil {
		if v <= 0 {
			return 0, "", fmt.Errorf("mosdepth: --by window size must be positive, got %d", v)
		}
		return v, "", nil
	}
	// Treat as a file path; verify it exists for a friendlier error.
	if _, err := os.Stat(arg); err != nil {
		return 0, "", fmt.Errorf("mosdepth: --by file %q not found: %w", arg, err)
	}
	return 0, arg, nil
}

// ParseReadGroups splits a comma-separated --read-groups argument, leaving
// an optional "OPS:" prefix on the first element so the downstream code can
// detect the OPS-tag mode.
func ParseReadGroups(spec string) []string {
	if spec == "" {
		return nil
	}
	parts := strings.Split(spec, ",")
	for i, p := range parts {
		parts[i] = strings.TrimSpace(p)
	}
	// Drop empties.
	out := parts[:0]
	for _, p := range parts {
		if p != "" {
			out = append(out, p)
		}
	}
	return out
}

// readerForFile opens path for reading and returns a buffered reader. The
// caller is responsible for closing the returned io.Closer.
func readerForFile(path string) (io.Reader, io.Closer, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	return bufio.NewReader(f), f, nil
}

// OpenAndRun is the convenience entry point used by the CLI: open path,
// stream it through Run, and ensure the file handle is closed.
//
// The input format (SAM, BAM, or CRAM) is auto-detected from its leading
// bytes. CRAM input is decoded through pkg/htsgo/alnio, honouring
// opts.Fasta (-f/--fasta) and the REF_CACHE environment variable for
// reference resolution; depth output is identical to the equivalent BAM.
//
// When opts.Threads >= 2 and the input is BGZF-compressed (the usual BAM
// case), the BGZF blocks are decompressed concurrently across opts.Threads
// worker goroutines via bgzf.NewMultiReader. The decoded byte stream — and
// therefore every output file — is byte-for-byte identical to the
// single-threaded path regardless of the thread count; threading only affects
// decode throughput. CRAM carries its own block framing and is read
// single-threaded.
func OpenAndRun(path string, opts Options) error {
	if err := validateOptions(opts); err != nil {
		return err
	}
	r, c, err := readerForFile(path)
	if err != nil {
		return err
	}
	defer c.Close()

	// Sniff the leading bytes to detect CRAM, which needs the reference-aware
	// alnio decoder rather than sam.NewReader. SAM/BAM fall through to the
	// existing (optionally BGZF-threaded) path unchanged. The bufio reader is
	// reused below so the peeked bytes are not lost.
	br := bufio.NewReader(r)
	head, _ := br.Peek(16)
	if looksLikeCRAM(head) {
		rd, rerr := alnio.NewReaderWithReference(br, opts.Fasta)
		if rerr != nil {
			return fmt.Errorf("mosdepth: open CRAM: %w", rerr)
		}
		// The CRAM reader owns the reference-FASTA handle opened by
		// NewReaderWithReference; close it so the descriptor is not leaked
		// (alnio returns it typed as sam.Reader, but the concrete CRAM reader
		// implements io.Closer).
		if rc, ok := rd.(io.Closer); ok {
			defer rc.Close()
		}
		return runWithReader(rd, opts)
	}
	r = br

	if opts.Threads >= 2 {
		if looksLikeBGZF(head) {
			mr, merr := bgzf.NewMultiReader(br, opts.Threads)
			if merr != nil {
				return fmt.Errorf("mosdepth: open parallel BGZF: %w", merr)
			}
			defer mr.Close()
			return Run(mr, opts)
		}
		// Not BGZF (raw BAM or SAM): nothing to parallelise; fall through to
		// the sequential path with the peeked bytes preserved.
		return Run(br, opts)
	}
	return Run(r, opts)
}

// looksLikeCRAM reports whether b begins with the four-byte "CRAM" file
// magic. CRAM has its own container framing (it is never BGZF-wrapped at the
// file level), so this sniff is sufficient to route the input to the
// reference-aware alnio decoder instead of sam.NewReader.
func looksLikeCRAM(b []byte) bool {
	return len(b) >= 4 && b[0] == 'C' && b[1] == 'R' && b[2] == 'A' && b[3] == 'M'
}

// looksLikeBGZF reports whether b begins with a BGZF gzip header (gzip magic +
// the BC subfield). It mirrors the sniff used by the sam and iohelper packages
// and gates whether OpenAndRun engages the parallel BGZF reader.
func looksLikeBGZF(b []byte) bool {
	if len(b) < 16 {
		return false
	}
	if b[0] != 0x1f || b[1] != 0x8b || b[2] != 0x08 || b[3]&0x04 == 0 {
		return false
	}
	xlen := uint16(b[10]) | uint16(b[11])<<8
	if xlen < 6 {
		return false
	}
	return b[12] == 'B' && b[13] == 'C' && b[14] == 0x02 && b[15] == 0x00
}
