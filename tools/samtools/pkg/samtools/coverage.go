package samtools

import (
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
	// computeCoverageStats rather than globally here.
	baseQSumAtPos map[int]uint64
	baseQCntAtPos map[int]uint64
	mapQSum       uint64
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

	states := make([]*coverageRefState, len(hdr.Refs))
	for i, ref := range hdr.Refs {
		states[i] = &coverageRefState{
			length:        ref.Length,
			posDelta:      map[int]int{},
			baseQSumAtPos: map[int]uint64{},
			baseQCntAtPos: map[int]uint64{},
		}
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

		// Walk CIGAR, accumulating per-base depth and baseQ stats. Bases
		// below MinBaseQ contribute neither to depth nor baseq sums
		// (matches upstream).
		refPos := int(rec.Pos) - 1
		queryPos := 0
		for _, op := range rec.Cigar {
			length := int(op.Length())
			switch op.Op() {
			case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
				for i := 0; i < length; i++ {
					var q uint8 = 0xff
					if queryPos+i < len(rec.Qual) {
						q = rec.Qual[queryPos+i]
					}
					if opts.MinBaseQ > 0 && q < opts.MinBaseQ {
						continue
					}
					pos := refPos + i
					if q != 0xff {
						st.baseQSumAtPos[pos] += uint64(q)
						st.baseQCntAtPos[pos]++
					}
					st.posDelta[pos]++
					st.posDelta[pos+1]--
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

	// Resolve regions list — when empty, every @SQ counts as the full
	// reference.
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
		for _, r := range opts.Regions {
			rg, perr := region.ParseRegion(r)
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
		}
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
	for _, s := range specs {
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

// computeCoverageStats resolves the per-position depths into the covered-base /
// summed-coverage / baseQ totals for the half-open 0-based window
// [start-1, end). Positions with depth < minDepth are not counted as covered,
// mirroring upstream's `depth >= mindepth` gate; baseQ sums are folded in only
// at covered positions to match upstream's per-column accounting.
func computeCoverageStats(start, end int32, st *coverageRefState, minDepth int) coverageStats {
	cs := coverageStats{covDepth: walkDepths(start, end, st, minDepth)}
	for p, depth := range cs.covDepth {
		cs.covBases++
		cs.summedCoverage += uint64(depth)
		cs.baseQSum += st.baseQSumAtPos[p]
		cs.baseQCnt += st.baseQCntAtPos[p]
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
	cs := computeCoverageStats(start, end, st, minDepth)
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

	cs := computeCoverageStats(start, end, st, minDepth)

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
