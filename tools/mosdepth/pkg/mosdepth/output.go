package mosdepth

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

// bedGzWriter wraps a BGZF writer plus a buffered writer on top, so callers
// can emit BED lines cheaply while still producing a valid bgzipped file
// indexable by tabix.
type bedGzWriter struct {
	f   *os.File
	bgz *bgzip.Writer
	bw  *bufio.Writer
	// path is retained so we can index it post-Close.
	path string
}

// newBedGzWriter opens path for writing and returns a wrapped writer that
// flushes its buffered text through a BGZF writer. The file is closed by
// Close.
func newBedGzWriter(path string) (*bedGzWriter, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	bgz := bgzip.NewWriter(f)
	return &bedGzWriter{f: f, bgz: bgz, bw: bufio.NewWriter(bgz), path: path}, nil
}

// writeBED emits one BED record: chrom\tstart\tend\textras...\n. The extras
// are joined by tabs and may be empty.
func (w *bedGzWriter) writeBED(chrom string, start, end int, extras ...string) error {
	if _, err := w.bw.WriteString(chrom); err != nil {
		return err
	}
	if err := w.bw.WriteByte('\t'); err != nil {
		return err
	}
	if _, err := w.bw.WriteString(strconv.Itoa(start)); err != nil {
		return err
	}
	if err := w.bw.WriteByte('\t'); err != nil {
		return err
	}
	if _, err := w.bw.WriteString(strconv.Itoa(end)); err != nil {
		return err
	}
	for _, e := range extras {
		if err := w.bw.WriteByte('\t'); err != nil {
			return err
		}
		if _, err := w.bw.WriteString(e); err != nil {
			return err
		}
	}
	return w.bw.WriteByte('\n')
}

// Close flushes the buffered writer, closes the BGZF stream (emitting the
// EOF block), and closes the underlying file.
func (w *bedGzWriter) Close() error {
	if err := w.bw.Flush(); err != nil {
		return err
	}
	if err := w.bgz.Close(); err != nil {
		return err
	}
	return w.f.Close()
}

// buildBedTbi runs the project's tabix builder against a freshly-closed
// bgzipped BED file at path, writing the .tbi alongside. Errors are
// returned to the caller; an empty file is allowed and produces an empty
// index.
func buildBedTbi(path string) error {
	cfg, err := tabix.PresetConfig(tabix.PresetBED)
	if err != nil {
		return err
	}
	idx, err := tabix.Build(path, cfg)
	if err != nil {
		return fmt.Errorf("tabix build %s: %w", path, err)
	}
	return idx.WriteFile(path + ".tbi")
}

// writeDistribution emits the cumulative depth-distribution file at path.
// histogram[d] is the number of bases observed at exactly depth d
// (length = max depth observed + 1).
//
// mosdepth's .global.dist.txt uses a "for each chrom plus 'total'" layout:
//
//	chrom<TAB>depth<TAB>proportion-of-bases-at-or-above-depth
//
// where proportion is the fraction of bases (across the entire chrom or
// genome) whose depth is >= the listed depth. The file always lists from
// the highest observed depth down to 0 so consumers can plot it directly.
func writeDistribution(path string, perChromHist map[string][]int64, chromOrder []string, forceEmit map[string]bool) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	bw := bufio.NewWriter(f)
	defer bw.Flush()

	// Aggregate "total" histogram. To match upstream mosdepth we only
	// roll up chromosomes that have at least one base of non-zero depth
	// — chroms whose histogram is `[length]` (no reads ever covered any
	// position) are excluded from the genome-wide row so the displayed
	// proportions reflect the covered subset of the reference.
	totalHist := []int64{}
	chromHasCoverage := func(h []int64) bool {
		for d := 1; d < len(h); d++ {
			if h[d] > 0 {
				return true
			}
		}
		return false
	}
	for _, name := range chromOrder {
		h := perChromHist[name]
		// Skip chromosomes with no coverage that weren't in
		// forceEmit (i.e. the BAM had no records there and the
		// user didn't target them via --chrom). Chroms in
		// forceEmit with no coverage still contribute their length
		// to `total` so the genome-wide proportions reflect them.
		if !chromHasCoverage(h) && !forceEmit[name] {
			continue
		}
		if len(h) > len(totalHist) {
			grown := make([]int64, len(h))
			copy(grown, totalHist)
			totalHist = grown
		}
		for i, c := range h {
			totalHist[i] += c
		}
	}
	emit := func(label string, hist []int64, force bool) error {
		if len(hist) == 0 {
			return nil
		}
		var totalBases int64
		for _, c := range hist {
			totalBases += c
		}
		if totalBases == 0 {
			return nil
		}
		// Match upstream: skip chroms whose only depth is zero (i.e.
		// no read ever covered any base), unless `force` says the
		// user explicitly wants this row (chrom had records pre-
		// filter or was named via --chrom). The `total` row is
		// similarly restricted: see chromHasCoverage above.
		var nonZero int64
		for d := 1; d < len(hist); d++ {
			nonZero += hist[d]
		}
		if nonZero == 0 && !force {
			return nil
		}
		// Cumulative-from-top: bases at depth >= d.
		var cum int64
		// Walk descending so each step adds hist[d] before printing the
		// row for d.
		props := make([]float64, len(hist))
		for d := len(hist) - 1; d >= 0; d-- {
			cum += hist[d]
			props[d] = float64(cum) / float64(totalBases)
		}
		for d := len(hist) - 1; d >= 0; d-- {
			if _, err := fmt.Fprintf(bw, "%s\t%d\t%s\n", label, d, formatProportion(props[d])); err != nil {
				return err
			}
		}
		return nil
	}
	for _, name := range chromOrder {
		if err := emit(name, perChromHist[name], forceEmit[name]); err != nil {
			return err
		}
	}
	// `total` row is only emitted if any chrom contributed to it.
	return emit("total", totalHist, len(totalHist) > 0)
}

// formatProportion renders a fractional proportion using the same width
// mosdepth upstream uses (4 decimal places, no trailing zeros trimmed) so
// downstream parsers see a stable layout.
func formatProportion(p float64) string {
	if p < 0 {
		p = 0
	}
	if p > 1 {
		p = 1
	}
	return strconv.FormatFloat(p, 'f', 2, 64)
}

// summaryRow is one row of the .mosdepth.summary.txt file.
type summaryRow struct {
	chrom  string
	length int64
	bases  int64
	mean   float64
	minD   int32
	maxD   int32
	// forceEmit, when true, makes writeSummary include this row even
	// if `bases == 0`. Used to mirror upstream mosdepth's behaviour of
	// emitting an all-zero row for chromosomes that did have records
	// in the BAM (which were then filtered out) or that the user
	// explicitly targeted via `--chrom`.
	forceEmit bool
}

// writeSummary emits the per-chromosome (plus total / total_region)
// summary file. Columns: chrom, length, bases, mean, min, max.
//
// Per-chromosome rows whose `bases == 0` and `!forceEmit` are skipped
// to match upstream mosdepth. The `total` row aggregates emitted
// non-`_region` rows; if `--by` was set so the input includes any
// `<chrom>_region` rows, a parallel `total_region` row is emitted
// right after `total` aggregating just those.
func writeSummary(path string, rows []summaryRow) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	bw := bufio.NewWriter(f)
	defer bw.Flush()
	if _, err := fmt.Fprintln(bw, "chrom\tlength\tbases\tmean\tmin\tmax"); err != nil {
		return err
	}
	var totLen, totBases int64
	var totMin, totMax int32
	totFirst := true
	var rtotLen, rtotBases int64
	var rtotMin, rtotMax int32
	rtotFirst := true
	hasRegion := false
	for _, r := range rows {
		if r.bases == 0 && !r.forceEmit {
			continue
		}
		isRegion := strings.HasSuffix(r.chrom, "_region")
		if isRegion {
			hasRegion = true
			rtotLen += r.length
			rtotBases += r.bases
			if rtotFirst {
				rtotMin = r.minD
				rtotMax = r.maxD
				rtotFirst = false
			} else {
				if r.minD < rtotMin {
					rtotMin = r.minD
				}
				if r.maxD > rtotMax {
					rtotMax = r.maxD
				}
			}
		} else {
			totLen += r.length
			totBases += r.bases
			if totFirst {
				totMin = r.minD
				totMax = r.maxD
				totFirst = false
			} else {
				if r.minD < totMin {
					totMin = r.minD
				}
				if r.maxD > totMax {
					totMax = r.maxD
				}
			}
		}
		mean := r.mean
		if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%s\t%d\t%d\n",
			r.chrom, r.length, r.bases, formatMean(mean), r.minD, r.maxD); err != nil {
			return err
		}
	}
	var totMean float64
	if totLen > 0 {
		totMean = float64(totBases) / float64(totLen)
	}
	if _, err := fmt.Fprintf(bw, "total\t%d\t%d\t%s\t%d\t%d\n",
		totLen, totBases, formatMean(totMean), totMin, totMax); err != nil {
		return err
	}
	if hasRegion {
		var rtotMean float64
		if rtotLen > 0 {
			rtotMean = float64(rtotBases) / float64(rtotLen)
		}
		if _, err := fmt.Fprintf(bw, "total_region\t%d\t%d\t%s\t%d\t%d\n",
			rtotLen, rtotBases, formatMean(rtotMean), rtotMin, rtotMax); err != nil {
			return err
		}
	}
	return nil
}

// formatMean renders a mean depth with two decimal places.
func formatMean(m float64) string {
	return strconv.FormatFloat(m, 'f', 2, 64)
}

// writeThresholdHeader writes the column header for a thresholds.bed.gz
// file: "#chrom\tstart\tend\tregion\t1X\t5X\t10X" given thresholds=[1,5,10].
func writeThresholdHeader(w *bedGzWriter, thresholds []int) error {
	parts := []string{"#chrom", "start", "end", "region"}
	for _, t := range thresholds {
		parts = append(parts, fmt.Sprintf("%dX", t))
	}
	if _, err := w.bw.WriteString(strings.Join(parts, "\t")); err != nil {
		return err
	}
	return w.bw.WriteByte('\n')
}

// parseThresholds parses a comma-separated list of non-negative integers
// from spec (e.g. "1,5,10,30"). An empty string yields a nil slice.
func parseThresholds(spec string) ([]int, error) {
	if spec == "" {
		return nil, nil
	}
	parts := strings.Split(spec, ",")
	out := make([]int, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p == "" {
			continue
		}
		v, err := strconv.Atoi(p)
		if err != nil {
			return nil, fmt.Errorf("mosdepth: bad threshold %q: %w", p, err)
		}
		if v < 0 {
			return nil, fmt.Errorf("mosdepth: threshold must be non-negative, got %d", v)
		}
		out = append(out, v)
	}
	// Keep them sorted ascending so 1X precedes 5X precedes 10X.
	sort.Ints(out)
	return out, nil
}

// ParseQuantize parses upstream mosdepth's `-q/--quantize` argument
// — a colon-separated, ascending list of integer cutoffs that defines
// the depth bins. For example `0:1:1000` defines three bins covering
// `[0,1)`, `[1,1000)`, and `[1000,+inf)`. Empty input yields a nil
// slice. Returns an error on non-integer entries or descending lists.
func ParseQuantize(spec string) ([]int, error) {
	spec = strings.TrimSpace(spec)
	if spec == "" {
		return nil, nil
	}
	parts := strings.Split(spec, ":")
	out := make([]int, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p == "" {
			continue
		}
		v, err := strconv.Atoi(p)
		if err != nil {
			return nil, fmt.Errorf("mosdepth: bad --quantize cutoff %q: %w", p, err)
		}
		out = append(out, v)
	}
	if len(out) < 2 {
		return nil, fmt.Errorf("mosdepth: --quantize needs at least two colon-separated cutoffs (got %q)", spec)
	}
	for i := 1; i < len(out); i++ {
		if out[i] <= out[i-1] {
			return nil, fmt.Errorf("mosdepth: --quantize cutoffs must be strictly ascending (got %v)", out)
		}
	}
	return out, nil
}

// resolveQuantizeLabels returns the label string for each closed bin
// `i` in `[0, len(cutoffs)-1)`. Bin `i` covers depth range
// `[cutoffs[i], cutoffs[i+1])`. The default label is the literal
// `"{cutoffs[i]}:{cutoffs[i+1]}"` string upstream mosdepth uses;
// the environment variable `MOSDEPTH_Q{i}` overrides it when set.
// Note that the open-ended top bin (`[cutoffs[N-1], +inf)`) is
// intentionally not assigned a label and is omitted from
// `.quantized.bed.gz` output to match upstream behaviour (the
// depth-too-low bin `[-inf, cutoffs[0])` is also dropped at the
// emission site in emitQuantized).
func resolveQuantizeLabels(cutoffs []int) []string {
	if len(cutoffs) < 2 {
		return nil
	}
	n := len(cutoffs) - 1
	labels := make([]string, n)
	for i := 0; i < n; i++ {
		envKey := fmt.Sprintf("MOSDEPTH_Q%d", i)
		if v, ok := osLookupEnv(envKey); ok {
			labels[i] = v
			continue
		}
		labels[i] = fmt.Sprintf("%d:%d", cutoffs[i], cutoffs[i+1])
	}
	return labels
}

// osLookupEnv is a tiny wrapper so tests can stub it out if needed.
func osLookupEnv(key string) (string, bool) { return os.LookupEnv(key) }

// quantizeBin returns the bin index for a depth value given an ascending
// cutoff list. Bin `i` covers `[cutoffs[i], cutoffs[i+1])`. Depths
// outside `[cutoffs[0], cutoffs[N-1])` return -1, signalling that the
// caller should skip emission (upstream mosdepth drops both the
// below-first-cutoff and at-or-above-last-cutoff runs from
// `.quantized.bed.gz`).
func quantizeBin(d int, cutoffs []int) int {
	if len(cutoffs) < 2 {
		return -1
	}
	if d < cutoffs[0] || d >= cutoffs[len(cutoffs)-1] {
		return -1
	}
	// Walk in order — bins are typically few. A binary search is
	// available but unnecessary at the expected scale.
	for i := 0; i < len(cutoffs)-1; i++ {
		if d >= cutoffs[i] && d < cutoffs[i+1] {
			return i
		}
	}
	return -1
}

// emitQuantized walks accum and emits one BED4 record per maximal run
// of consecutive bases that share the same bin index. Labels match the
// MOSDEPTH_Q{i} env-var override resolved by resolveQuantizeLabels.
// Runs whose bin is -1 (depth outside `[cutoffs[0], cutoffs[N-1])`)
// are skipped so the output matches upstream byte-for-byte.
func emitQuantized(w *bedGzWriter, chrom string, a *covAccum, cutoffs []int, labels []string) error {
	var emitErr error
	var runStart int = -1
	var runBin int = -1
	flush := func(end int) {
		if emitErr != nil || runStart < 0 || runBin < 0 {
			return
		}
		emitErr = w.writeBED(chrom, runStart, end, labels[runBin])
	}
	a.emit(func(pos int, depth int32) {
		if emitErr != nil {
			return
		}
		bin := quantizeBin(int(depth), cutoffs)
		if runStart < 0 {
			runStart = pos
			runBin = bin
			return
		}
		if bin != runBin {
			flush(pos)
			runStart = pos
			runBin = bin
		}
	})
	if emitErr != nil {
		return emitErr
	}
	if runStart >= 0 {
		end := a.refLen
		if end <= runStart {
			end = runStart + 1
		}
		flush(end)
	}
	return emitErr
}

// stringSliceUnique returns a stable de-duplicated copy of in.
func stringSliceUnique(in []string) []string {
	seen := map[string]struct{}{}
	out := make([]string, 0, len(in))
	for _, s := range in {
		if _, ok := seen[s]; ok {
			continue
		}
		seen[s] = struct{}{}
		out = append(out, s)
	}
	return out
}
