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
func writeDistribution(path string, perChromHist map[string][]int64, chromOrder []string) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	bw := bufio.NewWriter(f)
	defer bw.Flush()

	// Aggregate "total" histogram.
	totalHist := []int64{}
	for _, name := range chromOrder {
		h := perChromHist[name]
		if len(h) > len(totalHist) {
			grown := make([]int64, len(h))
			copy(grown, totalHist)
			totalHist = grown
		}
		for i, c := range h {
			totalHist[i] += c
		}
	}
	emit := func(label string, hist []int64) error {
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
		if err := emit(name, perChromHist[name]); err != nil {
			return err
		}
	}
	return emit("total", totalHist)
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
}

// writeSummary emits the per-chromosome (plus total_region) summary file.
// Columns: chrom, length, bases, mean, min, max.
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
	first := true
	for _, r := range rows {
		totLen += r.length
		totBases += r.bases
		if first {
			totMin = r.minD
			totMax = r.maxD
			first = false
		} else {
			if r.minD < totMin {
				totMin = r.minD
			}
			if r.maxD > totMax {
				totMax = r.maxD
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
	_, err = fmt.Fprintf(bw, "total\t%d\t%d\t%s\t%d\t%d\n",
		totLen, totBases, formatMean(totMean), totMin, totMax)
	return err
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
