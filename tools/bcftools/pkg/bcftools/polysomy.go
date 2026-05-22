// Package bcftools — see doc.go. This file implements `bcftools
// polysomy`, the chromosomal copy-number / aneuploidy estimator.
// Upstream source: `reference_code/bcftools/polysomy.c` plus the
// peak-fitting engine `reference_code/bcftools/peakfit.c`.
//
// The algorithm is a faithful port of upstream's Gaussian-mixture
// peak-fit, not a heuristic:
//
//  1. For each chromosome a BAF histogram is built over `NBins` bins
//     (default 150).
//  2. init() smooths the histogram, isolates the homozygous RR (BAF≈0)
//     and AA (BAF≈1) peaks from the heterozygous RA band, trims a
//     mis-centred AA peak, and per-segment normalises the three peaks
//     so they are comparable in height. A chromosome with too few hets
//     is short-circuited to CN1 (LOH / monosomy) or CN-unknown.
//  3. fit() runs three candidate Gaussian-mixture fits over the
//     heterozygous band:
//     - CN2: one bounded Gaussian centred near 0.5.
//     - CN3: two symmetric Gaussians near 1/3 and 2/3.
//     - CN4: three Gaussians (a central 0.5 peak plus two symmetric
//     side peaks).
//     Each fit's goodness is `Σ|model-y|`. The lowest copy number whose
//     fit passes `--fit-th` and (for CN3/CN4) the symmetry / peak-size
//     checks is chosen, with `--cn-penalty` acting as a tiebreaker: a
//     higher CN is only accepted when its fit beats
//     `(1-cn_penalty)·previous_fit`.
//
// The non-linear least-squares solver GSL provides upstream
// (`gsl_multifit_fdfsolver_lmsder`) is ported in-tree as pure Go in
// `peakfit_lm.go`; the peak models and Monte-Carlo restart driver are
// in `peakfit.go`. No third-party dependency is introduced.
//
// Validation: upstream ships no golden fixture for `polysomy` (its only
// output is per-chromosome PNG plots plus a dist.dat dump under
// `--output-dir`), so byte-for-byte parity cannot be demonstrated. The
// port is validated with unit tests of the LM solver against analytic
// curves, unit tests of every peak model, and hand-constructed BAF
// distributions for the canonical karyotypes (clean diploid → CN2,
// clear trisomy → CN3); see polysomy_test.go and peakfit_test.go.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Upstream default knob values (polysomy.c:main_polysomy).
const (
	polysomyDefaultNBins        = 150
	polysomyDefaultFitTh        = 3.3
	polysomyDefaultCnPenalty    = 0.7
	polysomyDefaultPeakSymmetry = 0.5
	polysomyDefaultMinPeakSize  = 0.1
	polysomyDefaultMinFraction  = 0.1
	polysomyDefaultSmooth       = -3
)

// PolysomyOptions controls a polysomy run. Zero-valued fields are
// replaced with upstream defaults by Polysomy / PolysomyFile.
type PolysomyOptions struct {
	// Sample restricts the analysis to a single sample name. Empty
	// means "all samples in the input". Upstream requires -s when the
	// VCF has more than one sample; the rule is applied by the CLI
	// runner.
	Sample string
	// SamplesFile is the upstream `-S/--samples-file` shortcut; a flat
	// list, one name per line.
	SamplesFile string
	// Regions / Targets / RegionsFile / TargetsFile match the rest of
	// bcftools: a chromosome-name post-filter.
	Regions     []string
	Targets     []string
	RegionsFile string
	TargetsFile string
	// IncludeExpr / ExcludeExpr are accepted at the CLI but applied
	// only at the record level (no per-sample filtering).
	IncludeExpr string
	ExcludeExpr string

	// CnPenalty is upstream `-c/--cn-penalty` (default 0.7). A higher
	// copy number is only accepted when its fit beats
	// (1-CnPenalty)·previous_fit.
	CnPenalty float64
	// FitTh is upstream `-f/--fit-th` (default 3.3): a candidate CN is
	// rejected when its Σ|model-y| fit exceeds this threshold.
	FitTh float64
	// PeakSymmetry is upstream `-p/--peak-symmetry` (default 0.5): the
	// minimum ratio of the smaller to the larger fitted peak area for a
	// CN3/CN4 call to be accepted.
	PeakSymmetry float64
	// MinPeakSize is upstream `-b/--peak-size` (default 0.1): the
	// minimum relative size of a CN4 side peak.
	MinPeakSize float64
	// MinFraction is upstream `-m/--min-fraction` (default 0.1): the
	// minimum distinguishable fraction of aberrant cells, controlling
	// how close to 0.5 the fitted side peaks may sit.
	MinFraction float64
	// IncludeAA is upstream `-i/--include-aa`: also fit the homozygous
	// AA peak when scoring CN2 and CN4.
	IncludeAA bool
	// NBins is the histogram bin count (hidden upstream `-n/--nbins`,
	// default 150).
	NBins int
	// Smooth is the smoothing half-window control (hidden upstream
	// `-S/--smooth`, default -3); a positive value also overwrites the
	// histogram with its smoothed form.
	Smooth int
	// RaRrScaling enables upstream's per-segment RA/RR/AA normalisation
	// (hidden upstream `--ra-rr-scaling`, default on). The CLI's
	// `--ra-rr-scaling` flag DISABLES it, matching upstream.
	RaRrScaling bool
	// ForceCN (upstream hidden `--force-cn`) tags every chromosome with
	// the requested copy number, bypassing the fit.
	ForceCN int
}

// PolysomyResult is the per-sample-per-chromosome CN call.
type PolysomyResult struct {
	Sample    string  // sample name
	Chrom     string  // chromosome / contig
	NHet      int     // number of heterozygous BAF values collected
	MeanBAF   float64 // arithmetic mean of BAF for the sample/chrom
	MedianBAF float64 // median of BAF for the sample/chrom
	CN        float64 // copy-number call (>=1, or -1 for unknown)
	Fit       float64 // Σ|model-y| of the winning fit (0 for short-circuits)
}

// PolysomySummary is what the file-level entry points return.
type PolysomySummary struct {
	Results []PolysomyResult
}

// applyPolysomyDefaults fills zero-valued knobs with upstream defaults.
// RaRrScaling defaults to enabled; the CLI flag of the same name turns
// it off, so callers that want it disabled must set the field after
// this call (the CLI runner does so explicitly).
func applyPolysomyDefaults(opts *PolysomyOptions) {
	if opts.CnPenalty == 0 {
		opts.CnPenalty = polysomyDefaultCnPenalty
	}
	if opts.FitTh == 0 {
		opts.FitTh = polysomyDefaultFitTh
	}
	if opts.PeakSymmetry == 0 {
		opts.PeakSymmetry = polysomyDefaultPeakSymmetry
	}
	if opts.MinPeakSize == 0 {
		opts.MinPeakSize = polysomyDefaultMinPeakSize
	}
	if opts.MinFraction == 0 {
		opts.MinFraction = polysomyDefaultMinFraction
	}
	if opts.NBins == 0 {
		opts.NBins = polysomyDefaultNBins
	}
	if opts.Smooth == 0 {
		// Upstream's default smoothing knob is -3, which yields a
		// 7-bin sliding window (win = |smooth|*2+1); a positive value
		// would additionally overwrite the histogram with its smoothed
		// form, so -3 keeps smoothing read-only. Defaulted here (rather
		// than re-defaulted inside polysomyInitDist) so the trap of two
		// places defining the same default is avoided.
		opts.Smooth = polysomyDefaultSmooth
	}
}

// PolysomyFile streams VCF/BCF from path and writes the per-sample ×
// per-chrom TSV to out.
func PolysomyFile(path string, out io.Writer, opts PolysomyOptions) (PolysomySummary, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return PolysomySummary{}, fmt.Errorf("bcftools polysomy: open %s: %w", path, err)
	}
	defer in.Close()
	return Polysomy(in, out, opts)
}

// Polysomy is the streaming entry point. It collects BAF values per
// (sample, chromosome), builds and fits each chromosome's BAF
// distribution, and emits CN calls.
func Polysomy(in io.Reader, out io.Writer, opts PolysomyOptions) (PolysomySummary, error) {
	applyPolysomyDefaults(&opts)

	hdr, variants, err := readAllVariants(in)
	if err != nil {
		return PolysomySummary{}, fmt.Errorf("bcftools polysomy: %w", err)
	}

	// Resolve sample subset.
	samples := hdr.Samples
	if opts.Sample != "" {
		samples = []string{opts.Sample}
	} else if opts.SamplesFile != "" {
		names, err := LoadSamplesFile(opts.SamplesFile)
		if err != nil {
			return PolysomySummary{}, fmt.Errorf("bcftools polysomy: %w", err)
		}
		samples = names
	}
	idxBySample := make(map[string]int, len(hdr.Samples))
	for i, s := range hdr.Samples {
		idxBySample[s] = i
	}
	for _, s := range samples {
		if _, ok := idxBySample[s]; !ok {
			return PolysomySummary{}, fmt.Errorf("bcftools polysomy: sample %q not in input", s)
		}
	}

	regionSet := buildPolysomyRegionSet(opts.Regions, opts.Targets, opts.RegionsFile, opts.TargetsFile)

	// Per-(sample, chrom) BAF accumulator.
	type key struct{ sample, chrom string }
	bafs := make(map[key][]float64)
	chromOrder := []string{}
	chromSeen := map[string]bool{}

	for _, v := range variants {
		if !polysomyKeepRecord(v, regionSet) {
			continue
		}
		if !chromSeen[v.Chrom] {
			chromSeen[v.Chrom] = true
			chromOrder = append(chromOrder, v.Chrom)
		}
		for _, sn := range samples {
			idx, ok := idxBySample[sn]
			if !ok {
				continue
			}
			baf, ok := readPolysomyBAF(v, idx)
			if !ok {
				continue
			}
			bafs[key{sample: sn, chrom: v.Chrom}] = append(bafs[key{sample: sn, chrom: v.Chrom}], baf)
		}
	}

	sum := PolysomySummary{}
	w := bufio.NewWriter(out)
	defer w.Flush()
	if _, err := fmt.Fprintln(w, "# sample\tchrom\tn_het\tmean_baf\tmedian_baf\tcn_call"); err != nil {
		return sum, err
	}
	for _, sn := range samples {
		for _, chrom := range chromOrder {
			vs := bafs[key{sample: sn, chrom: chrom}]
			cn, fit := polysomyCNCall(vs, opts)
			res := PolysomyResult{
				Sample:    sn,
				Chrom:     chrom,
				NHet:      len(vs),
				MeanBAF:   polysomyMean(vs),
				MedianBAF: polysomyMedian(vs),
				CN:        cn,
				Fit:       fit,
			}
			sum.Results = append(sum.Results, res)
			if _, err := fmt.Fprintf(w, "%s\t%s\t%d\t%.4f\t%.4f\t%s\n",
				res.Sample, res.Chrom, res.NHet, res.MeanBAF, res.MedianBAF, polysomyCNText(res.CN)); err != nil {
				return sum, err
			}
		}
	}
	return sum, nil
}

// polysomyKeepRecord applies the region / target post-filter. Empty set
// → keep everything.
func polysomyKeepRecord(v *vcf.Variant, regions map[string]bool) bool {
	if len(regions) == 0 {
		return true
	}
	return regions[v.Chrom]
}

// buildPolysomyRegionSet folds region / target lists into a CHROM-name
// allowlist. Only the chrom-name component is honoured (post-filter
// approximation; the per-base interval filter is a follow-up — see
// docs/PARITY_ROADMAP.md#bcftools).
func buildPolysomyRegionSet(regions, targets []string, regionsFile, targetsFile string) map[string]bool {
	out := map[string]bool{}
	addRegion := func(r string) {
		if i := strings.IndexByte(r, ':'); i >= 0 {
			r = r[:i]
		}
		r = strings.TrimSpace(r)
		if r != "" {
			out[r] = true
		}
	}
	for _, r := range regions {
		addRegion(r)
	}
	for _, r := range targets {
		addRegion(r)
	}
	if regionsFile != "" {
		if regs, err := LoadRegionsFile(regionsFile); err == nil {
			for _, r := range regs {
				addRegion(r)
			}
		}
	}
	if targetsFile != "" {
		if regs, err := LoadRegionsFile(targetsFile); err == nil {
			for _, r := range regs {
				addRegion(r)
			}
		}
	}
	return out
}

// readPolysomyBAF returns the heterozygous-site B-allele frequency for
// the sample at idx. An explicit FORMAT/BAF field is preferred
// (mirroring upstream's requirement); otherwise BAF is synthesised from
// FORMAT/AD = REF,ALT counts. Returns (baf, false) for non-het sites,
// missing data, or non-numeric fields.
func readPolysomyBAF(v *vcf.Variant, idx int) (float64, bool) {
	if idx < 0 || idx >= len(v.Samples) {
		return 0, false
	}
	data := v.Samples[idx].Data
	if data == nil {
		return 0, false
	}
	if gt, ok := data["GT"]; ok {
		if !polysomyIsHet(gt) {
			return 0, false
		}
	}
	if raw, ok := data["BAF"]; ok {
		if f, ok := polysomyParseFloat(raw); ok {
			return f, true
		}
	}
	if raw, ok := data["AD"]; ok {
		ad := strings.Split(raw, ",")
		if len(ad) >= 2 {
			refN, ok1 := polysomyParseFloat(ad[0])
			altN, ok2 := polysomyParseFloat(ad[1])
			if ok1 && ok2 && (refN+altN) > 0 {
				return altN / (refN + altN), true
			}
		}
	}
	return 0, false
}

// polysomyIsHet returns true if gt looks heterozygous. Both `/`- and
// `|`-separated diploid genotypes are accepted; everything else
// (missing, homozygous, polyploid) returns false.
func polysomyIsHet(gt string) bool {
	if gt == "" || gt == "." {
		return false
	}
	gt = strings.ReplaceAll(gt, "|", "/")
	parts := strings.Split(gt, "/")
	if len(parts) != 2 {
		return false
	}
	if parts[0] == "." || parts[1] == "." {
		return false
	}
	return parts[0] != parts[1]
}

// polysomyParseFloat parses a possibly-missing FORMAT field value.
func polysomyParseFloat(s string) (float64, bool) {
	s = strings.TrimSpace(s)
	if s == "" || s == "." {
		return 0, false
	}
	var f float64
	if _, err := fmt.Sscanf(s, "%f", &f); err != nil {
		return 0, false
	}
	return f, true
}

// polysomyMean is the arithmetic mean (0 for empty input).
func polysomyMean(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	var s float64
	for _, x := range xs {
		s += x
	}
	return s / float64(len(xs))
}

// polysomyMedian is the median (0 for empty input). For an even number
// of elements the average of the two middle values is returned.
func polysomyMedian(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	sorted := append([]float64(nil), xs...)
	sort.Float64s(sorted)
	n := len(sorted)
	if n%2 == 1 {
		return sorted[n/2]
	}
	return 0.5 * (sorted[n/2-1] + sorted[n/2])
}

// polysomyCNCall builds the BAF histogram for one chromosome's values
// and runs the upstream Gaussian-mixture peak fit, returning the CN
// call and the goodness of the winning fit. ForceCN, when non-zero,
// bypasses the fit. An empty BAF set is reported as CN1 (LOH / no
// informative hets).
func polysomyCNCall(bafs []float64, opts PolysomyOptions) (cn, fit float64) {
	if opts.ForceCN != 0 {
		return float64(opts.ForceCN), 0
	}
	if len(bafs) == 0 {
		return 1.0, 0
	}
	dist := buildBAFDist(bafs, opts)
	if dist.copyNumber != 0 {
		return float64(dist.copyNumber), 0
	}
	return fitBAFDist(dist, opts)
}

// polysomyCNText turns the float CN call into upstream's printable
// form. -1 is the "unknown" sentinel; everything else is printed with
// two decimals (matching upstream's `%.2f`).
func polysomyCNText(cn float64) string {
	if cn < 0 {
		return "?"
	}
	return fmt.Sprintf("%.2f", cn)
}

// bafDist is the per-chromosome BAF distribution after histogram
// binning, smoothing, RR/AA isolation and normalisation — the Go
// analogue of upstream's dist_t.
type bafDist struct {
	xvals      []float64 // bin centres in [0,1]
	yvals      []float64 // (smoothed, normalised) bin counts
	nvals      int       // valid prefix of xvals/yvals (AA peak may trim it)
	irr        int       // index where the heterozygous band starts
	ira        int       // index of the band centre (≈0.5)
	iaa        int       // index where the AA homozygous peak starts
	copyNumber int       // 0 => run the fit; 1 => CN1; -1 => unknown
}

// buildBAFDist bins the BAF values into a histogram and runs the
// upstream init_dist pipeline (smoothing, peak isolation, per-segment
// normalisation).
func buildBAFDist(bafs []float64, opts PolysomyOptions) *bafDist {
	nbins := opts.NBins
	if nbins < 2 {
		nbins = polysomyDefaultNBins
	}
	d := &bafDist{
		xvals: make([]float64, nbins),
		yvals: make([]float64, nbins),
		nvals: nbins,
	}
	for i := 0; i < nbins; i++ {
		d.xvals[i] = float64(i) / float64(nbins-1)
	}
	for _, b := range bafs {
		if b < 0 {
			b = 0
		} else if b > 1 {
			b = 1
		}
		bin := int(b * float64(nbins-1))
		if bin >= nbins {
			bin = nbins - 1
		}
		d.yvals[bin]++
	}
	polysomyInitDist(d, opts)
	return d
}

// polysomyInitDist ports polysomy.c:init_dist. It smooths the
// histogram, locates the RR/AA gaps, trims a mis-centred AA peak,
// decides whether the heterozygous peak exists at all (short-circuiting
// to CN1 / unknown), and per-segment normalises the three peaks.
func polysomyInitDist(d *bafDist, opts PolysomyOptions) {
	n := d.nvals
	// opts.Smooth is defaulted to -3 by applyPolysomyDefaults, so it is
	// never 0 here; abs(-3)*2+1 == 7 reproduces upstream's smooth==0
	// fallback window exactly.
	smooth := opts.Smooth
	win := 7
	if smooth != 0 {
		win = abs(smooth)*2 + 1 // must be odd
	}
	hwin := win / 2

	// Sliding-window mean, matching upstream's edge handling exactly.
	tmp := make([]float64, n)
	avg := d.yvals[0]
	tmp[0] = d.yvals[0]
	for i := 1; i < hwin; i++ {
		avg += d.yvals[2*i-1]
		tmp[i] = avg / float64(2*i+1)
	}
	avg = 0
	for i := 0; i < n; i++ {
		avg += d.yvals[i]
		if i >= win-1 {
			tmp[i-hwin] = avg / float64(win)
			avg -= d.yvals[i-win+1]
		}
	}
	hw := hwin
	for i := n - hw; i < n; i++ {
		avg -= d.yvals[i-hw]
		hw--
		tmp[i] = avg / float64(2*hw+1)
		avg -= d.yvals[i-hw]
	}

	// Find the RR gap in the left half and the AA gap in the right half.
	irr := 0
	for i := 0; i < n/2; i++ {
		if tmp[i] < tmp[irr] {
			irr = i
		}
	}
	iaa := n - 1
	for i := n - 1; i >= n/2; i-- {
		if tmp[i] < tmp[iaa] {
			iaa = i
		}
	}
	irr += int(float64(win) * 0.5)
	iaa += int(float64(win) * 0.5)
	if iaa >= n {
		iaa = n - 1
	}
	if irr >= iaa {
		// Upstream errors out here; we cannot, so fall back to
		// unknown-CN rather than aborting the whole run.
		d.copyNumber = -1
		d.irr, d.iaa, d.ira = 0, n-1, n/2
		return
	}
	if smooth > 0 {
		copy(d.yvals, tmp)
	}

	// Trim a mis-centred AA peak: the AA peak is occasionally closer to
	// the centre than to 1.0, so chop everything past its true maximum.
	imaxAA := iaa
	for i := iaa; i < n; i++ {
		if d.yvals[imaxAA] < d.yvals[i] {
			imaxAA = i
		}
	}
	d.nvals = imaxAA + 1
	n = d.nvals
	if iaa >= d.nvals {
		iaa = d.nvals - 1
	}

	// Per-segment maxima and sums for the RR / RA / AA bands.
	var maxRR, maxAA, maxRA, srr, saa, sra float64
	for i := 0; i < irr; i++ {
		srr += d.yvals[i]
		if maxRR < d.yvals[i] {
			maxRR = d.yvals[i]
		}
	}
	for i := irr; i <= iaa; i++ {
		sra += d.yvals[i]
		if maxRA < d.yvals[i] {
			maxRA = d.yvals[i]
		}
	}
	for i := iaa + 1; i < n; i++ {
		saa += d.yvals[i]
		if maxAA < d.yvals[i] {
			maxAA = d.yvals[i]
		}
	}

	if !opts.RaRrScaling {
		maxRA, maxAA = maxRR, maxRR
	}
	// Decide whether the heterozygous peak exists at all, matching
	// upstream's two-branch test. cDiv reproduces C's IEEE division so
	// that ratio-on-zero yields +Inf exactly as the C code relies on.
	raRr := cDiv(sra, srr)
	aaRa := cDiv(saa, sra)
	switch {
	case sra == 0 || (raRr < 0.1 && aaRa > 1.0):
		// Too few hets: CN1 (monosomy / LOH).
		maxRA = maxAA
		d.copyNumber = 1
	case raRr < 0.1 || aaRa > 1.0:
		maxRA = maxAA
		d.copyNumber = -1 // unknown copy number
	}
	if maxRR != 0 {
		for i := 0; i < irr; i++ {
			d.yvals[i] /= maxRR
		}
	}
	if maxRA != 0 {
		for i := irr; i <= iaa; i++ {
			d.yvals[i] /= maxRA
		}
	}
	if maxAA != 0 {
		for i := iaa + 1; i < n; i++ {
			d.yvals[i] /= maxAA
		}
	}

	d.irr = irr
	d.iaa = iaa
	d.ira = int(float64(n) * 0.5)
	if d.ira >= n {
		d.ira = n - 1
	}
}

// abs is an int absolute value (math.Abs works on float64 only).
func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// cDiv reproduces C's IEEE-754 division semantics for non-negative
// histogram sums: a finite numerator divided by zero yields +Inf
// (rather than panicking), and 0/0 yields NaN — both of which the
// upstream init_dist ratio tests rely on.
func cDiv(a, b float64) float64 {
	if b == 0 {
		if a == 0 {
			return math.NaN()
		}
		return math.Inf(1)
	}
	return a / b
}

// minRatio returns the smaller-over-larger ratio of two non-negative
// values, the symmetry measure upstream uses throughout fit_curves.
// When the larger value is zero the ratio is defined as zero.
func minRatio(a, b float64) float64 {
	if a > b {
		a, b = b, a
	}
	if b == 0 {
		return 0
	}
	return a / b
}

// fitBAFDist ports the per-chromosome body of polysomy.c:fit_curves: it
// runs the CN2, CN3 and CN4 candidate fits over the heterozygous band
// and picks the lowest copy number that passes the goodness, symmetry
// and peak-size checks, using cn_penalty as a tiebreaker.
func fitBAFDist(d *bafDist, opts PolysomyOptions) (cn, fit float64) {
	const nmc = 50
	pf := newPeakfit()

	irr, ira, iaa := d.irr, d.ira, d.iaa
	xrr := d.xvals[irr]
	xra := d.xvals[ira]
	xaa := d.xvals[iaa]
	xmax := d.xvals[d.nvals-1]

	nrrAA := iaa - irr + 1
	nrrRA := ira - irr + 1
	naaMax := d.nvals - iaa
	if nrrAA < 1 || nrrRA < 1 || naaMax < 1 {
		return -1, math.Inf(1)
	}

	// Sub-slices of the histogram, matching upstream's pointer offsets:
	//   xrr_vals / yrr_vals start at irr; xaa_vals / yaa_vals at iaa.
	xrrA := d.xvals[irr : irr+nrrAA]
	yrrA := d.yvals[irr : irr+nrrAA]
	xrrR := d.xvals[irr : irr+nrrRA]
	yrrR := d.yvals[irr : irr+nrrRA]
	xaaM := d.xvals[iaa : iaa+naaMax]
	yaaM := d.yvals[iaa : iaa+naaMax]

	// --- CN2: one bounded Gaussian near 0.5 (+ optional AA exp) ---
	var cn2aaFit float64
	if opts.IncludeAA {
		pf.reset()
		pf.addExp(1.0, 1.0, 0.2, 5)
		pf.setMC(0.01, 0.3, 2, nmc)
		pf.setMC(0.05, 1.0, 0, nmc)
		cn2aaFit = pf.run(xaaM, yaaM)
	}
	pf.reset()
	pf.addBoundedGaussian(1.0, 0.5, 0.03, 0.45, 0.55, 7)
	pf.setMC(0.01, 0.3, 2, nmc)
	pf.setMC(0.05, 1.0, 0, nmc)
	cn2raFit := pf.run(xrrA, yrrA)
	cn2Fit := cn2raFit + cn2aaFit

	// --- CN3: fit two bounded Gaussians, enforce symmetry, refit ---
	cn3aaFit := cn2aaFit
	minDx3 := 0.5 - 1.0/(opts.MinFraction+2)
	pf.reset()
	pf.addBoundedGaussian(1.0, 1.0/3.0, 0.03, xrr, xra-minDx3, 7)
	pf.setMC(xrr, xra-minDx3, 1, nmc)
	pf.addBoundedGaussian(1.0, 2.0/3.0, 0.03, xra+minDx3, xaa, 7)
	pf.setMC(xra+minDx3, xaa, 1, nmc)
	pf.run(xrrA, yrrA)
	s0, c0, sig0 := pf.getParams(0)
	s1, c1, sig1 := pf.getParams(1)
	cn3dx := (0.5 - c0 + c1 - 0.5) * 0.5
	if cn3dx > 0.5/3 {
		cn3dx = 0.5 / 3 // CN3 peaks must not be separated by more than 1/3
	}
	pf.reset()
	pf.addGaussian(s0, 0.5-cn3dx, sig0, 5)
	pf.addGaussian(s1, 0.5+cn3dx, sig1, 5)
	cn3raFit := pf.run(xrrA, yrrA)
	rraS, rraC, rraSig := pf.getParams(0)
	raaS, _, raaSig := pf.getParams(1)
	cn3dy := minRatio(rraS*rraS, raaS*raaS)
	cn3frac := (1 - 2*rraC) / rraC
	cn3Fit := cn3raFit + cn3aaFit
	// Exclude peaks far too broad or far too narrow.
	if rraSig > 0.3 || raaSig > 0.3 || rraSig < 1e-2 || raaSig < 1e-2 {
		cn3Fit = math.Inf(1)
	}

	// --- CN4: a central 0.5 peak plus two symmetric side peaks ---
	var cn4aaFit float64
	if opts.IncludeAA {
		pf.reset()
		pf.addExp(0.5, 1.0, 0.2, 5)
		pf.setMC(0.01, 0.3, 2, nmc)
		pf.addBoundedGaussian(0.4, (xaa+xmax)*0.5, 2e-2, xaa, xmax, 7)
		pf.setMC(xaa, xmax, 1, nmc)
		cn4aaFit = pf.run(xaaM, yaaM)
	}
	minDx4 := 0.25 * opts.MinFraction
	pf.reset()
	pf.addGaussian(1.0, 0.5, 0.03, 5)
	pf.addBoundedGaussian(0.6, 0.3, 0.03, xrr, xra-minDx4, 7)
	pf.setMC(xrr, xra-minDx4, 2, nmc)
	pf.run(xrrR, yrrR)
	raS, _, raSig := pf.getParams(0)
	rrS, rrC, rrSig := pf.getParams(1)
	cn4dx := 0.5 - rrC
	if cn4dx > 0.25 {
		cn4dx = 0.25 // CN4 side peaks must not be separated by more than 0.5
	}
	pf.reset()
	pf.addGaussian(raS, 0.5, raSig, 5)
	pf.addGaussian(rrS, 0.5-cn4dx, rrSig, 5)
	pf.addGaussian(rrS, 0.5+cn4dx, rrSig, 5)
	pf.setMC(0.1, raS, 0, nmc)
	pf.setMC(0.01, 0.1, 2, nmc)
	cn4raFit := pf.run(xrrA, yrrA)
	cn4raS, _, cn4raSig := pf.getParams(0)
	cn4rrS, cn4rrC, cn4rrSig := pf.getParams(1)
	cn4aaS, cn4aaC, cn4aaSig := pf.getParams(2)
	cn4RAraSize := math.Inf(1)
	if cn4raS != 0 {
		cn4RAraSize = cn4raS * cn4raS
	}
	cn4RArrSize := cn4rrS * cn4rrS
	cn4RAaaSize := cn4aaS * cn4aaS
	cn4rrDy := minRatio(cn4RArrSize, cn4RAraSize)
	cn4aaDy := minRatio(cn4RAaaSize, cn4RAraSize)
	cn4dy := minRatio(cn4rrDy, cn4aaDy)
	cn4ymin := cn4RAaaSize / cn4RAraSize
	if cn4RArrSize < cn4RAaaSize {
		cn4ymin = cn4RArrSize / cn4RAraSize
	}
	cn4dx = (cn4aaC - 0.5) - (0.5 - cn4rrC)
	cn4frac := cn4aaC - cn4rrC
	cn4Fit := cn4raFit + cn4aaFit
	if cn4raSig > 0.3 || cn4rrSig > 0.3 || cn4aaSig > 0.3 ||
		cn4raSig < 1e-2 || cn4rrSig < 1e-2 || cn4aaSig < 1e-2 {
		cn4Fit = math.Inf(1)
	}

	// --- Choose the best match (polysomy.c:fit_curves decision) ---
	cn2OK := cn2Fit <= opts.FitTh
	cn3OK := cn3Fit <= opts.FitTh && cn3dy >= opts.PeakSymmetry
	cn4OK := cn4Fit <= opts.FitTh && cn4ymin >= opts.MinPeakSize &&
		cn4dy >= opts.PeakSymmetry && cn4dx <= 0.1

	cn = -1
	fit = cn2Fit
	if cn2OK {
		cn = 2
		fit = cn2Fit
	}
	if cn3OK {
		// cn_penalty as a tiebreaker: a CN3 fit must beat
		// (1-cn_penalty)·cn2_fit to be preferred over CN2.
		if cn < 0 || cn3Fit < (1-opts.CnPenalty)*fit {
			cn = 2 + cn3frac
			fit = cn3Fit
		}
	}
	if cn4OK {
		if cn < 0 || cn4Fit < (1-opts.CnPenalty)*fit {
			cn = 3 + cn4frac
			fit = cn4Fit
		}
	}
	return cn, fit
}
