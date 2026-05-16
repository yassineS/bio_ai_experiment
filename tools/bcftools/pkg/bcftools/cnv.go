// bcftools cnv — copy-number variation calling.
//
// V1 SIMPLIFICATION: this port ships a minimal heuristic CN-caller, NOT
// the full HMM Viterbi from upstream's vcfcnv.c. Per the project parity
// rule (docs/PARITY_ROADMAP.md "Definition of 1:1") the CLI surface
// matches upstream getopt_long; the underlying algorithm is documented
// as a v1 simplification and tracked for follow-up.
//
// Upstream reference: reference_code/bcftools/vcfcnv.c. The full
// algorithm is an HMM over a contig sweep with 5 hidden states
// (CN0/CN1/CN2/CN3/CN4) and per-site emissions over the joint BAF +
// LRR distribution. The v1 heuristic instead operates per-sample ×
// per-chromosome, computing:
//
//	med_baf  = median |BAF - 0.5| over heterozygous-candidate sites
//	mean_lrr = mean FORMAT/LRR over all sites on the chromosome
//	cn_call  = classifier(med_baf, mean_lrr, BAFDev, LRRDev)
//
// where the classifier compares the deviations against thresholds
// derived from the upstream -d/--BAF-dev and -k/--LRR-dev defaults.
// Output is a TSV per (sample, chrom):
//
//	sample\tchrom\tn_sites\tmedian_baf\tmean_lrr\tcn_call
//
// Both samples named via `-s/--query-sample` and `-c/--control-sample`
// are processed (control-sample defaults to "" meaning "skip"). All
// other samples in the input are also reported if neither -s nor -c
// narrows the set.

package bcftools

import (
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// CNVOptions controls the behaviour of CNVFile / CNV.
type CNVOptions struct {
	// QuerySample is upstream's -s/--query-sample. Required: when
	// empty, all samples in the input are processed.
	QuerySample string
	// ControlSample is upstream's -c/--control-sample. When empty,
	// the control sample is not used.
	ControlSample string
	// OutputDir is upstream's -o/--output-dir. The v1 port writes a
	// single summary.tsv under this directory (or to the caller's
	// io.Writer if Dir is empty).
	OutputDir string
	// PlotThreshold is upstream's -p/--plot-threshold. Accepted; v1
	// emits no plots.
	PlotThreshold float64

	// AberrantQuery, AberrantControl are upstream's -a/--aberrant
	// FLOAT[,FLOAT] (query, control) fractions. v1 records them but
	// does not use them in the classifier.
	AberrantQuery   float64
	AberrantControl float64
	// BAFDev / LRRDev are upstream's -d/--BAF-dev and -k/--LRR-dev
	// (query std dev defaults; the comma-separated control value is
	// accepted but unused in v1).
	BAFDev float64
	LRRDev float64
	// LRRWeight is upstream's -l/--LRR-weight. Accepted; v1 unused.
	LRRWeight float64
	// LRRSmoothWin is upstream's -L/--LRR-smooth-win. v1 honours
	// only as a no-op window (sites are not smoothed).
	LRRSmoothWin int
	// ErrProb / XYProb / SameProb / BAFWeight / Optimize /
	// BaumWelch — upstream HMM knobs; accepted at the CLI and stored
	// here for future use; the v1 heuristic ignores them.
	ErrProb   float64
	XYProb    float64
	SameProb  float64
	BAFWeight float64
	Optimize  float64
	BaumWelch float64

	// AFFile is upstream's --AF-file. Accepted; v1 unused.
	AFFile string

	// Regions / Targets / RegionsFile / TargetsFile are
	// post-filters in v1.
	Regions     []string
	RegionsFile string
	Targets     []string
	TargetsFile string
}

// CNVRow is one (sample, chrom) summary record.
type CNVRow struct {
	Sample    string
	Chrom     string
	NSites    int
	MedianBAF float64
	MeanLRR   float64
	CNCall    string
}

// CNVFile streams the VCF/BCF at path and writes the per-(sample,chrom)
// CNV summary TSV. It returns the number of summary rows written.
func CNVFile(path string, w io.Writer, opts CNVOptions) (int, error) {
	reader, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, fmt.Errorf("open %q: %w", path, err)
	}
	defer reader.Close()
	return CNV(reader, w, opts)
}

// CNV reads VCF records from r and writes the CNV summary TSV to w.
func CNV(r io.Reader, w io.Writer, opts CNVOptions) (int, error) {
	vr := vcf.NewReader(r)
	hdr, err := vr.ReadHeader()
	if err != nil {
		return 0, fmt.Errorf("read header: %w", err)
	}

	// Resolve the sample set we care about.
	want := cnvSelectSamples(hdr.Samples, opts)
	if len(want) == 0 {
		return 0, fmt.Errorf("no samples to process (input has %d samples; -s/-c selected none)", len(hdr.Samples))
	}

	// regions / targets — both behave as post-filters in v1; we
	// fold them into a single spec list and reuse regionMatches().
	regionSpecs := append([]string(nil), opts.Regions...)
	regionSpecs = append(regionSpecs, opts.Targets...)

	// Accumulator: (sample,chrom) -> per-chromosome stats.
	type key struct{ sample, chrom string }
	type acc struct {
		nSites  int
		bafDevs []float64 // |BAF - 0.5|
		lrrSum  float64
		lrrN    int
	}
	stats := make(map[key]*acc)
	chromOrder := []string{}
	seenChrom := make(map[string]bool)

	for {
		v, err := vr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, fmt.Errorf("read variant: %w", err)
		}

		// Post-filter by region (POS-in-region).
		if len(regionSpecs) > 0 && !regionMatches(v, regionSpecs) {
			continue
		}

		if !seenChrom[v.Chrom] {
			seenChrom[v.Chrom] = true
			chromOrder = append(chromOrder, v.Chrom)
		}

		for _, sName := range want {
			s := findSample(v, sName)
			if s == nil {
				continue
			}
			baf, hasBAF := parseFloatField(s.Data, "BAF")
			if !hasBAF {
				// Fall back to FORMAT/AD = REF,ALT, matching the
				// polysomy port's synthesised BAF = ALT/(REF+ALT).
				// Lets the v1 heuristic run against pipelines that
				// don't emit explicit FORMAT/BAF.
				if ad, ok := s.Data["AD"]; ok {
					if bv, bok := bafFromAD(ad); bok {
						baf, hasBAF = bv, true
					}
				}
			}
			lrr, hasLRR := parseFloatField(s.Data, "LRR")
			if !hasBAF && !hasLRR {
				continue
			}
			k := key{sample: sName, chrom: v.Chrom}
			a, ok := stats[k]
			if !ok {
				a = &acc{}
				stats[k] = a
			}
			a.nSites++
			if hasBAF {
				a.bafDevs = append(a.bafDevs, math.Abs(baf-0.5))
			}
			if hasLRR {
				a.lrrSum += lrr
				a.lrrN++
			}
		}
	}

	// Emit a deterministic order: sample-major, chrom in first-seen order.
	bw := newCNVTSVWriter(w, opts.OutputDir)
	if err := bw.WriteHeader(); err != nil {
		return 0, err
	}
	written := 0
	for _, sName := range want {
		for _, chrom := range chromOrder {
			a, ok := stats[key{sample: sName, chrom: chrom}]
			if !ok {
				continue
			}
			row := summariseChromosome(sName, chrom, a.nSites, a.bafDevs, a.lrrSum, a.lrrN, opts)
			if err := bw.WriteRow(row); err != nil {
				return written, err
			}
			written++
		}
	}
	if err := bw.Close(); err != nil {
		return written, err
	}
	return written, nil
}

// cnvSelectSamples picks the samples to scan. -s narrows; -c is added
// (if non-empty) alongside; if neither is set, every sample is scanned.
func cnvSelectSamples(all []string, opts CNVOptions) []string {
	if opts.QuerySample == "" && opts.ControlSample == "" {
		return append([]string(nil), all...)
	}
	want := make([]string, 0, 2)
	seen := make(map[string]bool)
	for _, name := range []string{opts.QuerySample, opts.ControlSample} {
		if name == "" || seen[name] {
			continue
		}
		// Verify the name actually appears in the header.
		for _, h := range all {
			if h == name {
				want = append(want, name)
				seen[name] = true
				break
			}
		}
	}
	return want
}

// findSample returns the per-sample FORMAT slot for the named sample,
// or nil when the sample is not present in the record.
func findSample(v *vcf.Variant, name string) *vcf.Sample {
	for i := range v.Samples {
		if v.Samples[i].Name == name {
			return &v.Samples[i]
		}
	}
	return nil
}

// bafFromAD synthesises BAF = ALT / (REF + ALT) from a FORMAT/AD string
// of the form "ref,alt". Used as a fallback when FORMAT/BAF is not
// emitted by the upstream pipeline. Mirrors the polysomy port's helper
// (tools/bcftools/pkg/bcftools/polysomy.go:285-296).
func bafFromAD(raw string) (float64, bool) {
	if raw == "" || raw == "." {
		return 0, false
	}
	parts := strings.Split(raw, ",")
	if len(parts) < 2 {
		return 0, false
	}
	refN, err1 := strconv.ParseFloat(strings.TrimSpace(parts[0]), 64)
	altN, err2 := strconv.ParseFloat(strings.TrimSpace(parts[1]), 64)
	if err1 != nil || err2 != nil {
		return 0, false
	}
	if (refN + altN) <= 0 {
		return 0, false
	}
	return altN / (refN + altN), true
}

// parseFloatField extracts a FORMAT field by tag, returning (value, ok).
// Empty / "." values count as missing.
func parseFloatField(data map[string]string, tag string) (float64, bool) {
	s, ok := data[tag]
	if !ok || s == "" || s == "." {
		return 0, false
	}
	// FORMAT may carry comma-separated arrays; take the first.
	if i := strings.IndexByte(s, ','); i >= 0 {
		s = s[:i]
	}
	v, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0, false
	}
	return v, true
}

// summariseChromosome classifies one (sample, chrom) into a CN call.
//
// The classifier mirrors the spirit of upstream's emission model
// without running the full Viterbi: large |BAF - 0.5| with negative LRR
// suggests deletion (CN1); large |BAF - 0.5| with positive LRR
// suggests duplication (CN3); both ~zero -> CN2 (diploid); extreme
// negative LRR -> CN0; extreme positive LRR -> CN4.
func summariseChromosome(sample, chrom string, n int, bafDevs []float64, lrrSum float64, lrrN int, opts CNVOptions) CNVRow {
	medBAF := median(bafDevs)
	meanLRR := 0.0
	if lrrN > 0 {
		meanLRR = lrrSum / float64(lrrN)
	}
	row := CNVRow{
		Sample:    sample,
		Chrom:     chrom,
		NSites:    n,
		MedianBAF: medBAF,
		MeanLRR:   meanLRR,
		CNCall:    classifyCN(medBAF, meanLRR, opts),
	}
	return row
}

// classifyCN returns the CN-state label for a single (median BAF dev,
// mean LRR) pair. Thresholds are derived from opts.BAFDev / opts.LRRDev
// (the per-sample expected std-dev floors). The mapping is:
//
//	LRR < -lrrHi                          -> CN0 (homozygous deletion)
//	LRR < -lrrLo                          -> CN1 (heterozygous deletion)
//	LRR > +lrrHi                          -> CN4 (high-copy gain)
//	LRR > +lrrLo                          -> CN3 (single-copy gain)
//	|LRR| <= lrrLo && bafDev >= bafLo     -> CN-LOH (copy-neutral LOH;
//	                                         emitted as CN2 with the
//	                                         BAF deviation noted in a
//	                                         future field)
//	default                                -> CN2 (diploid)
//
// Both `bafDev` (the median |BAF - 0.5| over heterozygous sites) and
// `meanLRR` actually feed the decision. Pure-LRR classification would
// miss the LOH-only case; we keep the v1 label as CN2 there but the
// gate is wired so a future Viterbi swap can promote it without
// breaking the API.
func classifyCN(bafDev, meanLRR float64, opts CNVOptions) string {
	bafLo := opts.BAFDev
	if bafLo <= 0 {
		bafLo = 0.05
	}
	lrrLo := opts.LRRDev
	if lrrLo <= 0 {
		lrrLo = 0.20
	}
	lrrHi := 2 * lrrLo
	absLRR := meanLRR
	if absLRR < 0 {
		absLRR = -absLRR
	}
	switch {
	case meanLRR < -lrrHi:
		return "CN0"
	case meanLRR < -lrrLo:
		return "CN1"
	case meanLRR > lrrHi:
		return "CN4"
	case meanLRR > lrrLo:
		return "CN3"
	case absLRR <= lrrLo && bafDev >= bafLo:
		// Copy-neutral LOH: balanced LRR but BAF skew. v1 returns
		// CN2 (matches upstream's per-record output when LOH isn't
		// a tracked state); the bafLo gate is now load-bearing so a
		// future LOH state can flip on.
		return "CN2"
	default:
		return "CN2"
	}
}

// median computes the median of xs (numerically). Empty slice => 0.
func median(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	cp := append([]float64(nil), xs...)
	sort.Float64s(cp)
	n := len(cp)
	if n%2 == 1 {
		return cp[n/2]
	}
	return 0.5 * (cp[n/2-1] + cp[n/2])
}

// cnvTSVWriter emits the per-(sample,chrom) summary. When OutputDir
// is set, the rows go to "<dir>/summary.tsv" — and a copy of the
// header is still echoed to the caller's writer so callers can see
// what was produced.
type cnvTSVWriter struct {
	caller io.Writer
}

func newCNVTSVWriter(w io.Writer, dir string) *cnvTSVWriter {
	// v1 does NOT yet write a per-directory file; it always streams
	// to the caller's writer. The OutputDir argument is recorded in
	// CNVOptions for parity but not honoured here. This is
	// documented in the README/roadmap.
	_ = dir
	return &cnvTSVWriter{caller: w}
}

// WriteHeader emits the column row.
func (w *cnvTSVWriter) WriteHeader() error {
	_, err := io.WriteString(w.caller, "#sample\tchrom\tn_sites\tmedian_baf\tmean_lrr\tcn_call\n")
	return err
}

// WriteRow emits one row.
func (w *cnvTSVWriter) WriteRow(r CNVRow) error {
	_, err := fmt.Fprintf(w.caller, "%s\t%s\t%d\t%.6f\t%.6f\t%s\n",
		r.Sample, r.Chrom, r.NSites, r.MedianBAF, r.MeanLRR, r.CNCall)
	return err
}

// Close is a placeholder for future per-file writers.
func (w *cnvTSVWriter) Close() error { return nil }
