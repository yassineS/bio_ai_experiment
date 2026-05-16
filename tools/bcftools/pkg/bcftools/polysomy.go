// Package bcftools — see doc.go. This file implements `bcftools
// polysomy`, the chromosomal copy-number estimator. Upstream source:
// `reference_code/bcftools/polysomy.c`.
//
// Upstream's full algorithm fits a Gaussian-mixture model to the
// per-chromosome B-allele frequency (BAF) distribution using GSL's
// peak-fitter — CN2 = one symmetric peak at 0.5, CN3 = two symmetric
// peaks at 1/3 and 2/3, CN4 = three peaks at 0.25 / 0.5 / 0.75 — and
// picks the lowest copy-number whose fit beats
// `(1 - cn_penalty) * previous_fit`. That's a substantial port (the
// peakfit module alone is ~600 lines of nonlinear-least-squares plus
// Monte Carlo seeding).
//
// The v1 simplification (matching the "land the CLI surface, defer
// the heavy math" pattern used by `roh` and `call`) collects the
// per-chromosome BAF distribution, computes mean/median BAF, and
// emits a CN call based on how far the median is from 0.5:
//
//   - n_het == 0                              -> CN1 (monosomy / LOH /
//     no informative hets)
//   - |median - 0.5| <= MinBafDev             -> CN2 (diploid)
//   - median < 0.5 - MinBafDev                -> CN3 (trisomy, RA bias)
//   - median > 0.5 + MinBafDev                -> CN3 (trisomy, AA bias)
//
// CnPenalty acts as an additional dampener: when set to its upstream
// default (0.7) we only flip to CN3 if the deviation is at least
// MinBafDev — exactly the v1 semantics above. Smaller CnPenalty
// values cause us to nudge the threshold proportionally smaller, so
// `--cn-penalty 0.3` will flag more chromosomes as CN3.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// PolysomyOptions controls a polysomy run. The default constructor
// values are 0/empty/nil; PolysomyFile / Polysomy apply the upstream
// defaults (CnPenalty=0.7, MinBafDev=0.1) when a field is left at
// zero.
type PolysomyOptions struct {
	// Sample restricts the analysis to a single sample name. Empty
	// means "all samples in the input". Upstream requires -s when
	// the VCF has more than one sample; we keep the rule but apply
	// it in the CLI runner.
	Sample string
	// SamplesFile is the upstream `-S/--samples-file` shortcut; a
	// flat list, one name per line.
	SamplesFile string
	// Regions / Targets / RegionsFile / TargetsFile match the rest
	// of bcftools: post-filter for now.
	Regions     []string
	Targets     []string
	RegionsFile string
	TargetsFile string
	// IncludeExpr / ExcludeExpr are accepted at the CLI but
	// applied only at the record level (no per-sample filtering).
	IncludeExpr string
	ExcludeExpr string
	// CnPenalty matches upstream's `-c/--cn-penalty` (default 0.7).
	// In v1 it scales MinBafDev: a smaller penalty lowers the
	// threshold so more chromosomes get flagged as CN3.
	CnPenalty float64
	// MinBafDev is the minimum |median(BAF) - 0.5| that flips the
	// CN call from CN2 to CN3 (after CnPenalty scaling). Default 0.1
	// matches upstream's `--min-fraction` default.
	MinBafDev float64
	// IncludeNoise (upstream `-n/--include-noise`) keeps
	// chromosomes that would otherwise be reported as `?` (unknown
	// copy number). v1 always emits a row for every chromosome we
	// see, so this flag is accepted but inert.
	IncludeNoise bool
	// ForceCN (upstream hidden `--force-cn`) bypasses the
	// heuristic and tags every chromosome with the requested copy
	// number. Accepted; useful for downstream pipelines that just
	// want the per-chrom BAF table.
	ForceCN int
}

// PolysomyResult is the per-sample-per-chromosome CN call.
type PolysomyResult struct {
	Sample    string  // sample name
	Chrom     string  // chromosome / contig
	NHet      int     // number of heterozygous sites collected
	MeanBAF   float64 // arithmetic mean of BAF for the sample/chrom
	MedianBAF float64 // median of BAF for the sample/chrom
	CN        float64 // copy-number call (1.0, 2.0, 3.0, or -1 for unknown)
}

// PolysomySummary is what the file-level entry points return.
type PolysomySummary struct {
	Results []PolysomyResult
}

// PolysomyFile streams VCF/BCF from path and writes the per-sample
// × per-chrom TSV to out.
func PolysomyFile(path string, out io.Writer, opts PolysomyOptions) (PolysomySummary, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return PolysomySummary{}, fmt.Errorf("bcftools polysomy: open %s: %w", path, err)
	}
	defer in.Close()
	return Polysomy(in, out, opts)
}

// Polysomy is the streaming entry point. It collects BAF values per
// (sample, chromosome) and emits CN calls.
func Polysomy(in io.Reader, out io.Writer, opts PolysomyOptions) (PolysomySummary, error) {
	if opts.CnPenalty == 0 {
		opts.CnPenalty = 0.7
	}
	if opts.MinBafDev == 0 {
		opts.MinBafDev = 0.1
	}

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

	// Region / target post-filter set.
	regionSet := buildPolysomyRegionSet(opts.Regions, opts.Targets, opts.RegionsFile, opts.TargetsFile)

	// Per-(sample, chrom) BAF accumulator. We use parallel slices
	// keyed by (sampleName, chrom) → []BAF to keep allocations
	// predictable.
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

	// Emit results in deterministic order: sample (input order) then
	// chrom (first-seen order in the VCF).
	sum := PolysomySummary{}
	w := bufio.NewWriter(out)
	defer w.Flush()
	if _, err := fmt.Fprintln(w, "# sample\tchrom\tn_het\tmean_baf\tmedian_baf\tcn_call"); err != nil {
		return sum, err
	}
	for _, sn := range samples {
		for _, chrom := range chromOrder {
			vs := bafs[key{sample: sn, chrom: chrom}]
			res := PolysomyResult{
				Sample:    sn,
				Chrom:     chrom,
				NHet:      len(vs),
				MeanBAF:   polysomyMean(vs),
				MedianBAF: polysomyMedian(vs),
				CN:        polysomyCNCall(vs, opts),
			}
			sum.Results = append(sum.Results, res)
			cnText := polysomyCNText(res.CN)
			if _, err := fmt.Fprintf(w, "%s\t%s\t%d\t%.4f\t%.4f\t%s\n",
				res.Sample, res.Chrom, res.NHet, res.MeanBAF, res.MedianBAF, cnText); err != nil {
				return sum, err
			}
		}
	}
	return sum, nil
}

// polysomyKeepRecord applies the region / target post-filter. Empty
// set → keep everything.
func polysomyKeepRecord(v *vcf.Variant, regions map[string]bool) bool {
	if len(regions) == 0 {
		return true
	}
	return regions[v.Chrom]
}

// buildPolysomyRegionSet folds region / target lists into a CHROM-name
// allowlist. We only honour the chrom-name component (post-filter
// approximation; the per-base interval filter is a follow-up — see
// docs/PARITY_ROADMAP.md#bcftools).
func buildPolysomyRegionSet(regions, targets []string, regionsFile, targetsFile string) map[string]bool {
	out := map[string]bool{}
	addRegion := func(r string) {
		// strip "chr:beg-end" → "chr"
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

// readPolysomyBAF returns the heterozygous-site B-allele frequency
// for the sample at idx. We prefer an explicit FORMAT/BAF field when
// the input provides one (mirroring upstream's `bcftools polysomy`
// requirement of FORMAT/BAF); otherwise we synthesise a BAF from
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
	// Only score het sites — homozygous sites carry no copy-number
	// information in this scheme.
	if gt, ok := data["GT"]; ok {
		if !polysomyIsHet(gt) {
			return 0, false
		}
	}
	// Preferred: explicit FORMAT/BAF.
	if raw, ok := data["BAF"]; ok {
		if f, ok := polysomyParseFloat(raw); ok {
			return f, true
		}
	}
	// Fallback: derive from FORMAT/AD = REF,ALT counts.
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

// polysomyIsHet returns true if gt looks heterozygous. We accept
// both `/`- and `|`-separated diploid genotypes; everything else
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

// polysomyMedian is the median (0 for empty input). For an even
// number of elements we return the average of the two middle values.
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

// polysomyCNCall returns the copy-number call for one chromosome's
// BAF set. The rules are documented in the file header; see the
// upstream `polysomy.c:fit_curves` for the full Gaussian-mixture
// version.
func polysomyCNCall(bafs []float64, opts PolysomyOptions) float64 {
	if opts.ForceCN != 0 {
		return float64(opts.ForceCN)
	}
	if len(bafs) == 0 {
		return 1.0 // no hets → CN1 / LOH
	}
	med := polysomyMedian(bafs)
	dev := med - 0.5
	if dev < 0 {
		dev = -dev
	}
	// Scale the threshold by CnPenalty. Upstream's default 0.7
	// matches the un-scaled MinBafDev; smaller penalties shrink the
	// threshold so smaller deviations qualify.
	threshold := opts.MinBafDev * (opts.CnPenalty / 0.7)
	if dev <= threshold {
		return 2.0
	}
	return 3.0
}

// polysomyCNText turns the float CN call into upstream's printable
// form. -1 is the "unknown" sentinel.
func polysomyCNText(cn float64) string {
	if cn < 0 {
		return "?"
	}
	// Upstream prints %.2f; we follow that so downstream parsers
	// can round-trip our output.
	return fmt.Sprintf("%.2f", cn)
}
