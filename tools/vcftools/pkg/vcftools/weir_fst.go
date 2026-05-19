package vcftools

import (
	"fmt"
	"math"
	"os"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// popGenotypeCounts summarises the diploid genotypes of one population at a
// single biallelic SNP: number of individuals with non-missing GT, count of
// alternate alleles across those individuals, and count of heterozygous
// genotypes.
type popGenotypeCounts struct {
	n        int // number of individuals with non-missing diploid genotype
	altCount int // count of "1" alleles across the 2*n chromosomes
	hetCount int // count of heterozygous genotypes (one ref, one alt)
}

// weirFstSite holds the per-site Weir & Cockerham 1984 Fst components for a
// single SNP. The (a, b, c) triplet is the per-site decomposition; fst = a /
// (a + b + c). When the site is undefined (insufficient data, or denominator
// zero) defined is false and fst is NaN.
type weirFstSite struct {
	chrom   string
	pos     int
	a       float64
	b       float64
	c       float64
	fst     float64
	defined bool
}

// loadPopulationFiles reads the per-population sample lists for --weir-fst-pop.
// It returns one slice of sample names per file (in the same order as the input
// files) and an error if any file cannot be read, fewer than two populations
// are supplied, or a sample appears in more than one population.
func loadPopulationFiles(paths []string) ([][]string, error) {
	if len(paths) < 2 {
		return nil, fmt.Errorf("--weir-fst-pop requires at least 2 population files, got %d", len(paths))
	}
	pops := make([][]string, len(paths))
	seen := make(map[string]int) // sample -> index of first pop containing it
	for i, p := range paths {
		samples, err := loadSampleFile(p)
		if err != nil {
			return nil, fmt.Errorf("reading population file %s: %w", p, err)
		}
		// De-duplicate within a single file but preserve order.
		uniq := make([]string, 0, len(samples))
		inThis := make(map[string]bool)
		for _, s := range samples {
			if inThis[s] {
				continue
			}
			inThis[s] = true
			if prev, ok := seen[s]; ok {
				return nil, fmt.Errorf("sample %q appears in multiple population files (%s and %s)", s, paths[prev], p)
			}
			seen[s] = i
			uniq = append(uniq, s)
		}
		pops[i] = uniq
	}
	return pops, nil
}

// weirCockerhamFst implements the per-site Weir & Cockerham 1984 Fst estimator
// (eqs. 2-4 of Weir, B.S. & Cockerham, C.C., "Estimating F-Statistics for the
// Analysis of Population Structure", Evolution 38(6):1358-1370). It takes the
// per-population genotype summaries at a biallelic SNP and returns the (a, b,
// c) variance components together with Fst = a / (a + b + c).
//
// ok is false when n_bar <= 1 (any population has fewer than 2 individuals) or
// when (a + b + c) == 0; callers should treat the SNP as undefined in that
// case (emit nan in per-site output and skip in window aggregation).
func weirCockerhamFst(pops []popGenotypeCounts) (a, b, c, fst float64, ok bool) {
	r := len(pops)
	if r < 2 {
		return 0, 0, 0, 0, false
	}
	for _, p := range pops {
		if p.n < 2 {
			return 0, 0, 0, 0, false
		}
	}

	rf := float64(r)
	var sumN, sumNSq float64
	for _, p := range pops {
		nf := float64(p.n)
		sumN += nf
		sumNSq += nf * nf
	}
	nBar := sumN / rf

	// Allele frequency p_i and observed heterozygosity h_i in each pop.
	pI := make([]float64, r)
	hI := make([]float64, r)
	for i, p := range pops {
		pI[i] = float64(p.altCount) / (2.0 * float64(p.n))
		hI[i] = float64(p.hetCount) / float64(p.n)
	}

	// Weighted means.
	var pBar, hBar float64
	for i, p := range pops {
		nf := float64(p.n)
		pBar += nf * pI[i]
		hBar += nf * hI[i]
	}
	pBar /= sumN
	hBar /= sumN

	// Sample variance of allele frequencies.
	var sSq float64
	for i, p := range pops {
		nf := float64(p.n)
		d := pI[i] - pBar
		sSq += nf * d * d
	}
	sSq /= (rf - 1.0) * nBar

	nC := (sumN - sumNSq/sumN) / (rf - 1.0)

	if nBar <= 1.0 || nC == 0 {
		return 0, 0, 0, 0, false
	}

	pq := pBar * (1.0 - pBar)
	rTerm := ((rf - 1.0) / rf) * sSq

	a = (nBar / nC) * (sSq - (1.0/(nBar-1.0))*(pq-rTerm-hBar/4.0))
	b = (nBar / (nBar - 1.0)) * (pq - rTerm - ((2.0*nBar-1.0)/(4.0*nBar))*hBar)
	c = hBar / 2.0

	denom := a + b + c
	if denom == 0 {
		return a, b, c, math.NaN(), false
	}
	return a, b, c, a / denom, true
}

// weirFstAccumulator groups the per-site Weir & Cockerham Fst data collected
// while streaming a VCF and the population assignment derived from the
// --weir-fst-pop files.
type weirFstAccumulator struct {
	// popNames[i] is the set of sample names assigned to population i.
	popNames []map[string]bool
	// sites holds one entry per processed SNP, in input order.
	sites []weirFstSite
}

// newWeirFstAccumulator builds an accumulator from the loaded population sample
// lists.
func newWeirFstAccumulator(pops [][]string) *weirFstAccumulator {
	sets := make([]map[string]bool, len(pops))
	for i, samples := range pops {
		sets[i] = make(map[string]bool, len(samples))
		for _, s := range samples {
			sets[i][s] = true
		}
	}
	return &weirFstAccumulator{popNames: sets}
}

// addVariant computes the per-site Weir & Cockerham Fst for v and appends it
// to the accumulator. SNPs that are not biallelic single-base substitutions
// are silently skipped (matching upstream vcftools, which only computes Fst
// for biallelic SNPs).
func (acc *weirFstAccumulator) addVariant(v *vcf.Variant) {
	// Biallelic single-base SNPs only.
	if len(v.Alt) != 1 {
		return
	}
	if isIndelVariant(v) {
		return
	}
	if len(v.Ref) != 1 || len(v.Alt[0]) != 1 {
		return
	}

	pops := make([]popGenotypeCounts, len(acc.popNames))
	for _, s := range v.Samples {
		popIdx := -1
		for i, set := range acc.popNames {
			if set[s.Name] {
				popIdx = i
				break
			}
		}
		if popIdx < 0 {
			continue
		}
		gt, ok := s.Data["GT"]
		if !ok {
			continue
		}
		alleles := strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' })
		if len(alleles) != 2 {
			continue
		}
		if alleles[0] == "" || alleles[0] == "." || alleles[1] == "" || alleles[1] == "." {
			continue
		}
		// Only the reference (0) and the single alt (1) alleles count.
		if (alleles[0] != "0" && alleles[0] != "1") || (alleles[1] != "0" && alleles[1] != "1") {
			continue
		}
		pc := &pops[popIdx]
		pc.n++
		if alleles[0] == "1" {
			pc.altCount++
		}
		if alleles[1] == "1" {
			pc.altCount++
		}
		if alleles[0] != alleles[1] {
			pc.hetCount++
		}
	}

	// Skip the SNP entirely if any population has < 2 individuals.
	for _, p := range pops {
		if p.n < 2 {
			return
		}
	}

	a, b, c, fst, ok := weirCockerhamFst(pops)
	site := weirFstSite{
		chrom: v.Chrom,
		pos:   v.Pos,
		a:     a,
		b:     b,
		c:     c,
	}
	if ok {
		site.fst = fst
		site.defined = true
	} else {
		site.fst = math.NaN()
		site.defined = false
	}
	acc.sites = append(acc.sites, site)
}

// outputWeirFst writes the per-site Weir & Cockerham Fst table to
// <prefix>.weir.fst and prints the mean / weighted summary lines to stderr,
// matching upstream vcftools' summary output.
func (acc *weirFstAccumulator) outputWeirFst(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".weir.fst")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tPOS\tWEIR_AND_COCKERHAM_FST")

	var fstSum float64
	var fstCount int
	var sumA, sumABC float64

	for _, site := range acc.sites {
		if site.defined {
			fmt.Fprintf(f, "%s\t%d\t%.6f\n", site.chrom, site.pos, site.fst)
			fstSum += site.fst
			fstCount++
		} else {
			fmt.Fprintf(f, "%s\t%d\tnan\n", site.chrom, site.pos)
		}
		// Weighted Fst uses every processed SNP regardless of whether the
		// per-site ratio is defined; an undefined ratio just means a + b + c
		// is zero and contributes nothing to either sum.
		sumA += site.a
		sumABC += site.a + site.b + site.c
	}

	meanFst := math.NaN()
	if fstCount > 0 {
		meanFst = fstSum / float64(fstCount)
	}
	weightedFst := math.NaN()
	if sumABC != 0 {
		weightedFst = sumA / sumABC
	}

	fmt.Fprintf(os.Stderr, "Weir and Cockerham mean Fst estimate: %s\n", formatFloat(meanFst))
	fmt.Fprintf(os.Stderr, "Weir and Cockerham weighted Fst estimate: %s\n", formatFloat(weightedFst))

	return nil
}

// outputWindowedWeirFst writes the windowed Weir & Cockerham Fst table to
// <prefix>.windowed.weir.fst. Window starts are 1, 1+step, 1+2*step, ...; a
// SNP at 1-based position p belongs to every window [ws, ws+windowSize-1]
// containing it (the same scheme as outputWindowedPi). If stepSize is zero or
// larger than windowSize the windows are non-overlapping.
func (acc *weirFstAccumulator) outputWindowedWeirFst(prefix string, windowSize, stepSize int) error {
	if windowSize <= 0 {
		return nil
	}
	if stepSize <= 0 || stepSize > windowSize {
		stepSize = windowSize
	}

	f, err := iohelper.OpenWriter(prefix + ".windowed.weir.fst")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST")

	type winAcc struct {
		winStart  int
		nVariants int
		sumA      float64
		sumABC    float64
		fstSum    float64
		fstCount  int
	}

	var chromOrder []string
	windows := make(map[string]map[int]*winAcc)

	for _, site := range acc.sites {
		if _, seen := windows[site.chrom]; !seen {
			windows[site.chrom] = make(map[int]*winAcc)
			chromOrder = append(chromOrder, site.chrom)
		}
		p := site.pos
		kMax := (p - 1) / stepSize
		kMin := 0
		if p > windowSize {
			kMin = (p - windowSize + stepSize - 1) / stepSize
		}
		for k := kMin; k <= kMax; k++ {
			ws := 1 + k*stepSize
			acc := windows[site.chrom][ws]
			if acc == nil {
				acc = &winAcc{winStart: ws}
				windows[site.chrom][ws] = acc
			}
			acc.nVariants++
			acc.sumA += site.a
			acc.sumABC += site.a + site.b + site.c
			if site.defined {
				acc.fstSum += site.fst
				acc.fstCount++
			}
		}
	}

	for _, chrom := range chromOrder {
		var starts []int
		for ws := range windows[chrom] {
			starts = append(starts, ws)
		}
		sort.Ints(starts)
		for _, ws := range starts {
			w := windows[chrom][ws]
			weighted := math.NaN()
			if w.sumABC != 0 {
				weighted = w.sumA / w.sumABC
			}
			mean := math.NaN()
			if w.fstCount > 0 {
				mean = w.fstSum / float64(w.fstCount)
			}
			fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%s\t%s\n",
				chrom, ws, ws+windowSize-1, w.nVariants,
				formatFloat6(weighted), formatFloat6(mean))
		}
	}

	return nil
}

// formatFloat formats a float for the stderr summary lines, emitting "nan"
// when the value is NaN to match the per-site output.
func formatFloat(x float64) string {
	if math.IsNaN(x) {
		return "nan"
	}
	return fmt.Sprintf("%f", x)
}

// formatFloat6 is the %.6f counterpart of formatFloat used in the per-site /
// windowed table columns.
func formatFloat6(x float64) string {
	if math.IsNaN(x) {
		return "nan"
	}
	return fmt.Sprintf("%.6f", x)
}
