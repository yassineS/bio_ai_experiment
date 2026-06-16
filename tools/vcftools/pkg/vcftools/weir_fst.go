package vcftools

import (
	"fmt"
	"math"
	"os"
	"strconv"

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
// to the accumulator. It mirrors upstream output_weir_and_cockerham_fst
// (variant_file_output.cpp:3466-3637): all fully-diploid sites contribute
// (multi-allelic and monomorphic included); the Fst is the per-allele sum
// sum_a / sum_all, and per-allele a/b/c terms that are NaN are excluded from
// those sums. site.a / site.b / site.c carry the per-allele sums so the
// weighted Fst and windowed aggregation reuse the exact upstream numerators.
func (acc *weirFstAccumulator) addVariant(v *vcf.Variant) {
	// Diploid sites only.
	if !siteIsDiploid(v) {
		return
	}

	nPops := len(acc.popNames)
	nAlleles := len(siteAlleles(v))
	if nAlleles == 0 {
		return
	}

	// Per-population, per-allele homozygote / heterozygote counts.
	nHom := make([][]int, nPops)
	nHet := make([][]int, nPops)
	for i := range nHom {
		nHom[i] = make([]int, nAlleles)
		nHet[i] = make([]int, nAlleles)
	}
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
		first, second, _ := parseGTAlleles(s.Data["GT"])
		for j := 0; j < nAlleles; j++ {
			if first == j && second == j {
				nHom[popIdx][j]++
			} else if (first == j || second == j) && first != -1 && second != -1 {
				nHet[popIdx][j]++
			}
		}
	}

	r := float64(nPops)
	n := make([]float64, nPops)
	pMat := make([][]float64, nPops)
	pbar := make([]float64, nAlleles)
	hbar := make([]float64, nAlleles)
	var nbar, sumNsqr, nSum float64

	for i := 0; i < nPops; i++ {
		pMat[i] = make([]float64, nAlleles)
		for j := 0; j < nAlleles; j++ {
			n[i] += float64(nHom[i][j]) + 0.5*float64(nHet[i][j])
			pMat[i][j] = float64(nHet[i][j]) + 2.0*float64(nHom[i][j])
			nbar += n[i]
			pbar[j] += pMat[i][j]
			hbar[j] += float64(nHet[i][j])
		}
		for j := 0; j < nAlleles; j++ {
			pMat[i][j] /= 2.0 * n[i]
		}
		sumNsqr += n[i] * n[i]
	}
	for i := 0; i < nPops; i++ {
		nSum += n[i]
	}
	nbar = nSum / r
	for j := 0; j < nAlleles; j++ {
		pbar[j] /= nSum * 2.0
		hbar[j] /= nSum
	}

	ssqr := make([]float64, nAlleles)
	for j := 0; j < nAlleles; j++ {
		for i := 0; i < nPops; i++ {
			d := pMat[i][j] - pbar[j]
			ssqr[j] += n[i] * d * d
		}
		ssqr[j] /= (r - 1.0) * nbar
	}
	nc := (nSum - (sumNsqr / nSum)) / (r - 1.0)

	var sumA, sumAll float64
	for j := 0; j < nAlleles; j++ {
		a := (ssqr[j] - (pbar[j]*(1.0-pbar[j])-(((r-1.0)*ssqr[j])/r)-(hbar[j]/4.0))/(nbar-1.0)) * nbar / nc
		b := (pbar[j]*(1.0-pbar[j]) - (ssqr[j] * (r - 1.0) / r) - hbar[j]*(((2.0*nbar)-1.0)/(4.0*nbar))) * nbar / (nbar - 1.0)
		c := hbar[j] / 2.0
		if !math.IsNaN(a) && !math.IsNaN(b) && !math.IsNaN(c) {
			sumA += a
			sumAll += a + b + c
		}
	}
	fst := sumA / sumAll

	site := weirFstSite{
		chrom: v.Chrom,
		pos:   v.Pos,
		a:     sumA,
		b:     sumAll - sumA, // so a+b+c == sumAll
		c:     0,
		fst:   fst,
	}
	site.defined = !math.IsNaN(fst)
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
		// Upstream writes the per-site Fst straight to a default ostream
		// (precision 6); a NaN value (monomorphic / undefined) prints "-nan".
		fmt.Fprintf(f, "%s\t%d\t%s\n", site.chrom, site.pos, formatCppDefault(site.fst))
		if site.defined {
			fstSum += site.fst
			fstCount++
		}
		// Weighted Fst uses every processed site; upstream accumulates the
		// per-allele a/b/c sums (here a == sumA, a+b+c == sumABC).
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

	// Upstream prints these summary lines via dbl2str(_, 5) — defaultfloat,
	// precision 5 — to the log, not the .weir.fst file.
	fmt.Fprintf(os.Stderr, "Weir and Cockerham mean Fst estimate: %s\n", dbl2strPrec5(meanFst))
	fmt.Fprintf(os.Stderr, "Weir and Cockerham weighted Fst estimate: %s\n", dbl2strPrec5(weightedFst))

	return nil
}

// dbl2strPrec5 mirrors upstream output_log::dbl2str(n, 5): a default-format
// ostream with precision 5.
func dbl2strPrec5(x float64) string {
	switch {
	case math.IsNaN(x):
		return "-nan"
	case math.IsInf(x, 1):
		return "inf"
	case math.IsInf(x, -1):
		return "-inf"
	}
	return strconv.FormatFloat(x, 'g', 5, 64)
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

	// Bins are indexed by 0-based window index per chromosome, matching
	// upstream output_windowed_weir_and_cockerham_fst: only non-NaN sites
	// contribute, and a site at pos lands in bins [first, last) where
	// first = ceil((pos-win)/step) (clamped) and last = ceil(pos/step).
	type winAcc struct {
		sumA   float64 // sum of per-site sum_a
		sumABC float64 // sum of per-site sum_all
		fstSum float64
		count  int
	}

	var chromOrder []string
	bins := make(map[string][]*winAcc)

	for _, site := range acc.sites {
		if !site.defined {
			continue // upstream skips NaN-fst sites entirely
		}
		if _, seen := bins[site.chrom]; !seen {
			chromOrder = append(chromOrder, site.chrom)
		}
		first := int(math.Ceil(float64(site.pos-windowSize) / float64(stepSize)))
		if first < 0 {
			first = 0
		}
		last := int(math.Ceil(float64(site.pos) / float64(stepSize)))
		if last >= len(bins[site.chrom]) {
			grown := make([]*winAcc, last+1)
			copy(grown, bins[site.chrom])
			bins[site.chrom] = grown
		}
		for idx := first; idx < last; idx++ {
			if bins[site.chrom][idx] == nil {
				bins[site.chrom][idx] = &winAcc{}
			}
			w := bins[site.chrom][idx]
			w.sumA += site.a
			w.sumABC += site.a + site.b + site.c
			w.fstSum += site.fst
			w.count++
		}
	}

	for _, chrom := range chromOrder {
		for sIdx, w := range bins[chrom] {
			// Emit only bins with non-zero sum_all and at least one site.
			if w == nil || w.sumABC == 0 || w.count == 0 ||
				math.IsNaN(w.sumA) || math.IsNaN(w.sumABC) {
				continue
			}
			weighted := w.sumA / w.sumABC
			mean := w.fstSum / float64(w.count)
			fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%s\t%s\n",
				chrom, sIdx*stepSize+1, sIdx*stepSize+windowSize, w.count,
				formatCppDefault(weighted), formatCppDefault(mean))
		}
	}

	return nil
}
