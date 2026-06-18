package vcftools

import (
	"bufio"
	"fmt"
	"math"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// relatednessRunner streams variants and accumulates pairwise relatedness
// statistics using the Yang et al. (2010) unadjusted A_jk estimator
// (vcftools' --relatedness). The estimator for individuals j,k at biallelic
// SNP i with reference-allele frequency p_i and diploid ALT counts
// x_ij, x_ik in {0,1,2} is
//
//	A_jk^i = ((x_ij - 2 p_i)(x_ik - 2 p_i)) / (2 p_i (1 - p_i))
//
// and the per-pair statistic is the average of A_jk^i over the SNPs where
// both individuals are non-missing and where 0 < p_i < 1 (we restrict to
// biallelic SNPs).
type relatednessRunner struct {
	samples []string
	// sumA[i][j] is the running numerator (sum of A_jk^i) for pair (samples[i],
	// samples[j]) over SNPs where both individuals contributed a non-missing
	// diploid genotype.
	sumA [][]float64
	// nSites[i][j] counts the SNPs that contributed to the (i,j) pair.
	nSites [][]int
}

func newRelatednessRunner(samples []string) *relatednessRunner {
	n := len(samples)
	r := &relatednessRunner{samples: append([]string(nil), samples...)}
	r.sumA = make([][]float64, n)
	r.nSites = make([][]int, n)
	for i := range r.sumA {
		r.sumA[i] = make([]float64, n)
		r.nSites[i] = make([]int, n)
	}
	return r
}

// addVariant adds a single variant's contribution to all pairs. We restrict
// to biallelic SNPs (1 ALT, REF/ALT both single bases) so the genotype encoding
// 0/1/2 carries the usual meaning. Missing or multi-allelic genotypes are
// skipped per-individual.
func (r *relatednessRunner) addVariant(v *vcf.Variant) {
	if r == nil || len(r.samples) == 0 {
		return
	}
	// Biallelic only (REF + exactly one non-"." ALT) and fully diploid.
	if len(siteAlleles(v)) != 2 {
		return
	}
	if !siteIsDiploid(v) {
		return
	}

	alleleCounts, nChr := siteAlleleCountsIndexed(v)
	if nChr == 0 {
		return
	}
	// ALT-allele frequency, exactly as upstream (allele_counts[1]/N_chr).
	freq := float64(alleleCounts[1]) / float64(nChr)
	const eps = 2.220446049250313e-16
	if freq <= eps || freq >= 1.0-eps {
		return
	}

	// x[i] = sum of the two allele ids (0,1,2) for individual i, or -1 if
	// the genotype is missing (upstream initialises x to -1).
	x := make([]int, len(r.samples))
	for i := range r.samples {
		x[i] = -1
		if i >= len(v.Samples) {
			continue
		}
		first, second, _ := parseGTAlleles(v.Samples[i].Data["GT"])
		if first < 0 || second < 0 {
			continue
		}
		x[i] = first + second
	}

	div := 1.0 / (2.0 * freq * (1.0 - freq))
	for i := range r.samples {
		if x[i] < 0 {
			continue
		}
		xi := float64(x[i])
		// Diagonal term (Yang 2010 eq. 6 self-relatedness numerator).
		r.sumA[i][i] += (xi*xi - (1.0+2.0*freq)*xi + 2.0*freq*freq) * div
		r.nSites[i][i]++
		for j := i + 1; j < len(r.samples); j++ {
			if x[j] < 0 {
				continue
			}
			r.sumA[i][j] += (xi - 2.0*freq) * (float64(x[j]) - 2.0*freq) * div
			r.nSites[i][j]++
		}
	}
}

// relatedness2Runner accumulates the KING-robust kinship coefficient
// (Manichaikul et al. 2010; vcftools' --relatedness2). For each unordered
// pair of individuals (i,j) and biallelic SNPs where both individuals are
// non-missing, count:
//
//	N_Aa_i = # SNPs where i is heterozygous (genotype 1)
//	N_Aa_j = # SNPs where j is heterozygous (genotype 1)
//	N_Aa_Aa = # SNPs where i AND j are both heterozygous
//	N_AA_aa = # SNPs where one is hom-ref (0) and the other is hom-alt (2)
//
// The KING-robust kinship estimator is
//
//	kinship_ij = (2 * N_Aa_Aa - 4 * N_AA_aa - N_Aa_i - N_Aa_j + 2 * N_Aa_min) / (4 * N_Aa_min)
//
// where N_Aa_min = min(N_Aa_i, N_Aa_j). Equivalently the form used by
// vcftools / KING:
//
//	kinship_ij = 0.5 + (2 * N_Aa_Aa - 4 * N_AA_aa - N_Aa_i - N_Aa_j) / (4 * N_Aa_min)
//
// We emit <prefix>.relatedness2 with columns
// INDV1 INDV2 N_AaAa N_AAaa N1_Aa N2_Aa RELATEDNESS_PHI.
//
// References:
//   - Manichaikul et al. (2010) "Robust relationship inference in genome-wide
//     association studies" Bioinformatics 26(22):2867-73.
type relatedness2Runner struct {
	samples []string
	n       int
	nAaAa   [][]int // both heterozygous
	nAAaa   [][]int // both homozygous, for different alleles
	nAa     []int   // per-individual count of heterozygous sites
}

func newRelatedness2Runner(samples []string) *relatedness2Runner {
	n := len(samples)
	r := &relatedness2Runner{
		samples: append([]string(nil), samples...),
		n:       n,
		nAaAa:   make([][]int, n),
		nAAaa:   make([][]int, n),
		nAa:     make([]int, n),
	}
	for i := range samples {
		r.nAaAa[i] = make([]int, n)
		r.nAAaa[i] = make([]int, n)
	}
	return r
}

// addVariant accumulates the Manichaikul (2010) KING-robust tallies over a
// biallelic, fully-diploid site, mirroring upstream
// output_indv_relatedness_Manichaikul (variant_file_output.cpp:4706-4741):
// N_Aa per individual, and the full (ordered) N_AaAa / N_AAaa matrices.
func (r *relatedness2Runner) addVariant(v *vcf.Variant) {
	if r == nil || r.n == 0 {
		return
	}
	if len(siteAlleles(v)) != 2 {
		return
	}
	if !siteIsDiploid(v) {
		return
	}

	// Parse each individual's (first, second) allele ids once.
	type gp struct{ a, b int }
	g := make([]gp, r.n)
	for i := range r.samples {
		g[i] = gp{-1, -1}
		if i >= len(v.Samples) {
			continue
		}
		first, second, _ := parseGTAlleles(v.Samples[i].Data["GT"])
		g[i] = gp{first, second}
	}

	het := func(p gp) bool { return p.a != p.b && p.a >= 0 && p.b >= 0 }
	hom := func(p gp) bool { return p.a == p.b && p.a >= 0 && p.b >= 0 }

	for ui := 0; ui < r.n; ui++ {
		gi := g[ui]
		if het(gi) {
			r.nAa[ui]++
		}
		for uj := 0; uj < r.n; uj++ {
			gj := g[uj]
			if het(gi) && het(gj) {
				r.nAaAa[ui][uj]++
			}
			if hom(gi) && hom(gj) && gi.a != gj.a {
				r.nAAaa[ui][uj]++
			}
		}
	}
}

func (r *relatedness2Runner) writeOutput(prefix string) error {
	if r == nil {
		return nil
	}
	f, err := iohelper.OpenWriter(prefix + ".relatedness2")
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()
	if _, err := w.WriteString("INDV1\tINDV2\tN_AaAa\tN_AAaa\tN1_Aa\tN2_Aa\tRELATEDNESS_PHI\n"); err != nil {
		return err
	}
	// Upstream emits the FULL ordered N×N matrix (every (ui, uj) pair,
	// including ui==uj and both orderings). phi = (N_AaAa - 2*N_AAaa) /
	// (N_Aa[ui] + N_Aa[uj]); a zero denominator yields "-nan".
	for ui := 0; ui < r.n; ui++ {
		for uj := 0; uj < r.n; uj++ {
			phi := (float64(r.nAaAa[ui][uj]) - 2.0*float64(r.nAAaa[ui][uj])) /
				float64(r.nAa[ui]+r.nAa[uj])
			if _, err := fmt.Fprintf(w, "%s\t%s\t%d\t%d\t%d\t%d\t%s\n",
				r.samples[ui], r.samples[uj],
				r.nAaAa[ui][uj], r.nAAaa[ui][uj], r.nAa[ui], r.nAa[uj],
				formatCppDefault(phi)); err != nil {
				return err
			}
		}
	}
	return nil
}

// phasedBlockRunner detects contiguous runs of phased ("a|b") diploid
// genotypes within each sample and emits <prefix>.blocks with one row per
// run that spans at least 2 variants. Columns: CHROM BLOCK_START BLOCK_END
// N_VARIANTS INDV.
type phasedBlockRunner struct {
	samples  []string
	curChrom []string
	curStart []int
	curEnd   []int
	curN     []int
	blocks   []phasedBlock
}

type phasedBlock struct {
	chrom     string
	start     int
	end       int
	n         int
	sampleIdx int
}

func newPhasedBlockRunner(samples []string) *phasedBlockRunner {
	return &phasedBlockRunner{
		samples:  append([]string(nil), samples...),
		curChrom: make([]string, len(samples)),
		curStart: make([]int, len(samples)),
		curEnd:   make([]int, len(samples)),
		curN:     make([]int, len(samples)),
	}
}

// addVariant updates per-sample phased-block state. A sample is "phased" at
// a site when its GT uses '|' as the separator and both alleles parse to
// {0,1} (no missing).
func (r *phasedBlockRunner) addVariant(v *vcf.Variant) {
	if r == nil {
		return
	}
	for i := range r.samples {
		if i >= len(v.Samples) {
			r.closeRun(i)
			continue
		}
		gt, ok := v.Samples[i].Data["GT"]
		if !ok {
			r.closeRun(i)
			continue
		}
		phased := isPhasedDiploid(gt)
		if !phased {
			r.closeRun(i)
			continue
		}
		if r.curChrom[i] != v.Chrom || r.curN[i] == 0 {
			r.closeRun(i)
			r.curChrom[i] = v.Chrom
			r.curStart[i] = v.Pos
			r.curEnd[i] = v.Pos
			r.curN[i] = 1
			continue
		}
		r.curEnd[i] = v.Pos
		r.curN[i]++
	}
}

// isPhasedDiploid returns true for non-missing phased diploid genotypes
// using the '|' separator. Returns false for "0|.", ".|1", missing, or
// unphased ("a/b") genotypes.
func isPhasedDiploid(gt string) bool {
	if gt == "" {
		return false
	}
	if strings.ContainsRune(gt, '.') {
		return false
	}
	bar := strings.IndexByte(gt, '|')
	if bar < 0 {
		return false
	}
	left := gt[:bar]
	right := gt[bar+1:]
	if left == "" || right == "" {
		return false
	}
	return true
}

func (r *phasedBlockRunner) closeRun(i int) {
	if r.curN[i] >= 2 {
		r.blocks = append(r.blocks, phasedBlock{
			chrom:     r.curChrom[i],
			start:     r.curStart[i],
			end:       r.curEnd[i],
			n:         r.curN[i],
			sampleIdx: i,
		})
	}
	r.curChrom[i] = ""
	r.curStart[i] = 0
	r.curEnd[i] = 0
	r.curN[i] = 0
}

func (r *phasedBlockRunner) writeOutput(prefix string) error {
	if r == nil {
		return nil
	}
	for i := range r.samples {
		r.closeRun(i)
	}
	f, err := iohelper.OpenWriter(prefix + ".blocks")
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()
	if _, err := w.WriteString("CHROM\tBLOCK_START\tBLOCK_END\tN_VARIANTS\tINDV\n"); err != nil {
		return err
	}
	for _, b := range r.blocks {
		if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%s\n",
			b.chrom, b.start, b.end, b.n, r.samples[b.sampleIdx]); err != nil {
			return err
		}
	}
	return nil
}

// writeOutput emits <prefix>.relatedness with rows for every pair
// (INDV1, INDV2) including each individual's self-pair (INDV1 == INDV2), the
// latter being the per-individual diagonal of the GRM. RELATEDNESS_AJK is
// sumA / nSites for the pair, or 0 if nSites == 0.
//
// Layout: INDV1\tINDV2\tRELATEDNESS_AJK (one header row, tab-separated).
func (r *relatednessRunner) writeOutput(prefix string) error {
	if r == nil {
		return nil
	}
	f, err := iohelper.OpenWriter(prefix + ".relatedness")
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()
	if _, err := w.WriteString("INDV1\tINDV2\tRELATEDNESS_AJK\n"); err != nil {
		return err
	}
	for i := range r.samples {
		for j := i; j < len(r.samples); j++ {
			ajk := r.sumA[i][j] / float64(r.nSites[i][j])
			if i == j {
				// Upstream: Ajk[ui][ui] = 1.0 + (sum / N_sites).
				ajk = 1.0 + ajk
			}
			if _, err := fmt.Fprintf(w, "%s\t%s\t%s\n",
				r.samples[i], r.samples[j], formatCppDefault(ajk)); err != nil {
				return err
			}
		}
	}
	return nil
}

// lrohRunner implements --LROH: detect Long Runs of Homozygosity by the
// Boyko/Auton (Genome Research 2009) hidden-Markov method, using the
// forward-backward algorithm exactly as vcftools' output_LROH does. Per site
// it accumulates, for each individual, an emission probability pair
// (autozygous, non-autozygous) and a transition matrix derived from the
// inter-site genetic distance (assuming 1cM/Mb); after the scan it runs
// forward-backward per individual and reports each maximal run of sites whose
// posterior P(autozygous) exceeds the 0.99 threshold.
//
// Output: <prefix>.LROH with the 8-column header
// CHROM AUTO_START AUTO_END MIN_START MAX_END N_VARIANTS_BETWEEN_MAX_BOUNDARIES
// N_MISMATCHES INDV. vcftools requires a single --chr, so the runner assumes
// all sites it sees are on one chromosome.
type lrohRunner struct {
	samples []string
	nIndv   int
	chrom   string
	// per-individual accumulators (only sites the individual contributes to).
	emAuto  [][]float64    // emission given autozygous
	emNon   [][]float64    // emission given non-autozygous
	trans   [][][4]float64 // transition matrix A = [AtoA, AtoN, NtoA, NtoN]
	sPos    [][]int        // site positions
	isHet   [][]bool       // per-site het flag
	lastPOS []int          // previous contributing position per individual
}

// LROH model constants (vcftools output_LROH).
const (
	lrohNGen           = 4.0  // generations since common ancestry
	lrohGenoErr        = 0.01 // assumed genotype error rate
	lrohPAutoPrior     = 0.05 // prior probability of the autozygous state
	lrohPAutoThreshold = 0.99 // posterior threshold for reporting a region
	lrohMinSNPs        = 0    // minimum sites in a reported region
)

func newLROHRunner(samples []string, _ int) *lrohRunner {
	n := len(samples)
	last := make([]int, n)
	for i := range last {
		last[i] = -1
	}
	return &lrohRunner{
		samples: append([]string(nil), samples...),
		nIndv:   n,
		emAuto:  make([][]float64, n),
		emNon:   make([][]float64, n),
		trans:   make([][][4]float64, n),
		sPos:    make([][]int, n),
		isHet:   make([][]bool, n),
		lastPOS: last,
	}
}

// addVariant accumulates the per-individual HMM emission/transition data for
// one site. Sites where every genotype is homozygous-reference are skipped
// (matching upstream's has_non_ref guard); missing and non-diploid genotypes do
// not contribute.
func (r *lrohRunner) addVariant(v *vcf.Variant) {
	if r == nil || r.nIndv == 0 {
		return
	}
	type geno struct {
		het   bool
		valid bool
	}
	g := make([]geno, r.nIndv)
	nGeno, nHet := 0, 0
	hasNonRef := false
	for i := 0; i < r.nIndv; i++ {
		if i >= len(v.Samples) {
			continue
		}
		a, b, ok := parseLROHGenotype(v.Samples[i].Data["GT"])
		if !ok {
			continue
		}
		if a > 0 || b > 0 {
			hasNonRef = true
		}
		het := a != b
		g[i] = geno{het: het, valid: true}
		nGeno++
		if het {
			nHet++
		}
	}
	if !hasNonRef || nGeno == 0 {
		return
	}
	h := float64(nHet) / float64(nGeno) // site heterozygosity
	r.chrom = v.Chrom
	pos := v.Pos
	for i := 0; i < r.nIndv; i++ {
		if !g[i].valid {
			continue
		}
		var emA, emN float64
		if g[i].het {
			emA, emN = lrohGenoErr, h
		} else {
			emA, emN = 1.0-lrohGenoErr, 1.0-h
		}
		rDist := 0.0
		if r.lastPOS[i] > 0 {
			rDist = float64(pos-r.lastPOS[i]) / 1000000.0 / 100.0 // Morgans, 1cM/Mb
		}
		e := 1.0 - math.Exp(-2.0*lrohNGen*rDist)
		pAtoN := (1.0 - lrohPAutoPrior) * e
		pNtoA := lrohPAutoPrior * e
		r.emAuto[i] = append(r.emAuto[i], emA)
		r.emNon[i] = append(r.emNon[i], emN)
		r.trans[i] = append(r.trans[i], [4]float64{1.0 - pNtoA, pAtoN, pNtoA, 1.0 - pAtoN})
		r.sPos[i] = append(r.sPos[i], pos)
		r.isHet[i] = append(r.isHet[i], g[i].het)
		r.lastPOS[i] = pos
	}
}

// writeOutput runs forward-backward per individual and emits <prefix>.LROH.
func (r *lrohRunner) writeOutput(prefix string) error {
	if r == nil {
		return nil
	}
	f, err := iohelper.OpenWriter(prefix + ".LROH")
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()
	if _, err := w.WriteString("CHROM\tAUTO_START\tAUTO_END\tMIN_START\tMAX_END\tN_VARIANTS_BETWEEN_MAX_BOUNDARIES\tN_MISMATCHES\tINDV\n"); err != nil {
		return err
	}
	for ui := 0; ui < r.nIndv; ui++ {
		if err := r.emitIndividual(w, ui); err != nil {
			return err
		}
	}
	return nil
}

// emitIndividual runs the forward-backward HMM for one individual and writes
// its autozygous regions, mirroring vcftools output_LROH exactly (including the
// p_trans[i-1] indexing and the 1e-20/1e20 underflow renormalisation).
func (r *lrohRunner) emitIndividual(w *bufio.Writer, ui int) error {
	nObs := len(r.emAuto[ui])
	if nObs == 0 {
		return nil
	}
	emA, emN, trans := r.emAuto[ui], r.emNon[ui], r.trans[ui]
	sPos, isHet := r.sPos[ui], r.isHet[ui]

	alpha := make([][2]float64, nObs)
	beta := make([][2]float64, nObs)
	alpha[0][0] = emA[0]
	alpha[0][1] = emN[0]
	for i := 1; i < nObs; i++ {
		alpha[i][0] = alpha[i-1][0]*trans[i-1][0]*emA[i] + alpha[i-1][1]*trans[i-1][2]*emA[i]
		alpha[i][1] = alpha[i-1][1]*trans[i-1][3]*emN[i] + alpha[i-1][0]*trans[i-1][1]*emN[i]
		for alpha[i][0]+alpha[i][1] < 1e-20 {
			alpha[i][0] *= 1e20
			alpha[i][1] *= 1e20
		}
	}
	beta[nObs-1][0] = 1.0
	beta[nObs-1][1] = 1.0
	for i := nObs - 2; i >= 0; i-- {
		beta[i][0] = beta[i+1][0]*trans[i][0]*emA[i] + beta[i+1][1]*trans[i][2]*emA[i]
		beta[i][1] = beta[i+1][1]*trans[i][3]*emN[i] + beta[i+1][0]*trans[i][1]*emN[i]
		for beta[i][0]+beta[i][1] < 1e-20 {
			beta[i][0] *= 1e20
			beta[i][1] *= 1e20
		}
	}

	pAuto := make([]float64, nObs)
	for i := 0; i < nObs; i++ {
		pAuto[i] = alpha[i][0] * beta[i][0] / (alpha[i][0]*beta[i][0] + alpha[i][1]*beta[i][1])
	}

	inAuto := false
	startPos, endPos := 0, 0
	nSNPs, nSNPsBetweenHets, nHetsInRegion := 0, 0, 0
	lastHetPos := sPos[0]
	nextHetPos := -1
	for i := 0; i < nObs; i++ {
		if pAuto[i] > lrohPAutoThreshold {
			if !inAuto {
				startPos = sPos[i]
			}
			nSNPs++
			nSNPsBetweenHets++
			if isHet[i] {
				nHetsInRegion++
			}
			inAuto = true
		} else {
			if inAuto {
				nextHetPos = sPos[nObs-1]
				for j := i; j < nObs; j++ {
					if isHet[j] {
						nextHetPos = sPos[j]
						break
					}
					nSNPsBetweenHets++
				}
				endPos = sPos[i-1]
				if nSNPs >= lrohMinSNPs {
					if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",
						r.chrom, startPos, endPos, lastHetPos+1, nextHetPos-1, nSNPsBetweenHets, nHetsInRegion, r.samples[ui]); err != nil {
						return err
					}
				}
			}
			inAuto = false
			nSNPs = 0
			nHetsInRegion = 0
			if isHet[i] {
				lastHetPos = sPos[i]
				nSNPsBetweenHets = 0
			}
		}
	}
	if inAuto {
		endPos = sPos[nObs-1]
		nextHetPos = sPos[nObs-1]
		if nSNPs >= lrohMinSNPs {
			if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",
				r.chrom, startPos, endPos, lastHetPos+1, nextHetPos, nSNPsBetweenHets, nHetsInRegion, r.samples[ui]); err != nil {
				return err
			}
		}
	}
	return nil
}

// parseLROHGenotype parses a diploid GT into its two allele indices. It returns
// ok=false for missing, haploid, or unparseable genotypes (which do not
// contribute to the LROH HMM).
func parseLROHGenotype(gt string) (a, b int, ok bool) {
	if gt == "" || gt == "." {
		return 0, 0, false
	}
	sep := strings.IndexAny(gt, "/|")
	if sep < 0 {
		return 0, 0, false // haploid
	}
	left, right := gt[:sep], gt[sep+1:]
	// Trim any trailing per-sample fields should not appear here (GT only).
	an, aok := parseAlleleForLD(left)
	bn, bok := parseAlleleForLD(right)
	if !aok || !bok || an < 0 || bn < 0 {
		return 0, 0, false
	}
	return an, bn, true
}
