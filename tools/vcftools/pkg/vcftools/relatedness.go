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
	if len(v.Alt) != 1 {
		return
	}
	if len(v.Ref) != 1 || len(v.Alt[0]) != 1 {
		return
	}
	counts := make([]int, len(r.samples))
	var sumX, nIndv int
	for i := range r.samples {
		if i >= len(v.Samples) {
			counts[i] = -1
			continue
		}
		gt, ok := v.Samples[i].Data["GT"]
		if !ok {
			counts[i] = -1
			continue
		}
		gc, _, _ := parseGTForLD(gt)
		counts[i] = gc
		if gc < 0 {
			continue
		}
		sumX += gc
		nIndv++
	}
	if nIndv < 2 {
		return
	}
	// ALT-allele frequency p_alt; reference-allele frequency p = 1 - p_alt.
	// The Yang estimator is symmetric in the choice of reference allele
	// (the (x-2p)(x-2p) / (2p(1-p)) form is invariant under p -> 1-p with
	// x -> 2-x), but vcftools uses the ALT count and ALT freq; we follow
	// the same convention.
	p := float64(sumX) / float64(2*nIndv)
	if p <= 0 || p >= 1 {
		return
	}
	denom := 2 * p * (1 - p)
	for i := range r.samples {
		ci := counts[i]
		if ci < 0 {
			continue
		}
		ai := (float64(ci) - 2*p) * (float64(ci) - 2*p) / denom
		// Diagonal: A_jj contribution.
		r.sumA[i][i] += ai
		r.nSites[i][i]++
		for j := i + 1; j < len(r.samples); j++ {
			cj := counts[j]
			if cj < 0 {
				continue
			}
			a := (float64(ci) - 2*p) * (float64(cj) - 2*p) / denom
			r.sumA[i][j] += a
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
	nAaAa   [][]int // ordered N×N: count of SNPs where both i and j are het
	nAAaa   [][]int // ordered N×N: count of SNPs where i and j are both hom but with different homozygous alleles
	nAa     []int   // per-individual count of het SNPs across all biallelic sites where individual was non-missing
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

// addVariant updates the ordered N×N matrices following upstream's
// output_indv_relatedness_Manichaikul (variant_file_output.cpp:4706-4741).
// N_Aa[i] is a per-individual het count. N_AaAa[i][j] increments when
// both i and j are het. N_AAaa[i][j] increments when both are homozygous
// AND their homozygous allele identities differ.
func (r *relatedness2Runner) addVariant(v *vcf.Variant) {
	if r == nil || r.n == 0 {
		return
	}
	if len(v.Alt) != 1 {
		return
	}
	counts := make([]int, r.n)
	for i := range r.samples {
		counts[i] = -1
		if i >= len(v.Samples) {
			continue
		}
		gt, ok := v.Samples[i].Data["GT"]
		if !ok {
			continue
		}
		gc, _, _ := parseGTForLD(gt)
		counts[i] = gc
	}
	for i := 0; i < r.n; i++ {
		if counts[i] < 0 {
			continue
		}
		if counts[i] == 1 {
			r.nAa[i]++
		}
		for j := 0; j < r.n; j++ {
			if counts[j] < 0 {
				continue
			}
			if counts[i] == 1 && counts[j] == 1 {
				r.nAaAa[i][j]++
			}
			if counts[i] != 1 && counts[j] != 1 && counts[i] != counts[j] {
				r.nAAaa[i][j]++
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
	// Upstream emits the full ordered N×N matrix
	// (variant_file_output.cpp:4744-4757) with
	// phi = (N_AaAa - 2*N_AAaa) / (N_Aa[i] + N_Aa[j]).
	for i := 0; i < r.n; i++ {
		for j := 0; j < r.n; j++ {
			denom := float64(r.nAa[i] + r.nAa[j])
			var phi float64
			num := float64(r.nAaAa[i][j]) - 2.0*float64(r.nAAaa[i][j])
			switch {
			case denom == 0 && num == 0:
				phi = math.NaN()
			case denom == 0 && num < 0:
				phi = math.Inf(-1)
			case denom == 0:
				phi = math.Inf(1)
			default:
				phi = num / denom
			}
			if _, err := fmt.Fprintf(w, "%s\t%s\t%d\t%d\t%d\t%d\t%s\n",
				r.samples[i], r.samples[j], r.nAaAa[i][j], r.nAAaa[i][j],
				r.nAa[i], r.nAa[j], formatCppDouble(phi)); err != nil {
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
			var ajk float64
			if r.nSites[i][j] > 0 {
				ajk = r.sumA[i][j] / float64(r.nSites[i][j])
			}
			if _, err := fmt.Fprintf(w, "%s\t%s\t%g\n", r.samples[i], r.samples[j], ajk); err != nil {
				return err
			}
		}
	}
	return nil
}

// lrohRunner implements --LROH, a 1:1 port of vcftools' output_LROH
// (variant_file_output.cpp). It detects Long Runs Of Homozygosity per
// individual using the HMM of Auton et al. (Genome Research, 2009) decoded
// with the forward-backward algorithm.
//
// Two hidden states per site are modelled: autozygous (index 0) and
// non-autozygous (index 1). Each kept, diploid, non-missing genotype emits a
// homozygous/heterozygous observation whose emission probabilities depend on
// the per-site heterozygosity and a fixed genotype error rate. Transition
// probabilities are derived from the genetic distance between consecutive
// observed sites (assuming 1cM/Mb and a fixed number of generations since
// common ancestry). A site is called autozygous when its posterior
// probability exceeds p_auto_threshold.
//
// Output: <prefix>.LROH with the 8-column header
// CHROM AUTO_START AUTO_END MIN_START MAX_END
// N_VARIANTS_BETWEEN_MAX_BOUNDARIES N_MISMATCHES INDV.
type lrohRunner struct {
	samples []string
	// chrom is the CHROM of the most recently processed kept site. Upstream
	// keeps a single CHROM variable and uses its final value for every output
	// row, so we replicate that quirk for byte-for-byte parity.
	chrom string

	// Per-individual observation vectors, appended only at sites where the
	// individual has a usable diploid genotype.
	emissionAuto    [][]float64 // P(emission | autozygous)
	emissionNonAuto [][]float64 // P(emission | non-autozygous)
	isHet           [][]bool    // whether each observation is heterozygous
	posVec          [][]int     // POS of each observation
	transA          [][][4]float64
	lastPOS         []int

	// minVariants is the LROH-min-variants threshold (min_SNPs upstream).
	minVariants int
}

// LROH HMM constants, mirroring output_LROH exactly.
const (
	lrohNGen              = 4    // generations since common ancestry
	lrohGenotypeErrorRate = 0.01 // assumed genotype error rate
	lrohPAutoPrior        = 0.05 // prior probability of the autozygous state
	lrohPAutoThreshold    = 0.99 // posterior threshold for reporting autozygosity
)

const defaultLROHMinVariants = 0

func newLROHRunner(samples []string, minVariants int) *lrohRunner {
	if minVariants < 0 {
		minVariants = defaultLROHMinVariants
	}
	n := len(samples)
	return &lrohRunner{
		samples:         append([]string(nil), samples...),
		emissionAuto:    make([][]float64, n),
		emissionNonAuto: make([][]float64, n),
		isHet:           make([][]bool, n),
		posVec:          make([][]int, n),
		transA:          make([][][4]float64, n),
		lastPOS:         initIntSlice(n, -1),
		minVariants:     minVariants,
	}
}

func initIntSlice(n, v int) []int {
	s := make([]int, n)
	for i := range s {
		s[i] = v
	}
	return s
}

// lrohGenotype classifies a genotype string for LROH. It returns:
//
//	state = -1  missing / non-diploid / unusable (no observation)
//	state =  0  homozygous
//	state =  1  heterozygous
//
// and nonRef indicating whether any allele is non-reference. This mirrors
// upstream's get_indv_GENOTYPE_ids + ploidy check: any integer allele indices
// are accepted (not just 0/1), the genotype must be diploid, and a missing
// allele on either side skips the individual at that site.
func lrohGenotype(gt string) (state int, nonRef bool) {
	if gt == "" {
		return -1, false
	}
	// Find the allele separator ('/' or '|').
	sep := -1
	for i := 0; i < len(gt); i++ {
		if gt[i] == '/' || gt[i] == '|' {
			sep = i
			break
		}
	}
	if sep < 0 {
		// Haploid (or malformed): not diploid, skip.
		return -1, false
	}
	left := gt[:sep]
	right := gt[sep+1:]
	// A second separator would make this non-diploid (e.g. triploid).
	for i := 0; i < len(right); i++ {
		if right[i] == '/' || right[i] == '|' {
			return -1, false
		}
	}
	a, aOK := parseAlleleForLD(left)
	b, bOK := parseAlleleForLD(right)
	if !aOK || !bOK {
		return -1, false
	}
	if a < 0 || b < 0 {
		// Missing genotype on either allele.
		return -1, false
	}
	nonRef = a > 0 || b > 0
	if a != b {
		return 1, nonRef
	}
	return 0, nonRef
}

// addVariant processes one kept site, appending an observation for each
// individual with a usable diploid genotype. Sites at which no individual
// carries a non-reference allele are skipped entirely (upstream's
// has_non_ref == false short-circuit).
func (r *lrohRunner) addVariant(v *vcf.Variant) {
	if r == nil || len(r.samples) == 0 {
		return
	}

	states := make([]int, len(r.samples)) // -1 unusable, 0 hom, 1 het
	nGenotypes := 0
	nHets := 0
	hasNonRef := false
	for i := range r.samples {
		states[i] = -1
		if i >= len(v.Samples) {
			continue
		}
		gt, ok := v.Samples[i].Data["GT"]
		if !ok {
			continue
		}
		st, nonRef := lrohGenotype(gt)
		if st < 0 {
			continue
		}
		if nonRef {
			hasNonRef = true
		}
		nGenotypes++
		if st == 1 {
			nHets++
		}
		states[i] = st
	}

	if !hasNonRef || nGenotypes == 0 {
		return
	}

	r.chrom = v.Chrom
	pos := v.Pos
	h := float64(nHets) / float64(nGenotypes) // site heterozygosity

	for i := range r.samples {
		st := states[i]
		if st < 0 {
			continue
		}
		var pAuto, pNonAuto float64
		if st == 1 { // heterozygote
			pNonAuto = h
			pAuto = lrohGenotypeErrorRate
			r.isHet[i] = append(r.isHet[i], true)
		} else { // homozygote
			pNonAuto = 1.0 - h
			pAuto = 1.0 - lrohGenotypeErrorRate
			r.isHet[i] = append(r.isHet[i], false)
		}
		r.emissionAuto[i] = append(r.emissionAuto[i], pAuto)
		r.emissionNonAuto[i] = append(r.emissionNonAuto[i], pNonAuto)

		rec := 0.0
		if r.lastPOS[i] > 0 {
			// Assume 1cM/Mb. Distance in Morgans.
			rec = float64(pos-r.lastPOS[i]) / 1000000.0 / 100.0
		}
		eVal := 1.0 - math.Exp(-2.0*float64(lrohNGen)*rec)
		pAutoToNon := (1.0 - lrohPAutoPrior) * eVal
		pNonToAuto := lrohPAutoPrior * eVal
		var A [4]float64
		A[0] = 1.0 - pNonToAuto // auto -> auto
		A[1] = pAutoToNon       // auto -> non-auto
		A[2] = pNonToAuto       // non-auto -> auto
		A[3] = 1.0 - pAutoToNon // non-auto -> non-auto
		r.transA[i] = append(r.transA[i], A)
		r.posVec[i] = append(r.posVec[i], pos)
		r.lastPOS[i] = pos
	}
}

// writeOutput runs the forward-backward decode per individual and emits
// <prefix>.LROH.
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
	for i := range r.samples {
		if err := r.decodeIndividual(w, i); err != nil {
			return err
		}
	}
	return nil
}

// decodeIndividual runs forward-backward over individual i's observations and
// writes its autozygous runs.
func (r *lrohRunner) decodeIndividual(w *bufio.Writer, ui int) error {
	emA := r.emissionAuto[ui]
	emN := r.emissionNonAuto[ui]
	trans := r.transA[ui]
	pos := r.posVec[ui]
	isHet := r.isHet[ui]
	nObs := len(emA)
	if nObs == 0 {
		return nil
	}

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

	// Run extraction, mirroring output_LROH's loop exactly.
	inAuto := false
	startPos := 0
	endPos := 0
	nSNPs := 0
	nSNPsBetweenHets := 0
	nHetsInRegion := 0
	lastHetPos := pos[0]
	nextHetPos := -1
	for i := 0; i < nObs; i++ {
		if pAuto[i] > lrohPAutoThreshold {
			if !inAuto {
				startPos = pos[i]
			}
			nSNPs++
			nSNPsBetweenHets++
			if isHet[i] {
				nHetsInRegion++
			}
			inAuto = true
		} else {
			if inAuto {
				nextHetPos = pos[nObs-1]
				for j := i; j < nObs; j++ {
					if isHet[j] {
						nextHetPos = pos[j]
						break
					}
					nSNPsBetweenHets++
				}
				endPos = pos[i-1]
				if nSNPs >= r.minVariants {
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
				lastHetPos = pos[i]
				nSNPsBetweenHets = 0
			}
		}
	}
	if inAuto {
		endPos = pos[nObs-1]
		nextHetPos = pos[nObs-1]
		if nSNPs >= r.minVariants {
			if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",
				r.chrom, startPos, endPos, lastHetPos+1, nextHetPos, nSNPsBetweenHets, nHetsInRegion, r.samples[ui]); err != nil {
				return err
			}
		}
	}
	return nil
}
