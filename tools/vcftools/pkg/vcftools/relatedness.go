package vcftools

import (
	"bufio"
	"fmt"
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

// lrohRunner implements --LROH: for each individual, emit contiguous runs of
// homozygous diploid genotypes ("0/0" or "1/1") that span at least
// `lrohMinVariants` consecutive variants on the same chromosome.
//
// Output: <prefix>.LROH with header CHROM AUTO_START AUTO_END N_VARIANTS INDV
// where AUTO_START and AUTO_END are the 1-based VCF positions of the first
// and last homozygous variant in the run.
type lrohRunner struct {
	samples []string
	// running state per sample: current chrom, start pos, end pos, count.
	curChrom []string
	curStart []int
	curEnd   []int
	curN     []int
	// closed runs: chrom, start, end, n, sampleIdx.
	runs []lrohRun
	// minimum number of consecutive homozygous variants for a run to be emitted.
	minVariants int
}

type lrohRun struct {
	chrom     string
	start     int
	end       int
	n         int
	sampleIdx int
}

const defaultLROHMinVariants = 10

func newLROHRunner(samples []string, minVariants int) *lrohRunner {
	if minVariants <= 0 {
		minVariants = defaultLROHMinVariants
	}
	r := &lrohRunner{
		samples:     append([]string(nil), samples...),
		curChrom:    make([]string, len(samples)),
		curStart:    make([]int, len(samples)),
		curEnd:      make([]int, len(samples)),
		curN:        make([]int, len(samples)),
		minVariants: minVariants,
	}
	return r
}

// addVariant updates per-sample homozygous-run state. We treat ANY variant
// (not just biallelic SNPs) as a usable site for LROH; non-homozygous,
// missing, or chromosome-change events close the current run.
func (r *lrohRunner) addVariant(v *vcf.Variant) {
	if r == nil || len(r.samples) == 0 {
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
		gc, _, _ := parseGTForLD(gt)
		// Homozygous when gc == 0 (0/0) or gc == 2 (1/1). gc==1 is het,
		// gc<0 missing/skipped.
		if gc != 0 && gc != 2 {
			r.closeRun(i)
			continue
		}
		if r.curChrom[i] != v.Chrom || r.curN[i] == 0 {
			// Either first variant or chromosome change: start a fresh run.
			r.closeRun(i)
			r.curChrom[i] = v.Chrom
			r.curStart[i] = v.Pos
			r.curEnd[i] = v.Pos
			r.curN[i] = 1
			continue
		}
		// Extend the existing run.
		r.curEnd[i] = v.Pos
		r.curN[i]++
	}
}

// closeRun appends the current open run for sample i to runs (if long enough)
// and resets the per-sample state.
func (r *lrohRunner) closeRun(i int) {
	if r.curN[i] >= r.minVariants {
		r.runs = append(r.runs, lrohRun{
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

// writeOutput closes any open runs and emits <prefix>.LROH.
// Layout: CHROM AUTO_START AUTO_END N_VARIANTS INDV.
func (r *lrohRunner) writeOutput(prefix string) error {
	if r == nil {
		return nil
	}
	for i := range r.samples {
		r.closeRun(i)
	}
	f, err := iohelper.OpenWriter(prefix + ".LROH")
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()
	if _, err := w.WriteString("CHROM\tAUTO_START\tAUTO_END\tN_VARIANTS\tINDV\n"); err != nil {
		return err
	}
	for _, run := range r.runs {
		if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%s\n",
			run.chrom, run.start, run.end, run.n, r.samples[run.sampleIdx]); err != nil {
			return err
		}
	}
	return nil
}
