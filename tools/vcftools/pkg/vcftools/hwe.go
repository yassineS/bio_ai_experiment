package vcftools

import (
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// snpHWE implements Wigginton/Cao/Abecasis (2005) exact test of Hardy-Weinberg
// equilibrium for a single biallelic SNP. It mirrors upstream vcftools'
// `entry::SNPHWE` (reference_code/vcftools/src/cpp/entry.cpp:18) line-for-line
// so that `--hwe` produces byte-identical filtering decisions.
//
// Inputs are the observed counts of heterozygotes, homozygote-1 (REF/REF) and
// homozygote-2 (ALT/ALT) at the site. The function returns the two-sided
// exact p-value `pHWE` (sum of probabilities of every possible heterozygote
// count whose probability is <= that of the observed count). When the site
// has no called genotypes pHWE defaults to 1.0 (vacuously in equilibrium).
//
// Reference:
//
//	Wigginton, Cao & Abecasis (2005) "A note on exact tests of Hardy-Weinberg
//	equilibrium" Am. J. Hum. Genet. 76:887-883.
func snpHWE(obsHets, obsHom1, obsHom2 int) float64 {
	if obsHom1+obsHom2+obsHets == 0 {
		return 1.0
	}
	if obsHom1 < 0 || obsHom2 < 0 || obsHets < 0 {
		return 1.0
	}

	// homc = the more common homozygote count, homr = the rarer.
	obsHomC := obsHom1
	obsHomR := obsHom2
	if obsHom1 < obsHom2 {
		obsHomC = obsHom2
		obsHomR = obsHom1
	}

	rareCopies := 2*obsHomR + obsHets
	genotypes := obsHets + obsHomC + obsHomR
	if genotypes == 0 || rareCopies == 0 {
		return 1.0
	}

	hetProbs := make([]float64, rareCopies+1)

	// Start at midpoint (most-probable heterozygote count). Upstream uses
	// integer arithmetic: `mid = rare_copies * (2 * genotypes - rare_copies)
	// / (2 * genotypes)`.
	mid := rareCopies * (2*genotypes - rareCopies) / (2 * genotypes)

	// Ensure midpoint and rare alleles have same parity.
	if (rareCopies & 1) != (mid & 1) {
		mid++
	}
	if mid > rareCopies {
		mid = rareCopies
	}

	currHets := mid
	currHomR := (rareCopies - mid) / 2
	currHomC := genotypes - currHets - currHomR

	hetProbs[mid] = 1.0
	sum := hetProbs[mid]

	// Walk downward by 2 from the midpoint.
	for currHets = mid; currHets > 1; currHets -= 2 {
		hetProbs[currHets-2] = hetProbs[currHets] *
			float64(currHets) * float64(currHets-1) /
			(4.0 * float64(currHomR+1) * float64(currHomC+1))
		sum += hetProbs[currHets-2]
		// 2 fewer heterozygotes -> add one rare, one common homozygote.
		currHomR++
		currHomC++
	}

	currHets = mid
	currHomR = (rareCopies - mid) / 2
	currHomC = genotypes - currHets - currHomR

	// Walk upward by 2 from the midpoint.
	for currHets = mid; currHets <= rareCopies-2; currHets += 2 {
		hetProbs[currHets+2] = hetProbs[currHets] *
			4.0 * float64(currHomR) * float64(currHomC) /
			(float64(currHets+2) * float64(currHets+1))
		sum += hetProbs[currHets+2]
		// 2 more heterozygotes -> subtract one rare, one common homozygote.
		currHomR--
		currHomC--
	}

	if sum <= 0 {
		return 1.0
	}
	for i := range hetProbs {
		hetProbs[i] /= sum
	}

	// Two-sided p_hwe: sum probabilities of het counts whose probability
	// is no larger than the observed count's probability. Upstream guards
	// against the observed-index being out of range by relying on the
	// caller; we replicate that — `obs_hets` must be in [0, rare_copies]
	// for any sensible call.
	if obsHets < 0 || obsHets > rareCopies {
		return 1.0
	}
	obsP := hetProbs[obsHets]
	pHWE := 0.0
	for _, p := range hetProbs {
		if p > obsP {
			continue
		}
		pHWE += p
	}
	if pHWE > 1.0 {
		pHWE = 1.0
	}
	return pHWE
}

// countDiploidGenotypes counts the diploid hom1/het/hom2 categories at a
// biallelic site over all included samples. Missing genotypes, haploid calls,
// and calls referencing higher-index alleles are skipped (mirroring
// upstream's `entry::get_genotype_counts`, which only considers diploid
// REF/ALT pairs at biallelic sites).
//
// Returns (hom1, het, hom2, ok). `ok` is false if the site is not biallelic.
func countDiploidGenotypes(v *vcf.Variant) (hom1, het, hom2 int, ok bool) {
	if len(v.Alt) != 1 {
		return 0, 0, 0, false
	}
	for _, sample := range v.Samples {
		gt, present := sample.Data["GT"]
		if !present || gt == "" || gt == "." {
			continue
		}
		alleles := strings.FieldsFunc(gt, func(r rune) bool {
			return r == '/' || r == '|'
		})
		if len(alleles) != 2 {
			continue
		}
		if alleles[0] == "." || alleles[1] == "." {
			continue
		}
		// Restrict to {0,1} — out-of-range alleles can't appear at a
		// biallelic site by construction, but the guard mirrors upstream.
		if alleles[0] == "0" && alleles[1] == "0" {
			hom1++
		} else if alleles[0] == "1" && alleles[1] == "1" {
			hom2++
		} else if (alleles[0] == "0" && alleles[1] == "1") ||
			(alleles[0] == "1" && alleles[1] == "0") {
			het++
		}
	}
	return hom1, het, hom2, true
}

// countMissingChromosomes returns the count of missing-haploid alleles (in
// upstream terms, `N_chr - N_non_missing_chr`) summed across all samples at
// the site. A "./." diploid contributes 2, "0/." contributes 1, "0/0"
// contributes 0, "." haploid contributes 1, "0" haploid contributes 0. This
// mirrors `entry::get_allele_counts`'s accounting in
// `entry_filters.cpp:872-919` where each missing allele is counted as one
// chromosome short.
func countMissingChromosomes(v *vcf.Variant) int {
	missing := 0
	for _, sample := range v.Samples {
		gt, present := sample.Data["GT"]
		if !present || gt == "" {
			// No GT recorded at all: upstream treats this as fully
			// missing (both chromosomes), matching how
			// `entry::parse_genotype_entry` reports missing.
			missing += 2
			continue
		}
		// Find a single phasing/unphasing separator. If there is no
		// separator the call is haploid.
		sep := -1
		for i, r := range gt {
			if r == '/' || r == '|' {
				sep = i
				break
			}
		}
		if sep < 0 {
			// Haploid call. Only one chromosome here; missing if ".".
			if gt == "." {
				missing++
			}
			continue
		}
		left := gt[:sep]
		right := gt[sep+1:]
		if left == "." {
			missing++
		}
		if right == "." {
			missing++
		}
	}
	return missing
}
