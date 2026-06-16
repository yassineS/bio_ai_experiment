// Linkage-disequilibrium calculation for the native `prune` plugin, a faithful
// port of calc_ld / _calc_r2_ld and vcfbuf_ld in vcfbuf.c. It computes, for a
// pair of sites, the three LD measures upstream exposes:
//
//   - r2 (ldIdxR2): the squared genotype-dosage correlation,
//   - LD (ldIdxLD): Lewontin's normalised D' (capped at 1),
//   - RD (ldIdxRD): Ragsdale's unbiased \hat{D}.
//
// The arithmetic is `+ - * /` and sqrt only (all IEEE correctly-rounded), and
// the accumulation order matches upstream exactly so the float results are
// byte-identical after htslib's float32 narrowing on output.
package bcftools

import (
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// LD measure indices, matching VCFBUF_LD_IDX_* in vcfbuf.h. Note HD and RD
// share index 2 upstream (RD is the user-facing name for Ragsdale's D).
const (
	ldIdxR2 = 0
	ldIdxLD = 1
	ldIdxRD = 2
	ldN     = 3
)

// ldResult mirrors vcfbuf_ld_t: the three measure values and, for each, the
// position (1-based) of the partner record that produced the maximum.
type ldResult struct {
	val [ldN]float64
	pos [ldN]int // 1-based POS of the winning partner record; 0 if none
	set [ldN]bool
}

// dosageOf parses sample i's GT into a dosage following _calc_r2_ld's inner
// loop: it sums 1 for each ALT allele and, on a missing allele, either breaks
// (default) or, with rand-missing, draws a random replacement from the site
// allele frequency. It returns the dosage and the number of alleles consumed.
func dosageOf(g genotype, ok bool, randMissing bool, af float64, rng *drand48) (dsg, n int) {
	if !ok {
		return 0, 0
	}
	for _, a := range g.alleles {
		// There is no explicit vector-end in the textual model; a sample's
		// own ploidy is its allele count, so the loop simply ends.
		if a == missingAllele {
			if !randMissing {
				break
			}
			if rng.float64() >= af {
				dsg++
			}
			n++
			continue
		}
		if a > 0 {
			dsg++
		}
		n++
	}
	return dsg, n
}

// pruneEstimateAF ports _estimate_af: the alt-allele frequency over all called
// (non-missing) alleles of a site, used only when rand-missing fills in
// missing genotypes. It stops a sample's alleles at the first missing entry,
// matching the upstream `break`.
func pruneEstimateAF(v *vcf.Variant) float64 {
	nref, nalt := 0, 0
	for i := range v.Samples {
		g, ok := sampleGT(v, i)
		if !ok {
			continue
		}
		for _, a := range g.alleles {
			if a == missingAllele {
				break
			}
			if a > 0 {
				nalt++
			} else {
				nref++
			}
		}
	}
	if nref+nalt == 0 {
		return 0
	}
	return float64(nalt) / float64(nref+nalt)
}

// calcR2LD ports _calc_r2_ld for the pair (arec=buffered, brec=current). It
// returns false when the values cannot be determined (no common data), exactly
// as upstream returns -1. randMissing/aaf/baf drive the missing-genotype draws;
// rng must be the shared drand48 generator so the draw order matches upstream.
func calcR2LD(arec, brec *vcf.Variant, randMissing bool, aaf, baf float64, rng *drand48) (ldResult, bool) {
	var out ldResult

	var nhd [9]float64
	var ab, aa, bb, a, b float64
	nab := 0
	ndiff := 0
	anTot, bnTot := 0, 0

	n := len(arec.Samples)
	if len(brec.Samples) < n {
		n = len(brec.Samples)
	}
	for i := 0; i < n; i++ {
		ag, aok := sampleGT(arec, i)
		bg, bok := sampleGT(brec, i)
		adsg, an := dosageOf(ag, aok, randMissing, aaf, rng)
		bdsg, bn := dosageOf(bg, bok, randMissing, baf, rng)

		if an != 0 && bn != 0 {
			anTot += an
			aa += float64(adsg * adsg)
			a += float64(adsg)

			bnTot += bn
			bb += float64(bdsg * bdsg)
			b += float64(bdsg)

			if adsg != bdsg {
				ndiff++
			}
			ab += float64(adsg * bdsg)
			nab++
		}
		if an == 2 && bn == 2 { // diploid genotypes only
			nhd[bdsg*3+adsg]++
		}
	}
	if nab == 0 {
		return out, false // no data in common for the two sites
	}

	pa := a / float64(anTot)
	pb := b / float64(bnTot)
	var cor float64
	if ndiff == 0 {
		cor = 1
	} else {
		if aa == a*a/float64(nab) || bb == b*b/float64(nab) { // zero variance, add small noise
			aa += 1e-4
			bb += 1e-4
			ab += 1e-4
			a += 1e-2
			b += 1e-2
			nab++
		}
		cor = (ab - a*b/float64(nab)) / math.Sqrt(aa-a*a/float64(nab)) / math.Sqrt(bb-b*b/float64(nab))
	}

	out.val[ldIdxR2] = cor * cor

	// Lewontin's normalisation of D, capped at 1.
	ld := cor * math.Sqrt(pa*(1-pa)*pb*(1-pb))
	var norm float64
	if ld < 0 {
		norm = -pa * pb
		if -(1-pa)*(1-pb) > norm {
			norm = -(1 - pa) * (1 - pb)
		}
	} else {
		norm = pa * (1 - pb)
		if (1-pa)*pb > norm {
			norm = (1 - pa) * pb
		}
	}
	if norm != 0 {
		if math.Abs(norm) > math.Abs(ld) {
			ld = ld / norm
		} else {
			ld = 1
		}
	}
	if ld == 0 {
		ld = math.Abs(ld) // avoid "-0" on output
	}
	out.val[ldIdxLD] = ld

	hd := (nhd[0] + nhd[1]/2. + nhd[3]/2. + nhd[4]/4.) * (nhd[4]/4. + nhd[5]/2. + nhd[7]/2. + nhd[8])
	hd -= (nhd[1]/2. + nhd[2] + nhd[4]/4. + nhd[5]/2.) * (nhd[3]/2. + nhd[4]/4. + nhd[6] + nhd[7]/2.)
	hd /= float64(nab)
	hd /= float64(nab + 1)
	out.val[ldIdxRD] = hd

	return out, true
}

// calcAC computes (reference allele count, alt allele count) for a record,
// mirroring bcf_calc_ac with BCF_UN_INFO|BCF_UN_FMT: it prefers INFO/AN+AC
// (nref = AN - sum(AC)) and falls back to counting genotypes. This matches the
// ac[0]=ref, ac[1..]=alt layout _prune_sites relies on.
func calcAC(v *vcf.Variant) (nref, nalt int, ok bool) {
	an, anOK := v.Info["AN"]
	ac, acOK := v.Info["AC"]
	if anOK && acOK {
		if anv, err := strconv.Atoi(an); err == nil {
			sum := 0
			good := true
			for _, s := range strings.Split(ac, ",") {
				n, err := strconv.Atoi(s)
				if err != nil {
					good = false
					break
				}
				sum += n
			}
			if good {
				return anv - sum, sum, true
			}
		}
	}
	// Fall back to genotypes.
	ref, alt := 0, 0
	any := false
	for i := range v.Samples {
		gt, gok := sampleGT(v, i)
		if !gok {
			continue
		}
		for _, a := range gt.alleles {
			if a == missingAllele {
				continue
			}
			any = true
			if a > 0 {
				alt++
			} else {
				ref++
			}
		}
	}
	if !any {
		return 0, 0, false
	}
	return ref, alt, true
}
