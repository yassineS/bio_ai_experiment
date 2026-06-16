// Floating-point de-novo mutation models for the trio-dnm3 plugin: DMM
// (Dirichlet-multinomial over FORMAT/AD, the default), ALM (allele-likelihood
// over FORMAT/QS or AD/PL), and DNG (the original DeNovoGear model over
// FORMAT/PL). These port process_trio_DMM / process_trio_ALM / process_trio_DNG
// and their priors + noise helpers from plugins/trio-dnm3.c.
//
// Numerical note (the libm-tolerance boundary): the per-trio score is a
// log/exp/pow/lgamma reduction over the full father x mother x child genotype
// grid. Upstream links the platform libm and htslib's kfunc.c; this port uses
// Go's math package plus the in-tree AS245 kfLgamma / kfBetai transcribed from
// kfunc.c (callc.go). The incomplete-beta and lgamma kernels are bit-stable
// because they avoid libm (kfBetai/kfLgamma are pure +,-,*,/ over a series),
// but math.Log/Exp/Pow are only guaranteed to the last ULP, so after a long sum
// and the final narrowing to a 32-bit float printed with %g, the DNM/VA scores
// can differ from upstream in the last printed digit. They are therefore
// validated with the tolerance-aware proximity helper (numeric_parity_test.go),
// not byte-for-byte; the NAIVE model stays byte-exact. See the package doc and
// tools/bcftools/README.md.
package bcftools

import "math"

// =============================================================
//  Phred / log conversions (mirrors the inline helpers in trio-dnm3.c)
// =============================================================

// phred2num converts a phred-scaled value to a linear probability: 10^(-0.1*phred).
func phred2num(phred float64) float64 { return math.Pow(10, -0.1*phred) }

// log2phred converts a natural-log value to a phred score: |4.3429 * num|.
func log2phred(num float64) float64 { return math.Abs(4.3429 * num) }

// phred2log converts a phred-scaled value to a natural-log value: -phred/4.3429.
func phred2log(phred float64) float64 { return -phred / 4.3429 }

// subtractLog returns log(exp(aLog)-exp(bLog)), mirroring subtract_log().
func subtractLog(aLog, bLog float64) float64 {
	return aLog + math.Log(1-math.Exp(bLog-aLog))
}

// sumLog returns log(exp(a)+exp(b)), mirroring sum_log().
func sumLog(a, b float64) float64 {
	if math.IsInf(a, -1) && math.IsInf(b, -1) {
		return math.Inf(-1)
	}
	if a > b {
		return math.Log(1+math.Exp(b-a)) + a
	}
	return math.Log(1+math.Exp(a-b)) + b
}

// lfac is the log-space factorial log(n!) = lgamma(n+1), mirroring lfac().
func lfac(n float64) float64 { return kfLgamma(n + 1.0) }

// logMultinomCoeff returns the log multinomial coefficient for the counts,
// mirroring log_multinom_coeff_dbl().
func logMultinomCoeff(cnt []float64) float64 {
	sum, n := 0.0, 0.0
	for _, c := range cnt {
		n += c
		sum += lfac(c)
	}
	return lfac(n) - sum
}

// ldirichletMultinom is the Dirichlet-multinomial log PMF, mirroring
// ldirichlet_multinom_dbl(). phi is the overdispersion parameter.
func ldirichletMultinom(cnt, prob []float64, phi float64) float64 {
	n := 0.0
	for _, c := range cnt {
		n += c
	}
	sumAlpha := 0.0
	for _, p := range prob {
		sumAlpha += phi * p
	}
	ll := logMultinomCoeff(cnt) + kfLgamma(sumAlpha) - kfLgamma(n+sumAlpha)
	for i := range cnt {
		ll += kfLgamma(cnt[i]+phi*prob[i]) - kfLgamma(phi*prob[i])
	}
	return ll
}

// =============================================================
//  Priors (full pprob table) — mirrors init_priors() in trio-dnm3.c
// =============================================================

// countUnique counts the distinct alleles used by the genotype indices in gts,
// mirroring count_unique_alleles(). When includeRef is false the reference
// allele (0) is excluded.
func countUnique(gts []int, includeRef bool) int {
	var als [4]bool
	for _, igt := range gts {
		als[dnmSeq1[igt]] = true
		als[dnmSeq2[igt]] = true
	}
	nals := 0
	begin := 0
	if !includeRef {
		begin = 1
	}
	for i := begin; i < 4; i++ {
		if als[i] {
			nals++
		}
	}
	return nals
}

// newDNMPriorsFull builds the full denovo/denovo_allele/pprob tables for the
// given rule set, mirroring init_priors(). useDNGPriors selects the original
// (buggy-by-design) DeNovoGear prior assignment. mrate is the mutation rate.
func newDNMPriorsFull(strictlyNovel, useDNGPriors bool, mrate float64, typ dnmPriorsType) *dnmPriors {
	p := &dnmPriors{}
	for fi := 0; fi < 10; fi++ {
		for mi := 0; mi < 10; mi++ {
			for ci := 0; ci < 10; ci++ {
				var gtPrior float64
				switch {
				case useDNGPriors:
					gtPrior = initDNGMFPriors(fi, mi, ci)
				case typ == autosomalPriors:
					gtPrior = initMFPriors(fi, mi)
				case typ == chrXPriors:
					gtPrior = initMFPriorsChrX(mi)
				default: // chrXXPriors
					gtPrior = initMFPriorsChrXX(fi, mi)
				}

				var tprob, mprob float64
				var allele int
				switch {
				case useDNGPriors:
					tprob, mprob, allele = initDNGTprobMprob(mrate, fi, mi, ci)
				case typ == autosomalPriors:
					tprob, mprob, allele = initTprobMprob(strictlyNovel, mrate, fi, mi, ci)
				case typ == chrXPriors:
					tprob, mprob, allele = initTprobMprobChrX(strictlyNovel, mrate, fi, mi, ci)
				default:
					tprob, mprob, allele = initTprobMprobChrXX(strictlyNovel, mrate, fi, mi, ci)
				}

				if tprob == 0 {
					p.denovo[fi][mi][ci] = 1
					p.denovoAllele[fi][mi][ci] = allele
				} else {
					p.denovo[fi][mi][ci] = 0
					p.denovoAllele[fi][mi][ci] = 255
				}
				t := tprob
				if tprob == 0 {
					t = 1
				}
				p.pprob[fi][mi][ci] = math.Log(gtPrior * mprob * t)
			}
		}
	}
	return p
}

// initMFPriors returns the parent genotype prior L(GM,GF) for the autosomal rule
// set with DNG bugs fixed, mirroring init_mf_priors().
func initMFPriors(fi, mi int) float64 {
	fa, fb := dnmSeq1[fi], dnmSeq2[fi]
	ma, mb := dnmSeq1[mi], dnmSeq2[mi]
	naltMF := countUnique([]int{fi, mi}, false)
	nrefMF := b2i(fa == 0) + b2i(fb == 0) + b2i(ma == 0) + b2i(mb == 0)

	const pHomref = 0.998
	pPoly := (1 - pHomref) * (1 - pHomref)
	pNonref := 1 - pHomref - pPoly

	switch {
	case naltMF >= 3:
		return 1e-26
	case naltMF >= 2:
		return pPoly / 57.
	case nrefMF == 4:
		return pHomref
	case nrefMF == 3:
		return pNonref * (4.0 / 15.0) * (1.0 / 3.0)
	case nrefMF == 2 && ma == mb:
		return pNonref * (2.0 / 15.0) * (1.0 / 3.0)
	case nrefMF == 2:
		return pNonref * (4.0 / 15.0) * (1.0 / 3.0)
	case nrefMF == 1:
		return pNonref * (4.0 / 15.0) * (1.0 / 3.0)
	default: // nrefMF == 0
		return pNonref * (1.0 / 15.0) * (1.0 / 3.0)
	}
}

// initMFPriorsChrX returns the maternal genotype prior L(GM) for a male proband
// on chrX, mirroring init_mf_priors_chrX().
func initMFPriorsChrX(mi int) float64 {
	ma, mb := dnmSeq1[mi], dnmSeq2[mi]
	naltM := countUnique([]int{mi}, false)
	nrefM := b2i(ma == 0) + b2i(mb == 0)

	const pHomref = 0.999
	pPoly := (1 - pHomref) * (1 - pHomref)
	pNonref := 1 - pHomref - pPoly

	switch {
	case naltM >= 2:
		return pPoly / 3.
	case nrefM == 2:
		return pHomref
	case nrefM == 1:
		return pNonref * (2.0 / 3.0) * (1.0 / 3.0)
	default: // nrefM == 0
		return pNonref * (1.0 / 3.0) * (1.0 / 3.0)
	}
}

// initMFPriorsChrXX returns the parent genotype prior for a female proband on
// chrX, mirroring init_mf_priors_chrXX(). A heterozygous father yields prior 0.
func initMFPriorsChrXX(fi, mi int) float64 {
	fa, fb := dnmSeq1[fi], dnmSeq2[fi]
	ma, mb := dnmSeq1[mi], dnmSeq2[mi]
	naltMF := countUnique([]int{fi, mi}, false)
	nrefMF := b2i(fa == 0) + b2i(fb == 0) + b2i(ma == 0) + b2i(mb == 0)
	if fa != fb {
		return 0 // father can't be a het
	}
	if fa == 0 {
		nrefMF--
	} else {
		naltMF--
	}

	const pHomref = 0.998
	pPoly := (1 - pHomref) * (1 - pHomref)
	pNonref := 1 - pHomref - pPoly

	switch {
	case naltMF >= 3:
		return 1e-26
	case naltMF >= 2:
		return pPoly * (1.0 / 9.0) * (1.0 / 3.0)
	case nrefMF == 3:
		return pHomref
	case nrefMF == 2:
		return pNonref * (3.0 / 7.0) * (1.0 / 3.0)
	case nrefMF == 1:
		return pNonref * (3.0 / 7.0) * (1.0 / 3.0)
	default: // nrefMF == 0
		return pNonref * (1.0 / 7.0) * (1.0 / 3.0)
	}
}

// initDNGMFPriors returns the original DeNovoGear parent genotype prior,
// including the prior-assignment bugs, mirroring init_DNG_mf_priors().
func initDNGMFPriors(fi, mi, ci int) float64 {
	fa, fb := dnmSeq1[fi], dnmSeq2[fi]
	ma, mb := dnmSeq1[mi], dnmSeq2[mi]
	ca, cb := dnmSeq1[ci], dnmSeq2[ci]
	nalsMF := countUnique([]int{fi, mi}, true)
	nalsMFC := countUnique([]int{fi, mi, ci}, true)
	nrefMF := b2i(fa == 0) + b2i(fb == 0) + b2i(ma == 0) + b2i(mb == 0)

	switch {
	case nalsMFC > 3:
		return 1e-26
	case nalsMF >= 3:
		return 0.002 * 0.002 / 414
	case nalsMFC == 3:
		return 1e-3 * 1e-3
	case nrefMF == 4:
		return 0.995 * 0.998
	case nrefMF == 3:
		return 0.995 * 0.002 * (3.0 / 5.0) * (4.0 / 5.0) * 0.5
	case nrefMF == 2 && fa == fb && ma == mb:
		return 0.995 * 0.002 * (2.0 / 5.0) * (1.0 / 5.0) * 0.5
	case nrefMF == 2:
		return 0.995 * 0.002 * (2.0 / 5.0) * (2.0 / 5.0)
	case nrefMF == 1:
		return 0.995 * 0.002 * (2.0 / 5.0) * (2.0 / 5.0) * 0.5
	default: // nrefMF == 0
		if nalsMF == 1 {
			return 0.995 * 0.002 * (3.0 / 5.0) * (1.0 / 5.0)
		}
		// nalsMF == 2; ca,cb != 0 guaranteed by upstream's asserts
		_ = ca
		_ = cb
		return 0.002 * 0.002 / 414
	}
}

// initTprobMprob returns (tprob, mprob, denovoAllele) for the autosomal rule,
// mirroring init_tprob_mprob().
func initTprobMprob(strictlyNovel bool, mrate float64, fi, mi, ci int) (float64, float64, int) {
	fa, fb := dnmSeq1[fi], dnmSeq2[fi]
	ma, mb := dnmSeq1[mi], dnmSeq2[mi]
	ca, cb := dnmSeq1[ci], dnmSeq2[ci]

	var denovoAllele int
	if ca != fa && ca != fb && ca != ma && ca != mb {
		denovoAllele = ca
	} else {
		denovoAllele = cb
	}

	var isNovel bool
	if strictlyNovel {
		isNovel = (ca != fa && ca != fb && ca != ma && ca != mb) || (cb != fa && cb != fb && cb != ma && cb != mb)
		if isNovel && denovoAllele == 0 {
			isNovel = false
		}
	} else {
		isNovel = !(((ca == fa || ca == fb) && (cb == ma || cb == mb)) || ((ca == ma || ca == mb) && (cb == fa || cb == fb)))
	}

	if !isNovel {
		var tprob float64
		switch {
		case fa == fb && ma == mb:
			tprob = 1
		case fa == fb || ma == mb:
			tprob = 0.5
		default:
			tprob = 0.25
		}
		return tprob, 1 - mrate, denovoAllele
	}
	var mprob float64
	if (ca == fa || ca == fb) || (ca == ma || ca == mb) || (cb == fa || cb == fb) || (cb == ma || cb == mb) {
		mprob = mrate
	} else {
		mprob = mrate * mrate
	}
	return 0, mprob, denovoAllele
}

// initTprobMprobChrX returns (tprob, mprob, denovoAllele) for a male proband on
// chrX, mirroring init_tprob_mprob_chrX().
func initTprobMprobChrX(strictlyNovel bool, mrate float64, fi, mi, ci int) (float64, float64, int) {
	fa, fb := dnmSeq1[fi], dnmSeq2[fi]
	ma, mb := dnmSeq1[mi], dnmSeq2[mi]
	ca, cb := dnmSeq1[ci], dnmSeq2[ci]

	var denovoAllele int
	if ca != ma && ca != mb {
		denovoAllele = ca
	} else {
		denovoAllele = cb
	}

	if ca != cb { // male cannot be het in X, but can be mosaic
		var isNovel bool
		if strictlyNovel {
			isNovel = (ca != fa && ca != fb && ca != ma && ca != mb) || (cb != fa && cb != fb && cb != ma && cb != mb)
		} else {
			isNovel = (ca != ma && ca != mb) || (cb != ma && cb != mb)
		}
		if isNovel {
			return 0, mrate, denovoAllele
		}
		// genotype error: fall back to autosomal
		return initTprobMprob(strictlyNovel, mrate, fi, mi, ci)
	}
	if ca == ma || ca == mb { // inherited
		var tprob float64
		if ma == mb {
			tprob = 1
		} else {
			tprob = 0.5
		}
		return tprob, 1 - mrate, denovoAllele
	}
	// de novo
	return 0, mrate, denovoAllele
}

// initTprobMprobChrXX returns (tprob, mprob, denovoAllele) for a female proband
// on chrX, mirroring init_tprob_mprob_chrXX().
func initTprobMprobChrXX(strictlyNovel bool, mrate float64, fi, mi, ci int) (float64, float64, int) {
	fa, fb := dnmSeq1[fi], dnmSeq2[fi]
	ma, mb := dnmSeq1[mi], dnmSeq2[mi]
	ca, cb := dnmSeq1[ci], dnmSeq2[ci]

	var denovoAllele int
	if ca != fa && ca != fb && ca != ma && ca != mb {
		denovoAllele = ca
	} else {
		denovoAllele = cb
	}

	if fa != fb {
		// father can't be het in X; treat as genotype error, fall back to autosomal
		return initTprobMprob(strictlyNovel, mrate, fi, mi, ci)
	}
	if (ca == fa && (cb == ma || cb == mb)) || (cb == fa && (ca == ma || ca == mb)) {
		var tprob float64
		if ma == mb {
			tprob = 1
		} else {
			tprob = 0.5
		}
		return tprob, 1 - mrate, denovoAllele
	}
	var mprob float64
	if (ca == fa || (ca == ma || ca == mb)) || (cb == fa || (cb == ma || cb == mb)) {
		mprob = mrate
	} else {
		mprob = mrate * mrate
	}
	return 0, mprob, denovoAllele
}

// initDNGTprobMprob returns (tprob, mprob, denovoAllele) for the original
// DeNovoGear transmission model, mirroring init_DNG_tprob_mprob().
func initDNGTprobMprob(mrate float64, fi, mi, ci int) (float64, float64, int) {
	fa, fb := dnmSeq1[fi], dnmSeq2[fi]
	ma, mb := dnmSeq1[mi], dnmSeq2[mi]
	ca, cb := dnmSeq1[ci], dnmSeq2[ci]
	nalsMFC := countUnique([]int{fi, mi, ci}, true)

	tprob := 1.0
	mprob := 1 - mrate
	var denovoAllele int
	if ca != fa && ca != fb && ca != ma && ca != mb {
		denovoAllele = ca
	} else {
		denovoAllele = cb
	}

	switch {
	case nalsMFC == 4:
		tprob = 0
	case nalsMFC == 3:
		if ((ca == fa || ca == fb) && (cb == ma || cb == mb)) ||
			((cb == fa || cb == fb) && (ca == ma || ca == mb)) {
			switch {
			case ca == cb:
				tprob = 0.25
			case fa == fb || ma == mb:
				tprob = 0.5
			default:
				tprob = 0.25
			}
		} else {
			if ca != fa && ca != fb && ca != ma && ca != mb &&
				cb != fa && cb != fb && cb != ma && cb != mb {
				mprob = mrate * mrate
			} else {
				mprob = mrate
			}
			tprob = 0
		}
	case nalsMFC == 2:
		switch {
		case fa != fb && ma != mb:
			tprob = 0.25
		case fa == fb && ma == mb:
			switch {
			case fa == ma && ca == cb:
				tprob, mprob = 0, mrate*mrate
			case fa == ma:
				tprob, mprob = 0, mrate
			case ca == cb:
				tprob, mprob = 0, mrate
			}
		case ca == cb && ((fa == fb && fa != ca) || (ma == mb && ma != ca)):
			tprob, mprob = 0, mrate
		default:
			tprob = 0.5
		}
	}
	return tprob, mprob, denovoAllele
}
