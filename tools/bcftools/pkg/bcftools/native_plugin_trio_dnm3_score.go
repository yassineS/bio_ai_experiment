// Per-trio scoring for the trio-dnm3 float models and their noise-model
// helpers, porting process_trio_DMM / process_trio_ALM / process_trio_DNG and
// the site_noise / parental_emission / missed_inherited / mosaic_noise helpers
// from plugins/trio-dnm3.c. See native_plugin_trio_dnm3_models.go for the
// libm-tolerance note.
package bcftools

import "math"

// pnoise mirrors pnoise_t: the tolerated fraction/count of unexpected parental
// reads, with the *1 variants applied to alleles observed in a single parent.
type pnoise struct {
	abs, frac   float64
	abs1, frac1 float64
}

// dnmModelParams bundles the float-model knobs shared by the process_trio_*
// functions, mirroring the relevant args_t fields.
type dnmModelParams struct {
	phi         float64
	minQM       float64 // signed; the model uses fabs(minQM) as the per-read error
	minVAF      float64
	noisePrior  float64
	allelicDrop float64
	withCAD     bool
	withPPL     bool
	strictNovel bool
	pnCur       pnoise // current site (snv or indel)
}

// ldirichletMultinomWithSpurious mirrors ldirichlet_multinom_with_spurious():
// a Dirichlet-multinomial likelihood for a parental genotype that tolerates a
// bounded fraction/count of unexpected (spurious) alternate reads. als is the
// allele bitmask of the genotype; strict is the bitmask of alleles observed in
// both parents (for which the stricter pn is used).
func ldirichletMultinomWithSpurious(p *dnmModelParams, dp []int, err []float64, als, strict int) float64 {
	cnt := [3]float64{}
	prob := [3]float64{}
	fst := true
	for i := 0; i < len(dp); i++ {
		inAls := (als>>i)&1 == 1
		switch {
		case !inAls:
			cnt[2] += float64(dp[i])
			if err == nil {
				prob[2] = math.Abs(p.minQM)
			} else if prob[2] < err[i] {
				prob[2] = err[i]
			}
		case fst:
			cnt[0] = float64(dp[i])
			if err != nil {
				prob[0] = 1.0 - err[i]
			} else {
				prob[0] = 1 - math.Abs(p.minQM)
			}
			fst = false
		default:
			cnt[1] = float64(dp[i])
			if err != nil {
				prob[1] = 1.0 - err[i]
			} else {
				prob[1] = 1 - math.Abs(p.minQM)
			}
		}
	}
	if cnt[0]+cnt[1] == 0 {
		return math.Inf(-1)
	}
	const eps = 1e-12
	sum := 0.0
	for i := 0; i < 3; i++ {
		if prob[i] < eps {
			prob[i] = eps
		} else if prob[i] > 1-eps {
			prob[i] = 1 - eps
		}
		sum += prob[i]
	}
	for i := 0; i < 3; i++ {
		prob[i] /= sum
	}

	phi := p.phi
	var eta, kmax float64
	if (^als)&strict != 0 {
		eta, kmax = p.pnCur.frac, p.pnCur.abs
	} else {
		eta, kmax = p.pnCur.frac1, p.pnCur.abs1
	}
	ncnt := cnt[0] + cnt[1] + cnt[2]
	if kmax < eta*ncnt {
		kmax = eta * ncnt
	}
	if kmax > cnt[2] {
		kmax = cnt[2]
	}
	cnt[2] -= kmax
	return ldirichletMultinom(cnt[:], prob[:], phi)
}

// normProb3 clamps and renormalises a 3-element probability vector, mirroring
// norm_prob().
func normProb3(prob []float64) {
	const eps = 1e-12
	sum := 0.0
	for i := 0; i < 3; i++ {
		if prob[i] < eps {
			prob[i] = eps
		} else if prob[i] > 1-eps {
			prob[i] = 1 - eps
		}
		sum += prob[i]
	}
	for i := 0; i < 3; i++ {
		prob[i] /= sum
	}
}

// ldirichletMultinomAD is the AD-based child genotype likelihood, mirroring
// ldirichlet_multinom_AD() (used by DMM with --with-cAD).
func ldirichletMultinomAD(p *dnmModelParams, dp []int, als int) float64 {
	cnt := [3]float64{}
	prob := [3]float64{}
	nals := 0
	for i := 0; i < len(dp); i++ {
		inAls := (als>>i)&1 == 1
		switch {
		case !inAls:
			cnt[2] += float64(dp[i])
		case nals == 0:
			cnt[0] = float64(dp[i])
			nals++
		default:
			cnt[1] = float64(dp[i])
			nals++
		}
	}
	if cnt[0]+cnt[1] == 0 {
		return math.Inf(-1)
	}
	if nals == 1 {
		prob[0] = 1 - math.Abs(p.minQM)
		prob[1] = math.Abs(p.minQM)
		prob[2] = math.Abs(p.minQM)
		normProb3(prob[:])
		return ldirichletMultinom(cnt[:], prob[:], p.phi)
	}
	llMax := math.Inf(-1)
	if cnt[0] < cnt[1] {
		cnt[0], cnt[1] = cnt[1], cnt[0]
	}
	for vaf := 0.5; vaf >= p.minVAF; vaf -= 0.1 {
		prob[0] = (1 - math.Abs(p.minQM)) * (1 - vaf)
		prob[1] = (1 - math.Abs(p.minQM)) * vaf
		prob[2] = math.Abs(p.minQM)
		normProb3(prob[:])
		ll := ldirichletMultinom(cnt[:], prob[:], p.phi)
		if ll < llMax {
			return llMax
		}
		llMax = ll
	}
	return llMax
}

// siteNoise is the upper binomial-tail penalty for unexpected reads at a site,
// mirroring site_noise(). als is the called genotype bitmask.
func siteNoise(ad []int, als int, err, lprior float64) float64 {
	n, k := 0, 0
	for i := 0; i < len(ad); i++ {
		if (1<<i)&als == 0 {
			k += ad[i]
		}
		n += ad[i]
	}
	if n <= 0 || k <= 0 {
		return 0
	}
	if err <= 0 {
		return 0
	}
	if err > 0.5 {
		err = 0.5
	}
	pTail := kfBetai(float64(k), float64(n-k+1), err)
	if pTail < 1e-300 {
		pTail = 1e-300
	}
	lprob := math.Log(pTail) + lprior
	if lprob > 0 {
		lprob = 0
	}
	return lprob
}

// parentalEmission penalises a candidate de-novo allele that is emitted by the
// parents, mirroring parental_emission(). ad holds [father, mother, child]
// per-allele depths; only father/mother are used.
func parentalEmission(p *dnmModelParams, ad [3][]int, ial int) float64 {
	npe := b2i(ad[0][ial] != 0) + b2i(ad[1][ial] != 0)
	if npe == 0 {
		return 0
	}
	err := math.Abs(p.minQM)
	var pval [2]float64
	for i := 0; i < 2; i++ { // father and mother
		if ad[i][ial] == 0 {
			continue
		}
		k := ad[i][ial]
		n := 0
		for _, v := range ad[i] {
			n += v
		}
		var eta, kmax float64
		if npe == 2 {
			eta, kmax = p.pnCur.frac, p.pnCur.abs
		} else {
			eta, kmax = p.pnCur.frac1, p.pnCur.abs1
		}
		if kmax < eta*float64(n) {
			kmax = eta * float64(n)
		}
		if kmax > float64(k) {
			kmax = float64(k)
		}
		// Upstream stores n and k as int; `n -= kmax; k -= kmax` with a double
		// kmax promotes to double, subtracts, then truncates the RESULT back to
		// int (e.g. 18 - 0.81 -> 17), not kmax.
		n = int(float64(n) - kmax)
		k = int(float64(k) - kmax)
		if k == 0 {
			continue
		}
		pv := kfBetai(float64(k), float64(n-k+1), err)
		if pv < 1e-300 {
			pval[i] = math.Log(1e-300)
		} else {
			pval[i] = math.Log(pv)
		}
	}
	return pval[0] + pval[1]
}

// missedInherited mirrors missed_inherited(): a correction for an inherited
// allele potentially missed due to low depth (ALM --allelic-dropout).
func missedInherited(ad []int, idenovo int) float64 {
	n, k := 0, 0
	for i := 0; i < len(ad); i++ {
		if i == idenovo {
			k = ad[i]
		}
		n += ad[i]
	}
	pTail := kfBetai(float64(n-k), float64(k+1), 0.5)
	if pTail < 1e-300 {
		pTail = 1e-300
	}
	lprob := math.Log(pTail)
	if lprob > 0 {
		lprob = 0
	}
	return subtractLog(0, lprob)
}

// mosaicScoreCDFSmooth mirrors mosaic_score_cdf_smooth(): the Beta-Binomial
// posterior probability that the minor-allele fraction exceeds m0.
func mosaicScoreCDFSmooth(c0, c1 int, m0 float64, ncap int, a0, b0 float64) float64 {
	n := c0 + c1
	if n <= 0 {
		return 0.0
	}
	mhat := float64(c1) / float64(n)
	ne := float64(n)
	if n >= ncap {
		ne = float64(ncap)
	}
	ke := ne * mhat
	if ke < 0 {
		ke = 0
	}
	if ke > ne {
		ke = ne
	}
	a := a0 + ke
	b := b0 + (ne - ke)
	cdf := kfBetai(a, b, m0)
	score := 1.0 - cdf
	if score < 0.0 {
		score = 0.0
	}
	if score > 1.0 {
		score = 1.0
	}
	return score
}

// mosaicNoise mirrors mosaic_noise(): downplays de-novo calls whose child VAF is
// too low to be a confident mosaic, returning a log penalty.
func mosaicNoise(p *dnmModelParams, ad []int, als, idenovo int) float64 {
	var cnt [2]int
	for i := 0; i < len(ad); i++ {
		if (1<<i)&als == 0 {
			continue
		}
		idx := 0
		if i == idenovo {
			idx = 1
		}
		cnt[idx] = ad[i]
	}
	if cnt[1] == 0 {
		return 0
	}
	a0, b0 := 1.0, 1.0
	m0 := p.minVAF
	ncap := 50
	score := mosaicScoreCDFSmooth(cnt[0], cnt[1], m0, ncap, a0, b0)
	if score < 1e-300 {
		score = 1e-300
	}
	return math.Log(score)
}

// processTrioDMM mirrors process_trio_DMM(): the default Dirichlet-multinomial
// model. pl/ad/qm are the per-member values ordered [father, mother, child]
// (indices iFATHER=0, iMOTHER=1, iCHILD=2). It returns the log de-novo score and
// sets al0/al1 (the two child alleles, al1 being the de-novo allele).
func processTrioDMM(p *dnmModelParams, priors *dnmPriors, nals int, pl [3][]float64, ad [3][]int, qm [3][]float64) (score float64, al0, al1 int) {
	strict := 0
	for i := 0; i < nals; i++ {
		if ad[iFATHER][i] != 0 && ad[iMOTHER][i] != 0 {
			strict |= 1 << i
		}
	}

	sum := math.Inf(-1)
	maxDNM := math.Inf(-1)
	var maxGT [3]int
	ci := 0
	for ca := 0; ca < nals; ca++ {
		for cb := 0; cb <= ca; cb++ {
			cals := (1 << ca) | (1 << cb)
			var cpl float64
			if p.withCAD {
				cpl = ldirichletMultinomAD(p, ad[iCHILD], cals)
			} else {
				cpl = pl[iCHILD][ci]
			}
			fi := 0
			for fa := 0; fa < nals; fa++ {
				for fb := 0; fb <= fa; fb++ {
					fals := (1 << fa) | (1 << fb)
					fpl := ldirichletMultinomWithSpurious(p, ad[iFATHER], qm[iFATHER], fals, strict)
					mi := 0
					for ma := 0; ma < nals; ma++ {
						for mb := 0; mb <= ma; mb++ {
							mals := (1 << ma) | (1 << mb)
							mpl := ldirichletMultinomWithSpurious(p, ad[iMOTHER], qm[iMOTHER], mals, strict)
							val := cpl + mpl + fpl + priors.pprob[fi][mi][ci]
							sum = sumLog(sum, val)
							if priors.denovo[fi][mi][ci] != 0 && maxDNM < val {
								maxDNM = val
								maxGT[iCHILD] = cals
								maxGT[iMOTHER] = mals
								maxGT[iFATHER] = fals
								if priors.denovoAllele[fi][mi][ci] == ca {
									al0, al1 = cb, ca
								} else {
									al0, al1 = ca, cb
								}
							}
							mi++
						}
					}
					fi++
				}
			}
			ci++
		}
	}
	if math.IsInf(sum, 0) {
		return math.Inf(-1), al0, al1
	}

	lnoise := 0.0
	if p.noisePrior != 0 {
		err := math.Abs(p.minQM)
		lprior := math.Log(1e6 * math.Abs(p.noisePrior))
		lnoise += siteNoise(ad[iMOTHER], maxGT[iMOTHER], err, lprior)
		lnoise += siteNoise(ad[iFATHER], maxGT[iFATHER], err, lprior)
		if p.noisePrior > 0 {
			lnoise += siteNoise(ad[iCHILD], maxGT[iCHILD], err, lprior)
		}
	}
	lnoise += parentalEmission(p, ad, al1)
	if p.minVAF != 0 {
		lnoise += mosaicNoise(p, ad[iCHILD], maxGT[iCHILD], al1)
	}
	return maxDNM + lnoise - sum, al0, al1
}

// processTrioALM mirrors process_trio_ALM(): the allele-likelihood model over
// FORMAT/QS (or PL with --with-pPL). qs holds per-member log-prob QS.
func processTrioALM(p *dnmModelParams, priors *dnmPriors, nals int, pl [3][]float64, ad [3][]int, qs [3][]float64, isChrX bool) (score float64, al0, al1 int) {
	sum := math.Inf(-1)
	maxDNM := math.Inf(-1)
	var maxGT [3]int
	ci := 0
	for ca := 0; ca < nals; ca++ {
		for cb := 0; cb <= ca; cb++ {
			cals := (1 << ca) | (1 << cb)
			cpl := pl[iCHILD][ci]
			fi := 0
			for fa := 0; fa < nals; fa++ {
				for fb := 0; fb <= fa; fb++ {
					fals := (1 << fa) | (1 << fb)
					var fpl float64
					if p.withPPL {
						fpl = pl[iFATHER][fi]
					} else {
						for i := 0; i < nals; i++ {
							if fals&(1<<i) != 0 {
								fpl += subtractLog(0, qs[iFATHER][i])
							} else {
								fpl += qs[iFATHER][i]
							}
						}
					}
					mi := 0
					for ma := 0; ma < nals; ma++ {
						for mb := 0; mb <= ma; mb++ {
							mals := (1 << ma) | (1 << mb)
							var mpl float64
							if p.withPPL {
								mpl = pl[iMOTHER][mi]
							} else {
								for i := 0; i < nals; i++ {
									if mals&(1<<i) != 0 {
										mpl += subtractLog(0, qs[iMOTHER][i])
									} else {
										mpl += qs[iMOTHER][i]
									}
								}
							}
							val := cpl + fpl + mpl + priors.pprob[fi][mi][ci]
							sum = sumLog(sum, val)
							if priors.denovo[fi][mi][ci] != 0 && maxDNM < val {
								maxDNM = val
								maxGT[iCHILD] = cals
								maxGT[iMOTHER] = mals
								maxGT[iFATHER] = fals
								if priors.denovoAllele[fi][mi][ci] == ca {
									al0, al1 = cb, ca
								} else {
									al0, al1 = ca, cb
								}
							}
							mi++
						}
					}
					fi++
				}
			}
			ci++
		}
	}

	// ALM's noise helpers all iterate over nals (the trimmed allele count), unlike
	// DMM which iterates over the original nad. Slice the (possibly full-length,
	// stale-tailed) per-member AD arrays to nals so the stale tail is never read.
	adN := [3][]int{}
	for i := 0; i < 3; i++ {
		if ad[i] != nil {
			adN[i] = ad[i][:nals]
		}
	}
	lnoise := 0.0
	if p.noisePrior != 0 {
		err := math.Abs(p.minQM)
		lprior := math.Log(1e6 * math.Abs(p.noisePrior))
		lnoise += siteNoise(adN[iMOTHER], maxGT[iMOTHER], err, lprior)
		lnoise += siteNoise(adN[iFATHER], maxGT[iFATHER], err, lprior)
		if p.noisePrior > 0 {
			lnoise += siteNoise(adN[iCHILD], maxGT[iCHILD], err, lprior)
		}
	}
	if p.allelicDrop != 0 {
		// missed_inherited(args, nals, ad[iX], max_gt[iX], *al1): only the last
		// argument (idenovo == al1) is used by the function body.
		lnoise += missedInherited(adN[iMOTHER], al1)
		if !isChrX {
			lnoise += missedInherited(adN[iFATHER], al1)
		}
	}
	if p.minVAF != 0 {
		lnoise += mosaicNoise(p, adN[iCHILD], maxGT[iCHILD], al1)
	}
	if p.strictNovel {
		ial := al1
		if qs[iMOTHER][ial]+qs[iFATHER][ial] != 0 {
			tmp := 0.0
			if qs[iMOTHER][ial] != 0 {
				tmp += subtractLog(0, qs[iMOTHER][ial])
			}
			if qs[iFATHER][ial] != 0 {
				tmp += subtractLog(0, qs[iFATHER][ial])
			}
			sum = sumLog(sum, tmp)
			maxDNM += tmp
		}
	}
	return maxDNM + lnoise - sum, al0, al1
}

// processTrioDNG mirrors process_trio_DNG(): the original DeNovoGear model over
// FORMAT/PL only.
func processTrioDNG(priors *dnmPriors, nals int, pl [3][]float64) (score float64, al0, al1 int) {
	sum := math.Inf(-1)
	maxv := math.Inf(-1)
	ci := 0
	for ca := 0; ca < nals; ca++ {
		for cb := 0; cb <= ca; cb++ {
			fi := 0
			for fa := 0; fa < nals; fa++ {
				for fb := 0; fb <= fa; fb++ {
					mi := 0
					for ma := 0; ma < nals; ma++ {
						for mb := 0; mb <= ma; mb++ {
							val := pl[iCHILD][ci] + pl[iFATHER][fi] + pl[iMOTHER][mi] + priors.pprob[fi][mi][ci]
							sum = sumLog(val, sum)
							if priors.denovo[fi][mi][ci] != 0 && maxv < val {
								maxv = val
								if priors.denovoAllele[fi][mi][ci] == ca {
									al0, al1 = cb, ca
								} else {
									al0, al1 = ca, cb
								}
							}
							mi++
						}
					}
					fi++
				}
			}
			ci++
		}
	}
	return maxv - sum, al0, al1
}

// iFATHER, iMOTHER, iCHILD index the per-member arrays, matching the C macros.
const (
	iFATHER = 0
	iMOTHER = 1
	iCHILD  = 2
)
