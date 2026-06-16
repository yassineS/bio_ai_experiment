// Binomial-tail helpers shared by the native plugins, built on the existing
// in-tree kfunc port (kfBetai / kfLgamma in callc.go). These mirror
// calc_binom_one_sided() and calc_binom_two_sided() in bcftools.h.
//
// The key property is that kfBetai routes through the deterministic AS245
// kf_lgamma polynomial (callc.go::kfLgamma) and the Lentz continued fraction,
// NOT the platform libm lgamma/tgamma. So as long as the elementary log/exp/pow
// operations agree with C — which they do on linux/amd64 — these tails are
// byte-reproducible against upstream, which is what lets the parental-origin
// plugin match upstream's %e DBG probabilities exactly.
package bcftools

// calcBinomTwoSided reproduces calc_binom_two_sided() in bcftools.h: the
// two-sided binomial tail probability of observing na vs nb successes under
// success probability aprob.
func calcBinomTwoSided(na, nb int, aprob float64) float64 {
	if na == 0 && nb == 0 {
		return -1
	}
	if na == nb {
		return 1
	}
	var prob float64
	if na > nb {
		prob = 2 * kfBetai(float64(na), float64(nb)+1, aprob)
	} else {
		prob = 2 * kfBetai(float64(nb), float64(na)+1, aprob)
	}
	if prob > 1 {
		prob = 1
	}
	return prob
}

// calcBinomOneSided reproduces calc_binom_one_sided() in bcftools.h: the
// one-sided binomial tail. When ge is true it returns I_aprob(na, nb+1);
// otherwise I_{1-aprob}(nb, na+1).
func calcBinomOneSided(na, nb int, aprob float64, ge bool) float64 {
	if ge {
		return kfBetai(float64(na), float64(nb)+1, aprob)
	}
	return kfBetai(float64(nb), float64(na)+1, 1-aprob)
}
