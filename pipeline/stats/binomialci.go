// Package stats provides small, dependency-free statistical helpers used to
// attach confidence intervals to the project's parity rates.
//
// A parity rate is a binomial proportion: out of n compared cells, k matched
// the upstream oracle byte-for-byte (or within the floating-point tolerance).
// Reporting the point estimate k/n alone hides the sample size, so every rate
// in the manuscript is accompanied by a 95% confidence interval computed here.
//
// Two interval estimators are provided:
//
//   - WilsonCI: the Wilson score interval. A good default — well-behaved for
//     small n and for proportions near 0 or 1 (where the textbook
//     normal-approximation "Wald" interval is badly miscalibrated and can run
//     outside [0,1]).
//   - ClopperPearsonCI: the exact (Clopper-Pearson) interval derived from the
//     binomial CDF via its Beta-distribution relationship. It is guaranteed to
//     have at least the nominal coverage (it is conservative), which is the
//     right choice for the headline "k/k passed" cells where Wilson can look
//     optimistically narrow.
//
// Everything here is pure Go, standard library only.
package stats

import "math"

// WilsonCI returns the lower and upper bounds of the Wilson score confidence
// interval for a binomial proportion of successes out of n trials.
//
// z is the standard-normal quantile for the desired two-sided coverage: use
// z = 1.959963985 for a 95% interval, 2.575829304 for 99%. The interval is
// centred on a shrinkage-adjusted estimate rather than the raw k/n, which is
// why it stays inside [0,1] and remains sensible at the boundaries (e.g. it
// gives a non-degenerate lower bound for the all-success case k == n, where the
// Wald interval collapses to a point).
//
// For n == 0 it returns the whole unit interval (0, 1): with no data, no
// proportion can be excluded.
func WilsonCI(successes, n int, z float64) (lo, hi float64) {
	if n <= 0 {
		return 0, 1
	}
	nf := float64(n)
	phat := float64(successes) / nf
	z2 := z * z
	// denom = 1 + z^2/n
	denom := 1 + z2/nf
	// centre = (phat + z^2/(2n)) / denom
	centre := (phat + z2/(2*nf)) / denom
	// half-width = (z/denom) * sqrt( phat(1-phat)/n + z^2/(4n^2) )
	margin := (z / denom) * math.Sqrt(phat*(1-phat)/nf+z2/(4*nf*nf))
	lo = centre - margin
	hi = centre + margin
	if lo < 0 {
		lo = 0
	}
	if hi > 1 {
		hi = 1
	}
	return lo, hi
}

// ClopperPearsonCI returns the lower and upper bounds of the exact
// (Clopper-Pearson) confidence interval for a binomial proportion of successes
// out of n trials, at two-sided confidence level 1-alpha (pass alpha = 0.05 for
// a 95% interval).
//
// The interval is defined through the inverse of the regularized incomplete
// Beta function (the Beta quantile), using the standard identities:
//
//	lower = BetaQuantile(alpha/2,     k,   n-k+1)
//	upper = BetaQuantile(1-alpha/2,   k+1, n-k)
//
// with the boundary conventions lower = 0 when k == 0 and upper = 1 when
// k == n. These follow from the Beta–Binomial relationship
// F_Binomial(k-1; n, p) = I_{1-p}(n-k+1, k); inverting the binomial tail
// probabilities for p yields exactly the Beta quantiles above. The interval is
// conservative (coverage >= 1-alpha), which is what we want for the headline
// all-pass cells.
//
// For n == 0 it returns the whole unit interval (0, 1).
func ClopperPearsonCI(successes, n int, alpha float64) (lo, hi float64) {
	if n <= 0 {
		return 0, 1
	}
	k := successes
	if k <= 0 {
		lo = 0
	} else {
		lo = betaQuantile(alpha/2, float64(k), float64(n-k+1))
	}
	if k >= n {
		hi = 1
	} else {
		hi = betaQuantile(1-alpha/2, float64(k+1), float64(n-k))
	}
	return lo, hi
}

// betaQuantile returns x in [0,1] such that the regularized incomplete Beta
// function I_x(a, b) == p. It inverts incompleteBeta by bisection, which is
// robust over the full parameter range we use (a, b are positive integers up to
// a few hundred) and needs no derivative. ~60 iterations of bisection drive the
// bracket below 1e-12, far tighter than the three significant figures we report.
func betaQuantile(p, a, b float64) float64 {
	if p <= 0 {
		return 0
	}
	if p >= 1 {
		return 1
	}
	lo, hi := 0.0, 1.0
	for i := 0; i < 200; i++ {
		mid := (lo + hi) / 2
		if incompleteBeta(mid, a, b) < p {
			lo = mid
		} else {
			hi = mid
		}
		if hi-lo < 1e-12 {
			break
		}
	}
	return (lo + hi) / 2
}

// incompleteBeta returns the regularized incomplete Beta function I_x(a, b),
// i.e. the CDF of a Beta(a, b) distribution evaluated at x. It uses the
// Lentz continued-fraction expansion (Numerical Recipes §6.4), with the
// standard symmetry reflection I_x(a,b) = 1 - I_{1-x}(b,a) applied when x is on
// the slowly-converging side of the distribution mean (a+1)/(a+b+2).
func incompleteBeta(x, a, b float64) float64 {
	if x <= 0 {
		return 0
	}
	if x >= 1 {
		return 1
	}
	// ln of the leading factor x^a (1-x)^b / (a B(a,b)).
	lnFront := lgamma(a+b) - lgamma(a) - lgamma(b) +
		a*math.Log(x) + b*math.Log(1-x)
	front := math.Exp(lnFront)
	if x < (a+1)/(a+b+2) {
		return front * betaCF(x, a, b) / a
	}
	// Reflect to the fast-converging tail.
	return 1 - front*betaCF(1-x, b, a)/b
}

// betaCF evaluates the continued fraction for the incomplete Beta function via
// the modified Lentz algorithm. It converges for x < (a+1)/(a+b+2); callers
// reflect the argument otherwise (see incompleteBeta).
func betaCF(x, a, b float64) float64 {
	const (
		maxIter = 300
		eps     = 1e-14
		tiny    = 1e-300
	)
	qab := a + b
	qap := a + 1
	qam := a - 1
	c := 1.0
	d := 1 - qab*x/qap
	if math.Abs(d) < tiny {
		d = tiny
	}
	d = 1 / d
	h := d
	for m := 1; m <= maxIter; m++ {
		mf := float64(m)
		// even step
		aa := mf * (b - mf) * x / ((qam + 2*mf) * (a + 2*mf))
		d = 1 + aa*d
		if math.Abs(d) < tiny {
			d = tiny
		}
		c = 1 + aa/c
		if math.Abs(c) < tiny {
			c = tiny
		}
		d = 1 / d
		h *= d * c
		// odd step
		aa = -(a + mf) * (qab + mf) * x / ((a + 2*mf) * (qap + 2*mf))
		d = 1 + aa*d
		if math.Abs(d) < tiny {
			d = tiny
		}
		c = 1 + aa/c
		if math.Abs(c) < tiny {
			c = tiny
		}
		d = 1 / d
		del := d * c
		h *= del
		if math.Abs(del-1) < eps {
			break
		}
	}
	return h
}

// lgamma returns the natural log of the absolute value of the Gamma function,
// discarding the sign (always +1 for the positive arguments used here).
func lgamma(x float64) float64 {
	v, _ := math.Lgamma(x)
	return v
}
