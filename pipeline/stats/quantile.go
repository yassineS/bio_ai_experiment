// Quantile and ratio-confidence-interval helpers for the performance side of
// the parity pipeline. Where binomialci.go attaches confidence intervals to
// *parity rates* (a binomial proportion), this file attaches them to *timing
// measurements* (a set of per-repetition samples): the manuscript's speed claim
// (C3) needs the full distribution of each cell's runtime, not a single reduced
// number, plus a confidence interval on the ours/upstream speed ratio.
//
// Everything here is pure Go, standard library only.
package stats

import (
	"math"
	"math/rand"
	"sort"
)

// Quantile returns the p-quantile (0 <= p <= 1) of xs using linear
// interpolation between the two nearest order statistics. This is the
// "type 7" quantile of Hyndman & Fan (1996) — the default in R's quantile()
// and in NumPy's percentile() (method="linear") — so the values reproduce
// those reference implementations exactly.
//
// For a sorted sample x[0..n-1] the type-7 quantile is computed at the
// fractional rank h = (n-1)*p:
//
//	Q(p) = x[floor(h)] + (h - floor(h)) * (x[floor(h)+1] - x[floor(h)])
//
// xs need not be sorted on entry; a sorted copy is taken so the caller's slice
// is left untouched. An empty slice returns NaN; a single element returns that
// element for every p.
func Quantile(xs []float64, p float64) float64 {
	n := len(xs)
	if n == 0 {
		return math.NaN()
	}
	if n == 1 {
		return xs[0]
	}
	if p <= 0 {
		return minOf(xs)
	}
	if p >= 1 {
		return maxOf(xs)
	}
	s := make([]float64, n)
	copy(s, xs)
	sort.Float64s(s)
	h := float64(n-1) * p
	lo := int(math.Floor(h))
	frac := h - float64(lo)
	if lo+1 >= n {
		return s[n-1]
	}
	return s[lo] + frac*(s[lo+1]-s[lo])
}

// MedianIQR returns the median (Q2), the first quartile (Q1) and the third
// quartile (Q3) of xs, each computed with the type-7 interpolating quantile
// (see Quantile). The inter-quartile range is q3-q1; callers that want it can
// subtract. An empty slice yields NaN for all three.
func MedianIQR(xs []float64) (median, q1, q3 float64) {
	if len(xs) == 0 {
		return math.NaN(), math.NaN(), math.NaN()
	}
	// Sort once and reuse for all three quantiles.
	s := make([]float64, len(xs))
	copy(s, xs)
	sort.Float64s(s)
	return quantileSorted(s, 0.5), quantileSorted(s, 0.25), quantileSorted(s, 0.75)
}

// quantileSorted is Quantile for an already-sorted slice (no defensive copy).
func quantileSorted(s []float64, p float64) float64 {
	n := len(s)
	if n == 0 {
		return math.NaN()
	}
	if n == 1 {
		return s[0]
	}
	if p <= 0 {
		return s[0]
	}
	if p >= 1 {
		return s[n-1]
	}
	h := float64(n-1) * p
	lo := int(math.Floor(h))
	frac := h - float64(lo)
	if lo+1 >= n {
		return s[n-1]
	}
	return s[lo] + frac*(s[lo+1]-s[lo])
}

// BootstrapSeed is the fixed PRNG seed used by RatioCI so the reported interval
// is deterministic: the same samples always yield the same bounds, which keeps
// the manuscript numbers reproducible and the bench output diff-stable across
// runs. (math/rand with an explicit, constant source — not the global,
// process-time-seeded generator.)
const BootstrapSeed = 0x5eed1a1

// defaultBootstrapResamples is the number of bootstrap resamples RatioCI draws.
// 2000 is comfortably above the rule-of-thumb minimum (~1000) for a stable
// percentile interval while staying cheap for the handful of cells per report.
const defaultBootstrapResamples = 2000

// RatioCI returns a percentile bootstrap confidence interval for the ratio
//
//	r = median(ours) / median(upstream)
//
// at two-sided confidence level 1-alpha (pass alpha = 0.05 for a 95%
// interval), together with the point estimate r itself.
//
// Method: a paired/independent two-sample bootstrap. On each of B iterations it
// resamples the `ours` reps with replacement and the `upstream` reps with
// replacement (independently — the two sides are separate processes, not paired
// observations), forms the ratio of the two resampled medians, and records it.
// The CI is the empirical [alpha/2, 1-alpha/2] percentile range of those B
// ratios (Efron's percentile interval). Resampling uses a math/rand source
// seeded from the constant BootstrapSeed, so the result is fully deterministic.
//
// Edge cases:
//   - If either side is empty, or the upstream point median is <= 0 (so the
//     ratio is undefined), it returns point = NaN and an empty [NaN, NaN]
//     interval.
//   - With a single rep per side there is nothing to resample, so every
//     bootstrap ratio equals the point estimate and the interval collapses to
//     [r, r]; this is reported faithfully rather than hidden.
func RatioCI(ours, upstream []float64, alpha float64) (point, lo, hi float64) {
	return ratioCIWith(ours, upstream, alpha, defaultBootstrapResamples, BootstrapSeed)
}

// ratioCIWith is RatioCI with the resample count and seed exposed, so tests can
// pin both for exact, reproducible expectations.
func ratioCIWith(ours, upstream []float64, alpha float64, resamples int, seed int64) (point, lo, hi float64) {
	if len(ours) == 0 || len(upstream) == 0 {
		return math.NaN(), math.NaN(), math.NaN()
	}
	mo, _, _ := MedianIQR(ours)
	mu, _, _ := MedianIQR(upstream)
	if mu <= 0 {
		return math.NaN(), math.NaN(), math.NaN()
	}
	point = mo / mu

	if resamples < 1 {
		resamples = 1
	}
	rng := rand.New(rand.NewSource(seed))
	ratios := make([]float64, 0, resamples)
	bo := make([]float64, len(ours))
	bu := make([]float64, len(upstream))
	for i := 0; i < resamples; i++ {
		for j := range bo {
			bo[j] = ours[rng.Intn(len(ours))]
		}
		for j := range bu {
			bu[j] = upstream[rng.Intn(len(upstream))]
		}
		rmu, _, _ := MedianIQR(bu)
		if rmu <= 0 {
			continue
		}
		rmo, _, _ := MedianIQR(bo)
		ratios = append(ratios, rmo/rmu)
	}
	if len(ratios) == 0 {
		return point, math.NaN(), math.NaN()
	}
	lo = Quantile(ratios, alpha/2)
	hi = Quantile(ratios, 1-alpha/2)
	return point, lo, hi
}

func minOf(xs []float64) float64 {
	m := xs[0]
	for _, x := range xs[1:] {
		if x < m {
			m = x
		}
	}
	return m
}

func maxOf(xs []float64) float64 {
	m := xs[0]
	for _, x := range xs[1:] {
		if x > m {
			m = x
		}
	}
	return m
}
