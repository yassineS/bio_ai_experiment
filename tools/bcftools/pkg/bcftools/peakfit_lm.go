// Package bcftools — see doc.go. This file ports the numerical core
// that upstream `bcftools polysomy` borrows from GSL: a
// Levenberg-Marquardt non-linear least-squares solver and the glibc
// `rand()` pseudo-random generator used by `peakfit.c` for its
// Monte-Carlo restarts.
//
// Upstream's `peakfit.c` calls GSL's `gsl_multifit_fdfsolver_lmsder`.
// CLAUDE.md scopes the one sanctioned numerical dependency (gonum)
// narrowly and explicitly excludes stats-fitting tools like this, so
// the solver is ported in-tree as focused pure Go. It is the textbook
// LM damping loop over the normal equations
//
//	(JᵀJ + λ·diag(JᵀJ))·δ = Jᵀr
//
// with an analytic Jacobian supplied by the caller, a λ up/down
// schedule, and convergence tests on the residual and parameter
// deltas that mirror GSL's `gsl_multifit_test_delta` /
// `gsl_multifit_test_gradient` (both with tolerance 1e-8).
package bcftools

import "math"

// glibcRand reproduces glibc's TYPE_3 additive-feedback `rand()`
// (the default `random()` generator used by the C standard library
// on Linux). `peakfit.c` calls `srand(0)` before each fit and seeds
// its Monte-Carlo restarts from the resulting stream, so reproducing
// the exact sequence keeps the Go port's restart trajectory
// deterministic and faithful to upstream.
type glibcRand struct {
	table [344]int32
	idx   int
}

// newGlibcRand returns a generator seeded exactly as glibc's
// `srand(seed)` would seed `random()`.
func newGlibcRand(seed uint32) *glibcRand {
	r := &glibcRand{}
	if seed == 0 {
		seed = 1
	}
	r.table[0] = int32(seed)
	for i := 1; i < 31; i++ {
		// hi/lo split avoids 64-bit overflow, matching glibc.
		hi := r.table[i-1] / 127773
		lo := r.table[i-1] % 127773
		word := 16807*lo - 2836*hi
		if word < 0 {
			word += 2147483647
		}
		r.table[i] = word
	}
	for i := 31; i < 34; i++ {
		r.table[i] = r.table[i-31]
	}
	for i := 34; i < 344; i++ {
		r.table[i] = r.table[i-31] + r.table[i-3]
	}
	r.idx = 344
	return r
}

// next returns the next pseudo-random value in [0, 2^31).
func (r *glibcRand) next() int32 {
	if r.idx >= len(r.table) {
		// Shift the trailing 31-word state to the front and keep
		// generating; this mirrors glibc's circular buffer.
		copy(r.table[0:31], r.table[len(r.table)-31:])
		r.idx = 31
	}
	i := r.idx
	r.table[i] = r.table[i-31] + r.table[i-3]
	r.idx++
	return int32(uint32(r.table[i]) >> 1)
}

// randMax is glibc's RAND_MAX, the inclusive upper bound of rand().
const randMax = 2147483647

// float01 mirrors `rand()*(max-min)/RAND_MAX + min` from peakfit.c.
func (r *glibcRand) uniform(min, max float64) float64 {
	return float64(r.next())*(max-min)/randMax + min
}

// lmModel is the residual model the LM solver optimises. It returns,
// for the current parameter vector, the residual vector r (length m)
// and the m×n Jacobian J (J[i][j] = ∂r_i/∂p_j). The solver minimises
// Σ r_i².
type lmModel interface {
	// nResiduals reports m, the number of residual rows.
	nResiduals() int
	// nParams reports n, the number of free parameters.
	nParams() int
	// evaluate fills residuals (len m) and jacobian (m×n) for the
	// supplied parameter vector.
	evaluate(params []float64, residuals []float64, jacobian [][]float64)
}

// lmResult holds the outcome of a Levenberg-Marquardt run.
type lmResult struct {
	Params []float64 // best-fit parameter vector
	SSR    float64   // sum of squared residuals at Params
	Iters  int       // iterations performed
}

// levenbergMarquardt minimises Σ r_i² for model starting at params0.
// It follows GSL's lmsder behaviour closely enough for the polysomy
// peak fits: an analytic Jacobian, a λ damping schedule that grows on
// rejected steps and shrinks on accepted ones, and convergence when
// either the relative parameter delta or the gradient falls below
// 1e-8. maxIter caps the iteration count (GSL's polysomy loop uses
// 500).
func levenbergMarquardt(model lmModel, params0 []float64, maxIter int) lmResult {
	n := model.nParams()
	m := model.nResiduals()
	p := make([]float64, n)
	copy(p, params0)

	residuals := make([]float64, m)
	jac := make([][]float64, m)
	for i := range jac {
		jac[i] = make([]float64, n)
	}
	model.evaluate(p, residuals, jac)
	ssr := dotResiduals(residuals)

	// JᵀJ and Jᵀr scratch.
	jtj := make([][]float64, n)
	for i := range jtj {
		jtj[i] = make([]float64, n)
	}
	jtr := make([]float64, n)
	delta := make([]float64, n)
	trial := make([]float64, n)
	trialRes := make([]float64, m)
	trialJac := make([][]float64, m)
	for i := range trialJac {
		trialJac[i] = make([]float64, n)
	}

	lambda := 1e-3
	const (
		lambdaUp   = 10.0
		lambdaDown = 10.0
		tol        = 1e-8
	)
	iters := 0
	for ; iters < maxIter; iters++ {
		// Normal equations: JᵀJ and Jᵀr.
		for a := 0; a < n; a++ {
			jtr[a] = 0
			for b := 0; b < n; b++ {
				jtj[a][b] = 0
			}
		}
		for i := 0; i < m; i++ {
			ri := residuals[i]
			row := jac[i]
			for a := 0; a < n; a++ {
				jtr[a] += row[a] * ri
				for b := a; b < n; b++ {
					jtj[a][b] += row[a] * row[b]
				}
			}
		}
		for a := 0; a < n; a++ {
			for b := 0; b < a; b++ {
				jtj[a][b] = jtj[b][a]
			}
		}

		// Gradient convergence test (GSL test_gradient): gradient is
		// 2·Jᵀr; we test against tol directly on Jᵀr's max element.
		gradMax := 0.0
		for a := 0; a < n; a++ {
			if v := math.Abs(jtr[a]); v > gradMax {
				gradMax = v
			}
		}
		if gradMax < tol {
			break
		}

		// Inner loop: grow λ until a step reduces the SSR.
		accepted := false
		for inner := 0; inner < 30; inner++ {
			// (JᵀJ + λ·diag(JᵀJ))·δ = -Jᵀr
			aug := make([][]float64, n)
			for a := 0; a < n; a++ {
				aug[a] = make([]float64, n)
				copy(aug[a], jtj[a])
				diag := jtj[a][a]
				if diag == 0 {
					diag = 1
				}
				aug[a][a] += lambda * diag
				delta[a] = -jtr[a]
			}
			if !solveLinear(aug, delta) {
				lambda *= lambdaUp
				continue
			}
			for a := 0; a < n; a++ {
				trial[a] = p[a] + delta[a]
			}
			model.evaluate(trial, trialRes, trialJac)
			trialSSR := dotResiduals(trialRes)
			if trialSSR < ssr && !math.IsNaN(trialSSR) && !math.IsInf(trialSSR, 0) {
				// Accept the step.
				dxNorm, xNorm := 0.0, 0.0
				for a := 0; a < n; a++ {
					dxNorm += delta[a] * delta[a]
					xNorm += trial[a] * trial[a]
				}
				copy(p, trial)
				copy(residuals, trialRes)
				for i := 0; i < m; i++ {
					copy(jac[i], trialJac[i])
				}
				ssr = trialSSR
				lambda /= lambdaDown
				if lambda < 1e-12 {
					lambda = 1e-12
				}
				accepted = true
				// Parameter-delta convergence (GSL test_delta).
				if math.Sqrt(dxNorm) < tol*(math.Sqrt(xNorm)+tol) {
					iters++
					return lmResult{Params: p, SSR: ssr, Iters: iters}
				}
				break
			}
			lambda *= lambdaUp
			if lambda > 1e12 {
				accepted = false
				break
			}
		}
		if !accepted {
			break
		}
	}
	return lmResult{Params: p, SSR: ssr, Iters: iters}
}

// dotResiduals returns Σ r_i².
func dotResiduals(r []float64) float64 {
	var s float64
	for _, v := range r {
		s += v * v
	}
	return s
}

// solveLinear solves a·x = b in place (b becomes x) by Gaussian
// elimination with partial pivoting. It returns false if a is
// singular. a is destroyed.
func solveLinear(a [][]float64, b []float64) bool {
	n := len(b)
	for col := 0; col < n; col++ {
		// Partial pivot.
		pivot := col
		best := math.Abs(a[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(a[r][col]); v > best {
				best = v
				pivot = r
			}
		}
		if best < 1e-300 {
			return false
		}
		if pivot != col {
			a[col], a[pivot] = a[pivot], a[col]
			b[col], b[pivot] = b[pivot], b[col]
		}
		inv := 1.0 / a[col][col]
		for r := col + 1; r < n; r++ {
			f := a[r][col] * inv
			if f == 0 {
				continue
			}
			for c := col; c < n; c++ {
				a[r][c] -= f * a[col][c]
			}
			b[r] -= f * b[col]
		}
	}
	// Back-substitution.
	for r := n - 1; r >= 0; r-- {
		s := b[r]
		for c := r + 1; c < n; c++ {
			s -= a[r][c] * b[c]
		}
		b[r] = s / a[r][r]
	}
	return true
}
