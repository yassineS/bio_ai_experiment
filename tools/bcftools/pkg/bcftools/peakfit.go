// Package bcftools — see doc.go. This file ports upstream's
// `peakfit.c`: the peak-fitting engine that `bcftools polysomy` uses to
// fit Gaussian-mixture models to a chromosome's B-allele-frequency
// (BAF) histogram.
//
// Upstream's `peakfit.c` builds a sum of parametric peaks (Gaussian,
// centre-bounded Gaussian, exponential) and fits the free parameters to
// the observed histogram with GSL's Levenberg-Marquardt non-linear
// least-squares solver (`gsl_multifit_fdfsolver_lmsder`), wrapped in a
// Monte-Carlo restart loop that re-seeds selected parameters from a
// uniform range to escape local minima. The LM solver itself is ported
// in `peakfit_lm.go`; this file ports the peak models, the residual
// objective, the Monte-Carlo driver, and the public peakfit API.
//
// Faithful-port notes:
//
//   - The residual handed to the LM solver is `(model-y)/0.01`, exactly
//     as upstream's `peakfit_calc_f` scales it. The reported fit and the
//     value used to pick the best Monte-Carlo restart, however, use the
//     UNSCALED `Σ|model-y|` (upstream's `peakfit_evaluate`); both are
//     reproduced here verbatim so the CN decision thresholds line up.
//   - Monte-Carlo restarts are seeded from the in-tree glibc `rand()`
//     port (`peakfit_lm.go`) after `srand(0)`, matching upstream's
//     reproducibility guarantee.
//   - The bounded Gaussian re-parametrises its centre through a cosine
//     so the optimiser sees an unconstrained variable; the forward and
//     inverse transforms match `bounded_gaussian_convert_*`.
package bcftools

import (
	"fmt"
	"math"
	"strings"
)

// nPeakParams is the per-peak parameter slot count (upstream NPARAMS).
// Slots 0..2 are scale/centre/sigma; slots 3..4 hold the bounded
// Gaussian's interval endpoints.
const nPeakParams = 5

// peakKind enumerates the parametric peak shapes peakfit supports.
type peakKind int

const (
	peakGaussian        peakKind = iota // scale^2 * exp(-(x-c)^2/sigma^2)
	peakBoundedGaussian                 // Gaussian with centre bound to [d,e]
	peakExp                             // scale^2 * exp((x-c)/sigma^2)
)

// mcSetting records the Monte-Carlo restart range for one parameter
// slot: when scan is true the slot is re-seeded uniformly in [min,max]
// at the start of each restart.
type mcSetting struct {
	scan     bool
	min, max float64
	best     float64
}

// peak is one parametric component of a peakfit model.
type peak struct {
	kind    peakKind
	fitMask int                  // bit i set => params[i] is a free LM variable
	params  [nPeakParams]float64 // current working parameters
	ori     [nPeakParams]float64 // initial parameters (re-applied per restart)
	mc      [nPeakParams]mcSetting
}

// peakfit is the porcelain over a multi-peak fit. Construct one with
// newPeakfit, add peaks, optionally attach Monte-Carlo ranges, then call
// run.
type peakfit struct {
	peaks   []peak
	nparams int // total free parameters across all peaks
	nmcIter int // Monte-Carlo restart count (0 => single fit)
}

// newPeakfit returns an empty peakfit.
func newPeakfit() *peakfit { return &peakfit{} }

// reset clears the peak list so the same peakfit can be reused for the
// next candidate CN model, matching upstream's peakfit_reset.
func (pf *peakfit) reset() {
	pf.peaks = pf.peaks[:0]
	pf.nparams = 0
	pf.nmcIter = 0
}

// countFit returns the number of free parameters implied by a fit mask.
func countFit(mask int) int {
	n := 0
	for i := 0; i < nPeakParams; i++ {
		if mask&(1<<uint(i)) != 0 {
			n++
		}
	}
	return n
}

// addGaussian appends an unbounded Gaussian peak with initial
// scale=a, centre=b, sigma=c and the given free-parameter mask.
func (pf *peakfit) addGaussian(a, b, c float64, fitMask int) {
	p := peak{kind: peakGaussian, fitMask: fitMask}
	p.ori[0], p.ori[1], p.ori[2] = a, b, c
	pf.peaks = append(pf.peaks, p)
	pf.nparams += countFit(fitMask)
}

// addBoundedGaussian appends a Gaussian whose centre is constrained to
// the interval [d,e]. The constraint is enforced by a cosine
// re-parametrisation so the LM optimiser still sees an unconstrained
// variable, exactly as upstream's bounded_gaussian model.
func (pf *peakfit) addBoundedGaussian(a, b, c, d, e float64, fitMask int) {
	if !(d < e) {
		// Upstream asserts d<e; degenerate intervals collapse to a
		// point so the optimiser cannot move the centre.
		e = d + 1e-9
	}
	p := peak{kind: peakBoundedGaussian, fitMask: fitMask}
	p.ori[0], p.ori[2], p.ori[3], p.ori[4] = a, c, d, e
	p.ori[1] = boundedCenterSet(b, d, e)
	pf.peaks = append(pf.peaks, p)
	pf.nparams += countFit(fitMask)
}

// addExp appends an exponential peak with initial scale=a, centre=b,
// sigma=c. Upstream forbids fitting the centre of an exp peak, so bit 1
// is masked off defensively.
func (pf *peakfit) addExp(a, b, c float64, fitMask int) {
	// Deliberate divergence from upstream peakfit_add_exp, which does
	// `assert(!(fit_mask&2))` and aborts the process on a bad mask.
	// Silently clearing the bit is the right call for a library: a
	// caller cannot crash the host program with a stray flag.
	fitMask &^= 1 << 1
	p := peak{kind: peakExp, fitMask: fitMask}
	p.ori[0], p.ori[1], p.ori[2] = a, b, c
	pf.peaks = append(pf.peaks, p)
	pf.nparams += countFit(fitMask)
}

// setMC attaches a Monte-Carlo restart range to parameter slot iparam
// of the most recently added peak and records the restart count, in the
// same way as upstream's peakfit_set_mc.
func (pf *peakfit) setMC(xmin, xmax float64, iparam, niter int) {
	if len(pf.peaks) == 0 || iparam < 0 || iparam >= nPeakParams {
		return
	}
	p := &pf.peaks[len(pf.peaks)-1]
	p.mc[iparam].scan = true
	p.mc[iparam].min = xmin
	p.mc[iparam].max = xmax
	pf.nmcIter = niter
}

// boundedCenterSet maps a real centre value into the unconstrained
// cosine coordinate the optimiser works in (upstream
// bounded_gaussian_convert_set). The centre is clamped to [d,e] first.
func boundedCenterSet(value, d, e float64) float64 {
	if value < d {
		value = d
	} else if value > e {
		value = e
	}
	arg := 2*(value-d)/(e-d) - 1
	if arg < -1 {
		arg = -1
	} else if arg > 1 {
		arg = 1
	}
	return math.Acos(arg)
}

// boundedCenterGet maps the cosine coordinate back to a real centre in
// [d,e] (upstream bounded_gaussian_convert_get).
func boundedCenterGet(center, d, e float64) float64 {
	return 0.5*(math.Cos(center)+1)*(e-d) + d
}

// evalPeak adds peak p's contribution at every x into out (out is
// accumulated, not overwritten — peaks sum).
func evalPeak(p *peak, xvals, out []float64) {
	switch p.kind {
	case peakGaussian:
		scale2 := p.params[0] * p.params[0]
		center := p.params[1]
		sigma := p.params[2]
		for i, x := range xvals {
			t := (x - center) / sigma
			out[i] += scale2 * math.Exp(-t*t)
		}
	case peakBoundedGaussian:
		scale2 := p.params[0] * p.params[0]
		z := boundedCenterGet(p.params[1], p.params[3], p.params[4])
		sigma := p.params[2]
		for i, x := range xvals {
			t := (x - z) / sigma
			out[i] += scale2 * math.Exp(-t*t)
		}
	case peakExp:
		scale2 := p.params[0] * p.params[0]
		center := p.params[1]
		sigma := p.params[2]
		for i, x := range xvals {
			out[i] += scale2 * math.Exp((x-center)/sigma/sigma)
		}
	}
}

// derivPeak adds peak p's partial derivative ∂(peak)/∂(slot idf) at
// every x into out. Matches upstream's *_calc_df.
func derivPeak(p *peak, xvals, out []float64, idf int) {
	switch p.kind {
	case peakGaussian:
		scale := p.params[0]
		center := p.params[1]
		sigma := p.params[2]
		for i, x := range xvals {
			zi := x - center
			expv := math.Exp(-zi * zi / (sigma * sigma))
			switch idf {
			case 0:
				out[i] += 2 * scale * expv
			case 1:
				out[i] += 2 * scale * scale * zi * expv / (sigma * sigma)
			case 2:
				out[i] += 2 * scale * scale * zi * zi * expv / (sigma * sigma * sigma)
			}
		}
	case peakBoundedGaussian:
		scale := p.params[0]
		center := p.params[1]
		sigma := p.params[2]
		d, e := p.params[3], p.params[4]
		z := boundedCenterGet(center, d, e)
		for i, x := range xvals {
			zi := x - z
			expv := math.Exp(-zi * zi / sigma / sigma)
			switch idf {
			case 0:
				out[i] += 2 * scale * expv
			case 1:
				out[i] -= scale * scale * math.Sin(center) * (e - d) * zi * expv / sigma / sigma
			case 2:
				out[i] += 2 * scale * scale * zi * zi * expv / sigma / sigma / sigma
			}
		}
	case peakExp:
		scale := p.params[0]
		center := p.params[1]
		sigma := p.params[2]
		for i, x := range xvals {
			expv := math.Exp((x - center) / sigma / sigma)
			switch idf {
			case 0:
				out[i] += 2 * scale * expv
			case 2:
				out[i] -= 2 * scale * scale * (x - center) * expv / sigma / sigma / sigma
			}
		}
	}
}

// peakfitModel adapts a peakfit to the lmModel interface so the LM
// solver can drive it. The residual is upstream's (model-y)/0.01.
type peakfitModel struct {
	pf    *peakfit
	xvals []float64
	yvals []float64
	// freeMap[k] = (peakIdx, slot) for the k-th free LM variable, in
	// the same peak-major order GSL packs them.
	freeMap [][2]int
	scratch []float64
}

const residualScale = 0.01

// newPeakfitModel builds the lmModel adapter for a peakfit and data.
func newPeakfitModel(pf *peakfit, xvals, yvals []float64) *peakfitModel {
	m := &peakfitModel{pf: pf, xvals: xvals, yvals: yvals}
	for pi := range pf.peaks {
		for slot := 0; slot < nPeakParams; slot++ {
			if pf.peaks[pi].fitMask&(1<<uint(slot)) != 0 {
				m.freeMap = append(m.freeMap, [2]int{pi, slot})
			}
		}
	}
	m.scratch = make([]float64, len(xvals))
	return m
}

func (m *peakfitModel) nResiduals() int { return len(m.xvals) }
func (m *peakfitModel) nParams() int    { return len(m.freeMap) }

// evaluate fills residuals and the Jacobian for the given parameter
// vector. residuals[i] = (model(x_i)-y_i)/0.01; jacobian[i][k] is the
// derivative of that residual with respect to the k-th free variable.
func (m *peakfitModel) evaluate(params, residuals []float64, jacobian [][]float64) {
	// Push free params into the peaks.
	for k, fm := range m.freeMap {
		m.pf.peaks[fm[0]].params[fm[1]] = params[k]
	}
	// Model values.
	for i := range m.scratch {
		m.scratch[i] = 0
	}
	for pi := range m.pf.peaks {
		evalPeak(&m.pf.peaks[pi], m.xvals, m.scratch)
	}
	for i := range residuals {
		residuals[i] = (m.scratch[i] - m.yvals[i]) / residualScale
	}
	// Jacobian: column k is ∂residual/∂param_k = (1/0.01)·∂peak/∂slot.
	for k, fm := range m.freeMap {
		for i := range m.scratch {
			m.scratch[i] = 0
		}
		derivPeak(&m.pf.peaks[fm[0]], m.xvals, m.scratch, fm[1])
		for i := range jacobian {
			jacobian[i][k] = m.scratch[i] / residualScale
		}
	}
}

// evaluateFit returns upstream's peakfit_evaluate: the UNSCALED
// Σ|model(x_i)-y_i| for the current peak parameters. This is the value
// used both to rank Monte-Carlo restarts and as the fit reported to the
// CN decision logic.
func (pf *peakfit) evaluateFit(xvals, yvals []float64) float64 {
	vals := make([]float64, len(xvals))
	for pi := range pf.peaks {
		evalPeak(&pf.peaks[pi], xvals, vals)
	}
	var sum float64
	for i := range vals {
		sum += math.Abs(vals[i] - yvals[i])
	}
	return sum
}

// run fits the assembled peak model to (xvals,yvals) and returns the
// best UNSCALED Σ|model-y| fit. It mirrors upstream peakfit_run: srand(0),
// then nmcIter+1 Monte-Carlo restarts each seeding scanned parameter
// slots from their uniform range, running the LM solver to convergence,
// and keeping the parameter set with the smallest fit.
func (pf *peakfit) run(xvals, yvals []float64) float64 {
	if pf.nparams == 0 {
		// No free parameters: a pure evaluation, as upstream does.
		for pi := range pf.peaks {
			pf.applyOri(pi)
		}
		return pf.evaluateFit(xvals, yvals)
	}

	rng := newGlibcRand(0)
	model := newPeakfitModel(pf, xvals, yvals)
	bestFit := math.Inf(1)
	bestParams := make([][nPeakParams]float64, len(pf.peaks))

	for iter := 0; iter <= pf.nmcIter; iter++ {
		// Seed every peak: start from ori, then re-seed scanned slots.
		init := make([]float64, 0, pf.nparams)
		for pi := range pf.peaks {
			p := &pf.peaks[pi]
			for slot := 0; slot < nPeakParams; slot++ {
				p.params[slot] = p.ori[slot]
				if p.mc[slot].scan {
					v := rng.uniform(p.mc[slot].min, p.mc[slot].max)
					if p.kind == peakBoundedGaussian && slot == 1 {
						v = boundedCenterSet(v, p.ori[3], p.ori[4])
					}
					p.params[slot] = v
				}
				if p.fitMask&(1<<uint(slot)) != 0 {
					init = append(init, p.params[slot])
				}
			}
		}

		res := levenbergMarquardt(model, init, 500)
		// Push the converged free parameters back into the peaks.
		for k, fm := range model.freeMap {
			pf.peaks[fm[0]].params[fm[1]] = res.Params[k]
		}
		fit := pf.evaluateFit(xvals, yvals)
		// Intentional, more-correct deviation from upstream
		// peakfit_run (peakfit.c:584-600): the C code snapshots the
		// best parameters and updates best_fit in two SEPARATE `if
		// (fit<best_fit)` statements, so the second test sees the
		// stale best_fit and the snapshot/best_fit pair can drift out
		// of sync on a tie. Here a single branch snapshots
		// bestParams and bestFit together, so the recovered
		// parameters always belong to the recorded best fit. Do not
		// "restore" the split form.
		if fit < bestFit {
			bestFit = fit
			for pi := range pf.peaks {
				bestParams[pi] = pf.peaks[pi].params
			}
		}
	}

	for pi := range pf.peaks {
		pf.peaks[pi].params = bestParams[pi]
	}
	return bestFit
}

// applyOri copies a peak's initial parameters into its working slots,
// used for the no-free-parameter evaluation path.
func (pf *peakfit) applyOri(pi int) {
	pf.peaks[pi].params = pf.peaks[pi].ori
}

// getParams returns the converged (scale, centre, sigma) of peak ipk in
// the same human-readable convention as upstream's peakfit_get_params:
// scale and sigma are returned as absolute values, and the bounded
// Gaussian's centre is mapped back out of its cosine coordinate.
func (pf *peakfit) getParams(ipk int) (scale, center, sigma float64) {
	p := &pf.peaks[ipk]
	switch p.kind {
	case peakBoundedGaussian:
		scale = math.Abs(p.params[0])
		sigma = math.Abs(p.params[2])
		center = boundedCenterGet(p.params[1], p.params[3], p.params[4])
	default:
		scale = math.Abs(p.params[0])
		center = math.Abs(p.params[1])
		sigma = math.Abs(p.params[2])
	}
	return scale, center, sigma
}

// sprintFunc renders the fitted model as a human-readable expression,
// matching upstream's peakfit_sprint_func (used in the dist.dat dump).
func (pf *peakfit) sprintFunc() string {
	var parts []string
	for ipk := range pf.peaks {
		p := &pf.peaks[ipk]
		switch p.kind {
		case peakExp:
			parts = append(parts, fmt.Sprintf("%f**2 * exp((x-%f)/%f**2)",
				math.Abs(p.params[0]), math.Abs(p.params[1]), math.Abs(p.params[2])))
		case peakBoundedGaussian:
			z := boundedCenterGet(p.params[1], p.params[3], p.params[4])
			parts = append(parts, fmt.Sprintf("%f**2 * exp(-(x-%f)**2/%f**2)",
				math.Abs(p.params[0]), z, math.Abs(p.params[2])))
		default:
			parts = append(parts, fmt.Sprintf("%f**2 * exp(-(x-%f)**2/%f**2)",
				math.Abs(p.params[0]), p.params[1], math.Abs(p.params[2])))
		}
	}
	return strings.Join(parts, " + ")
}
