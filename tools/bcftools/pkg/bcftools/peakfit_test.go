package bcftools

import (
	"math"
	"testing"
)

// --- Levenberg-Marquardt solver ---

// quadModel fits y = p0 + p1*x + p2*x^2 to data, an exactly-solvable
// linear-in-parameters least-squares problem the LM solver must nail.
type quadModel struct {
	xs, ys []float64
}

func (m *quadModel) nResiduals() int { return len(m.xs) }
func (m *quadModel) nParams() int    { return 3 }
func (m *quadModel) evaluate(p, res []float64, jac [][]float64) {
	for i, x := range m.xs {
		model := p[0] + p[1]*x + p[2]*x*x
		res[i] = model - m.ys[i]
		jac[i][0] = 1
		jac[i][1] = x
		jac[i][2] = x * x
	}
}

func TestLM_QuadraticExactFit(t *testing.T) {
	// True parameters: 2 - 3x + 1.5x^2.
	want := []float64{2, -3, 1.5}
	xs := make([]float64, 41)
	ys := make([]float64, 41)
	for i := range xs {
		x := -2 + 0.1*float64(i)
		xs[i] = x
		ys[i] = want[0] + want[1]*x + want[2]*x*x
	}
	m := &quadModel{xs: xs, ys: ys}
	res := levenbergMarquardt(m, []float64{0, 0, 0}, 500)
	for j := range want {
		if math.Abs(res.Params[j]-want[j]) > 1e-6 {
			t.Errorf("param[%d] = %v, want %v", j, res.Params[j], want[j])
		}
	}
	if res.SSR > 1e-12 {
		t.Errorf("SSR = %v, want ~0", res.SSR)
	}
}

// gaussModel fits a single Gaussian a*exp(-((x-b)/c)^2) (a non-linear
// problem) to exercise the LM damping loop properly.
type gaussModel struct {
	xs, ys []float64
}

func (m *gaussModel) nResiduals() int { return len(m.xs) }
func (m *gaussModel) nParams() int    { return 3 }
func (m *gaussModel) evaluate(p, res []float64, jac [][]float64) {
	a, b, c := p[0], p[1], p[2]
	for i, x := range m.xs {
		t := (x - b) / c
		e := math.Exp(-t * t)
		res[i] = a*e - m.ys[i]
		jac[i][0] = e
		jac[i][1] = a * e * 2 * (x - b) / (c * c)
		jac[i][2] = a * e * 2 * (x - b) * (x - b) / (c * c * c)
	}
}

func TestLM_GaussianNonlinearFit(t *testing.T) {
	want := []float64{3.0, 0.4, 0.08}
	xs := make([]float64, 101)
	ys := make([]float64, 101)
	for i := range xs {
		x := 0.01 * float64(i)
		xs[i] = x
		tt := (x - want[1]) / want[2]
		ys[i] = want[0] * math.Exp(-tt*tt)
	}
	m := &gaussModel{xs: xs, ys: ys}
	// Start near, but not at, the optimum.
	res := levenbergMarquardt(m, []float64{1.0, 0.5, 0.05}, 500)
	for j := range want {
		if math.Abs(res.Params[j]-want[j]) > 1e-4 {
			t.Errorf("param[%d] = %v, want %v", j, res.Params[j], want[j])
		}
	}
}

func TestSolveLinear(t *testing.T) {
	// 2x + y = 5 ; x + 3y = 10  =>  x=1, y=3.
	a := [][]float64{{2, 1}, {1, 3}}
	b := []float64{5, 10}
	if !solveLinear(a, b) {
		t.Fatal("solveLinear reported singular for a solvable system")
	}
	if math.Abs(b[0]-1) > 1e-12 || math.Abs(b[1]-3) > 1e-12 {
		t.Errorf("solveLinear = %v, want [1 3]", b)
	}
}

func TestSolveLinearSingular(t *testing.T) {
	a := [][]float64{{1, 2}, {2, 4}} // rank deficient
	b := []float64{3, 6}
	if solveLinear(a, b) {
		t.Error("solveLinear should report a singular matrix")
	}
}

// --- glibc rand() ---

// TestGlibcRand checks the first values of srand(0) against the known
// glibc random() sequence (seed 0 is mapped to 1 by glibc).
func TestGlibcRand(t *testing.T) {
	r := newGlibcRand(0)
	// Reference values from glibc: srand(0); rand() x5.
	want := []int32{1804289383, 846930886, 1681692777, 1714636915, 1957747793}
	for i, w := range want {
		if got := r.next(); got != w {
			t.Errorf("rand()#%d = %d, want %d", i, got, w)
		}
	}
}

func TestGlibcRandUniform(t *testing.T) {
	r := newGlibcRand(0)
	for i := 0; i < 1000; i++ {
		v := r.uniform(0.2, 0.8)
		if v < 0.2 || v > 0.8 {
			t.Fatalf("uniform out of range: %v", v)
		}
	}
}

// --- peak models ---

// TestGaussianModel fits a peakfit Gaussian to a clean synthetic
// Gaussian sample and checks the recovered scale/centre/sigma.
func TestGaussianModel(t *testing.T) {
	// Target: scale^2=4 (scale=2), centre=0.5, sigma=0.06.
	xs := make([]float64, 101)
	ys := make([]float64, 101)
	for i := range xs {
		x := 0.01 * float64(i)
		xs[i] = x
		t := (x - 0.5) / 0.06
		ys[i] = 4 * math.Exp(-t*t)
	}
	pf := newPeakfit()
	pf.addGaussian(1.5, 0.45, 0.05, 7)
	fit := pf.run(xs, ys)
	scale, center, sigma := pf.getParams(0)
	if math.Abs(scale*scale-4) > 1e-2 {
		t.Errorf("scale^2 = %v, want 4", scale*scale)
	}
	if math.Abs(center-0.5) > 1e-3 {
		t.Errorf("center = %v, want 0.5", center)
	}
	if math.Abs(sigma-0.06) > 1e-3 {
		t.Errorf("sigma = %v, want 0.06", sigma)
	}
	if fit > 1e-2 {
		t.Errorf("fit residual = %v, want ~0", fit)
	}
}

// TestBoundedGaussianModel checks the cosine re-parametrisation keeps
// the fitted centre inside the [d,e] interval and still recovers it.
func TestBoundedGaussianModel(t *testing.T) {
	xs := make([]float64, 101)
	ys := make([]float64, 101)
	for i := range xs {
		x := 0.01 * float64(i)
		xs[i] = x
		t := (x - 0.52) / 0.05
		ys[i] = 2.25 * math.Exp(-t*t) // scale^2=2.25
	}
	pf := newPeakfit()
	pf.addBoundedGaussian(1.0, 0.5, 0.04, 0.45, 0.55, 7)
	pf.run(xs, ys)
	scale, center, sigma := pf.getParams(0)
	if center < 0.45 || center > 0.55 {
		t.Errorf("center %v escaped the [0.45,0.55] bound", center)
	}
	if math.Abs(center-0.52) > 5e-3 {
		t.Errorf("center = %v, want ~0.52", center)
	}
	if math.Abs(scale*scale-2.25) > 5e-2 {
		t.Errorf("scale^2 = %v, want ~2.25", scale*scale)
	}
	if sigma <= 0 {
		t.Errorf("sigma = %v, want positive", sigma)
	}
}

// TestBoundedGaussianClamp confirms a centre initialised outside the
// interval is clamped, not NaN-ed.
func TestBoundedGaussianClamp(t *testing.T) {
	c := boundedCenterSet(0.9, 0.45, 0.55)
	if math.IsNaN(c) {
		t.Fatal("boundedCenterSet produced NaN for an out-of-range centre")
	}
	back := boundedCenterGet(c, 0.45, 0.55)
	if back < 0.45 || back > 0.55 {
		t.Errorf("clamped centre %v escaped [0.45,0.55]", back)
	}
}

// TestExpModel fits the exponential peak model.
func TestExpModel(t *testing.T) {
	// y = scale^2 * exp((x-center)/sigma^2); scale^2=1, center=1, sigma=0.3.
	xs := make([]float64, 51)
	ys := make([]float64, 51)
	for i := range xs {
		x := 0.6 + 0.008*float64(i)
		xs[i] = x
		ys[i] = 1.0 * math.Exp((x-1.0)/(0.3*0.3))
	}
	pf := newPeakfit()
	pf.addExp(0.8, 1.0, 0.25, 5)
	fit := pf.run(xs, ys)
	if fit > 5e-2 {
		t.Errorf("exp fit residual = %v, want small", fit)
	}
}

// TestPeakfitEvaluateUnscaled confirms the reported fit is the
// UNSCALED Σ|model-y|, matching upstream peakfit_evaluate.
func TestPeakfitEvaluateUnscaled(t *testing.T) {
	xs := []float64{0, 1}
	ys := []float64{0, 0}
	pf := newPeakfit()
	// A fixed peak (no free params): scale^2=1 at centre 0, sigma huge
	// so the peak is ~1 everywhere. With no free params run() just
	// evaluates: Σ|model-y| ≈ |1-0| + |1-0| = 2.
	pf.addGaussian(1.0, 0.0, 1e6, 0)
	got := pf.run(xs, ys)
	if math.Abs(got-2) > 1e-3 {
		t.Errorf("evaluate fit = %v, want ~2 (unscaled abs sum)", got)
	}
}
