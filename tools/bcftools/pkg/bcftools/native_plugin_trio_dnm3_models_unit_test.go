// Binary-free unit tests for the trio-dnm3 float building blocks
// (native_plugin_trio_dnm3_models.go) and a pinned known-input -> known-output
// case for the DNG / DMM / ALM entry points (native_plugin_trio_dnm3_score.go).
// All values are derived from the math the source implements (log/exp identities
// and the documented model semantics), so the tests run with no upstream binary
// and no reference_code submodule. The libm-tolerance note in the source means
// the long log/exp reductions are validated within an epsilon, not byte-exact.
package bcftools

import (
	"math"
	"testing"
)

// unitEpsilon is the tolerance for the log/exp/lgamma reductions. kfLgamma and
// math.Log differ from a closed-form value only in the last few ULPs, so 1e-9
// is comfortably tight while staying robust across platforms.
const unitEpsilon = 1e-9

func closeTo(t *testing.T, name string, got, want, eps float64) {
	t.Helper()
	if math.Abs(got-want) > eps {
		t.Errorf("%s = %v, want %v (|diff| %g > %g)", name, got, want, math.Abs(got-want), eps)
	}
}

// TestUnitTrioDNM3PhredLog pins the phred/log conversion helpers against their
// closed-form definitions.
func TestUnitTrioDNM3PhredLog(t *testing.T) {
	closeTo(t, "phred2num(0)", phred2num(0), 1.0, unitEpsilon)
	closeTo(t, "phred2num(10)", phred2num(10), 0.1, unitEpsilon)
	closeTo(t, "phred2num(20)", phred2num(20), 0.01, unitEpsilon)
	closeTo(t, "phred2log(0)", phred2log(0), 0.0, unitEpsilon)
	closeTo(t, "phred2log(10)", phred2log(10), -10.0/4.3429, unitEpsilon)
	closeTo(t, "log2phred(0)", log2phred(0), 0.0, unitEpsilon)
	closeTo(t, "log2phred(-2)", log2phred(-2), math.Abs(4.3429*-2), unitEpsilon)
	// phred2num and phred2log are consistent: log(phred2num(p)) == phred2log(p)
	// only approximately (4.3429 vs 10/ln10), so just pin each independently.
}

// TestUnitTrioDNM3SumSubtractLog pins the log-space add/subtract helpers using
// exact log identities: sumLog(0,0)=log2, subtractLog(log3,0)=log2, and the
// -Inf absorbing behavior of sumLog.
func TestUnitTrioDNM3SumSubtractLog(t *testing.T) {
	closeTo(t, "sumLog(0,0)", sumLog(0, 0), math.Log(2), unitEpsilon)
	// sumLog(log a, log b) = log(a+b).
	closeTo(t, "sumLog(log2,log3)", sumLog(math.Log(2), math.Log(3)), math.Log(5), unitEpsilon)
	closeTo(t, "sumLog(log3,log2)", sumLog(math.Log(3), math.Log(2)), math.Log(5), unitEpsilon)
	// -Inf is the log-space zero: sumLog(-Inf, x) = x.
	closeTo(t, "sumLog(-Inf,log4)", sumLog(math.Inf(-1), math.Log(4)), math.Log(4), unitEpsilon)
	if !math.IsInf(sumLog(math.Inf(-1), math.Inf(-1)), -1) {
		t.Error("sumLog(-Inf,-Inf) should be -Inf")
	}
	// subtractLog(log a, log b) = log(a-b).
	closeTo(t, "subtractLog(log3,0)", subtractLog(math.Log(3), 0), math.Log(2), unitEpsilon)
	closeTo(t, "subtractLog(log5,log2)", subtractLog(math.Log(5), math.Log(2)), math.Log(3), unitEpsilon)
}

// TestUnitTrioDNM3MultinomCoeff pins the log multinomial coefficient against the
// exact factorial form: log C(n; cnt) where C is the multinomial coefficient.
func TestUnitTrioDNM3MultinomCoeff(t *testing.T) {
	// [1,1]: 2!/(1!1!) = 2 -> log 2.
	closeTo(t, "logMultinomCoeff[1,1]", logMultinomCoeff([]float64{1, 1}), math.Log(2), unitEpsilon)
	// [2,1]: 3!/(2!1!) = 3 -> log 3.
	closeTo(t, "logMultinomCoeff[2,1]", logMultinomCoeff([]float64{2, 1}), math.Log(3), unitEpsilon)
	// [2,2]: 4!/(2!2!) = 6 -> log 6.
	closeTo(t, "logMultinomCoeff[2,2]", logMultinomCoeff([]float64{2, 2}), math.Log(6), unitEpsilon)
	// A single category collapses to log 1 = 0.
	closeTo(t, "logMultinomCoeff[5]", logMultinomCoeff([]float64{5}), 0.0, unitEpsilon)
}

// TestUnitTrioDNM3DirichletMultinom pins ldirichletMultinom as a proper log-PMF:
// it must be finite for a sensible count/prob vector, and symmetric under
// swapping the two categories when their counts and probabilities are swapped
// together.
func TestUnitTrioDNM3DirichletMultinom(t *testing.T) {
	phi := 100.0
	a := ldirichletMultinom([]float64{2, 1}, []float64{0.5, 0.5}, phi)
	if math.IsNaN(a) || math.IsInf(a, 0) {
		t.Fatalf("ldirichletMultinom returned non-finite %v", a)
	}
	// A log-PMF is <= 0 (probabilities never exceed 1).
	if a > unitEpsilon {
		t.Errorf("ldirichletMultinom log-PMF = %v, want <= 0", a)
	}
	// Symmetry: swapping counts AND probs together leaves the value unchanged.
	b := ldirichletMultinom([]float64{1, 2}, []float64{0.5, 0.5}, phi)
	closeTo(t, "dirichlet symmetry", a, b, unitEpsilon)
}

// TestUnitTrioDNM3ProcessDNG pins the DNG entry point on a tiny fixed trio. A
// child that is strongly heterozygous for an allele neither parent carries is a
// near-certain de-novo: the best de-novo configuration dominates the total, so
// the log score is ~0, with the de-novo allele reported as al1==1.
func TestUnitTrioDNM3ProcessDNG(t *testing.T) {
	priors := newDNMPriorsFull(false, false, 1e-8, autosomalPriors)
	nals := 2 // 3 genotypes: idx0=0/0, idx1=0/1, idx2=1/1

	// Strong de-novo: child het, both parents hom-ref.
	pl := [3][]float64{}
	pl[iFATHER] = []float64{0, -100, -100}
	pl[iMOTHER] = []float64{0, -100, -100}
	pl[iCHILD] = []float64{-100, 0, -100}
	score, al0, al1 := processTrioDNG(priors, nals, pl)
	closeTo(t, "DNG strong de-novo score", score, 0.0, 1e-3)
	if al0 != 0 || al1 != 1 {
		t.Errorf("DNG alleles = (%d,%d), want (0,1)", al0, al1)
	}

	// A clearly inherited het (father het, mother hom-ref) is NOT de novo: the
	// best de-novo config is many log-units below the total, so the score is
	// strongly negative.
	plInh := [3][]float64{}
	plInh[iFATHER] = []float64{-100, 0, -100}
	plInh[iMOTHER] = []float64{0, -100, -100}
	plInh[iCHILD] = []float64{-100, 0, -100}
	inhScore, _, _ := processTrioDNG(priors, nals, plInh)
	if inhScore > -10 {
		t.Errorf("DNG inherited score = %v, want strongly negative (< -10)", inhScore)
	}
	// Determinism: repeated calls are bit-identical.
	again, _, _ := processTrioDNG(priors, nals, plInh)
	if again != inhScore {
		t.Errorf("DNG not deterministic: %v != %v", again, inhScore)
	}
}

// TestUnitTrioDNM3ProcessDMM pins the default Dirichlet-multinomial model on a
// fixed trio. With the noise priors and mosaic term disabled and clean AD
// (parents hom-ref with depth, child a balanced het), the de-novo configuration
// again dominates, giving a score near 0 with the de-novo allele al1==1.
func TestUnitTrioDNM3ProcessDMM(t *testing.T) {
	priors := newDNMPriorsFull(false, false, 1e-8, autosomalPriors)
	nals := 2
	p := &dnmModelParams{
		phi:        200,
		minQM:      0.01,
		minVAF:     0,
		noisePrior: 0,
		withCAD:    false,
		pnCur:      pnoise{}, // zero noise tolerance
	}
	ad := [3][]int{}
	ad[iFATHER] = []int{30, 0}
	ad[iMOTHER] = []int{30, 0}
	ad[iCHILD] = []int{15, 15}
	var qm [3][]float64 // nil -> the model uses |minQM| as the per-read error
	pl := [3][]float64{}
	pl[iFATHER] = []float64{0, -100, -100}
	pl[iMOTHER] = []float64{0, -100, -100}
	pl[iCHILD] = []float64{-100, 0, -100}

	score, al0, al1 := processTrioDMM(p, priors, nals, pl, ad, qm)
	closeTo(t, "DMM strong de-novo score", score, 0.0, 1e-2)
	if al0 != 0 || al1 != 1 {
		t.Errorf("DMM alleles = (%d,%d), want (0,1)", al0, al1)
	}
	again, _, _ := processTrioDMM(p, priors, nals, pl, ad, qm)
	if again != score {
		t.Errorf("DMM not deterministic: %v != %v", again, score)
	}
}

// TestUnitTrioDNM3ProcessALM pins the allele-likelihood model: with noise/dropout
// disabled it must return a finite, non-positive log score deterministically for
// a fixed QS/PL trio. (The exact value is config-dependent; the test pins the
// model's structural guarantees rather than a magic constant.)
func TestUnitTrioDNM3ProcessALM(t *testing.T) {
	priors := newDNMPriorsFull(false, false, 1e-8, autosomalPriors)
	nals := 2
	p := &dnmModelParams{
		phi:         200,
		minQM:       0.01,
		minVAF:      0,
		noisePrior:  0,
		allelicDrop: 0,
		withPPL:     false,
		strictNovel: false,
	}
	ad := [3][]int{}
	ad[iFATHER] = []int{30, 0}
	ad[iMOTHER] = []int{30, 0}
	ad[iCHILD] = []int{15, 15}
	qs := [3][]float64{}
	qs[iFATHER] = []float64{-0.0001, -100}
	qs[iMOTHER] = []float64{-0.0001, -100}
	qs[iCHILD] = []float64{-100, -100}
	pl := [3][]float64{}
	pl[iFATHER] = []float64{0, -100, -100}
	pl[iMOTHER] = []float64{0, -100, -100}
	pl[iCHILD] = []float64{-100, 0, -100}

	score, _, _ := processTrioALM(p, priors, nals, pl, ad, qs, false)
	if math.IsNaN(score) || math.IsInf(score, 0) {
		t.Fatalf("ALM score is non-finite: %v", score)
	}
	if score > 1e-6 {
		t.Errorf("ALM log score = %v, want <= 0", score)
	}
	again, _, _ := processTrioALM(p, priors, nals, pl, ad, qs, false)
	if again != score {
		t.Errorf("ALM not deterministic: %v != %v", again, score)
	}
}
