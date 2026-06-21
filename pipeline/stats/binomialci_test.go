package stats

import (
	"math"
	"testing"
)

// z95 is the two-sided 95% standard-normal quantile.
const z95 = 1.959963984540054

func approx(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Errorf("%s = %.6f, want %.6f (tol %.1g)", name, got, want, tol)
	}
}

// TestWilsonCI checks the Wilson score interval against values reproduced from
// the closed-form formula (and matching R's binom::binom.wilson and the
// Newcombe 1998 worked examples).
func TestWilsonCI(t *testing.T) {
	tests := []struct {
		name           string
		k, n           int
		wantLo, wantHi float64
		tol            float64
	}{
		// 398/400 = 0.995. Wilson 95% ≈ [0.9821, 0.9986].
		{"398of400", 398, 400, 0.982054, 0.998614, 1e-4},
		// 89/89 = 1.0 all-success: Wilson lower bound is non-degenerate.
		// 95% ≈ [0.9586, 1.0].
		{"89of89", 89, 89, 0.958612, 1.0, 1e-4},
		// Textbook midpoint sanity: 50/100. 95% ≈ [0.4038, 0.5962].
		{"50of100", 50, 100, 0.403820, 0.596180, 1e-4},
		// Small all-success 10/10: 95% ≈ [0.7225, 1.0].
		{"10of10", 10, 10, 0.722460, 1.0, 1e-4},
		// 223/314 ≈ 0.7102 (project audit pass rate). 95% ≈ [0.6573, 0.7581].
		{"223of314", 223, 314, 0.657290, 0.758101, 1e-3},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			lo, hi := WilsonCI(tc.k, tc.n, z95)
			approx(t, "lo", lo, tc.wantLo, tc.tol)
			approx(t, "hi", hi, tc.wantHi, tc.tol)
			if lo < 0 || hi > 1 || lo > hi {
				t.Errorf("interval out of order or out of [0,1]: [%f, %f]", lo, hi)
			}
		})
	}
}

// TestClopperPearsonCI checks the exact interval against values reproduced from
// R's binom.test (the canonical Clopper-Pearson reference).
func TestClopperPearsonCI(t *testing.T) {
	tests := []struct {
		name           string
		k, n           int
		wantLo, wantHi float64
		tol            float64
	}{
		// binom.test(398, 400)$conf.int -> [0.98203, 0.99940].
		{"398of400", 398, 400, 0.982029, 0.999395, 1e-3},
		// binom.test(89, 89)$conf.int -> [0.95937, 1.0]. All-success: exact
		// lower bound = (alpha/2)^(1/n).
		{"89of89", 89, 89, 0.959374, 1.0, 1e-3},
		// binom.test(50, 100)$conf.int -> [0.39833, 0.60167].
		{"50of100", 50, 100, 0.398328, 0.601672, 1e-3},
		// binom.test(10, 10)$conf.int -> [0.69150, 1.0].
		{"10of10", 10, 10, 0.691504, 1.0, 1e-3},
		// Zero-success boundary: binom.test(0, 20) -> [0.0, 0.16843].
		{"0of20", 0, 20, 0.0, 0.168157, 1e-3},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			lo, hi := ClopperPearsonCI(tc.k, tc.n, 0.05)
			approx(t, "lo", lo, tc.wantLo, tc.tol)
			approx(t, "hi", hi, tc.wantHi, tc.tol)
			if lo < 0 || hi > 1 || lo > hi {
				t.Errorf("interval out of order or out of [0,1]: [%f, %f]", lo, hi)
			}
		})
	}
}

// TestAllSuccessLowerBoundFormula verifies the exact all-success lower bound
// matches the closed form (alpha/2)^(1/n), an independent check on the Beta
// quantile inversion at the boundary.
func TestAllSuccessLowerBoundFormula(t *testing.T) {
	for _, n := range []int{5, 50, 89, 400} {
		lo, hi := ClopperPearsonCI(n, n, 0.05)
		want := math.Pow(0.025, 1.0/float64(n))
		approx(t, "lo", lo, want, 1e-6)
		if hi != 1.0 {
			t.Errorf("n=%d: upper bound = %f, want exactly 1.0", n, hi)
		}
	}
}

// TestZeroTrials documents the no-data convention: the whole unit interval.
func TestZeroTrials(t *testing.T) {
	if lo, hi := WilsonCI(0, 0, z95); lo != 0 || hi != 1 {
		t.Errorf("WilsonCI(0,0) = [%f,%f], want [0,1]", lo, hi)
	}
	if lo, hi := ClopperPearsonCI(0, 0, 0.05); lo != 0 || hi != 1 {
		t.Errorf("ClopperPearsonCI(0,0) = [%f,%f], want [0,1]", lo, hi)
	}
}

// TestIncompleteBetaIdentities sanity-checks the incomplete Beta engine that
// backs the Clopper-Pearson interval: I_0=0, I_1=1, the symmetric midpoint of a
// symmetric Beta is 1/2, and the reflection identity holds.
func TestIncompleteBetaIdentities(t *testing.T) {
	approx(t, "I_0(2,3)", incompleteBeta(0, 2, 3), 0, 1e-12)
	approx(t, "I_1(2,3)", incompleteBeta(1, 2, 3), 1, 1e-12)
	approx(t, "I_.5(3,3)", incompleteBeta(0.5, 3, 3), 0.5, 1e-9)
	// Reflection: I_x(a,b) = 1 - I_{1-x}(b,a).
	approx(t, "reflection", incompleteBeta(0.3, 2, 5), 1-incompleteBeta(0.7, 5, 2), 1e-9)
	// Round-trip: quantile then CDF recovers p.
	p := 0.025
	x := betaQuantile(p, 7.0, 84.0)
	approx(t, "roundtrip", incompleteBeta(x, 7, 84), p, 1e-6)
}
