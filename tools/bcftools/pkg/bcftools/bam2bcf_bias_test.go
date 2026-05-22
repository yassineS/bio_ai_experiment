package bcftools

import (
	"math"
	"testing"
)

// TestKfErfc checks the ported htslib kf_erfc against a few reference
// points. kf_erfc(x) evaluates erfc(x*sqrt2), so kf_erfc(0)==1.
func TestKfErfc(t *testing.T) {
	cases := []struct {
		x, want float64
	}{
		{0, 1},
		{37.1, 0},  // z>37 short-circuit, positive
		{-37.1, 2}, // z>37 short-circuit, negative
	}
	for _, c := range cases {
		got := kfErfc(c.x)
		if math.Abs(got-c.want) > 1e-9 {
			t.Errorf("kfErfc(%v) = %v, want %v", c.x, got, c.want)
		}
	}
	// Symmetry: kf_erfc(x) + kf_erfc(-x) == 2.
	for _, x := range []float64{0.1, 0.7, 1.5, 3.0} {
		if s := kfErfc(x) + kfErfc(-x); math.Abs(s-2) > 1e-9 {
			t.Errorf("kfErfc(%v)+kfErfc(-%v) = %v, want 2", x, x, s)
		}
	}
}

// TestCalcVDB covers the insufficient-depth guard, the dp==2 exact
// branch and a deeper distribution.
func TestCalcVDB(t *testing.T) {
	// Fewer than two variant reads: VDB is undefined.
	if _, ok := calcVDB(make([]int, b2bNpos)); ok {
		t.Error("calcVDB with no reads should be undefined")
	}
	one := make([]int, b2bNpos)
	one[10] = 1
	if _, ok := calcVDB(one); ok {
		t.Error("calcVDB with one read should be undefined")
	}
	// dp==2: two reads at distinct positions -> a defined value in (0,1].
	two := make([]int, b2bNpos)
	two[10], two[40] = 1, 1
	v, ok := calcVDB(two)
	if !ok {
		t.Fatal("calcVDB with two reads should be defined")
	}
	if v <= 0 || v > 1 {
		t.Errorf("calcVDB dp==2 = %v, want in (0,1]", v)
	}
	// Deeper, well-spread distribution -> defined, in [0,1].
	deep := make([]int, b2bNpos)
	for i := 20; i < 60; i++ {
		deep[i] = 1
	}
	v, ok = calcVDB(deep)
	if !ok || v < 0 || v > 1.0001 {
		t.Errorf("calcVDB deep = %v ok=%v, want defined in [0,1]", v, ok)
	}
}

// TestCalcMWUBiasZ covers the empty-histogram guard and the symmetry of
// identical distributions (z-score ~ 0).
func TestCalcMWUBiasZ(t *testing.T) {
	a := make([]int, b2bNqual)
	b := make([]int, b2bNqual)
	// One side empty -> undefined.
	a[10] = 5
	if _, ok := calcMWUBiasZ(a, b, false); ok {
		t.Error("calcMWUBiasZ with empty b should be undefined")
	}
	// Identical distributions -> z-score 0.
	for i := 10; i < 20; i++ {
		a[i], b[i] = 3, 3
	}
	a[10] = 3 // overwrite the lone spike so a and b match exactly
	z, ok := calcMWUBiasZ(a, b, false)
	if !ok {
		t.Fatal("calcMWUBiasZ identical should be defined")
	}
	if math.Abs(z) > 1e-9 {
		t.Errorf("calcMWUBiasZ identical = %v, want ~0", z)
	}
	// A shifted distribution yields a non-zero z-score.
	c := make([]int, b2bNqual)
	d := make([]int, b2bNqual)
	for i := 0; i < 10; i++ {
		c[i] = 4
	}
	for i := 40; i < 50; i++ {
		d[i] = 4
	}
	z, ok = calcMWUBiasZ(c, d, false)
	if !ok || math.Abs(z) < 1 {
		t.Errorf("calcMWUBiasZ shifted = %v ok=%v, want a sizable z-score", z, ok)
	}
}

// TestCalcSegBias covers the no-variant-read guard and a basic run.
func TestCalcSegBias(t *testing.T) {
	var call bcfCall
	// No non-reference reads -> undefined.
	calls := []bcfCallret{{}}
	if _, ok := calcSegBias(calls, &call); ok {
		t.Error("calcSegBias with no non-ref reads should be undefined")
	}
	// Some non-reference reads -> defined finite value.
	call.anno[0], call.anno[1] = 10, 10 // ref depth
	call.anno[2], call.anno[3] = 3, 2   // non-ref depth
	calls[0].anno[2], calls[0].anno[3] = 3, 2
	v, ok := calcSegBias(calls, &call)
	if !ok {
		t.Fatal("calcSegBias should be defined")
	}
	if math.IsInf(v, 0) || math.IsNaN(v) {
		t.Errorf("calcSegBias = %v, want finite", v)
	}
}

// TestGetPosition checks the soft-clip-aware read-position annotation.
func TestGetPosition(t *testing.T) {
	// Plain 100M read, base at qpos 49: no soft-clip, pos = 50.
	r := getPosition([]int{0}, []int{100}, 49, 100)
	if r.pos != 50 || r.length != 100 || r.scLen != 0 {
		t.Errorf("getPosition 100M = %+v, want pos=50 length=100 scLen=0", r)
	}
	// 10S90M: 10 bp leading soft-clip, base at qpos 14 (5th aligned
	// base). edist = 15 - 10 = 5; length = 100-10 = 90.
	r = getPosition([]int{4, 0}, []int{10, 90}, 14, 100)
	if r.pos != 5 || r.length != 90 || r.scLen != 10 {
		t.Errorf("getPosition 10S90M = %+v, want pos=5 length=90 scLen=10", r)
	}
	// 90M10S: trailing soft-clip; aligned base near the start.
	r = getPosition([]int{0, 4}, []int{90, 10}, 5, 100)
	if r.length != 90 || r.scLen != 10 {
		t.Errorf("getPosition 90M10S = %+v, want length=90 scLen=10", r)
	}
}

// TestWangHash pins the read-name hash machinery so the smart-overlap
// mate selection stays deterministic.
func TestWangHash(t *testing.T) {
	// Two different names must (almost surely) hash differently.
	if x31HashString("readA") == x31HashString("readB") {
		t.Error("x31HashString collided on readA/readB")
	}
	// wangHash is a bijection on uint32; distinct inputs differ.
	if wangHash(1) == wangHash(2) {
		t.Error("wangHash collided on 1/2")
	}
	// Empty string hashes to 0 (matches the C guard).
	if x31HashString("") != 0 {
		t.Error("x31HashString(\"\") should be 0")
	}
}
