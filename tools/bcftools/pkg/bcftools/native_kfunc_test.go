package bcftools

import (
	"math"
	"testing"
)

// TestKfBetaiKnownValues checks the in-tree regularized incomplete beta against
// well-known closed-form values. I_x(a,b) with a==b==1 equals x; the symmetric
// midpoint I_{0.5}(a,a)==0.5; and I_x(a,b)+I_{1-x}(b,a)==1.
func TestKfBetaiKnownValues(t *testing.T) {
	cases := []struct {
		a, b, x, want float64
	}{
		{1, 1, 0.25, 0.25},
		{1, 1, 0.5, 0.5},
		{1, 1, 0.9, 0.9},
		{2, 2, 0.5, 0.5},
		{5, 5, 0.5, 0.5},
		{2, 3, 0.5, 0.6875}, // I_{0.5}(2,3) = 11/16
	}
	for _, c := range cases {
		got := kfBetai(c.a, c.b, c.x)
		if math.Abs(got-c.want) > 1e-12 {
			t.Errorf("kfBetai(%g,%g,%g)=%.15g, want %.15g", c.a, c.b, c.x, got, c.want)
		}
	}
	// Reflection identity.
	for _, x := range []float64{0.1, 0.3, 0.5, 0.7, 0.9} {
		s := kfBetai(3, 5, x) + kfBetai(5, 3, 1-x)
		if math.Abs(s-1) > 1e-12 {
			t.Errorf("reflection identity failed at x=%g: sum=%.15g", x, s)
		}
	}
}

// TestCalcBinomTwoSided checks the two-sided binomial tail's documented edge
// cases (both zero -> -1; equal counts -> 1; symmetric fair split -> 1) and a
// strongly skewed case (<0.05).
func TestCalcBinomTwoSided(t *testing.T) {
	if got := calcBinomTwoSided(0, 0, 0.5); got != -1 {
		t.Errorf("calcBinomTwoSided(0,0)=%g, want -1", got)
	}
	if got := calcBinomTwoSided(5, 5, 0.5); got != 1 {
		t.Errorf("calcBinomTwoSided(5,5)=%g, want 1", got)
	}
	if got := calcBinomTwoSided(1, 1, 0.5); got != 1 {
		t.Errorf("calcBinomTwoSided(1,1)=%g, want 1", got)
	}
	// 18 vs 2 under p=0.5 is strongly skewed: two-sided p well below 0.05.
	if got := calcBinomTwoSided(18, 2, 0.5); got <= 0 || got >= 0.05 {
		t.Errorf("calcBinomTwoSided(18,2,0.5)=%g, want a small positive tail (<0.05)", got)
	}
	// Result is clamped to <=1.
	if got := calcBinomTwoSided(2, 1, 0.5); got > 1 {
		t.Errorf("calcBinomTwoSided result not clamped: %g", got)
	}
}

// TestCalcBinomOneSided checks the one-sided tail equals the regularized
// incomplete-beta closed forms it wraps (ge -> I_aprob(na,nb+1); !ge ->
// I_{1-aprob}(nb,na+1)) and stays within [0,1].
func TestCalcBinomOneSided(t *testing.T) {
	cases := []struct {
		na, nb int
		aprob  float64
		ge     bool
	}{
		{11, 9, 1.0 / 3, true},
		{11, 9, 2.0 / 3, false},
		{5, 5, 0.5, true},
		{5, 5, 0.5, false},
		{2, 18, 1.0 / 3, true},
	}
	for _, c := range cases {
		got := calcBinomOneSided(c.na, c.nb, c.aprob, c.ge)
		var want float64
		if c.ge {
			want = kfBetai(float64(c.na), float64(c.nb)+1, c.aprob)
		} else {
			want = kfBetai(float64(c.nb), float64(c.na)+1, 1-c.aprob)
		}
		if got != want {
			t.Errorf("calcBinomOneSided(%d,%d,%g,%v)=%.15g, want %.15g", c.na, c.nb, c.aprob, c.ge, got, want)
		}
		if got < 0 || got > 1 {
			t.Errorf("one-sided tail out of [0,1]: %g", got)
		}
	}
}

// TestDNMPriorsTrivialMendel checks the NAIVE de-novo priors on hand-verified
// trio genotypes: a clear de-novo (both parents hom-ref, child het) is flagged,
// and a normal transmission (parents het, child het) is not.
func TestDNMPriorsTrivialMendel(t *testing.T) {
	p := newDNMPriors(false, autosomalPriors)
	// Genotype indices via the seq mapping: mask (1<<a)|(1<<b) -> dnmSeq3.
	gt := func(a, b int) int { return dnmSeq3[(1<<a)|(1<<b)] }
	homRef := gt(0, 0) // 0/0
	het := gt(0, 1)    // 0/1
	homAlt := gt(1, 1) // 1/1

	// 0/0 x 0/0 -> 0/1 child: de novo, allele 1.
	if p.denovo[homRef][homRef][het] != 1 {
		t.Errorf("expected de-novo for 0/0 x 0/0 -> 0/1")
	}
	if p.denovoAllele[homRef][homRef][het] != 1 {
		t.Errorf("expected de-novo allele 1, got %d", p.denovoAllele[homRef][homRef][het])
	}
	// 0/1 x 0/1 -> 0/1 child: normal transmission, not de novo.
	if p.denovo[het][het][het] != 0 {
		t.Errorf("expected NOT de-novo for 0/1 x 0/1 -> 0/1")
	}
	// 0/0 x 1/1 -> 0/1 child: normal transmission, not de novo.
	if p.denovo[homRef][homAlt][het] != 0 {
		t.Errorf("expected NOT de-novo for 0/0 x 1/1 -> 0/1")
	}
	// 0/0 x 0/0 -> 1/1 child: de novo (two alleles mutated).
	if p.denovo[homRef][homRef][homAlt] != 1 {
		t.Errorf("expected de-novo for 0/0 x 0/0 -> 1/1")
	}
}

// TestChrXMatcher checks the default GRCh37 PAR-exclusion list and the GRCh38
// shortcut place positions in/out of the non-PAR chrX region.
func TestChrXMatcher(t *testing.T) {
	m := buildChrXMatcher("")
	if !m.overlaps("X", 100) {
		t.Errorf("expected X:100 inside the non-PAR region (PAR1 1-60000)")
	}
	if m.overlaps("X", 1000000) {
		t.Errorf("expected X:1000000 in the PAR gap (60000-2699521), not matched")
	}
	if !m.overlaps("X", 5000000) {
		t.Errorf("expected X:5000000 inside the non-PAR region")
	}
	if m.overlaps("1", 100) {
		t.Errorf("autosome must never match the chrX list")
	}
	m38 := buildChrXMatcher("GRCh38")
	if !m38.overlaps("chrX", 100) {
		t.Errorf("expected chrX:100 inside the GRCh38 non-PAR region (1-9999)")
	}
}
