package bcftools

import "testing"

// Unit tests for the proximity (tolerance-aware) parity comparison helper.

func TestProximityNumericClose(t *testing.T) {
	tol := defaultProximityTolerance
	cases := []struct {
		a, b float64
		want bool
	}{
		{-46.0521, -46.0522, true}, // last %g digit differs -> within 6 sig figs
		{-46.0521, -46.0521, true},
		{1.0, 1.0000001, true}, // within rel eps
		{0.0, 1e-7, true},      // within abs eps
		{-46.05, -46.5, false}, // genuinely different score
		{100.0, 101.0, false},  // off by a whole point
		{1e-300, 2e-300, true}, // both effectively zero; within abs eps
		{2.0, 3.0, false},      // 2x apart, well above abs eps
	}
	for _, c := range cases {
		if got := numericClose(c.a, c.b, tol); got != c.want {
			t.Errorf("numericClose(%g,%g)=%v, want %v", c.a, c.b, got, c.want)
		}
	}
}

func TestProximityNaNInf(t *testing.T) {
	tol := defaultProximityTolerance
	nan, ninf := parseNumericFieldMust(t, "nan"), parseNumericFieldMust(t, "-nan")
	if !numericClose(nan, ninf, tol) {
		t.Errorf("nan and -nan should compare equal")
	}
	pinf := parseNumericFieldMust(t, "inf")
	minf := parseNumericFieldMust(t, "-inf")
	if numericClose(pinf, minf, tol) {
		t.Errorf("+inf and -inf must differ in sign")
	}
	if !numericClose(pinf, pinf, tol) {
		t.Errorf("+inf should equal +inf")
	}
	if numericClose(nan, 1.0, tol) {
		t.Errorf("nan must not equal a finite number")
	}
}

func parseNumericFieldMust(t *testing.T, s string) float64 {
	t.Helper()
	f, ok := parseNumericField(s)
	if !ok {
		t.Fatalf("parseNumericField(%q) failed", s)
	}
	return f
}

func TestProximityFieldAware(t *testing.T) {
	// String fields must match exactly; numeric sub-fields within a composite
	// VCF column may differ in the last printed digit.
	want := "chr1\t100\t.\tA\tT\t50\tPASS\tAC=4\tGT:DNM:VA\t0|1:-46.0521:0"
	got := "chr1\t100\t.\tA\tT\t50\tPASS\tAC=4\tGT:DNM:VA\t0|1:-46.0522:0"
	if d := compareProximityDefault(want, got); d != nil {
		t.Errorf("expected proximity-equal, got diffs:\n%v", d)
	}

	// A string field difference (REF A->C) must be caught.
	got2 := "chr1\t100\t.\tC\tT\t50\tPASS\tAC=4\tGT:DNM:VA\t0|1:-46.0521:0"
	d := compareProximityDefault(want, got2)
	if d == nil {
		t.Fatalf("expected a diff for differing REF allele")
	}
	if d[0].Field != 3 {
		t.Errorf("expected mismatch at field 3 (REF), got field %d: %v", d[0].Field, d[0])
	}

	// A real numeric divergence (score off by a whole point) must be caught.
	got3 := "chr1\t100\t.\tA\tT\t50\tPASS\tAC=4\tGT:DNM:VA\t0|1:-47.0521:0"
	if compareProximityDefault(want, got3) == nil {
		t.Fatalf("expected a diff for a score off by a whole point")
	}
}

func TestProximityLineFieldCount(t *testing.T) {
	if compareProximityDefault("a\tb", "a\tb\tc") == nil {
		t.Errorf("expected a diff for differing field count")
	}
	if compareProximityDefault("a\nb", "a") == nil {
		t.Errorf("expected a diff for differing line count")
	}
	// Trailing newline must not register as an extra line.
	if d := compareProximityDefault("a\tb\n", "a\tb"); d != nil {
		t.Errorf("trailing newline should be ignored, got: %v", d)
	}
}

func TestRoundSig(t *testing.T) {
	cases := []struct {
		x    float64
		n    int
		want float64
	}{
		{-46.05213, 6, -46.0521},
		{123456.7, 6, 123457},
		{0.00012345678, 6, 0.000123457},
		{0, 6, 0},
	}
	for _, c := range cases {
		if got := roundSig(c.x, c.n); got != c.want {
			t.Errorf("roundSig(%g,%d)=%g, want %g", c.x, c.n, got, c.want)
		}
	}
}
