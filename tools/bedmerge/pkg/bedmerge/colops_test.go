package bedmerge

import (
	"strings"
	"testing"
)

// TestApplyOp_NewOps covers the column ops added for full-parity closure
// with upstream bedtools: stdev, sstdev, absmin, absmax, cat, cat_uniq.
// Expected values are hand-computed against a fixed 4-element input.
func TestApplyOp_NewOps(t *testing.T) {
	// Numeric ops fixture: {2, 4, 4, 6}. mean=4. ssq=(4+0+0+4)=8.
	// stdev = sqrt(8/4)  = sqrt(2)   ≈ 1.4142135623730951
	// sstdev = sqrt(8/3) = sqrt(2.6) ≈ 1.632993161855452
	numVals := []string{"2", "4", "4", "6"}

	cases := []struct {
		op       string
		vals     []string
		wantPfx  string // we accept any %g-formatted prefix of this length
		wantFull string // if non-empty, must match exactly
	}{
		{op: "stdev", vals: numVals, wantPfx: "1.4142135"},
		{op: "sstdev", vals: numVals, wantPfx: "1.6329931"},
		// Sign drop: {-3, 1, 2, -4} -> |.| = {3,1,2,4}; absmin=1, absmax=4.
		{op: "absmin", vals: []string{"-3", "1", "2", "-4"}, wantFull: "1"},
		{op: "absmax", vals: []string{"-3", "1", "2", "-4"}, wantFull: "4"},
		// cat: concat with no separator.
		{op: "cat", vals: []string{"a", "b", "a", "c"}, wantFull: "abac"},
		// cat_uniq: unique in first-appearance order, no separator.
		{op: "cat_uniq", vals: []string{"a", "b", "a", "c"}, wantFull: "abc"},
	}
	for _, c := range cases {
		got, err := ApplyOp(c.op, 4, c.vals)
		if err != nil {
			t.Errorf("%s: unexpected error: %v", c.op, err)
			continue
		}
		if c.wantFull != "" {
			if got != c.wantFull {
				t.Errorf("%s: got %q want %q", c.op, got, c.wantFull)
			}
			continue
		}
		if !strings.HasPrefix(got, c.wantPfx) {
			t.Errorf("%s: got %q, want prefix %q", c.op, got, c.wantPfx)
		}
	}
}

// TestApplyOp_SstdevSingle confirms that sample stdev on a single value
// returns "." (upstream prints its NullValue on NaN; with n=1 the divisor
// is zero so sqrt(0/0) = NaN).
func TestApplyOp_SstdevSingle(t *testing.T) {
	got, err := ApplyOp("sstdev", 4, []string{"5"})
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if got != "." {
		t.Errorf("sstdev[n=1]: got %q want %q", got, ".")
	}
}

// TestApplyOp_StdevSingle confirms population stdev on a single value is 0
// (the population variance of one point is 0).
func TestApplyOp_StdevSingle(t *testing.T) {
	got, err := ApplyOp("stdev", 4, []string{"5"})
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if got != "0" {
		t.Errorf("stdev[n=1]: got %q want %q", got, "0")
	}
}

// TestApplyOp_AbsMinSigned exercises absmin on already-positive values: the
// answer is just the regular min.
func TestApplyOp_AbsMinSigned(t *testing.T) {
	got, err := ApplyOp("absmin", 4, []string{"5", "3", "7"})
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if got != "3" {
		t.Errorf("absmin: got %q want %q", got, "3")
	}
}

// TestApplyOp_CatEmpty verifies cat / cat_uniq on a single value just echo it.
func TestApplyOp_CatEmpty(t *testing.T) {
	for _, op := range []string{"cat", "cat_uniq"} {
		got, err := ApplyOp(op, 4, []string{"only"})
		if err != nil {
			t.Errorf("%s: unexpected error: %v", op, err)
			continue
		}
		if got != "only" {
			t.Errorf("%s[singleton]: got %q want %q", op, got, "only")
		}
	}
}

// TestParseColumnOps_NewOps confirms the new op names are accepted by the
// flag parser.
func TestParseColumnOps_NewOps(t *testing.T) {
	for _, op := range []string{"stdev", "sstdev", "absmin", "absmax", "cat", "cat_uniq"} {
		if _, err := ParseColumnOps("5", op); err != nil {
			t.Errorf("ParseColumnOps(5, %q) failed: %v", op, err)
		}
	}
}
