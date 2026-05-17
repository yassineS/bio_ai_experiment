package vcftools

// Unit and parity tests for --diff-switch-error.
//
// Goldens generated from upstream vcftools (reference_code submodule) with the
// documented FORTIFY_SOURCE workaround. See diff_parity_test.go's header for
// the make incantation.

import (
	"path/filepath"
	"testing"
)

// TestParity_DiffSwitchError exercises --diff-switch-error end-to-end against
// upstream byte-for-byte. The two fixtures share sample names s1/s2 (no
// --diff-indv-map needed) and contain one switch event for s1 between sites
// 100→200, another for s1 between sites 300→400, and one for s2 between
// 300→400.
func TestParity_DiffSwitchError(t *testing.T) {
	fix := vcftoolsFixtureDir(t)
	prefix := runVcftoolsParity(t, "switch_f1.vcf", &Params{
		Diff:            filepath.Join(fix, "switch_f2.vcf"),
		DiffSwitchError: true,
	})
	gotEvents := readFileBytes(t, prefix+".diff.switch")
	wantEvents := readFileBytes(t, filepath.Join(fix, "switch.expected.diff.switch"))
	if string(gotEvents) != string(wantEvents) {
		t.Errorf(".diff.switch byte mismatch:\n--- got ---\n%s\n--- want ---\n%s",
			string(gotEvents), string(wantEvents))
	}
	gotIndv := readFileBytes(t, prefix+".diff.indv.switch")
	wantIndv := readFileBytes(t, filepath.Join(fix, "switch.expected.diff.indv.switch"))
	if string(gotIndv) != string(wantIndv) {
		t.Errorf(".diff.indv.switch byte mismatch:\n--- got ---\n%s\n--- want ---\n%s",
			string(gotIndv), string(wantIndv))
	}
}

// TestSplitPhasedAlleles exercises the GT parser used by diff_switch.go.
// Diff-switch's hot path only cares about diploid phased calls with two
// non-missing alleles; haploid and missing inputs must return phased=false /
// empty allele strings so the caller drops them.
func TestSplitPhasedAlleles(t *testing.T) {
	cases := []struct {
		gt string
		a1 string
		a2 string
		ph bool
	}{
		{"0|1", "0", "1", true},
		{"1|0", "1", "0", true},
		{"0/1", "0", "1", false},
		{"./.", ".", ".", false},
		{"0", "", "", false}, // haploid
		{"", "", "", false},
		{".|.", ".", ".", true},
	}
	for _, c := range cases {
		gotA1, gotA2, gotPh := splitPhasedAlleles(c.gt)
		if gotA1 != c.a1 || gotA2 != c.a2 || gotPh != c.ph {
			t.Errorf("splitPhasedAlleles(%q) = (%q, %q, %v); want (%q, %q, %v)",
				c.gt, gotA1, gotA2, gotPh, c.a1, c.a2, c.ph)
		}
	}
}

// TestFormatSwitchRate confirms our switch-rate formatter matches upstream's
// C++ default `cout << double` precision (6 significant digits, trimmed).
func TestFormatSwitchRate(t *testing.T) {
	cases := []struct {
		v    float64
		want string
	}{
		{0, "0"},
		{0.5, "0.5"},
		{0.25, "0.25"},
		{1.0 / 3.0, "0.333333"},
		{1, "1"},
	}
	for _, c := range cases {
		got := formatSwitchRate(c.v)
		if got != c.want {
			t.Errorf("formatSwitchRate(%v) = %q; want %q", c.v, got, c.want)
		}
	}
}
