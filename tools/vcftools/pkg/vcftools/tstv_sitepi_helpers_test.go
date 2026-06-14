package vcftools

import "testing"

// TestSiteIsDiploid exercises the upstream entry::is_diploid gate used by
// --site-pi: a site is diploid only when every non-empty GT has exactly two
// allele slots; a fully-missing "./." still counts as ploidy 2.
func TestSiteIsDiploid(t *testing.T) {
	cases := []struct {
		name string
		gt   []string
		want bool
	}{
		{"all diploid", []string{"0/0", "0|1", "1/1"}, true},
		{"missing diploid kept", []string{"./.", "0/1", "1|0"}, true},
		{"haploid present", []string{"0", "0/1", "1/1"}, false},
		{"haploid missing", []string{".", "0/0"}, false},
		{"empty GT ignored", []string{"", "0/0", "1/1"}, true},
		{"triploid", []string{"0/0/1", "0/0"}, false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			v := createTestVariant("1", 100, "A", []string{"G"}, 50, tc.gt)
			if got := siteIsDiploid(v); got != tc.want {
				t.Fatalf("siteIsDiploid(%v) = %v, want %v", tc.gt, got, tc.want)
			}
		})
	}
}

// TestIsNucleotide checks the single-base allele predicate that excludes ".",
// "N", and symbolic ALTs from the Ts/Tv tallies.
func TestIsNucleotide(t *testing.T) {
	cases := []struct {
		s    string
		want bool
	}{
		{"A", true}, {"C", true}, {"G", true}, {"T", true},
		{".", false}, {"N", false}, {"a", false}, {"AC", false},
		{"", false}, {"<DEL>", false},
	}
	for _, tc := range cases {
		if got := isNucleotide(tc.s); got != tc.want {
			t.Errorf("isNucleotide(%q) = %v, want %v", tc.s, got, tc.want)
		}
	}
}

// TestCppRatio pins the C++-ostream Ts/Tv ratio formatting: "-nan" for 0/0,
// "inf"/"-inf" for a finite-over-zero, and C++ defaultfloat for finite ratios.
func TestCppRatio(t *testing.T) {
	cases := []struct {
		num, den float64
		want     string
	}{
		{0, 0, "-nan"},
		{1, 0, "inf"},
		{-1, 0, "-inf"},
		{0, 1, "0"},
		{1, 2, "0.5"},
		{1, 1, "1"},
		{2, 3, "0.666667"},
	}
	for _, tc := range cases {
		if got := cppRatio(tc.num, tc.den); got != tc.want {
			t.Errorf("cppRatio(%v, %v) = %q, want %q", tc.num, tc.den, got, tc.want)
		}
	}
}
