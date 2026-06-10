package main

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

// TestCheckFilterDeferred locks in the upstream-flag-name surface
// that runFilter hard-rejects rather than silently accepting. Per
// the project parity rule (docs/PARITY_ROADMAP.md "Definition of 1:1")
// every documented upstream flag must be recognised — either
// implemented or gracefully rejected with a roadmap pointer. A
// regression that drops any of these from the rejection set without
// implementing the underlying behaviour is a parity bug.
func TestCheckFilterDeferred(t *testing.T) {
	if got := checkFilterDeferred(checkFilterDeferredInputs{}); got != "" {
		t.Fatalf("empty input: got deferred=%q, want \"\"", got)
	}
	cases := []struct {
		name string
		in   checkFilterDeferredInputs
		want string
	}{
		{"mask-region", checkFilterDeferredInputs{maskRegion: "chr1:100-200"}, "--mask"},
		{"mask-file", checkFilterDeferredInputs{maskFile: "mask.bed"}, "-M/--mask-file"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkFilterDeferred(tc.in); got != tc.want {
				t.Errorf("deferred(%s): got %q, want %q", tc.name, got, tc.want)
			}
		})
	}
}

// TestParseSnpGap accepts the upstream "INT[:TYPE]" form.
func TestParseSnpGap(t *testing.T) {
	cases := []struct {
		in   string
		want int
		err  bool
	}{
		{"", 0, false},
		{"5", 5, false},
		{"10:indel", 10, false},
		{"3:mnp,bnd,overlap", 3, false},
		{"abc", 0, true},
	}
	for _, tc := range cases {
		got, err := parseSnpGap(tc.in)
		if tc.err && err == nil {
			t.Errorf("parseSnpGap(%q) expected error, got %d", tc.in, got)
		}
		if !tc.err && err != nil {
			t.Errorf("parseSnpGap(%q) unexpected error: %v", tc.in, err)
		}
		if !tc.err && got != tc.want {
			t.Errorf("parseSnpGap(%q) = %d, want %d", tc.in, got, tc.want)
		}
	}
}

// TestParseHaplotypeNpIu checks the "NpIu" (phased-index / unphased-IUPAC)
// haplotype encoding parses into HapPhasedIUPAC with the right index, and
// that malformed variants are rejected.
func TestParseHaplotypeNpIu(t *testing.T) {
	cases := []struct {
		in      string
		wantSel bcftools.HaplotypeSelector
		wantIdx int
		wantErr bool
	}{
		{"2pIu", bcftools.HapPhasedIUPAC, 2, false},
		{"1pIu", bcftools.HapPhasedIUPAC, 1, false},
		{"10pIu", bcftools.HapPhasedIUPAC, 10, false},
		{"2PIU", bcftools.HapPhasedIUPAC, 2, false}, // case-insensitive suffix
		{"2", bcftools.HapIndex, 2, false},
		{"I", bcftools.HapIUPAC, 0, false},
		{"2pIuX", 0, 0, true},
		{"pIu", 0, 0, true},
		{"apIu", 0, 0, true},
	}
	for _, tc := range cases {
		sel, idx, err := bcftools.ParseHaplotypeSelector(tc.in)
		if tc.wantErr {
			if err == nil {
				t.Errorf("ParseHaplotypeSelector(%q): expected error, got sel=%v idx=%d", tc.in, sel, idx)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseHaplotypeSelector(%q): unexpected error %v", tc.in, err)
			continue
		}
		if sel != tc.wantSel || idx != tc.wantIdx {
			t.Errorf("ParseHaplotypeSelector(%q) = (%v,%d), want (%v,%d)", tc.in, sel, idx, tc.wantSel, tc.wantIdx)
		}
	}
}
