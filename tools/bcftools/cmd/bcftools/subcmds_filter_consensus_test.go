package main

import "testing"

// TestCheckFilterDeferred is now an empty-set sentinel: every
// documented filter flag is either implemented or gracefully
// accepted-and-ignored. The previous --mask / -M deferred branch
// landed when filter_mask.go shipped, so the set is empty. If a
// future regression removes a flag from runFilter, add a test here.
func TestCheckFilterDeferred(t *testing.T) {
	if got := checkFilterDeferred(checkFilterDeferredInputs{}); got != "" {
		t.Errorf("empty input: got deferred=%q, want \"\"", got)
	}
}

// TestCheckConsensusDeferred is the corresponding test for consensus.
func TestCheckConsensusDeferred(t *testing.T) {
	if got := checkConsensusDeferred(checkConsensusDeferredInputs{}); got != "" {
		t.Fatalf("empty input: got deferred=%q, want \"\"", got)
	}
	cases := []struct {
		name string
		in   checkConsensusDeferredInputs
		want string
	}{
		{"chain", checkConsensusDeferredInputs{chainFile: "out.chain"}, "-c/--chain (liftover chain output)"},
		{"NpIu", checkConsensusDeferredInputs{haplotype: "2pIu"}, "-H NpIu (phased-index / unphased-IUPAC)"},
		// Other haplotype codes must not be flagged as deferred.
		{"R", checkConsensusDeferredInputs{haplotype: "R"}, ""},
		{"LA", checkConsensusDeferredInputs{haplotype: "LA"}, ""},
		{"2", checkConsensusDeferredInputs{haplotype: "2"}, ""},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkConsensusDeferred(tc.in); got != tc.want {
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

// TestIsPhasedIUPAC matches upstream's "NpIu" haplotype encoding.
func TestIsPhasedIUPAC(t *testing.T) {
	cases := []struct {
		in   string
		want bool
	}{
		{"", false},
		{"R", false},
		{"2", false},
		{"2pIu", true},
		{"10pIu", true},
		{"2pIuX", false},
		{"pIu", false},
		{"apIu", false},
	}
	for _, tc := range cases {
		if got := isPhasedIUPAC(tc.in); got != tc.want {
			t.Errorf("isPhasedIUPAC(%q) = %v, want %v", tc.in, got, tc.want)
		}
	}
}
