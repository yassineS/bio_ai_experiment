package main

import "testing"

// TestCheckMendelian2Deferred locks in the upstream-flag-name surface
// that runMendelian2 hard-rejects rather than silently accepting.
// Per the project parity rule (docs/PARITY_ROADMAP.md "Definition of
// 1:1") every documented upstream flag must be recognised — either
// implemented or gracefully rejected with a roadmap pointer. A
// regression that drops any of these from the rejection set without
// implementing the underlying behaviour is a parity bug.
func TestCheckMendelian2Deferred(t *testing.T) {
	if got := checkMendelian2Deferred(checkMendelian2DeferredInputs{}); got != "" {
		t.Fatalf("empty inputs: got deferred=%q, want \"\"", got)
	}
	cases := []struct {
		name string
		in   checkMendelian2DeferredInputs
		want string
	}{
		{"rules", checkMendelian2DeferredInputs{rules: "GRCh38"}, "--rules"},
		{"rules-file", checkMendelian2DeferredInputs{rulesFile: "rules.txt"}, "--rules-file"},
		{"write-index", checkMendelian2DeferredInputs{writeIndex: "csi"}, "-W/--write-index"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkMendelian2Deferred(tc.in); got != tc.want {
				t.Errorf("deferred(%s): got %q, want %q", tc.name, got, tc.want)
			}
		})
	}
}
