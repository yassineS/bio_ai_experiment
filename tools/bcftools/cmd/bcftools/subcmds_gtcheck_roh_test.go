package main

import "testing"

// TestCheckGtcheckDeferred locks in the upstream-flag-name surface
// that runGtcheck hard-rejects rather than silently accepting. Per
// the "every documented upstream flag must be recognised — either
// implemented or gracefully rejected with a roadmap pointer" rule
// (docs/PARITY_ROADMAP.md#definition-of-11), a regression that drops
// any of these from the rejection set without implementing the
// underlying behaviour is a parity bug.
func TestCheckGtcheckDeferred(t *testing.T) {
	if got := checkGtcheckDeferred(checkGtcheckDeferredInputs{outputType: "t"}); got != "" {
		t.Fatalf("default outputType=t: got deferred=%q, want \"\"", got)
	}
	cases := []struct {
		name string
		in   checkGtcheckDeferredInputs
		want string
	}{
		{"cluster", checkGtcheckDeferredInputs{cluster: "2,4"}, "--cluster"},
		{"distinctive-sites", checkGtcheckDeferredInputs{distinctiveSites: "0.1"}, "--distinctive-sites"},
		{"n-matches", checkGtcheckDeferredInputs{nMatches: 5}, "--n-matches"},
		{"output-type-z", checkGtcheckDeferredInputs{outputType: "z"}, "-O z (compressed output)"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkGtcheckDeferred(tc.in); got != tc.want {
				t.Errorf("deferred(%s): got %q, want %q", tc.name, got, tc.want)
			}
		})
	}
}

// TestCheckRohDeferred is the corresponding test for roh.
func TestCheckRohDeferred(t *testing.T) {
	if got := checkRohDeferred(checkRohDeferredInputs{outputType: "sr"}); got != "" {
		t.Fatalf("default outputType=sr: got deferred=%q, want \"\"", got)
	}
	cases := []struct {
		name string
		in   checkRohDeferredInputs
		want string
	}{
		{"buffer-size", checkRohDeferredInputs{bufferSize: "100000"}, "-b/--buffer-size"},
		{"estimate-AF", checkRohDeferredInputs{estimateAF: "GT,-"}, "-e/--estimate-AF"},
		{"genetic-map", checkRohDeferredInputs{geneticMap: "map.txt"}, "-m/--genetic-map"},
		{"rec-rate", checkRohDeferredInputs{recRate: 1e-9}, "-M/--rec-rate"},
		{"viterbi-training", checkRohDeferredInputs{viterbiTraining: 1e-10}, "-V/--viterbi-training"},
		{"output-z", checkRohDeferredInputs{outputType: "srz"}, "-O z (compressed output)"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkRohDeferred(tc.in); got != tc.want {
				t.Errorf("deferred(%s): got %q, want %q", tc.name, got, tc.want)
			}
		})
	}
}

// TestSplitSamplesArg validates the upstream "qry:" / "gt:" prefix
// handling on `-s`.
func TestSplitSamplesArg(t *testing.T) {
	cases := []struct {
		in      string
		wantQry []string
		wantGT  []string
	}{
		{"", nil, nil},
		{"a,b,c", []string{"a", "b", "c"}, []string{"a", "b", "c"}},
		{"qry:a,b", []string{"a", "b"}, nil},
		{"gt:x,y,z", nil, []string{"x", "y", "z"}},
	}
	for _, c := range cases {
		gotQ, gotG := splitSamplesArg(c.in)
		if !stringsEqual(gotQ, c.wantQry) {
			t.Errorf("splitSamplesArg(%q) qry = %v, want %v", c.in, gotQ, c.wantQry)
		}
		if !stringsEqual(gotG, c.wantGT) {
			t.Errorf("splitSamplesArg(%q) gt = %v, want %v", c.in, gotG, c.wantGT)
		}
	}
}

func stringsEqual(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
