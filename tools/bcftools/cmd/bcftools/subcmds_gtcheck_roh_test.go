package main

import "testing"

// TestCheckGtcheckDeferred locks in the upstream-flag-name surface
// that runGtcheck hard-rejects rather than silently accepting. Per
// the project parity rule (docs/PARITY_ROADMAP.md "Definition of 1:1")
// every documented upstream flag must be recognised — either
// implemented or gracefully rejected with a roadmap pointer. A
// regression that drops any of these from the rejection set without
// implementing the underlying behaviour is a parity bug.
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
		// --n-matches and -O z are now real-ported (the deferred
		// checker no longer rejects them); only --cluster and
		// --distinctive-sites remain in the rejection set.
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkGtcheckDeferred(tc.in); got != tc.want {
				t.Errorf("deferred(%s): got %q, want %q", tc.name, got, tc.want)
			}
		})
	}
}

// TestCheckRohDeferred is the corresponding test for roh. With the
// genetic-map / rec-rate / estimate-AF / buffer-size / viterbi-training
// features now implemented, only gzip output (-O z) remains deferred.
func TestCheckRohDeferred(t *testing.T) {
	if got := checkRohDeferred(checkRohDeferredInputs{outputType: "sr"}); got != "" {
		t.Fatalf("default outputType=sr: got deferred=%q, want \"\"", got)
	}
	if got := checkRohDeferred(checkRohDeferredInputs{outputType: "srz"}); got != "-O z (compressed output)" {
		t.Errorf("output-z: got %q, want %q", got, "-O z (compressed output)")
	}
}

// TestSplitSamplesArg validates the upstream "qry:" / "gt:" prefix
// handling on `-s`. Un-prefixed lists are rejected to match upstream
// vcfgtcheck.c's "Which one? Query samples (qry:...) or genotype
// samples (gt:...)?" diagnostic.
func TestSplitSamplesArg(t *testing.T) {
	cases := []struct {
		in      string
		wantQry []string
		wantGT  []string
		wantErr bool
	}{
		{"", nil, nil, false},
		{"qry:a,b", []string{"a", "b"}, nil, false},
		{"gt:x,y,z", nil, []string{"x", "y", "z"}, false},
		{"a,b,c", nil, nil, true},
	}
	for _, c := range cases {
		gotQ, gotG, gotErr := splitSamplesArg(c.in)
		if (gotErr != nil) != c.wantErr {
			t.Errorf("splitSamplesArg(%q) err = %v, wantErr = %v", c.in, gotErr, c.wantErr)
		}
		if c.wantErr {
			continue
		}
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
