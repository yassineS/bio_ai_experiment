package main

import "testing"

// TestCheckGtcheckDeferred locks the upstream-flag-name surface that
// runGtcheck hard-rejects rather than silently accepting. Same parity
// rule as TestCheckConvertDeferred: dropping any entry from this list
// without implementing the underlying behaviour is a regression.
func TestCheckGtcheckDeferred(t *testing.T) {
	if got := checkGtcheckDeferred(checkGtcheckDeferredInputs{}); got != "" {
		t.Fatalf("empty inputs: got deferred=%q, want \"\"", got)
	}
	cases := []struct {
		name string
		in   checkGtcheckDeferredInputs
		want string
	}{
		{"all-sites", checkGtcheckDeferredInputs{allSites: true}, "--all-sites"},
		{"homs-only", checkGtcheckDeferredInputs{homsOnly: true}, "--homs-only"},
		{"no-HWE-prob", checkGtcheckDeferredInputs{noHWEProb: true}, "--no-HWE-prob"},
		{"GTs-only", checkGtcheckDeferredInputs{gtsOnly: "1"}, "--GTs-only"},
		{"pl-units", checkGtcheckDeferredInputs{plUnits: "PL"}, "--pl-units"},
		{"dosage", checkGtcheckDeferredInputs{dosage: true}, "--dosage"},
		{"tags", checkGtcheckDeferredInputs{tagsFlag: "GT"}, "--tags"},
		{"cluster", checkGtcheckDeferredInputs{cluster: "1"}, "--cluster"},
		{"normalize", checkGtcheckDeferredInputs{normalize: "1"}, "--normalize"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkGtcheckDeferred(tc.in); got != tc.want {
				t.Errorf("deferred(%s): got %q, want %q", tc.name, got, tc.want)
			}
		})
	}
}

// TestCheckRohDeferred is the analogous coverage for runRoh.
func TestCheckRohDeferred(t *testing.T) {
	if got := checkRohDeferred(checkRohDeferredInputs{}); got != "" {
		t.Fatalf("empty inputs: got deferred=%q, want \"\"", got)
	}
	cases := []struct {
		name string
		in   checkRohDeferredInputs
		want string
	}{
		{"rec-rate", checkRohDeferredInputs{recRate: 1e-8}, "-M/--rec-rate"},
		{"genetic-map", checkRohDeferredInputs{geneticMap: "map.tsv"}, "-V/--genetic-map"},
		{"buffer-size", checkRohDeferredInputs{bufferSize: "1000"}, "--buffer-size"},
		{"skip-indels", checkRohDeferredInputs{skipIndels: true}, "--skip-indels"},
		{"include-noalt", checkRohDeferredInputs{includeNo: true}, "--include-noalt"},
		{"estimate-AF", checkRohDeferredInputs{estimateAF: true}, "--estimate-AF"},
		{"af-file", checkRohDeferredInputs{afFile: "af.tsv"}, "--AF-file"},
		{"hw-to-az", checkRohDeferredInputs{hwToAz: 1e-4}, "-a/--hw-to-az"},
		{"az-to-hw", checkRohDeferredInputs{azToHw: 1e-3}, "-H/--az-to-hw"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkRohDeferred(tc.in); got != tc.want {
				t.Errorf("deferred(%s): got %q, want %q", tc.name, got, tc.want)
			}
		})
	}
}

// TestPhredToProb is a quick sanity check on the helper used to convert
// `-G INT` into a linear genotype-error rate.
func TestPhredToProb(t *testing.T) {
	cases := []struct {
		in   float64
		want float64
	}{
		{0, 1},
		{10, 0.1},
		{20, 0.01},
		{30, 0.001},
	}
	for _, c := range cases {
		got := phredToProb(c.in)
		if got > c.want*1.001 || got < c.want*0.999 {
			t.Errorf("phredToProb(%f) = %g, want %g", c.in, got, c.want)
		}
	}
}
