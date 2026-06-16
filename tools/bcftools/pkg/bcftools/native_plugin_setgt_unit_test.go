// Binary-free unit tests for setGT's pure helpers (native_plugin_setgt.go):
// parseBinomExpr (the "b:TAG CMP VAL" selector) and scanFloatPrefix. These
// exercise the parser directly on a zero-value plugin, needing no VCF header or
// upstream binary.
package bcftools

import "testing"

// TestUnitSetGTBinomExpr pins parseBinomExpr: it must extract the FORMAT tag,
// comparison operator, and threshold from a "b:TAG CMP VAL" string, set the
// sgtBinom bit, and reject malformed input.
func TestUnitSetGTBinomExpr(t *testing.T) {
	tests := []struct {
		name    string
		expr    string
		wantTag string
		wantVal float64
		// cmpProbe is fed to the parsed comparator as cmp(a,b); cmpWant is the
		// expected result, used to identify which operator was selected.
		cmpA, cmpB float64
		cmpWant    bool
		wantErr    bool
	}{
		{
			name: "greater-than scientific", expr: "b:AD>1e-3",
			wantTag: "AD", wantVal: 1e-3, cmpA: 2, cmpB: 1, cmpWant: true,
		},
		{
			name: "less-or-equal", expr: "b:AD<=0.05",
			wantTag: "AD", wantVal: 0.05, cmpA: 0.05, cmpB: 0.05, cmpWant: true,
		},
		{
			name: "greater-or-equal with spaces", expr: "b:DP >= 10",
			wantTag: "DP", wantVal: 10, cmpA: 10, cmpB: 10, cmpWant: true,
		},
		{
			name: "equality double-equals", expr: "b:AD==0.5",
			wantTag: "AD", wantVal: 0.5, cmpA: 0.5, cmpB: 0.5, cmpWant: true,
		},
		{
			name: "equality single-equals", expr: "b:AD=0.5",
			wantTag: "AD", wantVal: 0.5, cmpA: 0.4, cmpB: 0.5, cmpWant: false,
		},
		{
			name: "less-than negative value", expr: "b:AD<-0.1",
			wantTag: "AD", wantVal: -0.1, cmpA: -0.2, cmpB: -0.1, cmpWant: true,
		},
		{name: "missing colon", expr: "bAD>1", wantErr: true},
		{name: "missing operator", expr: "b:AD1", wantErr: true},
		{name: "missing tag", expr: "b:>1", wantErr: true},
		{name: "missing value", expr: "b:AD>", wantErr: true},
		{name: "non-numeric value", expr: "b:AD>abc", wantErr: true},
		{name: "trailing garbage after value", expr: "b:AD>1xyz", wantErr: true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			p := &setGTPlugin{}
			err := p.parseBinomExpr(tt.expr)
			if tt.wantErr {
				if err == nil {
					t.Fatalf("parseBinomExpr(%q) = nil error, want error", tt.expr)
				}
				return
			}
			if err != nil {
				t.Fatalf("parseBinomExpr(%q) unexpected error: %v", tt.expr, err)
			}
			if p.binomTag != tt.wantTag {
				t.Errorf("binomTag = %q, want %q", p.binomTag, tt.wantTag)
			}
			if p.binomVal != tt.wantVal {
				t.Errorf("binomVal = %v, want %v", p.binomVal, tt.wantVal)
			}
			if p.tgtMask&sgtBinom == 0 {
				t.Errorf("sgtBinom bit not set after parsing %q", tt.expr)
			}
			if p.binomCmp == nil {
				t.Fatal("binomCmp is nil")
			}
			if got := p.binomCmp(tt.cmpA, tt.cmpB); got != tt.cmpWant {
				t.Errorf("binomCmp(%v,%v) = %v, want %v", tt.cmpA, tt.cmpB, got, tt.cmpWant)
			}
		})
	}
}

// TestUnitSetGTScanFloatPrefix pins scanFloatPrefix: the leading run that strtod
// would consume, and the trailing-garbage / no-digit cases.
func TestUnitSetGTScanFloatPrefix(t *testing.T) {
	tests := []struct {
		in       string
		wantStr  string
		wantSpan int
	}{
		{"1.5", "1.5", 3},
		{"1e-3rest", "1e-3", 4},
		{"-0.25 ", "-0.25", 5},
		{"+10", "+10", 3},
		{"42abc", "42", 2},
		{".5", ".5", 2},
		{"abc", "", 0},
		{"", "", 0},
		{"1e", "1", 1}, // a bare exponent marker is not consumed
	}
	for _, tt := range tests {
		gotStr, gotSpan := scanFloatPrefix(tt.in)
		if gotStr != tt.wantStr || gotSpan != tt.wantSpan {
			t.Errorf("scanFloatPrefix(%q) = (%q,%d), want (%q,%d)", tt.in, gotStr, gotSpan, tt.wantStr, tt.wantSpan)
		}
	}
}
