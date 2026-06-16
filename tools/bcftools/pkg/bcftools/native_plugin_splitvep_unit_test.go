// Binary-free unit tests for the standalone split-vep pure helpers
// (native_plugin_splitvep*.go): the format tokenizer, the :TYPE / index-range /
// column-type-code parsers, the QUAL renderer, and the severity scale logic
// (which is built from the in-tree default text, requiring no VCF header or
// upstream binary).
package bcftools

import (
	"reflect"
	"testing"
)

// TestUnitSplitVepFields pins the "a|b(...)|c" tokenizer: '|' separates fields
// and '(' ends the current field name (the bracketed text up to the next '|' is
// discarded).
func TestUnitSplitVepFields(t *testing.T) {
	tests := []struct {
		in   string
		want []string
	}{
		{"Allele|Consequence|IMPACT", []string{"Allele", "Consequence", "IMPACT"}},
		{"Allele|cDNA_position(1-based)|Codons", []string{"Allele", "cDNA_position", "Codons"}},
		{"Single", []string{"Single"}},
		{"A||B", []string{"A", "", "B"}},
	}
	for _, tt := range tests {
		if got := splitVepFields(tt.in); !reflect.DeepEqual(got, tt.want) {
			t.Errorf("splitVepFields(%q) = %v, want %v", tt.in, got, tt.want)
		}
	}
}

// TestUnitSplitVepParseSvType pins the :TYPE suffix spellings (case-insensitive)
// and the unknown-type rejection.
func TestUnitSplitVepParseSvType(t *testing.T) {
	tests := []struct {
		in     string
		want   int
		wantOK bool
	}{
		{"string", svTypeStr, true},
		{"str", svTypeStr, true},
		{"STRING", svTypeStr, true},
		{"integer", svTypeInt, true},
		{"int", svTypeInt, true},
		{"Int", svTypeInt, true},
		{"float", svTypeReal, true},
		{"real", svTypeReal, true},
		{"REAL", svTypeReal, true},
		{"bogus", 0, false},
		{"", 0, false},
	}
	for _, tt := range tests {
		got, ok := parseSvType(tt.in)
		if ok != tt.wantOK || (ok && got != tt.want) {
			t.Errorf("parseSvType(%q) = (%d,%v), want (%d,%v)", tt.in, got, ok, tt.want, tt.wantOK)
		}
	}
}

// TestUnitSplitVepParseIndexRange pins "N" and "N-M" parsing and the invalid
// forms.
func TestUnitSplitVepParseIndexRange(t *testing.T) {
	tests := []struct {
		in     string
		lo, hi int
		ok     bool
	}{
		{"5", 5, 5, true},
		{"2-7", 2, 7, true},
		{"0", 0, 0, true},
		{"10-3", 10, 3, true}, // the parser does not enforce lo<=hi
		// A leading '-' is NOT a range separator (the parser requires dash>0), so
		// "-5" falls through to a single Atoi and parses as the integer -5.
		{"-5", -5, -5, true},
		{"abc", 0, 0, false},
		{"2-x", 0, 0, false},
		{"", 0, 0, false},
	}
	for _, tt := range tests {
		lo, hi, ok := parseIndexRange(tt.in)
		if ok != tt.ok || (ok && (lo != tt.lo || hi != tt.hi)) {
			t.Errorf("parseIndexRange(%q) = (%d,%d,%v), want (%d,%d,%v)", tt.in, lo, hi, ok, tt.lo, tt.hi, tt.ok)
		}
	}
}

// TestUnitSplitVepColumnTypeCode pins the column-type-string mapping, including
// Flag mapping to the String renderer.
func TestUnitSplitVepColumnTypeCode(t *testing.T) {
	tests := []struct {
		in     string
		want   int
		wantOK bool
	}{
		{"Float", svTypeReal, true},
		{"Integer", svTypeInt, true},
		{"String", svTypeStr, true},
		{"Flag", svTypeStr, true},
		{"float", 0, false}, // case-sensitive, unlike parseSvType
		{"Bogus", 0, false},
	}
	for _, tt := range tests {
		got, ok := svColumnTypeCode(tt.in)
		if ok != tt.wantOK || (ok && got != tt.want) {
			t.Errorf("svColumnTypeCode(%q) = (%d,%v), want (%d,%v)", tt.in, got, ok, tt.want, tt.wantOK)
		}
	}
}

// TestUnitSplitVepFormatQual pins the QUAL renderer: a negative (missing) QUAL
// renders ".", otherwise the %g-style float.
func TestUnitSplitVepFormatQual(t *testing.T) {
	tests := []struct {
		in   float64
		want string
	}{
		{-1, "."},
		{-0.5, "."},
		{0, "0"},
		{42, "42"},
		{3.5, "3.5"},
	}
	for _, tt := range tests {
		if got := formatQual(tt.in); got != tt.want {
			t.Errorf("formatQual(%v) = %q, want %q", tt.in, got, tt.want)
		}
	}
}

// TestUnitSplitVepSeverityOrdering pins the relative severity ordering produced
// by csqToSeverity from the in-tree default scale. The absolute tier numbers
// depend on the default text, so the test asserts the documented ordering
// (transcript_ablation is the most severe, intergenic the least) rather than
// hard-coding tier indices.
func TestUnitSplitVepSeverityOrdering(t *testing.T) {
	p := &splitVepPlugin{}
	if err := p.initSeverityScale(); err != nil {
		t.Fatalf("initSeverityScale: %v", err)
	}
	sev := func(term string) int {
		lo, hi := p.csqToSeverity(term, -1)
		if lo != hi {
			t.Fatalf("single term %q spans tiers [%d,%d]", term, lo, hi)
		}
		return hi
	}
	// Ascending severity per the default scale comment.
	ascending := []string{
		"intergenic_variant",
		"intron_variant",
		"missense_variant",
		"stop_gained",
		"transcript_ablation",
	}
	prev := -1 << 30
	for _, term := range ascending {
		s := sev(term)
		if s <= prev {
			t.Fatalf("severity(%q) = %d, expected strictly greater than previous %d", term, s, prev)
		}
		prev = s
	}
}

// TestUnitSplitVepSeverityMulti pins that an '&'-joined consequence spans the
// min and max tiers of its terms.
func TestUnitSplitVepSeverityMulti(t *testing.T) {
	p := &splitVepPlugin{}
	if err := p.initSeverityScale(); err != nil {
		t.Fatalf("initSeverityScale: %v", err)
	}
	loLow, _ := p.csqToSeverity("intron_variant", -1)
	_, hiHigh := p.csqToSeverity("stop_gained", -1)
	minSev, maxSev := p.csqToSeverity("intron_variant&stop_gained", -1)
	if minSev != loLow {
		t.Errorf("combined min tier = %d, want %d (intron)", minSev, loLow)
	}
	if maxSev != hiHigh {
		t.Errorf("combined max tier = %d, want %d (stop_gained)", maxSev, hiHigh)
	}
}

// TestUnitSplitVepSeverityPass pins csqSeverityPass for the "any" sentinel, an
// exact-tier match, and an open-ended minimum range.
func TestUnitSplitVepSeverityPass(t *testing.T) {
	p := &splitVepPlugin{}
	if err := p.initSeverityScale(); err != nil {
		t.Fatalf("initSeverityScale: %v", err)
	}
	intronSev, _ := p.csqToSeverity("intron_variant", -1)
	stopSev, _ := p.csqToSeverity("stop_gained", -1)

	// any:any (both svSelectAny) passes everything.
	p.minSev, p.maxSev = svSelectAny, svSelectAny
	if !p.csqSeverityPass("intron_variant") {
		t.Error("any severity should pass intron_variant")
	}

	// Exact tier (min==max) passes only that tier.
	p.minSev, p.maxSev = stopSev, stopSev
	if !p.csqSeverityPass("stop_gained") {
		t.Error("exact stop_gained tier should pass stop_gained")
	}
	if p.csqSeverityPass("intron_variant") {
		t.Error("exact stop_gained tier should not pass intron_variant")
	}

	// Open-ended minimum (max == a large tier): stop_gained passes, intron fails.
	p.minSev, p.maxSev = stopSev, 1<<30
	if !p.csqSeverityPass("stop_gained") {
		t.Error("min=stop should pass stop_gained")
	}
	if p.csqSeverityPass("intron_variant") {
		t.Errorf("min=stop (intronSev=%d) should not pass intron_variant", intronSev)
	}
}
