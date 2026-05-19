package sam

import (
	"testing"
)

func TestCigarParseAndString(t *testing.T) {
	tests := []struct {
		in       string
		query    int
		ref      int
		outAgain string
	}{
		{"10M", 10, 10, "10M"},
		{"5S10M2I3D6M", 23, 19, "5S10M2I3D6M"},
		{"100=", 100, 100, "100="},
		{"3H10M3H", 10, 10, "3H10M3H"},
		{"*", 0, 0, "*"},
	}
	for _, tc := range tests {
		t.Run(tc.in, func(t *testing.T) {
			cig, err := ParseCigar(tc.in)
			if err != nil {
				t.Fatalf("ParseCigar: %v", err)
			}
			if got := cig.QueryLength(); got != tc.query {
				t.Errorf("QueryLength: got %d, want %d", got, tc.query)
			}
			if got := cig.ReferenceLength(); got != tc.ref {
				t.Errorf("ReferenceLength: got %d, want %d", got, tc.ref)
			}
			if got := cig.String(); got != tc.outAgain {
				t.Errorf("String: got %q, want %q", got, tc.outAgain)
			}
		})
	}
}

func TestCigarParseErrors(t *testing.T) {
	cases := []string{"", "10Z", "abc", "10", "M"}
	for _, c := range cases {
		if _, err := ParseCigar(c); err == nil {
			t.Errorf("expected error for ParseCigar(%q)", c)
		}
	}
}

func TestAuxFormatSAM(t *testing.T) {
	tests := []struct {
		a    Aux
		want string
	}{
		{Aux{Tag: "NM", Type: 'i', Value: int64(3)}, "NM:i:3"},
		{Aux{Tag: "AS", Type: 'i', Value: int64(-12)}, "AS:i:-12"},
		{Aux{Tag: "XF", Type: 'f', Value: 1.5}, "XF:f:1.5"},
		{Aux{Tag: "RG", Type: 'Z', Value: "rg1"}, "RG:Z:rg1"},
		{Aux{Tag: "XA", Type: 'A', Value: "U"}, "XA:A:U"},
		{Aux{Tag: "BB", Type: 'B', ArrayType: 'i', ArrayValues: []interface{}{int64(1), int64(2), int64(3)}}, "BB:B:i,1,2,3"},
		{Aux{Tag: "BF", Type: 'B', ArrayType: 'f', ArrayValues: []interface{}{1.5, 2.5}}, "BF:B:f,1.5,2.5"},
		{Aux{Tag: "HX", Type: 'H', Value: "DEADBEEF"}, "HX:H:DEADBEEF"},
	}
	for _, tc := range tests {
		if got := tc.a.FormatSAM(); got != tc.want {
			t.Errorf("FormatSAM(%+v): got %q, want %q", tc.a, got, tc.want)
		}
	}
}

func TestParseAux(t *testing.T) {
	tests := []struct {
		in   string
		want Aux
	}{
		{"NM:i:3", Aux{Tag: "NM", Type: 'i', Value: int64(3)}},
		{"XF:f:1.5", Aux{Tag: "XF", Type: 'f', Value: 1.5}},
		{"RG:Z:rg1", Aux{Tag: "RG", Type: 'Z', Value: "rg1"}},
		{"XA:A:U", Aux{Tag: "XA", Type: 'A', Value: "U"}},
	}
	for _, tc := range tests {
		got, err := ParseAux(tc.in)
		if err != nil {
			t.Fatalf("ParseAux(%q): %v", tc.in, err)
		}
		if got.Tag != tc.want.Tag || got.Type != tc.want.Type || got.Value != tc.want.Value {
			t.Errorf("ParseAux(%q): got %+v, want %+v", tc.in, got, tc.want)
		}
	}
}

func TestParseAuxArray(t *testing.T) {
	got, err := ParseAux("ZB:B:i,1,2,3")
	if err != nil {
		t.Fatalf("ParseAux: %v", err)
	}
	if got.Type != 'B' || got.ArrayType != 'i' || len(got.ArrayValues) != 3 {
		t.Fatalf("ParseAux B:i unexpected: %+v", got)
	}
	if v, _ := got.ArrayValues[2].(int64); v != 3 {
		t.Errorf("array value [2]: got %v, want 3", got.ArrayValues[2])
	}

	gotf, err := ParseAux("ZF:B:f,0.5,1.5,2.5")
	if err != nil {
		t.Fatalf("ParseAux B:f: %v", err)
	}
	if gotf.ArrayType != 'f' || len(gotf.ArrayValues) != 3 {
		t.Fatalf("ParseAux B:f unexpected: %+v", gotf)
	}
}

func TestParseAuxErrors(t *testing.T) {
	bad := []string{"", "X", "NM:i", "NM:i:abc", "XF:f:notnum", "XA:A:two", "ZZ:Q:1"}
	for _, c := range bad {
		if _, err := ParseAux(c); err == nil {
			t.Errorf("expected error for ParseAux(%q)", c)
		}
	}
}

func TestRecordFlagHelpers(t *testing.T) {
	r := &Record{Flag: FlagPaired | FlagProperPair | FlagRead1 | FlagReverse}
	if !r.IsPaired() || !r.IsProperPair() || !r.IsRead1() {
		t.Errorf("missing expected flag predicates")
	}
	if r.IsRead2() || r.IsSecondary() || r.IsSupplementary() {
		t.Errorf("unexpected flag predicate true")
	}
	if !r.IsPrimary() || !r.IsMapped() {
		t.Errorf("primary/mapped predicates wrong")
	}

	rUnmapped := &Record{Flag: FlagUnmapped | FlagPaired | FlagMateUnmapped | FlagQCFail | FlagDuplicate}
	if !rUnmapped.IsUnmapped() || !rUnmapped.IsMateUnmapped() || !rUnmapped.IsQCFail() || !rUnmapped.IsDuplicate() {
		t.Errorf("unmapped flags wrong")
	}
	rSec := &Record{Flag: FlagSecondary}
	rSup := &Record{Flag: FlagSupplementary}
	if rSec.IsPrimary() || rSup.IsPrimary() {
		t.Errorf("secondary/supplementary should not be primary")
	}
	if !rSec.IsSecondary() || !rSup.IsSupplementary() {
		t.Errorf("secondary/supplementary predicates wrong")
	}
}

func TestRecordEndPosition(t *testing.T) {
	cig, _ := ParseCigar("5M2D3M")
	r := &Record{Pos: 100, Cigar: cig}
	if got := r.EndPosition(); got != 109 {
		t.Errorf("EndPosition: got %d, want 109", got)
	}
	rEmpty := &Record{Pos: 50}
	if got := rEmpty.EndPosition(); got != 50 {
		t.Errorf("empty cigar EndPosition: got %d, want 50", got)
	}
}

func TestRecordGetAux(t *testing.T) {
	r := &Record{Aux: []Aux{
		{Tag: "NM", Type: 'i', Value: int64(2)},
		{Tag: "RG", Type: 'Z', Value: "rg1"},
	}}
	if a, ok := r.GetAux("NM"); !ok || a.Type != 'i' {
		t.Errorf("GetAux NM: %+v %v", a, ok)
	}
	if v, ok := r.Aux[0].Int(); !ok || v != 2 {
		t.Errorf("Int() helper: %d %v", v, ok)
	}
	if s, ok := r.Aux[1].String(); !ok || s != "rg1" {
		t.Errorf("String() helper: %q %v", s, ok)
	}
	if _, ok := r.GetAux("XX"); ok {
		t.Error("GetAux of missing tag should be false")
	}
	// Second lookup uses the cached index path.
	if _, ok := r.GetAux("NM"); !ok {
		t.Error("second GetAux NM failed")
	}
}

func TestCigarOpAccessors(t *testing.T) {
	op := CigarOp(uint32(7)<<4 | CigarSoftClip)
	if op.Length() != 7 || op.Op() != CigarSoftClip || op.Char() != 'S' {
		t.Errorf("CigarOp accessors wrong: len=%d op=%d ch=%c", op.Length(), op.Op(), op.Char())
	}
	bad := CigarOp(uint32(1)<<4 | 0xf) // unused op
	if bad.Char() != '?' {
		t.Errorf("expected '?' for invalid op, got %c", bad.Char())
	}
}
