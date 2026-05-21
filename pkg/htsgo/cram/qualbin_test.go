package cram

import (
	"bytes"
	"testing"
)

// TestBinTableBoundaries asserts every binning table maps each input
// quality 0-255 to exactly the documented representative for its bin.
func TestBinTableBoundaries(t *testing.T) {
	cases := []struct {
		name   string
		scheme QualityBinning
		bins   []binBoundary
	}{
		{"illumina-8", BinningIllumina8, illumina8Bins},
		{"illumina-4", BinningIllumina4, illumina4Bins},
		{"illumina-2", BinningIllumina2, illumina2Bins},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			tbl := tc.scheme.table()
			// Every entry 0-255 must fall in exactly one bin and carry that
			// bin's representative.
			for q := 0; q <= 255; q++ {
				want := byte(0)
				found := false
				for _, b := range tc.bins {
					if byte(q) >= b.lo && byte(q) <= b.hi {
						if found {
							t.Fatalf("quality %d covered by more than one bin", q)
						}
						want, found = b.rep, true
					}
				}
				if !found {
					t.Fatalf("quality %d not covered by any bin", q)
				}
				if tbl[q] != want {
					t.Errorf("table[%d] = %d, want %d", q, tbl[q], want)
				}
			}
		})
	}
}

// TestIllumina8CanonicalTable pins the exact Illumina 8-level boundaries
// and representatives from the "Reducing Whole-Genome Data Storage
// Footprint" technical note, so an accidental table edit is caught.
func TestIllumina8CanonicalTable(t *testing.T) {
	tbl := BinningIllumina8.table()
	type span struct {
		lo, hi int
		rep    byte
	}
	want := []span{
		{0, 2, 0},
		{3, 9, 6},
		{10, 19, 15},
		{20, 24, 22},
		{25, 29, 27},
		{30, 34, 33},
		{35, 39, 37},
		{40, 255, 40},
	}
	for _, s := range want {
		for q := s.lo; q <= s.hi; q++ {
			if tbl[q] != s.rep {
				t.Errorf("illumina-8 table[%d] = %d, want %d", q, tbl[q], s.rep)
			}
		}
	}
}

// TestBinQualityIdempotent confirms binning a value already binned by the
// same scheme is a no-op: bin(bin(q)) == bin(q) for every scheme.
func TestBinQualityIdempotent(t *testing.T) {
	schemes := []QualityBinning{
		BinningNone, BinningIllumina8, BinningIllumina4, BinningIllumina2,
	}
	input := make([]byte, 64)
	for i := range input {
		input[i] = byte(i)
	}
	for _, s := range schemes {
		once := s.BinQuality(input)
		twice := s.BinQuality(once)
		if !bytes.Equal(once, twice) {
			t.Errorf("scheme %s: binning is not idempotent\n once: %v\ntwice: %v", s, once, twice)
		}
	}
}

// TestBinQualityNoMutation confirms BinQuality never modifies its input
// slice — the caller's record quality must stay untouched.
func TestBinQualityNoMutation(t *testing.T) {
	input := []byte{0, 5, 12, 30, 45}
	orig := append([]byte(nil), input...)
	got := BinningIllumina8.BinQuality(input)
	if !bytes.Equal(input, orig) {
		t.Errorf("BinQuality mutated its input: got %v, want %v", input, orig)
	}
	want := []byte{0, 6, 15, 33, 40}
	if !bytes.Equal(got, want) {
		t.Errorf("BinQuality = %v, want %v", got, want)
	}
}

// TestBinQualityNoneIsCopy confirms BinningNone returns an exact,
// independent copy of the input.
func TestBinQualityNoneIsCopy(t *testing.T) {
	input := []byte{1, 2, 3, 40, 41}
	got := BinningNone.BinQuality(input)
	if !bytes.Equal(got, input) {
		t.Errorf("BinningNone.BinQuality = %v, want exact copy %v", got, input)
	}
	if len(input) > 0 && &got[0] == &input[0] {
		t.Error("BinningNone.BinQuality returned the input slice, want an independent copy")
	}
}

// TestBinQualityPreservesNoQualSentinel confirms the SAM no-quality
// sentinel 0xff survives every binning scheme unchanged.
func TestBinQualityPreservesNoQualSentinel(t *testing.T) {
	for _, s := range []QualityBinning{BinningIllumina8, BinningIllumina4, BinningIllumina2} {
		got := s.BinQuality([]byte{0xff, 0xff, 0xff})
		for i, q := range got {
			if q != 0xff {
				t.Errorf("scheme %s: byte %d = %d, want 0xff preserved", s, i, q)
			}
		}
	}
}

// TestBinQualityEmpty confirms an empty input yields an empty result.
func TestBinQualityEmpty(t *testing.T) {
	if got := BinningIllumina8.BinQuality(nil); len(got) != 0 {
		t.Errorf("BinQuality(nil) = %v, want empty", got)
	}
}

// TestQualityBinningString checks the human-readable scheme names.
func TestQualityBinningString(t *testing.T) {
	cases := map[QualityBinning]string{
		BinningNone:      "none",
		BinningIllumina8: "illumina-8",
		BinningIllumina4: "illumina-4",
		BinningIllumina2: "illumina-2",
	}
	for s, want := range cases {
		if got := s.String(); got != want {
			t.Errorf("%d.String() = %q, want %q", int(s), got, want)
		}
	}
}

// TestQualityBinningValid confirms valid recognises the four schemes and
// rejects an out-of-range value.
func TestQualityBinningValid(t *testing.T) {
	for _, s := range []QualityBinning{BinningNone, BinningIllumina8, BinningIllumina4, BinningIllumina2} {
		if !s.valid() {
			t.Errorf("scheme %d should be valid", int(s))
		}
	}
	if QualityBinning(99).valid() {
		t.Error("scheme 99 should be invalid")
	}
	if QualityBinning(-1).valid() {
		t.Error("scheme -1 should be invalid")
	}
}
