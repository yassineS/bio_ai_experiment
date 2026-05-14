package tabix

import (
	"reflect"
	"testing"
)

func TestReg2bin(t *testing.T) {
	// Worked examples derived from the htslib paper (Li 2011).
	// At each level the smallest bin that fully contains [beg, end) is
	// returned. We verify a handful of cases that exercise every level.
	cases := []struct {
		beg, end int
		want     int
	}{
		// Level 5 (16 kbp tiles): beg and end in the same 16 kbp tile.
		{0, 1, 4681},         // first base
		{16383, 16384, 4681}, // last base of tile 0
		{16384, 16385, 4682}, // first base of tile 1
		// Level 4 (128 kbp): crosses a 16-kbp boundary but stays in 128 kbp.
		{0, 1 << 17, 585},
		// Level 3 (1 Mbp): straddles 128 kbp boundary.
		{0, 1 << 20, 73},
		// Level 2 (8 Mbp).
		{0, 1 << 23, 9},
		// Level 1 (64 Mbp).
		{0, 1 << 26, 1},
		// Level 0 (top): span larger than 64 Mbp.
		{0, 1 << 27, 0},
	}
	for _, tc := range cases {
		got := Reg2bin(tc.beg, tc.end)
		if got != tc.want {
			t.Errorf("Reg2bin(%d, %d) = %d, want %d", tc.beg, tc.end, got, tc.want)
		}
	}
}

func TestReg2binSamePosAtBoundaries(t *testing.T) {
	// Two adjacent 16-kbp tiles must yield distinct level-5 bins.
	a := Reg2bin(0, 16384)
	b := Reg2bin(16384, 32768)
	if a == b {
		t.Errorf("adjacent 16-kbp tiles share a bin: %d", a)
	}
	if a != 4681 || b != 4682 {
		t.Errorf("unexpected bin numbers: a=%d b=%d", a, b)
	}
}

func TestReg2bins(t *testing.T) {
	// A 1-bp region falls in only its level-5 tile plus the ancestor bins.
	got := Reg2bins(0, 1)
	want := []int{0, 1, 9, 73, 585, 4681}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("Reg2bins(0, 1) = %v, want %v", got, want)
	}
}

func TestReg2binsSpansTiles(t *testing.T) {
	// A range that straddles two 16-kbp tiles should yield both level-5
	// bins.
	got := Reg2bins(16380, 16400)
	want := []int{0, 1, 9, 73, 585, 4681, 4682}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("Reg2bins crossing tile = %v, want %v", got, want)
	}
}

func TestReg2binsEmpty(t *testing.T) {
	// An empty / inverted region should still yield bin 0.
	got := Reg2bins(100, 100)
	if len(got) != 1 || got[0] != 0 {
		t.Errorf("Reg2bins for empty region = %v, want [0]", got)
	}
}

func TestLinearTile(t *testing.T) {
	cases := []struct {
		pos  int
		want int
	}{
		{0, 0},
		{1, 0},
		{16383, 0},
		{16384, 1},
		{32768, 2},
		{-1, 0},
	}
	for _, tc := range cases {
		if got := LinearTile(tc.pos); got != tc.want {
			t.Errorf("LinearTile(%d) = %d, want %d", tc.pos, got, tc.want)
		}
	}
}
