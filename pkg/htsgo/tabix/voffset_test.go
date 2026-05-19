package tabix

import "testing"

func TestVOffsetRoundTrip(t *testing.T) {
	cases := []struct {
		coff int64
		uoff int
	}{
		{0, 0},
		{1, 0},
		{0, 1},
		{0xDEADBEEF, 0x1234},
		{(1 << 47) - 1, 0xFFFE},
	}
	for _, tc := range cases {
		v := MakeVOffset(tc.coff, tc.uoff)
		if v.Coff() != tc.coff {
			t.Errorf("Coff: got %d want %d", v.Coff(), tc.coff)
		}
		if v.Uoff() != tc.uoff {
			t.Errorf("Uoff: got %d want %d", v.Uoff(), tc.uoff)
		}
	}
}

func TestVOffsetClampsUoff(t *testing.T) {
	v := MakeVOffset(100, 0x1FFFF) // way past 16 bits
	if v.Uoff() != 0xFFFF {
		t.Errorf("Uoff not clamped: %d", v.Uoff())
	}
	if v.Coff() != 100 {
		t.Errorf("Coff corrupted by uoff clamp: %d", v.Coff())
	}
}

func TestVOffsetNegativeInputs(t *testing.T) {
	v := MakeVOffset(-5, -3)
	if v.Coff() != 0 || v.Uoff() != 0 {
		t.Errorf("negative inputs not normalised: coff=%d uoff=%d", v.Coff(), v.Uoff())
	}
}

func TestVOffsetOrdering(t *testing.T) {
	// Virtual offsets are comparable as uint64. A later block must produce
	// a strictly larger value than any in-block offset of an earlier block.
	earlier := MakeVOffset(100, 0xFFFF)
	later := MakeVOffset(101, 0)
	if earlier >= later {
		t.Errorf("ordering broken: earlier=%x later=%x", earlier, later)
	}
}
