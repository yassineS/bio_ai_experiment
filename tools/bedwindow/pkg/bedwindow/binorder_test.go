package bedwindow

import "testing"

// TestUnitBinLevelAndNumber pins the UCSC bin placement: records within one
// finest bin share a level/number, a record spanning a coarse range lands on a
// coarser level, and a wider span lands coarser still.
func TestUnitBinLevelAndNumber(t *testing.T) {
	// 100000>>14 == 105000>>14 == 6, so both fall in the finest level, bin
	// offset+6.
	l1, b1 := binLevelAndNumber(100000, 100100)
	l2, b2 := binLevelAndNumber(105000, 105100)
	if l1 != 0 || l2 != 0 {
		t.Fatalf("expected finest level 0, got %d and %d", l1, l2)
	}
	if b1 != b2 {
		t.Fatalf("expected same finest bin number, got %d and %d", b1, b2)
	}
	// A span crossing finest-bin boundaries must land on a coarser level.
	lw, _ := binLevelAndNumber(50000, 200000)
	if lw <= 0 {
		t.Fatalf("expected coarser level for wide span, got %d", lw)
	}
	// An even wider span must land on an even coarser (>=) level.
	lwider, _ := binLevelAndNumber(16000, 600000)
	if lwider < lw {
		t.Fatalf("expected wider span at level >= %d, got %d", lw, lwider)
	}
}

// TestUnitOrderHitsByBin pins the hit ordering: finest level first, then bin
// number ascending, then file order on ties.
func TestUnitOrderHitsByBin(t *testing.T) {
	// File order interleaves a wide record between finest-bin records.
	hits := []*rec{
		{line: "fine1", start: 100000, end: 100100, order: 0},
		{line: "wide1", start: 50000, end: 200000, order: 1},
		{line: "fine2", start: 105000, end: 105100, order: 2},
		{line: "fineDup", start: 100000, end: 100100, order: 3},
	}
	orderHitsByBin(hits)
	want := []string{"fine1", "fine2", "fineDup", "wide1"}
	for i, w := range want {
		if hits[i].line != w {
			t.Fatalf("position %d = %q, want %q (full order: %v)", i, hits[i].line, w, names(hits))
		}
	}
}

func names(hits []*rec) []string {
	out := make([]string, len(hits))
	for i, h := range hits {
		out[i] = h.line
	}
	return out
}
