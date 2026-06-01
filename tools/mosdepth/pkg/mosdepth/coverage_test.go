package mosdepth

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

func mustCigar(t *testing.T, s string) sam.Cigar {
	t.Helper()
	c, err := sam.ParseCigar(s)
	if err != nil {
		t.Fatalf("ParseCigar(%q): %v", s, err)
	}
	return c
}

// TestCovAccumAddRecord_BasicCIGAR adds two overlapping 5M reads and
// checks the depth at every position via emit().
func TestCovAccumAddRecord_BasicCIGAR(t *testing.T) {
	a := newCovAccum(30)
	a.addRecord(&sam.Record{Pos: 10, Cigar: mustCigar(t, "5M")}, false) // 9..14 0-based
	a.addRecord(&sam.Record{Pos: 12, Cigar: mustCigar(t, "5M")}, false) // 11..16
	gotDepths := map[int]int32{}
	a.emit(func(pos int, depth int32) { gotDepths[pos] = depth })
	want := map[int]int32{
		// 1-based 10..16 → 0-based 9..15
		9: 1, 10: 1, 11: 2, 12: 2, 13: 2, 14: 1, 15: 1,
		0: 0, 5: 0, 20: 0, 29: 0,
	}
	for p, w := range want {
		if gotDepths[p] != w {
			t.Errorf("pos %d: got %d, want %d", p, gotDepths[p], w)
		}
	}
}

// TestCovAccumAddRecord_IgnoresInsertionsAndClips makes sure I and S do not
// consume reference, while M does. A "2M2I2M" at POS=10 should cover 0-based
// 9..10 and 11..12 (4 bases, two adjacent runs that merge into one).
func TestCovAccumAddRecord_IgnoresInsertionsAndClips(t *testing.T) {
	a := newCovAccum(20)
	a.addRecord(&sam.Record{Pos: 10, Cigar: mustCigar(t, "1S2M2I2M1S")}, false)
	gotDepths := map[int]int32{}
	a.emit(func(pos int, depth int32) { gotDepths[pos] = depth })
	// Reference-consuming runs: 2M @ refPos 9..10, then 2M @ refPos 11..12.
	for _, p := range []int{9, 10, 11, 12} {
		if gotDepths[p] != 1 {
			t.Errorf("pos %d: depth %d, want 1", p, gotDepths[p])
		}
	}
	// Position 8 and 13 should be 0.
	for _, p := range []int{8, 13} {
		if gotDepths[p] != 0 {
			t.Errorf("pos %d: depth %d, want 0", p, gotDepths[p])
		}
	}
}

// TestCovAccumAddRecord_DeletionBreaksRun verifies that a CIGAR "2M2D2M"
// produces two distinct runs (positions on either side of the deletion).
func TestCovAccumAddRecord_DeletionBreaksRun(t *testing.T) {
	a := newCovAccum(20)
	a.addRecord(&sam.Record{Pos: 10, Cigar: mustCigar(t, "2M2D2M")}, false)
	gotDepths := map[int]int32{}
	a.emit(func(pos int, depth int32) { gotDepths[pos] = depth })
	// 2M @ 9..10, then 2D advances ref to 11..12 (no depth), then 2M @ 13..14.
	if gotDepths[9] != 1 || gotDepths[10] != 1 {
		t.Errorf("first run depths: got %v", gotDepths)
	}
	if gotDepths[11] != 0 || gotDepths[12] != 0 {
		t.Errorf("deletion region: got %v", gotDepths)
	}
	if gotDepths[13] != 1 || gotDepths[14] != 1 {
		t.Errorf("second run depths: got %v", gotDepths)
	}
}

// TestCovAccumFastMode treats the whole read as covered: a "2M2D2M" in fast
// mode covers refPos 9..14 (6 bases) instead of two 2-base runs.
func TestCovAccumFastMode(t *testing.T) {
	a := newCovAccum(20)
	a.addRecord(&sam.Record{Pos: 10, Cigar: mustCigar(t, "2M2D2M")}, true)
	gotDepths := map[int]int32{}
	a.emit(func(pos int, depth int32) { gotDepths[pos] = depth })
	for p := 9; p <= 14; p++ {
		if gotDepths[p] != 1 {
			t.Errorf("fast-mode pos %d: depth %d, want 1", p, gotDepths[p])
		}
	}
}

// TestCovAccumEmitRuns confirms that contiguous equal-depth positions are
// collapsed into a single [start, end, depth] tuple by emitRuns.
func TestCovAccumEmitRuns(t *testing.T) {
	a := newCovAccum(20)
	a.addRecord(&sam.Record{Pos: 1, Cigar: mustCigar(t, "5M")}, false) // 0..4
	a.addRecord(&sam.Record{Pos: 3, Cigar: mustCigar(t, "5M")}, false) // 2..6
	type run struct {
		s, e int
		d    int32
	}
	got := []run{}
	a.emitRuns(func(s, e int, d int32) { got = append(got, run{s, e, d}) })
	// Expected layout:
	//  0..2  depth 1
	//  2..5  depth 2
	//  5..7  depth 1
	//  7..20 depth 0
	want := []run{{0, 2, 1}, {2, 5, 2}, {5, 7, 1}, {7, 20, 0}}
	if len(got) != len(want) {
		t.Fatalf("runs: got %v, want %v", got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("run %d: got %v, want %v", i, got[i], want[i])
		}
	}
}

// TestCovAccumRegionStats checks per-region sum, thresholds, min/max.
func TestCovAccumRegionStats(t *testing.T) {
	a := newCovAccum(20)
	// Build a region with depths:
	//   0..2 -> 1, 2..5 -> 2, 5..7 -> 1, 7..20 -> 0
	a.addRecord(&sam.Record{Pos: 1, Cigar: mustCigar(t, "5M")}, false)
	a.addRecord(&sam.Record{Pos: 3, Cigar: mustCigar(t, "5M")}, false)
	sum, perTh, minD, maxD := a.regionStats(0, 7, []int{1, 2}, nil)
	// Bases at depth >= 1: positions 0..6 (7 bases). At >= 2: positions 2..4 (3 bases).
	if perTh[0] != 7 {
		t.Errorf(">=1 threshold: got %d, want 7", perTh[0])
	}
	if perTh[1] != 3 {
		t.Errorf(">=2 threshold: got %d, want 3", perTh[1])
	}
	// Sum across [0,7) = 1*2 + 2*3 + 1*2 = 10.
	if sum != 10 {
		t.Errorf("sum: got %d, want 10", sum)
	}
	if minD != 1 || maxD != 2 {
		t.Errorf("min/max: got %d/%d, want 1/2", minD, maxD)
	}
}

// TestCovAccumRegionMedian verifies the upstream CountStat median
// convention, including the lower-median rule for even base counts.
func TestCovAccumRegionMedian(t *testing.T) {
	// Region [0,10): 7 bases at depth 1, 3 bases at depth 0.
	// Sorted depths: 0,0,0,1,1,1,1,1,1,1. n=10 (even).
	// stop_n = int(0.5 + 5) = 5. cum: depth0 -> 3 (<5), depth1 -> 10 (>=5) => 1.
	a := newCovAccum(10)
	a.add(0, 7) // depth 1 across 0..6
	if m := a.regionMedian(0, 10); m != 1 {
		t.Errorf("median(0,10): got %v, want 1", m)
	}

	// Even split exercising the LOWER-median rule: 5 bases depth 0, 5 depth 2.
	// Sorted: 0,0,0,0,0,2,2,2,2,2. n=10, stop_n=5. cum: depth0 -> 5 (>=5) => 0.
	b := newCovAccum(10)
	b.addSigned(5, 10, 2) // depth 2 across 5..9
	if m := b.regionMedian(0, 10); m != 0 {
		t.Errorf("lower-median(0,10): got %v, want 0 (lower of the two middle values)", m)
	}

	// Odd count: 5 bases, depths 0,0,1,2,2. n=5, stop_n=int(0.5+2.5)=3.
	// cum: depth0 -> 2 (<3), depth1 -> 3 (>=3) => 1.
	c := newCovAccum(5)
	c.add(2, 3)          // depth 1 at pos 2
	c.addSigned(3, 5, 2) // depth 2 at pos 3,4
	if m := c.regionMedian(0, 5); m != 1 {
		t.Errorf("odd median(0,5): got %v, want 1", m)
	}

	// Empty region returns 0.
	if m := a.regionMedian(5, 5); m != 0 {
		t.Errorf("empty median: got %v, want 0", m)
	}
}

// TestCovAccumRegionStats_Empty handles the degenerate case where end <= beg.
func TestCovAccumRegionStats_Empty(t *testing.T) {
	a := newCovAccum(20)
	a.addRecord(&sam.Record{Pos: 1, Cigar: mustCigar(t, "5M")}, false)
	sum, perTh, _, _ := a.regionStats(5, 5, []int{1}, nil)
	if sum != 0 {
		t.Errorf("empty sum: got %d", sum)
	}
	if len(perTh) != 1 || perTh[0] != 0 {
		t.Errorf("empty thresholds: got %v", perTh)
	}
}

// TestCovAccumAddOutOfRange clamps events that fall outside the reference.
func TestCovAccumAddOutOfRange(t *testing.T) {
	a := newCovAccum(5)
	a.add(-2, 10) // becomes 0..5
	got := map[int]int32{}
	a.emit(func(pos int, depth int32) { got[pos] = depth })
	for p := 0; p < 5; p++ {
		if got[p] != 1 {
			t.Errorf("pos %d depth %d, want 1", p, got[p])
		}
	}
}

// TestCovAccumAddRecord_StarPosIgnored — Pos==0 indicates unmapped; no
// events should be inserted.
func TestCovAccumAddRecord_StarPosIgnored(t *testing.T) {
	a := newCovAccum(10)
	a.addRecord(&sam.Record{Pos: 0, Cigar: mustCigar(t, "5M")}, false)
	if len(a.events) != 0 {
		t.Errorf("expected no events, got %d", len(a.events))
	}
}

// TestCovAccumFastMode_NoCIGAR uses an empty CIGAR (the BAM "*") with a
// sequence; fast mode should fall back to len(SEQ).
func TestCovAccumFastMode_NoCIGAR(t *testing.T) {
	a := newCovAccum(20)
	a.addRecord(&sam.Record{Pos: 5, Cigar: nil, Seq: "ACGTAC"}, true)
	gotDepths := map[int]int32{}
	a.emit(func(pos int, depth int32) { gotDepths[pos] = depth })
	for p := 4; p <= 9; p++ {
		if gotDepths[p] != 1 {
			t.Errorf("fastmode no-cigar pos %d depth %d, want 1", p, gotDepths[p])
		}
	}
}

// TestCovAccumEmit_NoRefLen emits up to the largest event position when
// refLen is unknown (<=0).
func TestCovAccumEmit_NoRefLen(t *testing.T) {
	a := newCovAccum(0)
	a.add(2, 5)
	count := 0
	a.emit(func(pos int, depth int32) { count++ })
	if count == 0 {
		t.Errorf("expected emission for inferred ref length")
	}
}
