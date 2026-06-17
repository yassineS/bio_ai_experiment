package mosdepth

import (
	"strconv"
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
	sum, perTh, minD, maxD, _ := a.regionStats(0, 7, []int{1, 2}, nil, 0)
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

// TestCovAccumRegionStats_IMean verifies the regions mean is the per-base
// divide-then-accumulate Σ(depth_i/L) that upstream mosdepth's imean computes,
// NOT (Σ depth_i)/L. The two differ by float rounding on boundary values: a
// region summing to 945 over width 280 has exact quotient 3.375, but the
// per-base accumulation lands just below 3.375 (so %.2f prints 3.37, matching
// upstream, where sum/width would print 3.38). We reconstruct that exact case.
func TestCovAccumRegionStats_IMean(t *testing.T) {
	const width = 280
	a := newCovAccum(width)
	// Lay down coverage summing to 945 over [0,280): 105 bases at depth 4
	// (420) + 175 bases at depth 3 (525) = 945; 945/280 = 3.375 exactly.
	a.addRecord(&sam.Record{Pos: 1, Cigar: mustCigar(t, "280M")}, false) // depth 1 over all
	a.addRecord(&sam.Record{Pos: 1, Cigar: mustCigar(t, "280M")}, false) // depth 2
	a.addRecord(&sam.Record{Pos: 1, Cigar: mustCigar(t, "280M")}, false) // depth 3
	a.addRecord(&sam.Record{Pos: 1, Cigar: mustCigar(t, "105M")}, false) // +1 over first 105 -> depth 4
	sum, _, _, _, fmean := a.regionStats(0, width, nil, nil, float64(width))
	if sum != 945 {
		t.Fatalf("sum: got %d, want 945", sum)
	}
	// sum/width is exactly 3.375; the per-base imean must be strictly below it.
	if fmean >= float64(sum)/float64(width) {
		t.Errorf("imean fmean=%.17g must be < sum/width=%.17g (per-base accumulation)", fmean, float64(sum)/float64(width))
	}
	if got := strconv.FormatFloat(fmean, 'f', 2, 64); got != "3.37" {
		t.Errorf("formatted imean = %s, want 3.37 (sum/width would give 3.38)", got)
	}
}

// TestCovAccumRegionStats_Empty handles the degenerate case where end <= beg.
func TestCovAccumRegionStats_Empty(t *testing.T) {
	a := newCovAccum(20)
	a.addRecord(&sam.Record{Pos: 1, Cigar: mustCigar(t, "5M")}, false)
	sum, perTh, _, _, _ := a.regionStats(5, 5, []int{1}, nil, 0)
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
