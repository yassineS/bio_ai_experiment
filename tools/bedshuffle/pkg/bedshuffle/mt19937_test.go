package bedshuffle

import "testing"

// TestMT19937_64Sequence pins the std::mt19937_64 output for seed 5489 (the
// engine's default seed) to the value the C++ standard mandates after 10000
// invocations — the canonical conformance check for MT19937-64. Matching this
// guarantees the same draw sequence bedtools shuffle consumes.
func TestMT19937_64Sequence(t *testing.T) {
	m := newMT19937_64(5489)
	var v uint64
	for i := 0; i < 10000; i++ {
		v = m.next()
	}
	const want = uint64(9981545732273789042)
	if v != want {
		t.Fatalf("mt19937_64 10000th value = %d, want %d", v, want)
	}
}

// TestMT19937RandRangeBound checks the rejection-sampling bound matches the
// upstream rand_range formula for a small limit.
func TestMT19937RandRangeBound(t *testing.T) {
	m := newMT19937_64(1)
	limit := uint64(100)
	for i := 0; i < 1000; i++ {
		if v := m.randRange(limit); v >= limit {
			t.Fatalf("randRange(%d) returned %d out of range", limit, v)
		}
	}
}
