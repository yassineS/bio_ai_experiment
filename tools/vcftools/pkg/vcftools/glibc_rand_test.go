package vcftools

import "testing"

// TestGlibcRandSequence pins the glibc rand() output to the known reference
// sequence (the canonical srand(1) prefix and an srand(42) prefix), captured
// from the system C library. Any deviation breaks --max-indv byte parity.
func TestGlibcRandSequence(t *testing.T) {
	// srand(1); rand() x10 — the textbook glibc rand() sequence.
	wantSeed1 := []int32{
		1804289383, 846930886, 1681692777, 1714636915, 1957747793,
		424238335, 719885386, 1649760492, 596516649, 1189641421,
	}
	g := newGlibcRand(1)
	for i, want := range wantSeed1 {
		if got := g.rand(); got != want {
			t.Fatalf("srand(1) rand()[%d] = %d, want %d", i, got, want)
		}
	}

	// srand(42); rand() x5.
	wantSeed42 := []int32{71876166, 708592740, 1483128881, 907283241, 442951012}
	g = newGlibcRand(42)
	for i, want := range wantSeed42 {
		if got := g.rand(); got != want {
			t.Fatalf("srand(42) rand()[%d] = %d, want %d", i, got, want)
		}
	}
}

// TestGlibcSeedZeroIsOne pins glibc's mapping of seed 0 to seed 1.
func TestGlibcSeedZeroIsOne(t *testing.T) {
	g0 := newGlibcRand(0)
	g1 := newGlibcRand(1)
	for i := 0; i < 20; i++ {
		if a, b := g0.rand(), g1.rand(); a != b {
			t.Fatalf("seed 0 vs 1 diverge at %d: %d != %d", i, a, b)
		}
	}
}

// TestGlibcRandomShuffle pins the std::random_shuffle swap sequence on a small
// index list against the upstream-verified expectation (glibc rand() driving
// the j = rand() % (i+1) two-argument form).
func TestGlibcRandomShuffle(t *testing.T) {
	// For seed 42 over [0..7], the C harness (verified against the real
	// vcftools binary) shuffles to this order.
	idx := []int{0, 1, 2, 3, 4, 5, 6, 7}
	newGlibcRand(42).randomShuffle(idx)
	want := []int{2, 6, 1, 0, 5, 7, 4, 3}
	for i := range want {
		if idx[i] != want[i] {
			t.Fatalf("random_shuffle(seed=42) = %v, want %v", idx, want)
		}
	}
}
