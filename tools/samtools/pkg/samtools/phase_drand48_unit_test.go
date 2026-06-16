package samtools

// Binary-free unit tests for the drand48 port and the dump_aln per-read
// routing logic. These pass with the reference_code submodules entirely
// unpopulated (no upstream binary, no exec.Command) — they exercise the
// pure helpers directly.

import "testing"

// TestUnitDrand48Sequence pins the first values of glibc's default-seeded
// drand48() stream. These constants were captured from the system C
// library (cc drand.c; for(i<6) printf drand48()). Matching them is what
// makes the phase -b split agree with upstream record-for-record.
func TestUnitDrand48Sequence(t *testing.T) {
	d := newDrand48()
	want := []float64{
		3.907985046680551e-14,
		0.00098539467465030839,
		0.041631001594613082,
		0.17664264254291595,
		0.36460224839060729,
		0.091330612112294318,
	}
	for i, w := range want {
		got := d.Float64()
		// Exact float64 equality is required: the LCG is integer-exact and
		// the IEEE-754 division reproduces glibc's __erand48_r bit-for-bit.
		if got != w {
			t.Errorf("drand48 step %d = %.17g, want %.17g", i, got, w)
		}
	}
}

// TestUnitDrand48IsFlipBranch confirms the first draw is < 0.5 (so
// dump_aln's is_flip is true on the first call) and the routing parity it
// implies. The very first drand48 ≈ 3.9e-14 < 0.5, which is exactly why
// upstream's confident haplotype reads get their 0<->1 bucket swapped on
// the first dump_aln call.
func TestUnitDrand48IsFlipBranch(t *testing.T) {
	d := newDrand48()
	if first := d.Float64(); first >= 0.5 {
		t.Fatalf("first drand48 = %.17g, expected < 0.5 (is_flip true)", first)
	}
}

// TestUnitDrand48Deterministic verifies two fresh generators produce the
// same stream — i.e. the default seed is fixed (Xi = 0), matching an
// unseeded upstream samtools process.
func TestUnitDrand48Deterministic(t *testing.T) {
	a, b := newDrand48(), newDrand48()
	for i := 0; i < 100; i++ {
		if x, y := a.Float64(), b.Float64(); x != y {
			t.Fatalf("step %d: a=%.17g b=%.17g (generators diverged)", i, x, y)
		}
	}
}
