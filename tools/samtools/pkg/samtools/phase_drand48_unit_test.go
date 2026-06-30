package samtools

// Binary-free unit tests for the drand48 port and the dump_aln per-read
// routing logic. These pass with the reference_code submodules entirely
// unpopulated (no upstream binary, no exec.Command) — they exercise the
// pure helpers directly.

import "testing"

// TestUnitDrand48Sequence pins the first values of the platform's
// default-seeded drand48() stream. The expected values are supplied by
// drand48PlatformTestSequence (phase_drand48_linux_test.go /
// phase_drand48_other_test.go). Matching these exactly is what makes the
// phase -b split agree with upstream record-for-record.
func TestUnitDrand48Sequence(t *testing.T) {
	d := newDrand48()
	for i, w := range drand48PlatformTestSequence {
		got := d.Float64()
		// Exact float64 equality is required: the LCG is integer-exact and
		// the IEEE-754 division reproduces the platform C library bit-for-bit.
		if got != w {
			t.Errorf("drand48 step %d = %.17g, want %.17g", i, got, w)
		}
	}
}

// TestUnitDrand48IsFlipBranch confirms the first draw matches the platform
// drand48 default state and is < 0.5 (so dump_aln's is_flip is true on the
// first call, swapping haplotype buckets for confident reads).
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
