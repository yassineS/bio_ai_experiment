package samtools

// drand48 is a byte-exact Go port of the POSIX/glibc drand48(3) random
// number generator, used so that the `samtools phase -b` per-haplotype
// BAM split routes reads into the .0 / .1 / .chimera buckets identically
// to upstream samtools phase.c.
//
// Why this matters: phase.c's dump_aln (phase.c:361) draws drand48()
// twice — once for the per-call `is_flip` shuffle that may swap a
// confidently-phased read between bucket 0 and bucket 1 (phase.c:386),
// and once per evidence-less read to route it 50/50 (phase.c:388).
// Upstream never calls srand48(), so it runs from glibc's *default* seed
// state Xi = 0, which makes the whole output stream deterministic and
// therefore byte-reproducible. Matching that exact sequence is the only
// way the split BAMs agree with upstream record-for-record (not merely
// up to a 0<->1 relabelling).
//
// The generator is a 48-bit linear congruential generator:
//
//	Xi(n+1) = (a * Xi(n) + c) mod 2^48
//	drand48 = Xi(n+1) / 2^48
//
// with a = 0x5DEECE66D and c = 0xB, the constants mandated by POSIX.
// glibc's __erand48_r builds the IEEE-754 double straight from the 48
// state bits, which is exactly Xi/2^48; the float64 division below is
// bit-identical to that construction for all 2^48 states.
type drand48 struct {
	state uint64
}

// drand48 LCG constants (POSIX / glibc).
const (
	drand48Mult = 0x5DEECE66D
	drand48Add  = 0xB
	drand48Mask = (uint64(1) << 48) - 1
)

// newDrand48 returns a drand48 generator seeded to the platform's default
// startup state, matching an upstream samtools process that never calls
// srand48(). The initial state is platform-specific (see
// phase_drand48_linux.go and phase_drand48_other.go):
//
//   - Linux (glibc): Xi = 0; first Float64() ≈ 3.907985e-14.
//   - macOS/BSD: Xi = 0x1234abcd330e; first Float64() ≈ 0.3965.
func newDrand48() *drand48 { return &drand48{state: drand48DefaultState} }

// Float64 advances the generator one step and returns the next value in
// [0, 1), matching glibc drand48(). It is named Float64 so a *drand48
// satisfies the same one-method shape used by the math/rand-based paths.
func (d *drand48) Float64() float64 {
	d.state = (drand48Mult*d.state + drand48Add) & drand48Mask
	return float64(d.state) / float64(uint64(1)<<48)
}

// phaseRNG abstracts the single Float64() entry point used by the
// phase BAM-split routing so the upstream-faithful drand48 generator
// and the legacy math/rand generator are interchangeable.
type phaseRNG interface {
	Float64() float64
}
