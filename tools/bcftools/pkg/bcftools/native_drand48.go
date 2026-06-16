// Port of htslib's deterministic drand48 PRNG (hts_srand48 / hts_drand48).
//
// htslib routes all of its randomness through hts_srand48()/hts_drand48(),
// which on every platform produce the SAME POSIX-defined sequence: where the
// host libc lacks a deterministic drand48 it falls back to the FreeBSD
// implementation amalgamated in htslib's os/rand.c, and on glibc it uses
// srand48_deterministic(). All of these share the classic 48-bit linear
// congruential generator
//
//	X_{n+1} = (a * X_n + c) mod 2^48,   a = 0x5DEECE66D, c = 0xB
//
// seeded so that the low 16 bits are 0x330E and the upper 32 bits are the
// supplied seed. Because setGT's `-t r` mode seeds this generator from a fixed
// default (0) or the `-s INT` option and nothing else (no time()/getpid()),
// the sequence is fully reproducible and the native plugin matches upstream
// byte-for-byte. (This corrects the earlier code comment that claimed the RNG
// "cannot be matched byte-for-byte" — it can, and is, validated by the setGT
// random-mode oracle.)
package bcftools

// drand48 is the deterministic 48-bit LCG used by htslib's hts_drand48.
type drand48 struct {
	x uint64 // current 48-bit state
}

// drand48 LCG constants (POSIX / SVID, matching htslib's RAND48_MULT/ADD).
const (
	drand48Mult = 0x5DEECE66D
	drand48Add  = 0xB
	drand48Mask = (uint64(1) << 48) - 1
)

// newDrand48 returns a generator seeded exactly as htslib's hts_srand48(seed):
// the low 16 bits of the state are the constant 0x330E and the upper 32 bits
// are the (32-bit-truncated) seed.
func newDrand48(seed int64) *drand48 {
	lo := uint16(0x330E)
	mid := uint16(seed)
	hi := uint16(seed >> 16)
	state := uint64(lo) | uint64(mid)<<16 | uint64(hi)<<32
	return &drand48{x: state & drand48Mask}
}

// next advances the LCG one step and returns the new 48-bit state.
func (r *drand48) next() uint64 {
	r.x = (drand48Mult*r.x + drand48Add) & drand48Mask
	return r.x
}

// float64 returns the next pseudo-random double in [0, 1), identical to
// htslib's hts_drand48(): the 48-bit state scaled by 2^-48.
func (r *drand48) float64() float64 {
	return float64(r.next()) / float64(uint64(1)<<48)
}
