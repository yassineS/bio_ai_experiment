package seqtk

// In-tree port of glibc's drand48 / srand48. Used by Randbase to
// reproduce upstream `seqtk randbase`'s output byte-for-byte. Upstream
// calls drand48() without an explicit srand48, so it relies on glibc's
// default initial state, which is 0 (i.e. `X0 = 0`). See
// `man drand48` and glibc's stdlib/drand48-iter.c.
//
// The algorithm is a 48-bit linear congruential generator:
//
//	X_{n+1} = (a * X_n + c) mod 2^48
//	with a = 0x5DEECE66D, c = 0xB
//
// drand48 returns the upper 48 bits of X / 2^48 as a double in [0, 1).

const (
	drand48A    uint64 = 0x5DEECE66D
	drand48C    uint64 = 0xB
	drand48Mask uint64 = (1 << 48) - 1
)

// drand48State holds a glibc-compatible 48-bit LCG state. The zero
// value represents the same state upstream `seqtk randbase` sees on a
// fresh process — no srand48 call required.
type drand48State struct {
	x uint64
}

// next returns the next drand48 draw in [0, 1).
func (d *drand48State) next() float64 {
	d.x = (drand48A*d.x + drand48C) & drand48Mask
	return float64(d.x) / float64(uint64(1)<<48)
}
