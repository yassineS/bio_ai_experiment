package seqtk

// In-tree port of seqtk's `krand` Mersenne Twister 64 (MT19937-64). The
// upstream definition lives in `reference_code/seqtk/seqtk.c` lines
// 300-357 (struct `_krand_t`, `kr_srand0`, `kr_srand`, `kr_rand`, and
// the `kr_drand` macro). The algorithm is the standard Matsumoto /
// Nishimura MT19937-64; this Go port reproduces upstream's byte-for-byte
// output, which `seqtk sample` and `seqtk randbase` rely on (default
// seeds 11 and 0 respectively).
//
// We do not use math/rand because Go's PRNG differs from MT19937-64 and
// would prevent byte parity with upstream's seeded output.

const (
	krNN = 312
	krMM = 156
	krUM = 0xFFFFFFFF80000000 // most significant 33 bits
	krLM = 0x7FFFFFFF         // least significant 31 bits
)

// krand is the MT19937-64 state. The zero value is invalid; use newKrand
// to seed it.
type krand struct {
	mti int
	mt  [krNN]uint64
}

// newKrand returns a fresh MT19937-64 RNG seeded with `seed`, matching
// `kr_srand(seed)` in `seqtk.c`.
func newKrand(seed uint64) *krand {
	k := &krand{}
	k.seed(seed)
	return k
}

// seed mirrors `kr_srand0(seed, kr)` (seqtk.c:315-320).
func (k *krand) seed(seed uint64) {
	k.mt[0] = seed
	for k.mti = 1; k.mti < krNN; k.mti++ {
		prev := k.mt[k.mti-1]
		k.mt[k.mti] = 6364136223846793005*(prev^(prev>>62)) + uint64(k.mti)
	}
}

// rand mirrors `kr_rand(kr)` (seqtk.c:330-355) — one 64-bit MT step
// with the standard tempering and Matsumoto's twist.
func (k *krand) rand() uint64 {
	var mag01 = [2]uint64{0, 0xB5026F5AA96619E9}
	if k.mti >= krNN {
		if k.mti == krNN+1 {
			// Upstream re-seeds with 5489 if rand() is called on an
			// uninitialised generator. newKrand always seeds, so this
			// branch is unreachable in practice — preserved for
			// faithfulness to upstream.
			k.seed(5489)
		}
		var i int
		for i = 0; i < krNN-krMM; i++ {
			x := (k.mt[i] & krUM) | (k.mt[i+1] & krLM)
			k.mt[i] = k.mt[i+krMM] ^ (x >> 1) ^ mag01[x&1]
		}
		for ; i < krNN-1; i++ {
			x := (k.mt[i] & krUM) | (k.mt[i+1] & krLM)
			k.mt[i] = k.mt[i+(krMM-krNN)] ^ (x >> 1) ^ mag01[x&1]
		}
		x := (k.mt[krNN-1] & krUM) | (k.mt[0] & krLM)
		k.mt[krNN-1] = k.mt[krMM-1] ^ (x >> 1) ^ mag01[x&1]
		k.mti = 0
	}
	x := k.mt[k.mti]
	k.mti++
	x ^= (x >> 29) & 0x5555555555555555
	x ^= (x << 17) & 0x71D67FFFEDA60000
	x ^= (x << 37) & 0xFFF7EEE000000000
	x ^= x >> 43
	return x
}

// drand returns a double in [0,1), mirroring the `kr_drand` macro
// at seqtk.c:357: `(kr_rand(kr) >> 11) * (1.0 / 9007199254740992.0)`.
func (k *krand) drand() float64 {
	return float64(k.rand()>>11) * (1.0 / 9007199254740992.0)
}
