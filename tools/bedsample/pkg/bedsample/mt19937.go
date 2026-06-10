package bedsample

// mt19937_64 is a faithful Go port of the C++11 std::mt19937_64 engine, the
// 64-bit Mersenne Twister that upstream `bedtools sample` uses for its
// reservoir replacement decisions (see
// reference_code/bedtools/src/utils/general/Random.cpp, the non-USE_RAND
// branch). Porting the exact engine — rather than reusing Go's math/rand — is
// what makes this port's sampled output byte-for-byte identical to upstream's
// for a given seed.
//
// The constants are the standard MT19937-64 parameters (Matsumoto & Nishimura,
// 2004), which std::mt19937_64 is defined to use, so the generated sequence is
// identical across the standard library and this port.
type mt19937_64 struct {
	state [mtNN]uint64
	index int
}

const (
	mtNN      = 312
	mtMM      = 156
	mtMatrixA = 0xB5026F5AA96619E9
	mtUpper   = 0xFFFFFFFF80000000 // most significant 33 bits
	mtLower   = 0x7FFFFFFF         // least significant 31 bits
)

// newMT19937_64 returns an engine seeded exactly as std::mt19937_64::seed(s)
// does (the C++11 initialization recurrence with multiplier 6364136223846793005).
func newMT19937_64(seed uint64) *mt19937_64 {
	m := &mt19937_64{}
	m.state[0] = seed
	for i := uint64(1); i < mtNN; i++ {
		m.state[i] = 6364136223846793005*(m.state[i-1]^(m.state[i-1]>>62)) + i
	}
	m.index = mtNN
	return m
}

// next returns the next 64-bit value, matching std::mt19937_64::operator().
func (m *mt19937_64) next() uint64 {
	if m.index >= mtNN {
		m.generate()
	}
	x := m.state[m.index]
	m.index++

	// Tempering.
	x ^= (x >> 29) & 0x5555555555555555
	x ^= (x << 17) & 0x71D67FFFEDA60000
	x ^= (x << 37) & 0xFFF7EEE000000000
	x ^= x >> 43
	return x
}

// generate refills the state array (the twist transformation).
func (m *mt19937_64) generate() {
	for i := 0; i < mtNN-mtMM; i++ {
		x := (m.state[i] & mtUpper) | (m.state[i+1] & mtLower)
		m.state[i] = m.state[i+mtMM] ^ (x >> 1) ^ mag01(x)
	}
	for i := mtNN - mtMM; i < mtNN-1; i++ {
		x := (m.state[i] & mtUpper) | (m.state[i+1] & mtLower)
		m.state[i] = m.state[i+(mtMM-mtNN)] ^ (x >> 1) ^ mag01(x)
	}
	x := (m.state[mtNN-1] & mtUpper) | (m.state[0] & mtLower)
	m.state[mtNN-1] = m.state[mtMM-1] ^ (x >> 1) ^ mag01(x)
	m.index = 0
}

// mag01 returns mtMatrixA for odd x and 0 for even x — the conditional XOR
// term of the twist.
func mag01(x uint64) uint64 {
	if x&1 != 0 {
		return mtMatrixA
	}
	return 0
}

// max is the engine's maximum output, std::mt19937_64::max() == 2^64 - 1.
const mtMax = uint64(0xFFFFFFFFFFFFFFFF)

// randRange returns a uniform value in [0, limit) using upstream's exact
// rejection sampling from rand_range (the non-USE_RAND branch):
//
//	max = mt.max() - (mt.max() % limit)
//	do { n = mt(); } while (n >= max);
//	return n % limit
//
// limit must be > 0. Matching this rejection bound exactly is required for
// byte-for-byte parity, since it determines which draws are discarded.
func (m *mt19937_64) randRange(limit uint64) uint64 {
	max := mtMax - (mtMax % limit)
	var n uint64
	for {
		n = m.next()
		if n < max {
			break
		}
	}
	return n % limit
}
