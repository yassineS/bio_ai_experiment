package vcftools

// glibcRand is a faithful Go port of glibc's TYPE_3 additive-feedback
// pseudo-random generator — the algorithm behind the C library's rand() /
// random() functions (the default GNU/Linux PRNG). Upstream vcftools'
// --max-indv downsampling (variant_file_filters.cpp,
// filter_individuals_randomly) draws from srand()/rand() through
// std::random_shuffle, so reproducing that exact selection requires the exact
// generator rather than Go's math/rand.
//
// glibc's random() uses a 31-word state and the trinomial x^31 + x^3 + 1
// (degree 31, separation 3). srandom(seed) initialises the state with a
// minimal-standard LCG warm-up, sets the feedback pointers, then steps the
// generator 310 times (10*DEG) to discard the seeding transient. Each rand()
// returns the high 31 bits of an additive-feedback step. This implementation
// matches glibc bit-for-bit; see glibc stdlib/random_r.c and stdlib/rand.c.
type glibcRand struct {
	// table is the 31-word additive-feedback state (glibc's "r" array).
	table [glibcDEG3]int32
	// fptr and rptr are the front and rear feedback pointers.
	fptr int
	rptr int
}

const (
	// glibcDEG3 is the degree of the TYPE_3 trinomial (31 state words).
	glibcDEG3 = 31
	// glibcSEP3 is the separation (the x^3 term) of the TYPE_3 trinomial.
	glibcSEP3 = 3
)

// newGlibcRand returns a generator seeded exactly as glibc's srand(seed) /
// srandom(seed). A seed of 0 is treated as 1, mirroring glibc.
func newGlibcRand(seed uint32) *glibcRand {
	g := &glibcRand{}
	g.srandom(seed)
	return g
}

// srandom initialises the state, matching glibc's __srandom_r for a TYPE_3
// generator (the default state size rand()/random() use).
func (g *glibcRand) srandom(seed uint32) {
	if seed == 0 {
		seed = 1 // glibc maps seed 0 to 1.
	}
	g.table[0] = int32(seed)
	for i := 1; i < glibcDEG3; i++ {
		// glibc's seeding recurrence is the minimal-standard LCG
		//   state[i] = (16807 * state[i-1]) % 2147483647
		// implemented via Schrage's hi/lo split to avoid 32-bit overflow,
		// with the same negative-result correction glibc applies.
		hi := int64(g.table[i-1]) / 127773
		lo := int64(g.table[i-1]) % 127773
		word := 16807*lo - 2836*hi
		if word < 0 {
			word += 2147483647
		}
		g.table[i] = int32(word)
	}
	g.fptr = glibcSEP3
	g.rptr = 0
	// Discard the first 10*DEG outputs to flush the seeding transient.
	for i := 0; i < 10*glibcDEG3; i++ {
		g.next()
	}
}

// next performs one additive-feedback step and returns the 31-bit output,
// matching glibc's __random_r for a TYPE_3 generator.
func (g *glibcRand) next() int32 {
	// The addition wraps modulo 2^32 (glibc uses uint32 arithmetic).
	val := uint32(g.table[g.fptr]) + uint32(g.table[g.rptr])
	g.table[g.fptr] = int32(val)
	// glibc returns (val >> 1) & 0x7fffffff.
	result := int32((val >> 1) & 0x7fffffff)

	g.fptr++
	if g.fptr >= glibcDEG3 {
		g.fptr = 0
	}
	g.rptr++
	if g.rptr >= glibcDEG3 {
		g.rptr = 0
	}
	return result
}

// rand returns the next pseudo-random integer in [0, 2^31), exactly as the C
// library's rand() does on glibc-based systems.
func (g *glibcRand) rand() int32 {
	return g.next()
}

// randomShuffle reorders idx in place using glibc rand() exactly as
// libstdc++'s std::random_shuffle(first, last) does:
//
//	for i := 1; i < n; i++ {
//	    j := rand() % (i+1)
//	    swap(idx[i], idx[j])
//	}
//
// This is the two-argument form (no custom RNG functor), which is what
// upstream vcftools calls. Reproducing the exact swap sequence is what makes
// the downsample selection match a glibc/libstdc++ run for a given seed.
func (g *glibcRand) randomShuffle(idx []int) {
	for i := 1; i < len(idx); i++ {
		j := int(g.rand()) % (i + 1)
		if i != j {
			idx[i], idx[j] = idx[j], idx[i]
		}
	}
}
