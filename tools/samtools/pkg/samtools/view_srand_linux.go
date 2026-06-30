//go:build linux

package samtools

// platformSrandRandTestCases holds reference values for platformSrandRand
// verified against glibc's TYPE_3 generator (the C values from a Linux build).
var platformSrandRandTestCases = []struct {
	seed uint32
	want uint32
}{
	{1, 1804289383},
	{2, 1505335290},
	{3, 1205554746},
	{42, 71876166},
}

// platformSrandRand reproduces `srand(seed); return rand();` under glibc's
// default TYPE_3 additive-feedback generator (stdlib/random_r.c). Upstream
// samtools on Linux links against glibc, so this matches its subsample seed
// transform (sam_view.c:1307-1311) byte-for-byte.
func platformSrandRand(seed uint32) uint32 {
	if seed == 0 {
		seed = 1
	}
	var r [344]int32
	r[0] = int32(seed)
	for i := 1; i < 31; i++ {
		// r[i] = (16807 * r[i-1]) % 2147483647, via Schrage's method to
		// avoid 32-bit signed overflow (matches glibc's int arithmetic).
		hi := r[i-1] / 127773
		lo := r[i-1] % 127773
		w := 16807*lo - 2836*hi
		if w < 0 {
			w += 2147483647
		}
		r[i] = w
	}
	for i := 31; i < 34; i++ {
		r[i] = r[i-31]
	}
	for i := 34; i < 344; i++ {
		r[i] = r[i-31] + r[i-3]
	}
	val := r[344-31] + r[344-3]
	return (uint32(val) >> 1) & 0x7fffffff
}
