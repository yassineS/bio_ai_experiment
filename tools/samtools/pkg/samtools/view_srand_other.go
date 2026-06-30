//go:build !linux

package samtools

// platformSrandRandTestCases holds reference values for platformSrandRand
// verified against macOS libc's Park-Miller LCG (srand/rand on Darwin).
var platformSrandRandTestCases = []struct {
	seed uint32
	want uint32
}{
	{1, 16807},
	{2, 33614},
	{3, 50421},
	{42, 705894},
}

// platformSrandRand reproduces `srand(seed); return rand();` under the
// Park-Miller LCG used by macOS (and BSD-derived) libc. The formula is:
//
//	result = (seed * 16807) % (2^31 - 1)
//
// Upstream samtools compiled on macOS links against this libc, so the macOS
// ARM64 binary in reference_code/samtools/ uses this seed transform. Our
// port must match it for byte-identical subsample decisions when running
// against a macOS-built oracle (sam_view.c:1307-1311).
func platformSrandRand(seed uint32) uint32 {
	const modulus = uint64(1<<31 - 1) // 2147483647
	return uint32((uint64(seed) * 16807) % modulus)
}
