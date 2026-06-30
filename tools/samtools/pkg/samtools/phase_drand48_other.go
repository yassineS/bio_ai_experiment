//go:build !linux

package samtools

// drand48DefaultState is BSD libc's initial drand48 state when srand48() is
// never called. BSD (macOS, FreeBSD, etc.) initialises the three unsigned
// short words as {0x330E, 0xABCD, 0x1234}, which in 48-bit little-endian
// form is 0x1234ABCD330E. Upstream samtools compiled on macOS links against
// this libc, so the macOS ARM64 binary in reference_code/samtools/ uses this
// starting state. Our port must start from the same state for byte-identical
// phase -b routing when tested against the macOS oracle. The first Float64()
// returns ≈ 0.3965.
const drand48DefaultState = uint64(0x1234abcd330e)
