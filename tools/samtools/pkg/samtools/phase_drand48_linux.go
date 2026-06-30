//go:build linux

package samtools

// drand48DefaultState is glibc's initial drand48 state when srand48() is
// never called (Xi = 0). Upstream samtools on Linux never seeds drand48, so
// this matches its random stream from the first call. The first Float64()
// returns 0xB / 2^48 ≈ 3.907985e-14.
const drand48DefaultState = uint64(0)
