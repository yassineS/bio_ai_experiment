//go:build linux

package samtools

// drand48PlatformTestSequence holds the first values of glibc's default-seeded
// drand48() stream (Xi = 0). Captured from: cc drand.c && ./a.out (Linux).
var drand48PlatformTestSequence = []float64{
	3.907985046680551e-14,
	0.00098539467465030839,
	0.041631001594613082,
	0.17664264254291595,
	0.36460224839060729,
	0.091330612112294318,
}
