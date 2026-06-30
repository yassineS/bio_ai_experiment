//go:build !linux

package samtools

// drand48PlatformTestSequence holds the first values of macOS/BSD's
// default-seeded drand48() stream (Xi = 0x1234abcd330e). Captured from:
// cc drand.c && ./a.out (macOS arm64).
var drand48PlatformTestSequence = []float64{
	0.39646477376027534,
	0.84048536941142515,
	0.35333609724524351,
	0.44658343479654405,
	0.31869277231188065,
	0.88642843322303122,
}
