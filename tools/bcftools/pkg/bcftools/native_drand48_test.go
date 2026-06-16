package bcftools

import (
	"math"
	"testing"
)

// TestDrand48KnownVectors pins the ported drand48 generator to the canonical
// POSIX/glibc drand48 sequence. The expected values were produced by the system
// C library (srand48(seed); drand48()...) and are identical to htslib's
// hts_srand48/hts_drand48, which is what setGT's random mode consumes. If these
// drift, the random-mode oracle parity is invalid.
func TestDrand48KnownVectors(t *testing.T) {
	cases := []struct {
		seed int64
		want []float64
	}{
		{0, []float64{0.1708280361, 0.7499019805, 0.0963716556, 0.8704652270, 0.5773035068}},
		{1, []float64{0.0416303448, 0.4544924447, 0.8348172182, 0.3359860301, 0.5654894036}},
		{42, []float64{0.7445250001, 0.3427014787, 0.1110852824, 0.4223389580, 0.0811111712}},
	}
	for _, tc := range cases {
		r := newDrand48(tc.seed)
		for i, want := range tc.want {
			got := r.float64()
			if math.Abs(got-want) > 1e-9 {
				t.Errorf("seed=%d draw %d: got %.10f want %.10f", tc.seed, i, got, want)
			}
		}
	}
}

// TestDrand48Range checks the generator stays in [0,1) over many draws.
func TestDrand48Range(t *testing.T) {
	r := newDrand48(12345)
	for i := 0; i < 100000; i++ {
		v := r.float64()
		if v < 0 || v >= 1 {
			t.Fatalf("draw %d out of range: %v", i, v)
		}
	}
}
