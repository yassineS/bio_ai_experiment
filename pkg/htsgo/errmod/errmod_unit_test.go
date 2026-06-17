package errmod

import (
	"math"
	"testing"
)

// These binary-free unit tests pin the errmod/gl2cns genotype-likelihood
// computation to values produced by the REAL upstream htslib errmod.c
// (errmod_init(1-0.83) + errmod_cal), captured as exact IEEE-754 float32
// bit patterns. They let us verify the model byte-for-byte without the
// upstream binary, mirroring the live-oracle parity layer.
//
// How the goldens were produced: errmod.c was compiled standalone
// (errmod.c + hts_os.c) with a driver that packs `bases[i]` as
// (qual<<5|strand<<4|base) exactly as Cal expects, calls
// errmod_cal(em, n, m, bases, q), and prints math.Float32bits(q[i]) for
// each of the m*m entries. The model parameter is depcorr = 1-0.83 (the
// value bcftools mpileup and samtools targetcut both use). Reproduce with
// the standalone oracle in the package's parity harness.
//
// The historical low-`-q` het-set divergence (see PARITY_ROADMAP /
// UPSTREAM_BUGS) was caused by accumulating the per-genotype `tmp1` sum in
// float64 and rounding once, whereas upstream declares `tmp1` as a C
// `float` and re-rounds to single precision at every `tmp1 += bsum[k]`
// step. These goldens lock in the upstream-faithful float32 accumulation.

// errmodGolden is one synthetic pileup column with its upstream-exact
// float32 q-matrix (as raw bits). bases are pre-packed.
type errmodGolden struct {
	name  string
	m     int
	bases []uint16
	qbits []uint32
}

func pb(qual, strand, base int) uint16 { return uint16(qual<<5 | strand<<4 | base) }

// errmodGoldens were captured from upstream htslib errmod.c (see file
// doc). Each qbits slice is exactly m*m entries.
var errmodGoldens = []errmodGolden{
	{
		// 4 bases at Q30: 1 A, 1 C, 2 split across strand — a clean het.
		name: "q30_ACGT_one_each",
		m:    5,
		bases: []uint16{
			pb(30, 0, 0), pb(30, 1, 0), pb(30, 0, 1), pb(30, 1, 1),
		},
		qbits: []uint32{
			0x4250e5aa, 0x40884fcd, 0x4268fad6, 0x4268fad6, 0x4268fad6,
			0x40884fcd, 0x4250e5aa, 0x4268fad6, 0x4268fad6, 0x4268fad6,
			0x4268fad6, 0x4268fad6, 0x42d0e5aa, 0x42d0e5aa, 0x42d0e5aa,
			0x4268fad6, 0x4268fad6, 0x42d0e5aa, 0x42d0e5aa, 0x42d0e5aa,
			0x4268fad6, 0x4268fad6, 0x42d0e5aa, 0x42d0e5aa, 0x42d0e5aa,
		},
	},
	{
		// 8 bases at Q30, 4 A / 4 C — deeper het.
		name: "q30_AC_4each",
		m:    5,
		bases: []uint16{
			pb(30, 0, 0), pb(30, 1, 0), pb(30, 0, 0), pb(30, 1, 0),
			pb(30, 0, 1), pb(30, 1, 1), pb(30, 0, 1), pb(30, 1, 1),
		},
		qbits: []uint32{
			0x42b973ad, 0x40b4352c, 0x42d188d9, 0x42d188d9, 0x42d188d9,
			0x40b4352c, 0x42b973ad, 0x42d188d9, 0x42d188d9, 0x42d188d9,
			0x42d188d9, 0x42d188d9, 0x433973ad, 0x433973ad, 0x433973ad,
			0x42d188d9, 0x42d188d9, 0x433973ad, 0x433973ad, 0x433973ad,
			0x42d188d9, 0x42d188d9, 0x433973ad, 0x433973ad, 0x433973ad,
		},
	},
	{
		// Mixed marginal qualities (Q9-Q36) across A and C — exactly the
		// low-quality regime where the het set diverged before the fix.
		name: "mixed_marginal_AC",
		m:    5,
		bases: []uint16{
			pb(10, 0, 0), pb(12, 1, 0), pb(9, 0, 0),
			pb(36, 1, 1), pb(25, 0, 1), pb(33, 1, 1),
		},
		qbits: []uint32{
			0x4298672e, 0x40a1a669, 0x42aa770f, 0x42aa770f, 0x42aa770f,
			0x40a1a669, 0x418a9514, 0x41d2d499, 0x41d2d499, 0x41d2d499,
			0x42aa770f, 0x41d2d499, 0x42bb0c73, 0x42bb0c73, 0x42bb0c73,
			0x42aa770f, 0x41d2d499, 0x42bb0c73, 0x42bb0c73, 0x42bb0c73,
			0x42aa770f, 0x41d2d499, 0x42bb0c73, 0x42bb0c73, 0x42bb0c73,
		},
	},
	{
		// m=4 (samtools targetcut path), 10 bases at Q20, 4 A / 6 C.
		name: "m4_q20_AC",
		m:    4,
		bases: []uint16{
			pb(20, 0, 0), pb(20, 1, 0), pb(20, 0, 0), pb(20, 1, 0),
			pb(20, 0, 1), pb(20, 1, 1), pb(20, 0, 1), pb(20, 1, 1),
			pb(20, 0, 1), pb(20, 1, 1),
		},
		qbits: []uint32{
			0x42a10c3e, 0x40dc3049, 0x42b9216a, 0x42b9216a,
			0x40dc3049, 0x424f76b5, 0x428bdb1d, 0x428bdb1d,
			0x42b9216a, 0x428bdb1d, 0x430463cc, 0x430463cc,
			0x42b9216a, 0x428bdb1d, 0x430463cc, 0x430463cc,
		},
	},
	{
		// 3 bases, three distinct alleles at low quality (Q5/Q5/Q6) —
		// stresses the float32 accumulation across multiple bsum terms.
		name: "lowq_three_alleles",
		m:    5,
		bases: []uint16{
			pb(5, 0, 0), pb(5, 1, 1), pb(6, 0, 2),
		},
		qbits: []uint32{
			0x40813502, 0x40ac00d7, 0x4095dd8d, 0x40e189b3, 0x40e189b3,
			0x40ac00d7, 0x40813502, 0x4095dd8d, 0x40e189b3, 0x40e189b3,
			0x4095dd8d, 0x4095dd8d, 0x40562370, 0x40cb6669, 0x40cb6669,
			0x40e189b3, 0x40e189b3, 0x40cb6669, 0x40b6bdde, 0x40b6bdde,
			0x40e189b3, 0x40e189b3, 0x40cb6669, 0x40b6bdde, 0x40b6bdde,
		},
	},
	{
		// All five alleles present at varied marginal qualities. Several
		// q entries sum three or four bsum terms, so the float32-vs-float64
		// accumulation-width difference is observable: with the old float64
		// accumulation q[0] was 0x430c8688, but upstream's `float tmp1`
		// re-rounds each step and yields 0x430c8689 (locked in below). This
		// is the load-bearing column for the het-set fix.
		name: "all5_alleles_marginal",
		m:    5,
		bases: []uint16{
			pb(7, 0, 0), pb(11, 1, 0), pb(13, 0, 0),
			pb(9, 1, 1), pb(33, 0, 1), pb(28, 1, 1),
			pb(41, 0, 2), pb(19, 1, 2),
			pb(5, 0, 3), pb(50, 1, 3),
			pb(8, 0, 4), pb(22, 1, 4),
		},
		qbits: []uint32{
			0x430c8689, 0x42c6ef13, 0x42cf68dd, 0x42d43685, 0x4304812c,
			0x42c6ef13, 0x42cf5664, 0x4285b232, 0x428a7fd9, 0x42bf4bac,
			0x42cf68dd, 0x4285b232, 0x42d7d02f, 0x42916439, 0x42c6300d,
			0x42d43685, 0x428a7fd9, 0x42916439, 0x42dc9dd6, 0x42cafdb4,
			0x4304812c, 0x42bf4bac, 0x42c6300d, 0x42cafdb4, 0x4308b4d5,
		},
	},
}

// TestUnitCalUpstreamGoldens verifies Cal reproduces the upstream
// errmod_cal q-matrix bit-for-bit (no upstream binary needed).
func TestUnitCalUpstreamGoldens(t *testing.T) {
	em := Init(1.0 - 0.83)
	for _, g := range errmodGoldens {
		g := g
		t.Run(g.name, func(t *testing.T) {
			if len(g.qbits) != g.m*g.m {
				t.Fatalf("golden %s: have %d qbits, want m*m=%d", g.name, len(g.qbits), g.m*g.m)
			}
			bases := make([]uint16, len(g.bases))
			copy(bases, g.bases)
			q := make([]float32, g.m*g.m)
			em.Cal(bases, g.m, q)
			for i := range q {
				got := math.Float32bits(q[i])
				if got != g.qbits[i] {
					t.Errorf("q[%d] = %08x (%v), want upstream %08x (%v)",
						i, got, q[i], g.qbits[i], math.Float32frombits(g.qbits[i]))
				}
			}
		})
	}
}

// TestUnitPhredScaleMatchesUpstreamConstant pins the -10/ln(10)
// coefficient to the exact IEEE-754 double the C compiler folds for
// `-10. / M_LN10`. An untyped Go constant `-10.0 / math.Ln10` rounds in
// arbitrary precision and lands 1 ULP away (0xc0115f2ced384f29); the
// upstream-faithful value is 0xc0115f2ced384f28. This is what makes the
// double-precision beta table bit-identical to upstream.
func TestUnitPhredScaleMatchesUpstreamConstant(t *testing.T) {
	const want = uint64(0xc0115f2ced384f28)
	got := math.Float64bits(phredScale)
	if got != want {
		t.Fatalf("phredScale bits = %016x, want upstream %016x", got, want)
	}
	// And the naive untyped-constant form is the (wrong) 1-ULP neighbour,
	// confirming the fix is load-bearing.
	const naive = -10.0 / math.Ln10
	if math.Float64bits(naive) == want {
		t.Skip("toolchain folds the constant to the upstream value; runtime guard moot")
	}
}

// TestUnitTmp1AccumulationIsFloat32 guards the load-bearing accumulation
// width: the per-genotype likelihood sum must be accumulated in float32
// and re-rounded at every step to match upstream's `float tmp1`. The
// all5_alleles_marginal column is chosen because its hom(A,A) entry sums
// four bsum terms in an order where float64-accumulate-then-round-once
// (the old, divergent behaviour) gives 0x430c8688 while upstream's
// step-rounded float32 gives 0x430c8689. We assert Cal lands on the
// upstream value AND that the two accumulation strategies genuinely
// differ here (so the test cannot silently pass on a non-discriminating
// input).
func TestUnitTmp1AccumulationIsFloat32(t *testing.T) {
	em := Init(1.0 - 0.83)
	var g errmodGolden
	for _, c := range errmodGoldens {
		if c.name == "all5_alleles_marginal" {
			g = c
			break
		}
	}
	if g.name == "" {
		t.Fatal("missing all5_alleles_marginal golden")
	}

	// Recompute aux.bsum exactly as Cal does (sorted, reverse scan).
	n := len(g.bases)
	sorted := make([]uint16, n)
	copy(sorted, g.bases)
	for i := 1; i < n; i++ {
		for j := i; j > 0 && sorted[j] < sorted[j-1]; j-- {
			sorted[j], sorted[j-1] = sorted[j-1], sorted[j]
		}
	}
	var aux callAux
	var w [32]int
	for j := n - 1; j >= 0; j-- {
		b := sorted[j]
		qual := int(b >> 5)
		if qual < 4 {
			qual = 4
		} else if qual > 63 {
			qual = 63
		}
		bs := int(b & 0x1f)
		base := int(b & 0xf)
		aux.bsum[base] += em.fk[w[bs]] * em.beta[qual<<16|n<<8|int(aux.c[base])]
		aux.c[base]++
		w[bs]++
	}

	// Homozygous (A,A): upstream sums bsum over all k != A(=0).
	const jA = 0
	var f32 float32 // step-rounded, like upstream `float tmp1`
	var f64 float64 // round-once, the old divergent behaviour
	for k := 0; k < g.m; k++ {
		if k == jA {
			continue
		}
		f32 = float32(float64(f32) + aux.bsum[k])
		f64 += aux.bsum[k]
	}
	if math.Float32bits(f32) == math.Float32bits(float32(f64)) {
		t.Fatal("golden no longer separates float32 vs float64 accumulation; pick a new column")
	}

	bases := make([]uint16, n)
	copy(bases, g.bases)
	q := make([]float32, g.m*g.m)
	em.Cal(bases, g.m, q)
	got := math.Float32bits(q[jA*g.m+jA])
	if got != math.Float32bits(f32) {
		t.Errorf("hom(A,A) = %08x; want step-rounded float32 %08x (round-once float64 would be %08x)",
			got, math.Float32bits(f32), math.Float32bits(float32(f64)))
	}
	// And it must equal the captured upstream golden bits.
	if got != g.qbits[jA*g.m+jA] {
		t.Errorf("hom(A,A) = %08x, want upstream %08x", got, g.qbits[jA*g.m+jA])
	}
}
