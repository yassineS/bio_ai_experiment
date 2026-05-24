package errmod

import (
	"math"
	"testing"
)

// packBase builds a bases[] entry the way Cal expects: the 6-bit
// quality occupies bits 5..15, strand is bit 4, and the 4-bit base
// is bits 0..3.
func packBase(qual int, strand int, base int) uint16 {
	return uint16(qual<<5 | strand<<4 | base)
}

const (
	baseA = 0
	baseC = 1
	baseG = 2
	baseT = 3
	baseN = 4
)

func TestInitFkTable(t *testing.T) {
	depcorr := 1.0 - 0.83
	em := Init(depcorr)
	if em.fk[0] != 1.0 {
		t.Fatalf("fk[0] = %v, want 1.0", em.fk[0])
	}
	// fk[n] = (1-depcorr)^n * (1-eta) + eta.
	for _, n := range []int{1, 2, 5, 10, 100, 255} {
		want := math.Pow(1.0-depcorr, float64(n))*(1.0-Eta) + Eta
		if math.Abs(em.fk[n]-want) > 1e-12 {
			t.Errorf("fk[%d] = %v, want %v", n, em.fk[n], want)
		}
	}
	// fk decays monotonically from fk[1] toward eta and stays >= eta.
	for n := 2; n < 256; n++ {
		if em.fk[n] > em.fk[n-1]+1e-15 {
			t.Errorf("fk not monotone non-increasing at n=%d", n)
		}
		if em.fk[n] < Eta-1e-12 {
			t.Errorf("fk[%d] = %v below eta", n, em.fk[n])
		}
	}
}

func TestInitLhetTable(t *testing.T) {
	em := Init(1.0 - 0.83)
	// lhet[n<<8|k] = logBinom(n,k) - ln2*n.
	lbinom := func(n, k int) float64 {
		lf := func(x int) float64 {
			v, _ := math.Lgamma(float64(x) + 1)
			return v
		}
		return lf(n) - lf(k) - lf(n-k)
	}
	cases := []struct{ n, k int }{{0, 0}, {1, 0}, {1, 1}, {4, 2}, {10, 3}, {20, 10}}
	for _, c := range cases {
		want := lbinom(c.n, c.k) - math.Ln2*float64(c.n)
		got := em.lhet[c.n<<8|c.k]
		if math.Abs(got-want) > 1e-9 {
			t.Errorf("lhet[%d,%d] = %v, want %v", c.n, c.k, got, want)
		}
	}
	// A fair 50/50 split (k = n/2) is the most likely het outcome, so its
	// lhet entry is the largest over k for fixed n.
	n := 20
	best := math.Inf(-1)
	bestK := -1
	for k := 0; k <= n; k++ {
		if em.lhet[n<<8|k] > best {
			best = em.lhet[n<<8|k]
			bestK = k
		}
	}
	if bestK != n/2 {
		t.Errorf("lhet peak at k=%d, want %d", bestK, n/2)
	}
}

func TestInitBetaTable(t *testing.T) {
	em := Init(1.0 - 0.83)
	// beta[q<<16|n<<8|n] is +Inf for every q,n (the all-error endpoint).
	for _, q := range []int{1, 10, 40, 63} {
		for _, n := range []int{1, 5, 50, 255} {
			if !math.IsInf(em.beta[q<<16|n<<8|n], 1) {
				t.Errorf("beta[q=%d,n=%d,k=%d] not +Inf", q, n, n)
			}
		}
	}
	// beta is a phred-scaled positive penalty: for k < n it must be finite
	// and non-negative.
	for _, n := range []int{4, 20, 100} {
		for q := 1; q < 64; q++ {
			v := em.beta[q<<16|n<<8|0]
			if math.IsInf(v, 0) || math.IsNaN(v) || v < -1e-9 {
				t.Errorf("beta[q=%d,n=%d,k=0] = %v invalid", q, n, v)
			}
		}
		// beta[...|0] is the phred penalty for the first observed copy of a
		// base at quality q; a higher-quality base carries more weight, so
		// the term is non-decreasing in q.
		for q := 2; q < 64; q++ {
			if em.beta[q<<16|n<<8|0] < em.beta[(q-1)<<16|n<<8|0]-1e-9 {
				t.Errorf("beta[k=0] not non-decreasing in q at q=%d n=%d", q, n)
			}
		}
		// Within a fixed (q,n), each successive observed copy costs no less
		// than the previous one: beta is non-decreasing in k for k < n.
		for q := 1; q < 64; q++ {
			for k := 1; k < n; k++ {
				if em.beta[q<<16|n<<8|k] < em.beta[q<<16|n<<8|(k-1)]-1e-9 {
					t.Errorf("beta not non-decreasing in k at q=%d n=%d k=%d", q, n, k)
				}
			}
		}
	}
}

// bestGenotype returns the (i,j) index pair (i<=j) of the smallest
// score in the m*m phred-scaled likelihood matrix.
func bestGenotype(q []float32, m int) (int, int) {
	bi, bj := 0, 0
	best := float32(math.MaxFloat32)
	for i := 0; i < m; i++ {
		for j := i; j < m; j++ {
			if q[i*m+j] < best {
				best = q[i*m+j]
				bi, bj = i, j
			}
		}
	}
	return bi, bj
}

func TestCalAllRef(t *testing.T) {
	em := Init(1.0 - 0.83)
	m := 5
	q := make([]float32, m*m)
	bases := make([]uint16, 0, 30)
	for i := 0; i < 30; i++ {
		bases = append(bases, packBase(30, i&1, baseA))
	}
	em.Cal(bases, m, q)
	bi, bj := bestGenotype(q, m)
	if bi != baseA || bj != baseA {
		t.Fatalf("all-REF best genotype = (%d,%d), want (A,A)", bi, bj)
	}
	if q[baseA*m+baseA] != 0 {
		t.Errorf("hom-REF score = %v, want 0", q[baseA*m+baseA])
	}
	// Every non-REF genotype must be clearly worse.
	for j := 0; j < m; j++ {
		if j == baseA {
			continue
		}
		if q[baseA*m+j] <= 1.0 {
			t.Errorf("het (A,%d) score %v not clearly worse", j, q[baseA*m+j])
		}
		if q[j*m+j] <= 1.0 {
			t.Errorf("hom (%d,%d) score %v not clearly worse", j, j, q[j*m+j])
		}
	}
}

func TestCalAllAlt(t *testing.T) {
	em := Init(1.0 - 0.83)
	m := 5
	q := make([]float32, m*m)
	bases := make([]uint16, 0, 30)
	for i := 0; i < 30; i++ {
		bases = append(bases, packBase(30, i&1, baseG))
	}
	em.Cal(bases, m, q)
	bi, bj := bestGenotype(q, m)
	if bi != baseG || bj != baseG {
		t.Fatalf("all-ALT best genotype = (%d,%d), want (G,G)", bi, bj)
	}
}

func TestCalHet5050(t *testing.T) {
	em := Init(1.0 - 0.83)
	m := 5
	q := make([]float32, m*m)
	bases := make([]uint16, 0, 40)
	for i := 0; i < 40; i++ {
		b := baseA
		if i%2 == 1 {
			b = baseC
		}
		bases = append(bases, packBase(30, i&1, b))
	}
	em.Cal(bases, m, q)
	bi, bj := bestGenotype(q, m)
	if !(bi == baseA && bj == baseC) {
		t.Fatalf("50/50 best genotype = (%d,%d), want (A,C)", bi, bj)
	}
	// The het beats both relevant homozygous calls.
	if q[baseA*m+baseC] >= q[baseA*m+baseA] || q[baseA*m+baseC] >= q[baseC*m+baseC] {
		t.Errorf("het score %v not better than homs (%v,%v)",
			q[baseA*m+baseC], q[baseA*m+baseA], q[baseC*m+baseC])
	}
}

func TestCalQualityMonotone(t *testing.T) {
	em := Init(1.0 - 0.83)
	m := 5
	// An all-REF pile: higher base quality => more confident, so the
	// competing wrong genotypes get monotonically worse (larger) scores.
	prevHom := float32(-1)
	prevHet := float32(-1)
	for _, qual := range []int{10, 20, 30, 40} {
		q := make([]float32, m*m)
		bases := make([]uint16, 0, 20)
		for i := 0; i < 20; i++ {
			bases = append(bases, packBase(qual, i&1, baseA))
		}
		em.Cal(bases, m, q)
		hom := q[baseC*m+baseC]
		het := q[baseA*m+baseC]
		if prevHom >= 0 && hom < prevHom-1e-3 {
			t.Errorf("hom-ALT score not monotone in quality at Q%d", qual)
		}
		if prevHet >= 0 && het < prevHet-1e-3 {
			t.Errorf("het score not monotone in quality at Q%d", qual)
		}
		prevHom, prevHet = hom, het
	}
}

func TestCalEmptyPile(t *testing.T) {
	em := Init(1.0 - 0.83)
	m := 5
	q := make([]float32, m*m)
	for i := range q {
		q[i] = 7 // garbage that must be zeroed.
	}
	em.Cal(nil, m, q)
	for i, v := range q {
		if v != 0 {
			t.Fatalf("q[%d] = %v, want 0 for empty pile", i, v)
		}
	}
}

func TestCalSymmetric(t *testing.T) {
	em := Init(1.0 - 0.83)
	m := 5
	q := make([]float32, m*m)
	bases := []uint16{
		packBase(25, 0, baseA), packBase(25, 1, baseC),
		packBase(30, 0, baseA), packBase(20, 1, baseG),
		packBase(35, 0, baseC), packBase(28, 1, baseA),
	}
	em.Cal(bases, m, q)
	for i := 0; i < m; i++ {
		for j := 0; j < m; j++ {
			if q[i*m+j] != q[j*m+i] {
				t.Errorf("q not symmetric at (%d,%d): %v vs %v",
					i, j, q[i*m+j], q[j*m+i])
			}
			if q[i*m+j] < 0 {
				t.Errorf("q[%d,%d] = %v negative (clamp failed)", i, j, q[i*m+j])
			}
		}
	}
	_ = baseN
	_ = baseT
}

// TestCalDownsampleTruncates exercises the >MaxObs path: a 600-deep
// pile is deterministically truncated to the first 255 observations,
// the call still completes, and the dominant genotype still wins.
func TestCalDownsampleTruncates(t *testing.T) {
	em := Init(1.0 - 0.83)
	m := 5
	bases := make([]uint16, 0, 600)
	for i := 0; i < 600; i++ {
		bases = append(bases, packBase(30, i&1, baseT))
	}
	q := make([]float32, m*m)
	em.Cal(bases, m, q)
	bi, bj := bestGenotype(q, m)
	if bi != baseT || bj != baseT {
		t.Fatalf("downsampled all-T best genotype = (%d,%d), want (T,T)", bi, bj)
	}
	for i, v := range q {
		if math.IsNaN(float64(v)) || math.IsInf(float64(v), 0) {
			t.Fatalf("q[%d] = %v not finite after downsample", i, v)
		}
	}
}

// TestCalDownsampleDeterministic guards the documented divergence from
// upstream: truncation to the first MaxObs observations is fully
// deterministic, so two identical inputs yield identical outputs even
// when n > MaxObs. (Upstream's ks_shuffle + drand48 path is not byte-
// reproducible across calls; ours is.)
func TestCalDownsampleDeterministic(t *testing.T) {
	em := Init(1.0 - 0.83)
	m := 5
	mk := func() []uint16 {
		b := make([]uint16, 0, 400)
		for i := 0; i < 400; i++ {
			b = append(b, packBase(30, i&1, baseA))
		}
		return b
	}
	q1 := make([]float32, m*m)
	q2 := make([]float32, m*m)
	em.Cal(mk(), m, q1)
	em.Cal(mk(), m, q2)
	for i := range q1 {
		if q1[i] != q2[i] {
			t.Fatalf("downsample not deterministic at q[%d]: %v vs %v", i, q1[i], q2[i])
		}
	}
}

// TestCalEquivalentToFirstMaxObs documents the truncation contract: a
// pile of MaxObs+k bases produces the same result as just the first
// MaxObs of them. (The downstream consumers — gencns in samtools
// targetcut and bam2bcf in bcftools mpileup — both cap depth at 255
// in their output, so this is the documented divergence from
// upstream's randomised downsampling.)
func TestCalEquivalentToFirstMaxObs(t *testing.T) {
	em := Init(1.0 - 0.83)
	m := 4
	mk := func(n int) []uint16 {
		b := make([]uint16, 0, n)
		for i := 0; i < n; i++ {
			b = append(b, packBase(30, i&1, baseA))
		}
		return b
	}
	qLong := make([]float32, m*m)
	qShort := make([]float32, m*m)
	em.Cal(mk(MaxObs+50), m, qLong)
	em.Cal(mk(MaxObs), m, qShort)
	for i := range qLong {
		if qLong[i] != qShort[i] {
			t.Fatalf("truncate-first contract broken at q[%d]: %v vs %v",
				i, qLong[i], qShort[i])
		}
	}
}
