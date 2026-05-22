package bcftools

import (
	"math"
	"math/rand"
	"testing"
)

// packBase builds a bases[] entry the way errmod_cal expects: the 6-bit
// quality occupies bits 5..15, strand is bit 4, and the 4-bit base is
// bits 0..3.
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

func TestErrmodInitFkTable(t *testing.T) {
	depcorr := 1.0 - 0.83
	em := ErrmodInit(depcorr)
	if em.fk[0] != 1.0 {
		t.Fatalf("fk[0] = %v, want 1.0", em.fk[0])
	}
	// fk[n] = (1-depcorr)^n * (1-eta) + eta.
	for _, n := range []int{1, 2, 5, 10, 100, 255} {
		want := math.Pow(1.0-depcorr, float64(n))*(1.0-errmodEta) + errmodEta
		if math.Abs(em.fk[n]-want) > 1e-12 {
			t.Errorf("fk[%d] = %v, want %v", n, em.fk[n], want)
		}
	}
	// fk decays monotonically from fk[1] toward eta and stays >= eta.
	for n := 2; n < 256; n++ {
		if em.fk[n] > em.fk[n-1]+1e-15 {
			t.Errorf("fk not monotone non-increasing at n=%d", n)
		}
		if em.fk[n] < errmodEta-1e-12 {
			t.Errorf("fk[%d] = %v below eta", n, em.fk[n])
		}
	}
}

func TestErrmodInitLhetTable(t *testing.T) {
	em := ErrmodInit(1.0 - 0.83)
	// lhet[n<<8|k] = logBinom(n,k) - ln2*n. Compute logBinom by hand.
	lbinom := func(n, k int) float64 {
		lf := func(x int) float64 {
			v, _ := math.Lgamma(float64(x) + 1)
			return v
		}
		if k == 0 || n == 0 {
			return lf(n) - lf(k) - lf(n-k)
		}
		return lf(n) - lf(k) - lf(n-k)
	}
	cases := []struct{ n, k int }{{0, 0}, {1, 0}, {1, 1}, {4, 2}, {10, 3}, {20, 10}}
	for _, c := range cases {
		want := lbinom(c.n, c.k) - errmodLn2*float64(c.n)
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

func TestErrmodInitBetaTable(t *testing.T) {
	em := ErrmodInit(1.0 - 0.83)
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

// best returns the index of the smallest (best) phred score in the m*m
// genotype matrix, decoded as the (i,j) pair with i <= j.
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

func TestErrmodCalAllRef(t *testing.T) {
	em := ErrmodInit(1.0 - 0.83)
	m := 5
	q := make([]float32, m*m)
	bases := make([]uint16, 0, 30)
	for i := 0; i < 30; i++ {
		strand := i & 1
		bases = append(bases, packBase(30, strand, baseA))
	}
	em.ErrmodCal(len(bases), m, bases, q, nil)
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

func TestErrmodCalAllAlt(t *testing.T) {
	em := ErrmodInit(1.0 - 0.83)
	m := 5
	q := make([]float32, m*m)
	bases := make([]uint16, 0, 30)
	for i := 0; i < 30; i++ {
		bases = append(bases, packBase(30, i&1, baseG))
	}
	em.ErrmodCal(len(bases), m, bases, q, nil)
	bi, bj := bestGenotype(q, m)
	if bi != baseG || bj != baseG {
		t.Fatalf("all-ALT best genotype = (%d,%d), want (G,G)", bi, bj)
	}
}

func TestErrmodCalHet5050(t *testing.T) {
	em := ErrmodInit(1.0 - 0.83)
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
	em.ErrmodCal(len(bases), m, bases, q, nil)
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

func TestErrmodCalQualityMonotone(t *testing.T) {
	em := ErrmodInit(1.0 - 0.83)
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
		em.ErrmodCal(len(bases), m, bases, q, nil)
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

func TestErrmodCalEmptyPile(t *testing.T) {
	em := ErrmodInit(1.0 - 0.83)
	m := 5
	q := make([]float32, m*m)
	for i := range q {
		q[i] = 7 // garbage that must be zeroed.
	}
	em.ErrmodCal(0, m, nil, q, nil)
	for i, v := range q {
		if v != 0 {
			t.Fatalf("q[%d] = %v, want 0 for empty pile", i, v)
		}
	}
}

func TestErrmodCalSymmetric(t *testing.T) {
	em := ErrmodInit(1.0 - 0.83)
	m := 5
	q := make([]float32, m*m)
	bases := []uint16{
		packBase(25, 0, baseA), packBase(25, 1, baseC),
		packBase(30, 0, baseA), packBase(20, 1, baseG),
		packBase(35, 0, baseC), packBase(28, 1, baseA),
	}
	em.ErrmodCal(len(bases), m, bases, q, nil)
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
}

func TestErrmodCalDownsample(t *testing.T) {
	em := ErrmodInit(1.0 - 0.83)
	m := 5
	rng := rand.New(rand.NewSource(1))
	// 600 observations > 255 forces the shuffle/downsample path.
	bases := make([]uint16, 0, 600)
	for i := 0; i < 600; i++ {
		bases = append(bases, packBase(30, i&1, baseT))
	}
	q := make([]float32, m*m)
	em.ErrmodCal(len(bases), m, bases, q, rng)
	// Path must complete and still pick the dominant genotype.
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

func TestErrmodCalDownsampleStable(t *testing.T) {
	em := ErrmodInit(1.0 - 0.83)
	m := 5
	// A homogeneous pile downsampled with the same seed twice must give
	// identical results (deterministic given the rng).
	mk := func() []uint16 {
		b := make([]uint16, 0, 400)
		for i := 0; i < 400; i++ {
			b = append(b, packBase(30, i&1, baseA))
		}
		return b
	}
	q1 := make([]float32, m*m)
	q2 := make([]float32, m*m)
	em.ErrmodCal(400, m, mk(), q1, rand.New(rand.NewSource(42)))
	em.ErrmodCal(400, m, mk(), q2, rand.New(rand.NewSource(42)))
	for i := range q1 {
		if q1[i] != q2[i] {
			t.Fatalf("downsample not stable at q[%d]: %v vs %v", i, q1[i], q2[i])
		}
	}
}
