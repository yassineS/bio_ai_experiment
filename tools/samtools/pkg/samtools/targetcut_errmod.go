package samtools

import (
	"math"
	"sort"
)

// errModel is a Go port of htslib's errmod_t (errmod.c). It is the
// MAQ revised error model used by samtools targetcut (cut_target.c
// calls `errmod_cal` to score the per-position base likelihoods).
//
// Implementation notes:
//   - Only the bits the upstream cut_target.c caller exercises are
//     supported: m=4 (ACGT), n up to 255 (256-cell precomputed tables).
//   - The pre-computed coefficient tables (fk, beta, lhet) are sized
//     identically to upstream so byte-equivalent results are
//     guaranteed: fk[0..255], beta[64*256*256], lhet[256*256].
//   - All math is float64 internally; the public Calc method writes
//     the m*m likelihood matrix as float32 to match the upstream
//     `float *q` API and to keep the gencns sort byte-identical.
//
// The error model itself is described in Heng Li, R Durbin, et al.,
// "Mapping short DNA sequencing reads ...", 2008; the constants
// (depcorr, eta=0.03) come from the upstream defaults.
type errModel struct {
	depCorr float64
	fk      []float64 // [256]
	beta    []float64 // [64*256*256] indexed by (qual<<16 | n<<8 | k)
	lhet    []float64 // [256*256] indexed by (n<<8 | k)
}

// newErrModel constructs an errModel for the given depCorr. Upstream
// samtools targetcut calls `errmod_init(1. - 0.83)` so callers should
// pass errModDepCorr (= 0.17).
func newErrModel(depCorr float64) *errModel {
	em := &errModel{depCorr: depCorr}
	em.calcCoef(depCorr, 0.03)
	return em
}

// calcCoef precomputes fk, beta, and lhet. Faithful port of cal_coef
// in errmod.c.
func (em *errModel) calcCoef(depCorr, eta float64) {
	// fk[n] = (1-depcorr)^n * (1-eta) + eta, with fk[0] = 1.0.
	em.fk = make([]float64, 256)
	em.fk[0] = 1.0
	for n := 1; n < 256; n++ {
		em.fk[n] = math.Pow(1.0-depCorr, float64(n))*(1.0-eta) + eta
	}

	lC := logBinomialTable(256) // [256*256], lC[n<<8|k] = log(C(n,k))

	em.beta = make([]float64, 64*256*256)
	for q := 1; q < 64; q++ {
		e := math.Pow(10.0, -float64(q)/10.0)
		le := math.Log(e)
		le1 := math.Log(1.0 - e)
		for n := 1; n <= 255; n++ {
			base := q<<16 | n<<8
			// sum1 = lC[n<<8|n] + n*le; beta[base+n] = +Inf
			sum1 := lC[n<<8|n] + float64(n)*le
			em.beta[base+n] = math.Inf(1)
			var sum float64
			for k := n - 1; k >= 0; k-- {
				// sum = sum1 + log1p(exp(lC[n<<8|k] + k*le + (n-k)*le1 - sum1))
				inner := lC[n<<8|k] + float64(k)*le + float64(n-k)*le1 - sum1
				sum = sum1 + math.Log1p(math.Exp(inner))
				em.beta[base+k] = -10.0 / math.Ln10 * (sum1 - sum)
				sum1 = sum
			}
		}
	}

	em.lhet = make([]float64, 256*256)
	for n := 0; n < 256; n++ {
		for k := 0; k < 256; k++ {
			em.lhet[n<<8|k] = lC[n<<8|k] - math.Ln2*float64(n)
		}
	}
}

// logBinomialTable returns a (size*size) flat table where entry
// [n<<8|k] is log(C(n, k)) for 1<=n<size, 1<=k<=n; other cells are 0.
// Matches errmod.c's logbinomial_table.
func logBinomialTable(size int) []float64 {
	out := make([]float64, size*size)
	for n := 1; n < size; n++ {
		lfn := lFact(n)
		for k := 1; k <= n; k++ {
			out[n<<8|k] = lfn - lFact(k) - lFact(n-k)
		}
	}
	return out
}

// lFact(n) = log(n!) computed via lgamma(n+1). Mirrors the errmod.c
// `lfact(n) lgamma(n+1)` macro.
func lFact(n int) float64 {
	lg, _ := math.Lgamma(float64(n + 1))
	return lg
}

// calc fills q[] with the m*m phred-scaled likelihood matrix for the
// observed packed bases. q must have length >= m*m. m is the allele
// alphabet size (always 4 for nt ACGT in cut_target.c). Faithful port
// of errmod_cal in errmod.c.
//
// bases entries are packed as upstream: bits 0..3 = base, bit 4 =
// strand (1 = reverse), bits 5..10 = quality. Quality is clamped to
// [4, 63] inside this function (matching the upstream `qual = b>>5 <
// 4? 4 : b>>5`).
func (em *errModel) calc(bases []uint16, m int, q []float32) {
	// Zero q[0..m*m).
	for i := 0; i < m*m; i++ {
		q[i] = 0
	}
	n := len(bases)
	if n == 0 {
		return
	}

	// Upstream randomly downsamples to 255 (ks_shuffle then truncate)
	// when n > 255. We use a deterministic truncation here: the
	// downstream consumer (targetcut.gencns) caps depth at 255 in the
	// output anyway, and the shuffle's seed in upstream depends on
	// drand48 state which would make our output non-reproducible. For
	// pileups ≤255 deep (the overwhelming common case) we are
	// byte-equivalent to upstream.
	work := bases
	if n > 255 {
		work = make([]uint16, 255)
		copy(work, bases[:255])
		n = 255
	}
	// ks_introsort by uint16 value: upstream sorts ascending.
	sortable := make([]uint16, n)
	copy(sortable, work)
	sort.Slice(sortable, func(i, j int) bool { return sortable[i] < sortable[j] })

	var fsum, bsum [16]float64
	var c [16]int
	var w [32]int

	for j := n - 1; j >= 0; j-- {
		b := sortable[j]
		qual := int(b >> 5)
		if qual < 4 {
			qual = 4
		}
		if qual > 63 {
			qual = 63
		}
		baseStrand := int(b & 0x1f)
		base := int(b & 0xf)
		fk := em.fk[w[baseStrand]]
		fsum[base] += fk
		bsum[base] += fk * em.beta[qual<<16|n<<8|c[base]]
		c[base]++
		w[baseStrand]++
	}

	// Likelihood matrix.
	for j := 0; j < m; j++ {
		// homozygous q[j*m+j]
		var tmp1, tmp3 float64
		var tmp2 int
		for k := 0; k < m; k++ {
			if k == j {
				continue
			}
			tmp1 += bsum[k]
			tmp2 += c[k]
			tmp3 += fsum[k]
		}
		if tmp2 != 0 {
			q[j*m+j] = float32(tmp1)
		}
		_ = tmp3 // upstream computes tmp3 but never reads it in the homozygous branch
		// heterozygous q[j*m+k] for k > j
		for k := j + 1; k < m; k++ {
			cjk := c[j] + c[k]
			var h1, h3 float64
			var h2 int
			for i := 0; i < m; i++ {
				if i == j || i == k {
					continue
				}
				h1 += bsum[i]
				h2 += c[i]
				h3 += fsum[i]
			}
			_ = h3
			het := -4.343 * em.lhet[cjk<<8|c[k]]
			if h2 != 0 {
				v := float32(het + h1)
				q[j*m+k] = v
				q[k*m+j] = v
			} else {
				v := float32(het)
				q[j*m+k] = v
				q[k*m+j] = v
			}
		}
		// Clamp negatives.
		for k := 0; k < m; k++ {
			if q[j*m+k] < 0.0 {
				q[j*m+k] = 0.0
			}
		}
	}
}
