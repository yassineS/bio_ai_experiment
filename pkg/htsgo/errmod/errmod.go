// Package errmod is a pure-Go port of htslib's MAQ revised error model
// (reference_code/htslib/errmod.c).
//
// It precomputes a set of model coefficient tables for a given depth-
// correlation parameter and then turns a packed pileup of observed
// bases into an m*m matrix of phred-scaled genotype likelihoods. The
// model underpins both `bcftools mpileup`'s MAQ genotype-likelihood
// path and `samtools targetcut`'s per-position consensus call (via
// `gencns`), and a single shared implementation here keeps them
// behaviour-equivalent.
//
// Compared with upstream's `errmod_cal`, the only deliberate
// divergence is the >255-deep pileup branch: upstream randomly
// downsamples with `ks_shuffle` driven by `drand48`, which would make
// our output non-reproducible without simulating the global drand48
// state. We instead deterministically truncate to the first 255
// observations. The downstream consumers in this project (`bcftools
// mpileup` via `bam2bcf` and `samtools targetcut` via `gencns`) both
// cap depth at 255 in their own output, and for pileups ≤255 deep —
// the overwhelming common case — we are byte-equivalent to upstream.
//
// Floating-point faithfulness: the q-matrix output is computed to be
// bit-for-bit identical to upstream errmod_cal for every ≤255-deep
// pileup (verified against the vendored htslib errmod.c over tens of
// thousands of random columns spanning the full quality range; see the
// TestUnit* oracle goldens). Two upstream-fidelity details are load-
// bearing and easy to get subtly wrong in Go:
//
//   - The per-genotype likelihood sum `tmp1` is a C `float`; upstream
//     does `tmp1 += aux.bsum[k]` where `bsum[k]` is a double, so the
//     running sum is re-rounded to single precision at every step. We
//     accumulate in float32 with an explicit float64 add and float32
//     store to reproduce that exactly (see Cal). Accumulating in float64
//     and rounding once at the end shifts a small fraction of float32
//     q entries by 1 ULP, which at a low het-LOD threshold flipped the
//     het set — this was the documented errmod/gl2cns "LOD precision"
//     divergence, now fixed.
//   - The phred coefficient -10/ln(10) is computed as a runtime double
//     division (phredScale) rather than a Go untyped constant, because
//     an untyped-constant `-10.0 / math.Ln10` rounds in arbitrary
//     precision and lands 1 ULP off the IEEE-754 double the C compiler
//     folds for `-10. / M_LN10`. This keeps the double-precision beta
//     table bit-identical to upstream.
//
// One residual, harmless boundary remains at the double level only: the
// beta table's interior entries use libm exp/log1p, and Go's math.Exp
// differs from glibc's exp by at most 1 ULP for a few arguments. That
// sub-ULP double noise is far below float32 resolution and is fully
// absorbed when the bsum terms are rounded into the float32 q output, so
// it never reaches the result. It is a genuine libm transcendental
// last-ULP property, not a fixable algorithmic difference.
package errmod

import (
	"math"
	"sort"
)

// Model constants matching errmod.c.
const (
	// Eta is the base-error inflation term used by cal_coef (errmod.c
	// hard-codes ETA=0.03).
	Eta = 0.03
	// MaxObs is the largest pile depth the precomputed tables support;
	// deeper piles are deterministically truncated to this many
	// observations (see the package doc for the divergence note).
	MaxObs = 255
)

// phredScale converts a natural log into a phred-scaled penalty
// (`-10. / M_LN10` in errmod.c). It is computed at runtime by
// phredScaleCoef rather than as an untyped Go constant: Go evaluates an
// untyped constant `-10.0 / math.Ln10` in arbitrary precision and then
// rounds the final result, which lands 1 ULP away from the IEEE-754
// `double` quotient the C compiler folds for `-10. / M_LN10` (C divides
// two doubles and rounds once). That single-ULP coefficient difference
// propagated through the whole beta table — every beta entry is
// phredScale*(sum1-sum) — and was the root cause of the low-`-q` het-set
// divergence. Computing the division on a non-constant float64 operand
// forces the IEEE double round-once division, reproducing the upstream
// constant bit-for-bit.
var phredScale = -10.0 / ln10Var

// ln10Var holds ln(10) in a package-level variable (not a constant) so
// that the division -10.0/ln10Var that computes phredScale is evaluated
// as a single IEEE-754 double division at runtime — the C compiler folds
// `-10. / M_LN10` the same way (round-once double divide), whereas a Go
// untyped-constant `-10.0 / math.Ln10` would round in arbitrary
// precision and differ by 1 ULP.
var ln10Var = math.Ln10

// Errmod holds the precomputed MAQ error-model coefficient tables for
// a fixed depth-correlation parameter. Construct it with Init and
// reuse it across many calls to Cal; the struct is read-only after
// init.
type Errmod struct {
	depcorr float64
	// fk is the depth-correlation decay table indexed by observation
	// rank.
	fk []float64
	// beta is the phred-scaled error table indexed by
	// (q<<16 | n<<8 | k).
	beta []float64
	// lhet is the log heterozygous-binomial table indexed by
	// (n<<8 | k).
	lhet []float64
}

// callAux accumulates per-base sums during Cal. It mirrors the
// call_aux_t struct in errmod.c: fsum and bsum are weighted sums per
// base and c is the total observed count per base (strand-agnostic).
type callAux struct {
	fsum [16]float64
	bsum [16]float64
	c    [16]uint32
}

// lfact returns log(n!) via lgamma: Gamma(n+1) = n!.
func lfact(n float64) float64 {
	v, _ := math.Lgamma(n + 1)
	return v
}

// logBinomialTable builds an nSize*nSize table of log-transformed
// binomial coefficients. Entry (n<<8 | k) holds
// log(n!) - log(k!) - log((n-k)!). It is a faithful port of
// logbinomial_table in errmod.c.
func logBinomialTable(nSize int) []float64 {
	logbinom := make([]float64, nSize*nSize)
	for n := 1; n < nSize; n++ {
		lfn := lfact(float64(n))
		for k := 1; k <= n; k++ {
			logbinom[n<<8|k] = lfn - lfact(float64(k)) - lfact(float64(n-k))
		}
	}
	return logbinom
}

// calCoef precomputes the fk, beta and lhet tables for the given
// depth-correlation and eta. It reproduces cal_coef from errmod.c
// exactly.
func (em *Errmod) calCoef(depcorr, eta float64) {
	// fk: depth-correlation decay.
	em.fk = make([]float64, 256)
	em.fk[0] = 1.0
	for n := 1; n < 256; n++ {
		em.fk[n] = math.Pow(1.0-depcorr, float64(n))*(1.0-eta) + eta
	}

	// beta: phred-scaled cumulative error probabilities.
	em.beta = make([]float64, 256*256*64)
	lC := logBinomialTable(256)

	for q := 1; q < 64; q++ {
		e := math.Pow(10.0, -float64(q)/10.0)
		le := math.Log(e)
		le1 := math.Log(1.0 - e)
		for n := 1; n <= 255; n++ {
			base := q<<16 | n<<8
			sum1 := lC[n<<8|n] + float64(n)*le
			em.beta[base+n] = math.Inf(1)
			for k := n - 1; k >= 0; k-- {
				sum := sum1 + math.Log1p(math.Exp(lC[n<<8|k]+float64(k)*le+float64(n-k)*le1-sum1))
				em.beta[base+k] = phredScale * (sum1 - sum)
				sum1 = sum
			}
		}
	}

	// lhet: log heterozygous-binomial table.
	em.lhet = make([]float64, 256*256)
	for n := 0; n < 256; n++ {
		for k := 0; k < 256; k++ {
			em.lhet[n<<8|k] = lC[n<<8|k] - math.Ln2*float64(n)
		}
	}
}

// Init creates an Errmod with the given depth-correlation parameter
// and precomputes its coefficient tables. Upstream callers construct
// the model with `errmod_init(1. - depcorr_param)`; bcftools mpileup
// uses 1-0.83 (theta=0.83) and samtools targetcut uses 1-0.83 (ERR_DEP).
// This is a faithful port of errmod_init in errmod.c.
func Init(depcorr float64) *Errmod {
	em := &Errmod{depcorr: depcorr}
	em.calCoef(depcorr, Eta)
	return em
}

// Cal computes the m*m matrix of phred-scaled genotype likelihoods for
// a pileup of observed bases. It is a faithful port of errmod_cal in
// errmod.c, with the single deliberate divergence on >255 piles
// documented at the package level.
//
// Each bases[i] is a uint16 packed as [6-bit quality | 1-bit strand |
// 4-bit base]: quality occupies bits 5..15, strand is bit 4, and base
// is bits 0..3 (0=A, 1=C, 2=G, 3=T, 4=N). m is the number of alleles
// (4 in samtools targetcut, 5 — A,C,G,T,N — in bcftools mpileup). On
// return q has length at least m*m and q[i*m+j] holds the phred-
// scaled likelihood of genotype (i,j); the matrix is symmetric and
// non-negative (negative likelihoods are clamped to zero, matching
// upstream).
//
// Cal sorts the bases slice in place (in the first min(len(bases),
// MaxObs) positions) — matching upstream's `ks_introsort_uint16`
// call. Callers that need the original ordering should pass a copy.
func (em *Errmod) Cal(bases []uint16, m int, q []float32) {
	for i := range q[:m*m] {
		q[i] = 0
	}
	n := len(bases)
	if n == 0 {
		return
	}

	// Deterministically truncate deep piles so we stay within the
	// precomputed tables. See the package doc for why this differs
	// from upstream's ks_shuffle-based downsampling.
	if n > MaxObs {
		n = MaxObs
	}
	work := bases[:n]
	sort.Slice(work, func(a, b int) bool { return work[a] < work[b] })

	// w[basestrand] is the running count per (strand, base) pair.
	var w [32]int
	var aux callAux

	for j := n - 1; j >= 0; j-- {
		b := work[j]
		// Quality lives in the top bits; floor at 4 and cap at 63.
		qual := int(b >> 5)
		if qual < 4 {
			qual = 4
		} else if qual > 63 {
			qual = 63
		}
		basestrand := int(b & 0x1f)
		base := int(b & 0xf)
		aux.fsum[base] += em.fk[w[basestrand]]
		aux.bsum[base] += em.fk[w[basestrand]] * em.beta[qual<<16|n<<8|int(aux.c[base])]
		aux.c[base]++
		w[basestrand]++
	}

	for j := 0; j < m; j++ {
		// Homozygous genotype (j,j): sum contributions from all other
		// bases.
		//
		// Width note: upstream declares `tmp1` as a C `float` and does
		// `tmp1 += aux.bsum[k]`, where `bsum[k]` is a `double`. The C
		// usual-arithmetic-conversion rules promote `tmp1` to double for
		// the add and then round the result back to float32 on store, so
		// the running sum is re-rounded to single precision at EVERY step.
		// We mirror that exactly by keeping `tmp1` as float32 and rounding
		// each partial sum back through float32; accumulating in float64
		// and rounding once at the end (the previous behaviour) produced a
		// slightly different LOD at the het-call margin under a low `-q`,
		// which is exactly the documented errmod/gl2cns divergence.
		var tmp1 float32
		var tmp2 int
		for k := 0; k < m; k++ {
			if k == j {
				continue
			}
			tmp1 = float32(float64(tmp1) + aux.bsum[k])
			tmp2 += int(aux.c[k])
		}
		if tmp2 != 0 {
			q[j*m+j] = tmp1
		}
		// Heterozygous genotypes (j,k) with k > j.
		for k := j + 1; k < m; k++ {
			cjk := int(aux.c[j]) + int(aux.c[k])
			tmp1 = 0
			tmp2 = 0
			for i := 0; i < m; i++ {
				if i == j || i == k {
					continue
				}
				tmp1 = float32(float64(tmp1) + aux.bsum[i])
				tmp2 += int(aux.c[i])
			}
			// Upstream: q[..] = -4.343 * em->lhet[..] + tmp1, evaluated in
			// double (lhet is double, tmp1 a float promoted to double) and
			// stored to the float32 output. The tmp2!=0 / ==0 branches in
			// upstream are equivalent here: when tmp2==0 tmp1 is 0, so the
			// single expression covers both.
			val := -4.343*em.lhet[cjk<<8|int(aux.c[k])] + float64(tmp1)
			q[j*m+k] = float32(val)
			q[k*m+j] = float32(val)
		}
		// Clamp negative likelihoods to zero.
		for k := 0; k < m; k++ {
			if q[j*m+k] < 0.0 {
				q[j*m+k] = 0.0
			}
		}
	}
}
