// Revised MAQ error model, ported from htslib's errmod.c.
//
// errmod precomputes a set of model coefficient tables for a given
// depth-correlation parameter and then turns a packed pileup of observed
// bases into an m*m matrix of phred-scaled genotype likelihoods. It is the
// genotype-likelihood machinery underneath `bcftools mpileup`'s MAQ model.

package bcftools

import (
	"math"
	"math/rand"
	"sort"
)

// Model constants matching errmod.c.
const (
	errmodEta   = 0.03      // base-error inflation term used by cal_coef
	errmodLn2   = math.Ln2  // M_LN2
	errmodLn10  = math.Ln10 // M_LN10
	errmodPhred = -10.0 / errmodLn10
)

// Errmod holds the precomputed MAQ error-model coefficient tables for a
// fixed depth-correlation parameter. Construct it with ErrmodInit and reuse
// it across many calls to ErrmodCal; the struct is read-only after init.
type Errmod struct {
	depcorr float64
	// fk is the depth-correlation decay table indexed by observation rank.
	fk []float64
	// beta is the phred-scaled error table indexed by (q<<16 | n<<8 | k).
	beta []float64
	// lhet is the log heterozygous-binomial table indexed by (n<<8 | k).
	lhet []float64
}

// errmodCallAux accumulates the per-base sums during ErrmodCal. It mirrors
// the call_aux_t struct in errmod.c: fsum and bsum are weighted sums per
// base and c is the total observed count per base (strand-agnostic).
type errmodCallAux struct {
	fsum [16]float64
	bsum [16]float64
	c    [16]uint32
}

// lfact returns log(n!) via the log-gamma function: Gamma(n+1) = n!.
func lfact(n float64) float64 {
	v, _ := math.Lgamma(n + 1)
	return v
}

// logBinomialTable builds an nSize*nSize table of log-transformed binomial
// coefficients. Entry (n<<8 | k) holds log(n!) - log(k!) - log((n-k)!).
// It is a faithful port of logbinomial_table in errmod.c.
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
// depth-correlation and eta. It reproduces cal_coef from errmod.c exactly.
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
				em.beta[base+k] = errmodPhred * (sum1 - sum)
				sum1 = sum
			}
		}
	}

	// lhet: log heterozygous-binomial table.
	em.lhet = make([]float64, 256*256)
	for n := 0; n < 256; n++ {
		for k := 0; k < 256; k++ {
			em.lhet[n<<8|k] = lC[n<<8|k] - errmodLn2*float64(n)
		}
	}
}

// ErrmodInit creates an Errmod with the given depth-correlation parameter
// and precomputes its coefficient tables. bcftools mpileup constructs the
// model with depcorr = 1 - 0.83. It mirrors errmod_init in errmod.c.
func ErrmodInit(depcorr float64) *Errmod {
	em := &Errmod{depcorr: depcorr}
	em.calCoef(depcorr, errmodEta)
	return em
}

// errmodMaxObs is the largest pile depth the precomputed tables support;
// deeper piles are randomly downsampled to this many observations.
const errmodMaxObs = 255

// ErrmodCal computes the m*m matrix of phred-scaled genotype likelihoods for
// a pileup of n observed bases. It is a faithful port of errmod_cal.
//
// Each bases[i] is a uint16 packed as [6-bit quality | 1-bit strand |
// 4-bit base]: quality occupies bits 5..15, strand is bit 4, and base is
// bits 0..3 (0=A,1=C,2=G,3=T,4=N). m is the number of alleles (5 for
// A,C,G,T,N). On return q has length m*m and q[i*m+j] holds the
// phred-scaled likelihood of genotype (i,j); the matrix is symmetric.
//
// ErrmodCal may reorder and (when n > 255) shuffle the bases slice in
// place, matching upstream behaviour. rng supplies the downsampling
// randomness; pass nil to use the default global source.
func (em *Errmod) ErrmodCal(n, m int, bases []uint16, q []float32, rng *rand.Rand) {
	for i := range q[:m*m] {
		q[i] = 0
	}
	if n == 0 {
		return
	}

	// Downsample deep piles so we stay within the precomputed tables.
	if n > errmodMaxObs {
		shuffleUint16(bases[:n], rng)
		n = errmodMaxObs
	}
	sort.Slice(bases[:n], func(a, b int) bool { return bases[a] < bases[b] })

	// w[basestrand] is the running count per (strand, base) combination.
	var w [32]int
	var aux errmodCallAux

	for j := n - 1; j >= 0; j-- {
		b := bases[j]
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
		// Homozygous genotype (j,j): sum contributions from all other bases.
		var tmp1 float64
		var tmp2 int
		for k := 0; k < m; k++ {
			if k == j {
				continue
			}
			tmp1 += aux.bsum[k]
			tmp2 += int(aux.c[k])
		}
		if tmp2 != 0 {
			q[j*m+j] = float32(tmp1)
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
				tmp1 += aux.bsum[i]
				tmp2 += int(aux.c[i])
			}
			val := -4.343*em.lhet[cjk<<8|int(aux.c[k])] + tmp1
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

// shuffleUint16 performs an in-place Fisher-Yates shuffle, matching the
// ks_shuffle macro used by errmod.c to randomly downsample deep piles.
func shuffleUint16(a []uint16, rng *rand.Rand) {
	for i := len(a) - 1; i > 0; i-- {
		var j int
		if rng != nil {
			j = rng.Intn(i + 1)
		} else {
			j = rand.Intn(i + 1)
		}
		a[i], a[j] = a[j], a[i]
	}
}
