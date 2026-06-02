package bcftools

import (
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Faithful Go port of reference_code/bcftools/{ccall.c, em.c,
// prob1.c} + the kfunc.c primitives (kfLgamma, kfGammap/q, kfBetai)
// needed by the consensus caller. Everything is bundled in one
// file because the parallel-agent harness was racing per-file
// Writes during this session.
//
// Implements `bcftools call -c` for single-group (no -G, no -C
// trio) mpileup-style inputs with INFO/I16 + INFO/QS. The numerical
// kernels are byte-equivalent transcriptions of upstream so
// floating-point results match the C build on IEEE-754 hardware.
//
// Synthetic PL-only fixtures (no I16) fall through to the existing
// v1 heuristic in call.go::callVariant — no test regressions there.

// =============================================================
//  kfunc — kfLgamma / kfGammap / kfGammaq / kfBetai  (kfErfc lives
//  in bam2bcf.go and is shared with the bias-test path).
// =============================================================

func kfLgamma(z float64) float64 {
	x := 0.0
	x += 0.1659470187408462e-06 / (z + 7)
	x += 0.9934937113930748e-05 / (z + 6)
	x -= 0.1385710331296526 / (z + 5)
	x += 12.50734324009056 / (z + 4)
	x -= 176.6150291498386 / (z + 3)
	x += 771.3234287757674 / (z + 2)
	x -= 1259.139216722289 / (z + 1)
	x += 676.5203681218835 / z
	x += 0.9999999999995183
	return math.Log(x) - 5.58106146679532777 - z + (z-0.5)*math.Log(z+6.5)
}

const (
	kfGammaEPS = 1e-14
	kfTiny     = 1e-290
)

func kfGammapSeries(s, z float64) float64 {
	sum := 1.0
	x := 1.0
	for k := 1; k < 100; k++ {
		x *= z / (s + float64(k))
		sum += x
		if x/sum < kfGammaEPS {
			break
		}
	}
	return math.Exp(s*math.Log(z) - z - kfLgamma(s+1.) + math.Log(sum))
}

func kfGammaqCF(s, z float64) float64 {
	f := 1 + z - s
	C := f
	D := 0.0
	for j := 1; j < 100; j++ {
		a := float64(j) * (s - float64(j))
		b := float64((j<<1)+1) + z - s
		D = b + a*D
		if D < kfTiny {
			D = kfTiny
		}
		C = b + a/C
		if C < kfTiny {
			C = kfTiny
		}
		D = 1 / D
		d := C * D
		f *= d
		if math.Abs(d-1) < kfGammaEPS {
			break
		}
	}
	return math.Exp(s*math.Log(z) - z - kfLgamma(s) - math.Log(f))
}

func kfGammap(s, z float64) float64 {
	if z <= 1 || z < s {
		return kfGammapSeries(s, z)
	}
	return 1 - kfGammaqCF(s, z)
}

func kfGammaq(s, z float64) float64 {
	if z <= 1 || z < s {
		return 1 - kfGammapSeries(s, z)
	}
	return kfGammaqCF(s, z)
}

func kfBetaiAux(a, b, x float64) float64 {
	if x == 0 {
		return 0
	}
	if x == 1 {
		return 1
	}
	f := 1.0
	C := f
	D := 0.0
	for j := 1; j < 200; j++ {
		m := j >> 1
		var aa float64
		if j&1 != 0 {
			aa = -(a + float64(m)) * (a + b + float64(m)) * x / ((a + float64(2*m)) * (a + float64(2*m+1)))
		} else {
			aa = float64(m) * (b - float64(m)) * x / ((a + float64(2*m-1)) * (a + float64(2*m)))
		}
		D = 1 + aa*D
		if D < kfTiny {
			D = kfTiny
		}
		C = 1 + aa/C
		if C < kfTiny {
			C = kfTiny
		}
		D = 1 / D
		d := C * D
		f *= d
		if math.Abs(d-1) < kfGammaEPS {
			break
		}
	}
	return math.Exp(kfLgamma(a+b)-kfLgamma(a)-kfLgamma(b)+a*math.Log(x)+b*math.Log(1-x)) / a / f
}

func kfBetai(a, b, x float64) float64 {
	if x < (a+1)/(a+b+2) {
		return kfBetaiAux(a, b, x)
	}
	return 1 - kfBetaiAux(b, a, 1-x)
}

// =============================================================
//  em.c port — single-locus EM for AF, HWE chi^2
// =============================================================

const (
	emIterMax = 50
	emIterTry = 10
	emEPS     = 1e-5
)

func emEstFreq(n int, pdg []float64) float64 {
	var gcnt [3]int
	for i := 0; i < n; i++ {
		p := pdg[i*3:]
		if p[0] != 1 || p[1] != 1 || p[2] != 1 {
			which := 0
			if p[0] <= p[1] {
				which = 1
			}
			if p[which] <= p[2] {
				which = 2
			}
			gcnt[which]++
		}
	}
	tmp := gcnt[0] + gcnt[1] + gcnt[2]
	if tmp == 0 {
		return -1
	}
	return (0.5*float64(gcnt[1]) + float64(gcnt[2])) / float64(tmp)
}

func emFreqIter(f *float64, pdg []float64, beg, end int) float64 {
	f0 := *f
	var f3 [3]float64
	f3[0] = (1 - f0) * (1 - f0)
	f3[1] = 2 * f0 * (1 - f0)
	f3[2] = f0 * f0
	acc := 0.0
	for i := beg; i < end; i++ {
		p := pdg[i*3:]
		acc += (p[1]*f3[1] + 2*p[2]*f3[2]) / (p[0]*f3[0] + p[1]*f3[1] + p[2]*f3[2])
	}
	acc /= float64(end-beg) * 2
	err := math.Abs(acc - *f)
	*f = acc
	return err
}

func emProb1(f float64, pdg []float64, beg, end int) float64 {
	if f < 0 || f > 1 {
		return 1e300
	}
	p := 1.0
	l := 0.0
	f3 := [3]float64{(1 - f) * (1 - f), 2 * f * (1 - f), f * f}
	for i := beg; i < end; i++ {
		pd := pdg[i*3:]
		p *= pd[0]*f3[0] + pd[1]*f3[1] + pd[2]*f3[2]
		if p < 1e-200 {
			l -= math.Log(p)
			p = 1
		}
	}
	return l - math.Log(p)
}

func kminBrent(fn func(float64) float64, a, b float64, tol float64, xmin *float64) float64 {
	const (
		ITMAX = 100
		CGOLD = 0.3819660112501051
		ZEPS  = 1e-10
	)
	var d, e float64
	x := b
	w := b
	v := b
	fx := fn(x)
	fw := fx
	fv := fx
	aa := a
	bb := b
	if aa > bb {
		aa, bb = bb, aa
	}
	for iter := 0; iter < ITMAX; iter++ {
		xm := 0.5 * (aa + bb)
		tol1 := tol*math.Abs(x) + ZEPS
		tol2 := 2 * tol1
		if math.Abs(x-xm) <= tol2-0.5*(bb-aa) {
			*xmin = x
			return fx
		}
		useGold := true
		if math.Abs(e) > tol1 {
			r := (x - w) * (fx - fv)
			q := (x - v) * (fx - fw)
			p := (x-v)*q - (x-w)*r
			q = 2 * (q - r)
			if q > 0 {
				p = -p
			}
			q = math.Abs(q)
			etemp := e
			e = d
			if math.Abs(p) >= math.Abs(0.5*q*etemp) || p <= q*(aa-x) || p >= q*(bb-x) {
				useGold = true
			} else {
				d = p / q
				u := x + d
				if u-aa < tol2 || bb-u < tol2 {
					d = math.Copysign(tol1, xm-x)
				}
				useGold = false
			}
		}
		if useGold {
			if x >= xm {
				e = aa - x
			} else {
				e = bb - x
			}
			d = CGOLD * e
		}
		var u float64
		if math.Abs(d) >= tol1 {
			u = x + d
		} else {
			u = x + math.Copysign(tol1, d)
		}
		fu := fn(u)
		if fu <= fx {
			if u >= x {
				aa = x
			} else {
				bb = x
			}
			v = w
			fv = fw
			w = x
			fw = fx
			x = u
			fx = fu
		} else {
			if u < x {
				aa = u
			} else {
				bb = u
			}
			if fu <= fw || w == x {
				v = w
				fv = fw
				w = u
				fw = fu
			} else if fu <= fv || v == x || v == w {
				v = u
				fv = fu
			}
		}
	}
	*xmin = x
	return fx
}

func emFreqML(f0 float64, beg, end int, pdg []float64) float64 {
	f := f0
	i := 0
	for ; i < emIterTry; i++ {
		if emFreqIter(&f, pdg, beg, end) < emEPS {
			break
		}
	}
	if i == emIterTry {
		var x float64
		start := f
		if f0 == f {
			start = 0.5 * f0
		}
		kminBrent(func(v float64) float64 { return emProb1(v, pdg, beg, end) }, start, f, emEPS, &x)
		f = x
	}
	return f
}

func emG3Iter(g []float64, pdg []float64, beg, end int) float64 {
	var gg [3]float64
	for i := beg; i < end; i++ {
		var tmp [3]float64
		p := pdg[i*3:]
		tmp[0] = p[0] * g[0]
		tmp[1] = p[1] * g[1]
		tmp[2] = p[2] * g[2]
		sum := (tmp[0] + tmp[1] + tmp[2]) * float64(end-beg)
		gg[0] += tmp[0] / sum
		gg[1] += tmp[1] / sum
		gg[2] += tmp[2] / sum
	}
	err := math.Abs(gg[0] - g[0])
	if d := math.Abs(gg[1] - g[1]); d > err {
		err = d
	}
	if d := math.Abs(gg[2] - g[2]); d > err {
		err = d
	}
	g[0] = gg[0]
	g[1] = gg[1]
	g[2] = gg[2]
	return err
}

func bcfEM1(nAllele, nSample, n1 int, flag int, pdg []float64, x []float64) int {
	if nAllele < 2 {
		return -1
	}
	if n1 < 0 || n1 > nSample {
		n1 = 0
	}
	if flag&(1<<7) != 0 {
		flag |= 7 << 5
	}
	if flag&(0xf<<1) != 0 {
		flag |= 0xf << 1
	}
	n := nSample
	for i := 0; i < 10; i++ {
		x[i] = -1
	}
	x[0] = emEstFreq(n, pdg)
	if x[0] < 0 {
		return -1
	}
	x[0] = emFreqML(x[0], 0, n, pdg)

	if flag&(0xf<<1|3<<8) != 0 {
		g := x[1:4]
		var f3 [3]float64
		f3[0] = (1 - x[0]) * (1 - x[0])
		g[0] = f3[0]
		f3[1] = 2 * x[0] * (1 - x[0])
		g[1] = f3[1]
		f3[2] = x[0] * x[0]
		g[2] = f3[2]
		for i := 0; i < emIterMax; i++ {
			if emG3Iter(g, pdg, 0, n) < emEPS {
				break
			}
		}
		r := 1.0
		for i := 0; i < n; i++ {
			p := pdg[i*3:]
			r *= (p[0]*g[0] + p[1]*g[1] + p[2]*g[2]) /
				(p[0]*f3[0] + p[1]*f3[1] + p[2]*f3[2])
		}
		x[4] = kfGammaq(0.5, math.Log(r))
	}
	return 0
}

// =============================================================
//  prob1.c port — bcf_p1_* posterior machinery (single-group)
// =============================================================

const (
	probP1Tiny = 1e-20
)

type bcfP1Aux struct {
	n        int
	M        int
	n1       int
	isIndel  bool
	ploidy   []int
	q2p      [256]float64
	pdg      []float64
	phi      []float64
	phiIndel []float64
	z        []float64
	zswap    []float64
	z1       []float64
	z2       []float64
	afs      []float64
	afs1     []float64
	lf       []float64
	t        float64
	t1       float64
	PL       []int
	PLLen    int
}

type bcfP1RST struct {
	rank0      int
	ac         int
	fExp       float64
	fFlat      float64
	pRefFolded float64
	pRef       float64
	pVarFolded float64
	pVar       float64
}

func initPriorFull(theta float64, M int, phi []float64) {
	sum := 0.0
	for i := 0; i < M; i++ {
		phi[i] = theta / float64(M-i)
		sum += phi[i]
	}
	phi[M] = 1 - sum
}

func bcfP1Init(n int, ploidy []int) *bcfP1Aux {
	ma := &bcfP1Aux{n: n, n1: -1}
	ma.M = 2 * n
	if ploidy != nil {
		all2 := true
		mSum := 0
		for _, p := range ploidy {
			mSum += p
			if p != 2 {
				all2 = false
			}
		}
		if all2 {
			ma.ploidy = nil
		} else {
			ma.ploidy = append([]int(nil), ploidy...)
			ma.M = mSum
		}
	}
	for i := 0; i < 256; i++ {
		ma.q2p[i] = math.Pow(10, -float64(i)/10)
	}
	ma.pdg = make([]float64, 3*ma.n)
	ma.phi = make([]float64, ma.M+1)
	ma.phiIndel = make([]float64, ma.M+1)
	ma.z = make([]float64, ma.M+1)
	ma.zswap = make([]float64, ma.M+1)
	ma.z1 = make([]float64, ma.M+1)
	ma.z2 = make([]float64, ma.M+1)
	ma.afs = make([]float64, ma.M+1)
	ma.afs1 = make([]float64, ma.M+1)
	ma.lf = make([]float64, ma.M+1)
	for i := 0; i <= ma.M; i++ {
		v, _ := math.Lgamma(float64(i + 1))
		ma.lf[i] = v
	}
	initPriorFull(1e-3, ma.M, ma.phi)
	for i := 0; i < ma.M; i++ {
		ma.phiIndel[i] = ma.phi[i] * 0.15
	}
	ma.phiIndel[ma.M] = 1 - ma.phi[ma.M]*0.15
	return ma
}

func clampPL(v int) int {
	if v < 0 {
		return 0
	}
	if v > 255 {
		return 255
	}
	return v
}

func calPdg(nAllele int, PL []int, plLen int, ma *bcfP1Aux) int {
	p := make([]int64, nAllele)
	for j := 0; j < ma.n; j++ {
		pi := PL[j*plLen:]
		pdg := ma.pdg[j*3:]
		pdg[0] = ma.q2p[clampPL(pi[2])]
		pdg[1] = ma.q2p[clampPL(pi[1])]
		pdg[2] = ma.q2p[clampPL(pi[0])]
		for i := 0; i < nAllele; i++ {
			idx := (i+1)*(i+2)/2 - 1
			if idx < plLen {
				p[i] += int64(pi[idx])
			}
		}
	}
	for i := 0; i < nAllele; i++ {
		p[i] = p[i]<<4 | int64(i)
	}
	for i := 1; i < nAllele; i++ {
		for j := i; j > 0 && p[j] < p[j-1]; j-- {
			p[j], p[j-1] = p[j-1], p[j]
		}
	}
	i := nAllele - 1
	for ; i >= 0; i-- {
		if (p[i] & 0xf) == 0 {
			break
		}
	}
	return i
}

func mcCalYCore(ma *bcfP1Aux, beg int) {
	z0 := ma.z
	z1 := ma.zswap
	for i := range z0 {
		z0[i] = 0
	}
	for i := range z1 {
		z1[i] = 0
	}
	z0[0] = 1
	lastMin := 0
	lastMax := 0
	ma.t = 0
	if ma.M == ma.n*2 {
		M := 0
		for _j := beg; _j < ma.n; _j++ {
			j := _j - beg
			minK := lastMin
			maxK := lastMax
			M0 := M
			M += 2
			pdg := ma.pdg[_j*3:]
			p := [3]float64{pdg[0], 2 * pdg[1], pdg[2]}
			for minK < maxK && z0[minK] < probP1Tiny {
				z0[minK] = 0
				z1[minK] = 0
				minK++
			}
			for maxK > minK && z0[maxK] < probP1Tiny {
				z0[maxK] = 0
				z1[maxK] = 0
				maxK--
			}
			maxK += 2
			if minK == 0 {
				k := 0
				z1[k] = float64(M0-k+1) * float64(M0-k+2) * p[0] * z0[k]
			}
			if minK <= 1 {
				k := 1
				z1[k] = float64(M0-k+1)*float64(M0-k+2)*p[0]*z0[k] + float64(k)*float64(M0-k+2)*p[1]*z0[k-1]
			}
			kStart := minK
			if kStart < 2 {
				kStart = 2
			}
			for k := kStart; k <= maxK; k++ {
				z1[k] = float64(M0-k+1)*float64(M0-k+2)*p[0]*z0[k] +
					float64(k)*float64(M0-k+2)*p[1]*z0[k-1] +
					float64(k)*float64(k-1)*p[2]*z0[k-2]
			}
			sum := 0.0
			for k := minK; k <= maxK; k++ {
				sum += z1[k]
			}
			ma.t += math.Log(sum / (float64(M) * float64(M-1)))
			for k := minK; k <= maxK; k++ {
				z1[k] /= sum
			}
			if minK >= 1 {
				z1[minK-1] = 0
			}
			if minK >= 2 {
				z1[minK-2] = 0
			}
			if j < ma.n-1 {
				if maxK+1 <= ma.M {
					z1[maxK+1] = 0
				}
				if maxK+2 <= ma.M {
					z1[maxK+2] = 0
				}
			}
			z0, z1 = z1, z0
			lastMin = minK
			lastMax = maxK
		}
	} else {
		M := 0
		for j := 0; j < ma.n; j++ {
			minK := lastMin
			maxK := lastMax
			pdg := ma.pdg[j*3:]
			for minK < maxK && z0[minK] < probP1Tiny {
				z0[minK] = 0
				z1[minK] = 0
				minK++
			}
			for maxK > minK && z0[maxK] < probP1Tiny {
				z0[maxK] = 0
				z1[maxK] = 0
				maxK--
			}
			M0 := M
			pl := ma.ploidy[j]
			M += pl
			if pl == 1 {
				p := [2]float64{pdg[0], pdg[2]}
				maxK++
				if minK == 0 {
					k := 0
					z1[k] = float64(M0+1-k) * p[0] * z0[k]
				}
				kStart := minK
				if kStart < 1 {
					kStart = 1
				}
				for k := kStart; k <= maxK; k++ {
					z1[k] = float64(M0+1-k)*p[0]*z0[k] + float64(k)*p[1]*z0[k-1]
				}
				sum := 0.0
				for k := minK; k <= maxK; k++ {
					sum += z1[k]
				}
				ma.t += math.Log(sum / float64(M))
				for k := minK; k <= maxK; k++ {
					z1[k] /= sum
				}
				if minK >= 1 {
					z1[minK-1] = 0
				}
				if j < ma.n-1 && maxK+1 <= ma.M {
					z1[maxK+1] = 0
				}
			} else if pl == 2 {
				p := [3]float64{pdg[0], 2 * pdg[1], pdg[2]}
				maxK += 2
				if minK == 0 {
					k := 0
					z1[k] = float64(M0-k+1) * float64(M0-k+2) * p[0] * z0[k]
				}
				if minK <= 1 {
					k := 1
					z1[k] = float64(M0-k+1)*float64(M0-k+2)*p[0]*z0[k] + float64(k)*float64(M0-k+2)*p[1]*z0[k-1]
				}
				kStart := minK
				if kStart < 2 {
					kStart = 2
				}
				for k := kStart; k <= maxK; k++ {
					z1[k] = float64(M0-k+1)*float64(M0-k+2)*p[0]*z0[k] +
						float64(k)*float64(M0-k+2)*p[1]*z0[k-1] +
						float64(k)*float64(k-1)*p[2]*z0[k-2]
				}
				sum := 0.0
				for k := minK; k <= maxK; k++ {
					sum += z1[k]
				}
				ma.t += math.Log(sum / (float64(M) * float64(M-1)))
				for k := minK; k <= maxK; k++ {
					z1[k] /= sum
				}
				if minK >= 1 {
					z1[minK-1] = 0
				}
				if minK >= 2 {
					z1[minK-2] = 0
				}
				if j < ma.n-1 {
					if maxK+1 <= ma.M {
						z1[maxK+1] = 0
					}
					if maxK+2 <= ma.M {
						z1[maxK+2] = 0
					}
				}
			}
			z0, z1 = z1, z0
			lastMin = minK
			lastMax = maxK
		}
	}
	if len(z0) > 0 && len(ma.z) > 0 && &z0[0] != &ma.z[0] {
		copy(ma.z, z0)
	}
}

func mcCalAfs(ma *bcfP1Aux) (fExp, pRefFolded, pVarFolded float64) {
	phi := ma.phi
	if ma.isIndel {
		phi = ma.phiIndel
	}
	for i := range ma.afs1 {
		ma.afs1[i] = 0
	}
	mcCalYCore(ma, 0)
	sum := 0.0
	for k := 0; k <= ma.M; k++ {
		sum += phi[k] * ma.z[k]
	}
	for k := 0; k <= ma.M; k++ {
		ma.afs1[k] = phi[k] * ma.z[k] / sum
		if math.IsNaN(ma.afs1[k]) || math.IsInf(ma.afs1[k], 0) {
			return -1, 0, 0
		}
	}
	sumF := 0.0
	for k := 0; k <= ma.M; k++ {
		sumF += (phi[k] + phi[ma.M-k]) / 2 * ma.z[k]
	}
	sum2 := 0.0
	var k int
	for k = 1; k < ma.M; k++ {
		sum2 += (phi[k] + phi[ma.M-k]) / 2 * ma.z[k]
	}
	pVarFolded = sum2 / sumF
	pRefFolded = (phi[k] + phi[ma.M-k]) / 2 * (ma.z[ma.M] + ma.z[0]) / sumF
	sumK := 0.0
	for k := 0; k <= ma.M; k++ {
		ma.afs[k] += ma.afs1[k]
		sumK += float64(k) * ma.afs1[k]
	}
	fExp = sumK / float64(ma.M)
	return
}

func bcfP1Cal(ma *bcfP1Aux, nAllele int, isSNP bool, rst *bcfP1RST) int {
	if nAllele < 2 {
		return -1
	}
	ma.isIndel = !isSNP
	rst.rank0 = calPdg(nAllele, ma.PL, ma.PLLen, ma)
	rst.fExp, rst.pRefFolded, rst.pVarFolded = mcCalAfs(ma)
	rst.pRef = ma.afs1[ma.M]
	sumV := 0.0
	for k := 0; k < ma.M; k++ {
		sumV += ma.afs1[k]
	}
	rst.pVar = sumV
	maxV := -1.0
	rst.ac = -1
	for k := 0; k <= ma.M; k++ {
		if maxV < ma.z[k] {
			maxV = ma.z[k]
			rst.ac = k
		}
	}
	rst.ac = ma.M - rst.ac
	sumZ := 0.0
	for k := 0; k <= ma.M; k++ {
		sumZ += ma.z[k]
	}
	rst.fFlat = 0
	for k := 0; k <= ma.M; k++ {
		pk := ma.z[k] / sumZ
		rst.fFlat += float64(k) * pk
	}
	rst.fFlat /= float64(ma.M)
	return 0
}

// bcfP1CallGT ports upstream `bcf_p1_call_gt`. Returns (q<<2)|gt
// with gt 0=AA, 1=AB, 2=RR and q the Phred GQ capped at 99.
func bcfP1CallGT(ma *bcfP1Aux, f0 float64, k int, isVar bool) int {
	pdg := ma.pdg[k*3:]
	ploidy := 2
	if ma.ploidy != nil {
		ploidy = ma.ploidy[k]
	}
	var f3 [3]float64
	if ploidy == 2 {
		f3[0] = (1 - f0) * (1 - f0)
		f3[1] = 2 * f0 * (1 - f0)
		f3[2] = f0 * f0
	} else {
		f3[0] = 1 - f0
		f3[1] = 0
		f3[2] = f0
	}
	g := [3]float64{}
	sum := 0.0
	for i := 0; i < 3; i++ {
		g[i] = pdg[i] * f3[i]
		sum += g[i]
	}
	maxV := -1.0
	maxI := 0
	for i := 0; i < 3; i++ {
		g[i] /= sum
		if g[i] > maxV {
			maxV = g[i]
			maxI = i
		}
	}
	if !isVar {
		maxI = 2
		maxV = g[2]
	}
	rem := 1 - maxV
	if rem < 1e-308 {
		rem = 1e-308
	}
	q := int(-4.343*math.Log(rem) + 0.499)
	if q > 99 {
		q = 99
	}
	return q<<2 | maxI
}

// =============================================================
//  ccall.c driver
// =============================================================

// ccallSite is the per-record entry point invoked from callVariant.
// It returns (variant, keep, ok); ok=false signals the caller to
// fall through to the v1 heuristic (no I16 → synthetic fixture).
func ccallSite(v *vcf.Variant, opts CallOptions, samplePloidy []int) (*vcf.Variant, bool, bool) {
	if !hasI16(v) {
		return nil, false, false
	}
	nsmpl := len(v.Samples)
	if nsmpl == 0 {
		return nil, false, false
	}
	nalsOri := 1 + len(v.Alt)
	unseen := 0
	for i, a := range v.Alt {
		if a == "<*>" || a == "<X>" {
			unseen = i + 1
			break
		}
	}

	ngtsOri := nalsOri * (nalsOri + 1) / 2
	plBuf, plLen, plPloidy, ok := ccallDecodePLs(v, nalsOri, samplePloidy)
	if !ok {
		return nil, false, false
	}

	ploidy := make([]int, nsmpl)
	allDiploid := true
	for i := 0; i < nsmpl; i++ {
		p := 2
		if i < len(samplePloidy) {
			p = samplePloidy[i]
		}
		if p == 0 {
			p = 2
		}
		ploidy[i] = p
		if p != 2 {
			allDiploid = false
		}
	}
	var p1Ploidy []int
	if !allDiploid {
		p1Ploidy = ploidy
	}
	p1 := bcfP1Init(nsmpl, p1Ploidy)
	p1.PL = plBuf
	p1.PLLen = plLen

	setPdg3(plBuf, plLen, p1.pdg, nsmpl, ngtsOri)
	em := make([]float64, 10)
	for i := range em {
		em[i] = -1
	}
	bcfEM1(nalsOri, nsmpl, 0, 0x1ff, p1.pdg, em)

	var pr bcfP1RST
	isSNP := variantIsSNP(v)
	bcfP1Cal(p1, nalsOri, isSNP, &pr)

	pref := opts.PvalThreshold
	if pref == 0 {
		pref = 0.5
	}
	if pr.pRef >= pref && opts.VariantsOnly {
		return nil, false, true
	}

	a16 := parseI16(v)
	_, mq, _, isTested, dp4, pv := test16(a16)

	isVar := pr.pRef < pref
	var r float64
	if isVar {
		r = pr.pRef
	} else {
		r = pr.pVar
	}

	nalsOut := 1
	if isVar || opts.KeepAlts {
		if pr.rank0 < 2 {
			nalsOut = 2
		} else {
			nalsOut = pr.rank0 + 1
		}
	}
	if opts.KeepAlts && unseen == nalsOut-1 {
		nalsOut--
	}
	if nalsOut < 1 {
		nalsOut = 1
	}
	if nalsOut > nalsOri {
		nalsOut = nalsOri
	}

	out := *v
	out.Samples = make([]vcf.Sample, nsmpl)
	for i, s := range v.Samples {
		ns := vcf.Sample{Name: s.Name, Data: copyStringMap(s.Data)}
		var gt string
		if plPloidy[i] == 2 || plPloidy[i] == 0 {
			x := bcfP1CallGT(p1, pr.fExp, i, isVar)
			gtIdx := x & 3
			switch gtIdx {
			case 1:
				gt = "0/1"
			case 0:
				gt = "1/1"
			default:
				gt = "0/0"
			}
		} else {
			x := bcfP1CallGT(p1, pr.fExp, i, isVar)
			gtIdx := x & 3
			if gtIdx == 0 {
				gt = "1"
			} else {
				gt = "0"
			}
		}
		ns.Data["GT"] = gt
		if nalsOut < nalsOri {
			if pls, okPL := s.Data["PL"]; okPL {
				ns.Data["PL"] = ccallTrimPL(pls, nalsOri, nalsOut, plPloidy[i])
			}
			if ad, okAD := s.Data["AD"]; okAD {
				ns.Data["AD"] = ccallTrimAD(ad, nalsOri, nalsOut)
			}
		}
		out.Samples[i] = ns
	}

	if nalsOut <= 1 {
		out.Alt = []string{"."}
	} else {
		newAlt := make([]string, 0, nalsOut-1)
		for i := 1; i < nalsOut; i++ {
			newAlt = append(newAlt, v.Alt[i-1])
		}
		out.Alt = newAlt
	}
	out.Format = append([]string{"GT"}, dropFormatKey(v.Format, "GT")...)

	out.Info = copyStringMap(v.Info)
	out.InfoOrder = append([]string(nil), v.InfoOrder...)

	if em[0] >= 0 {
		setInfo(&out, "AF1", formatFloat32G(1-em[0]))
	}
	if em[4] >= 0 && em[4] <= 0.05 {
		setInfo(&out, "G3", strings.Join([]string{
			formatFloat32G(em[3]),
			formatFloat32G(em[2]),
			formatFloat32G(em[1]),
		}, ","))
		setInfo(&out, "HWE", formatFloat32G(em[4]))
	}
	setInfo(&out, "AC1", strconv.Itoa(pr.ac))
	setInfo(&out, "DP4", strings.Join([]string{
		strconv.Itoa(dp4[0]), strconv.Itoa(dp4[1]),
		strconv.Itoa(dp4[2]), strconv.Itoa(dp4[3]),
	}, ","))
	setInfo(&out, "MQ", strconv.Itoa(mq))
	fq := 0.0
	if pr.pRefFolded < 0.5 {
		fq = -4.343 * math.Log(pr.pRefFolded)
	} else {
		fq = 4.343 * math.Log(pr.pVarFolded)
	}
	if fq < -999 {
		fq = -999
	}
	if fq > 999 {
		fq = 999
	}
	setInfo(&out, "FQ", formatFloat32G(fq))
	if isTested {
		setInfo(&out, "PV4", strings.Join([]string{
			formatFloat32G(pv[0]),
			formatFloat32G(pv[1]),
			formatFloat32G(pv[2]),
			formatFloat32G(pv[3]),
		}, ","))
	}
	delInfo(&out, "I16")
	delInfo(&out, "QS")

	if r < 1e-100 {
		out.Qual = 999
	} else {
		q := -4.343 * math.Log(r)
		if q > 999 {
			q = 999
		}
		out.Qual = quantizeQual(q)
	}
	return &out, true, true
}

func ccallDecodePLs(v *vcf.Variant, nals int, samplePloidy []int) ([]int, int, []int, bool) {
	nsmpl := len(v.Samples)
	ngts := nals * (nals + 1) / 2
	buf := make([]int, nsmpl*ngts)
	plPloidy := make([]int, nsmpl)
	for i, s := range v.Samples {
		pl := 2
		if i < len(samplePloidy) {
			pl = samplePloidy[i]
		}
		plPloidy[i] = pl
		raw, okPL := s.Data["PL"]
		if !okPL || raw == "" || raw == "." {
			for k := 0; k < ngts; k++ {
				buf[i*ngts+k] = 255
			}
			continue
		}
		parts := strings.Split(raw, ",")
		if len(parts) == nals && pl == 1 {
			for k := 0; k < ngts; k++ {
				buf[i*ngts+k] = 255
			}
			for a := 0; a < nals; a++ {
				idx := (a+1)*(a+2)/2 - 1
				if a >= len(parts) {
					break
				}
				buf[i*ngts+idx] = ccallParsePL(parts[a])
			}
			continue
		}
		for k := 0; k < ngts; k++ {
			if k < len(parts) {
				buf[i*ngts+k] = ccallParsePL(parts[k])
			} else {
				buf[i*ngts+k] = 255
			}
		}
	}
	return buf, ngts, plPloidy, true
}

func ccallParsePL(s string) int {
	if s == "" || s == "." {
		return 255
	}
	n, err := strconv.Atoi(s)
	if err != nil {
		return 255
	}
	if n < 0 {
		return 0
	}
	if n > 255 {
		return 255
	}
	return n
}

func setPdg3(PLs []int, nGtPerSmpl int, pdg []float64, nsmpl, ngts int) {
	for i := 0; i < nsmpl; i++ {
		pls := PLs[i*nGtPerSmpl:]
		out := pdg[i*3:]
		out[2] = pl2pTable[clampPL(pls[0])]
		out[1] = pl2pTable[clampPL(pls[1])]
		out[0] = pl2pTable[clampPL(pls[2])]
	}
	_ = ngts
}

func hasI16(v *vcf.Variant) bool {
	_, ok := v.Info["I16"]
	return ok
}

func parseI16(v *vcf.Variant) [16]float64 {
	var out [16]float64
	raw := v.Info["I16"]
	if raw == "" {
		return out
	}
	parts := strings.Split(raw, ",")
	for i := 0; i < 16 && i < len(parts); i++ {
		f, err := strconv.ParseFloat(parts[i], 64)
		if err == nil {
			out[i] = f
		}
	}
	return out
}

func test16(anno [16]float64) (hasI16 bool, mq, depth int, isTested bool, dp4 [4]int, p [4]float64) {
	p[0], p[1], p[2], p[3] = 1, 1, 1, 1
	for i := 0; i < 4; i++ {
		dp4[i] = int(anno[i])
	}
	depthF := anno[0] + anno[1] + anno[2] + anno[3]
	depth = int(depthF)
	isTested = (anno[0]+anno[1] > 0 && anno[2]+anno[3] > 0)
	if depthF == 0 {
		hasI16 = false
		return
	}
	hasI16 = true
	mq = int(math.Sqrt((anno[9]+anno[11])/depthF) + 0.499)
	_, _, two := mpileupFisherExact(int64(anno[0]), int64(anno[1]), int64(anno[2]), int64(anno[3]))
	p[0] = two
	n1 := int(anno[0] + anno[1])
	n2 := int(anno[2] + anno[3])
	for i := 1; i < 4; i++ {
		a := [4]float64{anno[4*i], anno[4*i+1], anno[4*i+2], anno[4*i+3]}
		p[i] = ccallTTest(n1, n2, a)
	}
	return
}

func ccallTTest(n1, n2 int, a [4]float64) float64 {
	if n1 == 0 || n2 == 0 || n1+n2 < 3 {
		return 1.0
	}
	u1 := a[0] / float64(n1)
	u2 := a[2] / float64(n2)
	if u1 <= u2 {
		return 1.0
	}
	denom := ((a[1] - float64(n1)*u1*u1) + (a[3] - float64(n2)*u2*u2)) / float64(n1+n2-2) * (1.0/float64(n1) + 1.0/float64(n2))
	if denom <= 0 {
		return 1.0
	}
	t := (u1 - u2) / math.Sqrt(denom)
	v := float64(n1 + n2 - 2)
	if t < 0 {
		return 1
	}
	return 0.5 * kfBetai(0.5*v, 0.5, v/(v+t*t))
}

func ccallTrimPL(s string, nalsOri, nalsOut, ploidy int) string {
	parts := strings.Split(s, ",")
	if ploidy == 1 && len(parts) == nalsOri {
		out := make([]string, nalsOut)
		for i := 0; i < nalsOut; i++ {
			if i < len(parts) {
				out[i] = parts[i]
			} else {
				out[i] = "."
			}
		}
		return strings.Join(out, ",")
	}
	ngtsOut := nalsOut * (nalsOut + 1) / 2
	out := make([]string, ngtsOut)
	for i := 0; i < ngtsOut; i++ {
		if i < len(parts) {
			out[i] = parts[i]
		} else {
			out[i] = "."
		}
	}
	return strings.Join(out, ",")
}

func ccallTrimAD(s string, nalsOri, nalsOut int) string {
	parts := strings.Split(s, ",")
	out := make([]string, nalsOut)
	for i := 0; i < nalsOut; i++ {
		if i < len(parts) {
			out[i] = parts[i]
		} else {
			out[i] = "0"
		}
	}
	return strings.Join(out, ",")
}

func dropFormatKey(orig []string, key string) []string {
	out := make([]string, 0, len(orig))
	for _, k := range orig {
		if k != key {
			out = append(out, k)
		}
	}
	return out
}

func variantIsSNP(v *vcf.Variant) bool {
	if len(v.Ref) != 1 {
		return false
	}
	for _, a := range v.Alt {
		if a == "<*>" || a == "<X>" {
			continue
		}
		if len(a) != 1 {
			return false
		}
	}
	return true
}

func headerHasInfo(hdr *vcf.Header, tag string) bool {
	if hdr == nil {
		return false
	}
	prefix := `##INFO=<ID=` + tag + `,`
	for _, m := range hdr.MetaInfo {
		if strings.HasPrefix(m, prefix) {
			return true
		}
	}
	return false
}

func augmentCallHeaderConsensus(hdr *vcf.Header) *vcf.Header {
	if hdr == nil {
		return hdr
	}
	out := &vcf.Header{Samples: append([]string(nil), hdr.Samples...)}
	for _, m := range hdr.MetaInfo {
		if strings.HasPrefix(m, `##INFO=<ID=I16,`) || strings.HasPrefix(m, `##INFO=<ID=QS,`) {
			continue
		}
		out.MetaInfo = append(out.MetaInfo, m)
	}
	appends := []struct {
		marker string
		line   string
	}{
		{`##FORMAT=<ID=GT,`, `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`},
		{`##INFO=<ID=AF1,`, `##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">`},
		{`##INFO=<ID=AF2,`, `##INFO=<ID=AF2,Number=1,Type=Float,Description="Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)">`},
		{`##INFO=<ID=AC1,`, `##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">`},
		{`##INFO=<ID=MQ,`, `##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">`},
		{`##INFO=<ID=FQ,`, `##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">`},
		{`##INFO=<ID=PV4,`, `##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">`},
		{`##INFO=<ID=G3,`, `##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">`},
		{`##INFO=<ID=HWE,`, `##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">`},
		{`##INFO=<ID=DP4,`, `##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">`},
	}
	for _, d := range appends {
		found := false
		for _, m := range out.MetaInfo {
			if strings.HasPrefix(m, d.marker) {
				found = true
				break
			}
		}
		if !found {
			out.MetaInfo = append(out.MetaInfo, d.line)
		}
	}
	return out
}
