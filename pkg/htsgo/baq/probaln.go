// Package baq implements BAQ (Base Alignment Quality) realignment, a faithful
// pure-Go port of htslib's realn.c and probaln.c.
//
// BAQ measures the probability that a read base is misaligned. It is computed
// with a probabilistic banded "glocal" forward-backward HMM that aligns the
// read against the reference window around its mapped position. The result is
// written back as the BQ:Z: aux tag (or, when applied, the base qualities are
// lowered and a ZQ:Z: tag records the adjustment so it can be reversed).
//
// The two public entry points mirror the htslib functions of the same intent:
//
//   - ProbalnGlocal — the core HMM (htslib probaln_glocal).
//   - SamProbRealn  — the CIGAR-aware driver that extracts read/reference
//     windows, runs the HMM and edits the record (htslib sam_prob_realn).
//   - SamCapMapq    — caps mapping quality from observed mismatches
//     (htslib sam_cap_mapq).
//
// The numerical core is ported verbatim from probaln.c (the non-PROBALN_ORIG
// reference path), including the set_u banding macro, the EI/EM constants, the
// g_qual2prob table and the scaled forward/backward recurrences, so results are
// byte-identical to upstream for the same inputs.
package baq

import (
	"errors"
	"math"
)

// EI and EM are the emission constants from probaln.c: EI is the per-base
// insertion emission probability and EM the mismatch emission scale.
const (
	eI = 0.25
	eM = 0.33333333333
)

// Par holds the probaln_par_t tuning parameters: D is the gap-open
// probability, E the gap-extension probability and BW the band half-width.
type Par struct {
	D  float64 // gap-open probability
	E  float64 // gap-extension probability
	BW int     // band half-width
}

// ErrProbalnFailed indicates probaln_glocal could not complete (it maps to the
// upstream INT_MIN return).
var ErrProbalnFailed = errors.New("baq: probaln_glocal failed")

// qual2prob is the g_qual2prob table: index i holds 10^(-i/10), the error
// probability for a Phred quality of i.
var qual2prob = func() [256]float64 {
	var t [256]float64
	for i := range t {
		t[i] = math.Pow(10, -float64(i)/10.0)
	}
	return t
}()

// setU reproduces the C set_u(u,b,i,k) macro: it maps a (query i, band k)
// coordinate to a flat offset into a banded row of three scores.
func setU(b, i, k int) int {
	x := i - b
	if x < 0 {
		x = 0
	}
	return (k - x + 1) * 3
}

// ProbalnGlocal is a faithful port of htslib's probaln_glocal. It performs a
// probabilistic banded glocal forward-backward alignment of query against ref.
//
// ref and query are residue sequences encoded 0/1/2/3 for A/C/G/T and 4 for an
// ambiguous residue. iqual holds the per-query-base Phred qualities; passing
// nil makes every base quality 30.
//
// When state and q are both non-nil the backward pass and MAP decoding run:
// state[i] receives (refPos<<2 | s) where s is 0 for a match or 1 for an
// insertion, and q[i] the Phred-scaled posterior probability that state[i] is
// wrong (capped at 99). Passing nil for both skips the backward pass.
//
// It returns the Phred-scaled alignment likelihood, or ErrProbalnFailed.
func ProbalnGlocal(ref []byte, query []byte, iqual []byte, c Par, state []int, q []byte) (int, error) {
	lRef := len(ref)
	lQuery := len(query)
	if lRef == 0 || lQuery == 0 {
		return 0, nil
	}

	isBackward := state != nil && q != nil

	// Band width selection (probaln.c lines 95-99).
	bw := lRef
	if lQuery > bw {
		bw = lQuery
	}
	if bw > c.BW {
		bw = c.BW
	}
	if d := abs(lRef - lQuery); bw < d {
		bw = d
	}
	bw2 := bw*2 + 1
	iDim := bw2*3 + 6
	if lRef*3+6 < iDim {
		iDim = lRef*3 + 6
	}

	// Forward / backward matrices and the per-row scaling array.
	f := make([]float64, (lQuery+1)*iDim)
	var b []float64
	if isBackward {
		b = make([]float64, (lQuery+1)*iDim)
	}
	s := make([]float64, lQuery+2)

	// Per-base error probabilities (probaln.c qual[]).
	qual := make([]float64, lQuery)
	for i := 0; i < lQuery; i++ {
		if iqual != nil {
			qual[i] = qual2prob[iqual[i]]
		} else {
			qual[i] = qual2prob[30]
		}
	}

	// Transition probabilities m[9] (probaln.c lines 132-142).
	sM := 1.0 / float64(2*lQuery+2)
	sI := sM
	var m [9]float64
	m[0] = (1 - c.D - c.D) * (1 - sM)
	m[1] = c.D * (1 - sM)
	m[2] = c.D * (1 - sM)
	m[3] = (1 - c.E) * (1 - sI)
	m[4] = c.E * (1 - sI)
	m[5] = 0
	m[6] = 1 - c.E
	m[7] = 0
	m[8] = c.E
	bM := (1 - c.D) / float64(lRef)
	bI := c.D / float64(lRef)

	/*** forward ***/
	// f[0]
	f[0*iDim+setU(bw, 0, 0)] = 1.0
	s[0] = 1.0
	{ // f[1]
		fi := f[1*iDim:]
		end := bw + 1
		if lRef < end {
			end = lRef
		}
		sum := 0.0
		for k := 1; k <= end; k++ {
			var e float64
			switch {
			case ref[k-1] > 3 || query[0] > 3:
				e = 1.0
			case ref[k-1] == query[0]:
				e = 1.0 - qual[0]
			default:
				e = qual[0] * eM
			}
			u := setU(bw, 1, k)
			fi[u+0] = e * bM
			fi[u+1] = eI * bI
			sum += fi[u] + fi[u+1]
		}
		s[1] = sum
	}
	// f[2..lQuery]
	for i := 2; i <= lQuery; i++ {
		fi := f[i*iDim:]
		fi1 := f[(i-1)*iDim:]
		qli := qual[i-1]
		qyi := query[i-1]
		beg := 1
		if x := i - bw; beg < x {
			beg = x
		}
		end := lRef
		if x := i + bw; end > x {
			end = x
		}
		E := [4]float64{qli * eM, 1.0 - qli, 1.0, 1.0}
		M := 1.0 / s[i-1]

		var xm [5]float64
		xm[0] = M * m[0]
		xm[1] = M * m[3]
		xm[2] = M * m[6]
		xm[3] = eI * M * m[1]
		xm[4] = eI * M * m[4]

		u := setU(bw, i, beg)
		v11 := setU(bw, i-1, beg-1)
		lX0 := m[2] * fi[u+0]
		lX2 := m[8] * fi[u+2]
		sum := 0.0
		xi := u
		yi := v11
		for k := beg; k <= end; k++ {
			cond := 0
			if ref[k-1] > 3 || qyi > 3 {
				cond = 2
			} else if ref[k-1] == qyi {
				cond = 1
			}

			z0 := xm[0] * fi1[yi+0]
			z1 := xm[1] * fi1[yi+1]
			z2 := xm[2] * fi1[yi+2]
			z3 := xm[3] * fi1[yi+3]
			z4 := xm[4] * fi1[yi+4]

			fi[xi+0] = E[cond] * (z0 + z1 + z2)
			fi[xi+1] = z3 + z4
			fi[xi+2] = lX0 + lX2
			sum += fi[xi+0] + fi[xi+1] + fi[xi+2]

			lX0 = m[2] * fi[xi+0]
			lX2 = m[8] * fi[xi+2]
			xi += 3
			yi += 3
		}
		s[i] = sum
	}
	{ // f[lQuery+1]
		M := 1.0 / s[lQuery]
		sum := 0.0
		for k := 1; k <= lRef; k++ {
			u := setU(bw, lQuery, k)
			if u < 3 || u >= iDim {
				continue
			}
			sum += M*f[lQuery*iDim+u+0]*sM + M*f[lQuery*iDim+u+1]*sI
		}
		s[lQuery+1] = sum
	}

	// likelihood
	var pr int
	{
		p := 1.0
		pr1 := 0.0
		for i := 0; i <= lQuery+1; i++ {
			p *= s[i]
			if p < 1e-100 {
				pr1 += -4.343 * math.Log(p)
				p = 1.0
			}
		}
		pr1 += -4.343 * math.Log(p*float64(lRef)*float64(lQuery))
		pr = int(pr1 + 0.499)
		if !isBackward {
			return pr, nil
		}
	}

	/*** backward ***/
	// b[lQuery]
	for k := 1; k <= lRef; k++ {
		bi := b[lQuery*iDim:]
		u := setU(bw, lQuery, k)
		if u < 3 || u >= iDim {
			continue
		}
		bi[u+0] = sM / s[lQuery] / s[lQuery+1]
		bi[u+1] = sI / s[lQuery] / s[lQuery+1]
	}
	// b[lQuery-1..1]
	for i := lQuery - 1; i >= 1; i-- {
		bi := b[i*iDim:]
		bi1 := b[(i+1)*iDim:]
		var y float64
		if i > 1 {
			y = 1.0
		}
		qli1 := qual[i]
		qyi1 := query[i]
		beg := 1
		if x := i - bw; beg < x {
			beg = x
		}
		end := lRef
		if x := i + bw; end > x {
			end = x
		}
		E := [4]float64{qli1 * eM, 1.0 - qli1, 1.0, 1.0}

		u := setU(bw, i, end)
		v10 := setU(bw, i+1, end)
		xi := u
		yi := v10
		xi5 := bi[xi+5]
		e1 := eI * m[1]
		e4 := eI * m[4]
		n := 1.0 / s[i]
		for k := end; k >= beg; k-- {
			var e float64
			if k < lRef {
				cond := 0
				if ref[k] > 3 || qyi1 > 3 {
					cond = 2
				} else if ref[k] == qyi1 {
					cond = 1
				}
				e = E[cond] * bi1[yi+3]
			}

			bi[xi+1] = e*m[3] + e4*bi1[yi+1]
			bi[xi+0] = e*m[0] + e1*bi1[yi+1] + m[2]*xi5
			bi[xi+2] = (e*m[6] + m[8]*xi5) * y
			xi5 = bi[xi+2]

			bi[xi+1] *= n
			bi[xi+0] *= n
			bi[xi+2] *= n

			xi -= 3
			yi -= 3
		}
	}
	{ // b[0]
		end := bw + 1
		if lRef < end {
			end = lRef
		}
		sum := 0.0
		for k := end; k >= 1; k-- {
			var e float64
			switch {
			case ref[k-1] > 3 || query[0] > 3:
				e = 1.0
			case ref[k-1] == query[0]:
				e = 1.0 - qual[0]
			default:
				e = qual[0] * eM
			}
			u := setU(bw, 1, k)
			if u < 3 || u >= iDim {
				continue
			}
			sum += e*b[1*iDim+u+0]*bM + eI*b[1*iDim+u+1]*bI
		}
		b[0*iDim+setU(bw, 0, 0)] = sum / s[0]
	}

	/*** MAP ***/
	for i := 1; i <= lQuery; i++ {
		fi := f[i*iDim:]
		bi := b[i*iDim:]
		beg := 1
		if x := i - bw; beg < x {
			beg = x
		}
		end := lRef
		if x := i + bw; end > x {
			end = x
		}
		M := 1.0 / s[i]
		sum := 0.0
		max := 0.0
		maxK := -1
		u := setU(bw, i, beg)
		for k := beg; k <= end; k++ {
			z1 := M * fi[u+0] * bi[u+0]
			z2 := M * fi[u+1] * bi[u+1]
			which := 0
			zm := z1
			if z2 > z1 {
				which = 1
				zm = z2
			}
			if zm > max {
				max = zm
				maxK = (k-1)<<2 | which
			}
			sum += z1 + z2
			u += 3
		}
		max /= sum
		if state != nil {
			state[i-1] = maxK
		}
		if q != nil {
			k := int(-4.343*math.Log(1.0-max) + 0.499)
			if k > 100 {
				k = 99
			}
			q[i-1] = byte(k)
		}
	}

	return pr, nil
}

// abs returns the absolute value of an int.
func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}
