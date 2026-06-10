// Per-sample FORMAT tag helpers for bcftools mpileup. The bitset in
// MpileupOptions.FmtFlag selects optional per-sample columns beyond the
// always-on FORMAT/PL: DP (high-quality coverage), DV (high-quality
// non-reference coverage), DP4 (the four strand-split counts), and SP
// (Fisher's exact strand-bias phred score). The values are sourced from
// the per-sample I16 strand-split tallies that bcfCallCombine already
// stores on each bcfCallret (`anno[0..3]`), so this file only renders
// them — no recomputation. Mirrors bam2bcf.c:1341-1389.

package bcftools

import (
	"math"
	"strconv"
	"strings"
)

// perSampleDP4 returns the four strand-split depths for one sample at
// the current site: {ref-fwd, ref-rev, alt-fwd, alt-rev}, matching
// upstream's bc->DP4[4*i .. 4*i+3] (bam2bcf.c:1052-1058) which is sourced
// from calls[i].anno[0..3].
func perSampleDP4(r *bcfCallret) [4]int {
	return [4]int{
		int(r.anno[0]),
		int(r.anno[1]),
		int(r.anno[2]),
		int(r.anno[3]),
	}
}

// formatPerSampleDP renders FORMAT/DP for one sample: the sum of the
// four strand-split depths (bam2bcf.c:1341-1346).
func formatPerSampleDP(dp4 [4]int) string {
	return strconv.Itoa(dp4[0] + dp4[1] + dp4[2] + dp4[3])
}

// formatPerSampleDV renders FORMAT/DV for one sample: the count of
// high-quality non-reference bases (alt-fwd + alt-rev, bam2bcf.c:
// 1348-1353).
func formatPerSampleDV(dp4 [4]int) string {
	return strconv.Itoa(dp4[2] + dp4[3])
}

// formatPerSampleDP4 renders FORMAT/DP4 for one sample: the four
// strand-split counts as `fwdRef,revRef,fwdAlt,revAlt`
// (bam2bcf.c:1374-1375).
func formatPerSampleDP4(dp4 [4]int) string {
	var b strings.Builder
	for k, v := range dp4 {
		if k > 0 {
			b.WriteByte(',')
		}
		b.WriteString(strconv.Itoa(v))
	}
	return b.String()
}

// formatPerSampleSP renders FORMAT/SP for one sample: a phred-scaled
// Fisher's exact strand-bias score on the 2x2 contingency table of
// (ref-fwd, ref-rev, alt-fwd, alt-rev). Mirrors bam2bcf.c:1355-1372 —
// returns "0" when any margin is below 2, otherwise
// int(-4.343*ln(p_two)+.499), capped at 255.
func formatPerSampleSP(dp4 [4]int) string {
	fwdRef, revRef, fwdAlt, revAlt := dp4[0], dp4[1], dp4[2], dp4[3]
	if fwdRef+revRef < 2 || fwdAlt+revAlt < 2 ||
		fwdRef+fwdAlt < 2 || revRef+revAlt < 2 {
		return "0"
	}
	_, _, two := mpileupFisherExact(int64(fwdRef), int64(revRef),
		int64(fwdAlt), int64(revAlt))
	x := int(-4.343*math.Log(two) + 0.499)
	if x > 255 {
		x = 255
	}
	if x < 0 {
		x = 0
	}
	return strconv.Itoa(x)
}

// formatPerSampleDPR renders FORMAT/DPR (deprecated synonym for AD):
// the elementwise sum of ADF and ADR per allele. bam2bcf.c:1380-1386
// emits DPR from the same buffer as AD when either bit is set.
func formatPerSampleDPR(adf, adr []int) string {
	return formatPerAlleleSum(adf, adr)
}

// sumPerAllele adds per-sample per-allele counts into a single
// nAlleles-long slice. It mirrors the cross-sample sum that bam2bcf.c
// performs implicitly by accumulating into bca->ADF/ADR slots indexed
// by site allele (the first B2B_MAX_ALLELES slots of bc->ADF/ADR are
// the INFO-row totals).
func sumPerAllele(perSample [][]int, nAlleles int) []int {
	out := make([]int, nAlleles)
	for _, row := range perSample {
		for i := 0; i < nAlleles && i < len(row); i++ {
			out[i] += row[i]
		}
	}
	return out
}

// perAlleleToCSV joins a per-allele count slice with commas.
func perAlleleToCSV(vals []int) string {
	if len(vals) == 0 {
		return "0"
	}
	var b strings.Builder
	for i, v := range vals {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(strconv.Itoa(v))
	}
	return b.String()
}

// --- Fisher's exact test (pure-Go port of htslib's kt_fisher_exact) ---
//
// The same routine lives in tools/bedfisher; reproduce it here so the
// bcftools package has no cross-tool dependency. The implementation is
// byte-equivalent to upstream's htslib kfunc.c.

type mpileupHGAcc struct {
	n11, n1Row, n1Col, n int64
	p                    float64
}

func mpileupLgamma(z float64) float64 {
	v, _ := math.Lgamma(z)
	return v
}

func mpileupLbinom(n, k int64) float64 {
	if k == 0 || n == k {
		return 0
	}
	return mpileupLgamma(float64(n+1)) - mpileupLgamma(float64(k+1)) -
		mpileupLgamma(float64(n-k+1))
}

func mpileupHypergeo(n11, n1Row, n1Col, n int64) float64 {
	return math.Exp(mpileupLbinom(n1Row, n11) +
		mpileupLbinom(n-n1Row, n1Col-n11) - mpileupLbinom(n, n1Col))
}

func mpileupHypergeoAcc(n11, n1Row, n1Col, n int64, aux *mpileupHGAcc) float64 {
	if n1Row != 0 || n1Col != 0 || n != 0 {
		aux.n11 = n11
		aux.n1Row = n1Row
		aux.n1Col = n1Col
		aux.n = n
	} else {
		if n11%11 != 0 && n11+aux.n-aux.n1Row-aux.n1Col != 0 {
			if n11 == aux.n11+1 {
				aux.p *= (float64(aux.n1Row-aux.n11) / float64(n11)) *
					(float64(aux.n1Col-aux.n11) / float64(n11+aux.n-aux.n1Row-aux.n1Col))
				aux.n11 = n11
				return aux.p
			}
			if n11 == aux.n11-1 {
				aux.p *= (float64(aux.n11) / float64(aux.n1Row-n11)) *
					(float64(aux.n11+aux.n-aux.n1Row-aux.n1Col) / float64(aux.n1Col-n11))
				aux.n11 = n11
				return aux.p
			}
		}
		aux.n11 = n11
	}
	aux.p = mpileupHypergeo(aux.n11, aux.n1Row, aux.n1Col, aux.n)
	return aux.p
}

// mpileupFisherExact returns the (left, right, two-tail) p-values for
// the 2x2 contingency table {n11, n12, n21, n22}. Pure-Go port of
// htslib's kt_fisher_exact (kfunc.c).
func mpileupFisherExact(n11, n12, n21, n22 int64) (float64, float64, float64) {
	n1Row := n11 + n12
	n1Col := n11 + n21
	n := n11 + n12 + n21 + n22

	max := n1Col
	if n1Row < max {
		max = n1Row
	}
	min := n1Row + n1Col - n
	if min < 0 {
		min = 0
	}

	if min == max {
		return 1.0, 1.0, 1.0
	}

	var aux mpileupHGAcc
	q := mpileupHypergeoAcc(n11, n1Row, n1Col, n, &aux)
	if q == 0.0 {
		if int64(n11)*(n+2) < (n1Col+1)*(n1Row+1) {
			return 0.0, 1.0, 0.0
		}
		return 1.0, 0.0, 0.0
	}

	// Left tail.
	p := mpileupHypergeoAcc(min, 0, 0, 0, &aux)
	var left, right float64
	i := min + 1
	for p < 0.99999999*q && i <= max {
		left += p
		p = mpileupHypergeoAcc(i, 0, 0, 0, &aux)
		i++
	}
	i--
	if p < 1.00000001*q {
		left += p
	} else {
		i--
	}

	// Right tail.
	p = mpileupHypergeoAcc(max, 0, 0, 0, &aux)
	j := max - 1
	for p < 0.99999999*q && j >= 0 {
		right += p
		p = mpileupHypergeoAcc(j, 0, 0, 0, &aux)
		j--
	}
	j++
	if p < 1.00000001*q {
		right += p
	} else {
		j++
	}

	two := left + right
	if two > 1 {
		two = 1
	}

	if mpileupAbsI(i-n11) < mpileupAbsI(j-n11) {
		right = 1.0 - left + q
	} else {
		left = 1.0 - right + q
	}
	return left, right, two
}

func mpileupAbsI(x int64) int64 {
	if x < 0 {
		return -x
	}
	return x
}
