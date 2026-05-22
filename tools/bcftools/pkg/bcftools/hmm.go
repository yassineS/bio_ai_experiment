package bcftools

// Generic 2-plus-state hidden Markov model engine, ported from
// reference_code/bcftools/HMM.c. bcftools `roh` builds on this engine:
// it precomputes powers of the transition matrix so that the
// per-marker transition probabilities scale with the physical (or
// genetic) distance between consecutive sites.
//
// The port keeps the upstream maths verbatim — Viterbi decoding,
// forward-backward posteriors and Baum-Welch re-estimation are all
// linear-space (not log-space) with per-site renormalisation, exactly
// as in HMM.c — so that the numeric outputs are byte-faithful to
// upstream.

import "math"

// matAt indexes a row-major n*n matrix the way upstream's MAT macro
// does: MAT(ptr,n,i,j) == ptr[i*n+j].
func matAt(m []float64, n, i, j int) float64 { return m[i*n+j] }

// matSet writes element (i,j) of a row-major n*n matrix.
func matSet(m []float64, n, i, j int, v float64) { m[i*n+j] = v }

// multiplyMatrix computes dst = a*b for two n*n row-major matrices.
// It mirrors HMM.c's multiply_matrix, tolerating dst aliasing a or b.
func multiplyMatrix(n int, a, b, dst []float64) {
	out := dst
	tmp := make([]float64, n*n)
	if &a[0] == &dst[0] || &b[0] == &dst[0] {
		out = tmp
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			val := 0.0
			for k := 0; k < n; k++ {
				val += matAt(a, n, i, k) * matAt(b, n, k, j)
			}
			matSet(out, n, i, j, val)
		}
	}
	if &out[0] != &dst[0] {
		copy(dst, out)
	}
}

// setTprobFunc is an optional callback that adjusts the transition
// matrix at every site, given the previous and current positions. It
// matches upstream's set_tprob_f. The matrix passed in is the engine's
// scratch curr_tprob and may be mutated in place.
type setTprobFunc func(prevPos, pos uint32, tprob []float64)

// hmm is the HMM engine state. The exported wrapper for bcftools roh
// lives in roh.go; hmm itself is unexported.
type hmm struct {
	nstates int

	// tprobArr holds ntprobArr precomputed transition matrices. The
	// first is the base matrix; entry i is base^(i+1).
	tprobArr  []float64
	ntprobArr int

	currTprob []float64 // scratch matrix, valid for one Viterbi/fwd step

	setTprob setTprobFunc

	// initVit is the initial per-state Viterbi probability vector.
	initVit []float64
	// stateVit / stateFwd carry forward across buffer flushes; snapPos
	// records the position a snapshot was taken at (0 == none).
	stateVit []float64
	stateFwd []float64
	snapPos  uint32

	// Outputs of the last run.
	vpath []uint8   // Viterbi path, length nsites (states only)
	fwd   []float64 // forward-backward posteriors, (nsites)*nstates

	snapAtPos   uint32
	snapVitProb []float64
	snapFwdProb []float64
}

// newHMM creates a 2-or-more-state HMM. tprob is the base transition
// matrix (row-major, nstates*nstates); ntprob requests that many
// precomputed matrix powers (upstream roh uses 10000).
func newHMM(nstates int, tprob []float64, ntprob int) *hmm {
	h := &hmm{
		nstates:   nstates,
		currTprob: make([]float64, nstates*nstates),
	}
	h.setTprobMatrix(tprob, ntprob)
	h.initStates(nil)
	return h
}

// initStates sets the initial state probabilities; nil means uniform.
func (h *hmm) initStates(probs []float64) {
	n := h.nstates
	h.initVit = make([]float64, n)
	h.stateVit = make([]float64, n)
	h.stateFwd = make([]float64, n)
	h.snapPos = 0
	if probs != nil {
		copy(h.initVit, probs)
		sum := 0.0
		for _, p := range h.initVit {
			sum += p
		}
		for i := range h.initVit {
			h.initVit[i] /= sum
		}
	} else {
		for i := range h.initVit {
			h.initVit[i] = 1.0 / float64(n)
		}
	}
	copy(h.stateVit, h.initVit)
	copy(h.stateFwd, h.initVit)
}

// setTprobMatrix installs a new base transition matrix and precomputes
// ntprob matrix powers, mirroring HMM.c's hmm_set_tprob.
func (h *hmm) setTprobMatrix(tprob []float64, ntprob int) {
	ntprobLogic := ntprob
	if ntprob <= 0 {
		ntprob = 1
	}
	nmat := h.nstates * h.nstates
	if h.tprobArr == nil || len(h.tprobArr) < nmat*ntprob {
		h.tprobArr = make([]float64, nmat*ntprob)
	}
	h.ntprobArr = ntprobLogic
	copy(h.tprobArr, tprob[:nmat])
	for i := 1; i < ntprob; i++ {
		multiplyMatrix(h.nstates,
			h.tprobArr[:nmat],
			h.tprobArr[(i-1)*nmat:i*nmat],
			h.tprobArr[i*nmat:(i+1)*nmat])
	}
}

// tprobMatrix returns the base transition matrix slice.
func (h *hmm) tprobMatrix() []float64 {
	nmat := h.nstates * h.nstates
	return h.tprobArr[:nmat]
}

// reset clears the carried-forward state back to the initial vector.
func (h *hmm) reset() {
	h.snapPos = 0
	h.snapAtPos = 0
	copy(h.stateVit, h.initVit)
	copy(h.stateFwd, h.initVit)
}

// requestSnapshot asks the next Viterbi / forward-backward pass to
// capture the state vectors once it reaches position pos, mirroring
// HMM.c's hmm_snapshot.
func (h *hmm) requestSnapshot(pos uint32) {
	h.snapAtPos = pos
	if h.snapVitProb == nil {
		h.snapVitProb = make([]float64, h.nstates)
		h.snapFwdProb = make([]float64, h.nstates)
	}
}

// restoreSnapshot loads the carried-forward state from the last
// captured snapshot, mirroring HMM.c's hmm_restore. With no snapshot
// it falls back to the initial state.
func (h *hmm) restoreSnapshot() {
	if h.snapAtPos == 0 {
		h.snapPos = 0
		copy(h.stateVit, h.initVit)
		copy(h.stateFwd, h.initVit)
		return
	}
	h.snapPos = h.snapAtPos
	copy(h.stateVit, h.snapVitProb)
	copy(h.stateFwd, h.snapFwdProb)
}

// setTprobForDiff loads currTprob with the transition matrix
// appropriate to a position gap of posDiff, replicating HMM.c's
// _set_tprob: pick the (posDiff mod ntprob)-th precomputed power then
// matrix-multiply by full blocks.
func (h *hmm) setTprobForDiff(posDiff int) {
	nmat := h.nstates * h.nstates
	n := 0
	if h.ntprobArr != 0 {
		n = posDiff % h.ntprobArr
	}
	copy(h.currTprob, h.tprobArr[n*nmat:(n+1)*nmat])
	if h.ntprobArr > 0 {
		blocks := posDiff / h.ntprobArr
		last := h.tprobArr[(h.ntprobArr-1)*nmat : h.ntprobArr*nmat]
		for i := 0; i < blocks; i++ {
			multiplyMatrix(h.nstates, last, h.currTprob, h.currTprob)
		}
	}
}

// posGap mirrors HMM.c: a forward step from prev to cur uses a gap of
// cur-prev-1 (0 when the positions coincide).
func posGap(prev, cur uint32) int {
	if cur == prev {
		return 0
	}
	return int(cur - prev - 1)
}

// runViterbi decodes the most likely state path for nsites sites.
// eprobs is laid out [nsites*nstates]; sites are 0-based positions.
// The decoded path is left in h.vpath.
func (h *hmm) runViterbi(nsites int, eprobs []float64, sites []uint32) {
	n := h.nstates
	if cap(h.vpath) < nsites {
		h.vpath = make([]uint8, nsites)
	}
	h.vpath = h.vpath[:nsites]
	back := make([]uint8, nsites*n)

	vprob := make([]float64, n)
	vtmp := make([]float64, n)
	copy(vprob, h.stateVit)
	prevPos := sites[0]
	if h.snapPos != 0 {
		prevPos = h.snapPos
	}

	for i := 0; i < nsites; i++ {
		eprob := eprobs[i*n : (i+1)*n]
		posDiff := 0
		if sites[i] != prevPos {
			posDiff = int(sites[i] - prevPos - 1)
		}
		h.setTprobForDiff(posDiff)
		if h.setTprob != nil {
			h.setTprob(prevPos, sites[i], h.currTprob)
		}
		prevPos = sites[i]

		vnorm := 0.0
		for j := 0; j < n; j++ {
			vmax := 0.0
			kMax := 0
			for k := 0; k < n; k++ {
				pval := vprob[k] * matAt(h.currTprob, n, j, k)
				if vmax < pval {
					vmax = pval
					kMax = k
				}
			}
			back[i*n+j] = uint8(kMax)
			vtmp[j] = vmax * eprob[j]
			vnorm += vtmp[j]
		}
		for j := 0; j < n; j++ {
			vtmp[j] /= vnorm
		}
		vprob, vtmp = vtmp, vprob

		if h.snapAtPos != 0 && sites[i] == h.snapAtPos {
			copy(h.snapVitProb, vprob)
		}
	}

	iptr := 0
	for i := 1; i < n; i++ {
		if vprob[iptr] < vprob[i] {
			iptr = i
		}
	}
	for i := nsites - 1; i >= 0; i-- {
		prev := back[i*n+iptr]
		h.vpath[i] = uint8(iptr)
		iptr = int(prev)
	}
}

// runFwdBwd computes per-site forward-backward posteriors, leaving them
// in h.fwd as an [nsites*nstates] row-major array.
func (h *hmm) runFwdBwd(nsites int, eprobs []float64, sites []uint32) {
	n := h.nstates
	fwd := make([]float64, (nsites+1)*n)
	copy(fwd[:n], h.stateFwd)
	bwd := make([]float64, n)
	bwdTmp := make([]float64, n)
	for i := range bwd {
		bwd[i] = 1
	}

	prevPos := sites[0]
	if h.snapPos != 0 {
		prevPos = h.snapPos
	}
	for i := 0; i < nsites; i++ {
		fwdPrev := fwd[i*n : (i+1)*n]
		fwdCur := fwd[(i+1)*n : (i+2)*n]
		eprob := eprobs[i*n : (i+1)*n]

		posDiff := 0
		if sites[i] != prevPos {
			posDiff = int(sites[i] - prevPos - 1)
		}
		h.setTprobForDiff(posDiff)
		if h.setTprob != nil {
			h.setTprob(prevPos, sites[i], h.currTprob)
		}
		prevPos = sites[i]

		norm := 0.0
		for j := 0; j < n; j++ {
			pval := 0.0
			for k := 0; k < n; k++ {
				pval += fwdPrev[k] * matAt(h.currTprob, n, j, k)
			}
			fwdCur[j] = pval * eprob[j]
			norm += fwdCur[j]
		}
		for j := 0; j < n; j++ {
			fwdCur[j] /= norm
		}
		if h.snapAtPos != 0 && sites[i] == h.snapAtPos {
			copy(h.snapFwdProb, fwdCur)
		}
	}

	prevPos = sites[nsites-1]
	for i := 0; i < nsites; i++ {
		fwdCur := fwd[(nsites-i)*n : (nsites-i+1)*n]
		eprob := eprobs[(nsites-i-1)*n : (nsites-i)*n]

		posDiff := 0
		if sites[nsites-i-1] != prevPos {
			posDiff = int(prevPos - sites[nsites-i-1] - 1)
		}
		h.setTprobForDiff(posDiff)
		if h.setTprob != nil {
			h.setTprob(sites[nsites-i-1], prevPos, h.currTprob)
		}
		prevPos = sites[nsites-i-1]

		bwdNorm := 0.0
		for j := 0; j < n; j++ {
			pval := 0.0
			for k := 0; k < n; k++ {
				pval += bwd[k] * eprob[k] * matAt(h.currTprob, n, k, j)
			}
			bwdTmp[j] = pval
			bwdNorm += pval
		}
		norm := 0.0
		for j := 0; j < n; j++ {
			bwdTmp[j] /= bwdNorm
			fwdCur[j] *= bwd[j]
			norm += fwdCur[j]
		}
		for j := 0; j < n; j++ {
			fwdCur[j] /= norm
		}
		bwd, bwdTmp = bwdTmp, bwd
	}
	// Drop the leading prior row so h.fwd[i] is site i's posterior.
	h.fwd = fwd[n:]
}

// runBaumWelch performs a single Baum-Welch re-estimation pass over
// nsites sites and returns the re-estimated transition matrix. It
// mirrors HMM.c's hmm_run_baum_welch.
func (h *hmm) runBaumWelch(nsites int, eprobs []float64, sites []uint32) []float64 {
	n := h.nstates
	fwd := make([]float64, (nsites+1)*n)
	copy(fwd[:n], h.stateFwd)
	bwd := make([]float64, n)
	bwdTmp := make([]float64, n)
	for i := range bwd {
		bwd[i] = 1
	}

	tmpXi := make([]float64, n*n)
	tmpGamma := make([]float64, n)
	fwdBwd := make([]float64, n)

	prevPos := sites[0]
	if h.snapPos != 0 {
		prevPos = h.snapPos
	}
	for i := 0; i < nsites; i++ {
		fwdPrev := fwd[i*n : (i+1)*n]
		fwdCur := fwd[(i+1)*n : (i+2)*n]
		eprob := eprobs[i*n : (i+1)*n]

		posDiff := 0
		if sites[i] != prevPos {
			posDiff = int(sites[i] - prevPos - 1)
		}
		h.setTprobForDiff(posDiff)
		if h.setTprob != nil {
			h.setTprob(prevPos, sites[i], h.currTprob)
		}
		prevPos = sites[i]

		norm := 0.0
		for j := 0; j < n; j++ {
			pval := 0.0
			for k := 0; k < n; k++ {
				pval += fwdPrev[k] * matAt(h.currTprob, n, j, k)
			}
			fwdCur[j] = pval * eprob[j]
			norm += fwdCur[j]
		}
		for j := 0; j < n; j++ {
			fwdCur[j] /= norm
		}
	}

	prevPos = sites[nsites-1]
	for i := 0; i < nsites; i++ {
		fwdCur := fwd[(nsites-i)*n : (nsites-i+1)*n]
		eprob := eprobs[(nsites-i-1)*n : (nsites-i)*n]

		posDiff := 0
		if sites[nsites-i-1] != prevPos {
			posDiff = int(prevPos - sites[nsites-i-1] - 1)
		}
		h.setTprobForDiff(posDiff)
		if h.setTprob != nil {
			h.setTprob(sites[nsites-i-1], prevPos, h.currTprob)
		}
		prevPos = sites[nsites-i-1]

		bwdNorm := 0.0
		for j := 0; j < n; j++ {
			pval := 0.0
			for k := 0; k < n; k++ {
				pval += bwd[k] * eprob[k] * matAt(h.currTprob, n, k, j)
			}
			bwdTmp[j] = pval
			bwdNorm += pval
		}
		norm := 0.0
		for j := 0; j < n; j++ {
			bwdTmp[j] /= bwdNorm
			fwdBwd[j] = fwdCur[j] * bwdTmp[j]
			norm += fwdBwd[j]
		}
		for j := 0; j < n; j++ {
			fwdBwd[j] /= norm
			tmpGamma[j] += fwdBwd[j]
		}
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				tmpXi[k*n+j] += fwdCur[j] * bwd[k] *
					matAt(h.tprobMatrix(), n, k, j) * eprob[k] / norm
			}
		}
		for j := 0; j < n; j++ {
			fwdCur[j] = fwdBwd[j]
		}
		bwd, bwdTmp = bwdTmp, bwd
	}

	for j := 0; j < n; j++ {
		norm := 0.0
		for k := 0; k < n; k++ {
			matSet(h.currTprob, n, k, j, tmpXi[k*n+j]/tmpGamma[j])
			norm += matAt(h.currTprob, n, k, j)
		}
		for k := 0; k < n; k++ {
			matSet(h.currTprob, n, k, j, matAt(h.currTprob, n, k, j)/norm)
		}
	}
	return h.currTprob
}

// phredScore mirrors bcftools.h's phred_score: -4.3429*ln(prob),
// capped at 99, with prob==0 yielding 99.
func phredScore(prob float64) float64 {
	if prob == 0 {
		return 99
	}
	p := -4.3429 * math.Log(prob)
	if p > 99 {
		return 99
	}
	return p
}
