package edlib

// Needleman-Wunsch global alignment via straight dynamic programming, plus
// a traceback that reconstructs edlib-style opcodes (OpMatch / OpInsert /
// OpDelete / OpMismatch).
//
// The vendored bcftools edlib.c commented out the bit-parallel NW path
// (edlib.c:167-172) since --indels-cns only uses HW. We still need NW for:
// (a) generating the alignment path of HW/SHW results once we've located
//     the start/end window — the consensus caller's downstream consumers
//     walk this opcode stream the same way the C original does at
//     bam2bcf_edlib.c:866-885;
// (b) exposing ModeNW to direct callers and to our test suite, which uses
//     it to cross-check distances against textbook examples.
//
// Quadratic memory in min(len(q), len(t)) is acceptable for the read-vs-
// consensus length regime mpileup encounters (a few hundred bases each).

// nwAlign runs full DP for ModeNW with optional traceback. Inputs are
// already alphabet-transformed by edlib.go.
func nwAlign(q, t []byte, alpha int, cfg Config) (Result, error) {
	res := Result{EditDistance: -1, AlphabetLength: alpha}
	qLen := len(q)
	tLen := len(t)

	// Score matrix: (qLen+1) x (tLen+1). Row-major.
	stride := tLen + 1
	// We need the full matrix if we want a path; otherwise just rolling
	// rows would suffice, but to keep the code simple we always allocate.
	dp := make([]int32, (qLen+1)*stride)
	for i := 0; i <= qLen; i++ {
		dp[i*stride] = int32(i)
	}
	for j := 0; j <= tLen; j++ {
		dp[j] = int32(j)
	}
	for i := 1; i <= qLen; i++ {
		base := i * stride
		prev := (i - 1) * stride
		qi := q[i-1]
		for j := 1; j <= tLen; j++ {
			cost := int32(1)
			if qi == t[j-1] {
				cost = 0
			}
			diag := dp[prev+j-1] + cost
			up := dp[prev+j] + 1
			left := dp[base+j-1] + 1
			best := diag
			if up < best {
				best = up
			}
			if left < best {
				best = left
			}
			dp[base+j] = best
		}
	}
	dist := int(dp[qLen*stride+tLen])
	if cfg.K >= 0 && dist > cfg.K {
		// Behaviour matches edlib's "edit distance larger than k -> -1"
		// (edlib.c:104-106 result init, edlib.h:101-107 comment).
		return res, nil
	}
	res.EditDistance = dist
	res.EndLocations = []int{tLen - 1}
	if cfg.Task == TaskLoc || cfg.Task == TaskPath {
		res.StartLocations = []int{0}
	}
	if cfg.Task == TaskPath {
		res.Alignment = tracebackFromMatrix(dp, stride, q, t)
	}
	return res, nil
}

// tracebackFromMatrix walks an already-filled NW score matrix to produce
// the edlib opcode stream. Tie-break order matches what bcftools'
// downstream consumer assumes: prefer diagonal (match/mismatch) over
// vertical (insertion-to-query) over horizontal (insertion-to-target).
func tracebackFromMatrix(dp []int32, stride int, q, t []byte) []byte {
	i := len(q)
	j := len(t)
	// Worst case: i + j opcodes.
	ops := make([]byte, 0, i+j)
	for i > 0 || j > 0 {
		cur := dp[i*stride+j]
		if i > 0 && j > 0 {
			cost := int32(1)
			match := q[i-1] == t[j-1]
			if match {
				cost = 0
			}
			if dp[(i-1)*stride+(j-1)]+cost == cur {
				if match {
					ops = append(ops, OpMatch)
				} else {
					ops = append(ops, OpMismatch)
				}
				i--
				j--
				continue
			}
		}
		if i > 0 && dp[(i-1)*stride+j]+1 == cur {
			// Consume query but not target = deletion from target = OpDelete.
			ops = append(ops, OpDelete)
			i--
			continue
		}
		// Else: consume target but not query = insertion to target = OpInsert.
		ops = append(ops, OpInsert)
		j--
	}
	// Reverse in place.
	for l, r := 0, len(ops)-1; l < r; l, r = l+1, r-1 {
		ops[l], ops[r] = ops[r], ops[l]
	}
	return ops
}

// tracebackNW runs DP between an already-transformed query and a target
// window then reconstructs the opcode sequence. Used by the HW/SHW path
// once start/end locations are known.
func tracebackNW(q, t []byte) []byte {
	qLen := len(q)
	tLen := len(t)
	stride := tLen + 1
	dp := make([]int32, (qLen+1)*stride)
	for i := 0; i <= qLen; i++ {
		dp[i*stride] = int32(i)
	}
	for j := 0; j <= tLen; j++ {
		dp[j] = int32(j)
	}
	for i := 1; i <= qLen; i++ {
		base := i * stride
		prev := (i - 1) * stride
		qi := q[i-1]
		for j := 1; j <= tLen; j++ {
			cost := int32(1)
			if qi == t[j-1] {
				cost = 0
			}
			diag := dp[prev+j-1] + cost
			up := dp[prev+j] + 1
			left := dp[base+j-1] + 1
			best := diag
			if up < best {
				best = up
			}
			if left < best {
				best = left
			}
			dp[base+j] = best
		}
	}
	return tracebackFromMatrix(dp, stride, q, t)
}
