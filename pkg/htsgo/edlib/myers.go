package edlib

// Bit-parallel Myers' algorithm — semi-global variant. Direct port of
// myersCalcEditDistanceSemiGlobal at reference_code/bcftools/edlib.c:415-577.
// Variable names track the original almost identically; comments call out
// where the Go port deviates from the C.

const (
	wordSize = 64
	highBit  = uint64(1) << 63
	maxPos   = 100 // MAX_POS, edlib.c:451
	// strongReduceNum matches edlib.c:439. Every this-many columns we
	// run a more expensive band-shrinking pass.
	strongReduceNum = 2048
)

type block struct {
	P     uint64 // Pvin: +1 deltas
	M     uint64 // Mvin: -1 deltas
	score int    // score of last cell in this block
}

// ceilDiv mirrors edlib.c:351-353. Inputs must be non-negative.
func ceilDiv(x, y int) int {
	if x%y != 0 {
		return x/y + 1
	}
	return x / y
}

// buildPeq matches edlib.c:251-282. Peq is laid out as
// (alphabetLength + 1) * maxBlocks words; the trailing symbol acts as a
// wildcard (all-ones). Each block's bit i is 1 iff query[r] equals the
// symbol where r = b*WORD_SIZE + i (or r >= queryLength, which also yields 1
// to mimic the C's "padding equals everything" optimisation).
func buildPeq(alphabet int, query []byte) []uint64 {
	qLen := len(query)
	maxBlocks := ceilDiv(qLen, wordSize)
	peq := make([]uint64, (alphabet+1)*maxBlocks)
	for sym := 0; sym < alphabet; sym++ {
		for b := 0; b < maxBlocks; b++ {
			var w uint64
			// Build from top of block downwards so the bit shifts
			// match the C loop (edlib.c:264-269).
			for r := (b+1)*wordSize - 1; r >= b*wordSize; r-- {
				w <<= 1
				if r >= qLen || query[r] == byte(sym) {
					w |= 1
				}
			}
			peq[sym*maxBlocks+b] = w
		}
	}
	// Wildcard row, edlib.c:273-279.
	wildcardBase := alphabet * maxBlocks
	for b := 0; b < maxBlocks; b++ {
		peq[wildcardBase+b] = ^uint64(0)
	}
	return peq
}

// calculateBlock advances one block (64-cell column tile) by one column.
// Faithful port of edlib.c:310-345. Returns hout (-1, 0, or +1) and writes
// the new Pv/Mv via the pointer parameters in the C original; here we
// return them.
func calculateBlock(Pv, Mv, Eq uint64, hin int) (newPv, newMv uint64, hout int) {
	// hin encoding: 1 -> 0...01, 0 -> 0...00, -1 -> 1...11 (two's complement).
	// Replicate the C trick at edlib.c:317 to extract sign.
	hinIsNeg := uint64(hin>>2) & 1

	Xv := Eq | Mv
	Eq |= hinIsNeg
	Xh := (((Eq & Pv) + Pv) ^ Pv) | Eq

	Ph := Mv | ^(Xh | Pv)
	Mh := Pv & Xh

	hout = int((Ph & highBit) >> (wordSize - 1))
	hout -= int((Mh & highBit) >> (wordSize - 1))

	Ph <<= 1
	Mh <<= 1

	Mh |= hinIsNeg
	// (hin+1)>>1 -> 1 if hin>0, else 0. Matches edlib.c:339.
	Ph |= uint64((hin + 1) >> 1)

	newPv = Mh | ^(Xv | Ph)
	newMv = Ph & Xv
	return
}

// getBlockCellValues returns the 64 cell scores in a block, bottom cell
// first. Direct translation of edlib.c:364-376. Used by the band-reduction
// path (allBlockCellsLarger) and by the tail sweep after the column loop.
func getBlockCellValues(bl block) [wordSize]int {
	var scores [wordSize]int
	score := bl.score
	mask := highBit
	for i := 0; i < wordSize-1; i++ {
		scores[i] = score
		if bl.P&mask != 0 {
			score--
		}
		if bl.M&mask != 0 {
			score++
		}
		mask >>= 1
	}
	scores[wordSize-1] = score
	return scores
}

func allBlockCellsLarger(bl block, k int) bool {
	scores := getBlockCellValues(bl)
	for i := 0; i < wordSize; i++ {
		if scores[i] <= k {
			return false
		}
	}
	return true
}

// myersSemiGlobal mirrors myersCalcEditDistanceSemiGlobal at
// edlib.c:415-577. Returns bestScore (-1 if no path with cost <= k) and the
// list of end positions in target where bestScore was achieved.
func myersSemiGlobal(peq []uint64, w, maxBlocks, qLen int, target []byte, k int, mode Mode) (bestScore int, positions []int) {
	tLen := len(target)
	bestScore = -1

	// Ukkonen band: [firstBlock..lastBlock] inclusive.
	firstBlock := 0
	lastBlock := minInt(ceilDiv(k+1, wordSize), maxBlocks) - 1

	if mode == ModeHW {
		// For HW the answer is never larger than qLen (edlib.c:432-435).
		k = minInt(qLen, k)
	}

	blocks := make([]block, maxBlocks)
	// Initialise blocks 0..lastBlock with score=(b+1)*WORD_SIZE, P=all-ones,
	// M=0 (edlib.c:442-448).
	for b := 0; b <= lastBlock; b++ {
		blocks[b].score = (b + 1) * wordSize
		blocks[b].P = ^uint64(0)
		blocks[b].M = 0
	}

	// startHout encodes whether the leading gap into the query is free.
	// HW => 0 (free), SHW => +1 per column (edlib.c:454).
	startHout := 0
	if mode == ModeSHW {
		startHout = 1
	}

	posBuf := make([]int, 0, maxPos)

	for c := 0; c < tLen; c++ {
		peqC := peq[int(target[c])*maxBlocks:]

		// ---- Calculate column (edlib.c:459-468). ----
		hout := startHout
		for b := firstBlock; b <= lastBlock; b++ {
			var p, m uint64
			p, m, hout = calculateBlock(blocks[b].P, blocks[b].M, peqC[b], hout)
			blocks[b].P = p
			blocks[b].M = m
			blocks[b].score += hout
		}
		lastBlockScore := blocks[lastBlock].score

		// ---- Adjust Ukkonen band (edlib.c:472-483). ----
		if (lastBlock < maxBlocks-1) && (lastBlockScore-hout <= k) &&
			((peqC[lastBlock+1]&1) != 0 || hout < 0) {
			lastBlock++
			blocks[lastBlock].P = ^uint64(0)
			blocks[lastBlock].M = 0
			var p, m uint64
			p, m, ho := calculateBlock(blocks[lastBlock].P, blocks[lastBlock].M, peqC[lastBlock], hout)
			blocks[lastBlock].P = p
			blocks[lastBlock].M = m
			blocks[lastBlock].score = blocks[lastBlock-1].score - hout + wordSize + ho
		} else {
			for lastBlock >= firstBlock && blocks[lastBlock].score >= k+wordSize {
				lastBlock--
			}
		}

		// Strong reduction (edlib.c:489-493).
		if c%strongReduceNum == 0 {
			for lastBlock >= 0 && lastBlock >= firstBlock && allBlockCellsLarger(blocks[lastBlock], k) {
				lastBlock--
			}
		}

		// For HW, the first block always remains a candidate
		// (edlib.c:498-500).
		if mode == ModeHW && lastBlock == -1 {
			lastBlock++
		}

		// Increase firstBlock if possible (not applicable to HW)
		// (edlib.c:503-512).
		if mode != ModeHW {
			for firstBlock <= lastBlock && blocks[firstBlock].score >= k+wordSize {
				firstBlock++
			}
			if c%strongReduceNum == 0 {
				for firstBlock <= lastBlock && allBlockCellsLarger(blocks[firstBlock], k) {
					firstBlock++
				}
			}
		}

		// Band exhausted — finish (edlib.c:515-524).
		if lastBlock < firstBlock {
			if bestScore != -1 {
				positions = append([]int(nil), posBuf...)
			}
			return bestScore, positions
		}

		// ---- Update best score (edlib.c:528-544). ----
		if lastBlock == maxBlocks-1 {
			colScore := blocks[lastBlock].score
			if colScore <= k {
				if bestScore == -1 || colScore <= bestScore {
					if colScore != bestScore {
						posBuf = posBuf[:0]
						bestScore = colScore
						k = bestScore
					}
					if len(posBuf) < maxPos {
						posBuf = append(posBuf, c-w)
					}
				}
			}
		}
	}

	// Tail sweep: scan the last block's cells for any score <= k that
	// corresponds to a target position past the column loop's window
	// (edlib.c:552-566).
	if lastBlock == maxBlocks-1 {
		blockScores := getBlockCellValues(blocks[lastBlock])
		for i := 0; i < w; i++ {
			colScore := blockScores[i+1]
			if colScore <= k && (bestScore == -1 || colScore <= bestScore) {
				if colScore != bestScore {
					posBuf = posBuf[:0]
					bestScore = colScore
					k = bestScore
				}
				if len(posBuf) < maxPos {
					posBuf = append(posBuf, tLen-w+i)
				}
			}
		}
	}

	if bestScore != -1 {
		positions = append([]int(nil), posBuf...)
	}
	return bestScore, positions
}

// reverseCopy returns a fresh reversed copy of in. Mirrors
// createReverseCopy at edlib.c:289-295.
func reverseCopy(in []byte) []byte {
	n := len(in)
	out := make([]byte, n)
	for i := 0; i < n; i++ {
		out[i] = in[n-1-i]
	}
	return out
}

// semiGlobalAlign drives the HW/SHW Myers algorithm and (optionally) finds
// start locations and the alignment path. Mirrors the relevant slice of
// edlibAlign at edlib.c:142-230 for the semi-global branches.
func semiGlobalAlign(q, t []byte, alpha int, cfg Config) (Result, error) {
	res := Result{EditDistance: -1, AlphabetLength: alpha}
	qLen := len(q)

	maxBlocks := ceilDiv(qLen, wordSize)
	w := maxBlocks*wordSize - qLen
	peq := buildPeq(alpha, q)

	dynamicK := cfg.K < 0
	k := cfg.K
	if dynamicK {
		k = wordSize
	}

	var bestScore int
	var positions []int
	for {
		bestScore, positions = myersSemiGlobal(peq, w, maxBlocks, qLen, t, k, cfg.Mode)
		if !dynamicK || bestScore != -1 {
			break
		}
		k *= 2
	}

	res.EditDistance = bestScore
	if bestScore < 0 {
		return res, nil
	}
	res.EndLocations = positions

	if cfg.Task == TaskLoc || cfg.Task == TaskPath {
		res.StartLocations = make([]int, len(positions))
		if cfg.Mode == ModeHW {
			// Reverse-scan to find start positions
			// (edlib.c:187-223).
			rT := reverseCopy(t)
			rQ := reverseCopy(q)
			rPeq := buildPeq(alpha, rQ)
			for i, endLoc := range positions {
				if endLoc == -1 {
					res.StartLocations[i] = 0
					continue
				}
				tail := rT[len(t)-endLoc-1:][:endLoc+1]
				_, posSHW := myersSemiGlobal(rPeq, w, maxBlocks, qLen, tail, bestScore, ModeSHW)
				if len(posSHW) == 0 {
					res.StartLocations[i] = 0
				} else {
					// Take last position to mirror C's
					// "ensures alignment doesn't start with
					// insertions" comment at edlib.c:215-217.
					res.StartLocations[i] = endLoc - posSHW[len(posSHW)-1]
				}
			}
		} else {
			// SHW: start is always 0.
			for i := range res.StartLocations {
				res.StartLocations[i] = 0
			}
		}
	}

	if cfg.Task == TaskPath && len(res.StartLocations) > 0 {
		start := res.StartLocations[0]
		end := res.EndLocations[0]
		if end < start {
			// Empty alignment — query did not produce any opcodes
			// within target. Leave Alignment nil.
			return res, nil
		}
		sub := t[start : end+1]
		res.Alignment = tracebackNW(q, sub)
	}
	return res, nil
}
