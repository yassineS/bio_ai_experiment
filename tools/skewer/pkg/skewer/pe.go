package skewer

import (
	"math"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// This file ports skewer 0.2.2's paired-end overlap analysis with overlap-based
// base error-correction, from reference_code/skewer/src/matrix.cpp and main.cpp.
//
// The three upstream pieces reproduced here are:
//
//   - cAdapter::align with bBestAlign=false (matrix.cpp:297-435): instead of
//     returning a single best alignment, it accumulates the *set* of admissible
//     adapter alignments — one per candidate cut position — into a cElementSet
//     ordered ascending by position, keeping the highest-scoring element per
//     position.  This is the candidate generator the PE path consumes.
//   - cMatrix::findAdapterWithPE (matrix.cpp:726-851): walks the merged R1/R2
//     candidate sets in increasing position order, scoring the reverse-complement
//     overlap between the mates at each candidate insert length and choosing the
//     position that maximises the joint (overlap + adapter) score.  Its output is
//     the pair of cut positions index.pos / index2.pos.
//   - cMatrix::combinePairSeqs (matrix.cpp:1117-1151): within the trimmed
//     overlap, copies the higher-quality mate's base (and quality) onto the first
//     mate — the overlap-based error-correction that makes -trimmed-pair1 differ
//     from a plain per-mate trim.
//
// The driver in main.cpp:1390-1419 ties them together: trim both mates at
// index.pos / index2.pos, run combinePairSeqs over the overlap when both cut
// positions clear minLen, and drop the pair when either falls below minLen.

// peElement mirrors the upstream ELEMENT struct for the PE candidate set: a
// normalised alignment score together with the read position at which the
// adapter alignment starts (the prospective cut position).
type peElement struct {
	score float64 // normalised score: span*dMu - rawPenalty (higher is better)
	pos   int     // adapter start position in the read (cut position)
}

// peCandidateSet is the Go analogue of cElementSet: a set keyed on pos that
// keeps the best (highest-score) element per position and yields its members in
// ascending position order.  Upstream backs this with a std::set ordered by
// idx.pos; we keep an ascending-by-pos slice and an index map, matching the
// insert/replace semantics of cElementSet::insert (matrix.cpp:144-158).
type peCandidateSet struct {
	byPos map[int]int // pos -> index into elems
	elems []peElement // kept sorted ascending by pos
}

func newPECandidateSet() *peCandidateSet {
	return &peCandidateSet{byPos: make(map[int]int)}
}

// insert adds val, replacing an existing element at the same position only when
// the new score is strictly higher (cElementSet::insert: keep the larger score,
// drop the smaller). Positions stay sorted ascending.
func (s *peCandidateSet) insert(val peElement) {
	if idx, ok := s.byPos[val.pos]; ok {
		if val.score > s.elems[idx].score {
			s.elems[idx].score = val.score
		}
		return
	}
	// Binary-insert to keep ascending order by pos.
	lo, hi := 0, len(s.elems)
	for lo < hi {
		mid := (lo + hi) / 2
		if s.elems[mid].pos < val.pos {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	s.elems = append(s.elems, peElement{})
	copy(s.elems[lo+1:], s.elems[lo:])
	s.elems[lo] = val
	// Rebuild the position index for entries at/after the insertion point.
	for i := lo; i < len(s.elems); i++ {
		s.byPos[s.elems[i].pos] = i
	}
}

// alignCandidates ports cAdapter::align with bBestAlign=false and the
// (bc >= 0) branch used by findAdapterWithPE, for the TRIM_TAIL adapters that
// the default `-m pe` mode uses. It fills `result` with the admissible adapter
// alignments (one normalised score per candidate cut position).
//
// Parameters mirror alignTrimTail: read/qual are the mate and its qualities,
// adapter the 3' adapter, minOverlap = cMatrix::iMinOverlap, errorRate =
// cMatrix::dEpsilon, indelRate = cMatrix::dEpsilonIndel. The bBestAlign=false
// path always uses minK = 1 (matrix.cpp:303) and dMu = MEAN_PENALTY (bc >= 0).
func alignCandidates(read, adapter string, qual []byte, minOverlap int, errorRate, indelRate float64, result *peCandidateSet) {
	rLen := len(read)
	length := len(adapter)
	if length == 0 {
		return
	}
	if length > maxAdapterLen {
		length = maxAdapterLen
	}

	var matchBits [cdCnt]uint64
	for code := 0; code < cdCnt; code++ {
		var bits uint64
		for i := length - 1; i >= 0; i-- {
			code2 := codeMap[adapter[i]]
			bits = (bits << 1) | chrVadp[code][code2]
		}
		matchBits[code] = ^bits
	}

	dDelta := maxPenalty // TRIM (non-AP) default, matrix.cpp:540.
	dMu := meanPenalty   // bc >= 0 ⇒ dMu = MEAN_PENALTY (matrix.cpp:304).
	dPenaltyPerErr := errorRate * meanPenalty
	bSensitive := indelRate > 0

	dMaxPenalty := dPenaltyPerErr*float64(length) + 0.001
	iMaxIndel := int(math.Ceil(indelRate * float64(length)))
	minK := 1 // bBestAlign == false (matrix.cpp:303).

	queue := make([]alignElement, 0, length+2)
	var legalBits uint64

	// Leading-gap initialisation (TRIM_TAIL branch, matrix.cpp:321-330).
	score := dDelta
	for i := 1; i < length; i++ {
		if i > iMaxIndel {
			break
		}
		queue = append(queue, alignElement{score: score, nIndel: i, pos: -i})
		legalBits = (legalBits << 1) | 1
		score += dDelta
	}

	var unbits, dnbits uint64

	for j := 0; j < rLen; j++ {
		mbits := matchBits[codeMap[read[j]]]
		var penal float64
		if len(qual) > 0 {
			penal = mismatchPenalty(qual, j)
		} else {
			penal = dMu
		}

		var es float64
		if (mbits & 0x01) == 0 {
			es = penal
		}
		queue = append(queue, alignElement{})
		copy(queue[1:], queue[0:])
		queue[0] = alignElement{score: es, nIndel: 0, pos: j}

		xbits := mbits | unbits
		dnbits <<= 1
		unbits <<= 1
		d0bits := ((dnbits + (xbits & dnbits)) ^ dnbits) | xbits
		legalBits = (legalBits << 1) | 1

		legalBits = updateColumn(&queue, d0bits, legalBits, &unbits, &dnbits, penal, dMaxPenalty, iMaxIndel, dDelta, bSensitive)

		dnbits &= d0bits
		unbits &= d0bits

		if len(queue) == length {
			// Full-window alignment: insert normalised score (matrix.cpp:374-378).
			back := queue[len(queue)-1]
			result.insert(peElement{score: float64(length)*dMu - back.score, pos: back.pos})
			queue = queue[:len(queue)-1]
		}
	}

	// Tail handling for adapters overrunning the read end (matrix.cpp:417-426).
	dMaxPenalty = dPenaltyPerErr*float64(len(queue)) + 0.001
	for i := len(queue); i >= minK; i, dMaxPenalty = i-1, dMaxPenalty-dPenaltyPerErr {
		back := queue[len(queue)-1]
		if back.score < dMaxPenalty {
			result.insert(peElement{score: float64(i)*dMu - back.score, pos: back.pos})
		}
		queue = queue[:len(queue)-1]
	}
}

// calcRevCompScoreVal ports cMatrix::CalcRevCompScore (matrix.cpp:487-522). It
// scores the n bases of seq1 starting at off1 against the reverse-complement of
// the n bases of seq2 ending at off2+n-1, and reports whether the
// quality-weighted mismatch penalty stayed within dPenaltyPerErr*n; on success
// it also returns the normalised overlap score (n*dMu - penalty). off1/off2 are
// the iStart/iStart2 offsets from findAdapterWithPE. minQLen is the *minimum* of
// the two qualities' lengths, driving the qLen>0 quality-weighting branch.
func calcRevCompScoreVal(seq1, seq2 string, qual1, qual2 []byte, off1, off2, n int, errorRate float64, minQLen int) (bool, float64) {
	dPenaltyPerErr := errorRate * meanPenalty
	if n <= 0 {
		// matrix.cpp:492-496: prefer to detect empty overlaps.
		if minQLen > 0 {
			return true, meanPenalty * float64(minQLen) / 2
		}
		return true, 0.0
	}
	dMaxPenalty := dPenaltyPerErr * float64(n)
	var penalty float64
	for i := 0; i < n; i++ {
		code := codeMap[seq1[off1+i]]
		code2 := complementCode(codeMap[seq2[off2+n-1-i]])
		base := scoringMismatch(code, code2)
		if base > 0.0 {
			if minQLen > 0 {
				p1 := mismatchPenalty(qual1, off1+i)
				p2 := mismatchPenalty(qual2, off2+n-1-i)
				if p1 <= p2 {
					base *= p1
				} else {
					base *= p2
				}
			} else {
				base *= meanPenalty
			}
			penalty += base
			if penalty > dMaxPenalty {
				return false, 0.0
			}
		}
	}
	return true, float64(n)*meanPenalty - penalty
}

// findAdapterWithPE ports cMatrix::findAdapterWithPE for the default `-m pe`
// single-pair-adapter case (one forward adapter on R1, one reverse adapter on
// R2, the indices matrix never vetoing a rev-comp-validated overlap). It returns
// the cut positions (pos, pos2) and true when an adapter/overlap is detected, or
// (rLen, rLen2, false) when not.
//
// adapter1 is applied to read1, adapter2 to read2; when adapter2 is empty the
// first adapter is shared (upstream's bShareAdapter, used when -x is given but
// -y is not).
func findAdapterWithPE(read, read2 string, qual, qual2 []byte, adapter1, adapter2 string, minOverlap int, errorRate, indelRate float64) (int, int, bool) {
	rLen := len(read)
	rLen2 := len(read2)
	qLen := len(qual)
	qLen2 := len(qual2)

	result := newPECandidateSet()
	result2 := newPECandidateSet()
	alignCandidates(read, adapter1, qual, minOverlap, errorRate, indelRate, result)
	ad2 := adapter2
	if ad2 == "" {
		ad2 = adapter1
	}
	alignCandidates(read2, ad2, qual2, minOverlap, errorRate, indelRate, result2)

	if len(result.elems) == 0 && len(result2.elems) == 0 {
		return rLen, rLen2, false
	}

	minQLen := qLen
	if qLen2 < minQLen {
		minQLen = qLen2
	}

	const intMax = int(^uint(0) >> 1)
	maxScore := -1.0
	apos := rLen
	apos2 := rLen2
	matched := false

	i1 := 0 // walk index into result.elems
	i2 := 0 // walk index into result2.elems
	for {
		pos := intMax
		if i1 < len(result.elems) {
			pos = result.elems[i1].pos
		}
		pos2 := intMax
		if i2 < len(result2.elems) {
			pos2 = result2.elems[i2].pos
		}
		if pos == intMax && pos2 == intMax {
			break
		}

		var cpos, iStart, iStart2 int
		if pos <= pos2 {
			cpos = pos
			if pos <= rLen2 {
				iStart = 0
			} else {
				iStart = pos - rLen2
			}
			iStart2 = 0
		} else {
			cpos = pos2
			iStart = 0
			if pos2 <= rLen {
				iStart2 = 0
			} else {
				iStart2 = pos2 - rLen
			}
		}

		bRevComplement, score := calcRevCompScoreVal(read, read2,
			qual, qual2, iStart, iStart2, cpos, errorRate, minQLen)

		if pos < pos2 {
			if bRevComplement {
				for i1 < len(result.elems) && result.elems[i1].pos == cpos {
					// indices[bc][0] >= 0 always holds in the default case.
					if score+result.elems[i1].score > maxScore {
						maxScore = score + result.elems[i1].score
						apos = cpos
						if cpos <= rLen2 {
							apos2 = cpos
						} else {
							apos2 = rLen2
						}
						matched = true
					}
					i1++
				}
			} else {
				for i1 < len(result.elems) && result.elems[i1].pos == cpos {
					i1++
				}
			}
		} else if pos > pos2 {
			if bRevComplement {
				for i2 < len(result2.elems) && result2.elems[i2].pos == cpos {
					// indices[0][bc2] >= 0 always holds in the default case.
					if score+result2.elems[i2].score > maxScore {
						maxScore = score + result2.elems[i2].score
						if cpos <= rLen {
							apos = cpos
						} else {
							apos = rLen
						}
						apos2 = cpos
						matched = true
					}
					i2++
				}
			} else {
				for i2 < len(result2.elems) && result2.elems[i2].pos == cpos {
					i2++
				}
			}
		} else { // pos == pos2
			if bRevComplement {
				for k1 := i1; k1 < len(result.elems) && result.elems[k1].pos == cpos; k1++ {
					for k2 := i2; k2 < len(result2.elems) && result2.elems[k2].pos == cpos; k2++ {
						// indices[bc][bc2] >= 0 always holds in the default case.
						if score+result.elems[k1].score+result2.elems[k2].score <= maxScore {
							continue
						}
						maxScore = score + result.elems[k1].score + result2.elems[k2].score
						apos = cpos
						apos2 = cpos
						matched = true
					}
				}
			}
			for i1 < len(result.elems) && result.elems[i1].pos == cpos {
				i1++
			}
			for i2 < len(result2.elems) && result2.elems[i2].pos == cpos {
				i2++
			}
		}
	}

	if !matched {
		// index.bc < 0: upstream returns false with index.pos left at rLen.
		return rLen, rLen2, false
	}
	if apos <= 0 || apos2 <= 0 {
		return 0, 0, true
	}
	return apos, apos2, true
}

// combinePairSeqs ports cMatrix::combinePairSeqs (matrix.cpp:1117-1151): the
// overlap-based error-correction.  Over the trimmed overlap it copies the
// higher-quality mate's base (and quality) from read2 onto read1, in place. Only
// read1's sequence/quality slices are mutated; read2 stays untouched. seqs/quals
// are the full-length record buffers; len1/len2 are the cut positions (pos/pos2)
// and qLen1/qLen2 the quality lengths.
func combinePairSeqs(seq1, seq2 []byte, qual1, qual2 []byte, len1, len2, qLen1, qLen2 int) {
	off1, off2 := 0, 0
	if len1 != len2 {
		if len1 > len2 {
			offset := len1 - len2
			off1 += offset
			len1 -= offset
			qLen1 -= offset
		} else {
			offset := len2 - len1
			off2 += offset
			qLen2 -= offset
		}
	}
	minQLen := qLen1
	if qLen2 < minQLen {
		minQLen = qLen2
	}
	if minQLen < len1 {
		return
	}
	for i := 0; i < len1; i++ {
		code := codeMap[seq1[off1+i]]
		code2 := complementCode(codeMap[seq2[off2+len1-1-i]])
		if qual2[off2+len1-1-i] > qual1[off1+i] {
			qual1[off1+i] = qual2[off2+len1-1-i]
			if code != code2 {
				seq1[off1+i] = characterTable[code2]
			}
		}
	}
}

// trimPairWithPE runs skewer's PE overlap analysis + error-correction for one
// pair and returns the trimmed (and, for pair1, error-corrected) records, plus
// the cut positions so the caller can apply the minLen drop rule symmetrically.
// It mirrors the driver at main.cpp:1390-1419: detect the cut positions, and —
// when both clear minLen — error-correct read1 over the overlap before trimming.
func trimPairWithPE(record1, record2 *fastq.Record, opts TrimOptions, stats *TrimStats) (out1, out2 *fastq.Record, pos1, pos2 int) {
	adapter1 := opts.Adapter3
	if adapter1 == "" {
		// No adapter configured: no PE overlap detection (matrix.cpp returns the
		// untrimmed lengths). Fall back to plain per-record trimming so any
		// 5'/quality trimming still applies.
		t1 := trimRecord(record1, opts, stats)
		t2 := trimRecord(record2, opts, stats)
		return t1, t2, len(t1.Sequence), len(t2.Sequence)
	}

	// Work on copies of the sequence/quality buffers so error-correction does
	// not mutate the caller's record (it shares the reader's backing slice).
	seq1 := append([]byte(nil), record1.Sequence...)
	seq2 := append([]byte(nil), record2.Sequence...)
	qual1 := append([]byte(nil), record1.Quality...)
	qual2 := append([]byte(nil), record2.Quality...)

	p1, p2, found := findAdapterWithPE(string(record1.Sequence), string(record2.Sequence),
		record1.Quality, record2.Quality, adapter1, opts.Adapter3Pair2, opts.MinOverlap, opts.ErrorRate, opts.IndelRate)

	if found {
		if p1 >= opts.MinLength && p2 >= opts.MinLength {
			combinePairSeqs(seq1, seq2, qual1, qual2, p1, p2, len(qual1), len(qual2))
			if stats != nil {
				stats.AdapterFound3++
			}
		}
	} else {
		p1 = len(seq1)
		p2 = len(seq2)
	}

	out1 = buildPETrimmed(record1, seq1, qual1, p1)
	out2 = buildPETrimmed(record2, seq2, qual2, p2)
	return out1, out2, p1, p2
}

// buildPETrimmed assembles a record from the (possibly corrected) seq/qual
// buffers truncated to [0, pos). The caller's minLen drop rule is applied
// separately in TrimPairedEnd, so this never returns an empty record itself.
func buildPETrimmed(record *fastq.Record, seq, qual []byte, pos int) *fastq.Record {
	if pos < 0 {
		pos = 0
	}
	if pos > len(seq) {
		pos = len(seq)
	}
	return &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    seq[:pos],
		Quality:     qual[:pos],
	}
}

// complementCode returns the complement CODE for a nucleotide code, mirroring
// the complement[] table in matrix.cpp:84-86.
func complementCode(code int) int {
	switch code {
	case cdA:
		return cdT
	case cdC:
		return cdG
	case cdG:
		return cdC
	case cdT:
		return cdA
	case cdR:
		return cdY
	case cdY:
		return cdR
	case cdS:
		return cdW
	case cdW:
		return cdS
	case cdK:
		return cdM
	case cdM:
		return cdK
	case cdB:
		return cdV
	case cdD:
		return cdH
	case cdH:
		return cdD
	case cdV:
		return cdB
	case cdN:
		return cdN
	}
	return cdNone
}

// characterTable maps a CODE back to its representative ASCII base, mirroring
// the character[] table in matrix.cpp:88-90.
var characterTable = [cdCnt]byte{
	'N', 'A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N',
}

// scoringMismatch returns the scoring[code][code2] entry from matrix.cpp:92-113.
// It is 0 for a (degenerate) match and a positive weight (1, or a fractional
// IUPAC-overlap weight, or 0.05 for N) for a mismatch. Only the >0 distinction
// and the exact weight feed CalcRevCompScore.
func scoringMismatch(code, code2 int) float64 {
	return scoringTable[code][code2]
}

// scoringTable is the verbatim scoring[][] matrix from matrix.cpp:92-113.
var scoringTable = [cdCnt][cdCnt]float64{
	//   -    A    C    G    T  | R    Y    S    W    K    M |  B    D    H    V |  N
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.05},                                           // -
	{1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0.05},                                           // A
	{1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0.05},                                           // C
	{1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0.05},                                           // G
	{1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0.05},                                           // T
	{1, 0, 1, 0, 1, 0, 1, 0.75, 0.75, 0.75, 0.75, 0.5, 0, 0.5, 0, 0.05},                           // R
	{1, 1, 0, 1, 0, 1, 0, 0.75, 0.75, 0.75, 0.75, 0, 0.5, 0, 0.5, 0.05},                           // Y
	{1, 1, 0, 0, 1, 0.75, 0.75, 0, 1, 0.75, 0.75, 0, 0.5, 0.5, 0, 0.05},                           // S
	{1, 0, 1, 1, 0, 0.75, 0.75, 1, 0, 0.75, 0.75, 0.5, 0, 0, 0.5, 0.05},                           // W
	{1, 1, 1, 0, 0, 0.75, 0.75, 0.75, 0.75, 0, 1, 0, 0, 0.5, 0.5, 0.05},                           // K
	{1, 0, 0, 1, 1, 0.75, 0.75, 0.75, 0.75, 1, 0, 0.5, 0.5, 0, 0, 0.05},                           // M
	{1, 1, 0, 0, 0, 0.5, 0, 0, 0.5, 0, 0.5, 0, 0.4, 0.4, 0.4, 0.05},                               // B
	{1, 0, 1, 0, 0, 0, 0.5, 0.5, 0, 0, 0.5, 0.4, 0, 0.4, 0.4, 0.05},                               // D
	{1, 0, 0, 1, 0, 0.5, 0, 0.5, 0, 0.5, 0, 0.4, 0.4, 0, 0.4, 0.05},                               // H
	{1, 0, 0, 0, 1, 0, 0.5, 0, 0.5, 0.5, 0, 0.4, 0.4, 0.4, 0, 0.05},                               // V
	{0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0}, // N
}
