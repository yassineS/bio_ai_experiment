package skewer

import "math"

// This file is a faithful Go port of skewer's bit-parallel k-difference
// adapter aligner (cAdapter::align with bBestAlign=true, TRIM_TAIL mode) from
// reference_code/skewer/src/matrix.cpp:297-435, together with the supporting
// nucleotide-code tables (matrix.cpp:38-136) and the UPDATE_COLUMN inner loop
// (matrix.cpp:224-295).
//
// The earlier gap-free scanner (findAdapterWithQual) collapsed upstream's
// indel-aware engine to a no-indel approximation. That approximation diverged
// from upstream on reads whose 3' adapter is corrupted by interior N bases:
// upstream's bit-vector keeps propagating the best-scoring alignment rightward
// across indel moves, so the reported cut position (INDEX.pos) lands further
// into the read than a naive leftmost gap-free match would. Reproducing the
// engine verbatim is the only way to match upstream's cut byte-for-byte. See
// the package README / parity notes for the worked example (read232925 in the
// medium fixture).

// CODE values mirror the CODE enum in matrix.h:65-84.
const (
	cdNone = 0
	cdA    = 1
	cdC    = 2
	cdG    = 3
	cdT    = 4
	cdR    = 5
	cdY    = 6
	cdS    = 7
	cdW    = 8
	cdK    = 9
	cdM    = 10
	cdB    = 11
	cdD    = 12
	cdH    = 13
	cdV    = 14
	cdN    = 15
	cdCnt  = cdN + 1
)

// codeMap maps an ASCII byte to its CODE, mirroring matrix.cpp:38-59.
var codeMap [256]int

// chrVadp[readCode][adapterCode] is 1 when the read base mismatches the adapter
// base and 0 when it matches, verbatim from matrix.cpp:115-136.
var chrVadp = [cdCnt][cdCnt]uint64{
	//      -  A  C  G  T  R  Y  S  W  K  M  B  D  H  V  N
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0}, // -
	{1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0}, // A
	{1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0}, // C
	{1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0}, // G
	{1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0}, // T
	{1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0}, // R
	{1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0}, // Y
	{1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0}, // S
	{1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0}, // W
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0}, // K
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0}, // M
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0}, // B
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0}, // D
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0}, // H
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0}, // V
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0}, // N
}

func init() {
	for i := range codeMap {
		codeMap[i] = cdNone
	}
	// 0x41..0x5A (A..Z) and the lower-case aliases, mirroring matrix.cpp:45-49.
	upper := map[byte]int{
		'A': cdA, 'B': cdB, 'C': cdC, 'D': cdD, 'G': cdG, 'H': cdH,
		'K': cdK, 'M': cdM, 'N': cdN, 'R': cdR, 'S': cdS, 'T': cdT,
		'U': cdT, 'V': cdV, 'W': cdW, 'Y': cdY,
	}
	for b, c := range upper {
		codeMap[b] = c
		codeMap[b+32] = c
	}
}

// maxAdapterLen mirrors MAX_ADAPTER_LEN (common.h:48); the bit-vector also caps
// the packed adapter length at 64 bits.
const maxAdapterLen = 64

// epsilonHead is EPSILON from matrix.cpp:141 (MIN_PENALTY/10).
const epsilonHead = minPenalty / 10.0

// alignElement mirrors the ELEMENT struct (matrix.h:43-47); only the fields the
// aligner reads are kept.
type alignElement struct {
	score  float64
	nIndel int
	pos    int
}

// alignTrimTail ports cAdapter::align for the TRIM_TAIL / bBestAlign case used
// by 3'-adapter trimming. It returns the cut position (read index at which the
// adapter starts, i.e. INDEX.pos) and true when an adapter is detected, or
// (rLen, false) when none is found — matching cMatrix::findAdapter's
// initialisation of index.pos = rLen.
//
// Parameters:
//   - read, qual: the read sequence and its Phred-encoded qualities (qual may
//     be empty, mirroring the qLen==0 fallback to dMu).
//   - adapter: the 3' adapter sequence.
//   - minOverlap: cMatrix::iMinOverlap.
//   - errorRate: cMatrix::dEpsilon (the -r value).
//   - indelRate: cMatrix::dEpsilonIndel (the -d value).
func alignTrimTail(read, adapter string, qual []byte, minOverlap int, errorRate, indelRate float64) (int, bool) {
	rLen := len(read)
	length := len(adapter)
	if length == 0 {
		return rLen, false
	}
	if length > maxAdapterLen {
		length = maxAdapterLen
	}

	// Per-adapter match bits, mirroring cAdapter::Init (matrix.cpp:178-188):
	// matchBits[code] has bit b set when read code `code` MATCHES adapter[b]
	// (the table is inverted with ^bits, and chrVadp stores mismatches).
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
	dMu := meanPenalty   // cMatrix::dMu, matrix.cpp:473.
	dPenaltyPerErr := errorRate * meanPenalty
	bSensitive := indelRate > 0

	dMaxPenalty := dPenaltyPerErr*float64(length) + 0.001
	iMaxIndel := int(math.Ceil(indelRate * float64(length)))
	minK := minOverlap
	if minOverlap >= length-iMaxIndel+1 {
		minK = length - iMaxIndel + 1
	}

	determined := false
	var elem alignElement

	// queue is the sliding deque of partial alignments. queue[0] is the front
	// (most recently pushed); queue[len-1] is the back.
	queue := make([]alignElement, 0, length+2)
	var legalBits uint64

	// Initialise the deque with the leading gap states (TRIM_TAIL branch,
	// matrix.cpp:321-330).
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
		// push_front: prepend the new element.
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
			// Full-window alignment available (main loop, matrix.cpp:355-380).
			back := queue[len(queue)-1]
			elem = back
			dMaxPenalty = elem.score + epsilonHead // TRIM_TAIL adds EPSILON.
			elem.score = (float64(length)*dMu - elem.score) / (float64(length) + 1)
			determined = true
			if dMaxPenalty == 0 {
				break
			}
			queue = queue[:len(queue)-1]
		}
	}

	if dMaxPenalty > 0 {
		// Tail handling for adapters that overrun the read end
		// (matrix.cpp:384-397).
		dMaxPenalty = dPenaltyPerErr*float64(len(queue)) + 0.001
		for i := len(queue); i >= minK; i, dMaxPenalty = i-1, dMaxPenalty-dPenaltyPerErr {
			if dMaxPenalty <= 0 {
				break
			}
			back := queue[len(queue)-1]
			if back.score < dMaxPenalty {
				if !determined || (float64(i)*dMu-back.score) > elem.score*(float64(i)+1) {
					elem = back
					dMaxPenalty = elem.score
					elem.score = (float64(i)*dMu - elem.score) / (float64(i) + 1)
					determined = true
				}
			}
			queue = queue[:len(queue)-1]
		}
	}

	if !determined {
		return rLen, false
	}
	return elem.pos, true
}

// updateColumn ports cAdapter::UPDATE_COLUMN (matrix.cpp:224-295). It mutates
// the deque (*queue) and the unbits/dnbits accumulators in place, including the
// trailing pop_back of inadmissible elements, and returns the updated
// legalBits.
func updateColumn(queue *[]alignElement, d0bits, legalBits uint64, unbits, dnbits *uint64, penal, dMaxPenalty float64, iMaxIndel int, dDelta float64, bSensitive bool) uint64 {
	q := *queue
	bits := (^legalBits) | d0bits
	bits >>= 1
	i := 1
	for ; i < len(q)-1; i++ {
		if (bits & 0x01) == 0 {
			if bSensitive {
				s := q[i].score + (penal - dDelta)
				if q[i-1].score < s && q[i-1].nIndel < iMaxIndel {
					if q[i+1].score < s && q[i+1].nIndel < iMaxIndel {
						if q[i-1].score < q[i+1].score {
							q[i] = q[i-1]
							*dnbits |= 1 << uint(i-1)
						} else {
							q[i] = q[i+1]
							*unbits |= 1 << uint(i+1)
						}
					} else {
						q[i] = q[i-1]
						*dnbits |= 1 << uint(i-1)
					}
					q[i].nIndel++
				} else {
					if q[i+1].score < s && q[i+1].nIndel < iMaxIndel {
						q[i] = q[i+1]
						*unbits |= 1 << uint(i+1)
						q[i].nIndel++
					} else {
						q[i].score = s
					}
				}
				q[i].score += dDelta
			} else {
				q[i].score += penal
			}
			if q[i].score >= dMaxPenalty {
				legalBits &= ^(uint64(1) << uint(i))
			}
		}
		bits >>= 1
	}
	if len(q) > 1 {
		if (bits & 0x01) == 0 {
			if bSensitive {
				if q[i-1].nIndel < iMaxIndel && q[i-1].score+dDelta < q[i].score+penal {
					q[i] = q[i-1]
					*dnbits |= 1 << uint(i-1)
					q[i].score += dDelta
					q[i].nIndel++
				} else {
					q[i].score += penal
				}
			} else {
				q[i].score += penal
			}
			if q[i].score >= dMaxPenalty {
				legalBits &= ^(uint64(1) << uint(i))
			}
		}
		for ; i > 0; i-- {
			if q[len(q)-1].score < dMaxPenalty {
				break
			}
			q = q[:len(q)-1]
		}
	}
	*queue = q
	return legalBits
}
