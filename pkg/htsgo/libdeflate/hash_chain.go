package libdeflate

import "encoding/binary"

// Pure-Go port of libdeflate's hc_matchfinder
// (reference_code/libdeflate/lib/hc_matchfinder.h). The data
// representation mirrors the C struct byte-for-byte, including the use
// of a signed 16-bit position type (mf_pos_t) so that values "below" the
// sliding-window base sort as out-of-bounds via a simple comparison
// against `cutoff = curPos - windowSize`.

// Matchfinder constants. WINDOW_ORDER and HASH3/HASH4 orders mirror
// matchfinder_common.h:47 and hc_matchfinder.h:112-113.
const (
	hcWindowOrder       = 15
	hcWindowSize        = 1 << hcWindowOrder // 32768
	hcWindowMask        = hcWindowSize - 1
	hcHash3Order        = 15
	hcHash4Order        = 16
	hcHash3BucketCount  = 1 << hcHash3Order
	hcHash4BucketCount  = 1 << hcHash4Order
	hcMatchfinderInit   = int32(-hcWindowSize)
	hcUnalignedFastWord = 4 // we use 4-byte word compares (UNALIGNED_ACCESS_IS_FAST = true on amd64/arm64)
)

// hcMatchfinder is the hash-chain matchfinder. Positions are stored as
// signed int32 even though libdeflate uses int16; the extra width
// simplifies the Go port without changing observable behavior because
// every comparison treats values < -hcWindowSize as "out of range".
// When the sliding window advances we still saturate at
// -hcWindowSize, matching matchfinder_rebase().
type hcMatchfinder struct {
	hash3Tab [hcHash3BucketCount]int32
	hash4Tab [hcHash4BucketCount]int32
	nextTab  [hcWindowSize]int32
}

// init resets every hash/chain entry to the sentinel value used by
// libdeflate's matchfinder_init (matchfinder_common.h:108).
func (mf *hcMatchfinder) init() {
	for i := range mf.hash3Tab {
		mf.hash3Tab[i] = hcMatchfinderInit
	}
	for i := range mf.hash4Tab {
		mf.hash4Tab[i] = hcMatchfinderInit
	}
	for i := range mf.nextTab {
		mf.nextTab[i] = hcMatchfinderInit
	}
}

// slideWindow subtracts hcWindowSize from every stored position with
// signed saturation at -hcWindowSize. Mirrors matchfinder_rebase()
// (matchfinder_common.h:137). The 32768-byte fast path collapses to a
// simple branchless sub on int16; we use the general branched form.
func (mf *hcMatchfinder) slideWindow() {
	rebase := func(v int32) int32 {
		if v >= 0 {
			return v - hcWindowSize
		}
		return -hcWindowSize
	}
	for i, v := range mf.hash3Tab {
		mf.hash3Tab[i] = rebase(v)
	}
	for i, v := range mf.hash4Tab {
		mf.hash4Tab[i] = rebase(v)
	}
	for i, v := range mf.nextTab {
		mf.nextTab[i] = rebase(v)
	}
}

// lzHash mirrors lz_hash() in matchfinder_common.h:169. The product
// must be taken modulo 2^32 (Go's uint32 multiplication wraps for us).
func lzHash(seq uint32, numBits uint) uint32 {
	return (seq * 0x1E35A7BD) >> (32 - numBits)
}

// loadU24 loads three bytes at p as the low 24 bits of a uint32, in
// platform-LE order (matchfinder_common.h:35 loaded_u32_to_u24 +
// load_u32_unaligned). We assume LE everywhere.
func loadU24(p []byte) uint32 {
	return binary.LittleEndian.Uint32(p) & 0xFFFFFF
}

// loadU32 reads four bytes at p as little-endian.
func loadU32(p []byte) uint32 {
	return binary.LittleEndian.Uint32(p)
}

// lzExtend returns the length of the byte run that matches between
// in[strPos:] and in[matchPos:], starting from startLen and capped at
// maxLen. Mirrors lz_extend() (matchfinder_common.h:179) but in
// straight bytewise Go; the 4-word unrolled compare in C is a
// micro-optimization, not a correctness invariant.
func lzExtend(in []byte, strPos, matchPos int32, startLen, maxLen uint32) uint32 {
	s := int(strPos)
	m := int(matchPos)
	n := int(maxLen)
	i := int(startLen)
	for i < n && in[m+i] == in[s+i] {
		i++
	}
	return uint32(i)
}

// longestMatch finds the longest match for the sequence at in[inNext:]
// that is longer than bestLen. Mirrors hc_matchfinder_longest_match()
// (hc_matchfinder.h:182). The caller passes the precomputed
// nextHashes (hash codes for the *current* position); we update them in
// place to the hashes for the *next* position before returning.
//
// Returns (offset, length). When no match longer than bestLen exists
// the returned length equals bestLen and offset is 0 (callers always
// gate on `length > prevLen` so the offset is undefined in that case).
//
// inBase is the absolute index of cur_pos == 0 in the matchfinder's
// coordinate space; it is updated (and pointer-equivalent in the C
// code) when the window slides.
func (mf *hcMatchfinder) longestMatch(
	in []byte,
	inBase *int32,
	inNext int32,
	bestLen uint32,
	maxLen uint32,
	niceLen uint32,
	maxSearchDepth uint32,
	nextHashes *[2]uint32,
) (offset uint32, length uint32) {

	depthRemaining := maxSearchDepth
	bestMatchPtr := inNext
	curPos := inNext - *inBase

	if curPos == hcWindowSize {
		mf.slideWindow()
		*inBase += hcWindowSize
		curPos = 0
	}

	cutoff := curPos - hcWindowSize

	// hc_matchfinder.h:214 — we need 4 bytes at in_next+1 for the
	// next-position hash, so bail out for max_len < 5.
	if maxLen < 5 {
		offset = uint32(inNext - bestMatchPtr)
		return offset, bestLen
	}

	hash3 := nextHashes[0]
	hash4 := nextHashes[1]

	curNode3 := mf.hash3Tab[hash3]
	curNode4 := mf.hash4Tab[hash4]

	// Insert current sequence (hc_matchfinder.h:227-232).
	mf.hash3Tab[hash3] = curPos
	mf.hash4Tab[hash4] = curPos
	mf.nextTab[curPos&hcWindowMask] = curNode4

	// Precompute next-position hashes (hc_matchfinder.h:235-237).
	nextHashSeq := loadU32(in[inNext+1:])
	nextHashes[0] = lzHash(nextHashSeq&0xFFFFFF, hcHash3Order)
	nextHashes[1] = lzHash(nextHashSeq, hcHash4Order)

	inBaseVal := *inBase

	if bestLen < 4 {
		// Length-3 / length-4 entry path (hc_matchfinder.h:241).
		if curNode3 <= cutoff {
			offset = uint32(inNext - bestMatchPtr)
			return offset, bestLen
		}

		seq4 := loadU32(in[inNext:])

		if bestLen < 3 {
			matchPtr := inBaseVal + curNode3
			if loadU24(in[matchPtr:]) == (seq4 & 0xFFFFFF) {
				bestLen = 3
				bestMatchPtr = matchPtr
			}
		}

		if curNode4 <= cutoff {
			offset = uint32(inNext - bestMatchPtr)
			return offset, bestLen
		}

		// First length-4 search loop (hc_matchfinder.h:263).
		// Walk the chain until we find a node whose first 4 bytes
		// match seq4, then fall through to the length-5+ extension
		// loop below.
		for {
			matchPtr := inBaseVal + curNode4
			if loadU32(in[matchPtr:]) == seq4 {
				bestMatchPtr = matchPtr
				bestLen = lzExtend(in, inNext, matchPtr, 4, maxLen)
				if bestLen >= niceLen {
					offset = uint32(inNext - bestMatchPtr)
					return offset, bestLen
				}
				curNode4 = mf.nextTab[curNode4&hcWindowMask]
				if curNode4 <= cutoff {
					offset = uint32(inNext - bestMatchPtr)
					return offset, bestLen
				}
				depthRemaining--
				if depthRemaining == 0 {
					offset = uint32(inNext - bestMatchPtr)
					return offset, bestLen
				}
				break
			}
			curNode4 = mf.nextTab[curNode4&hcWindowMask]
			if curNode4 <= cutoff {
				offset = uint32(inNext - bestMatchPtr)
				return offset, bestLen
			}
			depthRemaining--
			if depthRemaining == 0 {
				offset = uint32(inNext - bestMatchPtr)
				return offset, bestLen
			}
		}
		// Fall through to the length-5+ extension loop.
	} else {
		// Caller already has a length >= 4 match (hc_matchfinder.h:284).
		if curNode4 <= cutoff || bestLen >= niceLen {
			offset = uint32(inNext - bestMatchPtr)
			return offset, bestLen
		}
	}

	// Length >= 5 extension loop (hc_matchfinder.h:291). Tight
	// inner loop: check (last 4 bytes) and (first 4 bytes) before
	// calling lz_extend.
	for {
		var matchPtr int32
		for {
			matchPtr = inBaseVal + curNode4
			lastWordStr := loadU32(in[int32(inNext)+int32(bestLen)-3:])
			lastWordMatch := loadU32(in[matchPtr+int32(bestLen)-3:])
			if lastWordMatch == lastWordStr && loadU32(in[matchPtr:]) == loadU32(in[inNext:]) {
				break
			}
			curNode4 = mf.nextTab[curNode4&hcWindowMask]
			if curNode4 <= cutoff {
				offset = uint32(inNext - bestMatchPtr)
				return offset, bestLen
			}
			depthRemaining--
			if depthRemaining == 0 {
				offset = uint32(inNext - bestMatchPtr)
				return offset, bestLen
			}
		}

		// UNALIGNED_ACCESS_IS_FAST path: we already verified 4
		// bytes match at offset 0, so start the extension at 4.
		length = lzExtend(in, inNext, matchPtr, 4, maxLen)
		if length > bestLen {
			bestLen = length
			bestMatchPtr = matchPtr
			if bestLen >= niceLen {
				offset = uint32(inNext - bestMatchPtr)
				return offset, bestLen
			}
		}
		curNode4 = mf.nextTab[curNode4&hcWindowMask]
		if curNode4 <= cutoff {
			offset = uint32(inNext - bestMatchPtr)
			return offset, bestLen
		}
		depthRemaining--
		if depthRemaining == 0 {
			offset = uint32(inNext - bestMatchPtr)
			return offset, bestLen
		}
	}
}

// skipBytes advances the matchfinder by count positions without
// searching for matches. Mirrors hc_matchfinder_skip_bytes()
// (hc_matchfinder.h:361). count must be > 0 and count+5 must not
// exceed the remaining input.
func (mf *hcMatchfinder) skipBytes(
	in []byte,
	inBase *int32,
	inNext int32,
	inEnd int32,
	count uint32,
	nextHashes *[2]uint32,
) {
	if count+5 > uint32(inEnd-inNext) {
		return
	}

	curPos := inNext - *inBase
	hash3 := nextHashes[0]
	hash4 := nextHashes[1]
	remaining := count
	for {
		if curPos == hcWindowSize {
			mf.slideWindow()
			*inBase += hcWindowSize
			curPos = 0
		}
		mf.hash3Tab[hash3] = curPos
		mf.nextTab[curPos&hcWindowMask] = mf.hash4Tab[hash4]
		mf.hash4Tab[hash4] = curPos

		inNext++
		nextHashSeq := loadU32(in[inNext:])
		hash3 = lzHash(nextHashSeq&0xFFFFFF, hcHash3Order)
		hash4 = lzHash(nextHashSeq, hcHash4Order)
		curPos++
		remaining--
		if remaining == 0 {
			break
		}
	}
	nextHashes[0] = hash3
	nextHashes[1] = hash4
}
