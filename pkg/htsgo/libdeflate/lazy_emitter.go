package libdeflate

// Port of libdeflate's deflate_compress_lazy_generic
// (reference_code/libdeflate/lib/deflate_compress.c:2604-2807) for the
// level-6 lazy parser. Slice 2 covers only the matchfinder/item-stream
// side; the dynamic-Huffman block emission and the cost-driven block
// chooser are handled in Slice 3/4. The resulting item list is the
// concatenation of the LIT / MATCH events that libdeflate's BTRACE
// instrumentation records, in the same order, so we can verify the
// matchfinder against the oracle traces directly.

import "math/bits"

// Block-split / sequence-store constants from deflate_compress.c.
const (
	softMaxBlockLength          = 300000 // deflate_compress.c:81
	minBlockLength              = 5000   // deflate_compress.c:66
	seqStoreLength              = 50000  // deflate_compress.c:93
	numLiteralObservationTypes  = 8      // deflate_compress.c:439
	numMatchObservationTypes    = 2      // deflate_compress.c:440
	numObservationTypes         = numLiteralObservationTypes + numMatchObservationTypes
	numObservationsPerBlkCheck  = 512   // deflate_compress.c:443
	recalcMinMatchLenChunkBytes = 10000 // deflate_compress.c:2625
)

// blockSplitStats mirrors struct block_split_stats (deflate_compress.c
// around the merge_new_observations definitions).
type blockSplitStats struct {
	newObservations    [numObservationTypes]uint32
	observations       [numObservationTypes]uint32
	numNewObservations uint32
	numObservations    uint32
}

func (s *blockSplitStats) init() { *s = blockSplitStats{} }

// observeLiteral mirrors observe_literal (deflate_compress.c:2110).
func (s *blockSplitStats) observeLiteral(lit byte) {
	s.newObservations[((lit>>5)&0x6)|(lit&1)]++
	s.numNewObservations++
}

// observeMatch mirrors observe_match (deflate_compress.c:2121).
func (s *blockSplitStats) observeMatch(length uint32) {
	idx := numLiteralObservationTypes
	if length >= 9 {
		idx++
	}
	s.newObservations[idx]++
	s.numNewObservations++
}

// merge folds new observations into the running totals; called when
// the split heuristic decides not to end the block yet
// (deflate_compress.c:2129).
func (s *blockSplitStats) merge() {
	for i := 0; i < numObservationTypes; i++ {
		s.observations[i] += s.newObservations[i]
		s.newObservations[i] = 0
	}
	s.numObservations += s.numNewObservations
	s.numNewObservations = 0
}

// doEndBlockCheck mirrors do_end_block_check (deflate_compress.c:2142).
// Returns true when the block should end.
func (s *blockSplitStats) doEndBlockCheck(blockLength uint32) bool {
	if s.numObservations > 0 {
		var totalDelta uint32
		for i := 0; i < numObservationTypes; i++ {
			expected := s.observations[i] * s.numNewObservations
			actual := s.newObservations[i] * s.numObservations
			var delta uint32
			if actual > expected {
				delta = actual - expected
			} else {
				delta = expected - actual
			}
			totalDelta += delta
		}
		numItems := s.numObservations + s.numNewObservations
		cutoff := s.numNewObservations * 200 / 512 * s.numObservations
		// Short-block penalty (deflate_compress.c:2187).
		if blockLength < 10000 && numItems < 8192 {
			cutoff += uint32(uint64(cutoff) * uint64(8192-numItems) / 8192)
		}
		if totalDelta+(blockLength/4096)*s.numObservations >= cutoff {
			return true
		}
	}
	s.merge()
	return false
}

// readyToCheckBlock mirrors ready_to_check_block (deflate_compress.c:2200).
func (s *blockSplitStats) readyToCheckBlock(blockBegin, next, end int32) bool {
	return s.numNewObservations >= numObservationsPerBlkCheck &&
		next-blockBegin >= minBlockLength &&
		end-next >= minBlockLength
}

// shouldEndBlock mirrors should_end_block (deflate_compress.c:2210).
func (s *blockSplitStats) shouldEndBlock(blockBegin, next, end int32) bool {
	if !s.readyToCheckBlock(blockBegin, next, end) {
		return false
	}
	return s.doEndBlockCheck(uint32(next - blockBegin))
}

// deflateFreqs tracks per-symbol frequencies. Only the literal half
// is needed for min_len recalculation; the length/offset frequencies
// will become relevant in Slice 3 when the dynamic Huffman code is
// built.
type deflateFreqs struct {
	litlen [numLitlenSyms]uint32
	offset [numOffsetSyms]uint32
}

func (f *deflateFreqs) reset() { *f = deflateFreqs{} }

// minLensTable mirrors the choose_min_match_len lookup table
// (deflate_compress.c:2299). It has 80 entries; entries past index 79
// are implicitly 3.
var minLensTable = [80]uint8{
	9, 9, 9, 9, 9, 9, 8, 8, 7, 7, 6, 6, 6, 6, 6, 6,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
}

// chooseMinMatchLen mirrors choose_min_match_len (deflate_compress.c:2296).
func chooseMinMatchLen(numUsedLiterals, maxSearchDepth uint32) uint32 {
	if numUsedLiterals >= uint32(len(minLensTable)) {
		return 3
	}
	minLen := uint32(minLensTable[numUsedLiterals])
	if maxSearchDepth < 16 {
		switch {
		case maxSearchDepth < 5:
			if minLen > 4 {
				minLen = 4
			}
		case maxSearchDepth < 10:
			if minLen > 5 {
				minLen = 5
			}
		default:
			if minLen > 7 {
				minLen = 7
			}
		}
	}
	return minLen
}

// calculateMinMatchLen mirrors calculate_min_match_len
// (deflate_compress.c:2330). It samples the first up-to-4096 bytes of
// the block to estimate literal diversity.
func calculateMinMatchLen(data []byte, maxSearchDepth uint32) uint32 {
	if len(data) < 512 {
		return minMatchLen
	}
	n := len(data)
	if n > 4096 {
		n = 4096
	}
	var used [256]bool
	for i := 0; i < n; i++ {
		used[data[i]] = true
	}
	var num uint32
	for _, u := range used {
		if u {
			num++
		}
	}
	return chooseMinMatchLen(num, maxSearchDepth)
}

// recalculateMinMatchLen mirrors recalculate_min_match_len
// (deflate_compress.c:2360). Uses the running literal frequencies.
func recalculateMinMatchLen(freqs *deflateFreqs, maxSearchDepth uint32) uint32 {
	var literalFreq uint32
	for i := 0; i < numLiterals; i++ {
		literalFreq += freqs.litlen[i]
	}
	cutoff := literalFreq >> 10
	var num uint32
	for i := 0; i < numLiterals; i++ {
		if freqs.litlen[i] > cutoff {
			num++
		}
	}
	return chooseMinMatchLen(num, maxSearchDepth)
}

// adjustMaxAndNiceLen mirrors adjust_max_and_nice_len
// (deflate_compress.c:2271).
func adjustMaxAndNiceLen(maxLen, niceLen *uint32, remaining int32) {
	if uint32(remaining) < maxMatchLen {
		*maxLen = uint32(remaining)
		if *niceLen > *maxLen {
			*niceLen = *maxLen
		}
	}
}

// chooseMaxBlockEnd mirrors choose_max_block_end (deflate_compress.c:2380).
func chooseMaxBlockEnd(blockBegin, end int32, softMax int32) int32 {
	if end-blockBegin < softMax+minBlockLength {
		return end
	}
	return blockBegin + softMax
}

// lazyBlock holds the matchfinder output for a single DEFLATE block.
// Slice 3 consumes this directly: `items` and `freqs` feed both the
// dynamic-Huffman code builder and the cost-based block chooser, while
// `length` is the uncompressed length used by the uncompressed-cost
// formula.
type lazyBlock struct {
	items  []item
	freqs  deflateFreqs
	length uint32 // uncompressed byte length
	begin  int32  // absolute start position in the input
}

// lazyEmit is a thin compatibility wrapper used by the Slice 2 oracle
// tests. It flattens lazyEmitBlocks back to a single item slice.
func lazyEmit(in []byte, level int) []item {
	blocks := lazyEmitBlocks(in, level)
	if blocks == nil {
		return nil
	}
	total := 0
	for _, b := range blocks {
		total += len(b.items)
	}
	out := make([]item, 0, total)
	for _, b := range blocks {
		out = append(out, b.items...)
	}
	return out
}

// lazyEmitBlocks runs the lazy parser at level `level` over `in` and
// returns the per-block matchfinder output: the literal/match item
// stream, per-symbol frequencies, and uncompressed block length, in
// the exact order the BTRACE harness records them. Only levels 5
// (lazy) and 6 (lazy, primary target) are covered today; lazy2 (8/9)
// is a small additional change planned for a follow-up slice.
func lazyEmitBlocks(in []byte, level int) []lazyBlock {
	if level < 2 || level > 9 {
		return nil
	}
	params := levelTable[level]
	if params.impl != implLazy && params.impl != implLazy2 {
		return nil
	}
	lazy2 := params.impl == implLazy2

	var (
		mf      hcMatchfinder
		freqs   deflateFreqs
		stats   blockSplitStats
		blocks  []lazyBlock
		items   []item
		inBase  int32
		inNext  int32
		inEnd   = int32(len(in))
		maxLen  = uint32(maxMatchLen)
		niceLen = params.niceMatchLength
	)
	if niceLen > maxLen {
		niceLen = maxLen
	}
	mf.init()
	var nextHashes [2]uint32

	for inNext != inEnd {
		// Start a new block.
		inBlockBegin := inNext
		inMaxBlockEnd := chooseMaxBlockEnd(inBlockBegin, inEnd, softMaxBlockLength)
		nextRecalc := inNext
		if d := inEnd - inNext; d < recalcMinMatchLenChunkBytes {
			nextRecalc += d
		} else {
			nextRecalc += recalcMinMatchLenChunkBytes
		}
		stats.init()
		freqs.reset()
		items = make([]item, 0, inMaxBlockEnd-inBlockBegin)

		minLen := calculateMinMatchLen(in[inNext:inMaxBlockEnd], params.maxSearchDepth)
		// Track sequence count to mirror the SEQ_STORE_LENGTH check.
		seqCount := 0

		for {
			// Periodically recalculate min_len from the running
			// literal frequencies (deflate_compress.c:2644).
			if inNext >= nextRecalc {
				minLen = recalculateMinMatchLen(&freqs, params.maxSearchDepth)
				delta := inNext - inBlockBegin
				rem := inEnd - nextRecalc
				if rem < delta {
					delta = rem
				}
				nextRecalc += delta
			}

			adjustMaxAndNiceLen(&maxLen, &niceLen, inEnd-inNext)

			// Find the longest match at the current position.
			curOffset, curLen := mf.longestMatch(
				in, &inBase, inNext,
				minLen-1, maxLen, niceLen,
				params.maxSearchDepth, &nextHashes,
			)
			if curLen < minLen || (curLen == minMatchLen && curOffset > 8192) {
				// No usable match. Emit a literal.
				items = append(items, litItem(in[inNext]))
				freqs.litlen[in[inNext]]++
				stats.observeLiteral(in[inNext])
				inNext++
			} else {
				inNext++
				// haveCurMatch loop: a match of length curLen at
				// position (inNext-1, offset=curOffset) is in hand.
				// The lazy heuristic may upgrade it to a better
				// deferred match before we commit, which is why
				// this is structured as a loop instead of a
				// straight branch.
				for {
					if curLen >= niceLen {
						items = append(items, matchItem(uint16(curLen), uint16(curOffset)))
						freqs.litlen[firstLenSym+int(lengthSlot[curLen])]++
						freqs.offset[offsetSlot(curOffset)]++
						stats.observeMatch(curLen)
						seqCount++
						mf.skipBytes(in, &inBase, inNext, inEnd, curLen-1, &nextHashes)
						inNext += int32(curLen - 1)
						break
					}

					// Lazy step: look one byte ahead with
					// half search depth.
					adjustMaxAndNiceLen(&maxLen, &niceLen, inEnd-inNext)
					nextOffset, nextLen := mf.longestMatch(
						in, &inBase, inNext,
						curLen-1, maxLen, niceLen,
						params.maxSearchDepth>>1, &nextHashes,
					)
					inNext++
					if nextLen >= curLen &&
						4*int(nextLen-curLen)+(bsr32(curOffset)-bsr32(nextOffset)) > 2 {
						lit := in[inNext-2]
						items = append(items, litItem(lit))
						freqs.litlen[lit]++
						stats.observeLiteral(lit)
						curLen = nextLen
						curOffset = nextOffset
						continue
					}

					if lazy2 {
						adjustMaxAndNiceLen(&maxLen, &niceLen, inEnd-inNext)
						next2Off, next2Len := mf.longestMatch(
							in, &inBase, inNext,
							curLen-1, maxLen, niceLen,
							params.maxSearchDepth>>2, &nextHashes,
						)
						inNext++
						if next2Len >= curLen &&
							4*int(next2Len-curLen)+(bsr32(curOffset)-bsr32(next2Off)) > 6 {
							l1 := in[inNext-3]
							l2 := in[inNext-2]
							items = append(items, litItem(l1), litItem(l2))
							freqs.litlen[l1]++
							stats.observeLiteral(l1)
							freqs.litlen[l2]++
							stats.observeLiteral(l2)
							curLen = next2Len
							curOffset = next2Off
							continue
						}
						items = append(items, matchItem(uint16(curLen), uint16(curOffset)))
						freqs.litlen[firstLenSym+int(lengthSlot[curLen])]++
						freqs.offset[offsetSlot(curOffset)]++
						stats.observeMatch(curLen)
						seqCount++
						if curLen > 3 {
							mf.skipBytes(in, &inBase, inNext, inEnd, curLen-3, &nextHashes)
							inNext += int32(curLen - 3)
						}
					} else {
						items = append(items, matchItem(uint16(curLen), uint16(curOffset)))
						freqs.litlen[firstLenSym+int(lengthSlot[curLen])]++
						freqs.offset[offsetSlot(curOffset)]++
						stats.observeMatch(curLen)
						seqCount++
						mf.skipBytes(in, &inBase, inNext, inEnd, curLen-2, &nextHashes)
						inNext += int32(curLen - 2)
					}
					break
				}
			}

			if inNext >= inMaxBlockEnd ||
				seqCount >= seqStoreLength ||
				stats.shouldEndBlock(inBlockBegin, inNext, inEnd) {
				break
			}
			if inNext == inEnd {
				break
			}
		}

		blocks = append(blocks, lazyBlock{
			items:  items,
			freqs:  freqs,
			length: uint32(inNext - inBlockBegin),
			begin:  inBlockBegin,
		})

		if inNext == inEnd {
			break
		}
	}

	return blocks
}

// bsr32 returns the 0-indexed position of the highest set bit of x.
// Mirrors libdeflate's bsr32(); the C function is undefined for x == 0
// and so is this one — callers ensure a non-zero offset.
func bsr32(x uint32) int {
	if x == 0 {
		return 0
	}
	return 31 - bits.LeadingZeros32(x)
}
