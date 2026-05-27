package libdeflate

// Dynamic-Huffman block emission and cost-based block-type chooser,
// ported from libdeflate's deflate_flush_block /
// deflate_precompute_huffman_header / deflate_compute_precode_items
// (reference_code/libdeflate/lib/deflate_compress.c around lines
// 1482-2038). The implementation is intentionally narrow: it covers
// only what's needed for the level-6 lazy parser to emit each block in
// the cheapest of {dynamic, static, uncompressed}, in
// byte-for-byte agreement with the reference.

const (
	numPrecodeSyms        = 19
	maxPrecodeCodewordLen = 7
)

// precodeLensPermutation matches deflate_precode_lens_permutation
// (deflate_compress.c:311). The dynamic block header emits the
// HCLEN+4 3-bit precode lengths in this special order so the most
// commonly nonzero entries appear first.
var precodeLensPermutation = [numPrecodeSyms]uint8{
	16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
}

// extraPrecodeBits matches deflate_extra_precode_bits
// (deflate_compress.c:316). Symbols 16/17/18 carry 2/3/7 extra bits
// describing run lengths; the rest are literal code-length values.
var extraPrecodeBits = [numPrecodeSyms]uint8{
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 7,
}

// dynamicCode bundles a fully-computed dynamic-Huffman code for one
// block. The litlen and offset code-lengths are stored contiguously
// in [lens] with offset entries starting at numLitlenSyms; this
// mirrors libdeflate's struct deflate_lens layout
// (deflate_compress.c:330) and lets the precode RLE pass treat them
// as a single concatenated array (deflate_compress.c:1592-1598).
type dynamicCode struct {
	litlenLens  [numLitlenSyms]uint8
	litlenCodes [numLitlenSyms]uint32
	offsetLens  [numOffsetSyms]uint8
	offsetCodes [numOffsetSyms]uint32

	numLitlenSyms   int                   // 257..286
	numOffsetSyms   int                   // 1..30
	numExplicitLens int                   // 4..19
	precodeLens     [numPrecodeSyms]uint8 // 0..7
	precodeCodes    [numPrecodeSyms]uint32
	precodeFreqs    [numPrecodeSyms]uint32
	precodeItems    []uint32 // RLE-encoded litlen+offset code lens
}

// buildDynamicCode constructs the dynamic litlen and offset Huffman
// codes from per-symbol frequencies, RLE-encodes the resulting code
// lengths using the precode alphabet, and builds the precode itself.
// The result is ready to be emitted by writeDynamicBlock.
//
// freqs.litlen[endOfBlock] is expected to have been incremented once
// by the caller (mirroring deflate_finish_block, deflate_compress.c:
// 2041), so the end-of-block symbol always has a codeword.
func buildDynamicCode(freqs *deflateFreqs) *dynamicCode {
	d := &dynamicCode{}

	makeHuffmanCode(numLitlenSyms, maxLitlenCodeLen,
		freqs.litlen[:], d.litlenLens[:], d.litlenCodes[:])
	makeHuffmanCode(numOffsetSyms, maxOffsetCodeLen,
		freqs.offset[:], d.offsetLens[:], d.offsetCodes[:])

	// Trim trailing zero entries to find the smallest legal HLIT /
	// HDIST values. The minima are 257 litlen syms and 1 offset sym
	// (deflate_compress.c:1574-1585).
	d.numLitlenSyms = numLitlenSyms
	for d.numLitlenSyms > 257 && d.litlenLens[d.numLitlenSyms-1] == 0 {
		d.numLitlenSyms--
	}
	d.numOffsetSyms = numOffsetSyms
	for d.numOffsetSyms > 1 && d.offsetLens[d.numOffsetSyms-1] == 0 {
		d.numOffsetSyms--
	}

	// Build a contiguous code-length array for the RLE pass: litlen
	// lengths up to numLitlenSyms followed immediately by offset
	// lengths (deflate_compress.c:1592-1598 does this in-place via
	// memmove).
	contig := make([]uint8, d.numLitlenSyms+d.numOffsetSyms)
	copy(contig[:d.numLitlenSyms], d.litlenLens[:d.numLitlenSyms])
	copy(contig[d.numLitlenSyms:], d.offsetLens[:d.numOffsetSyms])

	d.precodeItems = computePrecodeItems(contig, d.precodeFreqs[:])

	// Build the precode itself (deflate_compress.c:1612-1615).
	makeHuffmanCode(numPrecodeSyms, maxPrecodeCodewordLen,
		d.precodeFreqs[:], d.precodeLens[:], d.precodeCodes[:])

	// Count how many precode lengths we actually need to output
	// (deflate_compress.c:1618-1623). The minimum is 4 because the
	// HCLEN field stores num_explicit_lens - 4.
	d.numExplicitLens = numPrecodeSyms
	for d.numExplicitLens > 4 &&
		d.precodeLens[precodeLensPermutation[d.numExplicitLens-1]] == 0 {
		d.numExplicitLens--
	}

	return d
}

// computePrecodeItems mirrors deflate_compute_precode_items
// (deflate_compress.c:1482). It RLE-encodes the concatenated litlen
// +offset code-length array into precode items. Each returned item is
// packed: low 5 bits are the precode symbol, high bits are the extra
// bits for symbols 16/17/18.
func computePrecodeItems(lens []uint8, freqs []uint32) []uint32 {
	for i := range freqs[:numPrecodeSyms] {
		freqs[i] = 0
	}
	items := make([]uint32, 0, len(lens))
	runStart := 0
	for runStart != len(lens) {
		v := lens[runStart]
		runEnd := runStart + 1
		for runEnd != len(lens) && lens[runEnd] == v {
			runEnd++
		}

		if v == 0 {
			// Symbol 18: 11..138 zeroes.
			for (runEnd - runStart) >= 11 {
				extra := uint32((runEnd - runStart) - 11)
				if extra > 0x7F {
					extra = 0x7F
				}
				freqs[18]++
				items = append(items, 18|(extra<<5))
				runStart += 11 + int(extra)
			}
			// Symbol 17: 3..10 zeroes.
			if (runEnd - runStart) >= 3 {
				extra := uint32((runEnd - runStart) - 3)
				if extra > 0x7 {
					extra = 0x7
				}
				freqs[17]++
				items = append(items, 17|(extra<<5))
				runStart += 3 + int(extra)
			}
		} else {
			// Symbol 16: repeat previous length 3..6 more times,
			// only kicks in when we have >= 4 of the same length
			// (one literal + at least one 3-run).
			if (runEnd - runStart) >= 4 {
				freqs[v]++
				items = append(items, uint32(v))
				runStart++
				for (runEnd - runStart) >= 3 {
					extra := uint32((runEnd - runStart) - 3)
					if extra > 0x3 {
						extra = 0x3
					}
					freqs[16]++
					items = append(items, 16|(extra<<5))
					runStart += 3 + int(extra)
				}
			}
		}

		// Output any remaining lengths without RLE.
		for runStart != runEnd {
			freqs[v]++
			items = append(items, uint32(v))
			runStart++
		}
	}
	return items
}

// blockCosts holds the cost (in bits) of encoding a single block with
// each of the three block types, plus the dynamic code itself so it
// can be reused after the chooser commits.
type blockCosts struct {
	dynamic      uint64
	static       uint64
	uncompressed uint64
	dyn          *dynamicCode
}

// computeBlockCosts mirrors the cost-computation phase of
// deflate_flush_block (deflate_compress.c:1729-1801). All three costs
// include the 3-bit block header. The uncompressed cost also includes
// the padding to the next byte boundary plus the LEN/NLEN bytes for
// each 65535-byte chunk.
//
// bitcountBefore is the number of unflushed bits in the bitstream
// before the block starts; it only affects the uncompressed-cost
// padding term.
func computeBlockCosts(freqs *deflateFreqs, blockLength uint32, bitcountBefore uint8) blockCosts {
	d := buildDynamicCode(freqs)

	var c blockCosts
	c.dyn = d
	c.dynamic = 3
	c.static = 3
	c.uncompressed = 3

	// Cost of encoding the dynamic Huffman codes themselves
	// (deflate_compress.c:1751-1757). HLIT(5) + HDIST(5) + HCLEN(4)
	// + 3 bits per explicit precode length + each precode-item-emit.
	c.dynamic += 5 + 5 + 4 + 3*uint64(d.numExplicitLens)
	for sym := 0; sym < numPrecodeSyms; sym++ {
		extra := uint64(extraPrecodeBits[sym])
		c.dynamic += uint64(d.precodeFreqs[sym]) *
			(extra + uint64(d.precodeLens[sym]))
	}

	// Literal-symbol cost. Static lens are 8 (sym<144) or 9 (sym<256).
	for sym := 0; sym < 144; sym++ {
		c.dynamic += uint64(freqs.litlen[sym]) * uint64(d.litlenLens[sym])
		c.static += uint64(freqs.litlen[sym]) * 8
	}
	for sym := 144; sym < 256; sym++ {
		c.dynamic += uint64(freqs.litlen[sym]) * uint64(d.litlenLens[sym])
		c.static += uint64(freqs.litlen[sym]) * 9
	}

	// End-of-block symbol cost. Static EOB length is 7.
	c.dynamic += uint64(d.litlenLens[endOfBlock])
	c.static += 7

	// Length-symbol cost (litlen syms 257..285).
	for slot := 0; slot < len(extraLengthBits); slot++ {
		sym := firstLenSym + slot
		extra := uint64(extraLengthBits[slot])
		c.dynamic += uint64(freqs.litlen[sym]) *
			(extra + uint64(d.litlenLens[sym]))
		c.static += uint64(freqs.litlen[sym]) *
			(extra + uint64(staticLitlenLen[sym]))
	}

	// Offset-symbol cost. Static offset lens are all 5.
	for slot := 0; slot < len(extraOffsetBits); slot++ {
		extra := uint64(extraOffsetBits[slot])
		c.dynamic += uint64(freqs.offset[slot]) *
			(extra + uint64(d.offsetLens[slot]))
		c.static += uint64(freqs.offset[slot]) * (extra + 5)
	}

	// Uncompressed cost (deflate_compress.c:1797-1801): 3-bit header
	// (already counted) + padding to byte boundary + 32 bits of
	// LEN/NLEN + 40 extra bits per additional 65535-byte chunk +
	// 8 bits per source byte.
	padBits := uint64((-(int(bitcountBefore) + 3)) & 7)
	extraChunks := uint64(divRoundUp(blockLength, maxStoredBlockLen) - 1)
	c.uncompressed += padBits + 32 + 40*extraChunks + 8*uint64(blockLength)

	return c
}

// divRoundUp returns ceil(a/b) for unsigned a >= 0, b > 0.
func divRoundUp(a, b uint32) uint32 {
	return (a + b - 1) / b
}

// pickBlockType returns the type with the lowest cost. On ties the
// preference order is uncompressed > static > dynamic, matching
// deflate_compress.c:1804-1808 ("If there is a tie, prefer
// uncompressed, then static, then dynamic.").
func (c blockCosts) pickBlockType() int {
	best := c.dynamic
	t := blockTypeDynamic
	if c.static <= best {
		best = c.static
		t = blockTypeStatic
	}
	if c.uncompressed <= best {
		t = blockTypeStored
	}
	return t
}

// writeDynamicHeader emits the dynamic Huffman block header: BFINAL,
// BTYPE=10, HLIT, HDIST, HCLEN, the HCLEN+4 3-bit precode lengths in
// the special permutation, and the RLE-encoded litlen/offset code
// lengths (deflate_compress.c:1873-1925).
func writeDynamicHeader(bw *bitWriter, d *dynamicCode, last bool) {
	var bfinal uint64
	if last {
		bfinal = 1
	}
	bw.writeBits(bfinal, 1)
	bw.writeBits(uint64(blockTypeDynamic), 2)
	bw.writeBits(uint64(d.numLitlenSyms-257), 5)
	bw.writeBits(uint64(d.numOffsetSyms-1), 5)
	bw.writeBits(uint64(d.numExplicitLens-4), 4)

	for i := 0; i < d.numExplicitLens; i++ {
		sym := precodeLensPermutation[i]
		bw.writeBits(uint64(d.precodeLens[sym]), 3)
	}

	for _, it := range d.precodeItems {
		sym := it & 0x1F
		bw.writeBits(uint64(d.precodeCodes[sym]), d.precodeLens[sym])
		bw.writeBits(uint64(it>>5), extraPrecodeBits[sym])
	}
}

// writeDynamicBlock emits the full dynamic Huffman block: header,
// each literal/match item using the dynamic codes, then the
// end-of-block symbol.
func writeDynamicBlock(bw *bitWriter, items []item, d *dynamicCode, last bool) {
	writeDynamicHeader(bw, d, last)

	for _, it := range items {
		if it.isLiteral() {
			sym := uint32(it.literal)
			bw.writeBits(uint64(d.litlenCodes[sym]), d.litlenLens[sym])
			continue
		}
		writeDynamicMatch(bw, uint32(it.length), uint32(it.offset), d)
	}
	bw.writeBits(uint64(d.litlenCodes[endOfBlock]), d.litlenLens[endOfBlock])
}

// writeDynamicMatch emits a single match using the dynamic Huffman
// code: length symbol, length extra bits, offset symbol, offset
// extra bits.
func writeDynamicMatch(bw *bitWriter, length, offset uint32, d *dynamicCode) {
	lslot := lengthSlot[length]
	lsym := uint32(firstLenSym) + uint32(lslot)
	bw.writeBits(uint64(d.litlenCodes[lsym]), d.litlenLens[lsym])
	bw.writeBits(uint64(length-lengthSlotBase[lslot]), extraLengthBits[lslot])

	oslot := offsetSlot(offset)
	bw.writeBits(uint64(d.offsetCodes[oslot]), d.offsetLens[oslot])
	bw.writeBits(uint64(offset-offsetSlotBase[oslot]), extraOffsetBits[oslot])
}
