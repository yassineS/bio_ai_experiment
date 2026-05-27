package libdeflate

// item represents a single literal or back-reference produced by the
// LZ77 stage. When length is zero the item is a literal byte; otherwise
// it is a (length, offset) match with length in [minMatchLen, maxMatchLen]
// and offset in [1, maxMatchOffset].
type item struct {
	literal byte
	length  uint16
	offset  uint16
}

func litItem(b byte) item        { return item{literal: b} }
func matchItem(l, o uint16) item { return item{length: l, offset: o} }

func (it item) isLiteral() bool { return it.length == 0 }

// writeStaticBlock emits a final-or-not STATIC Huffman block: the 3-bit
// header (BFINAL | BTYPE=01), the static-coded items, and the
// end-of-block symbol.
func writeStaticBlock(bw *bitWriter, items []item, last bool) {
	var bfinal uint64
	if last {
		bfinal = 1
	}
	// 3-bit header: bit 0 = BFINAL, bits 1-2 = BTYPE (=01 for static).
	bw.writeBits(bfinal|(uint64(blockTypeStatic)<<1), 3)

	for _, it := range items {
		if it.isLiteral() {
			sym := uint32(it.literal)
			bw.writeBits(uint64(staticLitlenCode[sym]),
				staticLitlenLen[sym])
			continue
		}
		writeStaticMatch(bw, uint32(it.length), uint32(it.offset))
	}

	// End-of-block symbol.
	bw.writeBits(uint64(staticLitlenCode[endOfBlock]),
		staticLitlenLen[endOfBlock])
}

// writeStaticMatch emits a single match using the static Huffman code:
// the length symbol, its extra bits, the offset symbol, and its extra
// bits.
func writeStaticMatch(bw *bitWriter, length, offset uint32) {
	lslot := lengthSlot[length]
	lsym := uint32(firstLenSym) + uint32(lslot)
	bw.writeBits(uint64(staticLitlenCode[lsym]),
		staticLitlenLen[lsym])
	bw.writeBits(uint64(length-lengthSlotBase[lslot]),
		extraLengthBits[lslot])

	oslot := offsetSlot(offset)
	bw.writeBits(uint64(staticOffsetCode[oslot]),
		staticOffsetLen[oslot])
	bw.writeBits(uint64(offset-offsetSlotBase[oslot]),
		extraOffsetBits[oslot])
}
