package libdeflate

// Static Huffman code tables for DEFLATE (RFC 1951 §3.2.6). The
// literal/length code uses codeword lengths 8, 9, 7, 8 over the four
// symbol ranges; the offset code uses fixed 5-bit codewords. Codewords
// are bit-reversed at init time so they can be pushed straight into the
// LSB-first bitstream by the encoders.

var (
	staticLitlenCode [numLitlenSyms]uint32
	staticLitlenLen  [numLitlenSyms]uint8
	staticOffsetCode [numOffsetSyms]uint32
	staticOffsetLen  [numOffsetSyms]uint8
)

func init() {
	// Literal/length code: assign codeword lengths per RFC 1951 §3.2.6
	// and then fill in canonical codewords.
	for sym := 0; sym < 144; sym++ {
		staticLitlenLen[sym] = 8
	}
	for sym := 144; sym < 256; sym++ {
		staticLitlenLen[sym] = 9
	}
	for sym := 256; sym < 280; sym++ {
		staticLitlenLen[sym] = 7
	}
	for sym := 280; sym < 288; sym++ {
		staticLitlenLen[sym] = 8
	}
	assignCanonicalCodes(staticLitlenLen[:], staticLitlenCode[:])

	// Offset code: 30 valid 5-bit codewords (symbols 30 and 31 are
	// reserved and unused in static blocks; we still assign them a code
	// so the table is well-defined).
	for sym := 0; sym < numOffsetSyms; sym++ {
		staticOffsetLen[sym] = 5
	}
	assignCanonicalCodes(staticOffsetLen[:], staticOffsetCode[:])
}

// assignCanonicalCodes fills codes[sym] with the bit-reversed canonical
// codeword for each symbol, given the per-symbol code lengths lens.
// Codes are reversed so that writeBits can push them LSB-first.
func assignCanonicalCodes(lens []uint8, codes []uint32) {
	const maxLen = maxCodewordLen
	var lenCounts [maxLen + 1]uint32
	for _, l := range lens {
		lenCounts[l]++
	}
	var nextCode [maxLen + 1]uint32
	var code uint32
	lenCounts[0] = 0
	for bits := 1; bits <= maxLen; bits++ {
		code = (code + lenCounts[bits-1]) << 1
		nextCode[bits] = code
	}
	for sym, l := range lens {
		if l == 0 {
			continue
		}
		c := nextCode[l]
		nextCode[l]++
		codes[sym] = reverseBits(c, l)
	}
}

// reverseBits returns the low n bits of v in reverse order.
func reverseBits(v uint32, n uint8) uint32 {
	var r uint32
	for i := uint8(0); i < n; i++ {
		r = (r << 1) | (v & 1)
		v >>= 1
	}
	return r
}
