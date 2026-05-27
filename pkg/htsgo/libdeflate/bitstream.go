package libdeflate

// bitWriter is an LSB-first bitstream writer that mirrors libdeflate's
// deflate_output_bitstream. Bits are accumulated in a uint64 word and
// flushed to the output buffer one byte at a time once eight or more
// bits are pending. This matches the byte-for-byte output ordering that
// the reference encoder produces.
type bitWriter struct {
	out      []byte
	bitbuf   uint64
	bitcount uint8
}

// newBitWriter returns a writer that appends to dst.
func newBitWriter(dst []byte) *bitWriter {
	return &bitWriter{out: dst}
}

// writeBits appends the low n bits of value to the bitstream, LSB first.
// n must be in the range [0, 56] so that the accumulator never overflows
// between flushes (the same bound libdeflate uses on 64-bit targets).
func (b *bitWriter) writeBits(value uint64, n uint8) {
	if n == 0 {
		return
	}
	mask := (uint64(1) << n) - 1
	b.bitbuf |= (value & mask) << b.bitcount
	b.bitcount += n
	for b.bitcount >= 8 {
		b.out = append(b.out, byte(b.bitbuf))
		b.bitbuf >>= 8
		b.bitcount -= 8
	}
}

// alignToByte pads the current bit position with zero bits up to the
// next byte boundary. After this call, bitcount is zero.
func (b *bitWriter) alignToByte() {
	if b.bitcount != 0 {
		b.out = append(b.out, byte(b.bitbuf))
		b.bitbuf = 0
		b.bitcount = 0
	}
}

// writeBytes appends raw bytes to the underlying buffer. The bitstream
// must be byte-aligned before this is called.
func (b *bitWriter) writeBytes(p []byte) {
	if b.bitcount != 0 {
		panic("libdeflate: writeBytes called on unaligned bitstream")
	}
	b.out = append(b.out, p...)
}

// flush writes out any trailing partial byte, padding with zero bits,
// and returns the accumulated output buffer.
func (b *bitWriter) flush() []byte {
	b.alignToByte()
	return b.out
}
