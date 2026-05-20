package cram

import (
	"fmt"
	"io"
)

// bitReader reads a CRAM CORE-data-block bitstream most-significant-bit
// first. CRAM packs several data series — Huffman codes, the unary
// prefixes of the Golomb/Elias families, and fixed-width BETA values —
// as a contiguous MSB-first bitstream inside a slice's CORE block. This
// reader exposes that bitstream both bit-at-a-time (for the unary parts
// of the integer codes) and N-bits-at-a-time (for the fixed-width
// parts).
//
// A bitReader is not safe for concurrent use.
type bitReader struct {
	data []byte // the CORE block's uncompressed payload.
	pos  int    // index of the next unread byte in data.
	cur  byte   // the byte currently being consumed bit by bit.
	nb   uint   // number of unconsumed low bits remaining in cur.
}

// newBitReader returns a bitReader over the given uncompressed CORE
// block payload.
func newBitReader(data []byte) *bitReader {
	return &bitReader{data: data}
}

// readBit reads and returns the next single bit as a 0 or 1. It returns
// io.ErrUnexpectedEOF when the bitstream is exhausted.
func (br *bitReader) readBit() (uint32, error) {
	if br.nb == 0 {
		if br.pos >= len(br.data) {
			return 0, io.ErrUnexpectedEOF
		}
		br.cur = br.data[br.pos]
		br.pos++
		br.nb = 8
	}
	br.nb--
	return uint32((br.cur >> br.nb) & 1), nil
}

// readBits reads n bits (0 <= n <= 32) most-significant-bit first and
// returns them as an unsigned integer. Reading zero bits yields zero
// and consumes nothing. It returns io.ErrUnexpectedEOF if the
// bitstream is exhausted before n bits are available.
func (br *bitReader) readBits(n uint) (uint32, error) {
	if n > 32 {
		return 0, fmt.Errorf("cram: bit reader asked for %d bits (max 32)", n)
	}
	var v uint32
	for i := uint(0); i < n; i++ {
		b, err := br.readBit()
		if err != nil {
			return 0, err
		}
		v = (v << 1) | b
	}
	return v, nil
}

// readBits64 reads n bits (0 <= n <= 64) most-significant-bit first as
// an unsigned 64-bit integer. It is used by the longer integer codes
// whose value can exceed 32 bits.
func (br *bitReader) readBits64(n uint) (uint64, error) {
	if n > 64 {
		return 0, fmt.Errorf("cram: bit reader asked for %d bits (max 64)", n)
	}
	var v uint64
	for i := uint(0); i < n; i++ {
		b, err := br.readBit()
		if err != nil {
			return 0, err
		}
		v = (v << 1) | uint64(b)
	}
	return v, nil
}

// readUnary counts and returns the number of consecutive 1-bits before
// the next 0-bit (the 0-bit is consumed). It is the prefix used by the
// Elias gamma, sub-exponential and Golomb code families. The limit
// bounds the count so a corrupt all-ones bitstream cannot loop forever.
func (br *bitReader) readUnary(limit int) (int, error) {
	n := 0
	for {
		b, err := br.readBit()
		if err != nil {
			return 0, err
		}
		if b == 0 {
			return n, nil
		}
		n++
		if n > limit {
			return 0, fmt.Errorf("cram: unary code exceeds %d bits (corrupt bitstream)", limit)
		}
	}
}

// align discards any unconsumed bits of the current byte so the reader
// is positioned at a byte boundary.
func (br *bitReader) align() {
	br.nb = 0
}

// exhausted reports whether every bit of the bitstream has been
// consumed.
func (br *bitReader) exhausted() bool {
	return br.pos >= len(br.data) && br.nb == 0
}
