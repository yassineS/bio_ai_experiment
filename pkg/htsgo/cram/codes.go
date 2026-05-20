package cram

import (
	"fmt"
)

// unaryLimit bounds the unary prefix of any CRAM integer code. A
// well-formed CRAM file never approaches it; the cap exists so a corrupt
// all-ones bitstream errors out instead of looping forever.
const unaryLimit = 64

// readGamma decodes one Elias gamma code from the CORE bitstream. The
// CRAM gamma code is: a unary count of leading 1-bits giving the length
// L, then L further bits whose value, with an implicit leading 1
// prepended, is the result. A code of "0" (zero leading 1-bits, no
// trailing bits) decodes to 1.
func readGamma(br *bitReader) (int32, error) {
	length, err := br.readUnary(unaryLimit)
	if err != nil {
		return 0, fmt.Errorf("cram: gamma code prefix: %w", err)
	}
	if length > 31 {
		return 0, fmt.Errorf("cram: gamma code length %d exceeds 31 bits", length)
	}
	tail, err := br.readBits(uint(length))
	if err != nil {
		return 0, fmt.Errorf("cram: gamma code body: %w", err)
	}
	// The decoded value is the trailing bits with an implicit leading 1.
	return int32((uint32(1) << uint(length)) | tail), nil
}

// readSubexp decodes one sub-exponential code with order parameter k
// from the CORE bitstream. The sub-exponential code, per the CRAM spec,
// reads a unary count n of leading 1-bits; when n is zero the value is
// the next k bits, otherwise the value is formed from a (k+n-1)-bit
// tail with the implicit 1<<(k+n-1) bucket base added.
func readSubexp(br *bitReader, k int) (int32, error) {
	if k < 0 || k > 31 {
		return 0, fmt.Errorf("cram: subexp order k=%d out of range", k)
	}
	n, err := br.readUnary(unaryLimit)
	if err != nil {
		return 0, fmt.Errorf("cram: subexp code prefix: %w", err)
	}
	var nbits int
	var base uint32
	if n == 0 {
		nbits = k
		base = 0
	} else {
		nbits = k + n - 1
		base = uint32(1) << uint(nbits)
	}
	if nbits > 31 {
		return 0, fmt.Errorf("cram: subexp code tail %d bits exceeds 31", nbits)
	}
	tail, err := br.readBits(uint(nbits))
	if err != nil {
		return 0, fmt.Errorf("cram: subexp code body: %w", err)
	}
	return int32(base + tail), nil
}

// readGolomb decodes one Golomb code with divisor m from the CORE
// bitstream. The code is a unary quotient q (q leading 1-bits then a
// 0-bit) followed by a truncated-binary remainder r in [0, m). The
// decoded value is q*m + r. Golomb-Rice is the special case where m is
// a power of two; this routine handles both.
func readGolomb(br *bitReader, m int) (int32, error) {
	if m <= 0 {
		return 0, fmt.Errorf("cram: golomb divisor m=%d must be positive", m)
	}
	q, err := br.readUnary(unaryLimit << 4)
	if err != nil {
		return 0, fmt.Errorf("cram: golomb quotient: %w", err)
	}
	if m == 1 {
		// Divisor 1: the remainder is empty, the value is the quotient.
		return int32(q), nil
	}
	// Truncated binary coding of the remainder in [0, m). Let b be the
	// number of bits needed for m; the first 2^b - m remainders use
	// b-1 bits, the rest use b bits.
	b := 0
	for (1 << uint(b)) < m {
		b++
	}
	cutoff := (1 << uint(b)) - m
	r, err := br.readBits(uint(b - 1))
	if err != nil {
		return 0, fmt.Errorf("cram: golomb remainder prefix: %w", err)
	}
	rem := int(r)
	if rem >= cutoff {
		extra, err := br.readBit()
		if err != nil {
			return 0, fmt.Errorf("cram: golomb remainder suffix: %w", err)
		}
		rem = ((rem << 1) | int(extra)) - cutoff
	}
	return int32(q*m + rem), nil
}
