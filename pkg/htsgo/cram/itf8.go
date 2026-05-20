package cram

import (
	"fmt"
	"io"
)

// CRAM's two self-delimiting integer encodings. In both, the number of
// leading 1-bits in the first byte gives how many *extra* bytes follow;
// the remaining low bits of the first byte are the most-significant bits
// of the value, and each extra byte contributes 8 more bits (big-endian
// order). ITF-8 encodes a 32-bit value in 1-5 bytes; LTF-8 encodes a
// 64-bit value in 1-9 bytes.

// readITF8 reads one ITF-8-encoded 32-bit integer from r. The returned
// count is the number of bytes consumed (1-5). The value is interpreted
// as a signed 32-bit integer, matching htslib's int32_t usage.
func readITF8(r io.Reader) (value int32, n int, err error) {
	var b0 [1]byte
	if _, err = io.ReadFull(r, b0[:]); err != nil {
		return 0, 0, eofToUnexpected(err, "ITF-8 first byte")
	}
	first := b0[0]
	var extra int
	var v uint32
	switch {
	case first&0x80 == 0: // 0xxxxxxx — 0 extra bytes, 7 value bits.
		extra = 0
		v = uint32(first & 0x7f)
	case first&0x40 == 0: // 10xxxxxx — 1 extra byte, 6+8 value bits.
		extra = 1
		v = uint32(first & 0x3f)
	case first&0x20 == 0: // 110xxxxx — 2 extra bytes, 5+16 value bits.
		extra = 2
		v = uint32(first & 0x1f)
	case first&0x10 == 0: // 1110xxxx — 3 extra bytes, 4+24 value bits.
		extra = 3
		v = uint32(first & 0x0f)
	default: // 1111xxxx — 4 extra bytes; only the low 4 bits of the
		// last byte are used, the first byte's low 4 bits are the top.
		extra = 4
		v = uint32(first & 0x0f)
	}
	buf := make([]byte, extra)
	if extra > 0 {
		if _, err = io.ReadFull(r, buf); err != nil {
			return 0, 0, eofToUnexpected(err, "ITF-8 continuation bytes")
		}
	}
	for i := 0; i < extra; i++ {
		if extra == 4 && i == extra-1 {
			// The fifth (last) byte of a 5-byte ITF-8 contributes only
			// its low 4 bits: 4+8+8+8+4 = 32 value bits in total.
			v = (v << 4) | uint32(buf[i]&0x0f)
		} else {
			v = (v << 8) | uint32(buf[i])
		}
	}
	return int32(v), extra + 1, nil
}

// readLTF8 reads one LTF-8-encoded 64-bit integer from r. The returned
// count is the number of bytes consumed (1-9). The value is interpreted
// as a signed 64-bit integer, matching htslib's int64_t usage.
func readLTF8(r io.Reader) (value int64, n int, err error) {
	var b0 [1]byte
	if _, err = io.ReadFull(r, b0[:]); err != nil {
		return 0, 0, eofToUnexpected(err, "LTF-8 first byte")
	}
	first := b0[0]
	var extra int
	var v uint64
	switch {
	case first&0x80 == 0: // 0xxxxxxx
		extra, v = 0, uint64(first&0x7f)
	case first&0x40 == 0: // 10xxxxxx
		extra, v = 1, uint64(first&0x3f)
	case first&0x20 == 0: // 110xxxxx
		extra, v = 2, uint64(first&0x1f)
	case first&0x10 == 0: // 1110xxxx
		extra, v = 3, uint64(first&0x0f)
	case first&0x08 == 0: // 11110xxx
		extra, v = 4, uint64(first&0x07)
	case first&0x04 == 0: // 111110xx
		extra, v = 5, uint64(first&0x03)
	case first&0x02 == 0: // 1111110x
		extra, v = 6, uint64(first&0x01)
	case first&0x01 == 0: // 11111110
		extra, v = 7, 0
	default: // 11111111 — 8 extra bytes, the full 64-bit value follows.
		extra, v = 8, 0
	}
	buf := make([]byte, extra)
	if extra > 0 {
		if _, err = io.ReadFull(r, buf); err != nil {
			return 0, 0, eofToUnexpected(err, "LTF-8 continuation bytes")
		}
	}
	for i := 0; i < extra; i++ {
		v = (v << 8) | uint64(buf[i])
	}
	return int64(v), extra + 1, nil
}

// eofToUnexpected converts a plain io.EOF returned mid-structure into an
// io.ErrUnexpectedEOF-flavoured error, so truncated input is reported as
// a parse failure rather than a clean end-of-stream.
func eofToUnexpected(err error, what string) error {
	if err == io.EOF {
		return fmt.Errorf("cram: truncated input reading %s: %w", what, io.ErrUnexpectedEOF)
	}
	return fmt.Errorf("cram: reading %s: %w", what, err)
}
