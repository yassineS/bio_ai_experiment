package cram

import (
	"io"
)

// CRAM v4.0 replaces the ITF-8 / LTF-8 self-delimiting integers used by
// v2.x and v3.x with the htscodecs "uint7" variable-length encoding: a
// big-endian base-128 form in which each byte carries seven value bits in
// its low seven bits and a continuation flag (0x80) in its high bit. The
// bytes are most-significant first, so the decoder shifts the accumulator
// left by seven and ORs in each byte's low seven bits until a byte with
// the continuation bit clear is seen. A 32-bit value occupies at most five
// bytes; a 64-bit value at most ten.
//
// Signed values are folded into the unsigned form with zig-zag coding,
// exactly as htscodecs/varint.h's var_get_s32 / var_get_s64 do:
// (v >> 1) ^ -(v & 1). This places small-magnitude negatives near zero so
// they encode compactly, and lets the CRAM container's ref_seq_id (-1, -2)
// round-trip without the ITF-8 sign-extension trick v3 relied on.
//
// These mirror cram_io.c's uint7_get_32 / sint7_get_32 / uint7_get_64 /
// sint7_get_64 and varint.h's var_get_u32 / var_get_u64 (the BIG_END
// branch htscodecs compiles by default). The v2/v3 ITF-8 / LTF-8 path in
// itf8.go is left byte-for-byte unchanged; v4 dispatches here.

// uint7Max32 and uint7Max64 cap the byte count of a 32- and 64-bit uint7
// value. htslib stops the 32-bit decode after five bytes and the 64-bit
// decode after ten, ignoring any further continuation bits; matching that
// bound keeps a corrupt all-continuation stream from over-reading.
const (
	uint7Max32 = 5
	uint7Max64 = 10
)

// uint7At32 decodes one uint7-encoded value from p starting at off,
// truncated to 32 bits. It returns the value and the number of bytes
// consumed (1-5). It is the v4 counterpart of itf8At. It returns
// io.ErrUnexpectedEOF if p ends before the value is complete.
func uint7At32(p []byte, off int) (value int32, n int, err error) {
	v, n, err := uint7At64Limited(p, off, uint7Max32)
	if err != nil {
		return 0, 0, err
	}
	return int32(uint32(v)), n, nil
}

// sint7At32 decodes one signed (zig-zag) uint7 value from p starting at
// off, truncated to 32 bits. It returns the value and the number of bytes
// consumed (1-5). It is the v4 counterpart of itf8At for a signed field
// such as the container/slice ref_seq_id.
func sint7At32(p []byte, off int) (value int32, n int, err error) {
	u, n, err := uint7At64Limited(p, off, uint7Max32)
	if err != nil {
		return 0, 0, err
	}
	uu := uint32(u)
	return int32(uu>>1) ^ -int32(uu&1), n, nil
}

// uint7At64 decodes one uint7-encoded value from p starting at off as a
// 64-bit unsigned integer. It returns the value and the number of bytes
// consumed (1-10). It is the v4 counterpart of ltf8At.
func uint7At64(p []byte, off int) (value int64, n int, err error) {
	v, n, err := uint7At64Limited(p, off, uint7Max64)
	if err != nil {
		return 0, 0, err
	}
	return int64(v), n, nil
}

// sint7At64 decodes one signed (zig-zag) uint7 value from p starting at
// off as a 64-bit integer. It returns the value and the number of bytes
// consumed (1-10).
func sint7At64(p []byte, off int) (value int64, n int, err error) {
	u, n, err := uint7At64Limited(p, off, uint7Max64)
	if err != nil {
		return 0, 0, err
	}
	return int64(u>>1) ^ -int64(u&1), n, nil
}

// uint7At64Limited is the shared big-endian uint7 accumulator used by the
// byte-slice decoders. It reads at most maxBytes bytes (htslib stops the
// loop after the size-specific limit regardless of the continuation bit)
// and returns the accumulated unsigned value and the byte count consumed.
func uint7At64Limited(p []byte, off, maxBytes int) (uint64, int, error) {
	if off < 0 || off >= len(p) {
		return 0, 0, io.ErrUnexpectedEOF
	}
	var v uint64
	for i := 0; i < maxBytes; i++ {
		if off+i >= len(p) {
			return 0, 0, io.ErrUnexpectedEOF
		}
		c := p[off+i]
		v = (v << 7) | uint64(c&0x7f)
		if c&0x80 == 0 {
			return v, i + 1, nil
		}
	}
	// All maxBytes bytes had the continuation bit set; htslib stops here
	// and uses what it has, so the count is maxBytes.
	return v, maxBytes, nil
}

// readUint7_32 reads one uint7-encoded 32-bit value from r, returning the
// value and the byte count consumed (1-5). It is the streaming counterpart
// of uint7At32 used by the v4 container- and block-header parsers, which
// read directly from an io.Reader rather than a buffered payload.
func readUint7_32(r io.Reader) (value int32, n int, err error) {
	v, n, err := readUint7(r, uint7Max32)
	if err != nil {
		return 0, 0, err
	}
	return int32(uint32(v)), n, nil
}

// readSint7_32 reads one signed (zig-zag) uint7 32-bit value from r.
func readSint7_32(r io.Reader) (value int32, n int, err error) {
	v, n, err := readUint7(r, uint7Max32)
	if err != nil {
		return 0, 0, err
	}
	u := uint32(v)
	return int32(u>>1) ^ -int32(u&1), n, nil
}

// readUint7_64 reads one uint7-encoded 64-bit value from r, returning the
// value and the byte count consumed (1-10).
func readUint7_64(r io.Reader) (value int64, n int, err error) {
	v, n, err := readUint7(r, uint7Max64)
	if err != nil {
		return 0, 0, err
	}
	return int64(v), n, nil
}

// readUint7 is the shared streaming big-endian uint7 accumulator. It reads
// one byte at a time, stopping at the first byte without the continuation
// flag or after maxBytes bytes, whichever comes first.
func readUint7(r io.Reader, maxBytes int) (uint64, int, error) {
	var b [1]byte
	var v uint64
	for i := 0; i < maxBytes; i++ {
		if _, err := io.ReadFull(r, b[:]); err != nil {
			return 0, 0, eofToUnexpected(err, "uint7 varint")
		}
		c := b[0]
		v = (v << 7) | uint64(c&0x7f)
		if c&0x80 == 0 {
			return v, i + 1, nil
		}
	}
	return v, maxBytes, nil
}

// appendUint7 encodes v as a big-endian uint7 byte sequence and appends it
// to dst, returning the extended slice. It is the inverse of uint7At64 and
// exists so tests (and any future encode path) can round-trip values; it
// mirrors varint.h's var_put_u64 BIG_END branch.
func appendUint7(dst []byte, v uint64) []byte {
	// Count the 7-bit groups needed, most-significant first.
	shift := 0
	for x := v; ; x >>= 7 {
		shift += 7
		if x < 0x80 {
			break
		}
	}
	for shift > 7 {
		shift -= 7
		dst = append(dst, byte((v>>uint(shift))&0x7f)|0x80)
	}
	return append(dst, byte(v&0x7f))
}

// appendSint7 encodes a signed value as a zig-zag uint7 sequence.
func appendSint7(dst []byte, v int64) []byte {
	return appendUint7(dst, uint64(v<<1)^uint64(v>>63))
}
