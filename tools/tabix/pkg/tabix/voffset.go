package tabix

// VOffset is a BGZF "virtual offset": the 48 high bits are the compressed
// byte offset of the start of a BGZF block, and the 16 low bits are the byte
// offset of a record's first byte within that block's decompressed payload.
//
// Two virtual offsets compare as plain uint64s — a useful property that the
// linear/bin indices rely on for "is this record before that one" tests.
type VOffset uint64

// MakeVOffset packs a (compressed offset, uncompressed offset) pair into a
// VOffset. uoff must fit in 16 bits, which is guaranteed for BGZF since the
// uncompressed block size is capped at 64 KiB (the tabix MaxBlockSize is
// 65,280 bytes — still ≤ 65,535).
func MakeVOffset(coff int64, uoff int) VOffset {
	if coff < 0 {
		coff = 0
	}
	if uoff < 0 {
		uoff = 0
	}
	if uoff > 0xFFFF {
		// Truncate; callers should not pass uoff > MaxBlockSize. The cap
		// keeps the resulting offset well-formed instead of corrupting the
		// compressed-offset bits.
		uoff = 0xFFFF
	}
	return VOffset(uint64(coff)<<16 | uint64(uoff)&0xFFFF)
}

// Coff returns the compressed-file offset bits.
func (v VOffset) Coff() int64 { return int64(uint64(v) >> 16) }

// Uoff returns the in-block uncompressed offset bits.
func (v VOffset) Uoff() int { return int(uint64(v) & 0xFFFF) }
