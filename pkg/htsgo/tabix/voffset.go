package tabix

import bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"

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

// VOffsetAt converts an absolute uncompressed byte position into a BGZF
// virtual offset, given the BGZF block table for the file. It mirrors what
// htslib's bgzf_tell would return for that read position.
//
// The crucial detail is the block-boundary normalization: a position that
// lands exactly on the end of a block's payload (== the start of the next
// block) is represented as (next_block.CompressedOffset, 0) rather than
// (this_block.CompressedOffset, this_block.UncompressedSize). htslib advances
// to the next block as soon as the current block is fully consumed, so the
// virtual offset of the byte just past the last record in a block rolls over
// to the next block's start. Getting this wrong makes the meta/pseudo-bin and
// data-bin "end" offsets for the final record of each block disagree with
// upstream by exactly one block-rollover.
func VOffsetAt(offsets []bgzip.BlockOffset, pos int64) VOffset {
	if len(offsets) == 0 {
		return MakeVOffset(0, int(pos))
	}
	// Binary search for the last block whose UncompressedOffset <= pos.
	lo, hi := 0, len(offsets)
	for lo < hi {
		mid := (lo + hi) / 2
		if offsets[mid].UncompressedOffset <= pos {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	i := lo - 1
	if i < 0 {
		i = 0
	}
	blk := offsets[i]
	uoff := int(pos - blk.UncompressedOffset)
	// If pos is exactly at this block's payload end, roll over to the start
	// of the next block (uoff 0), matching htslib's bgzf_tell. The next
	// block's compressed offset is the current block's compressed offset
	// plus its compressed size — which is also correct for the final data
	// block, whose "next" is the trailing EOF marker block (Scan does not
	// surface that block, but CompressedSize lets us compute its start).
	if uoff == blk.UncompressedSize {
		next := blk.CompressedOffset + int64(blk.CompressedSize)
		return MakeVOffset(next, 0)
	}
	return MakeVOffset(blk.CompressedOffset, uoff)
}

// Coff returns the compressed-file offset bits.
func (v VOffset) Coff() int64 { return int64(uint64(v) >> 16) }

// Uoff returns the in-block uncompressed offset bits.
func (v VOffset) Uoff() int { return int(uint64(v) & 0xFFFF) }
