package libdeflate

import "fmt"

// minCompressionLevel and maxCompressionLevel mirror libdeflate's
// supported range. Slice 1 accepts any level in [1, 12] but does not
// honor the per-level matchfinder/cost tradeoffs beyond the
// passthrough-size threshold; subsequent slices will fill that in.
const (
	minCompressionLevel = 1
	maxCompressionLevel = 12
)

// GzipCompress compresses src at the given compression level and
// returns the gzip-wrapped DEFLATE byte stream. The output is
// byte-identical to libdeflate's libdeflate_gzip_compress on the
// fixtures in the oracle corpus (empty, single_byte, repeated_a,
// random_64k, bgzf_payload) at level 6. The Slice 3 path:
//
//  1. For inputs <= the per-level passthrough threshold, emit a
//     single all-STORED stream (mirrors deflate_compress_none in
//     reference_code/libdeflate/lib/deflate_compress.c).
//  2. Otherwise run the lazy matchfinder (Slice 2) to produce one or
//     more blocks of literal/match items, then for each block run the
//     cost-based chooser to pick dynamic / static / uncompressed and
//     emit accordingly.
//
// Levels other than 5-7 fall back to the Slice 1 trivial encoder for
// now; expanding the level coverage is a follow-up slice.
func GzipCompress(src []byte, level int) ([]byte, error) {
	if level < minCompressionLevel || level > maxCompressionLevel {
		return nil, fmt.Errorf("libdeflate: invalid compression level %d (want %d..%d)",
			level, minCompressionLevel, maxCompressionLevel)
	}

	out := make([]byte, 0, gzipMinOverhead+len(src))
	out = writeGzipHeader(out, level)

	bw := newBitWriter(out)
	writeDeflateStream(bw, src, level)
	out = bw.flush()

	out = writeGzipTrailer(out, src)
	return out, nil
}

// DeflateCompress compresses src at the given compression level and
// returns the raw DEFLATE byte stream (RFC 1951) — no gzip header or
// trailer. The output is byte-identical to libdeflate's
// libdeflate_deflate_compress on the same fixtures GzipCompress
// validates against. This is the entry point used by the BGZF writer
// (pkg/htsgo/bgzf), since BGZF wraps each block with its own
// gzip-with-BC-subfield header instead of the standard RFC 1952 one.
func DeflateCompress(src []byte, level int) ([]byte, error) {
	if level < minCompressionLevel || level > maxCompressionLevel {
		return nil, fmt.Errorf("libdeflate: invalid compression level %d (want %d..%d)",
			level, minCompressionLevel, maxCompressionLevel)
	}
	bw := newBitWriter(make([]byte, 0, len(src)/2+16))
	writeDeflateStream(bw, src, level)
	return bw.flush(), nil
}

// writeDeflateStream emits one or more DEFLATE blocks for src into bw,
// terminating with BFINAL=1 on the final block. The dispatch matches
// libdeflate's deflate_compress: short inputs go through the all-stored
// passthrough, mid-range levels (5–7) drive the lazy parser + cost
// chooser, and other levels fall back to the trivial encoder (Slice 1).
func writeDeflateStream(bw *bitWriter, src []byte, level int) {
	switch {
	case uint64(len(src)) <= maxPassthroughSize(level):
		writeStoredBlocks(bw, src, true)
	case level >= 2 && level <= 7:
		writeLazyBlocks(bw, src, level)
	default:
		items := trivialLZ77(src)
		writeStaticBlock(bw, items, true)
	}
}

// writeLazyBlocks runs the lazy matchfinder over src and emits each
// block using the cost-cheapest of {dynamic, static, uncompressed}.
// Mirrors the per-block dispatch in deflate_compress_lazy_generic ->
// deflate_finish_block -> deflate_flush_block.
func writeLazyBlocks(bw *bitWriter, src []byte, level int) {
	blocks := lazyEmitBlocks(src, level)
	for i := range blocks {
		blk := &blocks[i]
		last := i == len(blocks)-1
		// deflate_finish_block (deflate_compress.c:2041) increments
		// the EOB frequency before building the dynamic code, since
		// the EOB symbol is always emitted exactly once per block.
		blk.freqs.litlen[endOfBlock]++
		costs := computeBlockCosts(&blk.freqs, blk.length, bw.bitcount)
		blockBegin := int(blk.begin)
		blockData := src[blockBegin : blockBegin+int(blk.length)]
		switch costs.pickBlockType() {
		case blockTypeStored:
			writeStoredBlocks(bw, blockData, last)
		case blockTypeStatic:
			writeStaticBlockTail(bw, blk.items, last)
		default: // blockTypeDynamic
			writeDynamicBlock(bw, blk.items, costs.dyn, last)
		}
	}
}

// writeStaticBlockTail is a thin alias for writeStaticBlock; it exists
// purely so the dispatch site reads symmetrically with the dynamic /
// stored cases.
func writeStaticBlockTail(bw *bitWriter, items []item, last bool) {
	writeStaticBlock(bw, items, last)
}

// maxPassthroughSize matches libdeflate's per-level threshold below
// which the encoder switches to the all-STORED "passthrough" path.
// See libdeflate's deflate_compress.c:
//
//	c->max_passthrough_size = 55 - (compression_level * 4);
//	if (c->compression_level == 0) c->max_passthrough_size = SIZE_MAX;
//
// We never expose level 0 (the public API rejects it), so the
// SIZE_MAX clamp is not needed here.
func maxPassthroughSize(level int) uint64 {
	v := 55 - 4*level
	if v < 0 {
		return 0
	}
	return uint64(v)
}

// trivialLZ77 produces a deterministic literal/match item stream for
// the input. It mirrors libdeflate's hash-chain matchfinder closely
// enough to reproduce the level-6 match decisions on the Slice 1
// corpus inputs (only repeated_a exercises the matching path), but
// without the bells and whistles of the lazy heuristic, depth caps,
// nice-length cutoffs, or the precomputed `next_hashes` array. The
// full lazy matchfinder lands in Slice 2.
//
// Notably, libdeflate's matcher only registers a position in the hash
// table at the moment that position is examined — it does NOT
// pre-populate the table. Because the lookup happens before the
// insert, the very first occurrence of any 3-byte key always yields a
// miss. This matters for repeated_a: at pos=1 the hash bucket for
// "AAA" is still empty (pos=0 was indexed under the initial-hash
// bucket 0, not the real "AAA" bucket, because of the off-by-one in
// libdeflate's `next_hashes` priming), so libdeflate emits a literal
// rather than a length-3 match at offset 1. We reproduce that quirk
// here by skipping the matcher entirely at positions 0 and 1.
func trivialLZ77(src []byte) []item {
	items := make([]item, 0, len(src))
	const hashBits = 15
	const hashSize = 1 << hashBits
	// head[h] is the most recent position whose 3-byte hash equals h,
	// or -1 if there is none. prev[p] is the previous such position
	// for the chain through p.
	head := make([]int32, hashSize)
	for i := range head {
		head[i] = -1
	}
	prev := make([]int32, len(src))
	for i := range prev {
		prev[i] = -1
	}

	insert := func(pos int) {
		if pos+minMatchLen > len(src) {
			return
		}
		h := hash3(src, pos, hashBits)
		prev[pos] = head[h]
		head[h] = int32(pos)
	}

	i := 0
	for i < len(src) {
		// libdeflate's `next_hashes` initialization leaves the hash
		// for position 0 mis-bucketed (under bucket 0 instead of the
		// real 3-byte hash), and at position 1 the real bucket is
		// still empty. The net effect is that the first two
		// positions never produce matches. Replicate that here
		// rather than trying to model the priming exactly.
		if i < 2 {
			items = append(items, litItem(src[i]))
			insert(i)
			i++
			continue
		}
		bestLen, bestOff := findHashMatch(src, i, head, prev)
		insert(i)
		if bestLen >= minMatchLen {
			items = append(items, matchItem(uint16(bestLen), uint16(bestOff)))
			// Insert hashes for the positions skipped over by the
			// match so future positions can match against them.
			for k := 1; k < bestLen; k++ {
				insert(i + k)
			}
			i += bestLen
			continue
		}
		items = append(items, litItem(src[i]))
		i++
	}
	return items
}

// hash3 returns the lz_hash of the 3-byte sequence at src[pos] folded
// into hashBits bits. The constant matches libdeflate's lz_hash for a
// 24-bit input.
func hash3(src []byte, pos int, hashBits uint) uint32 {
	const golden32 = 0x1E35A7BD
	v := uint32(src[pos]) | uint32(src[pos+1])<<8 | uint32(src[pos+2])<<16
	return (v * golden32) >> (32 - hashBits)
}

// findHashMatch walks the hash chain for the 3-byte key at src[pos]
// and returns the longest match found. The chain is sorted by
// decreasing position, so the closest (smallest-offset) match is
// visited first.
func findHashMatch(src []byte, pos int, head, prev []int32) (length, offset int) {
	if pos+minMatchLen > len(src) {
		return 0, 0
	}
	maxLen := len(src) - pos
	if maxLen > maxMatchLen {
		maxLen = maxMatchLen
	}
	earliest := pos - maxMatchOffset
	if earliest < 0 {
		earliest = 0
	}
	h := hash3(src, pos, 15)
	bestLen := 0
	bestOff := 0
	for cur := int(head[h]); cur >= earliest; cur = int(prev[cur]) {
		if src[cur] != src[pos] {
			continue
		}
		n := 0
		for n < maxLen && src[cur+n] == src[pos+n] {
			n++
		}
		if n > bestLen {
			bestLen = n
			bestOff = pos - cur
			if bestLen >= maxLen {
				break
			}
		}
	}
	if bestLen < minMatchLen {
		return 0, 0
	}
	return bestLen, bestOff
}
