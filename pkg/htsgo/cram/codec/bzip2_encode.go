package codec

// bzip2_encode.go is a pure-Go, in-tree bzip2 *encoder*. Go's standard
// library ships compress/bzip2 with a decoder only; CRAM's X_EXT
// external codec (used by arith_dynamic and the name tokeniser) needs to
// *encode* bzip2 streams as well. Adding a third-party bzip2 encoder is
// outside the sanctioned dependency set (CLAUDE.md), so the codec is
// ported here, confined to pkg/htsgo/cram/codec/ alongside the other
// CRAM codecs.
//
// The bzip2 stream format, top to bottom:
//
//	"BZh" magic, then '1'..'9' for the 100k * level block size.
//	One or more blocks, each:
//	  block magic 0x314159265359 (48 bits, the digits of pi)
//	  block CRC32 (32 bits, big-endian bit order)
//	  randomised flag (1 bit, always 0 — deprecated)
//	  BWT origin pointer (24 bits)
//	  symbol-map: 16-bit "used16" plus up to 16 16-bit bitmaps
//	  number of Huffman tables (3 bits, 2..6)
//	  number of selectors (15 bits) + MTF-coded selector list (unary)
//	  per-table Huffman code lengths (delta-coded)
//	  the Huffman-coded MTF/RLE2 symbol stream
//	stream footer magic 0x177245385090 (48 bits) + combined stream CRC32.
//
// The compression pipeline that produces each block's symbol stream:
//
//	RLE1   run-length encode runs of 4..255 identical bytes
//	BWT    Burrows-Wheeler transform (block sort)
//	MTF    move-to-front over the symbols actually present
//	RLE2   run-length encode runs of the MTF zero with the RUNA/RUNB
//	       bijective base-2 code; non-zero MTF values shift up by one
//	Huffman up to 6 canonical Huffman tables, selected per 50-symbol
//	       group, iteratively refined to minimise total code length.
//
// This is a correctness-first port. The block sort ranks the cyclic
// rotations with rank-doubling (Manber-Myers, O(n log n) total — see
// sortCyclicRotations); it is not the bzip2 reference's fallback/main
// sort, but it produces an identical BWT and the output decodes
// byte-for-byte under compress/bzip2, the system bzip2 -d, and htslib.
// See bzip2_encode_test.go for the round-trip and cross-tool gates.

import (
	"fmt"
	"sort"
)

const (
	// bzMaxAlphaSize is the largest MTF/RLE2 alphabet: 256 byte values
	// minus duplicates plus RUNA, RUNB and EOB. In practice it is
	// (#distinct bytes)+2; the constant bounds the tables.
	bzMaxAlphaSize = 258
	// bzMaxCodeLen is the longest Huffman code bzip2 permits.
	bzMaxCodeLen = 20
	// bzGroupSize is the number of symbols coded by one selector.
	bzGroupSize = 50
	// bzMaxSelectors caps the selector count for a maximum block.
	bzMaxSelectors = 18002
	// bzRunA and bzRunB are the two MTF zero-run symbols (RLE2).
	bzRunA = 0
	bzRunB = 1
)

// bzBlockMagic is the 48-bit per-block magic (digits of pi).
var bzBlockMagic = []byte{0x31, 0x41, 0x59, 0x26, 0x53, 0x59}

// bzEOSMagic is the 48-bit end-of-stream magic.
var bzEOSMagic = []byte{0x17, 0x72, 0x45, 0x38, 0x50, 0x90}

// bitWriter accumulates bits MSB-first, exactly as bzip2 emits them.
type bitWriter struct {
	out  []byte
	cur  uint64 // pending bits, left-justified in the low (nbits) bits
	nbit uint   // number of valid bits currently buffered
}

// writeBits appends the low n bits of v, most-significant first.
func (w *bitWriter) writeBits(n uint, v uint32) {
	w.cur = (w.cur << n) | uint64(v&((1<<n)-1))
	w.nbit += n
	for w.nbit >= 8 {
		w.nbit -= 8
		w.out = append(w.out, byte(w.cur>>w.nbit))
	}
}

// writeBit appends a single bit.
func (w *bitWriter) writeBit(b uint32) { w.writeBits(1, b) }

// flush pads the final partial byte with zero bits.
func (w *bitWriter) flush() {
	if w.nbit > 0 {
		w.out = append(w.out, byte(w.cur<<(8-w.nbit)))
		w.nbit = 0
		w.cur = 0
	}
}

// Bzip2Encode compresses src into a complete bzip2 stream at the default
// block size (level 9, 900k — the bzip2 default). The result is a
// standard .bz2 stream that decodes byte-for-byte under compress/bzip2,
// the system bzip2 and htslib/libbz2. It is the in-tree encoder
// counterpart of compress/bzip2's decoder, used for CRAM bzip2 blocks
// (compression method 2) and the arith_dynamic / name-tokeniser X_EXT
// external codec.
func Bzip2Encode(src []byte) ([]byte, error) {
	return bzip2Encode(src, 9)
}

// bzip2Encode compresses src into a complete bzip2 stream at the given
// block-size level (1..9, * 100k bytes per block). Level 9 (900k, the
// bzip2 default) is recommended. The result decodes byte-for-byte under
// any conforming bzip2 decoder.
func bzip2Encode(src []byte, level int) ([]byte, error) {
	if level < 1 || level > 9 {
		return nil, fmt.Errorf("bzip2: block-size level %d out of range 1..9", level)
	}
	blockSize := level * 100000

	w := &bitWriter{}
	// Stream header: "BZh" + level digit.
	w.writeBits(8, 'B')
	w.writeBits(8, 'Z')
	w.writeBits(8, 'h')
	w.writeBits(8, uint32('0'+level))

	var combinedCRC uint32
	// Carve the input into blocks, applying RLE1 per block with a byte
	// budget. A block's RLE1 output must not exceed the working limit
	// (the reference uses blockSize - 19; we leave the same headroom), and
	// crucially RLE1 groups are never split across blocks: rle1EncodeBlock
	// stops on a whole-group boundary so each block round-trips
	// independently.
	limit := blockSize - 20
	if limit < 1 {
		limit = 1
	}
	// A block's symbol stream is at most one MTF/RLE2 symbol per RLE1 byte
	// plus the EOB marker, so its selector count is at most
	// ceil((limit+1)/bzGroupSize). Cap the per-block byte budget so that
	// count can never exceed the 15-bit selector field (bounded here by the
	// reference's bzMaxSelectors), independent of the block-size level.
	if maxByLimit := bzMaxSelectors*bzGroupSize - 1; limit > maxByLimit {
		limit = maxByLimit
	}
	srcOff := 0
	for srcOff < len(src) {
		blk, consumed := rle1EncodeBlock(src[srcOff:], limit)
		// The bzip2 block CRC is computed over the *original* input bytes
		// that make up the block, before the RLE1 transform.
		crc := bzCRC32(src[srcOff : srcOff+consumed])
		srcOff += consumed
		combinedCRC = (combinedCRC<<1 | combinedCRC>>31) ^ crc
		encodeBlock(w, blk, crc)
	}

	// Stream footer.
	for _, b := range bzEOSMagic {
		w.writeBits(8, uint32(b))
	}
	w.writeBits(32, combinedCRC)
	w.flush()
	return w.out, nil
}

// encodeBlock writes one compressed block: magic, CRC, BWT, MTF/RLE2 and
// the multi-table Huffman stream. blk is the post-RLE1 data (1..blockSize
// bytes) and crc is the bzip2 CRC32 of the block's *original* (pre-RLE1)
// bytes, computed by the caller.
func encodeBlock(w *bitWriter, blk []byte, crc uint32) {
	for _, b := range bzBlockMagic {
		w.writeBits(8, uint32(b))
	}
	w.writeBits(32, crc)
	w.writeBit(0) // not randomised

	bwt, origPtr := bwtEncode(blk)
	w.writeBits(24, uint32(origPtr))

	// Build the symbol map: which byte values occur in the BWT output.
	var inUse [256]bool
	for _, b := range bwt {
		inUse[b] = true
	}
	// unseqToSeq maps a present byte value to its 0-based MTF rank seed
	// (the initial MTF list is the sorted list of present bytes).
	var seqToUnseq []byte
	for i := 0; i < 256; i++ {
		if inUse[i] {
			seqToUnseq = append(seqToUnseq, byte(i))
		}
	}
	nInUse := len(seqToUnseq)
	// Alphabet: RUNA=0, RUNB=1, then literals 2..nInUse, EOB=nInUse+1.
	alphaSize := nInUse + 2

	// MTF + RLE2: turn the BWT bytes into the MTF/RLE2 symbol stream.
	mtfv := mtfAndRLE2(bwt, seqToUnseq, nInUse)
	// mtfv ends with the EOB symbol already appended.

	// Build the Huffman tables and per-group selectors.
	nGroups := chooseNumTables(len(mtfv))
	lengths, selectors := buildHuffmanTables(mtfv, alphaSize, nGroups)

	// Emit the symbol map.
	var inUse16 [16]bool
	for i := 0; i < 256; i++ {
		if inUse[i] {
			inUse16[i>>4] = true
		}
	}
	var map16 uint32
	for i := 0; i < 16; i++ {
		if inUse16[i] {
			map16 |= 1 << uint(15-i)
		}
	}
	w.writeBits(16, map16)
	for i := 0; i < 16; i++ {
		if !inUse16[i] {
			continue
		}
		var bits uint32
		for j := 0; j < 16; j++ {
			if inUse[i*16+j] {
				bits |= 1 << uint(15-j)
			}
		}
		w.writeBits(16, bits)
	}

	// Number of Huffman tables and selectors.
	w.writeBits(3, uint32(nGroups))
	w.writeBits(15, uint32(len(selectors)))

	// Selectors, MTF-coded then emitted in unary.
	var selMTF [6]byte
	for i := byte(0); int(i) < nGroups; i++ {
		selMTF[i] = i
	}
	for _, s := range selectors {
		// Find s's position in the MTF list, emit that position in unary,
		// then move s to the front.
		var j byte
		for selMTF[j] != s {
			j++
		}
		for k := j; k > 0; k-- {
			selMTF[k] = selMTF[k-1]
		}
		selMTF[0] = s
		for ; j > 0; j-- {
			w.writeBit(1)
		}
		w.writeBit(0)
	}

	// Per-table code lengths, delta-coded from a running value.
	for t := 0; t < nGroups; t++ {
		curr := int(lengths[t][0])
		w.writeBits(5, uint32(curr))
		for s := 0; s < alphaSize; s++ {
			target := int(lengths[t][s])
			for curr < target {
				w.writeBit(1) // 10 = increment
				w.writeBit(0)
				curr++
			}
			for curr > target {
				w.writeBit(1) // 11 = decrement
				w.writeBit(1)
				curr--
			}
			w.writeBit(0) // 0 = stop
		}
	}

	// Build canonical codes for each table.
	var codes [6][bzMaxAlphaSize]uint32
	for t := 0; t < nGroups; t++ {
		assignCanonicalCodes(lengths[t][:alphaSize], codes[t][:alphaSize])
	}

	// Emit the symbol stream: for each group of bzGroupSize symbols use
	// the selected table.
	groupNo := 0
	for i := 0; i < len(mtfv); {
		t := selectors[groupNo]
		groupNo++
		end := i + bzGroupSize
		if end > len(mtfv) {
			end = len(mtfv)
		}
		for ; i < end; i++ {
			sym := mtfv[i]
			w.writeBits(uint(lengths[t][sym]), codes[t][sym])
		}
	}
}

// rle1Encode applies bzip2's first run-length stage to the whole input
// with no size limit. A run of 4..255 identical bytes is written as four
// copies of the byte followed by a single length byte counting the
// *extra* repeats (0..251). Runs longer than 255 restart. This is the
// unbounded form used by the unit tests; the encoder proper uses
// rle1EncodeBlock to respect block boundaries.
func rle1Encode(src []byte) []byte {
	out, _ := rle1EncodeBlock(src, 1<<62)
	return out
}

// rle1EncodeBlock RLE1-encodes a prefix of src, stopping once a further
// group would push the output past limit bytes. It returns the encoded
// bytes and the number of *source* bytes consumed, always ending on a
// whole RLE1 group so the block round-trips independently. At least one
// group is always emitted so the loop makes progress even when limit is
// tiny.
func rle1EncodeBlock(src []byte, limit int) (out []byte, consumed int) {
	n := len(src)
	out = make([]byte, 0, min(n, limit)+16)
	i := 0
	for i < n {
		b := src[i]
		run := 1
		for i+run < n && run < 255 && src[i+run] == b {
			run++
		}
		var groupLen int
		if run >= 4 {
			groupLen = 5
		} else {
			groupLen = run
		}
		// Stop before exceeding the budget, but always emit at least one
		// group.
		if len(out)+groupLen > limit && len(out) > 0 {
			break
		}
		if run >= 4 {
			out = append(out, b, b, b, b, byte(run-4))
		} else {
			for k := 0; k < run; k++ {
				out = append(out, b)
			}
		}
		i += run
	}
	return out, i
}

// mtfAndRLE2 performs move-to-front over the BWT output restricted to the
// bytes actually present, then RLE2-encodes runs of the MTF zero using
// the RUNA/RUNB bijective base-2 code. Non-zero MTF ranks r are emitted
// as r+1 (shifting literals above the two run symbols). The EOB symbol is
// appended at the end. seqToUnseq lists the present byte values sorted
// ascending (the initial MTF list); nInUse is its length.
func mtfAndRLE2(bwt []byte, seqToUnseq []byte, nInUse int) []uint16 {
	// MTF list seeded with the present bytes in ascending order.
	mtf := make([]byte, nInUse)
	copy(mtf, seqToUnseq)
	// unseqToSeq maps a byte value to its index in the MTF list.
	var unseqToSeq [256]byte
	for i, v := range seqToUnseq {
		unseqToSeq[v] = byte(i)
	}

	out := make([]uint16, 0, len(bwt)+1)
	zeroRun := 0
	flushZeros := func() {
		// Emit the zero run in bijective base-2, mirroring the reference
		// generateMTFValues: decrement once, then loop emitting RUNA for
		// an even bit and RUNB for an odd bit, halving until done.
		if zeroRun == 0 {
			return
		}
		z := zeroRun - 1
		for {
			if z&1 != 0 {
				out = append(out, bzRunB)
			} else {
				out = append(out, bzRunA)
			}
			if z < 2 {
				break
			}
			z = (z - 2) / 2
		}
		zeroRun = 0
	}

	for _, b := range bwt {
		// Find b's current rank in the MTF list and move it to front.
		j := unseqToSeq[b]
		if j == 0 {
			zeroRun++
			continue
		}
		flushZeros()
		// Move-to-front: shift mtf[0..j] down, recompute unseqToSeq for
		// the shifted entries.
		tmp := mtf[j]
		for k := j; k > 0; k-- {
			mtf[k] = mtf[k-1]
			unseqToSeq[mtf[k]] = k
		}
		mtf[0] = tmp
		unseqToSeq[tmp] = 0
		// Literal rank j is emitted as j+1 (above RUNA/RUNB).
		out = append(out, uint16(j)+1)
	}
	flushZeros()
	out = append(out, uint16(nInUse+1)) // EOB
	return out
}

// chooseNumTables picks the number of Huffman tables for a symbol stream
// of the given length, following the reference heuristic.
func chooseNumTables(nMTF int) int {
	switch {
	case nMTF < 200:
		return 2
	case nMTF < 600:
		return 3
	case nMTF < 1200:
		return 4
	case nMTF < 2400:
		return 5
	default:
		return 6
	}
}

// buildHuffmanTables runs the bzip2 multi-table optimisation: it seeds
// nGroups tables by partitioning the alphabet by cumulative frequency,
// then iterates (selector assignment -> per-table cost -> rebuild
// lengths) a few times, mirroring BZ2_compressBlock. It returns the
// per-table code lengths and the per-group selector list.
func buildHuffmanTables(mtfv []uint16, alphaSize, nGroups int) ([6][bzMaxAlphaSize]uint8, []byte) {
	var lengths [6][bzMaxAlphaSize]uint8
	nMTF := len(mtfv)

	// Overall symbol frequencies.
	freq := make([]int, alphaSize)
	for _, s := range mtfv {
		freq[s]++
	}

	// Seed tables: split the alphabet into nGroups frequency bands so each
	// band carries roughly nMTF/nGroups symbols. Symbols inside a band get
	// length 0 (cheap), outside get a large length (expensive). This is the
	// reference's initial-table heuristic.
	nPart := nGroups
	remF := nMTF
	gs := 0
	ge := -1
	for nPart > 0 {
		tFreq := remF / nPart
		ge = gs - 1
		aFreq := 0
		for aFreq < tFreq && ge < alphaSize-1 {
			ge++
			aFreq += freq[ge]
		}
		if ge > gs && nPart != nGroups && nPart != 1 && ((nGroups-nPart)%2 == 1) {
			aFreq -= freq[ge]
			ge--
		}
		t := nGroups - nPart
		for v := 0; v < alphaSize; v++ {
			if v >= gs && v <= ge {
				lengths[t][v] = 0
			} else {
				lengths[t][v] = 15
			}
		}
		nPart--
		gs = ge + 1
		remF -= aFreq
	}

	nGroupsBlocks := (nMTF + bzGroupSize - 1) / bzGroupSize
	selectors := make([]byte, nGroupsBlocks)

	// Iterate to refine.
	const nIters = 4
	for iter := 0; iter < nIters; iter++ {
		// Per-table per-symbol frequencies accumulated over the groups
		// that select each table.
		var rfreq [6][bzMaxAlphaSize]int

		gi := 0
		for start := 0; start < nMTF; start += bzGroupSize {
			end := start + bzGroupSize
			if end > nMTF {
				end = nMTF
			}
			// Compute each table's bit cost for this group.
			var cost [6]int
			for t := 0; t < nGroups; t++ {
				c := 0
				for k := start; k < end; k++ {
					c += int(lengths[t][mtfv[k]])
				}
				cost[t] = c
			}
			// Pick the cheapest table.
			bt := 0
			bc := cost[0]
			for t := 1; t < nGroups; t++ {
				if cost[t] < bc {
					bc = cost[t]
					bt = t
				}
			}
			selectors[gi] = byte(bt)
			gi++
			for k := start; k < end; k++ {
				rfreq[bt][mtfv[k]]++
			}
		}

		// Rebuild each table's lengths from its accumulated frequencies.
		for t := 0; t < nGroups; t++ {
			buildCodeLengths(rfreq[t][:alphaSize], lengths[t][:alphaSize], bzMaxCodeLen)
		}
	}

	return lengths, selectors
}

// buildCodeLengths computes canonical Huffman code lengths for the given
// symbol frequencies, capped at maxLen, using the package's Huffman
// length builder. Zero-frequency symbols are given a small nonzero
// weight so every symbol receives a valid length (the reference does the
// same so all tables can encode any symbol).
func buildCodeLengths(freq []int, lengths []uint8, maxLen int) {
	n := len(freq)
	weights := make([]int, n)
	for i := 0; i < n; i++ {
		if freq[i] == 0 {
			weights[i] = 1
		} else {
			weights[i] = freq[i]
		}
	}
	bzHuffmanLengths(weights, lengths, maxLen)
}

// assignCanonicalCodes fills codes[] with canonical Huffman codes derived
// from the given per-symbol bit lengths, in bzip2's MSB-first canonical
// order (shortest codes first, ascending symbol index within a length).
func assignCanonicalCodes(lengths []uint8, codes []uint32) {
	minLen := uint8(32)
	maxLen := uint8(0)
	for _, l := range lengths {
		if l < minLen {
			minLen = l
		}
		if l > maxLen {
			maxLen = l
		}
	}
	var code uint32
	for bits := minLen; bits <= maxLen; bits++ {
		for i, l := range lengths {
			if l == bits {
				codes[i] = code
				code++
			}
		}
		code <<= 1
	}
}

// bzHuffmanLengths computes length-limited Huffman code lengths for the
// given positive weights, capped at maxLen. It mirrors
// BZ2_hbMakeCodeLengths: build a Huffman tree, and if any code exceeds
// maxLen, scale the weights down and retry. lengths is filled in place.
func bzHuffmanLengths(weights []int, lengths []uint8, maxLen int) {
	n := len(weights)
	if n == 0 {
		return
	}
	if n == 1 {
		lengths[0] = 1
		return
	}

	// Working weights, with frequency in the high bits and a tie-break in
	// the low bits (matching the reference's WEIGHTOF/DEPTHOF packing is
	// not required for correctness of *lengths*; we just need a valid
	// length-limited Huffman code, which decodes identically).
	const (
		minWeight = 1
	)
	w := make([]int, n)
	for i := range weights {
		if weights[i] < minWeight {
			w[i] = minWeight
		} else {
			w[i] = weights[i]
		}
	}

	for {
		ok := huffmanBuild(w, lengths, maxLen)
		if ok {
			return
		}
		// Scale weights down (reference: w = (w>>8 ... ) but we use a
		// gentler divide-and-bump) and retry until within maxLen.
		for i := range w {
			w[i] = (w[i] + 1) >> 1
			if w[i] < minWeight {
				w[i] = minWeight
			}
		}
	}
}

// huffNode is a node in the Huffman build heap/tree.
type huffNode struct {
	weight int
	parent int
}

// huffmanBuild constructs a Huffman tree from the weights and writes the
// resulting per-symbol code lengths. It returns false if any code length
// exceeds maxLen (so the caller can rescale and retry).
func huffmanBuild(weights []int, lengths []uint8, maxLen int) bool {
	n := len(weights)
	// Node array: leaves 0..n-1, internal nodes appended.
	nodes := make([]huffNode, n, 2*n)
	// Heap of node indices ordered by weight (ascending), tie-broken by
	// index for determinism.
	heap := make([]int, n)
	for i := 0; i < n; i++ {
		nodes[i] = huffNode{weight: weights[i], parent: -1}
		heap[i] = i
	}
	less := func(a, b int) bool {
		if nodes[a].weight != nodes[b].weight {
			return nodes[a].weight < nodes[b].weight
		}
		return a < b
	}
	sort.Slice(heap, func(i, j int) bool { return less(heap[i], heap[j]) })

	// Repeatedly combine the two lightest nodes.
	for len(heap) > 1 {
		a := heap[0]
		b := heap[1]
		heap = heap[2:]
		ni := len(nodes)
		nodes = append(nodes, huffNode{weight: nodes[a].weight + nodes[b].weight, parent: -1})
		nodes[a].parent = ni
		nodes[b].parent = ni
		// Insert ni into the (sorted) heap.
		pos := sort.Search(len(heap), func(i int) bool { return !less(heap[i], ni) })
		heap = append(heap, 0)
		copy(heap[pos+1:], heap[pos:])
		heap[pos] = ni
	}

	// Derive code lengths by walking each leaf to the root.
	for i := 0; i < n; i++ {
		l := 0
		p := nodes[i].parent
		for p != -1 {
			l++
			p = nodes[p].parent
		}
		if l == 0 {
			l = 1 // single-symbol edge case
		}
		if l > maxLen {
			return false
		}
		lengths[i] = uint8(l)
	}
	return true
}

// bwtEncode computes the Burrows-Wheeler transform of data: it sorts the
// n cyclic rotations and returns the last column plus the row index of
// the original (unrotated) string. Sorting uses rank-doubling
// (Manber-Myers): rotations are ranked by their first 1, 2, 4, ... 2^k
// characters until all ranks are distinct, in O(n log n) total — robust
// against highly repetitive input where a naive rotation comparison sort
// degrades to O(n^2 log n). The result is the exact bzip2 BWT and
// decodes byte-for-byte under any conforming decoder.
func bwtEncode(data []byte) (bwt []byte, origPtr int) {
	n := len(data)
	if n == 0 {
		return nil, 0
	}
	if n == 1 {
		return []byte{data[0]}, 0
	}

	sa := sortCyclicRotations(data)

	bwt = make([]byte, n)
	for i, r := range sa {
		bwt[i] = data[(r+n-1)%n]
		if r == 0 {
			origPtr = i
		}
	}
	return bwt, origPtr
}

// sortCyclicRotations returns the starting indices of data's cyclic
// rotations in lexicographic order, using rank-doubling. Equal rotations
// (possible only when the string is periodic) are ordered by start index
// to match a stable rotation sort; the inverse BWT recovers the original
// regardless of how ties among identical rotations are broken.
func sortCyclicRotations(data []byte) []int {
	n := len(data)
	sa := make([]int, n)
	rank := make([]int, n)
	tmp := make([]int, n)

	// Initial ranks from the single leading byte.
	for i := 0; i < n; i++ {
		sa[i] = i
		rank[i] = int(data[i])
	}

	for k := 1; ; k <<= 1 {
		// Comparator on (rank[i], rank[i+k mod n]) pairs.
		cmp := func(a, b int) int {
			if rank[a] != rank[b] {
				if rank[a] < rank[b] {
					return -1
				}
				return 1
			}
			ra := rank[(a+k)%n]
			rb := rank[(b+k)%n]
			if ra != rb {
				if ra < rb {
					return -1
				}
				return 1
			}
			return 0
		}
		sort.Slice(sa, func(i, j int) bool { return cmp(sa[i], sa[j]) < 0 })

		// Recompute ranks from the sorted order.
		tmp[sa[0]] = 0
		distinct := 0
		for i := 1; i < n; i++ {
			if cmp(sa[i-1], sa[i]) < 0 {
				distinct++
			}
			tmp[sa[i]] = distinct
		}
		copy(rank, tmp)
		if distinct == n-1 {
			// All ranks distinct: rotations are fully ordered.
			break
		}
		if k >= n {
			// Ranks have stopped changing yet ties remain (the string is
			// periodic, so some rotations are genuinely identical). The
			// current sa already orders distinct rotations correctly;
			// resolve remaining ties by start index for determinism.
			sort.SliceStable(sa, func(i, j int) bool {
				if rank[sa[i]] != rank[sa[j]] {
					return rank[sa[i]] < rank[sa[j]]
				}
				return sa[i] < sa[j]
			})
			break
		}
	}
	return sa
}

// bzCRC32Table is the bzip2 CRC32 lookup table (MSB-first, polynomial
// 0x04C11DB7), distinct from the IEEE CRC used elsewhere.
var bzCRC32Table = func() [256]uint32 {
	var t [256]uint32
	for i := uint32(0); i < 256; i++ {
		c := i << 24
		for k := 0; k < 8; k++ {
			if c&0x80000000 != 0 {
				c = (c << 1) ^ 0x04C11DB7
			} else {
				c <<= 1
			}
		}
		t[i] = c
	}
	return t
}()

// bzCRC32 computes the bzip2 block/stream CRC32 of data.
func bzCRC32(data []byte) uint32 {
	crc := uint32(0xFFFFFFFF)
	for _, b := range data {
		crc = (crc << 8) ^ bzCRC32Table[((crc>>24)^uint32(b))&0xFF]
	}
	return ^crc
}
