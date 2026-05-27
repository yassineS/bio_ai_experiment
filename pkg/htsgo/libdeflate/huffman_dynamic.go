package libdeflate

// Length-limited canonical Huffman code construction, ported from
// libdeflate's deflate_make_huffman_code family in
// reference_code/libdeflate/lib/deflate_compress.c:1318-1396 (the
// public entry point), supported by sort_symbols (846), build_tree
// (940), compute_length_counts (1023), and gen_codewords (1178). The
// algorithm is a Huffman-construction-then-cap variant: build the
// classical tree (using the standard sort-then-merge optimization),
// extract per-depth leaf counts, and clamp any depth > max_codeword_len
// down to the largest already-used depth. libdeflate's authors note
// the clamp is not strictly optimal, but it is good enough for DEFLATE
// (the inflation in code-length cost is well below 1% on real data)
// and — crucially for us — it is what produces the exact code lengths
// the BTRACE HUFFMAN_LEN events record.

// numSymbolBits matches NUM_SYMBOL_BITS (deflate_compress.c:815). The
// sort/tree-build phase packs a (frequency, symbol) pair into a single
// uint32 with the symbol in the low NUM_SYMBOL_BITS bits.
const (
	numSymbolBits = 10
	symbolMask    = (1 << numSymbolBits) - 1
	freqMask      = ^uint32(symbolMask)
)

// makeHuffmanCode builds a length-limited canonical Huffman code for
// the given symbol-frequency table. It writes per-symbol code lengths
// (0 for unused symbols) into lens and bit-reversed canonical
// codewords into codes. The output is byte-identical to libdeflate's
// deflate_make_huffman_code on the same input.
func makeHuffmanCode(numSyms, maxCodewordLen int, freqs []uint32, lens []uint8, codes []uint32) {
	// Sort symbols primarily by frequency, secondarily by symbol value
	// (deflate_compress.c:846 sort_symbols). The packed array shares
	// storage with the eventual codewords output.
	numUsed := sortSymbols(numSyms, freqs, lens, codes)

	// A complete Huffman code must contain at least 2 codewords; pad
	// when the input has fewer than 2 used symbols
	// (deflate_compress.c:1369-1378).
	if numUsed < 2 {
		var sym int
		if numUsed == 1 {
			sym = int(codes[0]) & symbolMask
		}
		nonzeroIdx := sym
		if sym == 0 {
			nonzeroIdx = 1
		}
		for i := range codes[:numSyms] {
			codes[i] = 0
		}
		for i := range lens[:numSyms] {
			lens[i] = 0
		}
		codes[0] = 0
		lens[0] = 1
		codes[nonzeroIdx] = 1
		lens[nonzeroIdx] = 1
		return
	}

	// Build the stripped-down Huffman tree in-place. After build_tree
	// the first numUsed-1 entries hold parent-pointer-encoded non-leaf
	// nodes, with the root at index numUsed-2
	// (deflate_compress.c:940).
	buildTree(codes, numUsed)

	// Extract per-length leaf counts from the tree, capping any
	// depth > maxCodewordLen to the largest used depth.
	lenCounts := make([]uint32, maxCodewordLen+1)
	computeLengthCounts(codes, numUsed-2, lenCounts, maxCodewordLen)

	// Emit code lengths and canonical bit-reversed codewords.
	genCodewords(codes, lens, lenCounts, maxCodewordLen, numSyms)
}

// sortSymbols mirrors sort_symbols (deflate_compress.c:846). It packs
// each used (sym, freq) into a uint32 with freq in the high bits and
// sym in the low NUM_SYMBOL_BITS bits, writes them into symout in
// ascending (freq, sym) order, and returns the number used. Symbols
// with zero frequency get lens[sym] = 0.
//
// libdeflate uses a count-sort-then-heapsort hybrid for performance.
// Since correctness is what matters here, we use a straight stable
// sort by (freq, sym). The output ordering is identical because the
// hybrid is also (freq, sym) ascending.
func sortSymbols(numSyms int, freqs []uint32, lens []uint8, symout []uint32) int {
	numUsed := 0
	for sym := 0; sym < numSyms; sym++ {
		f := freqs[sym]
		if f == 0 {
			lens[sym] = 0
			continue
		}
		symout[numUsed] = uint32(sym) | (f << numSymbolBits)
		numUsed++
	}
	// Insertion sort would be O(n^2); since DEFLATE_MAX_NUM_SYMS=288
	// is tiny, a simple insertion sort is plenty fast. But we use
	// stdlib sort.Slice's sort.Sort to be safe.
	used := symout[:numUsed]
	// We want ascending order by (freq, sym). Because (freq << 10 | sym)
	// already encodes that order, a numeric sort suffices.
	sortUint32Ascending(used)
	return numUsed
}

// sortUint32Ascending sorts a in ascending order. We avoid pulling in
// sort.Slice's reflection overhead by hand-rolling a quick sort over
// uint32. The slice length is bounded by numLitlenSyms=288.
func sortUint32Ascending(a []uint32) {
	if len(a) < 2 {
		return
	}
	quicksortU32(a, 0, len(a)-1)
}

func quicksortU32(a []uint32, lo, hi int) {
	for lo < hi {
		if hi-lo < 12 {
			// Insertion sort for small ranges.
			for i := lo + 1; i <= hi; i++ {
				v := a[i]
				j := i - 1
				for j >= lo && a[j] > v {
					a[j+1] = a[j]
					j--
				}
				a[j+1] = v
			}
			return
		}
		pivot := medianOfThree(a, lo, hi)
		i, j := lo, hi
		for i <= j {
			for a[i] < pivot {
				i++
			}
			for a[j] > pivot {
				j--
			}
			if i <= j {
				a[i], a[j] = a[j], a[i]
				i++
				j--
			}
		}
		// Recurse on the smaller partition, iterate the larger.
		if j-lo < hi-i {
			quicksortU32(a, lo, j)
			lo = i
		} else {
			quicksortU32(a, i, hi)
			hi = j
		}
	}
}

func medianOfThree(a []uint32, lo, hi int) uint32 {
	mid := lo + (hi-lo)/2
	x, y, z := a[lo], a[mid], a[hi]
	switch {
	case x <= y && y <= z, z <= y && y <= x:
		return y
	case y <= x && x <= z, z <= x && x <= y:
		return x
	default:
		return z
	}
}

// buildTree mirrors build_tree (deflate_compress.c:940). On entry A
// holds numUsed packed (sym, freq) values in ascending order. On
// exit A[0..numUsed-2] are non-leaf nodes; A[i] stores its parent's
// index in the high bits with the original symbol preserved in the
// low NUM_SYMBOL_BITS bits.
func buildTree(A []uint32, numUsed int) {
	if numUsed < 2 {
		return
	}
	lastIdx := numUsed - 1

	// i = next lowest-frequency leaf needing a parent.
	// b = next lowest-frequency non-leaf needing a parent (or e if none).
	// e = next spot for a new non-leaf (overwrites a leaf).
	i, b, e := 0, 0, 0

	for {
		var newFreq uint32

		// Pick the two lowest-frequency among leaves A[i], A[i+1]
		// and non-leaves A[b], A[b+1].
		leafLeaf := i+1 <= lastIdx &&
			(b == e || (A[i+1]&freqMask) <= (A[b]&freqMask))
		nonleafNonleaf := !leafLeaf && b+2 <= e &&
			(i > lastIdx || (A[b+1]&freqMask) < (A[i]&freqMask))

		switch {
		case leafLeaf:
			newFreq = (A[i] & freqMask) + (A[i+1] & freqMask)
			i += 2
		case nonleafNonleaf:
			newFreq = (A[b] & freqMask) + (A[b+1] & freqMask)
			A[b] = (uint32(e) << numSymbolBits) | (A[b] & symbolMask)
			A[b+1] = (uint32(e) << numSymbolBits) | (A[b+1] & symbolMask)
			b += 2
		default:
			// One leaf and one non-leaf.
			newFreq = (A[i] & freqMask) + (A[b] & freqMask)
			A[b] = (uint32(e) << numSymbolBits) | (A[b] & symbolMask)
			i++
			b++
		}
		A[e] = newFreq | (A[e] & symbolMask)
		e++
		if e >= lastIdx {
			break
		}
	}
}

// computeLengthCounts mirrors compute_length_counts
// (deflate_compress.c:1023). It walks the build_tree output in reverse
// order so parents are visited before children, accumulating per-depth
// leaf counts. The length-limited constraint is enforced by clamping
// any depth >= max_codeword_len down to the deepest non-zero existing
// depth (deflate_compress.c:1077).
func computeLengthCounts(A []uint32, rootIdx int, lenCounts []uint32, maxCodewordLen int) {
	for len := 0; len <= maxCodewordLen; len++ {
		lenCounts[len] = 0
	}
	lenCounts[1] = 2

	// Root has depth 0.
	A[rootIdx] &= symbolMask

	for node := rootIdx - 1; node >= 0; node-- {
		parent := int(A[node] >> numSymbolBits)
		parentDepth := int(A[parent] >> numSymbolBits)
		depth := parentDepth + 1

		A[node] = (A[node] & symbolMask) | (uint32(depth) << numSymbolBits)

		if depth >= maxCodewordLen {
			depth = maxCodewordLen
			for lenCounts[depth] == 0 {
				depth--
			}
		}

		lenCounts[depth]--
		lenCounts[depth+1] += 2
	}
}

// genCodewords mirrors gen_codewords (deflate_compress.c:1178). It
// assigns codeword lengths to symbols (in increasing-frequency,
// increasing-symbol order, with the longest lengths going to the
// least-frequent symbols), then derives canonical codewords and
// bit-reverses them so the LSB-first bitstream writer can emit them
// directly.
func genCodewords(A []uint32, lens []uint8, lenCounts []uint32, maxCodewordLen, numSyms int) {
	// Assign lengths in decreasing-length order to symbols sorted by
	// increasing (freq, sym). After build_tree the symbol payload is
	// still preserved in the low NUM_SYMBOL_BITS bits.
	i := 0
	for length := maxCodewordLen; length >= 1; length-- {
		count := lenCounts[length]
		for count > 0 {
			lens[int(A[i])&symbolMask] = uint8(length)
			i++
			count--
		}
	}

	// Build the lexicographically-first codeword of each length.
	nextCodewords := make([]uint32, maxCodewordLen+1)
	nextCodewords[0] = 0
	nextCodewords[1] = 0
	for length := 2; length <= maxCodewordLen; length++ {
		nextCodewords[length] = (nextCodewords[length-1] + lenCounts[length-1]) << 1
	}

	// Emit per-symbol bit-reversed canonical codewords.
	for sym := 0; sym < numSyms; sym++ {
		l := lens[sym]
		if l == 0 {
			A[sym] = 0
			continue
		}
		A[sym] = reverseBits(nextCodewords[l], l)
		nextCodewords[l]++
	}
}
