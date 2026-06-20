package tabix

// binindex.go reproduces htslib's on-disk bin ordering and bin-compression
// so that the .bai/.csi/.tbi files we write are byte-identical to those
// produced by `samtools index`, `tabix`, and `bcftools index`.
//
// htslib stores each reference's bins in a khash hash table
// (KHASH_MAP_INIT_INT, identity hash over the 32-bit bin number) and
// serialises them in raw hash-table iteration order — not sorted by bin
// number. It also runs hts_idx_finish -> compress_binning, which merges a
// bin's chunks into its parent bin whenever the bin's chunks span fewer
// than HTS_MIN_MARKER_DIST (0x10000) compressed BGZF bytes, then deletes the
// merged bin. Reproducing both behaviours requires a faithful khash port:
// the iteration order after deletions depends on the exact bucket layout,
// which in turn depends on the insertion order and the kick-out resize
// policy.
//
// This file ports the integer-keyed khash map (identity hash, 32-bit key)
// plus compress_binning. It is the single source of bin ordering for every
// index writer in the repo (BAI, BAM CSI, TBI, tabix CSI).

import "sort"

// htsMinMarkerDist mirrors htslib's HTS_MIN_MARKER_DIST: a bin whose chunk
// list spans fewer than this many compressed BGZF bytes is folded into its
// parent bin by compress_binning.
const htsMinMarkerDist = 0x10000

// binChunk is one (beg, end) virtual-offset pair within a bin. It mirrors
// htslib's hts_pair64_t.
type binChunk struct {
	Beg, End VOffset
}

// binEntry is one bin together with its chunk list and CSI loffset. The
// builders hand these to the bin index in first-appearance (insertion)
// order; the index reproduces htslib's khash layout from that order.
type binEntry struct {
	ID      uint32
	LOffset VOffset
	Chunks  []binChunk
}

// BinChunk is the exported form of a bin's (Beg, End) virtual-offset pair,
// used by the BAI/CSI builders in other packages that drive OrderBins.
type BinChunk = binChunk

// BinEntry is the exported form of one bin (ID, loffset, chunks) handed to
// OrderBins by an index builder in first-appearance order.
type BinEntry = binEntry

// OrderBins reproduces htslib's on-disk bin ordering and bin compression for
// one reference. inserted holds the reference's bins in first-appearance
// (insertion) order, including the meta pseudo-bin last; nLvls is the bin
// hierarchy depth; nBins is the regular-bin count (the boundary below which a
// bin is "regular" and eligible for compression); metaBin is the pseudo-bin
// number. The returned slice is in the exact order htslib would serialise.
func OrderBins(inserted []BinEntry, nLvls int, nBins, metaBin uint32) []BinEntry {
	return orderBins(inserted, nLvls, nBins, metaBin)
}

// htsBinFirst returns the first bin number at hierarchy level l, matching
// htslib's hts_bin_first(l) = ((1<<(3*l)) - 1) / 7.
func htsBinFirst(l int) uint32 { return uint32((uint64(1)<<uint(3*l) - 1) / 7) }

// htsBinParent returns the parent of bin b: (b-1) >> 3, matching htslib's
// hts_bin_parent.
func htsBinParent(b uint32) uint32 { return (b - 1) >> 3 }

// kroundup32 rounds x up to the next power of two, matching htslib's macro.
func kroundup32(x uint32) uint32 {
	x--
	x |= x >> 1
	x |= x >> 2
	x |= x >> 4
	x |= x >> 8
	x |= x >> 16
	x++
	return x
}

// binKhash is a faithful port of khash.h's KHASH_MAP_INIT_INT specialised to
// the bin index: 32-bit keys, identity hash, value = *binEntry. Its bucket
// layout and resize/kick-out behaviour match khash.h bit-for-bit so that the
// serialisation order (a plain walk of buckets 0..n_buckets-1) reproduces
// upstream byte-for-byte.
type binKhash struct {
	nBuckets   uint32
	size       uint32
	nOccupied  uint32
	upperBound uint32
	flags      []uint32
	keys       []uint32
	vals       []*binEntry
}

func khFsize(m uint32) uint32 {
	if m < 16 {
		return 1
	}
	return m >> 4
}

func khFlagIsEmpty(f []uint32, i uint32) bool  { return (f[i>>4]>>((i&0xf)<<1))&2 != 0 }
func khFlagIsDel(f []uint32, i uint32) bool    { return (f[i>>4]>>((i&0xf)<<1))&1 != 0 }
func khFlagIsEither(f []uint32, i uint32) bool { return (f[i>>4]>>((i&0xf)<<1))&3 != 0 }
func khSetBothFalse(f []uint32, i uint32)      { f[i>>4] &^= 3 << ((i & 0xf) << 1) }
func khSetEmptyFalse(f []uint32, i uint32)     { f[i>>4] &^= 2 << ((i & 0xf) << 1) }
func khSetDelTrue(f []uint32, i uint32)        { f[i>>4] |= 1 << ((i & 0xf) << 1) }

func (h *binKhash) resize(newBuckets uint32) {
	const hashUpper = 0.77
	newBuckets = kroundup32(newBuckets)
	if newBuckets < 4 {
		newBuckets = 4
	}
	if float64(h.size) >= float64(newBuckets)*hashUpper+0.5 {
		return
	}
	newFlags := make([]uint32, khFsize(newBuckets))
	for i := range newFlags {
		newFlags[i] = 0xaaaaaaaa
	}
	if h.nBuckets < newBuckets {
		newKeys := make([]uint32, newBuckets)
		copy(newKeys, h.keys)
		h.keys = newKeys
		newVals := make([]*binEntry, newBuckets)
		copy(newVals, h.vals)
		h.vals = newVals
	}
	mask := newBuckets - 1
	for j := uint32(0); j < h.nBuckets; j++ {
		if khFlagIsEither(h.flags, j) {
			continue
		}
		key := h.keys[j]
		val := h.vals[j]
		khSetDelTrue(h.flags, j)
		for {
			i := key & mask
			step := uint32(0)
			for !khFlagIsEmpty(newFlags, i) {
				step++
				i = (i + step) & mask
			}
			khSetEmptyFalse(newFlags, i)
			if i < h.nBuckets && !khFlagIsEither(h.flags, i) {
				key, h.keys[i] = h.keys[i], key
				val, h.vals[i] = h.vals[i], val
				khSetDelTrue(h.flags, i)
			} else {
				h.keys[i] = key
				h.vals[i] = val
				break
			}
		}
	}
	if h.nBuckets > newBuckets {
		h.keys = h.keys[:newBuckets]
		h.vals = h.vals[:newBuckets]
	}
	h.flags = newFlags
	h.nBuckets = newBuckets
	h.nOccupied = h.size
	h.upperBound = uint32(float64(newBuckets)*hashUpper + 0.5)
}

// put inserts key and returns its bucket index plus whether it was newly
// added. Matches kh_put_##name in khash.h (identity hash).
func (h *binKhash) put(key uint32) (uint32, bool) {
	if h.nOccupied >= h.upperBound {
		if h.nBuckets > h.size*2 {
			h.resize(h.nBuckets - 1)
		} else {
			h.resize(h.nBuckets + 1)
		}
	}
	mask := h.nBuckets - 1
	x := key & mask
	site := h.nBuckets
	if !khFlagIsEmpty(h.flags, x) {
		last := x
		step := uint32(0)
		for !khFlagIsEmpty(h.flags, x) && (khFlagIsDel(h.flags, x) || h.keys[x] != key) {
			if khFlagIsDel(h.flags, x) {
				site = x
			}
			step++
			x = (x + step) & mask
			if x == last {
				x = site
				break
			}
		}
		if khFlagIsEmpty(h.flags, x) && site != h.nBuckets {
			x = site
		}
	}
	if khFlagIsEmpty(h.flags, x) {
		khSetBothFalse(h.flags, x)
		h.keys[x] = key
		h.size++
		h.nOccupied++
		return x, true
	}
	if khFlagIsDel(h.flags, x) {
		khSetBothFalse(h.flags, x)
		h.keys[x] = key
		h.size++
		return x, true
	}
	return x, false
}

// get returns key's bucket index and true if present, else (nBuckets, false).
func (h *binKhash) get(key uint32) (uint32, bool) {
	if h.nBuckets == 0 {
		return 0, false
	}
	mask := h.nBuckets - 1
	x := key & mask
	last := x
	step := uint32(0)
	for !khFlagIsEmpty(h.flags, x) && (khFlagIsDel(h.flags, x) || h.keys[x] != key) {
		step++
		x = (x + step) & mask
		if x == last {
			return h.nBuckets, false
		}
	}
	if khFlagIsEither(h.flags, x) {
		return h.nBuckets, false
	}
	return x, true
}

// del marks bucket k deleted. Matches kh_del_##name.
func (h *binKhash) del(k uint32) {
	if k < h.nBuckets && !khFlagIsEither(h.flags, k) {
		khSetDelTrue(h.flags, k)
		h.size--
	}
}

// orderBins takes the per-reference bins in first-appearance (insertion)
// order, reproduces htslib's khash bucket layout, runs compress_binning for
// the given hierarchy depth (nLvls) and bin count (nBins, the number that
// separates regular bins from the meta pseudo-bin), and returns the bins in
// the exact order htslib's idx_save_core would serialise them.
//
// metaBin is the pseudo-bin number (nBins+1 for BAI/CSI); it participates in
// the hash layout but is never merged or chunk-coalesced by
// compress_binning, matching htslib's `kh_key(bidx, k) >= idx->n_bins`
// guard.
func orderBins(inserted []binEntry, nLvls int, nBins, metaBin uint32) []binEntry {
	h := &binKhash{}
	for i := range inserted {
		e := inserted[i]
		k, isNew := h.put(e.ID)
		if isNew {
			ent := e // copy
			h.vals[k] = &ent
		} else {
			// htslib appends chunks to an existing bin; with our
			// per-bin-aggregated input this only happens if a bin is handed
			// in twice, so append the chunks.
			h.vals[k].Chunks = append(h.vals[k].Chunks, e.Chunks...)
		}
	}
	compressBinning(h, nLvls, nBins)

	out := make([]binEntry, 0, h.size)
	for k := uint32(0); k < h.nBuckets; k++ {
		if khFlagIsEither(h.flags, k) {
			continue
		}
		out = append(out, *h.vals[k])
	}
	return out
}

// compressBinning ports htslib's compress_binning: for each level from the
// deepest up to 1, merge a bin into its parent when the bin's chunk list
// spans fewer than htsMinMarkerDist compressed BGZF bytes, then coalesce
// adjacent chunks within every surviving bin. The meta pseudo-bin (key >=
// nBins) is never touched.
func compressBinning(h *binKhash, nLvls int, nBins uint32) {
	for l := nLvls; l > 0; l-- {
		start := htsBinFirst(l)
		for k := uint32(0); k < h.nBuckets; k++ {
			if khFlagIsEither(h.flags, k) {
				continue
			}
			key := h.keys[k]
			if key >= nBins || key < start {
				continue
			}
			p := h.vals[k]
			if l < nLvls && len(p.Chunks) > 1 {
				sortChunks(p.Chunks)
			}
			if len(p.Chunks) == 0 {
				continue
			}
			span := (uint64(p.Chunks[len(p.Chunks)-1].End) >> 16) - (uint64(p.Chunks[0].Beg) >> 16)
			if span < htsMinMarkerDist {
				kp, ok := h.get(htsBinParent(key))
				if !ok {
					continue
				}
				q := h.vals[kp]
				q.Chunks = append(q.Chunks, p.Chunks...)
				h.del(k)
			}
		}
	}
	// htslib sorts bin 0's chunks after the level loop.
	if k, ok := h.get(0); ok {
		sortChunks(h.vals[k].Chunks)
	}
	// Coalesce adjacent chunks that begin in the same BGZF block, for every
	// regular bin (key < nBins). Mirrors the second loop of compress_binning.
	for k := uint32(0); k < h.nBuckets; k++ {
		if khFlagIsEither(h.flags, k) {
			continue
		}
		if h.keys[k] >= nBins {
			continue
		}
		p := h.vals[k]
		if len(p.Chunks) == 0 {
			continue
		}
		m := 0
		for li := 1; li < len(p.Chunks); li++ {
			if uint64(p.Chunks[m].End)>>16 >= uint64(p.Chunks[li].Beg)>>16 {
				if p.Chunks[m].End < p.Chunks[li].End {
					p.Chunks[m].End = p.Chunks[li].End
				}
			} else {
				m++
				p.Chunks[m] = p.Chunks[li]
			}
		}
		p.Chunks = p.Chunks[:m+1]
	}
}

// sortChunks sorts a chunk list by Beg ascending, matching htslib's
// ks_introsort_off comparator pair64_lt = (a.u < b.u), which orders on the
// chunk start (Beg) only.
func sortChunks(c []binChunk) {
	sort.SliceStable(c, func(i, j int) bool { return c[i].Beg < c[j].Beg })
}
