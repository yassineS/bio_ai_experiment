// UCSC-bin hit ordering for the default (non-sorted, non-sortout) intersect
// path. Upstream `bedtools intersect` indexes the B database in a UCSC binning
// tree and, for each A query, walks the bin levels finest-first, each level's
// bins in ascending number, and the records within a bin in insertion (file)
// order. The resulting hit order is NOT plain file order, so this file
// reproduces the bin placement so the emitted hits match upstream byte-for-byte.
package bedintersect

import "sort"

// Bin scheme constants, copied from upstream BinTree.h. NUM_BIN_LEVELS levels,
// the finest bin spanning 2^binFirstShift bases, each coarser level grouping
// 2^binNextShift finer bins.
const (
	numBinLevels  = 8
	binFirstShift = 14
	binNextShift  = 3
)

// binOffsetsExtended is the per-level starting bin number, copied verbatim from
// upstream BinTree.h's _binOffsetsExtended. Index 0 is the finest level.
var binOffsetsExtended = [numBinLevels]int64{
	262144 + 32678 + 4096 + 512 + 64 + 8 + 1,
	32678 + 4096 + 512 + 64 + 8 + 1,
	4096 + 512 + 64 + 8 + 1,
	512 + 64 + 8 + 1,
	64 + 8 + 1,
	8 + 1,
	1,
	0,
}

// binLevelAndNumber returns the UCSC bin level (0 = finest) and absolute bin
// number a [start,end) record is placed in, mirroring upstream BinTree::getBin.
// A record that does not fit any level (out of range) is placed past the last
// level so it still sorts deterministically after every real hit.
func binLevelAndNumber(start, end int) (level int, bin int64) {
	s := int64(start)
	e := int64(end) - 1
	if e < s {
		e = s
	}
	s >>= binFirstShift
	e >>= binFirstShift
	for i := 0; i < numBinLevels; i++ {
		if s == e {
			return i, binOffsetsExtended[i] + s
		}
		s >>= binNextShift
		e >>= binNextShift
	}
	return numBinLevels, 0
}

// orderHitsByBin reorders a query's hits into upstream's bin-traversal order:
// by bin level (finest first), then bin number ascending, then the record's
// insertion order within its chromosome (file order). It is used for every
// default-path output mode; the -sortout and -sorted paths use their own
// orderings and must not call this.
func orderHitsByBin(hits []rawHit) {
	if len(hits) < 2 {
		return
	}
	sort.SliceStable(hits, func(i, j int) bool {
		li, bi := binLevelAndNumber(hits[i].b.start, hits[i].b.end)
		lj, bj := binLevelAndNumber(hits[j].b.start, hits[j].b.end)
		if li != lj {
			return li < lj
		}
		if bi != bj {
			return bi < bj
		}
		return hits[i].b.order < hits[j].b.order
	})
}
