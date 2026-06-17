// UCSC-bin hit ordering for `bedtools window`. Upstream loads the B database
// into a UCSC binning tree (loadBedFileIntoMap) and, for each A query, walks the
// bin levels finest-first, each level's bins in ascending number, and the
// records within a bin in insertion (file) order (BedFile::allHits). The
// resulting per-A hit order is therefore NOT plain file order; this file
// reproduces the bin placement so the emitted hits match upstream byte-for-byte.
//
// The bin of a B record is computed from its ORIGINAL (unexpanded) coordinates,
// because upstream bins B by its loaded coordinates and only the A query window
// is fudged.
package bedwindow

import "sort"

// Bin scheme constants, copied from upstream BinTree.h. numBinLevels levels, the
// finest bin spanning 2^binFirstShift bases, each coarser level grouping
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
// number a [start, end) record is placed in, mirroring upstream
// BinTree::getBin. A record that does not fit any level is placed past the last
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

// orderHitsByBin reorders a query's hits into upstream's bin-traversal order: by
// bin level (finest first), then bin number ascending, then the record's
// in-chromosome insertion (file) order.
func orderHitsByBin(hits []*rec) {
	if len(hits) < 2 {
		return
	}
	sort.SliceStable(hits, func(i, j int) bool {
		li, bi := binLevelAndNumber(hits[i].start, hits[i].end)
		lj, bj := binLevelAndNumber(hits[j].start, hits[j].end)
		if li != lj {
			return li < lj
		}
		if bi != bj {
			return bi < bj
		}
		return hits[i].order < hits[j].order
	})
}
