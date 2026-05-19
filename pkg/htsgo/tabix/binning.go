package tabix

// The UCSC binning scheme as used by SAM/BAM, BAI, and TBI.
//
// The 1-based genome is partitioned into a fixed five-level hierarchy of
// fixed-size bins:
//
//	level 0:           bin 0 covers [0, 2^29)              — 1 bin
//	level 1: bins 1..8       cover 2^26 bp each             — 8 bins
//	level 2: bins 9..72      cover 2^23 bp each             — 64 bins
//	level 3: bins 73..584    cover 2^20 bp each             — 512 bins
//	level 4: bins 585..4680  cover 2^17 bp each             — 4 096 bins
//	level 5: bins 4681..37448 cover 2^14 bp each            — 32 768 bins
//
// A record is placed in the smallest bin that fully contains its [beg, end)
// half-open interval. A region query consults all bins at every level whose
// span overlaps the region.

// MaxBin is one past the largest legal bin number — bins are in [0, MaxBin).
// 4681 + 32768 = 37449.
const MaxBin = 37449

// TileSize is the linear index granularity: one entry per 16,384 bp.
const TileSize = 1 << 14

// Reg2bin returns the smallest bin number that fully contains the half-open
// interval [beg, end). beg and end are 0-based.
//
// The math follows the htslib paper (Li 2011, Bioinformatics 27): convert
// half-open to closed by decrementing end, then walk from the deepest level
// outward, returning the first level where beg and end share a tile.
func Reg2bin(beg, end int) int {
	end--
	if beg>>14 == end>>14 {
		return ((1<<15)-1)/7 + (beg >> 14)
	}
	if beg>>17 == end>>17 {
		return ((1<<12)-1)/7 + (beg >> 17)
	}
	if beg>>20 == end>>20 {
		return ((1<<9)-1)/7 + (beg >> 20)
	}
	if beg>>23 == end>>23 {
		return ((1<<6)-1)/7 + (beg >> 23)
	}
	if beg>>26 == end>>26 {
		return ((1<<3)-1)/7 + (beg >> 26)
	}
	return 0
}

// Reg2bins returns every bin number whose span overlaps the half-open
// interval [beg, end). The returned slice is in ascending bin order and
// always starts with bin 0.
func Reg2bins(beg, end int) []int {
	if end <= beg {
		return []int{0}
	}
	end--
	out := make([]int, 0, 16)
	out = append(out, 0)
	for k := 1 + (beg >> 26); k <= 1+(end>>26); k++ {
		out = append(out, k)
	}
	for k := 9 + (beg >> 23); k <= 9+(end>>23); k++ {
		out = append(out, k)
	}
	for k := 73 + (beg >> 20); k <= 73+(end>>20); k++ {
		out = append(out, k)
	}
	for k := 585 + (beg >> 17); k <= 585+(end>>17); k++ {
		out = append(out, k)
	}
	for k := 4681 + (beg >> 14); k <= 4681+(end>>14); k++ {
		out = append(out, k)
	}
	return out
}

// LinearTile returns the linear-index tile that contains 0-based position
// pos. Each tile covers TileSize bp.
func LinearTile(pos int) int {
	if pos < 0 {
		return 0
	}
	return pos >> 14
}
