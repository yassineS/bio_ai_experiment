// BAM CSI (Coordinate-Sorted Index) support. CSI is the index format
// needed for reference sequences longer than ~512 Mbp — the BAI bin
// scheme caps coordinates at 2^29. A BAM `.csi` shares the BAI chunk /
// linear-offset model but parameterises the bin hierarchy by
// (min_shift, depth) and stores no linear index (CSI replaces it with a
// per-bin "loffset" minimum virtual offset).
//
// The binning math (reg2bin / reg2bins for a configurable min_shift and
// depth) is identical to the tabix `.csi` used for VCF/BCF, so this file
// reuses pkg/htsgo/tabix's CSI primitives rather than re-deriving them.
// The only BAM-specific additions are the pseudo-bin metadata that BAI
// also carries (mapped/unmapped counts and the first/last virtual
// offsets) and an empty aux block (a BAM `.csi` has l_aux == 0, whereas
// a tabix `.csi` carries a column-config aux payload).

package bam

import (
	"fmt"
	"io"
	"os"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

// DefaultCSIMinShift is the htslib default min_shift for a BAM `.csi`.
// Combined with DefaultCSIDepth it yields the BAI-equivalent six-level
// bin hierarchy (level 0..5).
const DefaultCSIMinShift = 14

// DefaultCSIDepth is the htslib default CSI depth.
const DefaultCSIDepth = 5

// CSIIndex is the in-memory form of a BAM `.csi` index. It wraps a
// tabix.CSI — the format container, binning math, writer and reader are
// shared — and records the (min_shift, depth) the index was built at.
//
// The pseudo-bin metadata (mapped/unmapped counts and first/last virtual
// offsets) is stored as a regular bin whose ID is MetaBinCSI(depth);
// region queries skip it because that bin number is never returned by
// reg2bins.
type CSIIndex struct {
	// CSI is the underlying tabix CSI container. Its Aux field is empty
	// for a BAM `.csi`.
	CSI *tabix.CSI
}

// MetaBinCSI returns the pseudo-bin number for a CSI of the given depth.
// htslib places the meta pseudo-bin one past the last regular bin: with
// n_bins regular bins numbered 0..n_bins-1, META_BIN(idx) == n_bins + 1.
// BinLimit() is the regular-bin count (n_bins), so the meta bin is
// BinLimit() + 1.
func MetaBinCSI(minShift, depth int32) uint32 {
	c := tabix.NewCSI(minShift, depth)
	return c.BinLimit() + 1
}

// csiRefBuilder is the in-progress per-reference state for the CSI
// builder. It mirrors baiRefBuilder but addresses the configurable CSI
// bin space.
//
// linear is the CSI analogue of the BAI linear index: one slot per
// finest-level tile (tile width 2^min_shift), each holding the smallest
// virtual offset of any record whose span overlaps that tile. It is not
// serialised — CSI has no on-disk linear index — but it is the correct
// source for each bin's loffset (the linear-index value covering the
// first locus the bin spans), including large records that overlap a
// bin's left edge yet are assigned to a coarser parent bin.
type csiRefBuilder struct {
	bins map[uint32]*tabix.CSIBin
	// order records bin IDs in first-appearance (insertion) order so Finish
	// can replay htslib's khash insertion sequence and reproduce its on-disk
	// bin ordering byte-for-byte.
	order     []uint32
	linear    []uint64
	firstVOff uint64
	lastVOff  uint64
	mapped    uint64
	unmapped  uint64
	haveOff   bool
}

// CSIBuilder collects bin/meta state for an in-progress BAM `.csi`.
// Callers add records sequentially in BAM order via AddRecord, then call
// Finish to obtain the populated CSIIndex.
type CSIBuilder struct {
	csi       *tabix.CSI
	refs      []*csiRefBuilder
	noCoor    uint64
	totalRefs int
}

// NewCSIBuilder creates an empty CSI builder over numRefs reference
// sequences. minShift and depth select the bin hierarchy; non-positive
// values fall back to DefaultCSIMinShift / DefaultCSIDepth.
func NewCSIBuilder(numRefs int, minShift, depth int32) *CSIBuilder {
	if minShift <= 0 {
		minShift = DefaultCSIMinShift
	}
	if depth <= 0 {
		depth = DefaultCSIDepth
	}
	b := &CSIBuilder{
		csi:       tabix.NewCSI(minShift, depth),
		totalRefs: numRefs,
		refs:      make([]*csiRefBuilder, numRefs),
	}
	for i := range b.refs {
		b.refs[i] = &csiRefBuilder{bins: map[uint32]*tabix.CSIBin{}, firstVOff: ^uint64(0)}
	}
	return b
}

// AddRecord registers one BAM record with the index. refID is the
// 0-based reference index (or -1 for unmapped); beg/end are the 0-based
// half-open reference coordinates the record spans. vBeg is the virtual
// offset of the record's first byte and vEnd the virtual offset just
// past it. Records with refID < 0 bump the trailing n_no_coor counter.
func (b *CSIBuilder) AddRecord(refID int, beg, end int64, vBeg, vEnd uint64, mapped bool) error {
	if refID < 0 {
		b.noCoor++
		return nil
	}
	if refID >= b.totalRefs {
		return fmt.Errorf("bam: CSI refID %d >= %d", refID, b.totalRefs)
	}
	if end <= beg {
		end = beg + 1
	}
	ref := b.refs[refID]
	if !ref.haveOff || vBeg < ref.firstVOff {
		ref.firstVOff = vBeg
		ref.haveOff = true
	}
	if vEnd > ref.lastVOff {
		ref.lastVOff = vEnd
	}
	if mapped {
		ref.mapped++
	} else {
		ref.unmapped++
	}

	binID := b.csi.Reg2bin(beg, end)
	bin, ok := ref.bins[binID]
	if !ok {
		bin = &tabix.CSIBin{ID: binID}
		ref.bins[binID] = bin
		ref.order = append(ref.order, binID)
	}
	if n := len(bin.Chunks); n > 0 && uint64(bin.Chunks[n-1].End) >= vBeg {
		if tabix.VOffset(vEnd) > bin.Chunks[n-1].End {
			bin.Chunks[n-1].End = tabix.VOffset(vEnd)
		}
	} else {
		bin.Chunks = append(bin.Chunks, tabix.CSIChunk{Beg: tabix.VOffset(vBeg), End: tabix.VOffset(vEnd)})
	}

	// Linear-index update — every finest-level tile (tile width
	// 2^min_shift) the record's [beg, end) span touches records the
	// smallest virtual offset of any record overlapping it. The bin
	// loffsets are derived from this array in Finish, so a large record
	// assigned to a coarse parent bin still lowers the loffset of every
	// finer bin it overlaps.
	shift := uint32(b.csi.MinShift)
	first := beg >> shift
	last := (end - 1) >> shift
	if last < first {
		last = first
	}
	for int64(len(ref.linear)) <= last {
		ref.linear = append(ref.linear, ^uint64(0))
	}
	for t := first; t <= last; t++ {
		if ref.linear[t] == ^uint64(0) || vBeg < ref.linear[t] {
			ref.linear[t] = vBeg
		}
	}
	return nil
}

// binLeftmostTile returns the index of the leftmost finest-level CSI
// tile (tile width 2^min_shift) covered by the bin numbered binID, for a
// CSI of the given depth. It is the inverse of the CSI bin-hierarchy
// formula: a bin at level L has level-relative index i = binID - t_L
// (where t_L = (8^L-1)/7 is the count of bins at all shallower levels),
// spans 2^(3*(depth-L)) finest tiles, and its first such tile is
// i << (3*(depth-L)).
func binLeftmostTile(binID uint32, depth int32) int64 {
	var t uint32
	for l := int32(0); l <= depth; l++ {
		next := t + (1 << uint32(3*l))
		if binID < next {
			rel := int64(binID - t)
			return rel << uint32(3*(depth-l))
		}
		t = next
	}
	// binID is the meta pseudo-bin or out of range; its loffset is unused.
	return 0
}

// loffsetForBin returns the loffset for binID: the linear-index value of
// the first (leftmost) finest-level tile the bin covers. Per the CSI
// spec a bin's loffset is the linear-index value covering the first
// locus the bin spans, so a query trims chunks against the smallest
// virtual offset of any record overlapping that locus — including large
// records assigned to a coarser parent bin.
func loffsetForBin(binID uint32, depth int32, linear []uint64) tabix.VOffset {
	tile := binLeftmostTile(binID, depth)
	if tile < 0 || tile >= int64(len(linear)) {
		return 0
	}
	v := linear[tile]
	if v == ^uint64(0) {
		return 0
	}
	return tabix.VOffset(v)
}

// Finish assembles the CSIIndex. It injects the BAM pseudo-bin per
// reference (mapped/unmapped counts plus first/last virtual offsets),
// matching the metadata BAI carries, and returns the populated index.
//
// Each regular bin's LOffset is set to the linear-index value of the
// leftmost finest-level tile the bin covers (see loffsetForBin), which
// is the spec-defined loffset and never exceeds vBeg of any record whose
// span overlaps the bin's region.
func (b *CSIBuilder) Finish() *CSIIndex {
	nBins := b.csi.BinLimit()
	metaBin := nBins + 1
	b.csi.NoCoor = b.noCoor
	b.csi.Refs = make([]tabix.CSIRef, b.totalRefs)
	for i, ref := range b.refs {
		// Assemble bins in first-appearance (insertion) order, then append
		// the meta pseudo-bin last, mirroring htslib's khash insertion
		// order. OrderBins reproduces htslib's khash iteration order and
		// compress_binning so the serialised .csi is byte-identical to
		// `samtools index -c`. Each surviving bin's loffset is the
		// linear-index value of its leftmost finest-level tile (see
		// loffsetForBin), which compress_binning preserves on the parent
		// when a small child bin is folded in.
		inserted := make([]tabix.BinEntry, 0, len(ref.order)+1)
		for _, id := range ref.order {
			bb := ref.bins[id]
			chunks := make([]tabix.BinChunk, len(bb.Chunks))
			for j, c := range bb.Chunks {
				chunks[j] = tabix.BinChunk{Beg: c.Beg, End: c.End}
			}
			inserted = append(inserted, tabix.BinEntry{
				ID:      id,
				LOffset: loffsetForBin(id, b.csi.Depth, ref.linear),
				Chunks:  chunks,
			})
		}
		if ref.haveOff || ref.mapped > 0 || ref.unmapped > 0 {
			// The pseudo-bin LOffset is 0: htslib's update_loff sets loff = 0
			// for every bin whose key is >= n_bins (the meta pseudo-bin). Its
			// two chunks carry (firstVOff, lastVOff) and (mapped, unmapped)
			// exactly as the BAI meta pseudo-bin does.
			inserted = append(inserted, tabix.BinEntry{
				ID:      metaBin,
				LOffset: 0,
				Chunks: []tabix.BinChunk{
					{Beg: tabix.VOffset(ref.firstVOff), End: tabix.VOffset(ref.lastVOff)},
					{Beg: tabix.VOffset(ref.mapped), End: tabix.VOffset(ref.unmapped)},
				},
			})
		}
		ordered := tabix.OrderBins(inserted, int(b.csi.Depth), nBins, metaBin)
		bins := make([]tabix.CSIBin, len(ordered))
		for j, e := range ordered {
			chunks := make([]tabix.CSIChunk, len(e.Chunks))
			for k, c := range e.Chunks {
				chunks[k] = tabix.CSIChunk{Beg: c.Beg, End: c.End}
			}
			bins[j] = tabix.CSIBin{ID: e.ID, LOffset: e.LOffset, Chunks: chunks}
		}
		b.csi.Refs[i] = tabix.CSIRef{Bins: bins}
	}
	return &CSIIndex{CSI: b.csi}
}

// BuildCSI streams every record from a BAMReader and returns the
// assembled BAM `.csi` index. It does not validate sort order — callers
// must pass a coordinate-sorted BAM. minShift and depth parameterise the
// bin hierarchy; pass 0 for the htslib defaults.
func BuildCSI(br *sam.BAMReader, numRefs int, minShift, depth int32) (*CSIIndex, error) {
	bld := NewCSIBuilder(numRefs, minShift, depth)
	for {
		vBeg := br.VirtualOffset()
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		vEnd := br.VirtualOffset()

		refID := -1
		if rec.RName != "" && rec.RName != "*" {
			refID = br.Header().RefIndex(rec.RName)
			if refID < 0 {
				return nil, fmt.Errorf("bam: record references unknown @SQ %q", rec.RName)
			}
		}
		mapped := !rec.IsUnmapped()
		beg := int64(rec.Pos) - 1
		if beg < 0 {
			beg = 0
		}
		end := beg + int64(rec.Cigar.ReferenceLength())
		if err := bld.AddRecord(refID, beg, end, vBeg, vEnd, mapped); err != nil {
			return nil, err
		}
	}
	return bld.Finish(), nil
}

// WriteCSI serialises idx to w in the canonical BAM `.csi` byte layout
// and BGZF-compresses it. A BAM `.csi` has an empty aux block.
func WriteCSI(w io.Writer, idx *CSIIndex) error {
	bw := bgzip.NewWriter(w)
	if err := idx.CSI.Write(bw); err != nil {
		return err
	}
	return bw.Close()
}

// ReadCSI parses a BGZF-compressed BAM `.csi` stream into a CSIIndex.
func ReadCSI(r io.Reader) (*CSIIndex, error) {
	br, err := bgzip.NewReader(r)
	if err != nil {
		return nil, err
	}
	defer br.Close()
	csi, err := tabix.ReadCSI(br)
	if err != nil {
		return nil, err
	}
	return &CSIIndex{CSI: csi}, nil
}

// ReadCSIFile opens the BAM `.csi` index at path and parses it into a
// CSIIndex. It is the on-disk counterpart of ReadCSI.
func ReadCSIFile(path string) (*CSIIndex, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ReadCSI(f)
}

// MetaCounts returns the (mapped, unmapped) counts stored in the CSI
// pseudo-bin for the given reference, or (0, 0, false) when absent.
func (idx *CSIIndex) MetaCounts(refID int) (mapped, unmapped uint64, ok bool) {
	metaBin := idx.CSI.BinLimit() + 1
	if refID < 0 || refID >= len(idx.CSI.Refs) {
		return 0, 0, false
	}
	for _, bin := range idx.CSI.Refs[refID].Bins {
		if bin.ID == metaBin && len(bin.Chunks) >= 2 {
			return uint64(bin.Chunks[1].Beg), uint64(bin.Chunks[1].End), true
		}
	}
	return 0, 0, false
}

// MetaBounds returns the (firstVOff, lastVOff) virtual offsets recorded
// in the CSI pseudo-bin for the given reference, or (0, 0, false) when
// absent.
func (idx *CSIIndex) MetaBounds(refID int) (first, last uint64, ok bool) {
	metaBin := idx.CSI.BinLimit() + 1
	if refID < 0 || refID >= len(idx.CSI.Refs) {
		return 0, 0, false
	}
	for _, bin := range idx.CSI.Refs[refID].Bins {
		if bin.ID == metaBin && len(bin.Chunks) >= 1 {
			return uint64(bin.Chunks[0].Beg), uint64(bin.Chunks[0].End), true
		}
	}
	return 0, 0, false
}

// RegionChunks returns the merged BAIChunk list that may contain records
// overlapping [beg, end) on the given reference. beg and end are 0-based
// half-open. The result has the same shape as BAIIndex.RegionChunks so
// the samtools view seek-and-scan path can consume either index kind.
func (idx *CSIIndex) RegionChunks(refID, beg, end int) []BAIChunk {
	cc := idx.CSI.RegionChunks(refID, int64(beg), int64(end))
	if len(cc) == 0 {
		return nil
	}
	out := make([]BAIChunk, len(cc))
	for i, c := range cc {
		out[i] = BAIChunk{Beg: uint64(c.Beg), End: uint64(c.End)}
	}
	return out
}

// MaxPos returns the largest 0-based position addressable by this index
// at its (min_shift, depth). For the default CSI parameters this exceeds
// the BAI 2^29 ceiling, which is the whole point of CSI.
func (idx *CSIIndex) MaxPos() int64 { return idx.CSI.MaxPos() }
