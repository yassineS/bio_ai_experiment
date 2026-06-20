package bam

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

// vOff builds a virtual offset at compressed byte offset coff and in-block
// offset uoff. Spacing two records in the same bin by >= HTS_MIN_MARKER_DIST
// (0x10000) compressed bytes keeps the bin's chunk span above the threshold
// so compress_binning does not fold the bin into its parent. The bin span is
// measured on the compressed offset (the high 48 bits, i.e. value>>16).
func vOff(coff, uoff uint64) uint64 { return coff<<16 | (uoff & 0xffff) }

// TestCSILoffsetTrimsForOverlappingLargeRecord is the FIX 1 regression
// test. A large record assigned to a coarse parent bin overlaps the left
// edge of a finer bin that also holds a later, small record. The finer
// bin's loffset must reflect the large record's earlier virtual offset
// (the spec-defined linear-index value covering the bin's first locus),
// not the min vBeg of records merely assigned to that bin — otherwise a
// region query against the finer bin would wrongly trim the large
// record's chunk.
//
// The finest bin is given two records spanning well over a full BGZF block
// (compressed-offset delta >= 0x10000) so htslib's compress_binning leaves it
// in place rather than merging it into its parent — that is the regime in
// which the loffset trimming this test exercises matters.
func TestCSILoffsetTrimsForOverlappingLargeRecord(t *testing.T) {
	b := NewCSIBuilder(1, DefaultCSIMinShift, DefaultCSIDepth)
	tile := int64(1) << DefaultCSIMinShift // finest tile width (16 Kbp)
	const span = uint64(0x20000)           // > HTS_MIN_MARKER_DIST (0x10000)

	// Record A: large, spans finest tiles 0 and 1, so Reg2bin assigns it
	// to a coarse parent bin (a level-4 bin). Early virtual offset. A second
	// record A2 in the same coarse bin sits >0x10000 compressed bytes later
	// so the coarse bin survives compress_binning.
	vBegA, vEndA := vOff(0, 0), vOff(0, 100)
	if err := b.AddRecord(0, 0, 2*tile, vBegA, vEndA, true); err != nil {
		t.Fatalf("AddRecord A: %v", err)
	}
	vBegA2, vEndA2 := vOff(span, 0), vOff(span, 100)
	if err := b.AddRecord(0, 0, 2*tile, vBegA2, vEndA2, true); err != nil {
		t.Fatalf("AddRecord A2: %v", err)
	}
	// Records B/C: small, fully inside finest tile 1, so Reg2bin assigns
	// them to the finest-level bin for tile 1. They sit far enough apart in
	// compressed bytes that the bin's chunk span exceeds HTS_MIN_MARKER_DIST
	// and survives compression.
	vBegB, vEndB := vOff(span, 0), vOff(span, 100)
	if err := b.AddRecord(0, tile+10, tile+15, vBegB, vEndB, true); err != nil {
		t.Fatalf("AddRecord B: %v", err)
	}
	vBegC, vEndC := vOff(3*span, 0), vOff(3*span, 100)
	if err := b.AddRecord(0, tile+20, tile+25, vBegC, vEndC, true); err != nil {
		t.Fatalf("AddRecord C: %v", err)
	}
	idx := b.Finish()

	csi := idx.CSI
	finestBin := csi.Reg2bin(tile+10, tile+15)
	parentBin := csi.Reg2bin(0, 2*tile)
	if finestBin == parentBin {
		t.Fatalf("test setup: records share bin %d; expected distinct bins", finestBin)
	}

	// The finest bin for tile 1 must carry loffset == vBegA: record A
	// overlaps tile 1's left edge even though it is assigned to the
	// coarser parentBin. The old min-of-assigned-vBeg rule would have
	// produced vBegB here.
	var gotFinest tabix.VOffset
	var haveFinest bool
	for _, bin := range csi.Refs[0].Bins {
		if bin.ID == finestBin {
			gotFinest, haveFinest = bin.LOffset, true
		}
	}
	if !haveFinest {
		t.Fatalf("finest bin %d missing from index", finestBin)
	}
	if gotFinest != tabix.VOffset(vBegA) {
		t.Errorf("finest bin loffset: got %d, want %d (large record vBeg)", gotFinest, vBegA)
	}

	// Spec invariant: for every bin, loffset must not exceed vBeg of any
	// record whose span overlaps that bin's region.
	metaBin := csi.BinLimit() + 1
	for _, bin := range csi.Refs[0].Bins {
		if bin.ID == metaBin {
			continue // meta pseudo-bin
		}
		tl := binLeftmostTile(bin.ID, csi.Depth)
		// Record A overlaps tiles 0 and 1; records B/C overlap tile 1.
		if tl == 0 || tl == 1 {
			if bin.LOffset > tabix.VOffset(vBegA) {
				t.Errorf("bin %d (leftmost tile %d) loffset %d > record A vBeg %d",
					bin.ID, tl, bin.LOffset, vBegA)
			}
		}
	}

	// A region query against tile 1 must keep record A's chunk: with the
	// correct loffset (vBegA) it is not trimmed.
	chunks := idx.RegionChunks(0, int(tile+10), int(tile+15))
	if len(chunks) == 0 {
		t.Fatal("RegionChunks returned nothing for the populated tile")
	}
	keptA := false
	for _, c := range chunks {
		if c.Beg <= vBegA && c.End >= vEndA {
			keptA = true
		}
	}
	if !keptA {
		t.Errorf("region query trimmed the large record's chunk [%d,%d); chunks=%v",
			vBegA, vEndA, chunks)
	}
}

// TestCSILoffsetMonotone confirms loffsets stay consistent with the
// per-tile linear index: for a run of records at strictly increasing
// virtual offsets and increasing positions, each finest bin's loffset is
// the linear value of its own tile and the sequence is non-decreasing.
func TestCSILoffsetMonotone(t *testing.T) {
	b := NewCSIBuilder(1, DefaultCSIMinShift, DefaultCSIDepth)
	tile := int64(1) << DefaultCSIMinShift

	type rec struct {
		pos  int64
		vBeg uint64
	}
	// Each finest tile gets its own record at an early virtual offset, plus
	// a far-later trailing record in a high compressed block so the bin's
	// chunk span exceeds HTS_MIN_MARKER_DIST and survives compress_binning
	// (otherwise htslib folds each single-record finest bin into its
	// parent). The early record carries the bin's loffset.
	recs := []rec{
		{0 * tile, vOff(1, 0)},
		{1 * tile, vOff(2, 0)},
		{2 * tile, vOff(3, 0)},
		{3 * tile, vOff(4, 0)},
	}
	for _, r := range recs {
		if err := b.AddRecord(0, r.pos, r.pos+5, r.vBeg, r.vBeg+50, true); err != nil {
			t.Fatalf("AddRecord pos %d: %v", r.pos, err)
		}
		// trailing record in the same tile, many blocks later
		far := vOff(uint64((r.pos/tile)*4+20), 0)
		if err := b.AddRecord(0, r.pos+6, r.pos+11, far, far+50, true); err != nil {
			t.Fatalf("AddRecord far pos %d: %v", r.pos, err)
		}
	}
	idx := b.Finish()
	csi := idx.CSI

	loff := map[uint32]tabix.VOffset{}
	for _, bin := range csi.Refs[0].Bins {
		loff[bin.ID] = bin.LOffset
	}
	var prev tabix.VOffset
	for i, r := range recs {
		binID := csi.Reg2bin(r.pos, r.pos+5)
		got, ok := loff[binID]
		if !ok {
			t.Fatalf("bin %d for record %d missing", binID, i)
		}
		if got != tabix.VOffset(r.vBeg) {
			t.Errorf("record %d bin %d loffset: got %d, want %d", i, binID, got, r.vBeg)
		}
		if i > 0 && got < prev {
			t.Errorf("loffset not monotone at record %d: %d < %d", i, got, prev)
		}
		prev = got
	}
}
