package bam

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

// TestCSILoffsetTrimsForOverlappingLargeRecord is the FIX 1 regression
// test. A large record assigned to a coarse parent bin overlaps the left
// edge of a finer bin that also holds a later, small record. The finer
// bin's loffset must reflect the large record's earlier virtual offset
// (the spec-defined linear-index value covering the bin's first locus),
// not the min vBeg of records merely assigned to that bin — otherwise a
// region query against the finer bin would wrongly trim the large
// record's chunk.
func TestCSILoffsetTrimsForOverlappingLargeRecord(t *testing.T) {
	b := NewCSIBuilder(1, DefaultCSIMinShift, DefaultCSIDepth)
	tile := int64(1) << DefaultCSIMinShift // finest tile width (16 Kbp)

	// Record A: large, spans finest tiles 0 and 1, so Reg2bin assigns it
	// to a coarse parent bin (a level-4 bin). Early virtual offset.
	const vBegA, vEndA = 100, 200
	if err := b.AddRecord(0, 0, 2*tile, vBegA, vEndA, true); err != nil {
		t.Fatalf("AddRecord A: %v", err)
	}
	// Record B: small, fully inside finest tile 1, so Reg2bin assigns it
	// to the finest-level bin for tile 1. Later virtual offset.
	const vBegB, vEndB = 5000, 5100
	if err := b.AddRecord(0, tile+10, tile+15, vBegB, vEndB, true); err != nil {
		t.Fatalf("AddRecord B: %v", err)
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
	for _, bin := range csi.Refs[0].Bins {
		if bin.ID == csi.BinLimit() {
			continue // meta pseudo-bin
		}
		tl := binLeftmostTile(bin.ID, csi.Depth)
		// Record A overlaps tiles 0 and 1; record B overlaps tile 1.
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
	recs := []rec{
		{0 * tile, 100},
		{1 * tile, 500},
		{2 * tile, 900},
		{3 * tile, 1300},
	}
	for _, r := range recs {
		if err := b.AddRecord(0, r.pos, r.pos+5, r.vBeg, r.vBeg+50, true); err != nil {
			t.Fatalf("AddRecord pos %d: %v", r.pos, err)
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
