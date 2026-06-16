package bedcoverage

import (
	"reflect"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// These tests exercise the pure -split-over-blocked-query helpers (queryBlocks,
// splitDepth, coveredFromDepths) with synthetic records. They require NO
// upstream binary and NO populated submodule, mirroring the binary-free layer
// the other bed* tools keep alongside their live-parity tests.

// blockedQuery builds a BED12 query record [start,end) carrying the given
// blocks, expressed as (relativeStart, size) pairs.
func blockedQuery(chrom string, start, end int, blocks ...[2]int) *bed.Record {
	r := &bed.Record{
		Chrom:      chrom,
		ChromStart: start,
		ChromEnd:   end,
		BlockCount: len(blocks),
	}
	for _, blk := range blocks {
		r.BlockStarts = append(r.BlockStarts, blk[0])
		r.BlockSizes = append(r.BlockSizes, blk[1])
	}
	return r
}

func plainB(chrom string, start, end int) *bed.Record {
	return &bed.Record{Chrom: chrom, ChromStart: start, ChromEnd: end}
}

func TestUnitQueryBlocks(t *testing.T) {
	// Two-block record: chr1 100-300 with blocks 100-150 and 250-300.
	q := blockedQuery("chr1", 100, 300, [2]int{0, 50}, [2]int{150, 50})
	got := queryBlocks(q)
	want := []block{{100, 150}, {250, 300}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("queryBlocks = %v, want %v", got, want)
	}

	// Unblocked record yields a single whole-span block.
	plain := plainB("chr1", 10, 40)
	if got := queryBlocks(plain); !reflect.DeepEqual(got, []block{{10, 40}}) {
		t.Fatalf("queryBlocks(unblocked) = %v, want [{10 40}]", got)
	}
}

func TestUnitSplitDepthExcludesIntrons(t *testing.T) {
	// Query chr1 100-300, blocks 100-150 and 250-300; intron 150-250.
	q := blockedQuery("chr1", 100, 300, [2]int{0, 50}, [2]int{150, 50})
	bs := []*bed.Record{
		plainB("chr1", 100, 150), // covers block1 fully (50 bp)
		plainB("chr1", 160, 240), // intron-only: MUST NOT count
		plainB("chr1", 270, 300), // covers 30 bp of block2
	}
	hitCount, depths := splitDepth(q, bs)

	// Depth vector spans the full record (200 bases), introns stay 0.
	if len(depths) != 200 {
		t.Fatalf("depth length = %d, want 200 (full span)", len(depths))
	}
	// Intron-only B excluded => only b1 and b3 contribute (one block each).
	if hitCount != 2 {
		t.Fatalf("hitCount = %d, want 2 (intron-only B excluded)", hitCount)
	}
	// Covered = 50 (block1) + 30 (block2 partial) = 80.
	if cov := coveredFromDepths(depths); cov != 80 {
		t.Fatalf("covered = %d, want 80", cov)
	}
	// Spot-check: base 0 (abs 100) depth 1, base 60 (abs 160, intron) depth 0,
	// base 170 (abs 270) depth 1, base 149 (abs 249, intron end) depth 0.
	for i, want := range map[int]int{0: 1, 49: 1, 60: 0, 149: 0, 170: 1, 199: 1} {
		if depths[i] != want {
			t.Fatalf("depths[%d] = %d, want %d", i, depths[i], want)
		}
	}
}

func TestUnitSplitDepthBStraddlingIntronCountedPerBlock(t *testing.T) {
	// A single B straddling the intron and overlapping BOTH query blocks is
	// counted once per block it touches (upstream findBlockedOverlaps pushes one
	// overlap sub-interval per query-block x hit-block pair).
	q := blockedQuery("chr1", 100, 300, [2]int{0, 50}, [2]int{150, 50})
	bs := []*bed.Record{plainB("chr1", 120, 280)} // 120-150 (30bp) + 250-280 (30bp)
	hitCount, depths := splitDepth(q, bs)
	if hitCount != 2 {
		t.Fatalf("hitCount = %d, want 2 (counted per overlapping block)", hitCount)
	}
	if cov := coveredFromDepths(depths); cov != 60 {
		t.Fatalf("covered = %d, want 60", cov)
	}
}

func TestUnitSplitDepthNoBlockOverlap(t *testing.T) {
	// B falls entirely in the intron: no hits, no covered bases, full-span len.
	q := blockedQuery("chr1", 100, 300, [2]int{0, 50}, [2]int{150, 50})
	bs := []*bed.Record{plainB("chr1", 170, 230)}
	hitCount, depths := splitDepth(q, bs)
	if hitCount != 0 {
		t.Fatalf("hitCount = %d, want 0", hitCount)
	}
	if cov := coveredFromDepths(depths); cov != 0 {
		t.Fatalf("covered = %d, want 0", cov)
	}
	if len(depths) != 200 {
		t.Fatalf("depth length = %d, want 200", len(depths))
	}
}

func TestUnitCoveredFromDepths(t *testing.T) {
	if got := coveredFromDepths([]int{0, 1, 2, 0, 3}); got != 3 {
		t.Fatalf("coveredFromDepths = %d, want 3", got)
	}
	if got := coveredFromDepths(nil); got != 0 {
		t.Fatalf("coveredFromDepths(nil) = %d, want 0", got)
	}
}
