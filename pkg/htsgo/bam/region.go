// BAI-specific region helpers: UnionChunks aggregates per-region BAI
// chunks into a single sorted, merged list — the seek-and-scan entry
// point for region-query reads. Lives in pkg/htsgo/bam because it
// depends on BAIIndex/BAIChunk and is meaningless without them; the
// region package itself stays free of BAI references so non-BAM
// region consumers can use it cleanly.

package bam

import (
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
)

// UnionChunks combines the BAI region-chunks across every ResolvedRegion
// into a single sorted, merged slice. It is the entry point used by the
// seek-and-scan region-query path in samtools view.
func UnionChunks(idx *BAIIndex, regions []region.ResolvedRegion) []BAIChunk {
	var all []BAIChunk
	for _, r := range regions {
		all = append(all, idx.RegionChunks(r.RefID, r.Beg0, r.End0)...)
	}
	if len(all) == 0 {
		return nil
	}
	sort.Slice(all, func(a, b int) bool { return all[a].Beg < all[b].Beg })
	merged := all[:0]
	cur := all[0]
	for i := 1; i < len(all); i++ {
		if all[i].Beg <= cur.End {
			if all[i].End > cur.End {
				cur.End = all[i].End
			}
		} else {
			merged = append(merged, cur)
			cur = all[i]
		}
	}
	merged = append(merged, cur)
	return merged
}

// UnionChunksCSI is the CSI counterpart of UnionChunks: it combines the
// region-chunks of a BAM `.csi` index across every ResolvedRegion into a
// single sorted, merged slice. The result has the same BAIChunk shape so
// the samtools view seek-and-scan path is index-kind agnostic.
func UnionChunksCSI(idx *CSIIndex, regions []region.ResolvedRegion) []BAIChunk {
	var all []BAIChunk
	for _, r := range regions {
		all = append(all, idx.RegionChunks(r.RefID, r.Beg0, r.End0)...)
	}
	if len(all) == 0 {
		return nil
	}
	sort.Slice(all, func(a, b int) bool { return all[a].Beg < all[b].Beg })
	merged := all[:0]
	cur := all[0]
	for i := 1; i < len(all); i++ {
		if all[i].Beg <= cur.End {
			if all[i].End > cur.End {
				cur.End = all[i].End
			}
		} else {
			merged = append(merged, cur)
			cur = all[i]
		}
	}
	merged = append(merged, cur)
	return merged
}
