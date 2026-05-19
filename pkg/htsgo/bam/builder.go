// BAI builder helpers that consume a BAMReader and emit a *BAIIndex.
// These live in the bam package (alongside the BAI format primitives in
// bai.go) so the BAM reader and the BAI builder share a single import
// edge — the htslib precedent for the same coupling. The runtime
// `samtools index` orchestration (option parsing, file I/O) stays in
// `tools/samtools` because it's CLI-level, not format-level.

package bam

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// BuildBAI streams every record from a BAMReader and returns the assembled
// BAI index. It does not validate sort order — callers must pass a
// coordinate-sorted BAM (where coordinates means (refID, 0-based pos)).
func BuildBAI(br *sam.BAMReader, numRefs int) (*BAIIndex, error) {
	bld := NewBAIBuilder(numRefs)
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

		// Resolve refID from RName via header order.
		refID := -1
		if rec.RName != "" && rec.RName != "*" {
			refID = br.Header().RefIndex(rec.RName)
			if refID < 0 {
				return nil, fmt.Errorf("bam: record references unknown @SQ %q", rec.RName)
			}
		}
		mapped := !rec.IsUnmapped()
		// Unmapped records that nevertheless carry a refID + pos still go
		// into the regular bin/linear index per the SAM spec: htslib treats
		// them as "placed but unmapped". Records with refID == -1 are the
		// truly unplaced ones that bump n_no_coor.
		beg := int(rec.Pos) - 1
		if beg < 0 {
			beg = 0
		}
		end := beg + rec.Cigar.ReferenceLength()
		if err := bld.AddRecord(refID, beg, end, vBeg, vEnd, mapped); err != nil {
			return nil, err
		}
	}
	return bld.Finish(), nil
}
