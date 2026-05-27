package bamtobed

import (
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// DecodeCRAMToBED reads a CRAM stream from r, using refFASTA for the
// reference-compressed base reconstruction, and returns (bedText, refs,
// nil). Output format mirrors DecodeBAMToBED: one BED6 row per primary
// alignment, refs is the @SQ list from the embedded SAM header. The
// reference path may be empty for reference-free CRAM streams (e.g. the
// empty CRAM fixture has no slices and never needs the reference).
func DecodeCRAMToBED(r io.Reader, refFASTA string) ([]byte, []BAMRef, error) {
	rr, err := cram.NewRecordReader(r)
	if err != nil {
		return nil, nil, fmt.Errorf("bamtobed: open CRAM: %w", err)
	}
	defer rr.Close()
	if refFASTA != "" {
		if err := rr.SetReferenceFASTA(refFASTA); err != nil {
			return nil, nil, fmt.Errorf("bamtobed: open CRAM reference: %w", err)
		}
	}
	hdr := rr.Header()
	var refs []BAMRef
	if hdr != nil {
		refs = make([]BAMRef, len(hdr.Refs))
		for i, ref := range hdr.Refs {
			refs[i] = BAMRef{Name: ref.Name, Length: int(ref.Length)}
		}
	}
	var out strings.Builder
	for {
		rec, err := rr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, nil, fmt.Errorf("bamtobed: read CRAM: %w", err)
		}
		if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
			rec.IsDuplicate() || rec.IsQCFail() {
			continue
		}
		refLen := rec.Cigar.ReferenceLength()
		if refLen <= 0 {
			continue
		}
		start := int(rec.Pos) - 1
		if start < 0 {
			continue
		}
		end := start + refLen
		strand := "+"
		if rec.Flag&sam.FlagReverse != 0 {
			strand = "-"
		}
		name := rec.QName
		if name == "" {
			name = "."
		}
		fmt.Fprintf(&out, "%s\t%d\t%d\t%s\t%d\t%s\n",
			rec.RName, start, end, name, rec.MapQ, strand)
	}
	return []byte(out.String()), refs, nil
}
