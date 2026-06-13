package bedintersect

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestReadInRecords_CRAM verifies that a CRAM stream on -a/-b decodes to the
// same intervals as the equivalent BAM: bedtools only needs the alignment
// coordinates, which CRAM stores directly (no decode reference required).
func TestReadInRecords_CRAM(t *testing.T) {
	// Build the header from SAM text so the @SQ lines (which the CRAM writer
	// serialises and the reader needs to resolve reference ids) are present.
	hr, err := sam.NewSAMReader(strings.NewReader("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n"))
	if err != nil {
		t.Fatalf("NewSAMReader(header): %v", err)
	}
	hdr := hr.Header()
	mk := func(qname string, pos int32, cig string) *sam.Record {
		c, err := sam.ParseCigar(cig)
		if err != nil {
			t.Fatalf("ParseCigar(%q): %v", cig, err)
		}
		seqLen := c.QueryLength()
		return &sam.Record{
			QName: qname, Flag: 0, RName: "chr1", Pos: pos, MapQ: 60,
			Cigar: c, Seq: string(bytes.Repeat([]byte("A"), seqLen)),
			Qual: bytes.Repeat([]byte{30}, seqLen), RNext: "*", PNext: 0,
		}
	}
	recs := []*sam.Record{
		mk("r1", 100, "50M"),
		mk("r2", 200, "20M30N20M"), // spliced -> two BED12 blocks
		mk("r3", 500, "10M2D10M"),
	}
	// Build the CRAM bytes (reference-free writer).
	var cramBuf bytes.Buffer
	if err := cram.WriteCRAM(&cramBuf, hdr, recs); err != nil {
		t.Fatalf("WriteCRAM: %v", err)
	}

	got, err := readInRecords(bytes.NewReader(cramBuf.Bytes()))
	if err != nil {
		t.Fatalf("readInRecords(CRAM): %v", err)
	}
	// Oracle: run the same records through the BAM->BED12 converter directly.
	var want []*inRecord
	for _, r := range recs {
		want = append(want, bamToBED12(r))
	}
	if len(got) != len(want) {
		t.Fatalf("CRAM intervals: got %d, want %d", len(got), len(want))
	}
	for i := range want {
		if got[i].chrom != want[i].chrom || got[i].start != want[i].start || got[i].end != want[i].end {
			t.Errorf("record %d: got %s:%d-%d, want %s:%d-%d", i,
				got[i].chrom, got[i].start, got[i].end,
				want[i].chrom, want[i].start, want[i].end)
		}
	}
}
