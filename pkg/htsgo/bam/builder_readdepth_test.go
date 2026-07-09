package bam

import (
	"bytes"
	"fmt"
	"io"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// buildBAIFull rebuilds the BAI using the historical full-record decode path
// (br.Read, which materialises SEQ/QUAL/QNAME/aux) rather than the
// allocation-free ReadDepthInto path BuildBAI now uses. Keeping this reference
// implementation in the test lets us assert byte-for-byte that the perf change
// did not alter the emitted index.
func buildBAIFull(br *sam.BAMReader, numRefs int) (*BAIIndex, error) {
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
		refID := -1
		if rec.RName != "" && rec.RName != "*" {
			refID = br.Header().RefIndex(rec.RName)
			if refID < 0 {
				return nil, fmt.Errorf("unknown ref %q", rec.RName)
			}
		}
		mapped := !rec.IsUnmapped()
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

// writeTestBAM writes a coordinate-sorted BAM carrying records with rich
// variable-length fields (SEQ, QUAL, QNAME, aux) — precisely the data
// ReadDepthInto skips — so a byte-identical BAI proves the skip is correct.
func writeTestBAM(t *testing.T) ([]byte, *sam.Header) {
	t.Helper()
	h := &sam.Header{Refs: []sam.Reference{
		{Name: "chr1", Length: 100000},
		{Name: "chr2", Length: 50000},
	}}
	mkCigar := func(s string) sam.Cigar {
		c, err := sam.ParseCigar(s)
		if err != nil {
			t.Fatalf("ParseCigar(%q): %v", s, err)
		}
		return c
	}
	recs := []*sam.Record{
		{QName: "read1", Flag: 0, RName: "chr1", Pos: 100, MapQ: 60, Cigar: mkCigar("10M"),
			Seq: "ACGTACGTAC", Qual: []byte{30, 31, 32, 33, 34, 35, 36, 37, 38, 39},
			Aux: []sam.Aux{{Tag: "NM", Type: 'i', Value: int64(0)}, {Tag: "MD", Type: 'Z', Value: "10"}}},
		{QName: "read2", Flag: 0, RName: "chr1", Pos: 5000, MapQ: 40, Cigar: mkCigar("5M2I5M"),
			Seq: "ACGTAACGTA", Qual: []byte{20, 20, 20, 20, 20, 20, 20, 20, 20, 20},
			Aux: []sam.Aux{{Tag: "AS", Type: 'i', Value: int64(90)}}},
		{QName: "read3_longer_name", Flag: 16, RName: "chr1", Pos: 16384, MapQ: 55, Cigar: mkCigar("8M2D2M"),
			Seq: "TTTTGGGGCC", Qual: []byte{25, 25, 25, 25, 25, 25, 25, 25, 25, 25}},
		{QName: "read4", Flag: 0, RName: "chr2", Pos: 200, MapQ: 60, Cigar: mkCigar("10M"),
			Seq: "GGGGCCCCAA", Qual: []byte{40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
			Aux: []sam.Aux{{Tag: "RG", Type: 'Z', Value: "grp1"}}},
		// A placed-but-unmapped record (has ref+pos, FlagUnmapped set).
		{QName: "read5", Flag: sam.FlagUnmapped, RName: "chr2", Pos: 300, MapQ: 0, Cigar: mkCigar("*"),
			Seq: "NNNNNNNNNN", Qual: []byte{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
		// A truly unplaced record (no ref) — bumps n_no_coor.
		{QName: "read6", Flag: sam.FlagUnmapped, RName: "*", Pos: 0, MapQ: 0, Cigar: nil,
			Seq: "AAAACCCCGG", Qual: []byte{10, 10, 10, 10, 10, 10, 10, 10, 10, 10}},
	}

	var buf bytes.Buffer
	bw := sam.NewBAMWriter(&buf)
	if err := bw.WriteHeader(h); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, r := range recs {
		if err := bw.Write(r); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return buf.Bytes(), h
}

// TestBuildBAIByteIdenticalToFullDecode indexes the same BAM two ways — via the
// allocation-free ReadDepthInto path (BuildBAI) and via the full-record decode
// (buildBAIFull) — and asserts the serialised BAI bytes are identical. This is
// the unit-level guarantee behind the perf change: skipping SEQ/QUAL/QNAME/aux
// during indexing cannot alter the index.
func TestBuildBAIByteIdenticalToFullDecode(t *testing.T) {
	bamBytes, h := writeTestBAM(t)

	brFast, err := sam.NewBAMReader(bytes.NewReader(bamBytes))
	if err != nil {
		t.Fatalf("NewBAMReader (fast): %v", err)
	}
	fastIdx, err := BuildBAI(brFast, len(h.Refs))
	if err != nil {
		t.Fatalf("BuildBAI: %v", err)
	}

	brFull, err := sam.NewBAMReader(bytes.NewReader(bamBytes))
	if err != nil {
		t.Fatalf("NewBAMReader (full): %v", err)
	}
	fullIdx, err := buildBAIFull(brFull, len(h.Refs))
	if err != nil {
		t.Fatalf("buildBAIFull: %v", err)
	}

	var fastBuf, fullBuf bytes.Buffer
	if err := WriteBAI(&fastBuf, fastIdx); err != nil {
		t.Fatalf("write fast BAI: %v", err)
	}
	if err := WriteBAI(&fullBuf, fullIdx); err != nil {
		t.Fatalf("write full BAI: %v", err)
	}
	if !bytes.Equal(fastBuf.Bytes(), fullBuf.Bytes()) {
		t.Fatalf("BAI mismatch: ReadDepthInto path produced a different index than the full-decode path\n fast=%d bytes full=%d bytes",
			fastBuf.Len(), fullBuf.Len())
	}
}
