package sam

import (
	"bufio"
	"bytes"
	"io"
	"testing"
)

// fastPathRecords exercises every SAM field and every BAM aux type the
// direct BAM->SAM serialiser must reproduce: all integer widths (c/C/s/S/i/I),
// float (f), single char (A), string (Z), hex (H), and B arrays of each
// subtype, plus the SEQ/QUAL "*" conventions and the '=' RNEXT shorthand.
func fastPathRecords() []*Record {
	return []*Record{
		// Mapped read with a representative mix of aux tags.
		{
			QName: "read.one", Flag: 99, RName: "chr1", Pos: 100, MapQ: 60,
			Cigar: mustCigar("5M2I3M"), RNext: "=", PNext: 250, TLen: 200,
			Seq: "ACGTAGCTA", Qual: []byte{30, 31, 32, 33, 34, 35, 36, 37, 38},
			Aux: []Aux{
				{Tag: "NM", Type: 'c', Value: int64(2)},
				{Tag: "MC", Type: 'C', Value: int64(200)},
				{Tag: "ms", Type: 's', Value: int64(-1234)},
				{Tag: "MS", Type: 'S', Value: int64(40000)},
				{Tag: "AS", Type: 'i', Value: int64(-2000000000)},
				{Tag: "XI", Type: 'I', Value: int64(3000000000)},
				{Tag: "XF", Type: 'f', Value: float64(float32(1.5))},
				{Tag: "XA", Type: 'A', Value: "U"},
				{Tag: "RG", Type: 'Z', Value: "group-1"},
				{Tag: "HX", Type: 'H', Value: "DEADBEEF"},
				{Tag: "Bc", Type: 'B', ArrayType: 'c', ArrayValues: []interface{}{int64(-1), int64(2), int64(-3)}},
				{Tag: "BC", Type: 'B', ArrayType: 'C', ArrayValues: []interface{}{int64(1), int64(255)}},
				{Tag: "Bs", Type: 'B', ArrayType: 's', ArrayValues: []interface{}{int64(-100), int64(100)}},
				{Tag: "BS", Type: 'B', ArrayType: 'S', ArrayValues: []interface{}{int64(0), int64(65535)}},
				{Tag: "Bi", Type: 'B', ArrayType: 'i', ArrayValues: []interface{}{int64(-7), int64(7)}},
				{Tag: "BI", Type: 'B', ArrayType: 'I', ArrayValues: []interface{}{int64(1), int64(4000000000)}},
				{Tag: "Bf", Type: 'B', ArrayType: 'f', ArrayValues: []interface{}{float64(float32(0.25)), float64(float32(-3.5))}},
			},
		},
		// Reverse-strand mate-reference differs from this ref (RNEXT not '=').
		{
			QName: "read.two", Flag: 147, RName: "chr2", Pos: 555, MapQ: 0,
			Cigar: mustCigar("10M"), RNext: "chr1", PNext: 99, TLen: -200,
			Seq: "TTTTTAAAAA", Qual: []byte{2, 2, 2, 2, 2, 40, 40, 40, 40, 40},
		},
		// Unmapped read: SEQ present, QUAL all-0xFF -> "*", no CIGAR.
		{
			QName: "read.unmapped", Flag: 77, RName: "", Pos: 0, MapQ: 0,
			Seq: "ACGT", Qual: []byte{0xff, 0xff, 0xff, 0xff},
		},
		// SEQ "*" / QUAL "*" (both absent).
		{
			QName: "read.noseq", Flag: 141, RName: "", Pos: 0, MapQ: 0,
		},
	}
}

// fastPathTestHeader names the two references the test records use.
func fastPathTestHeader() *Header {
	h := &Header{Refs: []Reference{{Name: "chr1", Length: 1000}, {Name: "chr2", Length: 2000}}}
	h.Lines = []HeaderLine{
		{Tag: "HD", Fields: []HeaderField{{Tag: "VN", Value: "1.6"}}},
		{Tag: "SQ", Fields: []HeaderField{{Tag: "SN", Value: "chr1"}, {Tag: "LN", Value: "1000"}}},
		{Tag: "SQ", Fields: []HeaderField{{Tag: "SN", Value: "chr2"}, {Tag: "LN", Value: "2000"}}},
	}
	return h
}

// TestWriteSAMBodyMatchesSAMWriter is the core correctness gate for the
// BAM->SAM fast path: a record serialised straight from its raw BAM bytes via
// ReadSAMInto + WriteSAMBody must be byte-identical to decoding it into a full
// Record and formatting it through SAMWriter.Write. It covers every aux type
// and the SEQ/QUAL "*" conventions, since the fast path bypasses the
// intermediate Record entirely.
func TestWriteSAMBodyMatchesSAMWriter(t *testing.T) {
	hdr := fastPathTestHeader()
	recs := fastPathRecords()

	// Encode the records to an in-memory BAM stream.
	var bamBuf bytes.Buffer
	bw := NewBAMWriter(&bamBuf)
	if err := bw.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, r := range recs {
		if err := bw.Write(r); err != nil {
			t.Fatalf("BAM Write(%s): %v", r.QName, err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("BAM Close: %v", err)
	}

	// Reference path: decode each record and format via SAMWriter.
	wantR, err := NewBAMReader(bytes.NewReader(bamBuf.Bytes()))
	if err != nil {
		t.Fatalf("NewBAMReader (ref): %v", err)
	}
	var wantBuf bytes.Buffer
	sw := NewSAMWriter(&wantBuf)
	if err := sw.WriteHeader(nil); err != nil {
		t.Fatalf("SAMWriter.WriteHeader: %v", err)
	}
	for {
		rec, rerr := wantR.Read()
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			t.Fatalf("ref Read: %v", rerr)
		}
		if werr := sw.Write(rec); werr != nil {
			t.Fatalf("SAMWriter.Write: %v", werr)
		}
	}
	if err := sw.Close(); err != nil {
		t.Fatalf("SAMWriter.Close: %v", err)
	}

	// Fast path: serialise each record's raw bytes via WriteSAMBody.
	fastR, err := NewBAMReader(bytes.NewReader(bamBuf.Bytes()))
	if err != nil {
		t.Fatalf("NewBAMReader (fast): %v", err)
	}
	var gotBuf bytes.Buffer
	fbw := bufio.NewWriter(&gotBuf)
	var ff FastFields
	for {
		ferr := fastR.ReadSAMInto(&ff)
		if ferr == io.EOF {
			break
		}
		if ferr != nil {
			t.Fatalf("ReadSAMInto: %v", ferr)
		}
		if werr := fastR.WriteSAMBody(fbw, &ff); werr != nil {
			t.Fatalf("WriteSAMBody: %v", werr)
		}
	}
	if err := fbw.Flush(); err != nil {
		t.Fatalf("flush: %v", err)
	}

	if !bytes.Equal(wantBuf.Bytes(), gotBuf.Bytes()) {
		t.Fatalf("fast path output differs from SAMWriter:\nwant:\n%s\ngot:\n%s",
			wantBuf.String(), gotBuf.String())
	}
}

// TestReadSAMIntoFields checks the cheap fixed-prefix decode populates the
// fields the view filter pipeline relies on (flag, MAPQ, RName, Pos, RefSpan).
func TestReadSAMIntoFields(t *testing.T) {
	hdr := fastPathTestHeader()
	recs := fastPathRecords()
	var bamBuf bytes.Buffer
	bw := NewBAMWriter(&bamBuf)
	if err := bw.WriteHeader(hdr); err != nil {
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

	r, err := NewBAMReader(bytes.NewReader(bamBuf.Bytes()))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	var ff FastFields
	// First record: chr1, pos 100, flag 99, mapq 60, CIGAR 5M2I3M -> ref span 8.
	if err := r.ReadSAMInto(&ff); err != nil {
		t.Fatalf("ReadSAMInto: %v", err)
	}
	if ff.Flag != 99 || ff.MapQ != 60 || ff.RName != "chr1" || ff.Pos != 100 || ff.RefSpan != 8 {
		t.Fatalf("rec1 fields: flag=%d mapq=%d rname=%q pos=%d span=%d",
			ff.Flag, ff.MapQ, ff.RName, ff.Pos, ff.RefSpan)
	}
	// Second record: chr2, pos 555, CIGAR 10M -> span 10.
	if err := r.ReadSAMInto(&ff); err != nil {
		t.Fatalf("ReadSAMInto: %v", err)
	}
	if ff.RName != "chr2" || ff.Pos != 555 || ff.RefSpan != 10 {
		t.Fatalf("rec2 fields: rname=%q pos=%d span=%d", ff.RName, ff.Pos, ff.RefSpan)
	}
	// Third record: unmapped (no RName) -> RName "", Pos 0, no CIGAR -> span 0.
	if err := r.ReadSAMInto(&ff); err != nil {
		t.Fatalf("ReadSAMInto: %v", err)
	}
	if ff.RName != "" || ff.Pos != 0 || ff.RefSpan != 0 {
		t.Fatalf("rec3 fields: rname=%q pos=%d span=%d", ff.RName, ff.Pos, ff.RefSpan)
	}
}
