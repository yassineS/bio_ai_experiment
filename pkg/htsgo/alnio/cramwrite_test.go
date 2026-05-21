package alnio

import (
	"bytes"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// cramWriterTestHeader builds a small SAM header for the CRAM writer
// adapter tests.
func cramWriterTestHeader(t *testing.T) *sam.Header {
	t.Helper()
	h, err := sam.ParseHeaderText("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n")
	if err != nil {
		t.Fatalf("ParseHeaderText: %v", err)
	}
	return h
}

// TestCRAMWriterRoundTrip drives the alnio CRAM writer adapter through
// the sam.Writer interface and confirms the output is a CRAM stream the
// CRAM reader decodes back to the same records.
func TestCRAMWriterRoundTrip(t *testing.T) {
	h := cramWriterTestHeader(t)
	rec := &sam.Record{
		QName: "r1", Flag: 0, RName: "chr1", Pos: 100, MapQ: 30,
		RNext: "*", Seq: "ACGTAC",
	}
	rec.Cigar, _ = sam.ParseCigar("6M")

	var buf bytes.Buffer
	w := NewCRAMWriter(&buf)
	if err := w.WriteHeader(h); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := w.Write(rec); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	rr, err := cram.NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(recs) != 1 {
		t.Fatalf("decoded %d records, want 1", len(recs))
	}
	if recs[0].QName != "r1" || recs[0].Pos != 100 || recs[0].RName != "chr1" {
		t.Errorf("decoded record = %+v, want r1/chr1/100", recs[0])
	}
}

// TestCRAMWriterNilHeader confirms the adapter rejects a nil header,
// since CRAM cannot be written without one.
func TestCRAMWriterNilHeader(t *testing.T) {
	w := NewCRAMWriter(&bytes.Buffer{})
	if err := w.WriteHeader(nil); err == nil {
		t.Fatal("WriteHeader(nil) succeeded, want an error")
	}
}

// TestCRAMWriterWriteBeforeHeader confirms a Write before WriteHeader is
// an error rather than a panic.
func TestCRAMWriterWriteBeforeHeader(t *testing.T) {
	w := NewCRAMWriter(&bytes.Buffer{})
	if err := w.Write(&sam.Record{QName: "x"}); err == nil {
		t.Fatal("Write before WriteHeader succeeded, want an error")
	}
}
