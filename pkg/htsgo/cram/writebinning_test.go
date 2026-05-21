package cram

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// binningTestRecords returns a small set of mapped and unmapped reads
// carrying varied raw quality strings, for the binning round-trip tests.
func binningTestRecords() []*sam.Record {
	mapped := &sam.Record{
		QName: "r1", Flag: 0, RName: "chr1", Pos: 100, MapQ: 30,
		RNext: "*", Seq: "ACGTAC", Qual: []byte{2, 8, 12, 24, 35, 41},
	}
	mapped.Cigar, _ = sam.ParseCigar("6M")
	rev := &sam.Record{
		QName: "r2", Flag: sam.FlagReverse, RName: "chr2", Pos: 200, MapQ: 20,
		RNext: "*", Seq: "TTGG", Qual: []byte{0, 19, 29, 40},
	}
	rev.Cigar, _ = sam.ParseCigar("4M")
	unmapped := &sam.Record{
		QName: "r3", Flag: sam.FlagUnmapped, RName: "*", Pos: 0,
		RNext: "*", Seq: "AAAA", Qual: []byte{5, 15, 25, 38},
	}
	return []*sam.Record{mapped, rev, unmapped}
}

// roundTripBinned writes records through a RecordWriter configured with
// the given binning scheme and reads them back, returning the decoded
// records.
func roundTripBinned(t *testing.T, h *sam.Header, records []*sam.Record, b QualityBinning) []*sam.Record {
	t.Helper()
	var buf bytes.Buffer
	rw, err := NewRecordWriterOpts(&buf, h, WriterOptions{Binning: b})
	if err != nil {
		t.Fatalf("NewRecordWriterOpts: %v", err)
	}
	for i, rec := range records {
		if err := rw.Write(rec); err != nil {
			t.Fatalf("Write record %d: %v", i, err)
		}
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	out, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	return out
}

// TestWriteCRAMBinningIllumina8 confirms that with BinningIllumina8 the
// decoded QUAL of every record equals the input quality mapped through
// the 8-level table — lossy, but deterministic.
func TestWriteCRAMBinningIllumina8(t *testing.T) {
	h := writerTestHeader()
	records := binningTestRecords()
	got := roundTripBinned(t, h, records, BinningIllumina8)
	if len(got) != len(records) {
		t.Fatalf("decoded %d records, want %d", len(got), len(records))
	}
	for i, rec := range records {
		want := BinningIllumina8.BinQuality(rec.Qual)
		if !bytes.Equal(got[i].Qual, want) {
			t.Errorf("record %d QUAL = %v, want binned %v (raw %v)",
				i, got[i].Qual, want, rec.Qual)
		}
		// The caller's record must not have been mutated.
		if bytes.Equal(rec.Qual, want) && !qualAllSame(rec.Qual) {
			t.Errorf("record %d caller QUAL was mutated to the binned values", i)
		}
	}
}

// TestWriteCRAMBinningNoneIsLossless confirms BinningNone stores quality
// verbatim — the decoded QUAL is exactly the input.
func TestWriteCRAMBinningNoneIsLossless(t *testing.T) {
	h := writerTestHeader()
	records := binningTestRecords()
	got := roundTripBinned(t, h, records, BinningNone)
	for i, rec := range records {
		if !bytes.Equal(got[i].Qual, rec.Qual) {
			t.Errorf("record %d QUAL = %v, want lossless %v", i, got[i].Qual, rec.Qual)
		}
	}
}

// TestWriteCRAMBinningDefaultBytesUnchanged confirms the default writer
// (WriterOptions zero value) produces byte-for-byte the same CRAM as
// WriteCRAM, so no existing caller is affected by the new option field.
func TestWriteCRAMBinningDefaultBytesUnchanged(t *testing.T) {
	h := writerTestHeader()
	records := binningTestRecords()

	var legacy bytes.Buffer
	if err := WriteCRAM(&legacy, h, records); err != nil {
		t.Fatalf("WriteCRAM: %v", err)
	}

	var opted bytes.Buffer
	rw, err := NewRecordWriterOpts(&opted, h, WriterOptions{})
	if err != nil {
		t.Fatalf("NewRecordWriterOpts: %v", err)
	}
	for _, rec := range records {
		if err := rw.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if !bytes.Equal(legacy.Bytes(), opted.Bytes()) {
		t.Errorf("zero-value WriterOptions changed the CRAM bytes (%d vs %d bytes)",
			legacy.Len(), opted.Len())
	}
}

// TestWriteCRAMBinningProvenance confirms a real binning scheme appends a
// @CO provenance line to the embedded SAM header, while BinningNone leaves
// the header untouched and the caller's header is never mutated.
func TestWriteCRAMBinningProvenance(t *testing.T) {
	h := writerTestHeader()
	origText := h.Text()

	var buf bytes.Buffer
	if _, err := func() (*RecordWriter, error) {
		return NewRecordWriterOpts(&buf, h, WriterOptions{Binning: BinningIllumina8})
	}(); err != nil {
		t.Fatalf("NewRecordWriterOpts: %v", err)
	}
	// The caller's header must not have been modified in place.
	if h.Text() != origText {
		t.Error("the caller's SAM header was mutated by the binning provenance line")
	}

	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	if !strings.Contains(rr.Header().Text(), "lossy quality-score binning") {
		t.Errorf("embedded header missing binning provenance:\n%s", rr.Header().Text())
	}

	// BinningNone must not add any provenance line.
	var plain bytes.Buffer
	if _, err := NewRecordWriterOpts(&plain, h, WriterOptions{Binning: BinningNone}); err != nil {
		t.Fatalf("NewRecordWriterOpts(none): %v", err)
	}
	pr, err := NewRecordReader(bytes.NewReader(plain.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader(none): %v", err)
	}
	if strings.Contains(pr.Header().Text(), "lossy quality-score binning") {
		t.Error("BinningNone added an unexpected provenance line")
	}
}

// TestNewRecordWriterOptsRejectsUnknownScheme confirms an out-of-range
// binning scheme is a clear construction error.
func TestNewRecordWriterOptsRejectsUnknownScheme(t *testing.T) {
	h := writerTestHeader()
	var buf bytes.Buffer
	if _, err := NewRecordWriterOpts(&buf, h, WriterOptions{Binning: QualityBinning(99)}); err == nil {
		t.Fatal("NewRecordWriterOpts accepted an unknown binning scheme")
	}
}

// qualAllSame reports whether every byte of q is identical (so binning it
// is a no-op and cannot evidence mutation).
func qualAllSame(q []byte) bool {
	for i := 1; i < len(q); i++ {
		if q[i] != q[0] {
			return false
		}
	}
	return true
}
