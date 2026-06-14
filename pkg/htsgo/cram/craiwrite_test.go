package cram

import (
	"bytes"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestWriteCRAIRoundTrip confirms WriteCRAI and ReadCRAI are inverses:
// a hand-built entry list survives a serialise/parse cycle unchanged.
func TestWriteCRAIRoundTrip(t *testing.T) {
	want := []CRAIEntry{
		{RefID: 0, AlignmentStart: 100, AlignmentSpan: 50, ContainerOffset: 26, SliceOffset: 12, SliceSize: 300},
		{RefID: 1, AlignmentStart: 1, AlignmentSpan: 10, ContainerOffset: 400, SliceOffset: 14, SliceSize: 120},
		{RefID: -1, AlignmentStart: 0, AlignmentSpan: 0, ContainerOffset: 600, SliceOffset: 14, SliceSize: 80},
	}
	var buf bytes.Buffer
	if err := WriteCRAI(&buf, want); err != nil {
		t.Fatalf("WriteCRAI: %v", err)
	}
	idx, err := ReadCRAI(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("ReadCRAI: %v", err)
	}
	if len(idx.Entries) != len(want) {
		t.Fatalf("got %d entries, want %d", len(idx.Entries), len(want))
	}
	for i, e := range idx.Entries {
		if e != want[i] {
			t.Errorf("entry %d = %+v, want %+v", i, e, want[i])
		}
	}
}

// TestBuildCRAIFromWrittenCRAM builds a .crai from a CRAM produced by the
// in-tree writer and asserts the index is sane: one entry per slice,
// offsets within the file, and an overlap query returns the right slice.
func TestBuildCRAIFromWrittenCRAM(t *testing.T) {
	h := writerTestHeader()
	var records []*sam.Record
	for i := 0; i < 25; i++ {
		records = append(records, mkRec("read", "chr1", int32(100+i*10), "6M", "ACGTAC"))
	}
	var cram bytes.Buffer
	rw, err := NewRecordWriter(&cram, h)
	if err != nil {
		t.Fatalf("NewRecordWriter: %v", err)
	}
	rw.recordsPerSlice = 10 // three containers => three slices.
	for i, rec := range records {
		if err := rw.Write(rec); err != nil {
			t.Fatalf("Write %d: %v", i, err)
		}
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	entries, err := BuildCRAI(bytes.NewReader(cram.Bytes()))
	if err != nil {
		t.Fatalf("BuildCRAI: %v", err)
	}
	// Writer emits one slice per container; 25 records at 10/slice => 3.
	if len(entries) != 3 {
		t.Fatalf("got %d index entries, want 3", len(entries))
	}
	fileLen := int64(cram.Len())
	for i, e := range entries {
		if e.RefID != 0 {
			t.Errorf("entry %d RefID = %d, want 0 (chr1)", i, e.RefID)
		}
		if e.ContainerOffset <= 0 || e.ContainerOffset >= fileLen {
			t.Errorf("entry %d ContainerOffset = %d, out of file bounds [0,%d)", i, e.ContainerOffset, fileLen)
		}
		if e.SliceOffset <= 0 {
			t.Errorf("entry %d SliceOffset = %d, want > 0", i, e.SliceOffset)
		}
		if e.SliceSize <= 0 || e.ContainerOffset+e.SliceOffset+e.SliceSize > fileLen {
			t.Errorf("entry %d slice [%d+%d+%d] runs past file end %d",
				i, e.ContainerOffset, e.SliceOffset, e.SliceSize, fileLen)
		}
		if e.AlignmentSpan <= 0 {
			t.Errorf("entry %d AlignmentSpan = %d, want > 0", i, e.AlignmentSpan)
		}
	}
	// Containers must be emitted in increasing file order.
	for i := 1; i < len(entries); i++ {
		if entries[i].ContainerOffset <= entries[i-1].ContainerOffset {
			t.Errorf("entry %d container offset %d not after entry %d offset %d",
				i, entries[i].ContainerOffset, i-1, entries[i-1].ContainerOffset)
		}
	}

	// Round-trip the entries through WriteCRAI/ReadCRAI and run a query.
	var craiBuf bytes.Buffer
	if err := WriteCRAI(&craiBuf, entries); err != nil {
		t.Fatalf("WriteCRAI: %v", err)
	}
	idx, err := ReadCRAI(bytes.NewReader(craiBuf.Bytes()))
	if err != nil {
		t.Fatalf("ReadCRAI: %v", err)
	}
	// The first slice covers reads at 100..150; a query on [100,106)
	// (0-based) overlaps only the first slice's span.
	hits := idx.Query(0, 100, 106)
	if len(hits) != 1 {
		t.Fatalf("Query(0,100,106) returned %d hits, want 1", len(hits))
	}
	if hits[0] != entries[0] {
		t.Errorf("Query hit = %+v, want first entry %+v", hits[0], entries[0])
	}
	// A whole-reference query overlaps every slice.
	all := idx.Query(0, 0, 0)
	if len(all) != 3 {
		t.Errorf("whole-reference query returned %d hits, want 3", len(all))
	}
}

// TestBuildCRAIReferenceCRAM builds a .crai from the upstream samtools
// fixture and asserts the entries are sane. The samtools fixture is a hard
// requirement for the CRAM index parity rig, so a missing submodule is a
// fatal error (with an init hint), not a skip.
func TestBuildCRAIReferenceCRAM(t *testing.T) {
	const path = "../../../reference_code/samtools/test/dat/test_input_1_a.cram"
	entries, err := BuildCRAIFile(path)
	if err != nil {
		t.Fatalf("fixture %s unavailable: %v; run `git submodule update --init reference_code/samtools`", path, err)
	}
	if len(entries) == 0 {
		t.Fatal("BuildCRAI returned no entries for a non-empty CRAM")
	}
	for i, e := range entries {
		if e.ContainerOffset < 0 || e.SliceOffset < 0 || e.SliceSize <= 0 {
			t.Errorf("entry %d has a bad offset/size: %+v", i, e)
		}
		if e.AlignmentSpan < 0 {
			t.Errorf("entry %d has a negative span: %+v", i, e)
		}
	}
	// The serialised index must parse straight back.
	var buf bytes.Buffer
	if err := WriteCRAI(&buf, entries); err != nil {
		t.Fatalf("WriteCRAI: %v", err)
	}
	idx, err := ReadCRAI(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("ReadCRAI: %v", err)
	}
	if len(idx.Entries) != len(entries) {
		t.Errorf("re-parsed %d entries, built %d", len(idx.Entries), len(entries))
	}
}

// FuzzCRAIWriteRoundTrip checks that any entry list WriteCRAI accepts is
// recovered byte-identically by ReadCRAI. The .crai fields are bounded to
// the ranges parseCRAILine accepts so the fuzzer exercises the
// serialise/parse pair rather than the parser's rejection paths.
func FuzzCRAIWriteRoundTrip(f *testing.F) {
	f.Add(int32(0), int32(100), int32(50), int64(26), int64(12), int64(300))
	f.Add(int32(-1), int32(0), int32(0), int64(0), int64(0), int64(1))
	f.Fuzz(func(t *testing.T, ref, start, span int32, cOff, sOff, sSize int64) {
		// Keep the inputs within parseCRAILine's accepted ranges.
		if span < 0 {
			span = -span
		}
		if cOff < 0 {
			cOff = -cOff
		}
		if sOff < 0 {
			sOff = -sOff
		}
		if sSize < 0 {
			sSize = -sSize
		}
		want := CRAIEntry{
			RefID:           ref,
			AlignmentStart:  start,
			AlignmentSpan:   span,
			ContainerOffset: cOff,
			SliceOffset:     sOff,
			SliceSize:       sSize,
		}
		var buf bytes.Buffer
		if err := WriteCRAI(&buf, []CRAIEntry{want}); err != nil {
			t.Fatalf("WriteCRAI: %v", err)
		}
		idx, err := ReadCRAI(bytes.NewReader(buf.Bytes()))
		if err != nil {
			t.Fatalf("ReadCRAI: %v", err)
		}
		if len(idx.Entries) != 1 || idx.Entries[0] != want {
			t.Fatalf("round-trip mismatch: got %+v, want %+v", idx.Entries, want)
		}
	})
}
