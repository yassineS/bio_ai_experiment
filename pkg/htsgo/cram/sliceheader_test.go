package cram

import (
	"bytes"
	"testing"
)

// buildSliceHeader assembles a CRAM slice-header block payload from the
// given field values.
func buildSliceHeader(s SliceHeader) []byte {
	var b bytes.Buffer
	b.Write(encITF8(s.RefSeqID))
	b.Write(encITF8(s.AlignmentStart))
	b.Write(encITF8(s.AlignmentSpan))
	b.Write(encITF8(s.NumRecords))
	b.Write(encLTF8(s.RecordCounter))
	b.Write(encITF8(s.NumBlocks))
	b.Write(encITF8(int32(len(s.BlockContentIDs))))
	for _, id := range s.BlockContentIDs {
		b.Write(encITF8(id))
	}
	b.Write(encITF8(s.EmbeddedRefID))
	b.Write(s.ReferenceMD5[:])
	b.Write(s.Tags)
	return b.Bytes()
}

// TestParseSliceHeaderRoundTrip builds a slice header, parses it back
// and checks every field.
func TestParseSliceHeaderRoundTrip(t *testing.T) {
	want := SliceHeader{
		RefSeqID:        3,
		AlignmentStart:  100,
		AlignmentSpan:   250,
		NumRecords:      42,
		RecordCounter:   1000,
		NumBlocks:       5,
		BlockContentIDs: []int32{1, 2, 11, 37, 42},
		EmbeddedRefID:   -1,
	}
	for i := range want.ReferenceMD5 {
		want.ReferenceMD5[i] = byte(i + 1)
	}
	wire := buildSliceHeader(want)
	got, err := parseSliceHeader(wire)
	if err != nil {
		t.Fatalf("parseSliceHeader: %v", err)
	}
	if got.RefSeqID != want.RefSeqID || got.AlignmentStart != want.AlignmentStart ||
		got.AlignmentSpan != want.AlignmentSpan || got.NumRecords != want.NumRecords ||
		got.RecordCounter != want.RecordCounter || got.NumBlocks != want.NumBlocks ||
		got.EmbeddedRefID != want.EmbeddedRefID {
		t.Errorf("scalar mismatch: got %+v want %+v", got, want)
	}
	if len(got.BlockContentIDs) != len(want.BlockContentIDs) {
		t.Fatalf("block id count: got %d want %d", len(got.BlockContentIDs), len(want.BlockContentIDs))
	}
	for i := range want.BlockContentIDs {
		if got.BlockContentIDs[i] != want.BlockContentIDs[i] {
			t.Errorf("block id %d: got %d want %d", i, got.BlockContentIDs[i], want.BlockContentIDs[i])
		}
	}
	if got.ReferenceMD5 != want.ReferenceMD5 {
		t.Errorf("MD5 mismatch: got %x want %x", got.ReferenceMD5, want.ReferenceMD5)
	}
}

// TestParseSliceHeaderWithTags checks trailing optional tag bytes are
// retained verbatim.
func TestParseSliceHeaderWithTags(t *testing.T) {
	s := SliceHeader{NumBlocks: 1, BlockContentIDs: []int32{5}, EmbeddedRefID: -1,
		Tags: []byte{0xde, 0xad, 0xbe, 0xef}}
	got, err := parseSliceHeader(buildSliceHeader(s))
	if err != nil {
		t.Fatalf("parseSliceHeader: %v", err)
	}
	if !bytes.Equal(got.Tags, []byte{0xde, 0xad, 0xbe, 0xef}) {
		t.Errorf("Tags = %v, want [de ad be ef]", got.Tags)
	}
}

// TestParseSliceHeaderTruncated checks the parser errors — never panics —
// on input truncated at every byte offset.
func TestParseSliceHeaderTruncated(t *testing.T) {
	full := buildSliceHeader(SliceHeader{
		RefSeqID: 1, NumBlocks: 2, BlockContentIDs: []int32{1, 2}, EmbeddedRefID: -1,
	})
	for n := 0; n < len(full); n++ {
		func() {
			defer func() {
				if r := recover(); r != nil {
					t.Errorf("parseSliceHeader panicked at truncation %d: %v", n, r)
				}
			}()
			if _, err := parseSliceHeader(full[:n]); err == nil {
				t.Errorf("expected error for slice header truncated to %d bytes", n)
			}
		}()
	}
}

// TestParseSliceHeaderNegativeCounts rejects negative block and content
// id counts.
func TestParseSliceHeaderNegativeCounts(t *testing.T) {
	var b bytes.Buffer
	b.Write(encITF8(0)) // ref id
	b.Write(encITF8(0)) // start
	b.Write(encITF8(0)) // span
	b.Write(encITF8(0)) // records
	b.Write(encLTF8(0)) // counter
	b.Write(encITF8(-1))
	if _, err := parseSliceHeader(b.Bytes()); err == nil {
		t.Errorf("negative block count should be rejected")
	}
}
