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
	got, err := parseSliceHeader(wire, 3)
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
	got, err := parseSliceHeader(buildSliceHeader(s), 3)
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
			if _, err := parseSliceHeader(full[:n], 3); err == nil {
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
	if _, err := parseSliceHeader(b.Bytes(), 3); err == nil {
		t.Errorf("negative block count should be rejected")
	}
}

// TestParseSliceHeaderRecordCounterWidth proves the record-counter field
// is read as ITF-8 (32-bit) for CRAM v2 and LTF-8 (64-bit) for v3+,
// matching htslib's cram_decode_slice_header. The two encodings diverge
// for a counter >= 2^28: the same wire bytes parse to a different value
// and, crucially, leave the cursor at a different offset so every later
// field (block count, content ids, MD5) stays aligned only when the
// width matches the version.
func TestParseSliceHeaderRecordCounterWidth(t *testing.T) {
	// A counter at the 2^28 boundary, where ITF-8 and LTF-8 first differ:
	// ITF-8 encodes 32-bit values in up to 5 bytes, LTF-8 up to 9, and a
	// value needing the high ITF-8 bits is laid out differently from the
	// LTF-8 form of the same number.
	const counter = int64(1) << 28

	build := func(enc []byte) []byte {
		var b bytes.Buffer
		b.Write(encITF8(2))  // ref id
		b.Write(encITF8(50)) // start
		b.Write(encITF8(10)) // span
		b.Write(encITF8(7))  // records
		b.Write(enc)         // record counter (width varies)
		b.Write(encITF8(1))  // block count
		b.Write(encITF8(1))  // content id count
		b.Write(encITF8(99)) // content id
		b.Write(encITF8(-1)) // embedded ref id
		b.Write(make([]byte, 16))
		return b.Bytes()
	}

	// v2: counter is ITF-8. Parsed with major=2, the value round-trips and
	// the trailing fields line up.
	v2wire := build(encITF8(int32(counter)))
	gotV2, err := parseSliceHeader(v2wire, 2)
	if err != nil {
		t.Fatalf("v2 parseSliceHeader: %v", err)
	}
	if gotV2.RecordCounter != counter {
		t.Errorf("v2 RecordCounter = %d, want %d", gotV2.RecordCounter, counter)
	}
	if gotV2.NumBlocks != 1 || len(gotV2.BlockContentIDs) != 1 || gotV2.BlockContentIDs[0] != 99 {
		t.Errorf("v2 trailing fields misaligned: blocks=%d ids=%v", gotV2.NumBlocks, gotV2.BlockContentIDs)
	}

	// v3: counter is LTF-8. Parsed with major=3 it round-trips identically.
	v3wire := build(encLTF8(counter))
	gotV3, err := parseSliceHeader(v3wire, 3)
	if err != nil {
		t.Fatalf("v3 parseSliceHeader: %v", err)
	}
	if gotV3.RecordCounter != counter {
		t.Errorf("v3 RecordCounter = %d, want %d", gotV3.RecordCounter, counter)
	}
	if gotV3.NumBlocks != 1 || len(gotV3.BlockContentIDs) != 1 || gotV3.BlockContentIDs[0] != 99 {
		t.Errorf("v3 trailing fields misaligned: blocks=%d ids=%v", gotV3.NumBlocks, gotV3.BlockContentIDs)
	}

	// Reading a v2 (ITF-8) counter with the LTF-8 path — the old
	// always-LTF-8 behaviour — misparses: at >= 2^28 the ITF-8 and LTF-8
	// byte layouts differ, so the counter value and/or the subsequent
	// field alignment is wrong. This is exactly the bug the version
	// threading fixes.
	if wrong, werr := parseSliceHeader(v2wire, 3); werr == nil {
		if wrong.RecordCounter == counter &&
			wrong.NumBlocks == 1 && len(wrong.BlockContentIDs) == 1 &&
			wrong.BlockContentIDs[0] == 99 {
			t.Errorf("expected v2 ITF-8 counter to misparse when read as LTF-8 at value %d, but it matched", counter)
		}
	}

	// Below 2^28 the two encodings coincide, so version makes no
	// difference — the realistic-file invariant the roadmap relied on.
	const small = int64(1234567)
	smallWire := build(encITF8(int32(small)))
	for _, major := range []uint8{2, 3} {
		got, err := parseSliceHeader(smallWire, major)
		if err != nil {
			t.Fatalf("small counter major=%d: %v", major, err)
		}
		if got.RecordCounter != small {
			t.Errorf("small counter major=%d: got %d want %d", major, got.RecordCounter, small)
		}
	}
}
