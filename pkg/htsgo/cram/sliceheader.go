package cram

import (
	"fmt"
)

// SliceHeader is the fully-parsed data of a CRAM slice-header block (a
// block of content type MAPPED_SLICE / type 2). It locates a slice
// within its container's reference span and lists the data blocks the
// slice owns.
type SliceHeader struct {
	// RefSeqID is the reference sequence the slice's records align to:
	// >= 0 indexes the SAM @SQ lines, -1 means unmapped reads and -2
	// means the slice spans multiple references.
	RefSeqID int32
	// AlignmentStart is the 1-based reference coordinate of the slice's
	// first record. CRAM v4 stores it as a 64-bit varint; it is decoded
	// 64-bit-clean and narrowed to this int32 field (within int32 range a
	// fixture is unaffected). The 64-bit decode means a v4 slice whose
	// coordinate exceeds 2^31 is read without misaligning the fields that
	// follow it; only this header field narrows.
	AlignmentStart int32
	// AlignmentSpan is the number of reference bases the slice covers. As
	// with AlignmentStart, CRAM v4 stores it as a 64-bit varint, narrowed
	// to int32 here.
	AlignmentSpan int32
	// NumRecords is the count of alignment records in the slice.
	NumRecords int32
	// RecordCounter is the running record total across all slices that
	// precede this one in the file.
	RecordCounter int64
	// NumBlocks is the count of data blocks the slice owns (the CORE
	// block plus the external blocks), not counting the slice-header
	// block itself.
	NumBlocks int32
	// BlockContentIDs lists the content identifiers of the slice's data
	// blocks, in file order.
	BlockContentIDs []int32
	// EmbeddedRefID is the content id of an embedded reference block, or
	// -1 when the slice carries no embedded reference.
	EmbeddedRefID int32
	// ReferenceMD5 is the 16-byte MD5 digest of the slice's reference
	// region, used to verify the decoder is using the right reference.
	ReferenceMD5 [16]byte
	// Tags holds the optional trailing tag bytes some writers append
	// after the MD5; it is nil when none are present.
	Tags []byte
}

// parseSliceHeader parses the data of a CRAM slice-header block. The
// payload is a sequence of variable-length integer fields: ref id,
// alignment start, alignment span, record count, record counter, block
// count, an array of block content ids, the embedded-reference block id,
// the 16-byte reference MD5, and any optional trailing tag bytes.
//
// The major argument is the container's CRAM major version, which selects
// the integer encoding throughout (cram/cram_decode.c,
// cram_decode_slice_header):
//
//   - CRAM v2/v3 use ITF-8 (32-bit) and LTF-8 (64-bit). The record
//     counter is ITF-8 for v2 and LTF-8 for v3+; the two coincide below
//     2^28, so the distinction only matters for a v2 slice whose
//     preceding records number 2^28 or more.
//   - CRAM v4 uses uint7 varints: ref_seq_id is a signed varint, the
//     alignment start and span are 64-bit varints (decoded 64-bit-clean
//     and narrowed to the int32 header fields), the record counter is a
//     64-bit varint, and the remaining counts and ids are unsigned
//     varints.
func parseSliceHeader(p []byte, major uint8) (*SliceHeader, error) {
	r := newIntReader(major)
	s := &SliceHeader{}
	off := 0
	var err error
	// ref_seq_id: signed. The v3 ITF-8 reader already yields the correct
	// signed value (it sign-extends); v4 uses the zig-zag signed varint.
	if s.RefSeqID, off, err = sliceSignedInt(r, p, off, "slice ref seq id"); err != nil {
		return nil, err
	}
	// Alignment start/span: 64-bit for v4, 32-bit for v2/v3. Decode wide
	// and narrow to the int32 header fields.
	var start64, span64 int64
	if start64, off, err = sliceWide(r, p, off, major, "slice alignment start"); err != nil {
		return nil, err
	}
	if span64, off, err = sliceWide(r, p, off, major, "slice alignment span"); err != nil {
		return nil, err
	}
	s.AlignmentStart = int32(start64)
	s.AlignmentSpan = int32(span64)
	if s.NumRecords, off, err = sliceUnsignedInt(r, p, off, "slice record count"); err != nil {
		return nil, err
	}
	// CRAM v2 stores the record counter as a 32-bit varint; v3+ (including
	// v4) use a 64-bit varint. Reading the wrong width would misalign
	// every field that follows once the counter reaches 2^28.
	if major <= 2 {
		var rc int32
		if rc, off, err = sliceUnsignedInt(r, p, off, "slice record counter"); err != nil {
			return nil, err
		}
		s.RecordCounter = int64(rc)
	} else if s.RecordCounter, off, err = sliceUnsigned64(r, p, off, "slice record counter"); err != nil {
		return nil, err
	}
	if s.NumBlocks, off, err = sliceUnsignedInt(r, p, off, "slice block count"); err != nil {
		return nil, err
	}
	if s.NumBlocks < 0 {
		return nil, fmt.Errorf("cram: slice header declares negative block count %d", s.NumBlocks)
	}
	var nIDs int32
	if nIDs, off, err = sliceUnsignedInt(r, p, off, "slice block content id count"); err != nil {
		return nil, err
	}
	if nIDs < 0 {
		return nil, fmt.Errorf("cram: slice header declares negative block content id count %d", nIDs)
	}
	// Each content id is at least one byte; a count larger than the
	// remaining payload is corrupt and bounds the slice growth.
	if int(nIDs) > len(p)-off {
		return nil, fmt.Errorf("cram: slice header declares %d block content ids but only %d bytes remain", nIDs, len(p)-off)
	}
	for i := int32(0); i < nIDs; i++ {
		var id int32
		if id, off, err = sliceUnsignedInt(r, p, off, "slice block content id"); err != nil {
			return nil, err
		}
		s.BlockContentIDs = append(s.BlockContentIDs, id)
	}
	if s.EmbeddedRefID, off, err = sliceUnsignedInt(r, p, off, "slice embedded reference id"); err != nil {
		return nil, err
	}
	if off+16 > len(p) {
		return nil, fmt.Errorf("cram: slice header truncated reading the 16-byte reference MD5")
	}
	copy(s.ReferenceMD5[:], p[off:off+16])
	off += 16
	// Any bytes past the MD5 are optional writer-appended tags; retain
	// them verbatim for C4b without interpreting them.
	if off < len(p) {
		s.Tags = append([]byte(nil), p[off:]...)
	}
	return s, nil
}

// sliceUnsignedInt reads a version-aware unsigned 32-bit field at off
// within p, advancing off and wrapping any error with what for context.
func sliceUnsignedInt(r intReader, p []byte, off int, what string) (int32, int, error) {
	v, n, err := r.u32(p, off)
	if err != nil {
		return 0, off, fmt.Errorf("cram: %s: %w", what, err)
	}
	return v, off + n, nil
}

// sliceSignedInt reads a version-aware signed 32-bit field at off within
// p, advancing off and wrapping any error with what for context.
func sliceSignedInt(r intReader, p []byte, off int, what string) (int32, int, error) {
	v, n, err := r.s32(p, off)
	if err != nil {
		return 0, off, fmt.Errorf("cram: %s: %w", what, err)
	}
	return v, off + n, nil
}

// sliceUnsigned64 reads a version-aware unsigned 64-bit field at off
// within p, advancing off and wrapping any error with what for context.
func sliceUnsigned64(r intReader, p []byte, off int, what string) (int64, int, error) {
	v, n, err := r.u64(p, off)
	if err != nil {
		return 0, off, fmt.Errorf("cram: %s: %w", what, err)
	}
	return v, off + n, nil
}

// sliceWide reads a field that is 64-bit for CRAM v4 and 32-bit for v2/v3,
// returning it as an int64. The slice alignment start and span widened in
// v4; this keeps the v2/v3 32-bit decode unchanged while reading the v4
// field at its full width.
func sliceWide(r intReader, p []byte, off int, major uint8, what string) (int64, int, error) {
	if major >= 4 {
		return sliceUnsigned64(r, p, off, what)
	}
	v, n, err := sliceUnsignedInt(r, p, off, what)
	return int64(v), n, err
}
