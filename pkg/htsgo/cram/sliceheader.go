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
	// first record.
	AlignmentStart int32
	// AlignmentSpan is the number of reference bases the slice covers.
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
// payload is a sequence of ITF-8/LTF-8 fields: ref id, alignment start,
// alignment span, record count, record counter, block count, an ITF-8
// array of block content ids, the embedded-reference block id, the
// 16-byte reference MD5, and any optional trailing tag bytes.
func parseSliceHeader(p []byte) (*SliceHeader, error) {
	s := &SliceHeader{}
	off := 0
	var err error
	if s.RefSeqID, off, err = sliceITF8(p, off, "slice ref seq id"); err != nil {
		return nil, err
	}
	if s.AlignmentStart, off, err = sliceITF8(p, off, "slice alignment start"); err != nil {
		return nil, err
	}
	if s.AlignmentSpan, off, err = sliceITF8(p, off, "slice alignment span"); err != nil {
		return nil, err
	}
	if s.NumRecords, off, err = sliceITF8(p, off, "slice record count"); err != nil {
		return nil, err
	}
	if s.RecordCounter, off, err = sliceLTF8(p, off, "slice record counter"); err != nil {
		return nil, err
	}
	if s.NumBlocks, off, err = sliceITF8(p, off, "slice block count"); err != nil {
		return nil, err
	}
	if s.NumBlocks < 0 {
		return nil, fmt.Errorf("cram: slice header declares negative block count %d", s.NumBlocks)
	}
	var nIDs int32
	if nIDs, off, err = sliceITF8(p, off, "slice block content id count"); err != nil {
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
		if id, off, err = sliceITF8(p, off, "slice block content id"); err != nil {
			return nil, err
		}
		s.BlockContentIDs = append(s.BlockContentIDs, id)
	}
	if s.EmbeddedRefID, off, err = sliceITF8(p, off, "slice embedded reference id"); err != nil {
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

// sliceITF8 reads an ITF-8 field at off within p, advancing off and
// wrapping any error with what for context.
func sliceITF8(p []byte, off int, what string) (int32, int, error) {
	v, n, err := itf8At(p, off)
	if err != nil {
		return 0, off, fmt.Errorf("cram: %s: %w", what, err)
	}
	return v, off + n, nil
}

// sliceLTF8 reads an LTF-8 field at off within p, advancing off and
// wrapping any error with what for context.
func sliceLTF8(p []byte, off int, what string) (int64, int, error) {
	v, n, err := ltf8At(p, off)
	if err != nil {
		return 0, off, fmt.Errorf("cram: %s: %w", what, err)
	}
	return v, off + n, nil
}
