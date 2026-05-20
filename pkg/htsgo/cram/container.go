package cram

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"fmt"
	"hash/crc32"
	"io"
)

// eofMarkerV3 is the 38-byte sentinel container that terminates a
// well-formed CRAM v3 file. htslib writes it verbatim and a reader
// recognises it byte-for-byte rather than parsing it as a real
// container.
var eofMarkerV3 = []byte{
	0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0,
	0x45, 0x4f, 0x46, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x05,
	0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00,
	0x01, 0x00, 0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
}

// eofMarkerV2 is the 30-byte sentinel container that terminates a
// well-formed CRAM v2 file. It lacks the trailing CRC32 that v3 adds.
var eofMarkerV2 = []byte{
	0x0b, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0,
	0x45, 0x4f, 0x46, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
	0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00,
}

// ContainerHeader holds the parsed fields of a CRAM container header.
// The header is followed by the container's blocks; the first block of
// the first container carries the SAM text header, and every later
// container holds a compression-header block followed by one or more
// slices.
type ContainerHeader struct {
	// Length is the total byte size of all blocks that follow this
	// header (it does not include the header itself).
	Length int32
	// RefSeqID is the reference sequence identifier the container's
	// records align to: a value >= 0 indexes the SAM header @SQ lines,
	// -1 means unmapped/unplaced reads and -2 means the container holds
	// reads spanning multiple references.
	RefSeqID int32
	// StartPos is the alignment start (1-based) of the container's first
	// record on RefSeqID.
	StartPos int32
	// AlignmentSpan is the number of reference bases the container's
	// records span.
	AlignmentSpan int32
	// NumRecords is the count of alignment records across all the
	// container's slices.
	NumRecords int32
	// RecordCounter is the running total of records in all containers
	// before this one.
	RecordCounter int64
	// NumBases is the total number of read bases in the container.
	NumBases int64
	// NumBlocks is the count of blocks that follow this header.
	NumBlocks int32
	// Landmarks holds the byte offsets, relative to the end of the
	// container header, of each slice within the container.
	Landmarks []int32
	// CRC is the 4-byte little-endian CRC32 stored after the header in
	// CRAM v3+; it is zero and unused for v2.
	CRC uint32
	// IsEOF reports whether this header is the CRAM EOF marker container
	// rather than a container with real data.
	IsEOF bool
}

// crcReader wraps an io.Reader and accumulates an IEEE CRC32 over every
// byte it yields. It is used to checksum a CRAM container header or
// block header in a single streaming pass, so the trailing CRC32 field
// can be validated without buffering the whole structure twice.
type crcReader struct {
	r     io.Reader
	crc   uint32
	count int64
}

// newCRCReader returns a crcReader over r.
func newCRCReader(r io.Reader) *crcReader {
	return &crcReader{r: r}
}

// Read implements io.Reader, updating the running CRC32 with the bytes
// it returns.
func (c *crcReader) Read(p []byte) (int, error) {
	n, err := c.r.Read(p)
	if n > 0 {
		c.crc = crc32.Update(c.crc, crc32.IEEETable, p[:n])
		c.count += int64(n)
	}
	return n, err
}

// sum returns the CRC32 accumulated so far.
func (c *crcReader) sum() uint32 { return c.crc }

// readContainerHeader parses one CRAM container header from br. When the
// CRAM version embeds CRC32 fields it reads and validates the trailing
// checksum. A header whose bytes match the CRAM EOF marker is reported
// via ContainerHeader.IsEOF. It returns io.EOF when br is exhausted at a
// container boundary.
//
// The EOF-marker check peeks rather than consumes, so a real (and
// possibly much shorter) container header is not over-read; the parse
// then proceeds directly from the same buffered reader.
func readContainerHeader(br *bufio.Reader, def FileDefinition) (ContainerHeader, error) {
	marker := eofMarkerV3
	if !def.hasCRC() {
		marker = eofMarkerV2
	}
	// A zero-byte peek distinguishes a clean end-of-stream from a
	// truncated container.
	if _, err := br.Peek(1); err == io.EOF {
		return ContainerHeader{}, io.EOF
	}
	if head, err := br.Peek(len(marker)); err == nil && bytes.Equal(head, marker) {
		if _, derr := br.Discard(len(marker)); derr != nil {
			return ContainerHeader{}, fmt.Errorf("cram: consuming EOF marker: %w", derr)
		}
		return ContainerHeader{IsEOF: true}, nil
	}
	return parseContainerHeader(br, def)
}

// parseContainerHeader reads the variable-length container header fields
// from r and, for CRAM v3+, validates the trailing CRC32.
//
// Note on the length field: the CRAM specification text describes the
// container length as ITF-8, but htslib — and therefore every CRAM file
// in the wild — writes it as a fixed 4-byte little-endian int32 so the
// length can be back-patched once the container body is assembled. This
// parser follows htslib's on-disk reality. The CRAM v3 EOF marker, whose
// first four bytes are 0f 00 00 00 (= 15), confirms the fixed-width
// little-endian encoding.
func parseContainerHeader(r io.Reader, def FileDefinition) (ContainerHeader, error) {
	cr := newCRCReader(r)
	var h ContainerHeader
	var err error
	var lenBuf [4]byte
	if _, err = io.ReadFull(cr, lenBuf[:]); err != nil {
		return h, fmt.Errorf("cram: container length: %w",
			eofToUnexpected(err, "container length"))
	}
	h.Length = int32(binary.LittleEndian.Uint32(lenBuf[:]))
	if h.Length < 0 {
		return h, fmt.Errorf("cram: container declares negative length %d", h.Length)
	}
	if h.RefSeqID, _, err = readITF8(cr); err != nil {
		return h, fmt.Errorf("cram: container ref seq id: %w", err)
	}
	if h.StartPos, _, err = readITF8(cr); err != nil {
		return h, fmt.Errorf("cram: container start pos: %w", err)
	}
	if h.AlignmentSpan, _, err = readITF8(cr); err != nil {
		return h, fmt.Errorf("cram: container alignment span: %w", err)
	}
	if h.NumRecords, _, err = readITF8(cr); err != nil {
		return h, fmt.Errorf("cram: container record count: %w", err)
	}
	if h.RecordCounter, _, err = readLTF8(cr); err != nil {
		return h, fmt.Errorf("cram: container record counter: %w", err)
	}
	if h.NumBases, _, err = readLTF8(cr); err != nil {
		return h, fmt.Errorf("cram: container base count: %w", err)
	}
	if h.NumBlocks, _, err = readITF8(cr); err != nil {
		return h, fmt.Errorf("cram: container block count: %w", err)
	}
	if h.NumBlocks < 0 {
		return h, fmt.Errorf("cram: container declares negative block count %d", h.NumBlocks)
	}
	var nLandmarks int32
	if nLandmarks, _, err = readITF8(cr); err != nil {
		return h, fmt.Errorf("cram: container landmark count: %w", err)
	}
	if nLandmarks < 0 {
		return h, fmt.Errorf("cram: container declares negative landmark count %d", nLandmarks)
	}
	// Each landmark is at least one byte on the wire, so a landmark
	// count larger than the container's block region is corrupt. This
	// also bounds the slice growth against a malformed count on a
	// CRC-less CRAM v2 header.
	if int64(nLandmarks) > int64(h.Length) {
		return h, fmt.Errorf("cram: container declares %d landmarks but only %d body bytes",
			nLandmarks, h.Length)
	}
	// Landmarks grows by append rather than being pre-sized: nLandmarks
	// is untrusted and a pre-sized make would risk a huge allocation.
	for i := int32(0); i < nLandmarks; i++ {
		var lm int32
		if lm, _, err = readITF8(cr); err != nil {
			return h, fmt.Errorf("cram: container landmark %d: %w", i, err)
		}
		h.Landmarks = append(h.Landmarks, lm)
	}
	if def.hasCRC() {
		want := cr.sum()
		var crcBuf [4]byte
		if _, err = io.ReadFull(r, crcBuf[:]); err != nil {
			return h, fmt.Errorf("cram: container header CRC32: %w",
				eofToUnexpected(err, "container CRC"))
		}
		h.CRC = binary.LittleEndian.Uint32(crcBuf[:])
		if h.CRC != want {
			return h, fmt.Errorf("cram: container header CRC32 mismatch: stored %#08x, computed %#08x",
				h.CRC, want)
		}
	}
	return h, nil
}
