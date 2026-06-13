package cram

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"hash/crc32"
	"io"
	"testing"
)

// buildV4ContainerHeader assembles a CRAM v4 container header (uint7
// fields plus the trailing CRC32) from the given fields, mirroring the
// bytes htslib's cram_write_container writes for major>=4.
func buildV4ContainerHeader(length int32, refID int32, start, span int64, nrec int32, counter, nbases int64, nblocks int32, landmarks []int32) []byte {
	var b []byte
	b = appendUint7(b, uint64(uint32(length)))
	b = appendSint7(b, int64(refID))
	b = appendUint7(b, uint64(start))
	b = appendUint7(b, uint64(span))
	b = appendUint7(b, uint64(uint32(nrec)))
	b = appendUint7(b, uint64(counter))
	b = appendUint7(b, uint64(nbases))
	b = appendUint7(b, uint64(uint32(nblocks)))
	b = appendUint7(b, uint64(len(landmarks)))
	for _, lm := range landmarks {
		b = appendUint7(b, uint64(uint32(lm)))
	}
	crc := crc32.ChecksumIEEE(b)
	var crcBytes [4]byte
	binary.LittleEndian.PutUint32(crcBytes[:], crc)
	return append(b, crcBytes[:]...)
}

// TestParseV4ContainerHeader parses a synthetic v4 container header and
// checks every field round-trips, including the signed ref id and 64-bit
// alignment coordinates.
func TestParseV4ContainerHeader(t *testing.T) {
	hdr := buildV4ContainerHeader(509, 0, 12345, 678, 7, 100, 4000, 2, []int32{0, 332})
	def := FileDefinition{Major: 4, Minor: 0}
	h, err := parseContainerHeader(bytes.NewReader(hdr), def)
	if err != nil {
		t.Fatalf("parseContainerHeader (v4): %v", err)
	}
	if h.Length != 509 || h.RefSeqID != 0 || h.StartPos != 12345 || h.AlignmentSpan != 678 ||
		h.NumRecords != 7 || h.RecordCounter != 100 || h.NumBases != 4000 || h.NumBlocks != 2 {
		t.Fatalf("v4 header parsed wrong: %+v", h)
	}
	if len(h.Landmarks) != 2 || h.Landmarks[0] != 0 || h.Landmarks[1] != 332 {
		t.Fatalf("v4 landmarks = %v, want [0 332]", h.Landmarks)
	}
}

// TestParseV4ContainerHeaderNegativeRef confirms a v4 container with a
// signed ref id of -1 (unmapped) round-trips via the zig-zag varint, where
// the v3 ITF-8 reader would mis-interpret the sign.
func TestParseV4ContainerHeaderNegativeRef(t *testing.T) {
	hdr := buildV4ContainerHeader(15, -1, 0x454f46, 0, 0, 0, 0, 1, nil)
	def := FileDefinition{Major: 4, Minor: 0}
	h, err := parseContainerHeader(bytes.NewReader(hdr), def)
	if err != nil {
		t.Fatalf("parseContainerHeader (v4, ref -1): %v", err)
	}
	if h.RefSeqID != -1 {
		t.Errorf("ref seq id = %d, want -1", h.RefSeqID)
	}
	if h.StartPos != 0x454f46 {
		t.Errorf("start pos = %d, want 0x454f46", h.StartPos)
	}
}

// TestV4EOFMarkerRecognised confirms the captured 31-byte v4 EOF marker is
// recognised as the end-of-stream sentinel and reported via io.EOF, not
// parsed as a real container. The bytes are the marker upstream htslib
// writes for v4.0 (see container.go eofMarkerV4).
func TestV4EOFMarkerRecognised(t *testing.T) {
	if len(eofMarkerV4) != 31 {
		t.Fatalf("v4 EOF marker is %d bytes, want 31", len(eofMarkerV4))
	}
	def := FileDefinition{Major: 4, Minor: 0}
	br := bufio.NewReader(bytes.NewReader(eofMarkerV4))
	h, err := readContainerHeader(br, def)
	if err != nil {
		t.Fatalf("readContainerHeader on v4 EOF marker: %v", err)
	}
	if !h.IsEOF {
		t.Fatal("v4 EOF marker was not recognised as the EOF container")
	}
}

// TestV4EOFMarkerHeaderValid confirms the v4 EOF marker's container-header
// portion parses to the documented empty-container values (ref id -1,
// start 0x454f46 "EOF"), so the marker is structurally a real—if empty—
// container as htslib intends.
func TestV4EOFMarkerHeaderValid(t *testing.T) {
	def := FileDefinition{Major: 4, Minor: 0}
	h, err := parseContainerHeader(bytes.NewReader(eofMarkerV4), def)
	if err != nil {
		t.Fatalf("parsing the v4 EOF marker as a container header: %v", err)
	}
	if h.Length != 15 || h.RefSeqID != -1 || h.StartPos != 0x454f46 || h.NumBlocks != 1 {
		t.Fatalf("v4 EOF marker header = %+v; want length 15, ref -1, start 0x454f46, 1 block", h)
	}
}

// TestV4ReaderStopsAtEOF confirms a minimal v4 stream (file definition then
// EOF marker) is walked to a clean io.EOF without spurious errors.
func TestV4ReaderStopsAtEOF(t *testing.T) {
	var buf bytes.Buffer
	def := make([]byte, fileDefSize)
	copy(def, "CRAM")
	def[4] = 4
	def[5] = 0
	buf.Write(def)
	buf.Write(eofMarkerV4)

	rd, err := NewReader(&buf)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	if _, err := rd.Next(); err != io.EOF {
		t.Fatalf("Next on empty v4 stream = %v, want io.EOF", err)
	}
}
