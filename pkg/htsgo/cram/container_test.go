package cram

import (
	"bytes"
	"encoding/binary"
	"hash/crc32"
	"testing"
)

// encITF8 encodes v as an ITF-8 byte sequence. It is a test helper that
// mirrors the readITF8 decoder, used to build synthetic CRAM structures.
func encITF8(v int32) []byte {
	u := uint32(v)
	switch {
	case u < 1<<7:
		return []byte{byte(u)}
	case u < 1<<14:
		return []byte{byte(u>>8) | 0x80, byte(u)}
	case u < 1<<21:
		return []byte{byte(u>>16) | 0xc0, byte(u >> 8), byte(u)}
	case u < 1<<28:
		return []byte{byte(u>>24) | 0xe0, byte(u >> 16), byte(u >> 8), byte(u)}
	default:
		return []byte{byte(u>>28) | 0xf0, byte(u >> 20), byte(u >> 12), byte(u >> 4), byte(u) & 0x0f}
	}
}

// encLTF8 encodes v as an LTF-8 byte sequence. It is a test helper that
// mirrors the readLTF8 decoder.
func encLTF8(v int64) []byte {
	u := uint64(v)
	switch {
	case u < 1<<7:
		return []byte{byte(u)}
	case u < 1<<14:
		return []byte{byte(u>>8) | 0x80, byte(u)}
	case u < 1<<21:
		return []byte{byte(u>>16) | 0xc0, byte(u >> 8), byte(u)}
	case u < 1<<28:
		return []byte{byte(u>>24) | 0xe0, byte(u >> 16), byte(u >> 8), byte(u)}
	case u < 1<<35:
		return []byte{byte(u>>32) | 0xf0, byte(u >> 24), byte(u >> 16), byte(u >> 8), byte(u)}
	case u < 1<<42:
		return []byte{byte(u>>40) | 0xf8, byte(u >> 32), byte(u >> 24), byte(u >> 16), byte(u >> 8), byte(u)}
	case u < 1<<49:
		return []byte{byte(u>>48) | 0xfc, byte(u >> 40), byte(u >> 32), byte(u >> 24), byte(u >> 16), byte(u >> 8), byte(u)}
	case u < 1<<56:
		return []byte{0xfe, byte(u >> 48), byte(u >> 40), byte(u >> 32), byte(u >> 24), byte(u >> 16), byte(u >> 8), byte(u)}
	default:
		return []byte{0xff, byte(u >> 56), byte(u >> 48), byte(u >> 40), byte(u >> 32), byte(u >> 24), byte(u >> 16), byte(u >> 8), byte(u)}
	}
}

// appendCRC appends the little-endian IEEE CRC32 of b to b.
func appendCRC(b []byte) []byte {
	var crc [4]byte
	binary.LittleEndian.PutUint32(crc[:], crc32.ChecksumIEEE(b))
	return append(b, crc[:]...)
}

// buildContainerHeader assembles a CRAM v3 container header on the wire
// from the given field values, with the trailing CRC32.
func buildContainerHeader(h ContainerHeader) []byte {
	var buf bytes.Buffer
	var lenBuf [4]byte
	binary.LittleEndian.PutUint32(lenBuf[:], uint32(h.Length))
	buf.Write(lenBuf[:])
	buf.Write(encITF8(h.RefSeqID))
	buf.Write(encITF8(h.StartPos))
	buf.Write(encITF8(h.AlignmentSpan))
	buf.Write(encITF8(h.NumRecords))
	buf.Write(encLTF8(h.RecordCounter))
	buf.Write(encLTF8(h.NumBases))
	buf.Write(encITF8(h.NumBlocks))
	buf.Write(encITF8(int32(len(h.Landmarks))))
	for _, lm := range h.Landmarks {
		buf.Write(encITF8(lm))
	}
	return appendCRC(buf.Bytes())
}

// TestContainerHeaderRoundTrip builds a container header, reads it back
// and checks every field plus the CRC32.
func TestContainerHeaderRoundTrip(t *testing.T) {
	want := ContainerHeader{
		Length:        12345,
		RefSeqID:      -2,
		StartPos:      0,
		AlignmentSpan: 0,
		NumRecords:    7,
		RecordCounter: 8,
		NumBases:      162,
		NumBlocks:     18,
		Landmarks:     []int32{0, 189, 4096},
	}
	wire := buildContainerHeader(want)
	got, err := parseContainerHeader(bytes.NewReader(wire), v3Def)
	if err != nil {
		t.Fatalf("parseContainerHeader: %v", err)
	}
	if got.Length != want.Length || got.RefSeqID != want.RefSeqID ||
		got.NumRecords != want.NumRecords || got.RecordCounter != want.RecordCounter ||
		got.NumBases != want.NumBases || got.NumBlocks != want.NumBlocks {
		t.Errorf("scalar field mismatch: got %+v want %+v", got, want)
	}
	if len(got.Landmarks) != len(want.Landmarks) {
		t.Fatalf("landmark count: got %d want %d", len(got.Landmarks), len(want.Landmarks))
	}
	for i := range want.Landmarks {
		if got.Landmarks[i] != want.Landmarks[i] {
			t.Errorf("landmark %d: got %d want %d", i, got.Landmarks[i], want.Landmarks[i])
		}
	}
}

// TestContainerHeaderCRCMismatch corrupts a container header's CRC and
// checks it is rejected.
func TestContainerHeaderCRCMismatch(t *testing.T) {
	wire := buildContainerHeader(ContainerHeader{Length: 10, NumBlocks: 1})
	wire[len(wire)-2] ^= 0xff
	if _, err := parseContainerHeader(bytes.NewReader(wire), v3Def); err == nil {
		t.Errorf("expected CRC32 mismatch error")
	}
}

// TestContainerHeaderTruncated checks parseContainerHeader errors —
// never panics — on input truncated at each byte offset.
func TestContainerHeaderTruncated(t *testing.T) {
	full := buildContainerHeader(ContainerHeader{
		Length: 5, NumBlocks: 2, Landmarks: []int32{0, 1},
	})
	for n := 0; n < len(full); n++ {
		if _, err := parseContainerHeader(bytes.NewReader(full[:n]), v3Def); err == nil {
			t.Errorf("expected error for container header truncated to %d bytes", n)
		}
	}
}

// TestCRCReaderCount checks the crcReader tracks the byte count and CRC
// of what it forwards.
func TestCRCReaderCount(t *testing.T) {
	data := []byte("the quick brown fox")
	cr := newCRCReader(bytes.NewReader(data))
	buf := make([]byte, len(data))
	if _, err := readFull(cr, buf); err != nil {
		t.Fatalf("read: %v", err)
	}
	if cr.count != int64(len(data)) {
		t.Errorf("count = %d, want %d", cr.count, len(data))
	}
	if cr.sum() != crc32.ChecksumIEEE(data) {
		t.Errorf("crc = %#x, want %#x", cr.sum(), crc32.ChecksumIEEE(data))
	}
}

// readFull is a tiny io.ReadFull stand-in to avoid an import in the test
// above.
func readFull(r interface{ Read([]byte) (int, error) }, buf []byte) (int, error) {
	total := 0
	for total < len(buf) {
		n, err := r.Read(buf[total:])
		total += n
		if err != nil {
			return total, err
		}
	}
	return total, nil
}
