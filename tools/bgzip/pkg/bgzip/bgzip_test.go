package bgzip

import (
	"bytes"
	"compress/flate"
	"crypto/rand"
	"encoding/binary"
	"errors"
	"hash/crc32"
	"io"
	"strings"
	"testing"
)

// roundTrip compresses payload with NewWriter, then decompresses the result
// with NewReader and returns the decoded bytes.
func roundTrip(t *testing.T, payload []byte) []byte {
	t.Helper()
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write(payload); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	r, err := NewReader(&buf)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	out, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if err := r.Close(); err != nil {
		t.Fatalf("Reader.Close: %v", err)
	}
	return out
}

func TestRoundTrip(t *testing.T) {
	cases := []struct {
		name string
		size int
	}{
		{"empty", 0},
		{"tiny", 5},
		{"sub-block", 1024},
		{"exact-block", MaxBlockSize},
		{"just-over-block", MaxBlockSize + 1},
		{"two-blocks", 2 * MaxBlockSize},
		{"multi-block-uneven", 3*MaxBlockSize + 17},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			payload := make([]byte, tc.size)
			if _, err := rand.Read(payload); err != nil {
				t.Fatal(err)
			}
			out := roundTrip(t, payload)
			if !bytes.Equal(out, payload) {
				t.Fatalf("round trip mismatch: got %d bytes, want %d bytes", len(out), len(payload))
			}
		})
	}
}

func TestRoundTripTextual(t *testing.T) {
	payload := []byte(strings.Repeat("chr1\t1000\t2000\tfeatureA\n", 5000))
	out := roundTrip(t, payload)
	if !bytes.Equal(out, payload) {
		t.Fatalf("round trip mismatch for textual data")
	}
}

func TestEOFBlockAlwaysWritten(t *testing.T) {
	cases := [][]byte{
		nil,
		[]byte("hi"),
		bytes.Repeat([]byte{0xAA}, MaxBlockSize+10),
	}
	for i, payload := range cases {
		var buf bytes.Buffer
		w := NewWriter(&buf)
		if len(payload) > 0 {
			if _, err := w.Write(payload); err != nil {
				t.Fatalf("case %d: Write: %v", i, err)
			}
		}
		if err := w.Close(); err != nil {
			t.Fatalf("case %d: Close: %v", i, err)
		}
		got := buf.Bytes()
		if len(got) < len(EOFBlock) {
			t.Fatalf("case %d: stream shorter than EOF block (%d < %d)", i, len(got), len(EOFBlock))
		}
		tail := got[len(got)-len(EOFBlock):]
		if !bytes.Equal(tail, EOFBlock) {
			t.Fatalf("case %d: stream does not end with EOF block. tail=%x", i, tail)
		}
	}
}

func TestEmptyStreamHasOnlyEOFBlock(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(buf.Bytes(), EOFBlock) {
		t.Fatalf("empty stream != EOFBlock: %x", buf.Bytes())
	}
}

func TestWriterLevels(t *testing.T) {
	payload := []byte(strings.Repeat("ACGT", 2000))
	levels := []int{
		flate.NoCompression,
		flate.BestSpeed,
		flate.DefaultCompression,
		flate.BestCompression,
	}
	for _, lvl := range levels {
		var buf bytes.Buffer
		w, err := NewWriterLevel(&buf, lvl)
		if err != nil {
			t.Fatalf("NewWriterLevel(%d): %v", lvl, err)
		}
		if _, err := w.Write(payload); err != nil {
			t.Fatal(err)
		}
		if err := w.Close(); err != nil {
			t.Fatal(err)
		}
		r, _ := NewReader(&buf)
		out, err := io.ReadAll(r)
		if err != nil {
			t.Fatalf("level %d read: %v", lvl, err)
		}
		if !bytes.Equal(out, payload) {
			t.Fatalf("level %d round trip mismatch", lvl)
		}
	}
}

func TestWriterInvalidLevel(t *testing.T) {
	if _, err := NewWriterLevel(io.Discard, 99); err == nil {
		t.Fatal("expected error for invalid level")
	}
}

func TestWriterWriteAfterClose(t *testing.T) {
	w := NewWriter(io.Discard)
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	if _, err := w.Write([]byte("x")); err == nil {
		t.Fatal("expected error writing to closed writer")
	}
}

func TestWriterFlush(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write([]byte("hello")); err != nil {
		t.Fatal(err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}
	beforeClose := buf.Len()
	if beforeClose == 0 {
		t.Fatal("Flush did not emit any bytes")
	}
	// A second Flush with no buffered data is a no-op.
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}
	if buf.Len() != beforeClose {
		t.Fatal("Flush emitted bytes for empty buffer")
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}

func TestCloseIdempotent(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("second Close returned %v", err)
	}
}

func TestBlocksHaveBCSubfieldAndCorrectBSIZE(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	// Force three blocks: one full, one full, one partial.
	payload := make([]byte, 2*MaxBlockSize+100)
	for i := range payload {
		payload[i] = byte(i & 0xff)
	}
	if _, err := w.Write(payload); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	offsets, err := Scan(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("Scan: %v", err)
	}
	if len(offsets) != 3 {
		t.Fatalf("expected 3 blocks, got %d", len(offsets))
	}
	if offsets[0].UncompressedSize != MaxBlockSize {
		t.Fatalf("block 0 size = %d, want %d", offsets[0].UncompressedSize, MaxBlockSize)
	}
	if offsets[1].UncompressedSize != MaxBlockSize {
		t.Fatalf("block 1 size = %d, want %d", offsets[1].UncompressedSize, MaxBlockSize)
	}
	if offsets[2].UncompressedSize != 100 {
		t.Fatalf("block 2 size = %d, want 100", offsets[2].UncompressedSize)
	}

	// Cumulative compressed offsets must be monotonic and match the BSIZE of
	// each block, and each block must carry the BC subfield (verified by Scan
	// returning successfully — readBlockHeader requires it).
	total := int64(0)
	for i, b := range offsets {
		if b.CompressedOffset != total {
			t.Fatalf("block %d: CompressedOffset = %d, want %d", i, b.CompressedOffset, total)
		}
		total += int64(b.CompressedSize)
	}
}

func TestDecompressedSize(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	payload := []byte(strings.Repeat("X", 3*MaxBlockSize+42))
	if _, err := w.Write(payload); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	n, err := DecompressedSize(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("DecompressedSize: %v", err)
	}
	if n != int64(len(payload)) {
		t.Fatalf("DecompressedSize = %d, want %d", n, len(payload))
	}
}

func TestUncompressedOffsetAt(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	payload := make([]byte, 3*MaxBlockSize+42)
	if _, err := w.Write(payload); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	data := buf.Bytes()

	// 0 → 0
	got, err := UncompressedOffsetAt(bytes.NewReader(data), 0)
	if err != nil || got != 0 {
		t.Fatalf("offset 0 → %d, %v", got, err)
	}

	// Read all blocks to know boundary positions.
	offsets, err := Scan(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("Scan: %v", err)
	}
	for _, b := range offsets {
		u, err := UncompressedOffsetAt(bytes.NewReader(data), b.CompressedOffset)
		if err != nil {
			t.Fatalf("UncompressedOffsetAt(%d): %v", b.CompressedOffset, err)
		}
		if u != b.UncompressedOffset {
			t.Fatalf("UncompressedOffsetAt(%d) = %d, want %d", b.CompressedOffset, u, b.UncompressedOffset)
		}
	}

	// Past the last block → total decompressed size.
	last := offsets[len(offsets)-1]
	farOff := last.CompressedOffset + int64(last.CompressedSize) + 999
	u, err := UncompressedOffsetAt(bytes.NewReader(data), farOff)
	if err != nil {
		t.Fatalf("past-end offset: %v", err)
	}
	if u != int64(len(payload)) {
		t.Fatalf("past-end uncompressed = %d, want %d", u, len(payload))
	}

	// Offset that does not align to a block boundary → error.
	if _, err := UncompressedOffsetAt(bytes.NewReader(data), 1); err == nil {
		t.Fatal("expected error for misaligned offset")
	}

	// Negative offset → error.
	if _, err := UncompressedOffsetAt(bytes.NewReader(data), -1); err == nil {
		t.Fatal("expected error for negative offset")
	}
}

func TestGZIRoundTrip(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write(bytes.Repeat([]byte("ACGT"), 70000)); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	offsets, err := Scan(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("Scan: %v", err)
	}

	var gzi bytes.Buffer
	if err := WriteGZI(&gzi, offsets); err != nil {
		t.Fatalf("WriteGZI: %v", err)
	}
	// The expected count drops the leading (0,0) entry.
	wantCount := uint64(len(offsets) - 1)
	gotCount := binary.LittleEndian.Uint64(gzi.Bytes()[:8])
	if gotCount != wantCount {
		t.Fatalf("gzi count = %d, want %d", gotCount, wantCount)
	}
	if gzi.Len() != int(8+wantCount*16) {
		t.Fatalf("gzi length = %d, want %d", gzi.Len(), 8+wantCount*16)
	}

	parsed, err := ReadGZI(bytes.NewReader(gzi.Bytes()))
	if err != nil {
		t.Fatalf("ReadGZI: %v", err)
	}
	if len(parsed) != int(wantCount) {
		t.Fatalf("ReadGZI entries = %d, want %d", len(parsed), wantCount)
	}
	for i, p := range parsed {
		want := offsets[i+1]
		if p.CompressedOffset != want.CompressedOffset || p.UncompressedOffset != want.UncompressedOffset {
			t.Fatalf("entry %d: got (%d,%d) want (%d,%d)",
				i, p.CompressedOffset, p.UncompressedOffset,
				want.CompressedOffset, want.UncompressedOffset)
		}
	}
}

func TestReadGZIShortInput(t *testing.T) {
	if _, err := ReadGZI(bytes.NewReader([]byte{0x01})); err == nil {
		t.Fatal("expected error for truncated gzi header")
	}
	// Header announces 2 entries but only one entry follows.
	var b bytes.Buffer
	var hdr [8]byte
	binary.LittleEndian.PutUint64(hdr[:], 2)
	b.Write(hdr[:])
	b.Write(make([]byte, 16))
	if _, err := ReadGZI(&b); err == nil {
		t.Fatal("expected error for truncated gzi body")
	}
}

func TestDecodeKnownGoodFixture(t *testing.T) {
	// Build a 2-block fixture by hand and decode it, exercising the
	// "decompress a known-good external BGZF file" path. The first block
	// contains "hello"; the second contains "world".
	fixture := buildFixture(t, [][]byte{[]byte("hello"), []byte("world")})
	r, err := NewReader(bytes.NewReader(fixture))
	if err != nil {
		t.Fatal(err)
	}
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if string(got) != "helloworld" {
		t.Fatalf("got %q, want %q", got, "helloworld")
	}
}

func TestReadTruncatedMissingEOF(t *testing.T) {
	// Build a valid stream then chop off the EOF block.
	var buf bytes.Buffer
	w := NewWriter(&buf)
	w.Write([]byte("payload"))
	w.Close()
	truncated := buf.Bytes()[:buf.Len()-len(EOFBlock)]
	r, _ := NewReader(bytes.NewReader(truncated))
	_, err := io.ReadAll(r)
	if !errors.Is(err, ErrTruncated) {
		t.Fatalf("got err=%v, want ErrTruncated", err)
	}
}

func TestReadTruncatedMidBlock(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	w.Write([]byte("payload"))
	w.Close()
	// Chop in the middle of the first block (after the header).
	truncated := buf.Bytes()[:14]
	r, _ := NewReader(bytes.NewReader(truncated))
	_, err := io.ReadAll(r)
	if err == nil {
		t.Fatal("expected error reading truncated header")
	}
}

func TestReadBadMagic(t *testing.T) {
	// 18 bytes of garbage — enough to fool readBlockHeader into trying.
	garbage := make([]byte, 18)
	r, _ := NewReader(bytes.NewReader(garbage))
	_, err := io.ReadAll(r)
	if !errors.Is(err, ErrBadMagic) {
		t.Fatalf("got %v, want ErrBadMagic", err)
	}
}

func TestReadMissingBCSubfield(t *testing.T) {
	// Hand-build a gzip block with FEXTRA set but with a non-BC subfield.
	var b bytes.Buffer
	b.Write([]byte{0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff})
	// XLEN = 6, subfield SI=AA, SLEN=2, value=0
	binary.Write(&b, binary.LittleEndian, uint16(6))
	b.Write([]byte{'A', 'A', 0x02, 0x00, 0x00, 0x00})
	// Minimal empty deflate stream + crc + isize=0
	emptyDeflate := emptyDeflateBytes()
	b.Write(emptyDeflate)
	b.Write(make([]byte, 8))

	r, _ := NewReader(bytes.NewReader(b.Bytes()))
	_, err := io.ReadAll(r)
	if !errors.Is(err, ErrNoBCSubfield) {
		t.Fatalf("got %v, want ErrNoBCSubfield", err)
	}
}

func TestReadMissingFEXTRA(t *testing.T) {
	// Real gzip member with FEXTRA flag cleared.
	var b bytes.Buffer
	b.Write([]byte{0x1f, 0x8b, 0x08, 0x00, 0, 0, 0, 0, 0, 0xff, 0, 0})
	r, _ := NewReader(bytes.NewReader(b.Bytes()))
	_, err := io.ReadAll(r)
	if !errors.Is(err, ErrNoExtra) {
		t.Fatalf("got %v, want ErrNoExtra", err)
	}
}

func TestReadCorruptedCRC(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	w.Write([]byte("hello"))
	w.Close()
	data := buf.Bytes()
	// Find the trailing 8-byte footer of the first block (before the 28-byte
	// EOF block). Flip a bit in the CRC.
	footerStart := len(data) - len(EOFBlock) - 8
	data[footerStart] ^= 0xff
	r, _ := NewReader(bytes.NewReader(data))
	_, err := io.ReadAll(r)
	if !errors.Is(err, ErrChecksum) {
		t.Fatalf("got %v, want ErrChecksum", err)
	}
}

func TestReaderWriterReuse(t *testing.T) {
	// Multi-block read where the underlying Reader is exercised across more
	// than one block transition. Ensures flate.Resetter path actually runs.
	payload := make([]byte, 4*MaxBlockSize+5)
	for i := range payload {
		payload[i] = byte((i * 7) & 0xff)
	}
	got := roundTrip(t, payload)
	if !bytes.Equal(got, payload) {
		t.Fatal("multi-block round trip mismatch")
	}
}

func TestScanReturnsTruncated(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	w.Write([]byte("payload"))
	w.Close()
	truncated := buf.Bytes()[:buf.Len()-len(EOFBlock)]
	_, err := Scan(bytes.NewReader(truncated))
	if !errors.Is(err, ErrTruncated) {
		t.Fatalf("got %v, want ErrTruncated", err)
	}
}

// emptyDeflateBytes returns the 2-byte fixed-Huffman empty deflate stream that
// htslib uses for the BGZF EOF block (0x03 0x00).
func emptyDeflateBytes() []byte {
	return []byte{0x03, 0x00}
}

// buildFixture hand-builds a multi-block BGZF stream by compressing each
// payload with flate and wrapping it in the BGZF framing. It does NOT use the
// Writer in this package, so it can be used to verify decoding independently.
func buildFixture(t *testing.T, payloads [][]byte) []byte {
	t.Helper()
	var out bytes.Buffer
	for _, payload := range payloads {
		var defl bytes.Buffer
		fw, err := flate.NewWriter(&defl, flate.DefaultCompression)
		if err != nil {
			t.Fatal(err)
		}
		if _, err := fw.Write(payload); err != nil {
			t.Fatal(err)
		}
		if err := fw.Close(); err != nil {
			t.Fatal(err)
		}
		blockLen := 12 + 6 + defl.Len() + 8

		var hdr [18]byte
		hdr[0] = 0x1f
		hdr[1] = 0x8b
		hdr[2] = 8
		hdr[3] = 0x04
		hdr[9] = 0xff
		binary.LittleEndian.PutUint16(hdr[10:12], 6)
		hdr[12] = 'B'
		hdr[13] = 'C'
		binary.LittleEndian.PutUint16(hdr[14:16], 2)
		binary.LittleEndian.PutUint16(hdr[16:18], uint16(blockLen-1))
		out.Write(hdr[:])
		out.Write(defl.Bytes())

		var footer [8]byte
		binary.LittleEndian.PutUint32(footer[0:4], crc32.ChecksumIEEE(payload))
		binary.LittleEndian.PutUint32(footer[4:8], uint32(len(payload)))
		out.Write(footer[:])
	}
	out.Write(EOFBlock)
	return out.Bytes()
}

func TestFixtureMatchesWriter(t *testing.T) {
	// Both routes should accept each other's output.
	payload := []byte("BGZF interop test 1234567890")
	fixture := buildFixture(t, [][]byte{payload})

	r, _ := NewReader(bytes.NewReader(fixture))
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(got, payload) {
		t.Fatalf("got %q, want %q", got, payload)
	}
}
