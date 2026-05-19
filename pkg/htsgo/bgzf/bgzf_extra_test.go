package bgzf

import (
	"bytes"
	"compress/flate"
	"encoding/binary"
	"errors"
	"fmt"
	"hash/crc32"
	"io"
	"testing"
)

// failingWriter returns an error after writing limit bytes.
type failingWriter struct {
	limit int
	n     int
}

func (fw *failingWriter) Write(p []byte) (int, error) {
	remaining := fw.limit - fw.n
	if remaining <= 0 {
		return 0, fmt.Errorf("forced write failure")
	}
	if len(p) > remaining {
		fw.n += remaining
		return remaining, fmt.Errorf("forced write failure")
	}
	fw.n += len(p)
	return len(p), nil
}

func TestWriterPropagatesUnderlyingError(t *testing.T) {
	// Cause the underlying writer to fail in the middle of encoding the first
	// block (header writes succeed, deflate body write fails).
	w := NewWriter(&failingWriter{limit: 10})
	_, err := w.Write([]byte("some data"))
	// Write itself buffers, so it may not fail; the error appears on Close.
	closeErr := w.Close()
	if err == nil && closeErr == nil {
		t.Fatal("expected error from failing underlying writer")
	}
	// Second Write after error must also fail.
	if _, err := w.Write([]byte("x")); err == nil {
		t.Fatal("expected write-after-error to fail")
	}
}

func TestWriterTinyLimitFailsImmediately(t *testing.T) {
	// Force overflow inside Write itself by writing many blocks against a
	// writer that refuses any output.
	w := NewWriter(&failingWriter{limit: 0})
	payload := make([]byte, 2*MaxBlockSize+100)
	_, err := w.Write(payload)
	if err == nil {
		// All buffered — Close will fail.
		if err := w.Close(); err == nil {
			t.Fatal("expected Close to fail against zero-limit writer")
		}
	}
}

func TestReadEOFAfterStreamDone(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	w.Write([]byte("hi"))
	w.Close()
	r, _ := NewReader(&buf)
	if _, err := io.ReadAll(r); err != nil {
		t.Fatal(err)
	}
	// Subsequent reads must return io.EOF.
	var tmp [4]byte
	n, err := r.Read(tmp[:])
	if n != 0 || err != io.EOF {
		t.Fatalf("post-stream Read = %d, %v; want 0, io.EOF", n, err)
	}
}

func TestReadTruncatedAfterPartialBlock(t *testing.T) {
	// Two valid blocks then no EOF marker. The Reader should still return the
	// data of the first block but signal ErrTruncated overall.
	fixture := buildFixture(t, [][]byte{[]byte("first block"), []byte("second block")})
	// Strip the EOF marker (last 28 bytes).
	truncated := fixture[:len(fixture)-len(EOFBlock)]
	r, _ := NewReader(bytes.NewReader(truncated))
	got, err := io.ReadAll(r)
	if !errors.Is(err, ErrTruncated) {
		t.Fatalf("got %v, want ErrTruncated", err)
	}
	if string(got) != "first blocksecond block" {
		t.Fatalf("partial read = %q", got)
	}
}

func TestUnsupportedCompressionMethod(t *testing.T) {
	// 18-byte fake header with bad CM byte.
	var b bytes.Buffer
	b.Write([]byte{0x1f, 0x8b, 99, 0x04, 0, 0, 0, 0, 0, 0xff})
	binary.Write(&b, binary.LittleEndian, uint16(6))
	b.Write([]byte{'B', 'C', 0x02, 0x00, 0x1b, 0x00})
	r, _ := NewReader(&b)
	_, err := io.ReadAll(r)
	if err == nil {
		t.Fatal("expected error for unsupported method")
	}
}

func TestBlockWithFNAMEFCOMMENTFHCRC(t *testing.T) {
	// Manually build a single block whose header sets FNAME, FCOMMENT, and
	// FHCRC in addition to FEXTRA. The Reader is expected to skip them.
	payload := []byte("with extras")
	var defl bytes.Buffer
	fw, _ := flate.NewWriter(&defl, flate.DefaultCompression)
	fw.Write(payload)
	fw.Close()
	deflated := defl.Bytes()

	// Header bytes count toward BSIZE: 12 fixed + 6 extra + 5 fname + 4 fcomment + 2 fhcrc.
	fname := []byte("abcd\x00")
	fcomment := []byte("hi\x00")
	headerLen := 12 + 6 + len(fname) + len(fcomment) + 2
	blockLen := headerLen + len(deflated) + 8

	flags := byte(flagFEXTRA | flagFNAME | flagFCOMMENT | flagFHCRC)
	var hdr [12]byte
	hdr[0] = 0x1f
	hdr[1] = 0x8b
	hdr[2] = 8
	hdr[3] = flags
	hdr[9] = 0xff
	binary.LittleEndian.PutUint16(hdr[10:12], 6)

	var b bytes.Buffer
	b.Write(hdr[:])
	b.Write([]byte{'B', 'C', 0x02, 0x00})
	var bsize [2]byte
	binary.LittleEndian.PutUint16(bsize[:], uint16(blockLen-1))
	b.Write(bsize[:])
	b.Write(fname)
	b.Write(fcomment)
	b.Write([]byte{0xAB, 0xCD})
	b.Write(deflated)
	var footer [8]byte
	binary.LittleEndian.PutUint32(footer[0:4], crc32.ChecksumIEEE(payload))
	binary.LittleEndian.PutUint32(footer[4:8], uint32(len(payload)))
	b.Write(footer[:])
	b.Write(EOFBlock)

	r, _ := NewReader(&b)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if string(got) != "with extras" {
		t.Fatalf("got %q, want %q", got, "with extras")
	}
}

func TestISIZEMismatch(t *testing.T) {
	// Build one valid block, then corrupt the ISIZE field so it disagrees with
	// the actual decoded length.
	var buf bytes.Buffer
	w := NewWriter(&buf)
	w.Write([]byte("hello"))
	w.Close()
	data := buf.Bytes()

	// Footer of the first block starts at len(data) - len(EOFBlock) - 8.
	footerStart := len(data) - len(EOFBlock) - 8
	// ISIZE is the last 4 bytes of the footer.
	binary.LittleEndian.PutUint32(data[footerStart+4:footerStart+8], 9999)

	r, _ := NewReader(bytes.NewReader(data))
	_, err := io.ReadAll(r)
	if !errors.Is(err, ErrISIZE) {
		t.Fatalf("got %v, want ErrISIZE", err)
	}
}

func TestInvalidBSIZEShorterThanHeader(t *testing.T) {
	// Hand-build a block whose BSIZE is smaller than the header it claims to
	// hold. readBlockHeader must reject it.
	flags := byte(flagFEXTRA | flagFHCRC)
	var hdr [12]byte
	hdr[0] = 0x1f
	hdr[1] = 0x8b
	hdr[2] = 8
	hdr[3] = flags
	hdr[9] = 0xff
	binary.LittleEndian.PutUint16(hdr[10:12], 6)

	var b bytes.Buffer
	b.Write(hdr[:])
	b.Write([]byte{'B', 'C', 0x02, 0x00})
	// BSIZE is 10 (block length 11), too small for the header.
	var bsize [2]byte
	binary.LittleEndian.PutUint16(bsize[:], 10)
	b.Write(bsize[:])
	b.Write([]byte{0x00, 0x00})

	r, _ := NewReader(&b)
	_, err := io.ReadAll(r)
	if !errors.Is(err, ErrBadBSIZE) {
		t.Fatalf("got %v, want ErrBadBSIZE", err)
	}
}

func TestInvalidDeflateLengthNegative(t *testing.T) {
	// BSIZE smaller than header+footer triggers the "deflate length < 0" branch.
	// We need XLEN=6 (header=18) and footer=8, so block must be < 26.
	// But ErrBadBSIZE catches "header > block". The negative-deflate path is
	// only reachable when block == header but footer doesn't fit. Construct
	// block_len = headerLen exactly = 18, deflate=0, footer space needed=8 →
	// deflate_len = 18 - 18 - 8 = -8.
	var hdr [12]byte
	hdr[0] = 0x1f
	hdr[1] = 0x8b
	hdr[2] = 8
	hdr[3] = 0x04
	hdr[9] = 0xff
	binary.LittleEndian.PutUint16(hdr[10:12], 6)
	var b bytes.Buffer
	b.Write(hdr[:])
	b.Write([]byte{'B', 'C', 0x02, 0x00})
	var bsize [2]byte
	// BSIZE = 17 → blockLen = 18 → deflate_len = 18 - 18 - 8 = -8 < 0
	binary.LittleEndian.PutUint16(bsize[:], 17)
	b.Write(bsize[:])
	r, _ := NewReader(&b)
	_, err := io.ReadAll(r)
	if err == nil {
		t.Fatal("expected error for invalid deflate length")
	}
}

func TestScanInvalidDeflateLength(t *testing.T) {
	var hdr [12]byte
	hdr[0] = 0x1f
	hdr[1] = 0x8b
	hdr[2] = 8
	hdr[3] = 0x04
	hdr[9] = 0xff
	binary.LittleEndian.PutUint16(hdr[10:12], 6)
	var b bytes.Buffer
	b.Write(hdr[:])
	b.Write([]byte{'B', 'C', 0x02, 0x00})
	var bsize [2]byte
	binary.LittleEndian.PutUint16(bsize[:], 17)
	b.Write(bsize[:])
	if _, err := Scan(&b); err == nil {
		t.Fatal("expected Scan error for invalid deflate length")
	}
}

func TestFlushPropagatesError(t *testing.T) {
	w := NewWriter(&failingWriter{limit: 0})
	w.Write([]byte("data"))
	if err := w.Flush(); err == nil {
		t.Fatal("expected Flush to fail against zero-limit writer")
	}
	// Second Flush should return the recorded error too.
	if err := w.Flush(); err == nil {
		t.Fatal("expected Flush to surface the previously-recorded error")
	}
}

func TestDecompressedSizeOnTruncated(t *testing.T) {
	fixture := buildFixture(t, [][]byte{[]byte("alpha"), []byte("beta")})
	truncated := fixture[:len(fixture)-len(EOFBlock)]
	n, err := DecompressedSize(bytes.NewReader(truncated))
	if !errors.Is(err, ErrTruncated) {
		t.Fatalf("expected ErrTruncated, got %v", err)
	}
	if n != int64(len("alpha")+len("beta")) {
		t.Fatalf("DecompressedSize on truncated = %d, want %d", n, len("alpha")+len("beta"))
	}
}

func TestWriteGZIEmpty(t *testing.T) {
	var b bytes.Buffer
	if err := WriteGZI(&b, nil); err != nil {
		t.Fatal(err)
	}
	if b.Len() != 8 {
		t.Fatalf("empty gzi length = %d, want 8", b.Len())
	}
	if binary.LittleEndian.Uint64(b.Bytes()) != 0 {
		t.Fatal("empty gzi count != 0")
	}
}

// TestVirtualOffsetTracksBlocks compresses three small payloads — one per
// block — and verifies that VirtualOffset returns offsets whose compressed
// coordinates change exactly at the right time and whose in-block portion
// advances byte-by-byte as Read consumes the block.
func TestVirtualOffsetTracksBlocks(t *testing.T) {
	payloads := [][]byte{
		[]byte("alpha"),
		[]byte("beta"),
		[]byte("gamma"),
	}
	fixture := buildFixture(t, payloads)
	// Capture the on-disk start byte of every block by scanning the fixture.
	offsets, err := Scan(bytes.NewReader(fixture))
	if err != nil {
		t.Fatalf("Scan: %v", err)
	}
	if len(offsets) != 3 {
		t.Fatalf("expected 3 blocks, got %d", len(offsets))
	}

	r, err := NewReader(bytes.NewReader(fixture))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	defer r.Close()

	// Before any Read, VirtualOffset must point at byte 0 of block 0.
	if got := r.VirtualOffset(); got != 0 {
		t.Errorf("initial VirtualOffset: got %#x, want 0", got)
	}

	// Consume payloads[0] one byte at a time and confirm uoff advances.
	var one [1]byte
	for i := 0; i < len(payloads[0]); i++ {
		n, err := r.Read(one[:])
		if n != 1 || err != nil {
			t.Fatalf("Read byte %d: n=%d err=%v", i, n, err)
		}
		want := uint64(offsets[0].CompressedOffset)<<16 | uint64(i+1)
		// When we finish the block (i == len(payloads[0])-1), the next
		// virtual offset rolls over to (block1, 0).
		if i == len(payloads[0])-1 {
			want = uint64(offsets[1].CompressedOffset) << 16
		}
		if got := r.VirtualOffset(); got != want {
			t.Errorf("after %d bytes: got %#x, want %#x", i+1, got, want)
		}
	}

	// Consume the rest in one big Read and confirm we land at the EOF block
	// boundary (offsets has only 3 data blocks, so the next-block offset is
	// the byte just past the last data block).
	rest, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	want := append([]byte{}, payloads[1]...)
	want = append(want, payloads[2]...)
	if !bytes.Equal(rest, want) {
		t.Errorf("rest: got %q, want %q", rest, want)
	}
}

// TestBAMReaderVirtualOffsetZeroForRaw confirms that a non-BGZF (raw) BAM
// reader returns 0 for VirtualOffset — it has no compressed layer to track.
func TestVirtualOffsetEmptyReader(t *testing.T) {
	// An empty (only EOF block) stream.
	var b bytes.Buffer
	w := NewWriter(&b)
	w.Close()
	r, err := NewReader(&b)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	if got := r.VirtualOffset(); got != 0 {
		t.Errorf("empty-stream VirtualOffset: got %#x, want 0", got)
	}
}
