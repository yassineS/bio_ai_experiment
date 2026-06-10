package bgzf

import (
	"bytes"
	"crypto/rand"
	"io"
	"testing"
)

// multiRoundTrip compresses payload with a MultiWriter using the given thread
// count, then decompresses with NewReader and returns the framed bytes and the
// decoded plaintext.
func multiRoundTrip(t *testing.T, payload []byte, threads int) (framed, decoded []byte) {
	t.Helper()
	var buf bytes.Buffer
	w, err := NewMultiWriter(&buf, DefaultCompression, threads)
	if err != nil {
		t.Fatalf("NewMultiWriter: %v", err)
	}
	if _, err := w.Write(payload); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	framed = append([]byte(nil), buf.Bytes()...)

	r, err := NewReader(bytes.NewReader(framed))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	decoded, err = io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if err := r.Close(); err != nil {
		t.Fatalf("Reader.Close: %v", err)
	}
	return framed, decoded
}

func TestMultiWriterRoundTrip(t *testing.T) {
	big := make([]byte, MaxBlockSize*5+123)
	if _, err := rand.Read(big); err != nil {
		t.Fatalf("rand: %v", err)
	}
	repetitive := bytes.Repeat([]byte("ACGTACGTNNNN\n"), 100000)

	cases := []struct {
		name    string
		payload []byte
	}{
		{"empty", nil},
		{"tiny", []byte("hello world")},
		{"exact-one-block", bytes.Repeat([]byte("x"), MaxBlockSize)},
		{"just-over-one-block", bytes.Repeat([]byte("y"), MaxBlockSize+1)},
		{"multi-block-random", big},
		{"multi-block-repetitive", repetitive},
	}

	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			for _, threads := range []int{1, 2, 4, 8} {
				_, decoded := multiRoundTrip(t, tc.payload, threads)
				if !bytes.Equal(decoded, tc.payload) {
					t.Fatalf("threads=%d: decoded %d bytes, want %d (mismatch)", threads, len(decoded), len(tc.payload))
				}
			}
		})
	}
}

// TestMultiWriterStructuralValidity confirms the produced stream parses block
// by block, ends in the canonical EOF marker, and that Scan agrees with the
// expected block count and sizes.
func TestMultiWriterStructuralValidity(t *testing.T) {
	payload := bytes.Repeat([]byte("Z"), MaxBlockSize*3+7)
	framed, _ := multiRoundTrip(t, payload, 4)

	if !bytes.HasSuffix(framed, EOFBlock) {
		t.Fatalf("stream does not end with the canonical EOF block")
	}

	offsets, err := Scan(bytes.NewReader(framed))
	if err != nil {
		t.Fatalf("Scan: %v", err)
	}
	// 3 full blocks + 1 partial (the +7 tail) = 4 blocks.
	if len(offsets) != 4 {
		t.Fatalf("Scan found %d blocks, want 4", len(offsets))
	}
	var total int64
	for i, off := range offsets {
		if i < 3 && off.UncompressedSize != MaxBlockSize {
			t.Errorf("block %d uncompressed size = %d, want %d", i, off.UncompressedSize, MaxBlockSize)
		}
		total += int64(off.UncompressedSize)
	}
	if total != int64(len(payload)) {
		t.Fatalf("sum of block sizes = %d, want %d", total, len(payload))
	}
	// Compressed offsets must be strictly increasing and in order.
	for i := 1; i < len(offsets); i++ {
		if offsets[i].CompressedOffset <= offsets[i-1].CompressedOffset {
			t.Fatalf("block offsets out of order at %d", i)
		}
	}
}

// TestMultiWriterMatchesSingleThreadPlaintext checks that single- and
// multi-threaded output decompress to identical plaintext for the same input.
// The compressed bytes need not be identical, but the recovered data must be.
func TestMultiWriterMatchesSingleThreadPlaintext(t *testing.T) {
	payload := make([]byte, MaxBlockSize*4+555)
	if _, err := rand.Read(payload); err != nil {
		t.Fatalf("rand: %v", err)
	}

	var single bytes.Buffer
	sw := NewWriter(&single)
	if _, err := sw.Write(payload); err != nil {
		t.Fatalf("single Write: %v", err)
	}
	if err := sw.Close(); err != nil {
		t.Fatalf("single Close: %v", err)
	}

	for _, threads := range []int{2, 3, 7} {
		_, decoded := multiRoundTrip(t, payload, threads)
		if !bytes.Equal(decoded, payload) {
			t.Fatalf("threads=%d plaintext mismatch", threads)
		}
	}

	// The single-threaded output (one block boundary per MaxBlockSize) is the
	// canonical block layout. At equal block boundaries the MultiWriter must
	// produce byte-identical output, since each block is compressed the same
	// way; verify that for the same single block-sized boundaries.
	var multi bytes.Buffer
	mw, err := NewMultiWriter(&multi, DefaultCompression, 4)
	if err != nil {
		t.Fatalf("NewMultiWriter: %v", err)
	}
	if _, err := mw.Write(payload); err != nil {
		t.Fatalf("multi Write: %v", err)
	}
	if err := mw.Close(); err != nil {
		t.Fatalf("multi Close: %v", err)
	}
	if !bytes.Equal(single.Bytes(), multi.Bytes()) {
		t.Fatalf("single and multi compressed bytes differ; block boundaries should match")
	}
}

// TestMultiWriterManySmallWrites stresses the buffering path with many small
// writes that straddle block boundaries.
func TestMultiWriterManySmallWrites(t *testing.T) {
	var want bytes.Buffer
	var buf bytes.Buffer
	w, err := NewMultiWriter(&buf, DefaultCompression, 4)
	if err != nil {
		t.Fatalf("NewMultiWriter: %v", err)
	}
	chunk := bytes.Repeat([]byte("abc"), 1000) // 3000 bytes, not a divisor of MaxBlockSize
	for i := 0; i < 200; i++ {
		if _, err := w.Write(chunk); err != nil {
			t.Fatalf("Write %d: %v", i, err)
		}
		want.Write(chunk)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	r, err := NewReader(&buf)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, want.Bytes()) {
		t.Fatalf("roundtrip mismatch: got %d bytes, want %d", len(got), want.Len())
	}
}

func TestMultiWriterFlush(t *testing.T) {
	var buf bytes.Buffer
	w, err := NewMultiWriter(&buf, DefaultCompression, 2)
	if err != nil {
		t.Fatalf("NewMultiWriter: %v", err)
	}
	if _, err := w.Write([]byte("first")); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	if _, err := w.Write([]byte("second")); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	framed := append([]byte(nil), buf.Bytes()...)
	r, _ := NewReader(bytes.NewReader(framed))
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if string(got) != "firstsecond" {
		t.Fatalf("got %q, want %q", got, "firstsecond")
	}
	// Flushing produced two separate blocks before the trailing close.
	offsets, err := Scan(bytes.NewReader(framed))
	if err != nil {
		t.Fatalf("Scan: %v", err)
	}
	if len(offsets) != 2 {
		t.Fatalf("expected 2 blocks after Flush, got %d", len(offsets))
	}
}

func TestMultiWriterCloseIdempotent(t *testing.T) {
	var buf bytes.Buffer
	w, err := NewMultiWriter(&buf, DefaultCompression, 2)
	if err != nil {
		t.Fatalf("NewMultiWriter: %v", err)
	}
	w.Write([]byte("data"))
	if err := w.Close(); err != nil {
		t.Fatalf("first Close: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("second Close: %v", err)
	}
	if _, err := w.Write([]byte("x")); err == nil {
		t.Fatalf("Write after Close should error")
	}
}

func TestMultiWriterThreadsClamped(t *testing.T) {
	var buf bytes.Buffer
	w, err := NewMultiWriter(&buf, DefaultCompression, 0)
	if err != nil {
		t.Fatalf("NewMultiWriter(threads=0): %v", err)
	}
	w.Write([]byte("clamped"))
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	r, _ := NewReader(&buf)
	got, _ := io.ReadAll(r)
	if string(got) != "clamped" {
		t.Fatalf("got %q", got)
	}
}

func TestMultiWriterBadLevel(t *testing.T) {
	if _, err := NewMultiWriter(io.Discard, 42, 4); err == nil {
		t.Fatalf("expected error for invalid level")
	}
}

// errAfterWriter fails its Nth write to exercise the collector error path.
type errAfterWriter struct {
	calls int
	fail  int
}

func (e *errAfterWriter) Write(p []byte) (int, error) {
	e.calls++
	if e.calls >= e.fail {
		return 0, io.ErrClosedPipe
	}
	return len(p), nil
}

func TestMultiWriterDownstreamError(t *testing.T) {
	ew := &errAfterWriter{fail: 2}
	w, err := NewMultiWriter(ew, DefaultCompression, 4)
	if err != nil {
		t.Fatalf("NewMultiWriter: %v", err)
	}
	// Write several blocks so the collector attempts multiple downstream writes.
	payload := bytes.Repeat([]byte("q"), MaxBlockSize*5)
	// Write may or may not observe the error depending on timing; Close must.
	_, _ = w.Write(payload)
	if err := w.Close(); err == nil {
		t.Fatalf("expected Close to report downstream write error")
	}
}

func BenchmarkMultiWriter(b *testing.B) {
	payload := make([]byte, MaxBlockSize*16)
	rand.Read(payload)
	for _, threads := range []int{1, 2, 4, 8} {
		b.Run(threadName(threads), func(b *testing.B) {
			b.SetBytes(int64(len(payload)))
			for i := 0; i < b.N; i++ {
				w, err := NewMultiWriter(io.Discard, DefaultCompression, threads)
				if err != nil {
					b.Fatal(err)
				}
				if _, err := w.Write(payload); err != nil {
					b.Fatal(err)
				}
				if err := w.Close(); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

func threadName(n int) string {
	switch n {
	case 1:
		return "threads=1"
	case 2:
		return "threads=2"
	case 4:
		return "threads=4"
	case 8:
		return "threads=8"
	default:
		return "threads=N"
	}
}
