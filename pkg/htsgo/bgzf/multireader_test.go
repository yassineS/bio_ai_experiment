package bgzf

import (
	"bytes"
	"crypto/rand"
	"fmt"
	"io"
	"testing"
)

// makeBGZF compresses payload into a complete BGZF stream (with EOF block).
func makeBGZF(t *testing.T, payload []byte) []byte {
	t.Helper()
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write(payload); err != nil {
		t.Fatalf("write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("close: %v", err)
	}
	return buf.Bytes()
}

// TestMultiReader_MatchesSequential proves the parallel reader decodes a
// multi-block BGZF stream to exactly the same bytes as the sequential Reader,
// for a range of worker counts.
func TestMultiReader_MatchesSequential(t *testing.T) {
	// A payload several blocks long so the parallel path is exercised.
	payload := make([]byte, 5*MaxBlockSize+1234)
	if _, err := rand.Read(payload); err != nil {
		t.Fatalf("rand: %v", err)
	}
	compressed := makeBGZF(t, payload)

	// Sequential reference decode.
	seq, err := NewReader(bytes.NewReader(compressed))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	want, err := io.ReadAll(seq)
	if err != nil {
		t.Fatalf("sequential ReadAll: %v", err)
	}
	if !bytes.Equal(want, payload) {
		t.Fatalf("sequential decode mismatch: %d vs %d bytes", len(want), len(payload))
	}

	for _, threads := range []int{1, 2, 3, 8} {
		mr, err := NewMultiReader(bytes.NewReader(compressed), threads)
		if err != nil {
			t.Fatalf("NewMultiReader(%d): %v", threads, err)
		}
		got, err := io.ReadAll(mr)
		if err != nil {
			t.Fatalf("threads=%d ReadAll: %v", threads, err)
		}
		_ = mr.Close()
		if !bytes.Equal(got, want) {
			t.Fatalf("threads=%d decode differs from sequential (%d vs %d bytes)", threads, len(got), len(want))
		}
	}
}

// TestMultiReader_Truncated reports ErrTruncated when the EOF block is absent,
// matching the sequential reader's contract.
func TestMultiReader_Truncated(t *testing.T) {
	payload := bytes.Repeat([]byte("ACGT"), 10000)
	compressed := makeBGZF(t, payload)
	// Drop the trailing 28-byte EOF block.
	truncated := compressed[:len(compressed)-len(EOFBlock)]

	mr, err := NewMultiReader(bytes.NewReader(truncated), 4)
	if err != nil {
		t.Fatalf("NewMultiReader: %v", err)
	}
	_, err = io.ReadAll(mr)
	if err != ErrTruncated {
		t.Fatalf("want ErrTruncated, got %v", err)
	}
}

// TestMultiReader_EarlyClose closes a MultiReader after reading only a prefix of
// a multi-block stream, then confirms the goroutines unwind (no leak/hang). Run
// under -race to also assert there is no data race on the teardown path.
func TestMultiReader_EarlyClose(t *testing.T) {
	payload := make([]byte, 8*MaxBlockSize)
	if _, err := rand.Read(payload); err != nil {
		t.Fatalf("rand: %v", err)
	}
	compressed := makeBGZF(t, payload)
	mr, err := NewMultiReader(bytes.NewReader(compressed), 4)
	if err != nil {
		t.Fatalf("NewMultiReader: %v", err)
	}
	// Read only a small prefix, then abandon the stream.
	prefix := make([]byte, 100)
	if _, err := io.ReadFull(mr, prefix); err != nil {
		t.Fatalf("ReadFull prefix: %v", err)
	}
	if err := mr.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	// Close is idempotent.
	if err := mr.Close(); err != nil {
		t.Fatalf("second Close: %v", err)
	}
}

// TestMultiReader_Empty decodes an empty (EOF-only) BGZF stream to zero bytes.
func TestMultiReader_Empty(t *testing.T) {
	compressed := makeBGZF(t, nil)
	mr, err := NewMultiReader(bytes.NewReader(compressed), 4)
	if err != nil {
		t.Fatalf("NewMultiReader: %v", err)
	}
	got, err := io.ReadAll(mr)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(got) != 0 {
		t.Fatalf("want 0 bytes, got %d", len(got))
	}
}

// BenchmarkMultiReader measures pure BGZF block-decode throughput at 1 vs
// several worker goroutines, isolating the inflate work (the part `-@` input
// threading parallelises) from any downstream record parsing. The payload is
// compressible (so deflate does real work per block) and spans many blocks. Run
// with `go test -run=^$ -bench=BenchmarkMultiReader ./pkg/htsgo/bgzf/...`.
func BenchmarkMultiReader(b *testing.B) {
	// ~64 MiB of moderately compressible data => hundreds of BGZF blocks.
	payload := bytes.Repeat([]byte("ACGTACGTNNNNacgtACGTTTTTGGGGCCCCAAAA\n"), 1800000)
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write(payload); err != nil {
		b.Fatalf("write: %v", err)
	}
	if err := w.Close(); err != nil {
		b.Fatalf("close: %v", err)
	}
	compressed := buf.Bytes()
	for _, threads := range []int{1, 2, 4, 8} {
		b.Run(fmt.Sprintf("threads=%d", threads), func(b *testing.B) {
			b.SetBytes(int64(len(payload)))
			for i := 0; i < b.N; i++ {
				mr, err := NewMultiReader(bytes.NewReader(compressed), threads)
				if err != nil {
					b.Fatalf("NewMultiReader: %v", err)
				}
				if _, err := io.Copy(io.Discard, mr); err != nil {
					b.Fatalf("Copy: %v", err)
				}
				mr.Close()
			}
		})
	}
}
