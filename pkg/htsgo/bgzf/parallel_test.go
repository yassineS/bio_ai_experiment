package bgzf

import (
	"bytes"
	"crypto/rand"
	"io"
	"testing"
)

// TestParallelWriterByteIdentical compares the on-disk bytes produced by
// the parallel BGZF writer against the serial Writer across a range of
// payload sizes and worker counts. All outputs must be byte-identical to
// the serial path so the BGZF on-disk format remains canonical regardless
// of -@.
func TestParallelWriterByteIdentical(t *testing.T) {
	sizes := []int{
		0,
		1,
		MaxBlockSize - 1,
		MaxBlockSize,
		MaxBlockSize + 1,
		3*MaxBlockSize + 17,
		17*MaxBlockSize + 1,
	}
	workersList := []int{1, 2, 4, 8}
	for _, n := range sizes {
		payload := make([]byte, n)
		// Deterministic but non-trivial: use a simple LCG so each test
		// run sees the same bytes (matters for the byte-identity check).
		for i := range payload {
			payload[i] = byte((i*1103515245 + 12345) >> 8)
		}
		var serial bytes.Buffer
		sw := NewWriter(&serial)
		if _, err := sw.Write(payload); err != nil {
			t.Fatalf("serial write %d: %v", n, err)
		}
		if err := sw.Close(); err != nil {
			t.Fatalf("serial close %d: %v", n, err)
		}
		for _, workers := range workersList {
			var par bytes.Buffer
			pw, err := NewParallelWriterLevel(&par, DefaultCompression, workers)
			if err != nil {
				t.Fatalf("new parallel %d w%d: %v", n, workers, err)
			}
			if _, err := pw.Write(payload); err != nil {
				t.Fatalf("parallel write %d w%d: %v", n, workers, err)
			}
			if err := pw.Close(); err != nil {
				t.Fatalf("parallel close %d w%d: %v", n, workers, err)
			}
			if !bytes.Equal(serial.Bytes(), par.Bytes()) {
				t.Fatalf("parallel(workers=%d) output differs from serial for size=%d (serial=%d bytes, parallel=%d bytes)",
					workers, n, serial.Len(), par.Len())
			}
		}
	}
}

// TestParallelWriterRandomChunks exercises the Write boundary path: a
// large random payload is fed through both writers in many small randomly
// sized Write() calls, and the on-disk bytes must still match exactly.
func TestParallelWriterRandomChunks(t *testing.T) {
	payload := make([]byte, 7*MaxBlockSize+123)
	if _, err := rand.Read(payload); err != nil {
		t.Fatalf("rand: %v", err)
	}
	chunks := func(rng *bytesReader) [][]byte {
		var out [][]byte
		p := payload
		for len(p) > 0 {
			// Deterministic chunk sizes via the rng so both runs slice
			// identically.
			n := 1 + int(rng.next()%uint32(MaxBlockSize/3))
			if n > len(p) {
				n = len(p)
			}
			out = append(out, p[:n])
			p = p[n:]
		}
		return out
	}

	var serial bytes.Buffer
	sw := NewWriter(&serial)
	for _, c := range chunks(newBytesReader(1)) {
		if _, err := sw.Write(c); err != nil {
			t.Fatalf("serial: %v", err)
		}
	}
	if err := sw.Close(); err != nil {
		t.Fatalf("serial close: %v", err)
	}
	var par bytes.Buffer
	pw, err := NewParallelWriterLevel(&par, DefaultCompression, 4)
	if err != nil {
		t.Fatalf("parallel: %v", err)
	}
	for _, c := range chunks(newBytesReader(1)) {
		if _, err := pw.Write(c); err != nil {
			t.Fatalf("parallel: %v", err)
		}
	}
	if err := pw.Close(); err != nil {
		t.Fatalf("parallel close: %v", err)
	}
	if !bytes.Equal(serial.Bytes(), par.Bytes()) {
		t.Fatalf("parallel output differs from serial (serial=%d, parallel=%d)", serial.Len(), par.Len())
	}
	// Round-trip decode for good measure.
	r, err := NewReader(bytes.NewReader(par.Bytes()))
	if err != nil {
		t.Fatalf("reader: %v", err)
	}
	defer r.Close()
	decoded, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("decode: %v", err)
	}
	if !bytes.Equal(decoded, payload) {
		t.Fatalf("round-trip mismatch (decoded=%d, want=%d)", len(decoded), len(payload))
	}
}

// bytesReader is a tiny deterministic LCG used to make chunk slicing
// reproducible across the two runs of the random-chunks test.
type bytesReader struct{ state uint32 }

func newBytesReader(seed uint32) *bytesReader { return &bytesReader{state: seed | 1} }

func (b *bytesReader) next() uint32 {
	b.state = b.state*1664525 + 1013904223
	return b.state
}
