package hfile

import (
	"bytes"
	"io"
	"net/http"
	"net/http/httptest"
	"sync/atomic"
	"testing"
	"time"
)

// TestOpenSeekableReadAtAndSeek checks the buffered remote seekable handle
// against os.File semantics: ReadAt at arbitrary offsets, sequential Read via
// Seek, and EOF behaviour all match the served payload.
func TestOpenSeekableReadAtAndSeek(t *testing.T) {
	body := make([]byte, 300000)
	for i := range body {
		body[i] = byte(i * 7)
	}
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.ServeContent(w, r, "data", time.Unix(0, 0), bytes.NewReader(body))
	}))
	defer srv.Close()

	h, err := OpenSeekable(srv.URL)
	if err != nil {
		t.Fatalf("OpenSeekable: %v", err)
	}
	defer h.Close()

	// Random ReadAt windows, including one straddling EOF.
	for _, c := range []struct{ off, n int }{{0, 100}, {150000, 4096}, {len(body) - 10, 64}} {
		buf := make([]byte, c.n)
		got, err := h.ReadAt(buf, int64(c.off))
		want := len(body) - c.off
		if want > c.n {
			want = c.n
		}
		if got != want {
			t.Fatalf("ReadAt(%d,%d): n=%d want %d (err=%v)", c.off, c.n, got, want, err)
		}
		if !bytes.Equal(buf[:got], body[c.off:c.off+got]) {
			t.Fatalf("ReadAt(%d,%d): wrong bytes", c.off, c.n)
		}
	}

	// Seek + sequential Read reconstructs the whole tail.
	if _, err := h.Seek(100000, io.SeekStart); err != nil {
		t.Fatalf("Seek: %v", err)
	}
	rest, err := io.ReadAll(h)
	if err != nil {
		t.Fatalf("ReadAll after seek: %v", err)
	}
	if !bytes.Equal(rest, body[100000:]) {
		t.Fatalf("sequential read after seek mismatch: got %d bytes", len(rest))
	}
}

// TestOpenSeekableCoalescesReads proves the read-ahead window turns many small
// sequential reads into a small number of ranged GETs — the property that
// makes remote BGZF/BAM decoding practical.
func TestOpenSeekableCoalescesReads(t *testing.T) {
	body := make([]byte, 2<<20) // 2 MiB, fits in a single 4 MiB window
	for i := range body {
		body[i] = byte(i)
	}
	var gets int64
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.Method == http.MethodGet {
			atomic.AddInt64(&gets, 1)
		}
		http.ServeContent(w, r, "data", time.Unix(0, 0), bytes.NewReader(body))
	}))
	defer srv.Close()

	h, err := OpenSeekable(srv.URL)
	if err != nil {
		t.Fatalf("OpenSeekable: %v", err)
	}
	defer h.Close()

	// Read the whole body in tiny 13-byte chunks (mimicking a BGZF header/field
	// read pattern). Without buffering this would be ~160k GETs.
	small := make([]byte, 13)
	total := 0
	for {
		n, err := h.Read(small)
		total += n
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read: %v", err)
		}
	}
	if total != len(body) {
		t.Fatalf("read %d bytes, want %d", total, len(body))
	}
	if g := atomic.LoadInt64(&gets); g > 4 {
		t.Fatalf("buffered reads issued %d GETs for a 2 MiB body; want <= 4", g)
	}
}
