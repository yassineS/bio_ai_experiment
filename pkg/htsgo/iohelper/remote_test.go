package iohelper

import (
	"bytes"
	"io"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
	"time"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// serveBytes starts an httptest server that serves payload with full HTTP
// Range support (via http.ServeContent), mirroring the ranged-GET access the
// hfile backend performs. The returned URL points at the served object.
func serveBytes(t *testing.T, payload []byte) string {
	t.Helper()
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.ServeContent(w, r, "data", time.Unix(0, 0), bytes.NewReader(payload))
	}))
	t.Cleanup(srv.Close)
	return srv.URL + "/data"
}

// TestOpenReaderRemotePlain verifies that OpenReader transparently reads an
// uncompressed object served over HTTP, using the hfile backend.
func TestOpenReaderRemotePlain(t *testing.T) {
	want := []byte("plain remote payload\nsecond line\n")
	url := serveBytes(t, want)

	rc, err := OpenReader(url)
	if err != nil {
		t.Fatalf("OpenReader: %v", err)
	}
	defer rc.Close()

	got, err := io.ReadAll(rc)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("remote plain mismatch:\n got=%q\nwant=%q", got, want)
	}
}

// TestOpenReaderRemoteBGZF verifies that a BGZF-compressed object served over
// HTTP is transparently decompressed by OpenReader, exactly as a local .gz
// would be.
func TestOpenReaderRemoteBGZF(t *testing.T) {
	plain := []byte(strings.Repeat("ACGTACGTNN\n", 5000))

	var buf bytes.Buffer
	bw := bgzip.NewWriter(&buf)
	if _, err := bw.Write(plain); err != nil {
		t.Fatalf("bgzf write: %v", err)
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("bgzf close: %v", err)
	}
	url := serveBytes(t, buf.Bytes())

	rc, err := OpenReader(url)
	if err != nil {
		t.Fatalf("OpenReader: %v", err)
	}
	defer rc.Close()

	got, err := io.ReadAll(rc)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, plain) {
		t.Fatalf("remote BGZF roundtrip mismatch: got %d bytes, want %d", len(got), len(plain))
	}
}
