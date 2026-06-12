package alnio

import (
	"bytes"
	"io"
	"net/http"
	"net/http/httptest"
	"testing"
	"time"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// serveBytesRange starts an httptest server serving payload with full HTTP
// Range support, matching the ranged-GET access pattern the hfile backend
// uses. It returns the URL of the served object.
func serveBytesRange(t *testing.T, payload []byte) string {
	t.Helper()
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.ServeContent(w, r, "aln", time.Unix(0, 0), bytes.NewReader(payload))
	}))
	t.Cleanup(srv.Close)
	return srv.URL + "/aln.bam"
}

// buildBAM encodes a one-record BAM (BGZF-wrapped) in memory.
func buildBAM(t *testing.T) []byte {
	t.Helper()
	var buf bytes.Buffer
	bw := sam.NewBAMWriter(&buf)
	hdr := &sam.Header{
		Refs: []sam.Reference{{Name: "ref", Length: 100}},
	}
	if err := bw.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	cig, err := sam.ParseCigar("4M")
	if err != nil {
		t.Fatalf("ParseCigar: %v", err)
	}
	rec := &sam.Record{
		QName: "r1", Flag: 0, RName: "ref", Pos: 5, MapQ: 60,
		Cigar: cig,
		Seq:   "ACGT", Qual: []byte{40, 40, 40, 40},
		RNext: "*", PNext: 0, TLen: 0,
	}
	if err := bw.Write(rec); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return buf.Bytes()
}

// TestOpenReaderRemoteBAM proves a BAM object served over HTTP is opened and
// decoded end-to-end through OpenReader's hfile-backed remote path: format
// detection, BGZF decode and BAM record parsing all run over ranged GETs.
func TestOpenReaderRemoteBAM(t *testing.T) {
	url := serveBytesRange(t, buildBAM(t))

	rd, err := OpenReader(url, "")
	if err != nil {
		t.Fatalf("OpenReader(remote BAM): %v", err)
	}
	defer rd.Close()

	if got := len(rd.Header().Refs); got != 1 {
		t.Fatalf("header has %d refs, want 1", got)
	}
	rec, err := rd.Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if rec.QName != "r1" {
		t.Errorf("QName = %q, want r1", rec.QName)
	}
	if _, err := rd.Read(); err != io.EOF {
		t.Errorf("second Read = %v, want io.EOF", err)
	}
}
