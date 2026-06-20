package alnio

import (
	"bytes"
	"io"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
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

// TestOpenReaderRemoteCRAM proves a real CRAM file served over HTTP (with full
// Range support, as S3/GCS provide) decodes through the hfile-backed remote
// path identically to the same file read locally. It uses the vendored
// htslib/samtools CRAM corpus fixture.
func TestOpenReaderRemoteCRAM(t *testing.T) {
	path := filepath.Join(samtoolsTestDir, "dat/test_input_1_a.cram")
	data, err := os.ReadFile(path)
	if err != nil {
		t.Skipf("samtools submodule not checked out; run `git submodule update --init reference_code/samtools` to provide %s: %v", path, err)
	}

	// Local decode for the expected record set.
	lf, err := os.Open(path)
	if err != nil {
		t.Fatalf("open local CRAM: %v", err)
	}
	defer lf.Close()
	localRd, err := NewReader(lf)
	if err != nil {
		t.Fatalf("NewReader(local CRAM): %v", err)
	}
	wantRecs, err := readAll(localRd)
	if err != nil {
		t.Fatalf("local CRAM ReadAll: %v", err)
	}

	// Serve the same bytes over HTTP with Range support and decode via the URL.
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.ServeContent(w, r, "test.cram", time.Unix(0, 0), bytes.NewReader(data))
	}))
	defer srv.Close()

	rd, err := OpenReader(srv.URL+"/test.cram", "")
	if err != nil {
		t.Fatalf("OpenReader(remote CRAM): %v", err)
	}
	defer rd.Close()
	gotRecs, err := readAll(rd)
	if err != nil {
		t.Fatalf("remote CRAM ReadAll: %v", err)
	}

	if len(gotRecs) != len(wantRecs) {
		t.Fatalf("remote CRAM decoded %d records, local decoded %d", len(gotRecs), len(wantRecs))
	}
	if len(gotRecs) == 0 {
		t.Fatal("CRAM fixture decoded to zero records")
	}
	for i := range gotRecs {
		if got, w := recordText(t, gotRecs[i]), recordText(t, wantRecs[i]); got != w {
			t.Errorf("remote CRAM record %d mismatch:\n got: %s\nwant: %s", i, got, w)
		}
	}
	t.Logf("remote CRAM: decoded %d records identically to local", len(gotRecs))
}
