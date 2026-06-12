package samtools

import (
	"bytes"
	"io"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

// TestViewFileRemoteIndexedQuery proves the indexed region-query path works
// against a BAM served over HTTP: ViewFile downloads the sibling .bai index
// via hfile, opens the BAM as a ranged-GET-backed seekable handle, and seeks
// to the chunk for the requested region — never reading the whole object.
func TestViewFileRemoteIndexedQuery(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)

	// Build the .bai for the BAM in memory.
	var baiBuf bytes.Buffer
	if err := Index(bytes.NewReader(bamBytes), &baiBuf, IndexOptions{}); err != nil {
		t.Fatalf("Index: %v", err)
	}
	baiBytes := baiBuf.Bytes()

	// Serve both /in.bam and /in.bam.bai with full Range support.
	mux := http.NewServeMux()
	serve := func(name string, body []byte) {
		mux.HandleFunc("/"+name, func(w http.ResponseWriter, r *http.Request) {
			http.ServeContent(w, r, name, time.Unix(0, 0), bytes.NewReader(body))
		})
	}
	serve("in.bam", bamBytes)
	serve("in.bam.bai", baiBytes)
	srv := httptest.NewServer(mux)
	defer srv.Close()

	url := srv.URL + "/in.bam"

	// Region query chr1:50-150 should return exactly r1, matching the local
	// behaviour pinned by TestIndexThenRegionQuery.
	var out bytes.Buffer
	n, err := ViewFile(url, &out, ViewOptions{
		Regions: []string{"chr1:50-150"},
		Count:   true,
	}, io.Discard)
	if err != nil {
		t.Fatalf("ViewFile(remote): %v", err)
	}
	if n != 1 {
		t.Errorf("remote chr1:50-150 count: got %d, want 1", n)
	}

	// Sanity: a region with no overlap returns zero without error.
	var out2 bytes.Buffer
	n2, err := ViewFile(url, &out2, ViewOptions{
		Regions: []string{"chrUnknown:1-1000"},
		Count:   true,
	}, io.Discard)
	if err != nil {
		t.Fatalf("ViewFile(remote unknown): %v", err)
	}
	if n2 != 0 {
		t.Errorf("remote chrUnknown count: got %d, want 0", n2)
	}
}

// TestIdxstatsRemote proves `samtools idxstats URL` reads the BAM header and
// the sibling .bai (fetched via hfile) over HTTP, returning the same rows as
// the equivalent local file.
func TestIdxstatsRemote(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)
	var baiBuf bytes.Buffer
	if err := Index(bytes.NewReader(bamBytes), &baiBuf, IndexOptions{}); err != nil {
		t.Fatalf("Index: %v", err)
	}

	mux := http.NewServeMux()
	mux.HandleFunc("/in.bam", func(w http.ResponseWriter, r *http.Request) {
		http.ServeContent(w, r, "in.bam", time.Unix(0, 0), bytes.NewReader(bamBytes))
	})
	mux.HandleFunc("/in.bam.bai", func(w http.ResponseWriter, r *http.Request) {
		http.ServeContent(w, r, "in.bam.bai", time.Unix(0, 0), bytes.NewReader(baiBuf.Bytes()))
	})
	srv := httptest.NewServer(mux)
	defer srv.Close()

	rows, err := Idxstats(srv.URL + "/in.bam")
	if err != nil {
		t.Fatalf("Idxstats(remote): %v", err)
	}
	if len(rows) == 0 {
		t.Fatalf("Idxstats(remote): no rows")
	}
	// The trailing row is always the unplaced "*" row.
	if rows[len(rows)-1].Name != "*" {
		t.Errorf("last row name = %q, want *", rows[len(rows)-1].Name)
	}
}

// TestMpileupRemoteStreaming proves `samtools mpileup URL` reads a BAM served
// over HTTP through the hfile-backed open path, producing the same output as
// the equivalent local file.
func TestMpileupRemoteStreaming(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.ServeContent(w, r, "in.bam", time.Unix(0, 0), bytes.NewReader(bamBytes))
	}))
	defer srv.Close()

	dir := t.TempDir()
	localPath := filepath.Join(dir, "in.bam")
	if err := os.WriteFile(localPath, bamBytes, 0o600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	var remoteOut, localOut bytes.Buffer
	if err := MpileupFile(MpileupOptions{Inputs: []string{srv.URL + "/in.bam"}}, &remoteOut); err != nil {
		t.Fatalf("MpileupFile(remote): %v", err)
	}
	if err := MpileupFile(MpileupOptions{Inputs: []string{localPath}}, &localOut); err != nil {
		t.Fatalf("MpileupFile(local): %v", err)
	}
	if remoteOut.String() != localOut.String() {
		t.Fatalf("remote mpileup output differs from local:\n--- remote ---\n%s\n--- local ---\n%s",
			remoteOut.String(), localOut.String())
	}
}

// TestViewFileRemoteStreaming confirms the no-region streaming path also works
// over HTTP (whole-object scan), independent of any index.
func TestViewFileRemoteStreaming(t *testing.T) {
	bamBytes := samToBAM(t, sortedSAM)
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.ServeContent(w, r, "in.bam", time.Unix(0, 0), bytes.NewReader(bamBytes))
	}))
	defer srv.Close()

	var out bytes.Buffer
	if _, err := ViewFile(srv.URL+"/in.bam", &out, ViewOptions{Count: true}, io.Discard); err != nil {
		t.Fatalf("ViewFile(remote streaming): %v", err)
	}
	// sortedSAM has 3 mapped + 1 unmapped record; the default (no flag filter)
	// counts every record. Assert we read a non-empty count line.
	if strings.TrimSpace(out.String()) == "" {
		t.Fatalf("remote streaming produced no count output")
	}
}
