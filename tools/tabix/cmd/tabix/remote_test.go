package main

import (
	"bytes"
	"net/http"
	"net/http/httptest"
	"os"
	"testing"
	"time"
)

// serveBytesMux registers name → body handlers with full HTTP Range support
// (via http.ServeContent) and returns a started httptest.Server. The caller
// must Close it.
func serveBytesMux(t *testing.T, files map[string][]byte) *httptest.Server {
	t.Helper()
	mux := http.NewServeMux()
	for name, body := range files {
		body := body
		mux.HandleFunc("/"+name, func(w http.ResponseWriter, r *http.Request) {
			http.ServeContent(w, r, name, time.Unix(0, 0), bytes.NewReader(body))
		})
	}
	return httptest.NewServer(mux)
}

// TestRunRemoteRegionQuery proves the indexed region-query path works against a
// bgzipped+`.tbi`-indexed file served over HTTP: the tool downloads the sibling
// .tbi via hfile, opens the data file as a ranged-GET-backed seekable handle,
// and seeks to the chunk for the requested region — never reading the whole
// object. The result must match the local indexed query exactly.
func TestRunRemoteRegionQuery(t *testing.T) {
	dir := t.TempDir()
	gz := indexSample(t, dir)

	dataBytes, err := os.ReadFile(gz)
	if err != nil {
		t.Fatalf("read data: %v", err)
	}
	tbiBytes, err := os.ReadFile(gz + ".tbi")
	if err != nil {
		t.Fatalf("read tbi: %v", err)
	}

	srv := serveBytesMux(t, map[string][]byte{
		"in.vcf.gz":     dataBytes,
		"in.vcf.gz.tbi": tbiBytes,
	})
	defer srv.Close()
	url := srv.URL + "/in.vcf.gz"

	// Local reference for the same region.
	var local bytes.Buffer
	if code := run([]string{"-p", "vcf", gz, "chr1:140-210"}, nil, &local, &bytes.Buffer{}); code != 0 {
		t.Fatalf("local run exit=%d", code)
	}

	var remote bytes.Buffer
	if code := run([]string{"-p", "vcf", url, "chr1:140-210"}, nil, &remote, &bytes.Buffer{}); code != 0 {
		t.Fatalf("remote run exit=%d", code)
	}

	wantBody := "chr1\t150\t.\tC\tG\t.\t.\t.\nchr1\t200\t.\tG\tA\t.\t.\t.\n"
	if local.String() != wantBody {
		t.Fatalf("local region query: got %q want %q", local.String(), wantBody)
	}
	if remote.String() != local.String() {
		t.Errorf("remote region query mismatch:\n got %q\nwant %q", remote.String(), local.String())
	}

	// A region with no overlapping records returns nothing without error.
	var none bytes.Buffer
	if code := run([]string{"-p", "vcf", url, "chrUnknown:1-1000"}, nil, &none, &bytes.Buffer{}); code != 0 {
		t.Fatalf("remote unknown-chrom exit=%d", code)
	}
	if none.String() != "" {
		t.Errorf("remote unknown chrom: got %q, want empty", none.String())
	}
}

// TestRunRemoteListChroms exercises the no-seek index path (-l) over HTTP: the
// tool downloads only the .tbi and lists the chromosome names from it.
func TestRunRemoteListChroms(t *testing.T) {
	dir := t.TempDir()
	gz := indexSample(t, dir)

	dataBytes, _ := os.ReadFile(gz)
	tbiBytes, _ := os.ReadFile(gz + ".tbi")
	srv := serveBytesMux(t, map[string][]byte{
		"in.vcf.gz":     dataBytes,
		"in.vcf.gz.tbi": tbiBytes,
	})
	defer srv.Close()

	var out bytes.Buffer
	if code := run([]string{"-l", srv.URL + "/in.vcf.gz"}, nil, &out, &bytes.Buffer{}); code != 0 {
		t.Fatalf("remote -l exit=%d", code)
	}
	if got := out.String(); got != "chr1\nchr2\n" {
		t.Errorf("remote list-chroms: got %q, want %q", got, "chr1\nchr2\n")
	}
}

// TestRunRemoteHeaderQuery confirms the streaming header emit (-h with
// --only-header) also routes its data-file open through the remote opener.
func TestRunRemoteHeaderQuery(t *testing.T) {
	dir := t.TempDir()
	gz := indexSample(t, dir)

	dataBytes, _ := os.ReadFile(gz)
	tbiBytes, _ := os.ReadFile(gz + ".tbi")
	srv := serveBytesMux(t, map[string][]byte{
		"in.vcf.gz":     dataBytes,
		"in.vcf.gz.tbi": tbiBytes,
	})
	defer srv.Close()

	var out bytes.Buffer
	if code := run([]string{"-p", "vcf", "-H", srv.URL + "/in.vcf.gz"}, nil, &out, &bytes.Buffer{}); code != 0 {
		t.Fatalf("remote -H exit=%d", code)
	}
	want := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	if out.String() != want {
		t.Errorf("remote header: got %q, want %q", out.String(), want)
	}
}
