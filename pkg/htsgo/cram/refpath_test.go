package cram

import (
	"bytes"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

// TestSplitRefPath checks the scheme-aware REF_PATH splitter: ':' separates
// entries except inside a "://" URL scheme.
func TestSplitRefPath(t *testing.T) {
	cases := []struct {
		in   string
		want []string
	}{
		{"https://www.ebi.ac.uk/ena/cram/md5/%s", []string{"https://www.ebi.ac.uk/ena/cram/md5/%s"}},
		{"/local/cache:https://host/%s", []string{"/local/cache", "https://host/%s"}},
		{"https://a/%s:https://b/%s", []string{"https://a/%s", "https://b/%s"}},
		{"", nil},
	}
	for _, c := range cases {
		got := splitRefPath(c.in)
		if strings.Join(got, "|") != strings.Join(c.want, "|") {
			t.Errorf("splitRefPath(%q) = %v, want %v", c.in, got, c.want)
		}
	}
}

// TestRefPathFromEnv verifies REF_PATH parsing: only URL entries become network
// templates, a URL without "%s" gets "/%s" appended, and an unset or
// directory-only REF_PATH yields no network source (opt-in semantics).
func TestRefPathFromEnv(t *testing.T) {
	t.Setenv("REF_PATH", "")
	if _, ok := RefPathFromEnv(); ok {
		t.Error("unset REF_PATH must not enable network fetch")
	}
	t.Setenv("REF_PATH", "/only/a/dir")
	if _, ok := RefPathFromEnv(); ok {
		t.Error("a directory-only REF_PATH must not enable network fetch")
	}
	t.Setenv("REF_PATH", "https://host/ref")
	p, ok := RefPathFromEnv()
	if !ok || p == nil {
		t.Fatal("a URL REF_PATH must enable network fetch")
	}
	if len(p.urls) != 1 || p.urls[0] != "https://host/ref/%s" {
		t.Errorf("a %%s-less URL must get /%%s appended, got %v", p.urls)
	}
}

// TestReferenceBackedDecodeViaRefPath proves a reference-backed CRAM decodes
// when its reference is fetched over the network by MD5 (REF_PATH), with no
// local FASTA and no REF_CACHE. An httptest server stands in for the EBI ENA
// endpoint, serving the contig bases at /<md5> exactly as the registry does.
func TestReferenceBackedDecodeViaRefPath(t *testing.T) {
	cramPath := filepath.Join(samtoolsTestDir, referenceBackedFixture.cram)
	faPath := filepath.Join(samtoolsTestDir, referenceBackedFixture.fasta)
	if _, err := os.Stat(cramPath); err != nil {
		t.Fatalf("samtools submodule not initialised — CRAM fixture unavailable; run `git submodule update --init reference_code/samtools`")
	}
	if _, err := os.Stat(faPath); err != nil {
		t.Fatalf("samtools submodule not initialised — reference FASTA unavailable; run `git submodule update --init reference_code/samtools`")
	}

	probe, err := OpenRecords(cramPath)
	if err != nil {
		t.Fatalf("OpenRecords (probe): %v", err)
	}
	m5 := probe.contigMD5(0)
	probe.Close()
	if m5 == "" {
		t.Fatalf("the CRAM @SQ entry carries no M5 tag — REF_PATH keying not exercised; the %s fixture must carry an M5 tag", referenceBackedFixture.cram)
	}

	seq := readContigSequence(t, faPath, referenceBackedFixture.contig)

	// Serve the contig under its MD5, mirroring the EBI ENA md5 endpoint.
	// http.ServeContent honours Range requests, as the real registry does.
	var hits int
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.URL.Path == "/"+m5 {
			hits++
			http.ServeContent(w, r, m5, time.Unix(0, 0), bytes.NewReader(seq))
			return
		}
		http.NotFound(w, r)
	}))
	defer srv.Close()

	t.Setenv("REF_PATH", srv.URL+"/%s")
	t.Setenv("REF_CACHE", "") // ensure only the network source can resolve

	rr, err := OpenRecords(cramPath)
	if err != nil {
		t.Fatalf("OpenRecords: %v", err)
	}
	defer rr.Close()
	if !rr.UseRefPathFromEnv() {
		t.Fatal("UseRefPathFromEnv did not attach a network source for a URL REF_PATH")
	}
	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("REF_PATH-backed decode: %v", err)
	}
	if rr.NeedsReference() {
		t.Error("a REF_PATH-resolved decode must not report NeedsReference")
	}
	if len(recs) == 0 {
		t.Error("REF_PATH decode produced no records")
	}
	if hits == 0 {
		t.Error("the network reference endpoint was never queried")
	}
	if n := countN(recs); n > len(recs) {
		t.Errorf("REF_PATH decode left %d 'N' bases across %d records — resolution incomplete", n, len(recs))
	}
}
