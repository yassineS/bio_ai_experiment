package bcftools

import (
	"bytes"
	"io"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"testing"
	"time"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

// serveFiles registers name → body handlers with full HTTP Range support (via
// http.ServeContent) and returns a started httptest.Server; the caller must
// Close it.
func serveFiles(t *testing.T, files map[string][]byte) *httptest.Server {
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

// buildVCFGzFixture writes sampleVCF as a bgzipped file plus its sibling .tbi
// and returns the data-file path.
func buildVCFGzFixture(t *testing.T) string {
	t.Helper()
	dir := t.TempDir()
	bgzPath := filepath.Join(dir, "x.vcf.gz")
	f, err := os.Create(bgzPath)
	if err != nil {
		t.Fatal(err)
	}
	bw := bgzip.NewWriter(f)
	if _, err := bw.Write([]byte(sampleVCF)); err != nil {
		t.Fatal(err)
	}
	if err := bw.Close(); err != nil {
		t.Fatal(err)
	}
	f.Close()
	cfg, _ := tabix.PresetConfig(tabix.PresetVCF)
	idx, err := tabix.Build(bgzPath, cfg)
	if err != nil {
		t.Fatalf("tabix.Build: %v", err)
	}
	if err := idx.WriteFile(bgzPath + ".tbi"); err != nil {
		t.Fatalf("tabix.WriteFile: %v", err)
	}
	return bgzPath
}

// TestViewFileRemoteVCFGzRegion proves the .tbi-backed region-query path works
// against a bgzipped VCF served over HTTP: ViewFile downloads the sibling .tbi
// via hfile, opens the data file as a ranged-GET-backed seekable handle, and
// seeks to the chunk for the requested region. The result must match the local
// indexed query exactly.
func TestViewFileRemoteVCFGzRegion(t *testing.T) {
	bgzPath := buildVCFGzFixture(t)
	dataBytes, _ := os.ReadFile(bgzPath)
	tbiBytes, _ := os.ReadFile(bgzPath + ".tbi")

	srv := serveFiles(t, map[string][]byte{
		"x.vcf.gz":     dataBytes,
		"x.vcf.gz.tbi": tbiBytes,
	})
	defer srv.Close()
	url := srv.URL + "/x.vcf.gz"

	opts := ViewOptions{Regions: []string{"chr1:90-150"}}

	var local bytes.Buffer
	if _, err := ViewFile(bgzPath, &local, opts, io.Discard); err != nil {
		t.Fatalf("ViewFile(local): %v", err)
	}
	var remote bytes.Buffer
	if _, err := ViewFile(url, &remote, opts, io.Discard); err != nil {
		t.Fatalf("ViewFile(remote): %v", err)
	}

	if remote.String() != local.String() {
		t.Errorf("remote region query mismatch:\nremote:\n%s\nlocal:\n%s", remote.String(), local.String())
	}
	if recordsOf(remote.String()) != 1 {
		t.Errorf("remote region records: got %d, want 1", recordsOf(remote.String()))
	}
}

// TestViewFileRemoteBCFCSIRegion is the BCF + .csi counterpart: a BGZF-wrapped
// BCF served over HTTP is region-queried through its sibling .csi.
func TestViewFileRemoteBCFCSIRegion(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	DP=10
chr1	500	b	C	G	30	PASS	DP=20
chr2	200	c	G	A	30	PASS	DP=30
`
	bcfPath := writeBCFForIndex(t, vcfText)
	if _, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI, Force: true}); err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	dataBytes, _ := os.ReadFile(bcfPath)
	csiBytes, _ := os.ReadFile(bcfPath + ".csi")

	srv := serveFiles(t, map[string][]byte{
		"x.bcf":     dataBytes,
		"x.bcf.csi": csiBytes,
	})
	defer srv.Close()
	url := srv.URL + "/x.bcf"

	opts := ViewOptions{Regions: []string{"chr1:50-200"}}

	var local bytes.Buffer
	if _, err := ViewFile(bcfPath, &local, opts, io.Discard); err != nil {
		t.Fatalf("ViewFile(local bcf): %v", err)
	}
	var remote bytes.Buffer
	if _, err := ViewFile(url, &remote, opts, io.Discard); err != nil {
		t.Fatalf("ViewFile(remote bcf): %v", err)
	}

	if remote.String() != local.String() {
		t.Errorf("remote BCF region mismatch:\nremote:\n%s\nlocal:\n%s", remote.String(), local.String())
	}
	if recordsOf(remote.String()) != 1 {
		t.Errorf("remote BCF region records: got %d, want 1", recordsOf(remote.String()))
	}
}

// TestViewFileRemoteStreaming confirms the no-region streaming path works over
// HTTP (whole-object scan via iohelper), independent of any index.
func TestViewFileRemoteStreaming(t *testing.T) {
	bgzPath := buildVCFGzFixture(t)
	dataBytes, _ := os.ReadFile(bgzPath)

	// Serve only the data file (no index): ViewFile must fall back to a
	// streaming scan over the remote object.
	srv := serveFiles(t, map[string][]byte{"x.vcf.gz": dataBytes})
	defer srv.Close()

	var remote bytes.Buffer
	if _, err := ViewFile(srv.URL+"/x.vcf.gz", &remote, ViewOptions{}, io.Discard); err != nil {
		t.Fatalf("ViewFile(remote streaming): %v", err)
	}
	if recordsOf(remote.String()) != 3 {
		t.Errorf("remote streaming records: got %d, want 3", recordsOf(remote.String()))
	}
}

// TestQueryFileRemoteVCFGzRegion proves `bcftools query -r URL` runs the
// .tbi-backed region query over a bgzipped VCF served via HTTP, matching the
// local indexed query byte-for-byte.
func TestQueryFileRemoteVCFGzRegion(t *testing.T) {
	bgzPath := buildVCFGzFixture(t)
	dataBytes, _ := os.ReadFile(bgzPath)
	tbiBytes, _ := os.ReadFile(bgzPath + ".tbi")

	srv := serveFiles(t, map[string][]byte{
		"x.vcf.gz":     dataBytes,
		"x.vcf.gz.tbi": tbiBytes,
	})
	defer srv.Close()
	url := srv.URL + "/x.vcf.gz"

	opts := QueryOptions{Format: "%CHROM\t%POS\t%REF\t%ALT\n", Regions: []string{"chr1:90-150"}}

	var local bytes.Buffer
	if _, err := QueryFile(bgzPath, &local, opts, io.Discard); err != nil {
		t.Fatalf("QueryFile(local): %v", err)
	}
	var remote bytes.Buffer
	if _, err := QueryFile(url, &remote, opts, io.Discard); err != nil {
		t.Fatalf("QueryFile(remote): %v", err)
	}
	if remote.String() != local.String() {
		t.Errorf("remote query -r mismatch:\nremote:\n%s\nlocal:\n%s", remote.String(), local.String())
	}
	if remote.Len() == 0 {
		t.Error("remote query -r produced no output")
	}
}
