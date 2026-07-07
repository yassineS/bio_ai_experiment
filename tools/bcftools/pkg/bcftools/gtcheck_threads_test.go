package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// TestGtcheck_ThreadsByteIdentical verifies --threads only parallelises the
// BGZF I/O: the gtcheck report is byte-for-byte identical for any worker count.
// It covers the input-decode path (GtcheckFile over a BGZF-framed VCF) and the
// -O z output-deflate path (Gtcheck with OutputType "z").
func TestGtcheck_ThreadsByteIdentical(t *testing.T) {
	// --- input decode: GtcheckFile over a BGZF .vcf.gz ---
	dir := t.TempDir()
	gzPath := filepath.Join(dir, "q.vcf.gz")
	f, err := os.Create(gzPath)
	if err != nil {
		t.Fatalf("create: %v", err)
	}
	bw := bgzf.NewWriter(f)
	if _, err := bw.Write([]byte(fixtureGtcheck)); err != nil {
		t.Fatalf("bgzf write: %v", err)
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("bgzf close: %v", err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("file close: %v", err)
	}
	runFile := func(threads int) []byte {
		var out bytes.Buffer
		if _, err := GtcheckFile(gzPath, &out, GtcheckOptions{Threads: threads}); err != nil {
			t.Fatalf("GtcheckFile -@%d: %v", threads, err)
		}
		return out.Bytes()
	}
	if one, many := runFile(1), runFile(4); !bytes.Equal(one, many) {
		t.Errorf("GtcheckFile -@1 (%d bytes) vs -@4 (%d bytes) differ", len(one), len(many))
	}

	// --- output deflate: -O z framed output ---
	runZ := func(threads int) []byte {
		var out bytes.Buffer
		if _, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{OutputType: "z", Threads: threads}); err != nil {
			t.Fatalf("Gtcheck -O z -@%d: %v", threads, err)
		}
		return out.Bytes()
	}
	if one, many := runZ(1), runZ(4); !bytes.Equal(one, many) {
		t.Errorf("Gtcheck -O z -@1 (%d bytes) vs -@4 (%d bytes) differ", len(one), len(many))
	}
}
