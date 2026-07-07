package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// TestRoh_ThreadsByteIdentical verifies --threads only parallelises the BGZF
// I/O: the roh output is byte-for-byte identical for any worker count. It
// covers the input-decode path (RohFile over a BGZF-framed VCF) and the -O z
// output-deflate path (Roh with a 'z' output type).
func TestRoh_ThreadsByteIdentical(t *testing.T) {
	// --- input decode: RohFile over a BGZF .vcf.gz ---
	dir := t.TempDir()
	gzPath := filepath.Join(dir, "in.vcf.gz")
	f, err := os.Create(gzPath)
	if err != nil {
		t.Fatalf("create: %v", err)
	}
	bw := bgzf.NewWriter(f)
	if _, err := bw.Write([]byte(fixtureRoh)); err != nil {
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
		if _, err := RohFile(gzPath, &out, RohOptions{GTsOnly: 30, Threads: threads}); err != nil {
			t.Fatalf("RohFile -@%d: %v", threads, err)
		}
		return out.Bytes()
	}
	if one, many := runFile(1), runFile(4); !bytes.Equal(one, many) {
		t.Errorf("RohFile -@1 (%d bytes) vs -@4 (%d bytes) differ", len(one), len(many))
	}

	// --- output deflate: -O z framed output ---
	runZ := func(threads int) []byte {
		var out bytes.Buffer
		if _, err := Roh(strings.NewReader(fixtureRoh), &out, RohOptions{GTsOnly: 30, OutputTypes: "srz", Threads: threads}); err != nil {
			t.Fatalf("Roh -Oz -@%d: %v", threads, err)
		}
		return out.Bytes()
	}
	if one, many := runZ(1), runZ(4); !bytes.Equal(one, many) {
		t.Errorf("Roh -O z -@1 (%d bytes) vs -@4 (%d bytes) differ", len(one), len(many))
	}
}
