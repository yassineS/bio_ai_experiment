package tabix

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"testing"
)

// TestTabixIndexUpstreamOracle builds a .tbi from a real bgzipped VCF and
// asserts the index payload is byte-identical to the one a freshly-built
// upstream `tabix -p vcf` produces. The oracle inputs live under
// testdata/upstream_tbi/:
//
//	in.vcf        - the source VCF (2 contigs, 5 records)
//	in.vcf.gz     - bgzipped by genuine libdeflate-linked bgzip
//	expected.tbi  - genuine upstream `tabix -p vcf` output (BGZF-wrapped)
//
// The decompressed expected.tbi payload is 214 bytes and carries, per
// reference, both the data bin 4681 and the meta/pseudo-bin 37450 (with the
// {voffset span} and {n_mapped, n_unmapped} chunks), plus a trailing
// n_no_coor field — exactly the structures this package emits.
func TestTabixIndexUpstreamOracle(t *testing.T) {
	cfg, err := PresetConfig(PresetVCF)
	if err != nil {
		t.Fatal(err)
	}
	idx, err := Build("testdata/upstream_tbi/in.vcf.gz", cfg)
	if err != nil {
		t.Fatalf("Build: %v", err)
	}
	var got bytes.Buffer
	if err := idx.Write(&got); err != nil {
		t.Fatalf("Write: %v", err)
	}

	want := decompressTBI(t, "testdata/upstream_tbi/expected.tbi")
	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("index payload mismatch\n got (%d bytes) %x\nwant (%d bytes) %x",
			got.Len(), got.Bytes(), len(want), want)
	}
}

// decompressTBI reads a BGZF/gzip-wrapped .tbi file and returns its raw index
// payload. BGZF is gzip-compatible, so compress/gzip decodes it.
func decompressTBI(t *testing.T, path string) []byte {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	gr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatalf("gzip reader: %v", err)
	}
	defer gr.Close()
	gr.Multistream(true)
	out, err := io.ReadAll(gr)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return out
}
