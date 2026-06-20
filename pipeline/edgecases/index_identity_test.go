package edgecases

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// TestIndexByteIdentity builds .bai, .csi and .tbi indexes with our tools and
// with upstream on the *same* bgzipped input, and asserts the index payloads
// match. Identity depends on (a) identical BGZF block boundaries in the shared
// input (guaranteed here — both index the same file) and (b) an identical
// bin/interval layout, including htslib's khash bin ordering and
// compress_binning (the formerly-divergent "extra empty bin 4696" came from
// not reproducing compress_binning). A divergence is a real finding: it
// signals a bin-list or interval-layout discrepancy.
//
// The .bai index is stored uncompressed, so it is compared byte-for-byte. The
// .csi and .tbi indexes are BGZF-compressed; the DEFLATE backend differs
// between our writer (klauspost/compress) and upstream htslib (libdeflate/
// zlib), so the *compressed* bytes are not expected to match even when the
// index is identical (the same caveat the repo documents for all BGZF output,
// see docs/PARITY_ROADMAP.md). For those two we therefore compare the
// BGZF-decompressed index payload, which is the bin/interval layout this test
// exists to pin.
func TestIndexByteIdentity(t *testing.T) {
	fix := smallFixtureDir(t)
	if fix == "" {
		t.Skip("pipeline/.fixtures/small not present; run the fixture generator")
	}

	t.Run("bai", func(t *testing.T) {
		our := ourBin(t, "samtools")
		up := upBin(t, "samtools")
		bam := copyInto(t, filepath.Join(fix, "reads.bam"))
		if bam == "" {
			t.Skip("reads.bam fixture missing")
		}
		dir := filepath.Dir(bam)
		ourIdx := filepath.Join(dir, "our.bai")
		upIdx := filepath.Join(dir, "up.bai")
		mustRun(t, our, "index", "-b", "-o", ourIdx, bam)
		mustRun(t, up, "index", "-b", "-o", upIdx, bam)
		assertByteIdentical(t, "BAI", ourIdx, upIdx)
	})

	t.Run("csi", func(t *testing.T) {
		our := ourBin(t, "samtools")
		up := upBin(t, "samtools")
		bam := copyInto(t, filepath.Join(fix, "reads.bam"))
		if bam == "" {
			t.Skip("reads.bam fixture missing")
		}
		dir := filepath.Dir(bam)
		ourIdx := filepath.Join(dir, "our.csi")
		upIdx := filepath.Join(dir, "up.csi")
		mustRun(t, our, "index", "-c", "-o", ourIdx, bam)
		mustRun(t, up, "index", "-c", "-o", upIdx, bam)
		assertPayloadIdentical(t, "CSI", ourIdx, upIdx)
	})

	t.Run("tbi", func(t *testing.T) {
		ourTabix := ourBin(t, "tabix")
		upTabix := upBin(t, "tabix")
		gz := copyInto(t, filepath.Join(fix, "variants.vcf.gz"))
		if gz == "" {
			t.Skip("variants.vcf.gz fixture missing")
		}
		// tabix writes <input>.tbi next to the input; index two copies so the
		// outputs don't collide.
		dir := filepath.Dir(gz)
		ourGz := filepath.Join(dir, "our.vcf.gz")
		upGz := filepath.Join(dir, "up.vcf.gz")
		copyFileTo(t, gz, ourGz)
		copyFileTo(t, gz, upGz)
		mustRun(t, ourTabix, "-p", "vcf", ourGz)
		mustRun(t, upTabix, "-p", "vcf", upGz)
		assertPayloadIdentical(t, "TBI", ourGz+".tbi", upGz+".tbi")
	})
}

// assertPayloadIdentical BGZF/gzip-decompresses both index files and compares
// the decompressed bytes. It is used for the .csi and .tbi indexes, whose
// compressed framing legitimately differs between DEFLATE backends while the
// decoded index payload must be byte-identical.
func assertPayloadIdentical(t testing.TB, kind, ours, up string) {
	t.Helper()
	ob := gunzipFile(t, ours)
	ub := gunzipFile(t, up)
	if !bytes.Equal(ob, ub) {
		t.Errorf("PARITY GAP: %s index payload is not byte-identical to upstream (ours %d bytes, upstream %d bytes; first diff at byte %d).\n"+
			"  Both indexed the same input, so this is a bin/interval-layout difference, not a BGZF-offset one.",
			kind, len(ob), len(ub), firstDiff(ob, ub))
	}
}

// gunzipFile reads a BGZF/gzip file and returns its decompressed bytes.
func gunzipFile(t testing.TB, path string) []byte {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	gr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatalf("gzip %s: %v", path, err)
	}
	defer gr.Close()
	b, err := io.ReadAll(gr)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return b
}

// copyInto copies src into a fresh t.TempDir() and returns the new path, or ""
// if src does not exist.
func copyInto(t testing.TB, src string) string {
	t.Helper()
	if _, err := os.Stat(src); err != nil {
		return ""
	}
	dir := t.TempDir()
	dst := filepath.Join(dir, filepath.Base(src))
	copyFileTo(t, src, dst)
	return dst
}

func copyFileTo(t testing.TB, src, dst string) {
	t.Helper()
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatalf("read %s: %v", src, err)
	}
	if err := os.WriteFile(dst, b, 0o644); err != nil {
		t.Fatalf("write %s: %v", dst, err)
	}
}

func assertByteIdentical(t testing.TB, kind, ours, up string) {
	t.Helper()
	ob, err := os.ReadFile(ours)
	if err != nil {
		t.Fatalf("read our %s: %v", kind, err)
	}
	ub, err := os.ReadFile(up)
	if err != nil {
		t.Fatalf("read upstream %s: %v", kind, err)
	}
	if !bytes.Equal(ob, ub) {
		t.Errorf("PARITY GAP: %s index is not byte-identical to upstream (ours %d bytes, upstream %d bytes; first diff at byte %d).\n"+
			"  Both indexed the same input, so this is a bin/interval-layout difference, not a BGZF-offset one.",
			kind, len(ob), len(ub), firstDiff(ob, ub))
	}
}

func firstDiff(a, b []byte) int {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	for i := 0; i < n; i++ {
		if a[i] != b[i] {
			return i
		}
	}
	if len(a) != len(b) {
		return n
	}
	return -1
}
