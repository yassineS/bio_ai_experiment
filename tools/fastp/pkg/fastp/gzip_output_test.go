package fastp

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// gunzip decompresses a .gz file with the standard library reader, asserting the
// stream is a valid, standard gzip member (i.e. what the klauspost encoder used
// by the split writer must emit).
func gunzip(t *testing.T, path string) []byte {
	t.Helper()
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	if len(raw) < 2 || raw[0] != 0x1f || raw[1] != 0x8b {
		t.Fatalf("%s is not a gzip stream", path)
	}
	gr, err := gzip.NewReader(bytes.NewReader(raw))
	if err != nil {
		t.Fatalf("gzip.NewReader %s: %v", path, err)
	}
	defer gr.Close()
	out, err := io.ReadAll(gr)
	if err != nil {
		t.Fatalf("gunzip %s: %v", path, err)
	}
	return out
}

// TestSplitWriterGzipOutput drives the split writer to a .gz base path and
// verifies each split file is a standard gzip stream whose decompressed content
// is exactly the FASTQ the writer was fed. This guards the fastp fast-gzip
// output path: the encoder changed (klauspost) but the decompressed records
// must be byte-identical, and the level knob must not corrupt the stream.
func TestSplitWriterGzipOutput(t *testing.T) {
	dir := t.TempDir()
	base := filepath.Join(dir, "out.fq.gz")

	const total = 400
	// Compress the same input at two levels; decompressed output must match
	// regardless of the compression level chosen.
	for _, level := range []int{1, 4, 9} {
		sub := filepath.Join(dir, "l")
		if err := os.MkdirAll(sub, 0o755); err != nil {
			t.Fatal(err)
		}
		base = filepath.Join(sub, "out.fq.gz")
		sw := newSplitWriter(base, SplitConfig{Size: 300, Digits: 4, Threads: 1, ByFileLines: true, CompressLevel: level}, fastq.Phred33)

		var want bytes.Buffer
		w := fastq.NewWriter(&want, fastq.Phred33)
		for i := 0; i < total; i++ {
			sw.SetInputPos(i)
			rec := &fastq.Record{ID: "read", Sequence: []byte("ACGTACGTNN"), Quality: []byte("IIIIIIIIII")}
			if err := sw.Write(rec); err != nil {
				t.Fatalf("level %d write %d: %v", level, i, err)
			}
			if err := w.Write(rec); err != nil {
				t.Fatalf("reference write %d: %v", i, err)
			}
		}
		if err := w.Flush(); err != nil {
			t.Fatal(err)
		}
		if err := sw.Close(); err != nil {
			t.Fatalf("level %d close: %v", level, err)
		}

		// Two split files at Size 300 / pack 256: 512-quantized boundary means
		// file 1 gets 512-capped (but only 400 exist), so all 400 land in file 1.
		var got bytes.Buffer
		for i := 0; ; i++ {
			p := splitFileName(base, i, 4)
			if _, err := os.Stat(p); err != nil {
				break
			}
			got.Write(gunzip(t, p))
		}
		if !bytes.Equal(got.Bytes(), want.Bytes()) {
			t.Fatalf("level %d: decompressed split output differs from input records", level)
		}
		os.RemoveAll(sub)
	}
}
