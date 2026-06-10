package fastp

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

func TestSplitFileName(t *testing.T) {
	cases := []struct {
		base   string
		index  int
		digits int
		want   string
	}{
		{"out.fq", 0, 4, "0001.out.fq"},
		{"out.fq", 9, 4, "0010.out.fq"},
		{"out.fq", 0, 0, "1.out.fq"},
		{"out.fq", 0, 2, "01.out.fq"},
		{filepath.Join("sub", "out.fq"), 1, 4, filepath.Join("sub", "0002.out.fq")},
	}
	for _, c := range cases {
		if got := splitFileName(c.base, c.index, c.digits); got != c.want {
			t.Errorf("splitFileName(%q,%d,%d) = %q, want %q", c.base, c.index, c.digits, got, c.want)
		}
	}
}

func TestResolveSplitConfig(t *testing.T) {
	// --split N divides the exact count.
	opts := DefaultProcessOptions()
	opts.SplitNumber = 4
	cfg := resolveSplitConfig(opts, 5000)
	if cfg.Size != 1250 || cfg.Number != 4 || cfg.Digits != 4 {
		t.Fatalf("split-by-number cfg = %+v", cfg)
	}
	// Fewer reads than files: size clamps to 1.
	if c := resolveSplitConfig(opts, 2); c.Size != 1 {
		t.Fatalf("clamp size = %d, want 1", c.Size)
	}
	// --split_by_lines uses lines/4 and no file cap.
	opts2 := DefaultProcessOptions()
	opts2.SplitByLines = 4000
	cfg2 := resolveSplitConfig(opts2, 9999)
	if cfg2.Size != 1000 || cfg2.Number != 0 {
		t.Fatalf("split-by-lines cfg = %+v", cfg2)
	}
}

// TestSplitWriterPackBoundary checks that rollover happens only at pack
// (splitPackSize) boundaries, so a Size that isn't a multiple of the pack
// size still produces files quantized to the pack size.
func TestSplitWriterPackBoundary(t *testing.T) {
	dir := t.TempDir()
	base := filepath.Join(dir, "out.fq")
	// Size 300 with pack 256: the first rollover happens after 512 reads
	// (the first pack boundary at or beyond 300).
	sw := newSplitWriter(base, SplitConfig{Size: 300, Digits: 4, Number: 0}, fastq.Phred33)
	total := 600
	for i := 0; i < total; i++ {
		rec := &fastq.Record{
			ID:       "r",
			Sequence: []byte("ACGT"),
			Quality:  []byte("IIII"),
		}
		if err := sw.Write(rec); err != nil {
			t.Fatalf("write %d: %v", i, err)
		}
	}
	if err := sw.Close(); err != nil {
		t.Fatalf("close: %v", err)
	}
	f1 := countRecords(t, filepath.Join(dir, "0001.out.fq"))
	f2 := countRecords(t, filepath.Join(dir, "0002.out.fq"))
	if f1 != 512 {
		t.Fatalf("file1 records = %d, want 512 (pack-quantized)", f1)
	}
	if f2 != 88 {
		t.Fatalf("file2 records = %d, want 88", f2)
	}
}

// TestSplitWriterEmpty verifies a zero-record run still creates 0001.<base>.
func TestSplitWriterEmpty(t *testing.T) {
	dir := t.TempDir()
	base := filepath.Join(dir, "out.fq")
	sw := newSplitWriter(base, SplitConfig{Size: 100, Digits: 4, Number: 0}, fastq.Phred33)
	if err := sw.Close(); err != nil {
		t.Fatalf("close: %v", err)
	}
	if _, err := os.Stat(filepath.Join(dir, "0001.out.fq")); err != nil {
		t.Fatalf("expected empty first split file: %v", err)
	}
}

func countRecords(t *testing.T, path string) int {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return bytes.Count(b, []byte("\n")) / 4
}
