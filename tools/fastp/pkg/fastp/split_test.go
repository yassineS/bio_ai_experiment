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
	if cfg.Size != 1250 || cfg.Number != 4 || cfg.Digits != 4 || cfg.ByFileLines {
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
	if cfg2.Size != 1000 || cfg2.Number != 0 || !cfg2.ByFileLines {
		t.Fatalf("split-by-lines cfg = %+v", cfg2)
	}

	// Worker count is clamped to split.number in byFileNumber mode
	// (options.cpp: thread cannot exceed the number of split files).
	opts3 := DefaultProcessOptions()
	opts3.SplitNumber = 4
	opts3.Threads = 8
	if c := resolveSplitConfig(opts3, 5000); c.Threads != 4 {
		t.Fatalf("threads clamp = %d, want 4", c.Threads)
	}
	// byFileLines does not clamp the worker count.
	opts4 := DefaultProcessOptions()
	opts4.SplitByLines = 4000
	opts4.Threads = 8
	if c := resolveSplitConfig(opts4, 5000); c.Threads != 8 {
		t.Fatalf("byFileLines threads = %d, want 8", c.Threads)
	}
}

// writeSplitAll feeds total all-passing records to sw, announcing each input
// position, then closes it. It mirrors how the processing loops drive the
// splitWriter (SetInputPos before each input read).
func writeSplitAll(t *testing.T, sw *splitWriter, total int) {
	t.Helper()
	for i := 0; i < total; i++ {
		sw.SetInputPos(i)
		rec := &fastq.Record{ID: "r", Sequence: []byte("ACGT"), Quality: []byte("IIII")}
		if err := sw.Write(rec); err != nil {
			t.Fatalf("write %d: %v", i, err)
		}
	}
	if err := sw.Close(); err != nil {
		t.Fatalf("close: %v", err)
	}
}

// TestSplitWriterPackBoundary checks that, in single-thread byFileLines mode,
// rollover happens only at pack (splitPackSize) boundaries, so a Size that
// isn't a multiple of the pack size still produces files quantized to it.
func TestSplitWriterPackBoundary(t *testing.T) {
	dir := t.TempDir()
	base := filepath.Join(dir, "out.fq")
	// Size 300 with pack 256: the first rollover happens after 512 reads
	// (the first pack boundary at or beyond 300).
	sw := newSplitWriter(base, SplitConfig{Size: 300, Digits: 4, Threads: 1, ByFileLines: true}, fastq.Phred33)
	writeSplitAll(t, sw, 600)
	f1 := countRecords(t, filepath.Join(dir, "0001.out.fq"))
	f2 := countRecords(t, filepath.Join(dir, "0002.out.fq"))
	if f1 != 512 {
		t.Fatalf("file1 records = %d, want 512 (pack-quantized)", f1)
	}
	if f2 != 88 {
		t.Fatalf("file2 records = %d, want 88", f2)
	}
}

// TestSplitWriterMultiThread checks the upstream pack/thread round-robin
// distribution: with 1000 reads, byFileNumber size 250 (N=4), and 2 threads,
// thread 0 owns packs 0,2 (files 0001 then 0003) and thread 1 owns packs 1,3
// (files 0002 then 0004). The pack size (256) quantizes the boundaries.
func TestSplitWriterMultiThread(t *testing.T) {
	dir := t.TempDir()
	base := filepath.Join(dir, "out.fq")
	sw := newSplitWriter(base, SplitConfig{Size: 250, Digits: 4, Number: 4, Threads: 2}, fastq.Phred33)
	writeSplitAll(t, sw, 1000)
	// Packs: 256,256,256,232 (positions 0..999). thread0 gets packs 0 and 2,
	// thread1 gets packs 1 and 3. Each thread rolls after its first pack
	// (256 >= 250). So: file0001=pack0=256, file0002=pack1=256,
	// file0003=pack2=256, file0004=pack3=232.
	want := []int{256, 256, 256, 232}
	for i, w := range want {
		got := countRecords(t, splitFileName(base, i, 4))
		if got != w {
			t.Fatalf("file %04d records = %d, want %d", i+1, got, w)
		}
	}
}

// TestSplitWriterEmptyFilesByNumber verifies that when the input has fewer
// reads than the requested file count, the writer still emits split.number
// files, with the trailing ones empty (writeEmptyFilesForSplitting).
func TestSplitWriterEmptyFilesByNumber(t *testing.T) {
	dir := t.TempDir()
	base := filepath.Join(dir, "out.fq")
	sw := newSplitWriter(base, SplitConfig{Size: 1, Digits: 4, Number: 4, Threads: 2}, fastq.Phred33)
	writeSplitAll(t, sw, 3)
	if got := countRecords(t, splitFileName(base, 0, 4)); got != 3 {
		t.Fatalf("file 0001 records = %d, want 3", got)
	}
	for i := 1; i < 4; i++ {
		name := splitFileName(base, i, 4)
		if _, err := os.Stat(name); err != nil {
			t.Fatalf("expected empty split file %s: %v", name, err)
		}
		if got := countRecords(t, name); got != 0 {
			t.Fatalf("file %04d records = %d, want 0", i+1, got)
		}
	}
}

// TestSplitWriterEmpty verifies a zero-record run still creates 0001.<base>.
func TestSplitWriterEmpty(t *testing.T) {
	dir := t.TempDir()
	base := filepath.Join(dir, "out.fq")
	sw := newSplitWriter(base, SplitConfig{Size: 100, Digits: 4, Number: 0, Threads: 1, ByFileLines: true}, fastq.Phred33)
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
