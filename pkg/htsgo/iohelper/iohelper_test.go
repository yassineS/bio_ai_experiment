package iohelper

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/tools/bgzip/pkg/bgzip"
)

func TestOpenReaderUncompressed(t *testing.T) {
	// Create a temporary uncompressed file
	tmpDir := t.TempDir()
	tmpFile := filepath.Join(tmpDir, "test.txt")

	content := []byte("Hello, World!")
	if err := os.WriteFile(tmpFile, content, 0644); err != nil {
		t.Fatalf("Failed to create test file: %v", err)
	}

	// Open and read
	r, err := OpenReader(tmpFile)
	if err != nil {
		t.Fatalf("OpenReader failed: %v", err)
	}
	defer r.Close()

	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll failed: %v", err)
	}

	if !bytes.Equal(got, content) {
		t.Errorf("Content mismatch: got %q, want %q", got, content)
	}
}

func TestOpenReaderGzip(t *testing.T) {
	// Create a temporary gzip file
	tmpDir := t.TempDir()
	tmpFile := filepath.Join(tmpDir, "test.txt.gz")

	content := []byte("Hello, Gzip World!")

	// Write gzipped content
	f, err := os.Create(tmpFile)
	if err != nil {
		t.Fatalf("Failed to create test file: %v", err)
	}
	gzw := gzip.NewWriter(f)
	if _, err := gzw.Write(content); err != nil {
		t.Fatalf("Failed to write gzipped content: %v", err)
	}
	gzw.Close()
	f.Close()

	// Open and read
	r, err := OpenReader(tmpFile)
	if err != nil {
		t.Fatalf("OpenReader failed: %v", err)
	}
	defer r.Close()

	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll failed: %v", err)
	}

	if !bytes.Equal(got, content) {
		t.Errorf("Content mismatch: got %q, want %q", got, content)
	}
}

// TestOpenReaderGzipNoExtension verifies that a plain-gzip file is still
// decoded correctly even when the filename does not end in .gz — the helper
// detects format by magic bytes, not extension.
func TestOpenReaderGzipNoExtension(t *testing.T) {
	tmpDir := t.TempDir()
	tmpFile := filepath.Join(tmpDir, "test.dat")

	content := []byte("payload without .gz extension")

	f, err := os.Create(tmpFile)
	if err != nil {
		t.Fatalf("Failed to create test file: %v", err)
	}
	gzw := gzip.NewWriter(f)
	if _, err := gzw.Write(content); err != nil {
		t.Fatalf("Failed to write gzipped content: %v", err)
	}
	gzw.Close()
	f.Close()

	r, err := OpenReader(tmpFile)
	if err != nil {
		t.Fatalf("OpenReader failed: %v", err)
	}
	defer r.Close()

	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll failed: %v", err)
	}
	if !bytes.Equal(got, content) {
		t.Errorf("Content mismatch: got %q, want %q", got, content)
	}
}

// TestOpenReaderBGZF verifies transparent BGZF detection: a file written with
// bgzip.NewWriter (with the BGZF EOF block) round-trips correctly through
// OpenReader without the caller needing to know the format.
func TestOpenReaderBGZF(t *testing.T) {
	tmpDir := t.TempDir()
	tmpFile := filepath.Join(tmpDir, "test.vcf.gz")

	content := []byte("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n1\t100\t.\tA\tG\n")

	f, err := os.Create(tmpFile)
	if err != nil {
		t.Fatalf("Failed to create test file: %v", err)
	}
	bgw := bgzip.NewWriter(f)
	if _, err := bgw.Write(content); err != nil {
		t.Fatalf("Failed to write bgzipped content: %v", err)
	}
	if err := bgw.Close(); err != nil {
		t.Fatalf("bgzip.Writer.Close failed: %v", err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("file.Close failed: %v", err)
	}

	r, err := OpenReader(tmpFile)
	if err != nil {
		t.Fatalf("OpenReader failed: %v", err)
	}

	// Confirm sniffing routed us through the BGZF decoder rather than the
	// plain gzip path. This protects against a future regression where a
	// reordered switch silently falls back to compress/gzip (which would
	// still produce the right bytes but lose the random-access foothold).
	if _, ok := r.(*bgzfReadCloser); !ok {
		t.Errorf("expected *bgzfReadCloser for BGZF input, got %T", r)
	}

	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll failed: %v", err)
	}
	if err := r.Close(); err != nil {
		t.Fatalf("Close failed: %v", err)
	}

	if !bytes.Equal(got, content) {
		t.Errorf("Content mismatch: got %q, want %q", got, content)
	}
}

// TestOpenReaderBGZFMultiBlock exercises the BGZF path with input large
// enough that bgzip emits more than one block, ensuring the wrapper streams
// across block boundaries.
func TestOpenReaderBGZFMultiBlock(t *testing.T) {
	tmpDir := t.TempDir()
	tmpFile := filepath.Join(tmpDir, "big.bgz")

	// Build payload larger than MaxBlockSize so the writer emits >= 2 blocks.
	chunk := bytes.Repeat([]byte("ACGTN"), 1024)
	var content []byte
	for len(content) < bgzip.MaxBlockSize*2+1234 {
		content = append(content, chunk...)
	}

	f, err := os.Create(tmpFile)
	if err != nil {
		t.Fatalf("create: %v", err)
	}
	bgw := bgzip.NewWriter(f)
	if _, err := bgw.Write(content); err != nil {
		t.Fatalf("bgzip write: %v", err)
	}
	if err := bgw.Close(); err != nil {
		t.Fatalf("bgzip close: %v", err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("file close: %v", err)
	}

	r, err := OpenReader(tmpFile)
	if err != nil {
		t.Fatalf("OpenReader: %v", err)
	}
	defer r.Close()
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, content) {
		t.Errorf("multi-block content mismatch: got %d bytes, want %d", len(got), len(content))
	}
}

func TestOpenReaderStdinSentinel(t *testing.T) {
	// "-" and "" both route to stdin without decompression. We can't easily
	// read from real stdin in a test, but we can confirm OpenReader returns
	// without touching the filesystem and yields a non-nil ReadCloser.
	for _, name := range []string{"-", ""} {
		r, err := OpenReader(name)
		if err != nil {
			t.Fatalf("OpenReader(%q) failed: %v", name, err)
		}
		if r == nil {
			t.Fatalf("OpenReader(%q) returned nil reader", name)
		}
		// Close must be a no-op for stdin.
		if err := r.Close(); err != nil {
			t.Errorf("Close on stdin sentinel returned error: %v", err)
		}
	}
}

func TestOpenWriterUncompressed(t *testing.T) {
	// Create a temporary file
	tmpDir := t.TempDir()
	tmpFile := filepath.Join(tmpDir, "test.txt")

	content := []byte("Test content")

	// Write
	w, err := OpenWriter(tmpFile)
	if err != nil {
		t.Fatalf("OpenWriter failed: %v", err)
	}
	if _, err := w.Write(content); err != nil {
		t.Fatalf("Write failed: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close failed: %v", err)
	}

	// Read back
	got, err := os.ReadFile(tmpFile)
	if err != nil {
		t.Fatalf("ReadFile failed: %v", err)
	}

	if !bytes.Equal(got, content) {
		t.Errorf("Content mismatch: got %q, want %q", got, content)
	}
}

func TestOpenWriterGzip(t *testing.T) {
	// Create a temporary gzip file
	tmpDir := t.TempDir()
	tmpFile := filepath.Join(tmpDir, "test.txt.gz")

	content := []byte("Test gzip content")

	// Write
	w, err := OpenWriter(tmpFile)
	if err != nil {
		t.Fatalf("OpenWriter failed: %v", err)
	}
	if _, err := w.Write(content); err != nil {
		t.Fatalf("Write failed: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close failed: %v", err)
	}

	// Read back and decompress
	f, err := os.Open(tmpFile)
	if err != nil {
		t.Fatalf("Open failed: %v", err)
	}
	defer f.Close()

	gzr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatalf("gzip.NewReader failed: %v", err)
	}
	defer gzr.Close()

	got, err := io.ReadAll(gzr)
	if err != nil {
		t.Fatalf("ReadAll failed: %v", err)
	}

	if !bytes.Equal(got, content) {
		t.Errorf("Content mismatch: got %q, want %q", got, content)
	}
}

// TestOpenWriterStdoutSentinel covers the "-"/"" routes to stdout, where
// Close must not close the underlying file descriptor.
func TestOpenWriterStdoutSentinel(t *testing.T) {
	for _, name := range []string{"-", ""} {
		w, err := OpenWriter(name)
		if err != nil {
			t.Fatalf("OpenWriter(%q) failed: %v", name, err)
		}
		if w == nil {
			t.Fatalf("OpenWriter(%q) returned nil writer", name)
		}
		if err := w.Close(); err != nil {
			t.Errorf("Close on stdout sentinel returned error: %v", err)
		}
	}
}

func TestOpenReaderNonExistent(t *testing.T) {
	_, err := OpenReader("/nonexistent/file.txt")
	if err == nil {
		t.Error("Expected error for non-existent file, got nil")
	}
}

// TestOpenWriterNonExistentDir confirms that errors creating the output file
// propagate to the caller.
func TestOpenWriterNonExistentDir(t *testing.T) {
	if _, err := OpenWriter("/nonexistent/dir/output.txt"); err == nil {
		t.Error("Expected error when output directory does not exist, got nil")
	}
}

// TestOpenReaderEmptyFile confirms that a zero-byte file is handled as
// uncompressed (no sniffable magic) and yields an empty read.
func TestOpenReaderEmptyFile(t *testing.T) {
	tmpDir := t.TempDir()
	tmpFile := filepath.Join(tmpDir, "empty.txt")
	if err := os.WriteFile(tmpFile, nil, 0644); err != nil {
		t.Fatalf("create: %v", err)
	}
	r, err := OpenReader(tmpFile)
	if err != nil {
		t.Fatalf("OpenReader: %v", err)
	}
	defer r.Close()
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(got) != 0 {
		t.Errorf("expected empty read, got %d bytes", len(got))
	}
}

// TestOpenReaderCorruptGzipExtension ensures that a non-gzip payload with a
// .gz filename is read as plain bytes (since detection is now magic-based)
// rather than reported as a gzip decode error.
func TestOpenReaderCorruptGzipExtension(t *testing.T) {
	tmpDir := t.TempDir()
	tmpFile := filepath.Join(tmpDir, "fake.gz")
	content := []byte("not actually gzipped\n")
	if err := os.WriteFile(tmpFile, content, 0644); err != nil {
		t.Fatalf("write: %v", err)
	}
	r, err := OpenReader(tmpFile)
	if err != nil {
		t.Fatalf("OpenReader: %v", err)
	}
	defer r.Close()
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, content) {
		t.Errorf("content mismatch: got %q want %q", got, content)
	}
}

func TestBgzfSniff(t *testing.T) {
	tests := []struct {
		name string
		in   []byte
		want bool
	}{
		{"empty", nil, false},
		{"too short", []byte{0x1f, 0x8b, 0x08, 0x04}, false},
		{
			"plain gzip (no BC)",
			[]byte{0x1f, 0x8b, 0x08, 0x00, 0, 0, 0, 0, 0, 0xff, 0, 0, 0, 0, 0, 0},
			false,
		},
		{
			"bad magic",
			[]byte{0x00, 0x00, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 'B', 'C', 0x02, 0x00},
			false,
		},
		{
			"FEXTRA but xlen<6",
			[]byte{0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x02, 0, 'B', 'C', 0x02, 0x00},
			false,
		},
		{
			"FEXTRA but wrong SI",
			[]byte{0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 'X', 'Y', 0x02, 0x00},
			false,
		},
		{
			"canonical BGZF prefix",
			[]byte{0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 'B', 'C', 0x02, 0x00},
			true,
		},
		{
			"non-deflate CM",
			[]byte{0x1f, 0x8b, 0x09, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 'B', 'C', 0x02, 0x00},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := bgzfSniff(tt.in); got != tt.want {
				t.Errorf("bgzfSniff(%x) = %v, want %v", tt.in, got, tt.want)
			}
		})
	}
}

// TestNoCloseReader and TestNoCloseWriter exercise the stdin/stdout wrappers
// directly — they're trivially correct but coverage tooling only registers
// them when called.
func TestNoCloseReader(t *testing.T) {
	src := bytes.NewBufferString("abc")
	r := &noCloseReader{r: src}
	buf := make([]byte, 8)
	n, err := r.Read(buf)
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if string(buf[:n]) != "abc" {
		t.Errorf("read %q want %q", buf[:n], "abc")
	}
	if err := r.Close(); err != nil {
		t.Errorf("Close: %v", err)
	}
}

func TestNoCloseWriter(t *testing.T) {
	var sink bytes.Buffer
	w := &noCloseWriter{w: &sink}
	if _, err := w.Write([]byte("xyz")); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if sink.String() != "xyz" {
		t.Errorf("wrote %q want %q", sink.String(), "xyz")
	}
	if err := w.Close(); err != nil {
		t.Errorf("Close: %v", err)
	}
}

// TestWrappersIdempotentClose ensures Close on each wrapper does not panic
// when invoked after fields have already been cleared (defensive paths in
// each Close hit when err == nil and one of the inner closers returns nil).
func TestWrappersIdempotentClose(t *testing.T) {
	if err := (&readCloserWrapper{}).Close(); err != nil {
		t.Errorf("readCloserWrapper zero-value Close: %v", err)
	}
	if err := (&bgzfReadCloser{}).Close(); err != nil {
		t.Errorf("bgzfReadCloser zero-value Close: %v", err)
	}
	if err := (&plainReadCloser{}).Close(); err != nil {
		t.Errorf("plainReadCloser zero-value Close: %v", err)
	}
	if err := (&writeCloserWrapper{}).Close(); err != nil {
		t.Errorf("writeCloserWrapper zero-value Close: %v", err)
	}
}

func TestGzipSniff(t *testing.T) {
	if gzipSniff(nil) {
		t.Error("gzipSniff(nil) = true, want false")
	}
	if gzipSniff([]byte{0x1f}) {
		t.Error("gzipSniff(short) = true, want false")
	}
	if !gzipSniff([]byte{0x1f, 0x8b, 0x08}) {
		t.Error("gzipSniff(gzip) = false, want true")
	}
	if gzipSniff([]byte{0x00, 0x8b}) {
		t.Error("gzipSniff(non-gzip) = true, want false")
	}
}
