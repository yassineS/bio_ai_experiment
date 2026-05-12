package iohelper

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"testing"
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

func TestOpenReaderNonExistent(t *testing.T) {
	_, err := OpenReader("/nonexistent/file.txt")
	if err == nil {
		t.Error("Expected error for non-existent file, got nil")
	}
}
