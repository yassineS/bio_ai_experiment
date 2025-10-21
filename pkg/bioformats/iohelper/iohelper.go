// Package iohelper provides utilities for transparent handling of compressed and uncompressed files.
package iohelper

import (
	"compress/gzip"
	"io"
	"os"
	"strings"
)

// OpenReader opens a file for reading, automatically detecting and handling gzip compression.
// If the file has a .gz extension, it returns a gzip reader wrapped around the file.
// The caller is responsible for closing the returned ReadCloser.
func OpenReader(filename string) (io.ReadCloser, error) {
	if filename == "-" || filename == "" {
		// Read from stdin - no gzip detection for stdin
		return &noCloseReader{os.Stdin}, nil
	}

	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	// Check if file is gzip-compressed (by extension)
	if strings.HasSuffix(filename, ".gz") {
		gzr, err := gzip.NewReader(file)
		if err != nil {
			file.Close()
			return nil, err
		}
		return &readCloserWrapper{gzr, file}, nil
	}

	return file, nil
}

// OpenWriter opens a file for writing, automatically handling gzip compression.
// If the file has a .gz extension, it returns a gzip writer wrapped around the file.
// The caller is responsible for closing the returned WriteCloser.
func OpenWriter(filename string) (io.WriteCloser, error) {
	if filename == "-" || filename == "" {
		// Write to stdout - no gzip for stdout
		return &noCloseWriter{os.Stdout}, nil
	}

	file, err := os.Create(filename)
	if err != nil {
		return nil, err
	}

	// Check if file should be gzip-compressed (by extension)
	if strings.HasSuffix(filename, ".gz") {
		gzw := gzip.NewWriter(file)
		return &writeCloserWrapper{gzw, file}, nil
	}

	return file, nil
}

// noCloseReader wraps stdin/stdout to prevent closing them
type noCloseReader struct {
	r io.Reader
}

func (r *noCloseReader) Read(p []byte) (n int, err error) {
	return r.r.Read(p)
}

func (r *noCloseReader) Close() error {
	return nil // Don't close stdin
}

type noCloseWriter struct {
	w io.Writer
}

func (w *noCloseWriter) Write(p []byte) (n int, err error) {
	return w.w.Write(p)
}

func (w *noCloseWriter) Close() error {
	return nil // Don't close stdout
}

// readCloserWrapper wraps a gzip.Reader and its underlying file for proper cleanup
type readCloserWrapper struct {
	gzr  *gzip.Reader
	file *os.File
}

func (r *readCloserWrapper) Read(p []byte) (n int, err error) {
	return r.gzr.Read(p)
}

func (r *readCloserWrapper) Close() error {
	var err error
	if r.gzr != nil {
		err = r.gzr.Close()
	}
	if r.file != nil {
		if ferr := r.file.Close(); ferr != nil && err == nil {
			err = ferr
		}
	}
	return err
}

// writeCloserWrapper wraps a gzip.Writer and its underlying file for proper cleanup
type writeCloserWrapper struct {
	gzw  *gzip.Writer
	file *os.File
}

func (w *writeCloserWrapper) Write(p []byte) (n int, err error) {
	return w.gzw.Write(p)
}

func (w *writeCloserWrapper) Close() error {
	var err error
	if w.gzw != nil {
		err = w.gzw.Close()
	}
	if w.file != nil {
		if ferr := w.file.Close(); ferr != nil && err == nil {
			err = ferr
		}
	}
	return err
}
