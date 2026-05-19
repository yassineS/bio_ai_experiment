// Package iohelper provides utilities for transparent handling of compressed and uncompressed files.
//
// On read, OpenReader sniffs the first bytes of the stream and returns a
// decoder appropriate for the data:
//
//   - BGZF (RFC 1952 gzip with an `BC` extra subfield, as produced by htslib's
//     bgzip and used by VCF/BCF, BAM, and tabix-indexed files) is decoded by
//     the project's pure-Go BGZF reader in
//     github.com/yassineS/bio_ai_experiment/tools/bgzip/pkg/bgzip. This makes
//     bgzipped inputs round-trip correctly and lays the groundwork for future
//     random-access (tabix) integration.
//   - Plain gzip is decoded by compress/gzip.
//   - Anything else is read raw.
//
// On write, OpenWriter still chooses by extension: .gz inputs get a
// compress/gzip writer; everything else is written raw. Producing BGZF on
// write is not yet wired through this helper (callers that need it can use
// the bgzip package directly).
package iohelper

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/tools/bgzip/pkg/bgzip"
)

// sniffSize is the number of bytes peeked from the head of a file to detect
// the on-disk format. 16 bytes is enough to cover the fixed gzip header plus
// the start of the FEXTRA field where the `BC` subfield lives.
const sniffSize = 16

// OpenReader opens a file for reading, transparently decompressing gzip and
// BGZF inputs. The format is detected by sniffing the first bytes of the
// stream, not by file extension, so a .vcf.gz that is actually BGZF is read
// through the BGZF decoder and a plain .gz is read through compress/gzip.
//
// A filename of "-" or "" reads from standard input without any decompression.
// The caller is responsible for closing the returned ReadCloser; Close releases
// the decompressor (if any) and the underlying file.
func OpenReader(filename string) (io.ReadCloser, error) {
	if filename == "-" || filename == "" {
		// Read from stdin - no compression detection for stdin.
		return &noCloseReader{os.Stdin}, nil
	}

	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	br := bufio.NewReader(file)
	// Peek may return fewer than sniffSize bytes for very short files; that's
	// fine — bgzfSniff and the plain-gzip check both handle short input.
	head, _ := br.Peek(sniffSize)

	switch {
	case bgzfSniff(head):
		bgr, err := bgzip.NewReader(br)
		if err != nil {
			file.Close()
			return nil, err
		}
		return &bgzfReadCloser{bgr: bgr, file: file}, nil
	case gzipSniff(head):
		gzr, err := gzip.NewReader(br)
		if err != nil {
			file.Close()
			return nil, err
		}
		return &readCloserWrapper{gzr: gzr, file: file}, nil
	default:
		// Not compressed — return the buffered reader so peeked bytes are not
		// lost. The caller's Close will release the underlying file.
		return &plainReadCloser{r: br, file: file}, nil
	}
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

// gzipSniff reports whether b begins with the gzip magic bytes (1f 8b).
// It is true for both plain gzip and BGZF; use bgzfSniff first to distinguish
// the two.
func gzipSniff(b []byte) bool {
	return len(b) >= 2 && b[0] == 0x1f && b[1] == 0x8b
}

// bgzfSniff reports whether b is the start of a BGZF gzip member. A BGZF
// block is a gzip member with CM=deflate (08), FLG=FEXTRA (04), and a `BC`
// subfield (SI='B','C', SLEN=2) as the first record in the extra field. That
// makes the first 16 bytes of every block:
//
//	1f 8b 08 04 <mtime:4> <xfl> <os> <xlen:2> 'B' 'C' 02 00 ...
//
// We require xlen ≥ 6 because the BC subfield is 6 bytes long. xlen larger
// than 6 is legal (other subfields may follow BC), but htslib-produced BGZF
// always emits xlen=6 with BC as the only subfield, and we are strict to avoid
// false positives on hand-crafted gzip members.
func bgzfSniff(b []byte) bool {
	if len(b) < 16 {
		return false
	}
	if b[0] != 0x1f || b[1] != 0x8b {
		return false
	}
	if b[2] != 0x08 { // CM = deflate
		return false
	}
	if b[3]&0x04 == 0 { // FLG must have FEXTRA set
		return false
	}
	xlen := uint16(b[10]) | uint16(b[11])<<8
	if xlen < 6 {
		return false
	}
	// First subfield must be SI='B','C', SLEN=2.
	return b[12] == 'B' && b[13] == 'C' && b[14] == 0x02 && b[15] == 0x00
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

// bgzfReadCloser wraps a bgzip.Reader and its underlying file. The bgzip
// reader does not own its source, so Close releases both the decoder and the
// file descriptor.
type bgzfReadCloser struct {
	bgr  *bgzip.Reader
	file *os.File
}

func (r *bgzfReadCloser) Read(p []byte) (int, error) {
	return r.bgr.Read(p)
}

func (r *bgzfReadCloser) Close() error {
	var err error
	if r.bgr != nil {
		err = r.bgr.Close()
	}
	if r.file != nil {
		if ferr := r.file.Close(); ferr != nil && err == nil {
			err = ferr
		}
	}
	return err
}

// plainReadCloser wraps a buffered reader and its underlying file so that
// uncompressed inputs still close the file descriptor when the caller is done.
// We need a separate type (rather than returning the *os.File directly) because
// peeking into the bufio.Reader during format sniffing means the first bytes
// have already been consumed from the file and the buffered reader owns them.
type plainReadCloser struct {
	r    *bufio.Reader
	file *os.File
}

func (r *plainReadCloser) Read(p []byte) (int, error) {
	return r.r.Read(p)
}

func (r *plainReadCloser) Close() error {
	if r.file != nil {
		return r.file.Close()
	}
	return nil
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
