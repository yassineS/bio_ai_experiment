// Package iohelper provides utilities for transparent handling of compressed and uncompressed files.
//
// On read, OpenReader sniffs the first bytes of the stream and returns a
// decoder appropriate for the data:
//
//   - BGZF (RFC 1952 gzip with an `BC` extra subfield, as produced by htslib's
//     bgzip and used by VCF/BCF, BAM, and tabix-indexed files) is decoded by
//     the project's pure-Go BGZF reader in
//     github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf. This makes
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

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/hfile"
)

// sniffSize is the number of bytes peeked from the head of a file to detect
// the on-disk format. 16 bytes is enough to cover the fixed gzip header plus
// the start of the FEXTRA field where the `BC` subfield lives.
const sniffSize = 16

// Format classifies the on-disk format of an alignment stream as detected
// from its leading bytes. It distinguishes the three SAM-family formats so
// callers can route a stream to the right decoder without relying on the
// file extension.
type Format int

const (
	// FormatUnknown is returned when the leading bytes match none of the
	// recognised alignment formats. A plain-text SAM stream that does not
	// start with an `@` header line also classifies as FormatUnknown; the
	// SAM reader is the natural fallback for it.
	FormatUnknown Format = iota
	// FormatSAM is a text SAM stream, recognised by a leading `@` header
	// line.
	FormatSAM
	// FormatBAM is a BAM stream: either BGZF-wrapped (the on-disk form) or
	// an already-decompressed raw `BAM\1` body.
	FormatBAM
	// FormatCRAM is a CRAM stream, recognised by the four-byte `CRAM`
	// file-definition magic.
	FormatCRAM
)

// String returns the format name in lower case ("sam", "bam", "cram" or
// "unknown").
func (f Format) String() string {
	switch f {
	case FormatSAM:
		return "sam"
	case FormatBAM:
		return "bam"
	case FormatCRAM:
		return "cram"
	default:
		return "unknown"
	}
}

// cramMagic is the four-byte signature at the start of every CRAM file.
var cramMagic = [4]byte{'C', 'R', 'A', 'M'}

// classifyHead classifies leading bytes of a stream as SAM, BAM, CRAM or
// unknown. It looks only at magic bytes; it never consumes input. BGZF is
// reported as FormatBAM because every BGZF-wrapped alignment stream this
// project reads is BAM (a BGZF-wrapped CRAM is not a thing — CRAM has its
// own container framing and is never BGZF-wrapped at the file level).
func classifyHead(head []byte) Format {
	if len(head) >= 4 && head[0] == cramMagic[0] && head[1] == cramMagic[1] &&
		head[2] == cramMagic[2] && head[3] == cramMagic[3] {
		return FormatCRAM
	}
	if bgzfSniff(head) {
		return FormatBAM
	}
	if len(head) >= 4 && head[0] == 'B' && head[1] == 'A' && head[2] == 'M' && head[3] == 0x01 {
		return FormatBAM
	}
	if len(head) >= 1 && head[0] == '@' {
		return FormatSAM
	}
	return FormatUnknown
}

// DetectFormat classifies the alignment format of r from its leading bytes
// and returns the classification together with a reader that still yields
// those bytes — peeking is non-destructive, so the returned reader is a
// faithful replacement for r and must be used in its place.
//
// DetectFormat performs no decompression: a CRAM stream is reported as
// FormatCRAM and a BGZF-wrapped BAM as FormatBAM, leaving the actual
// decoding to the caller. It deliberately does not import the sam or cram
// packages, so it stays a low-level building block free of decoder
// dependencies.
func DetectFormat(r io.Reader) (Format, io.Reader, error) {
	br := bufio.NewReader(r)
	head, err := br.Peek(sniffSize)
	if err != nil && err != io.EOF && err != bufio.ErrBufferFull {
		return FormatUnknown, br, err
	}
	return classifyHead(head), br, nil
}

// OpenReader opens a file for reading, transparently decompressing gzip and
// BGZF inputs. The format is detected by sniffing the first bytes of the
// stream, not by file extension, so a .vcf.gz that is actually BGZF is read
// through the BGZF decoder and a plain .gz is read through compress/gzip.
//
// A filename of "-" or "" reads from standard input without any decompression.
// A remote URL (http://, https://, s3:// or gs://) is opened through the
// hfile package, so a bgzipped or gzipped object served from S3, GCS or an
// HTTP server is decompressed transparently just like a local file.
//
// The caller is responsible for closing the returned ReadCloser; Close releases
// the decompressor (if any) and the underlying file or remote handle.
func OpenReader(filename string) (io.ReadCloser, error) {
	if filename == "-" || filename == "" {
		// Read from stdin - no compression detection for stdin.
		return &noCloseReader{os.Stdin}, nil
	}

	src, err := openSource(filename)
	if err != nil {
		return nil, err
	}

	br := bufio.NewReader(src)
	// Peek may return fewer than sniffSize bytes for very short files; that's
	// fine — bgzfSniff and the plain-gzip check both handle short input.
	head, _ := br.Peek(sniffSize)

	switch {
	case bgzfSniff(head):
		bgr, err := bgzip.NewReader(br)
		if err != nil {
			src.Close()
			return nil, err
		}
		return &bgzfReadCloser{bgr: bgr, src: src}, nil
	case gzipSniff(head):
		gzr, err := gzip.NewReader(br)
		if err != nil {
			src.Close()
			return nil, err
		}
		return &readCloserWrapper{gzr: gzr, src: src}, nil
	default:
		// Not compressed — return the buffered reader so peeked bytes are not
		// lost. The caller's Close will release the underlying source.
		return &plainReadCloser{r: br, src: src}, nil
	}
}

// openSource opens filename and returns a raw (undecoded) ReadCloser. Remote
// URLs are routed through the hfile package; everything else is opened as a
// local file. Stdin ("-"/"") is handled by the caller and never reaches here.
func openSource(filename string) (io.ReadCloser, error) {
	if hfile.IsRemote(filename) {
		h, err := hfile.Open(filename)
		if err != nil {
			return nil, err
		}
		// hfile.Handle is an io.ReadCloser (Read + Close); the extra ReaderAt
		// and Size methods are unused on this sequential path.
		return h, nil
	}
	return os.Open(filename)
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

// readCloserWrapper wraps a gzip.Reader and its underlying source for proper
// cleanup. The source is an *os.File for local inputs or an hfile remote
// handle for URLs.
type readCloserWrapper struct {
	gzr *gzip.Reader
	src io.Closer
}

func (r *readCloserWrapper) Read(p []byte) (n int, err error) {
	return r.gzr.Read(p)
}

func (r *readCloserWrapper) Close() error {
	var err error
	if r.gzr != nil {
		err = r.gzr.Close()
	}
	if r.src != nil {
		if ferr := r.src.Close(); ferr != nil && err == nil {
			err = ferr
		}
	}
	return err
}

// bgzfReadCloser wraps a bgzip.Reader and its underlying source. The bgzip
// reader does not own its source, so Close releases both the decoder and the
// source (a local file descriptor or a remote hfile handle).
type bgzfReadCloser struct {
	bgr *bgzip.Reader
	src io.Closer
}

func (r *bgzfReadCloser) Read(p []byte) (int, error) {
	return r.bgr.Read(p)
}

func (r *bgzfReadCloser) Close() error {
	var err error
	if r.bgr != nil {
		err = r.bgr.Close()
	}
	if r.src != nil {
		if ferr := r.src.Close(); ferr != nil && err == nil {
			err = ferr
		}
	}
	return err
}

// plainReadCloser wraps a buffered reader and its underlying source so that
// uncompressed inputs still close the source when the caller is done. We need
// a separate type (rather than returning the source directly) because peeking
// into the bufio.Reader during format sniffing means the first bytes have
// already been consumed from the source and the buffered reader owns them.
type plainReadCloser struct {
	r   *bufio.Reader
	src io.Closer
}

func (r *plainReadCloser) Read(p []byte) (int, error) {
	return r.r.Read(p)
}

func (r *plainReadCloser) Close() error {
	if r.src != nil {
		return r.src.Close()
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
