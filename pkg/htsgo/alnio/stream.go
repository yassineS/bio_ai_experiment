package alnio

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/hfile"
)

// isStdin reports whether path names standard input. The empty string and
// "-" are both treated as stdin, matching the convention used across the
// samtools subcommands and pkg/htsgo/iohelper.
func isStdin(path string) bool {
	return path == "" || path == "-"
}

// stdinReader returns the process's standard input as an io.Reader.
func stdinReader() io.Reader { return os.Stdin }

// osOpen opens the named file for reading. It is a thin wrapper over
// os.Open kept as a package function so the file-opening seam is easy to
// stub in tests.
func osOpen(path string) (*os.File, error) { return os.Open(path) }

// openAlnSource opens an alignment file by path or URL for sequential reading.
// A remote URL (http(s)://, s3://, gs://) is opened through the hfile package
// so that `samtools`/`bcftools`-style tools can read a BAM/CRAM/SAM object
// straight from cloud storage; any other path is opened as a local file.
func openAlnSource(path string) (io.ReadCloser, error) {
	if hfile.IsRemote(path) {
		// OpenSeekable adds read-ahead buffering so the sequential alignment
		// scan pulls the object in a few large ranged GETs. The returned
		// SeekHandle satisfies io.ReadCloser.
		return hfile.OpenSeekable(path)
	}
	return osOpen(path)
}

// decompressStream returns a reader that transparently decompresses a
// plain-gzip-compressed SAM stream. It sniffs the leading bytes without
// consuming them: a plain gzip member is wrapped in a gzip decoder, and
// anything else (raw SAM text, or a raw or BGZF-wrapped BAM body) is
// returned unchanged.
//
// BGZF-wrapped BAM is deliberately left compressed: sam.NewReader detects
// the BGZF header itself and routes to its BAM reader, which owns the BGZF
// decode. Only plain gzip — which sam.NewReader does not handle — is
// stripped here.
func decompressStream(r io.Reader) (io.Reader, error) {
	br := bufio.NewReader(r)
	head, err := br.Peek(16)
	if err != nil && err != io.EOF && err != bufio.ErrBufferFull {
		return nil, err
	}
	if isPlainGzip(head) {
		return gzip.NewReader(br)
	}
	return br, nil
}

// isPlainGzip reports whether head begins a plain (non-BGZF) gzip member.
// A BGZF member also begins with the gzip magic, so the BGZF FEXTRA `BC`
// subfield is checked to exclude it.
func isPlainGzip(head []byte) bool {
	if len(head) < 2 || head[0] != 0x1f || head[1] != 0x8b {
		return false
	}
	return !looksLikeBGZF(head)
}

// looksLikeBGZF reports whether head is the start of a BGZF gzip member —
// a deflate gzip member carrying a `BC` extra subfield. It mirrors the
// detection in pkg/htsgo/iohelper and pkg/htsgo/sam.
func looksLikeBGZF(head []byte) bool {
	if len(head) < 16 {
		return false
	}
	if head[0] != 0x1f || head[1] != 0x8b || head[2] != 0x08 {
		return false
	}
	if head[3]&0x04 == 0 {
		return false
	}
	xlen := uint16(head[10]) | uint16(head[11])<<8
	if xlen < 6 {
		return false
	}
	return head[12] == 'B' && head[13] == 'C' && head[14] == 0x02 && head[15] == 0x00
}
