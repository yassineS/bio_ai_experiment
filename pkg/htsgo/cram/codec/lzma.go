package codec

import (
	"bytes"
	"fmt"
	"io"

	"github.com/ulikunitz/xz"
)

// LZMA — the optional per-block CRAM compression method 3.
//
// CRAM does not use raw "LZMA-alone" framing for method-3 blocks. htslib's
// CRAM implementation (cram/cram_io.c) compresses with liblzma's
// lzma_easy_buffer_encode(level, LZMA_CHECK_CRC32, ...) and decompresses
// with lzma_stream_decoder — both of which operate on the standard .xz
// container stream (the one beginning with the 6-byte "\xFD7zXZ\x00"
// magic). The htscodecs source carries the same call in a commented-out
// experiment (arith_dynamic.c). A method-3 CRAM block payload is therefore
// a complete .xz stream, which ulikunitz/xz reads and writes directly.
//
// This file is the sole place in the repository permitted to import
// github.com/ulikunitz/xz; the dependency is confined here per CLAUDE.md
// and docs/CRAM_ROADMAP.md.

// maxLZMARawSize is a defensive ceiling on the decompressed size accepted
// from a single LZMA block. LZMA can legitimately expand highly redundant
// data by a large factor, so the bound cannot be derived from the
// compressed size — it has to be absolute. 1 GiB sits orders of magnitude
// above any real CRAM block (slices are typically <=10 MB) while still
// rejecting a malicious or corrupt stream that would otherwise drive an
// unbounded allocation. It mirrors maxRANSRawSize used by the rANS codecs.
const maxLZMARawSize = 1 << 30

// LZMADecode decompresses a CRAM method-3 LZMA block. The input is a
// complete .xz container stream, as produced by htslib's CRAM writer; the
// returned slice is its uncompressed contents.
//
// A truncated, corrupt or non-.xz input yields an error rather than a
// panic. Decoding stops with an error once more than maxLZMARawSize bytes
// have been produced, so a hostile stream cannot exhaust memory.
func LZMADecode(in []byte) ([]byte, error) {
	zr, err := xz.NewReader(bytes.NewReader(in))
	if err != nil {
		return nil, fmt.Errorf("lzma: opening xz stream: %w", err)
	}
	// LimitReader caps the work an adversarial stream can force. Read one
	// byte past the ceiling so an over-size stream is detected rather than
	// silently truncated.
	out, err := io.ReadAll(io.LimitReader(zr, maxLZMARawSize+1))
	if err != nil {
		return nil, fmt.Errorf("lzma: decompressing xz stream: %w", err)
	}
	if len(out) > maxLZMARawSize {
		return nil, fmt.Errorf("lzma: decompressed size exceeds the %d-byte safety ceiling", maxLZMARawSize)
	}
	return out, nil
}

// LZMAEncode compresses data into a CRAM method-3 LZMA block: a complete
// .xz container stream matching the framing htslib's CRAM writer emits.
// It is provided so callers (and tests) can produce blocks a CRAM reader
// will accept; decoding is the codec's primary responsibility.
func LZMAEncode(in []byte) ([]byte, error) {
	var buf bytes.Buffer
	zw, err := xz.NewWriter(&buf)
	if err != nil {
		return nil, fmt.Errorf("lzma: opening xz writer: %w", err)
	}
	if _, err := zw.Write(in); err != nil {
		return nil, fmt.Errorf("lzma: writing xz stream: %w", err)
	}
	if err := zw.Close(); err != nil {
		return nil, fmt.Errorf("lzma: closing xz stream: %w", err)
	}
	return buf.Bytes(), nil
}
