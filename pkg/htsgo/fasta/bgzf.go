// BGZF-compressed FASTA support.
//
// Bedtools `getfasta -fi <ref.fa.gz>` (and samtools faidx on a bgzipped
// FASTA) work by combining two side-files:
//
//   - <ref>.fa.gz.fai records the same five columns as a plain `.fai` —
//     name, length, offset, line bases, line bytes — but the "offset"
//     here is into the *uncompressed* virtual stream the BGZF blocks
//     represent.
//   - <ref>.fa.gz.gzi maps every BGZF block's uncompressed-stream offset
//     back to its compressed-stream offset, enabling O(1) seek without
//     full decompression.
//
// For the small reference genomes we care about in tests and most
// bedtools use cases, the BGZF stream comfortably fits in memory.
// OpenRandomAccessBGZF therefore takes the pragmatic path: decompress
// the whole BGZF stream into a byte slice and reuse BuildIndex on it.
// This avoids 200+ lines of partial-decompression seek logic while
// remaining byte-for-byte compatible with upstream's `getfasta -fi
// ref.fa.gz` output. Streaming/.gzi-based seeking can be layered on
// later without API churn — see the `loadGZIBlockOffsets` helper which
// is wired up so callers can verify the sidecar `.gzi` exists.
package fasta

import (
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/tools/bgzip/pkg/bgzip"
)

// bgzfMagic is the first four bytes of every BGZF block. RFC 1952 says a
// gzip member starts with 1f 8b 08 04 when an FEXTRA field is present —
// BGZF wraps that with a BC-tagged extra field. We only need the first
// four bytes to distinguish BGZF from a plain gzip stream.
var bgzfMagic = []byte{0x1f, 0x8b, 0x08, 0x04}

// isBGZFFile reports whether path's first four bytes match the BGZF
// magic (0x1f 0x8b 0x08 0x04). A negative answer is not an error —
// callers fall back to plain-FASTA handling.
func isBGZFFile(path string) (bool, error) {
	f, err := os.Open(path)
	if err != nil {
		return false, err
	}
	defer f.Close()
	var buf [4]byte
	n, err := io.ReadFull(f, buf[:])
	if err != nil && !errors.Is(err, io.EOF) && !errors.Is(err, io.ErrUnexpectedEOF) {
		return false, err
	}
	if n < 4 {
		return false, nil
	}
	return bytes.Equal(buf[:], bgzfMagic), nil
}

// gziEntry mirrors `bgzip.BlockOffset` for the subset of fields we need
// (compressed and uncompressed offsets). Defined here so the fasta
// package does not leak the bgzip type into its public API.
type gziEntry struct {
	CompressedOffset   int64
	UncompressedOffset int64
}

// loadGZIBlockOffsets reads a `.gzi` file produced by `bgzip --reindex`
// or htslib's `bgzf_index_dump`. The on-disk format is:
//
//	uint64 N (little-endian) — number of entries
//	N * (uint64 compressed_offset, uint64 uncompressed_offset)
//
// The implicit leading (0,0) entry is NOT stored on disk and is NOT
// returned here either — callers must treat block 0 as starting at
// compressed offset 0 / uncompressed offset 0. Returning the entries
// (even though OpenRandomAccessBGZF doesn't yet seek with them) keeps
// the function reusable and the wire format covered by tests.
func loadGZIBlockOffsets(path string) ([]gziEntry, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return readGZIBlockOffsets(f)
}

// readGZIBlockOffsets is the io.Reader counterpart of loadGZIBlockOffsets.
func readGZIBlockOffsets(r io.Reader) ([]gziEntry, error) {
	var hdr [8]byte
	if _, err := io.ReadFull(r, hdr[:]); err != nil {
		return nil, fmt.Errorf("fasta: bad .gzi header: %w", err)
	}
	n := binary.LittleEndian.Uint64(hdr[:])
	if n > 1<<30 {
		return nil, fmt.Errorf("fasta: implausible .gzi entry count %d", n)
	}
	out := make([]gziEntry, 0, n)
	var entry [16]byte
	for i := uint64(0); i < n; i++ {
		if _, err := io.ReadFull(r, entry[:]); err != nil {
			return nil, fmt.Errorf("fasta: short .gzi at entry %d: %w", i, err)
		}
		out = append(out, gziEntry{
			CompressedOffset:   int64(binary.LittleEndian.Uint64(entry[0:8])),
			UncompressedOffset: int64(binary.LittleEndian.Uint64(entry[8:16])),
		})
	}
	return out, nil
}

// decompressBGZF reads the BGZF stream at path and returns the fully
// decompressed payload. Empty BGZF files (just the EOF marker) yield a
// nil slice and no error.
func decompressBGZF(path string) ([]byte, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	br, err := bgzip.NewReader(f)
	if err != nil {
		return nil, fmt.Errorf("fasta: bgzf header: %w", err)
	}
	data, err := io.ReadAll(br)
	if err != nil {
		return nil, fmt.Errorf("fasta: bgzf decompress: %w", err)
	}
	return data, nil
}

// OpenRandomAccessBGZF opens a BGZF-compressed FASTA (`.fa.gz`) for
// random-access fetches. It uses the same on-disk conventions as
// `samtools faidx` on a bgzipped FASTA:
//
//   - <path>.fai (samtools-style index over the *uncompressed* virtual
//     stream) is consulted if present; otherwise the index is built on
//     the fly from the decompressed payload.
//   - <path>.gzi (the BGZF block index produced by `bgzip --reindex` /
//     htslib's `bgzf_index_dump`) is consulted if present; it is not
//     strictly required for the current in-memory implementation, but
//     its presence is verified so callers get an early signal if the
//     sidecar was renamed or lost.
//
// The returned RandomAccess is backed by the fully-decompressed payload
// (held in memory) and is safe to use across concurrent Fetch calls
// because Fetch only reads from the underlying bytes.Reader through
// ReadAt. For large multi-GB references where in-memory decompression
// becomes prohibitive, see the package doc for the partial-decompression
// roadmap.
func OpenRandomAccessBGZF(path string) (*RandomAccess, error) {
	bgzf, err := isBGZFFile(path)
	if err != nil {
		return nil, err
	}
	if !bgzf {
		return nil, fmt.Errorf("fasta: %q is not BGZF-compressed", path)
	}
	payload, err := decompressBGZF(path)
	if err != nil {
		return nil, err
	}
	// Best-effort: probe for a sibling .gzi so misnamed sidecars surface
	// early. Absence is not fatal — we already have the full payload.
	if _, err := os.Stat(path + ".gzi"); err == nil {
		if _, err := loadGZIBlockOffsets(path + ".gzi"); err != nil {
			return nil, fmt.Errorf("fasta: reading .gzi: %w", err)
		}
	}

	// Try the on-disk .fai first (samtools-compatible — its offsets are
	// into the uncompressed stream, which exactly matches our payload).
	var idx *Index
	if fi, err := LoadIndex(path + ".fai"); err == nil {
		idx = fi
	} else if !os.IsNotExist(err) {
		return nil, err
	} else {
		idx, err = buildIndexFromBytes(payload, false)
		if err != nil {
			return nil, err
		}
	}
	return NewRandomAccess(bytes.NewReader(payload), idx), nil
}

// OpenRandomAccessBGZFFullHeader is the `-fullHeader` analogue of
// OpenRandomAccessBGZF: contigs are keyed by the full header line
// (whitespace included) instead of the first-token name. As with
// OpenRandomAccessFullHeader, any on-disk `.fai` is intentionally
// ignored and the index is rebuilt from the decompressed payload so
// the in-memory keys match upstream's lookup convention.
func OpenRandomAccessBGZFFullHeader(path string) (*RandomAccess, error) {
	bgzf, err := isBGZFFile(path)
	if err != nil {
		return nil, err
	}
	if !bgzf {
		return nil, fmt.Errorf("fasta: %q is not BGZF-compressed", path)
	}
	payload, err := decompressBGZF(path)
	if err != nil {
		return nil, err
	}
	idx, err := buildIndexFromBytes(payload, true)
	if err != nil {
		return nil, err
	}
	return NewRandomAccess(bytes.NewReader(payload), idx), nil
}
