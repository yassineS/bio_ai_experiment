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
// OpenRandomAccessBGZF chooses between two backends:
//
//   - When both a samtools-style `.fai` (offsets into the uncompressed
//     stream) and a `.gzi` (BGZF block index) sidecar are present, it
//     serves Fetch through a bgzf.SeekReader that inflates only the
//     blocks overlapping each request — true partial decompression,
//     matching `samtools faidx ref.fa.gz` without holding the whole
//     reference in memory.
//   - Otherwise it falls back to the pragmatic path: decompress the
//     whole BGZF stream into a byte slice and reuse BuildIndex on it.
//     For the small reference genomes we care about in tests and most
//     bedtools use cases, the BGZF stream comfortably fits in memory.
//
// Both paths are byte-for-byte compatible with upstream's
// `getfasta -fi ref.fa.gz` output.
package fasta

import (
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"os"
	"sync"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
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

// gziReaderAt is an io.ReaderAt over a BGZF stream that decompresses only the
// blocks overlapping each ReadAt request, using a .gzi block index. It is the
// partial-decompression backend for OpenRandomAccessBGZF, replacing the
// whole-file in-memory decompress when a .gzi sidecar is present. ReadAt
// requests are serialised with a mutex because the underlying bgzf.SeekReader
// keeps a single decoded-block cursor and is not safe for concurrent use; the
// FASTA Fetch path makes one small ReadAt per call, so the contention cost is
// negligible compared with re-decompressing whole references.
type gziReaderAt struct {
	mu  sync.Mutex
	sr  *bgzip.SeekReader
	src io.Closer
}

// ReadAt satisfies io.ReaderAt by seeking the BGZF stream to the requested
// uncompressed offset and inflating exactly enough blocks to fill p. It mirrors
// os.File.ReadAt semantics: a short read at end of stream returns io.EOF.
func (g *gziReaderAt) ReadAt(p []byte, off int64) (int, error) {
	g.mu.Lock()
	defer g.mu.Unlock()
	if err := g.sr.SeekTo(off); err != nil {
		if errors.Is(err, bgzip.ErrSeekPastEnd) {
			return 0, io.EOF
		}
		return 0, err
	}
	n, err := io.ReadFull(g.sr, p)
	// Normalise the partial-read sentinel to io.EOF so callers that key on
	// os.File.ReadAt semantics (fasta.RandomAccess.Fetch tolerates io.EOF on a
	// final short slice) behave identically on the .gzi-backed path.
	if err == io.ErrUnexpectedEOF {
		err = io.EOF
	}
	return n, err
}

// Close releases the underlying file handle.
func (g *gziReaderAt) Close() error {
	if g.src != nil {
		return g.src.Close()
	}
	return nil
}

// newGZIReaderAt builds a gziReaderAt for the BGZF file at path using the
// sidecar <path>.gzi. It returns (nil, nil) when no .gzi is present so the
// caller can fall back to the whole-file in-memory path.
func newGZIReaderAt(path string) (*gziReaderAt, error) {
	gziPath := path + ".gzi"
	gziData, err := os.ReadFile(gziPath)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, nil
		}
		return nil, err
	}
	index, err := bgzip.ReadGZI(bytes.NewReader(gziData))
	if err != nil {
		return nil, fmt.Errorf("fasta: reading .gzi: %w", err)
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	return &gziReaderAt{sr: bgzip.NewSeekReader(f, index), src: f}, nil
}

// OpenPayloadReaderAt opens path and returns an io.ReaderAt over its
// *decompressed* payload, together with a closer for the underlying handle.
// It never holds the whole payload in memory:
//
//   - A plain (non-BGZF) file is returned as its own *os.File (ReadAt reads the
//     bytes directly — the payload is the file).
//   - A BGZF file is served through a partial-decompression SeekReader that
//     inflates only the blocks overlapping each ReadAt. When a sibling .gzi is
//     present it is used as the block index; otherwise the block index is built
//     once by scanning the compressed stream (cheap — Scan never decodes the
//     deflate payloads), so genome-scale FASTQ quality fetches stay bounded.
//
// It is the streaming counterpart of decompressBGZF for callers (e.g. the
// fqidx quality accessor) that need byte-offset random access into the
// uncompressed stream without slurping the whole file. ReadAt follows
// os.File.ReadAt semantics: a short read at end of stream returns io.EOF.
func OpenPayloadReaderAt(path string) (io.ReaderAt, io.Closer, error) {
	isBgzf, err := isBGZFFile(path)
	if err != nil {
		return nil, nil, err
	}
	if !isBgzf {
		f, err := os.Open(path)
		if err != nil {
			return nil, nil, err
		}
		return f, f, nil
	}
	// BGZF: prefer the sidecar .gzi block index (no scan needed); otherwise
	// build the block index by scanning the compressed stream once. Neither
	// path decodes the payload eagerly.
	if ra, gerr := newGZIReaderAt(path); gerr != nil {
		return nil, nil, gerr
	} else if ra != nil {
		return ra, ra, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	offsets, err := bgzip.Scan(f)
	if err != nil {
		f.Close()
		return nil, nil, fmt.Errorf("fasta: scanning BGZF blocks: %w", err)
	}
	sr := bgzip.NewSeekReader(f, offsets)
	return &gziReaderAt{sr: sr, src: f}, f, nil
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

	// Fast path: when BOTH the samtools-style .fai (offsets into the
	// uncompressed stream) AND the .gzi (BGZF block index) are present, we can
	// serve Fetch without ever decompressing the whole file — the .fai gives
	// the uncompressed byte range and the .gzi lets the SeekReader inflate
	// only the overlapping blocks. This matches `samtools faidx ref.fa.gz`
	// behaviour and is the partial-decompression seek the roadmap called for.
	if fi, ferr := LoadIndex(path + ".fai"); ferr == nil {
		ra, gerr := newGZIReaderAt(path)
		if gerr != nil {
			return nil, gerr
		}
		if ra != nil {
			return newRandomAccessWithCloser(ra, fi, ra.Close), nil
		}
	} else if !os.IsNotExist(ferr) {
		return nil, ferr
	}

	// Fallback: decompress the whole BGZF stream into memory and index it.
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
