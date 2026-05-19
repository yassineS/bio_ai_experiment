package bgzf

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
)

// BlockOffset describes one BGZF block on disk: its byte position in the
// compressed stream and the corresponding byte position in the decompressed
// data. It mirrors one entry in a .gzi index.
type BlockOffset struct {
	// CompressedOffset is the byte offset of the start of the block in the
	// compressed BGZF stream.
	CompressedOffset int64
	// UncompressedOffset is the byte offset of the block's payload in the
	// virtual decompressed stream.
	UncompressedOffset int64
	// CompressedSize is BSIZE+1 — the block's length on disk.
	CompressedSize int
	// UncompressedSize is the block's ISIZE — the payload length.
	UncompressedSize int
}

// Scan walks every block in a BGZF stream, recording its compressed offset,
// compressed size, and uncompressed size. The deflate payload itself is not
// decoded; we read it as raw bytes and trust the ISIZE field in the footer.
// Scan accepts a possibly truncated stream — it returns ErrTruncated if the
// stream ends without the BGZF EOF block, but only after returning all blocks
// it did manage to parse.
//
// The returned offsets follow the .gzi convention: the first entry has
// CompressedOffset = 0 and UncompressedOffset = 0.
func Scan(r io.Reader) ([]BlockOffset, error) {
	var (
		offsets   []BlockOffset
		compOff   int64
		uncompOff int64
		sawEOF    bool
		blockBuf  = make([]byte, 0, 4096)
	)
	for {
		hdr, err := readBlockHeader(r)
		if err == io.EOF {
			break
		}
		if err != nil {
			return offsets, err
		}
		deflatedLen := hdr.compressedSize - hdr.headerLen - 8
		if deflatedLen < 0 {
			return offsets, fmt.Errorf("bgzip: invalid block layout (deflate length %d)", deflatedLen)
		}
		// Read and discard the deflate body and the 8-byte footer, but parse
		// ISIZE so we know the uncompressed size.
		need := deflatedLen + 8
		if cap(blockBuf) < need {
			blockBuf = make([]byte, need)
		} else {
			blockBuf = blockBuf[:need]
		}
		if _, err := io.ReadFull(r, blockBuf); err != nil {
			return offsets, ioErrUnexpected(err)
		}
		isize := binary.LittleEndian.Uint32(blockBuf[need-4 : need])

		// The EOF marker is an empty block (ISIZE = 0, deflate body = 2 bytes).
		if isize == 0 && deflatedLen == 2 {
			sawEOF = true
			break
		}

		offsets = append(offsets, BlockOffset{
			CompressedOffset:   compOff,
			UncompressedOffset: uncompOff,
			CompressedSize:     hdr.compressedSize,
			UncompressedSize:   int(isize),
		})
		compOff += int64(hdr.compressedSize)
		uncompOff += int64(isize)
	}
	if !sawEOF {
		return offsets, ErrTruncated
	}
	return offsets, nil
}

// DecompressedSize sums the uncompressed sizes of every block in r and returns
// the total. It is the implementation of `bgzip -s`.
func DecompressedSize(r io.Reader) (int64, error) {
	offsets, err := Scan(r)
	if err != nil && !errors.Is(err, ErrTruncated) {
		return 0, err
	}
	var total int64
	for _, b := range offsets {
		total += int64(b.UncompressedSize)
	}
	if err != nil {
		return total, err
	}
	return total, nil
}

// UncompressedOffsetAt returns the uncompressed byte offset that corresponds
// to compressed byte offset compOff. compOff must point to the start of a
// block. It is the implementation of `bgzip -b N`.
//
// Per upstream bgzip's semantics, when compOff lands exactly at the start of
// the i-th block, the returned uncompressed offset is the cumulative size of
// blocks [0, i). If compOff is beyond the last block in the stream, the
// total decompressed size is returned.
func UncompressedOffsetAt(r io.Reader, compOff int64) (int64, error) {
	if compOff < 0 {
		return 0, errors.New("bgzip: negative offset")
	}
	offsets, err := Scan(r)
	if err != nil && !errors.Is(err, ErrTruncated) {
		return 0, err
	}
	if compOff == 0 {
		return 0, nil
	}
	for _, b := range offsets {
		if b.CompressedOffset == compOff {
			return b.UncompressedOffset, nil
		}
		if b.CompressedOffset > compOff {
			return 0, fmt.Errorf("bgzip: offset %d does not start a block", compOff)
		}
	}
	// Past the last block — return the total decompressed size, which is the
	// uncompressed offset just after the final byte.
	if n := len(offsets); n > 0 {
		last := offsets[n-1]
		return last.UncompressedOffset + int64(last.UncompressedSize), nil
	}
	return 0, nil
}

// gziMagicSize is the size in bytes of the count prefix in a .gzi file. The
// .gzi format is not, strictly, magic-prefixed — htslib just writes a uint64
// entry count followed by 2*uint64 per entry. We keep the constant here as
// documentation.
const gziMagicSize = 8

// WriteGZI writes the .gzi index format used by htslib/tabix to w. The format
// is: a little-endian uint64 N giving the number of entries, then for each
// entry two little-endian uint64s — the compressed offset and the
// uncompressed offset of a block.
//
// Per htslib's bgzf_index_dump, the very first block (offset 0/0) is NOT
// written to the .gzi file; only blocks 1..N are. This function applies the
// same convention.
func WriteGZI(w io.Writer, offsets []BlockOffset) error {
	// Drop the leading zero-offset block, if present.
	if len(offsets) > 0 && offsets[0].CompressedOffset == 0 && offsets[0].UncompressedOffset == 0 {
		offsets = offsets[1:]
	}
	var hdr [8]byte
	binary.LittleEndian.PutUint64(hdr[:], uint64(len(offsets)))
	if _, err := w.Write(hdr[:]); err != nil {
		return err
	}
	var entry [16]byte
	for _, off := range offsets {
		binary.LittleEndian.PutUint64(entry[0:8], uint64(off.CompressedOffset))
		binary.LittleEndian.PutUint64(entry[8:16], uint64(off.UncompressedOffset))
		if _, err := w.Write(entry[:]); err != nil {
			return err
		}
	}
	return nil
}

// ReadGZI reads a .gzi file produced by htslib/WriteGZI and returns the
// recorded block offsets, NOT including the implicit leading (0, 0) entry.
func ReadGZI(r io.Reader) ([]BlockOffset, error) {
	var hdr [8]byte
	if _, err := io.ReadFull(r, hdr[:]); err != nil {
		return nil, err
	}
	n := binary.LittleEndian.Uint64(hdr[:])
	out := make([]BlockOffset, 0, n)
	var entry [16]byte
	for i := uint64(0); i < n; i++ {
		if _, err := io.ReadFull(r, entry[:]); err != nil {
			return nil, ioErrUnexpected(err)
		}
		out = append(out, BlockOffset{
			CompressedOffset:   int64(binary.LittleEndian.Uint64(entry[0:8])),
			UncompressedOffset: int64(binary.LittleEndian.Uint64(entry[8:16])),
		})
	}
	return out, nil
}
