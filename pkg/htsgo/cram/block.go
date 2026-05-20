package cram

import (
	"bytes"
	"compress/bzip2"
	"compress/gzip"
	"encoding/binary"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram/codec"
)

// CompressionMethod identifies how a CRAM block's data is compressed.
type CompressionMethod byte

// The CRAM block compression methods. Methods 0-2, 4 and 5 are handled
// by this package; 3, 6, 7 and 8 are out of scope for the C3 structural
// parser and decompressing them returns an error.
const (
	CompRaw      CompressionMethod = 0 // Uncompressed.
	CompGzip     CompressionMethod = 1 // RFC 1952 gzip (compress/gzip).
	CompBzip2    CompressionMethod = 2 // bzip2 (compress/bzip2).
	CompLZMA     CompressionMethod = 3 // LZMA — out of scope (PR C-LZMA).
	CompRANS4x8  CompressionMethod = 4 // rANS 4x8 (CRAM v3.0 codec).
	CompRANS4x16 CompressionMethod = 5 // rANS 4x16 (CRAM v3.1 codec).
	CompArith    CompressionMethod = 6 // Range/arithmetic coder — out of scope.
	CompFQZComp  CompressionMethod = 7 // fqzcomp quality codec — out of scope.
	CompNameTok  CompressionMethod = 8 // Name tokeniser — out of scope.
)

// String returns the spec name of the compression method.
func (m CompressionMethod) String() string {
	switch m {
	case CompRaw:
		return "raw"
	case CompGzip:
		return "gzip"
	case CompBzip2:
		return "bzip2"
	case CompLZMA:
		return "lzma"
	case CompRANS4x8:
		return "rans4x8"
	case CompRANS4x16:
		return "rans4x16"
	case CompArith:
		return "arith"
	case CompFQZComp:
		return "fqzcomp"
	case CompNameTok:
		return "name-tokeniser"
	default:
		return fmt.Sprintf("unknown(%d)", byte(m))
	}
}

// BlockContentType identifies the role a CRAM block plays in the file
// structure.
type BlockContentType byte

// The CRAM block content types, per the v3 specification.
const (
	ContentFileHeader        BlockContentType = 0 // SAM header block.
	ContentCompressionHeader BlockContentType = 1 // Per-container compression header.
	ContentMappedSlice       BlockContentType = 2 // Slice header block.
	ContentReserved          BlockContentType = 3 // Reserved / unused.
	ContentExternal          BlockContentType = 4 // External data series block.
	ContentCoreData          BlockContentType = 5 // Core (bit-packed) data block.
)

// String returns a human-readable name for the block content type.
func (t BlockContentType) String() string {
	switch t {
	case ContentFileHeader:
		return "file-header"
	case ContentCompressionHeader:
		return "compression-header"
	case ContentMappedSlice:
		return "slice-header"
	case ContentReserved:
		return "reserved"
	case ContentExternal:
		return "external"
	case ContentCoreData:
		return "core-data"
	default:
		return fmt.Sprintf("unknown(%d)", byte(t))
	}
}

// Block is a single parsed CRAM block: its header fields plus the raw
// (still-compressed) payload bytes. Call Decompress to obtain the
// uncompressed data.
type Block struct {
	// Method is the block's compression method.
	Method CompressionMethod
	// ContentType is the structural role of the block.
	ContentType BlockContentType
	// ContentID is the block's content identifier; for external data
	// blocks it keys the data series the block carries.
	ContentID int32
	// CompressedSize is the declared length of the on-disk payload.
	CompressedSize int32
	// UncompressedSize is the declared length after decompression.
	UncompressedSize int32
	// Data is the raw, still-compressed payload (CompressedSize bytes).
	Data []byte
	// CRC is the 4-byte little-endian CRC32 stored after the block in
	// CRAM v3+; it is zero and unused for v2.
	CRC uint32
}

// Decompress returns the block's uncompressed payload. For a raw block
// it returns the data unchanged. For gzip, bzip2, rANS 4x8 and rANS
// 4x16 blocks it decompresses via the standard library or the codec
// sub-package. For any other method it returns an "unsupported
// compression method" error. The returned length is verified against
// the block's declared UncompressedSize.
func (b *Block) Decompress() ([]byte, error) {
	var out []byte
	var err error
	switch b.Method {
	case CompRaw:
		out = b.Data
	case CompGzip:
		out, err = gunzip(b.Data)
	case CompBzip2:
		out, err = io.ReadAll(bzip2.NewReader(bytes.NewReader(b.Data)))
		if err != nil {
			err = fmt.Errorf("bzip2: %w", err)
		}
	case CompRANS4x8:
		out, err = codec.RANS4x8Decode(b.Data)
	case CompRANS4x16:
		out, err = codec.RANS4x16Decode(b.Data)
	case CompLZMA, CompArith, CompFQZComp, CompNameTok:
		return nil, fmt.Errorf("cram: unsupported compression method %d (%s)",
			byte(b.Method), b.Method)
	default:
		return nil, fmt.Errorf("cram: unsupported compression method %d", byte(b.Method))
	}
	if err != nil {
		return nil, fmt.Errorf("cram: decompressing %s block (content id %d): %w",
			b.Method, b.ContentID, err)
	}
	if int32(len(out)) != b.UncompressedSize {
		return nil, fmt.Errorf("cram: %s block (content id %d) decompressed to %d bytes, header declared %d",
			b.Method, b.ContentID, len(out), b.UncompressedSize)
	}
	return out, nil
}

// SupportedMethod reports whether Decompress can handle the block's
// compression method without returning an unsupported-method error.
func (b *Block) SupportedMethod() bool {
	switch b.Method {
	case CompRaw, CompGzip, CompBzip2, CompRANS4x8, CompRANS4x16:
		return true
	default:
		return false
	}
}

// readBlock parses one CRAM block from r: the 1-byte compression
// method, 1-byte content type, ITF-8 content id, ITF-8 compressed size,
// ITF-8 uncompressed size, the compressed data bytes, and — for CRAM
// v3+ — the trailing 4-byte little-endian CRC32. The CRC32 covers the
// block header and its data, and is validated here; a mismatch means
// the block was mis-delineated and is reported as an error.
//
// maxData is an upper bound — the number of bytes still available in
// the enclosing container — on the block's compressed payload. A block
// declaring a compressed size larger than maxData is rejected before
// any buffer is allocated, so a corrupt size field on a CRC-less CRAM
// v2 file cannot trigger an unbounded allocation.
func readBlock(r io.Reader, def FileDefinition, maxData int64) (Block, error) {
	cr := newCRCReader(r)
	var b Block
	var hdr [2]byte
	if _, err := io.ReadFull(cr, hdr[:]); err != nil {
		return b, fmt.Errorf("cram: block header: %w",
			eofToUnexpected(err, "block method/type"))
	}
	b.Method = CompressionMethod(hdr[0])
	b.ContentType = BlockContentType(hdr[1])
	var err error
	if b.ContentID, _, err = readITF8(cr); err != nil {
		return b, fmt.Errorf("cram: block content id: %w", err)
	}
	if b.CompressedSize, _, err = readITF8(cr); err != nil {
		return b, fmt.Errorf("cram: block compressed size: %w", err)
	}
	if b.UncompressedSize, _, err = readITF8(cr); err != nil {
		return b, fmt.Errorf("cram: block uncompressed size: %w", err)
	}
	if b.CompressedSize < 0 {
		return b, fmt.Errorf("cram: block declares negative compressed size %d", b.CompressedSize)
	}
	if b.UncompressedSize < 0 {
		return b, fmt.Errorf("cram: block declares negative uncompressed size %d", b.UncompressedSize)
	}
	if int64(b.CompressedSize) > maxData {
		return b, fmt.Errorf("cram: block declares compressed size %d exceeding the %d bytes left in the container",
			b.CompressedSize, maxData)
	}
	// Read the payload through io.CopyN into a growable buffer rather
	// than pre-allocating CompressedSize bytes: a corrupt size field
	// would otherwise drive an unbounded up-front allocation. CopyN
	// returns a non-nil error if the stream ends before CompressedSize
	// bytes, so truncation is reported cleanly.
	var data bytes.Buffer
	if _, cErr := io.CopyN(&data, cr, int64(b.CompressedSize)); cErr != nil {
		return b, fmt.Errorf("cram: block data (%d bytes): %w", b.CompressedSize,
			eofToUnexpected(cErr, "block data"))
	}
	b.Data = data.Bytes()
	if def.hasCRC() {
		want := cr.sum()
		var crcBuf [4]byte
		if _, err = io.ReadFull(r, crcBuf[:]); err != nil {
			return b, fmt.Errorf("cram: block CRC32: %w",
				eofToUnexpected(err, "block CRC"))
		}
		b.CRC = binary.LittleEndian.Uint32(crcBuf[:])
		if b.CRC != want {
			return b, fmt.Errorf("cram: block CRC32 mismatch (content id %d): stored %#08x, computed %#08x",
				b.ContentID, b.CRC, want)
		}
	}
	return b, nil
}

// gunzip decompresses a single gzip member.
func gunzip(in []byte) ([]byte, error) {
	zr, err := gzip.NewReader(bytes.NewReader(in))
	if err != nil {
		return nil, fmt.Errorf("gzip: %w", err)
	}
	defer zr.Close()
	out, err := io.ReadAll(zr)
	if err != nil {
		return nil, fmt.Errorf("gzip: %w", err)
	}
	return out, nil
}
