package bgzf

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"hash/crc32"
	"io"

	kflate "github.com/klauspost/compress/flate"
)

// readBufSize is the size of the buffer interposed between the caller's
// io.Reader (typically a raw *os.File) and the BGZF block parser. Each BGZF
// block is decoded by several small io.ReadFull calls (a 12-byte fixed header,
// the XLEN extra field, the ~64 KiB deflate payload and an 8-byte footer). On an
// unbuffered *os.File every one of those becomes its own read() syscall — and
// with tens of thousands of blocks that syscall churn dominates the profile
// (>80% of `samtools view` region/BED wall time was in syscall.rawsyscalln,
// almost none in the actual inflate). Interposing a 256 KiB buffered reader
// collapses those tiny reads into a handful of large sequential reads per
// buffer-fill, cutting the syscall count by orders of magnitude. 256 KiB spans
// ~4 whole BGZF blocks so a single fill usually serves several blocks. The
// buffer only ever holds bytes we are entitled to read (a chunk-bounded scan
// wraps this reader at a higher layer and stops on virtual offset, not on the
// file position), so read-ahead never over-reads across a chunk boundary in a
// way that reaches the caller.
const readBufSize = 256 * 1024

// Both the BGZF compression (writer) and decompression (reader) sides use
// klauspost/compress's pure-Go flate implementation (imported as kflate). On
// the writer side it is faster and produces a slightly better ratio than the
// standard library's compress/flate while emitting standard DEFLATE bit
// streams. On the reader side its inflater is also faster than the standard
// library's, and since BGZF blocks are ordinary DEFLATE the decoded bytes are
// identical — only the decode speed changes. BGZF read is on the hot path of
// every BAM/CRAM/.vcf.gz-consuming tool (view, flagstat, stats, depth,
// mpileup, bcftools view/query/...), so using the faster inflater there moves
// all of them toward upstream's libdeflate-backed throughput. The klauspost
// reader API mirrors the stdlib's (NewReader + Resetter), so the BGZF framing,
// block-size bounds and level mapping are unchanged.

// MaxBlockSize is the maximum number of uncompressed bytes a single BGZF block
// may carry. htslib uses 64 KiB minus a 256-byte safety margin so that the
// resulting compressed block (with the worst-case deflate expansion) still fits
// inside the 16-bit BSIZE field.
const MaxBlockSize = 65280

// MaxCompressedBlockSize is the largest legal compressed BGZF block, derived
// from the 16-bit BSIZE field plus one. A reader can refuse to allocate larger
// frames.
const MaxCompressedBlockSize = 1 << 16

// DefaultCompression matches the historical bgzip default and gzip's "level 6".
const DefaultCompression = 6

// EOFBlock is the canonical 28-byte empty BGZF block that terminates every
// well-formed BGZF stream. Decoders use its presence to distinguish a complete
// stream from a silently truncated one.
var EOFBlock = []byte{
	0x1f, 0x8b, 0x08, 0x04,
	0x00, 0x00, 0x00, 0x00,
	0x00, 0xff,
	0x06, 0x00,
	0x42, 0x43, 0x02, 0x00,
	0x1b, 0x00,
	0x03, 0x00,
	0x00, 0x00, 0x00, 0x00,
	0x00, 0x00, 0x00, 0x00,
}

// Errors returned by the BGZF reader.
var (
	// ErrBadMagic indicates the gzip magic bytes are missing.
	ErrBadMagic = errors.New("bgzf: not a gzip member (bad magic)")
	// ErrNoExtra indicates the FEXTRA flag is not set on a block header.
	ErrNoExtra = errors.New("bgzf: gzip block is missing the FEXTRA flag")
	// ErrNoBCSubfield indicates the BC subfield is absent from the extra field.
	ErrNoBCSubfield = errors.New("bgzf: gzip extra field is missing the BC subfield")
	// ErrBadBSIZE indicates BSIZE is shorter than the header it sits in.
	ErrBadBSIZE = errors.New("bgzf: BSIZE is shorter than the gzip header")
	// ErrTruncated indicates the stream ended before the BGZF EOF block.
	ErrTruncated = errors.New("bgzf: stream is truncated (missing EOF block)")
	// ErrChecksum indicates a block's CRC32 did not match its decoded payload.
	ErrChecksum = errors.New("bgzf: block CRC32 mismatch")
	// ErrISIZE indicates a block's ISIZE footer did not match the decoded length.
	ErrISIZE = errors.New("bgzf: block ISIZE mismatch")
)

// gzip header flag bits.
const (
	flagFTEXT    = 1 << 0
	flagFHCRC    = 1 << 1
	flagFEXTRA   = 1 << 2
	flagFNAME    = 1 << 3
	flagFCOMMENT = 1 << 4
)

// Writer compresses bytes into a BGZF stream. Writes are buffered up to
// MaxBlockSize bytes; each buffered chunk is emitted as one gzip member with
// the BGZF BC subfield. Close emits the final BGZF EOF block.
type Writer struct {
	w     io.Writer
	level int

	buf [MaxBlockSize]byte
	n   int

	// scratch buffer for the deflate output of one block. flate may write
	// slightly more than the input size on incompressible data, so size this
	// generously.
	deflated bytes.Buffer
	fw       *kflate.Writer

	err    error
	closed bool
}

// NewWriter returns a Writer using DefaultCompression.
func NewWriter(w io.Writer) *Writer {
	wr, _ := NewWriterLevel(w, DefaultCompression)
	return wr
}

// NewWriterLevel returns a Writer that compresses to w at the given level.
// Valid levels are flate.HuffmanOnly, flate.NoCompression (0),
// flate.BestSpeed (1) through flate.BestCompression (9), and
// flate.DefaultCompression (-1); klauspost/compress accepts the same constants.
func NewWriterLevel(w io.Writer, level int) (*Writer, error) {
	fw, err := kflate.NewWriter(io.Discard, level)
	if err != nil {
		return nil, err
	}
	return &Writer{w: w, level: level, fw: fw}, nil
}

// Write appends p to the writer's current block. Once the block fills it is
// flushed as a complete gzip member.
func (w *Writer) Write(p []byte) (int, error) {
	if w.err != nil {
		return 0, w.err
	}
	if w.closed {
		return 0, errors.New("bgzf: write on closed Writer")
	}
	total := 0
	for len(p) > 0 {
		space := MaxBlockSize - w.n
		if space == 0 {
			if err := w.flushBlock(); err != nil {
				return total, err
			}
			space = MaxBlockSize
		}
		n := copy(w.buf[w.n:], p[:min(len(p), space)])
		w.n += n
		total += n
		p = p[n:]
	}
	return total, nil
}

// Flush emits the buffered bytes as a BGZF block (even if the block is not
// full). It is safe to call Flush on an empty buffer; in that case it does
// nothing. Flush does not emit the EOF block.
func (w *Writer) Flush() error {
	if w.err != nil {
		return w.err
	}
	if w.n == 0 {
		return nil
	}
	return w.flushBlock()
}

// Close flushes any buffered bytes, writes the BGZF EOF block, and releases
// internal resources. Close does not close the underlying writer.
func (w *Writer) Close() error {
	if w.closed {
		return nil
	}
	w.closed = true
	if w.err != nil {
		return w.err
	}
	if w.n > 0 {
		if err := w.flushBlock(); err != nil {
			return err
		}
	}
	// Always emit the EOF block.
	if _, err := w.w.Write(EOFBlock); err != nil {
		w.err = err
		return err
	}
	return nil
}

// flushBlock encodes w.buf[:w.n] as a single BGZF gzip member and writes it to
// the underlying writer.
func (w *Writer) flushBlock() error {
	payload := w.buf[:w.n]
	if err := w.encodeBlock(payload); err != nil {
		w.err = err
		return err
	}
	w.n = 0
	return nil
}

// encodeBlock builds and writes one gzip member with the BC subfield for the
// given payload (which may be empty — that's how the EOF block is built, though
// EOFBlock is hard-coded for cross-implementation byte-identity).
func (w *Writer) encodeBlock(payload []byte) error {
	// Run deflate on the payload via the reusable per-Writer scratch buffer and
	// flate.Writer, then serialise the framed block to the underlying writer.
	w.deflated.Reset()
	w.fw.Reset(&w.deflated)
	if _, err := w.fw.Write(payload); err != nil {
		return err
	}
	if err := w.fw.Close(); err != nil {
		return err
	}
	frame, err := frameBlock(payload, w.deflated.Bytes(), nil)
	if err != nil {
		return err
	}
	_, err = w.w.Write(frame)
	return err
}

// frameBlock wraps an already-deflated payload in the BGZF gzip-member framing
// (fixed gzip header + BC subfield + deflate body + CRC32/ISIZE footer) and
// returns the complete on-disk block. The deflated argument must be the raw
// deflate (no gzip wrapping) of payload. dst, when non-nil and large enough, is
// reused as the destination buffer to avoid an allocation per block; otherwise
// a fresh slice is allocated.
func frameBlock(payload, deflated, dst []byte) ([]byte, error) {
	// Total compressed block size:
	//   12 bytes fixed header + 6 bytes BC subfield + len(deflated)
	//   + 8 bytes footer (CRC32 + ISIZE) = 26 + len(deflated).
	blockLen := 12 + 6 + len(deflated) + 8
	if blockLen > MaxCompressedBlockSize {
		return nil, fmt.Errorf("bgzf: compressed block size %d exceeds %d", blockLen, MaxCompressedBlockSize)
	}

	if cap(dst) >= blockLen {
		dst = dst[:blockLen]
	} else {
		dst = make([]byte, blockLen)
	}

	dst[0] = 0x1f
	dst[1] = 0x8b
	dst[2] = 8                                  // CM = deflate
	dst[3] = 0x04                               // FLG = FEXTRA
	dst[4], dst[5], dst[6], dst[7] = 0, 0, 0, 0 // MTIME
	dst[8] = 0                                  // XFL
	dst[9] = 0xff                               // OS = unknown
	// XLEN = 6 (one BC subfield: 4 bytes header + 2 bytes BSIZE)
	binary.LittleEndian.PutUint16(dst[10:12], 6)
	dst[12] = 'B'
	dst[13] = 'C'
	binary.LittleEndian.PutUint16(dst[14:16], 2)                  // SLEN = 2
	binary.LittleEndian.PutUint16(dst[16:18], uint16(blockLen-1)) // BSIZE = total-1

	copy(dst[18:18+len(deflated)], deflated)

	foot := dst[18+len(deflated):]
	binary.LittleEndian.PutUint32(foot[0:4], crc32.ChecksumIEEE(payload))
	binary.LittleEndian.PutUint32(foot[4:8], uint32(len(payload)))
	return dst, nil
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// Block is the decoded form of a single BGZF gzip member.
type Block struct {
	// CompressedSize is the on-disk size of the block (BSIZE+1).
	CompressedSize int
	// UncompressedSize is the size of the decoded payload (matches ISIZE).
	UncompressedSize int
	// Data is the decompressed payload. It is nil when only the header was
	// parsed (e.g. by Scan), and populated when the caller asks for it.
	Data []byte
}

// Reader streams a BGZF input. It decodes blocks lazily on Read and verifies
// per-block CRC32/ISIZE as it goes. A successful read of the EOF block ends
// the stream with io.EOF; if the underlying stream ends *before* the EOF block
// the reader returns ErrTruncated.
//
// VirtualOffset returns the BGZF virtual offset of the next byte that Read
// will deliver, which is suitable for use as a BAI/CSI/TBI virtual offset
// pointing at the start of a still-to-be-read record.
type Reader struct {
	// counted wraps the caller's reader and accumulates the number of
	// compressed bytes consumed so that VirtualOffset can compute the
	// compressed offset of the currently-active block.
	counted *countingReader

	// fr is reused across blocks via Reset on each new gzip member.
	fr io.ReadCloser

	// current block data and read cursor within it.
	block []byte
	off   int

	// deflated and decBuf are reused across blocks to avoid allocating a
	// fresh compressed-payload slice and decompressed-output buffer per BGZF
	// block. Both are internal: Read copies bytes out of br.block before the
	// next block is decoded, so reusing their backing arrays is safe.
	deflated []byte
	decBuf   bytes.Buffer

	// blockCoff is the compressed-stream byte offset of the current block —
	// i.e. the byte at which the gzip member containing br.block begins.
	blockCoff int64
	// nextBlockCoff is the compressed-stream byte offset of the *next* block
	// to be parsed. After nextBlock is called this becomes blockCoff for the
	// freshly decoded block.
	nextBlockCoff int64

	sawEOFBlock bool
	streamDone  bool
}

// NewReader returns a Reader that decodes BGZF bytes from r, treating the first
// byte of r as compressed-stream offset 0.
func NewReader(r io.Reader) (*Reader, error) {
	return NewReaderAt(r, 0)
}

// NewReaderAt returns a Reader that decodes BGZF bytes from r, treating the
// first byte of r as compressed-stream offset baseCoff. Use it after seeking the
// underlying stream to a BGZF block boundary so VirtualOffset reports *absolute*
// virtual offsets (block coffset << 16 | in-block uoffset). Without the base a
// reader opened mid-file numbers blocks from zero, so virtual offsets compared
// against absolute BAI/CSI chunk bounds (e.g. a chunk-bounded region scan) would
// be off by baseCoff<<16 and the bound would never fire correctly — making the
// scan over-read into later chunks and emit duplicate records.
func NewReaderAt(r io.Reader, baseCoff int64) (*Reader, error) {
	// Interpose a buffered reader so the per-block header/payload/footer
	// io.ReadFull calls are served from memory instead of one read() syscall
	// each (see readBufSize). A reader that is already a *bufio.Reader (or any
	// source that buffers internally, e.g. bytes.Reader/bytes.Buffer used by
	// in-memory tests) is left untouched to avoid a redundant copy — the
	// countingReader still tracks exactly the bytes consumed either way, so
	// virtual offsets are unchanged.
	src := r
	if bufferedReadSource(r) {
		// already cheap to read in small chunks; no wrapping benefit.
	} else {
		src = bufio.NewReaderSize(r, readBufSize)
	}
	return &Reader{counted: &countingReader{r: src}, nextBlockCoff: baseCoff}, nil
}

// bufferedReadSource reports whether r already serves small reads without a
// per-read syscall, so wrapping it in another bufio.Reader would only add a
// redundant copy. It recognises *bufio.Reader and the standard in-memory
// readers (*bytes.Reader, *bytes.Buffer, *strings.Reader) used throughout the
// tests; anything else (notably *os.File) is treated as syscall-backed and
// gets buffered.
func bufferedReadSource(r io.Reader) bool {
	switch r.(type) {
	case *bufio.Reader, *bytes.Reader, *bytes.Buffer:
		return true
	default:
		return false
	}
}

// countingReader wraps an io.Reader and tracks how many bytes have been
// returned through it. The bgzip.Reader uses this to know where in the
// compressed stream it currently is, which is what virtual offsets need.
type countingReader struct {
	r io.Reader
	n int64
}

func (c *countingReader) Read(p []byte) (int, error) {
	n, err := c.r.Read(p)
	c.n += int64(n)
	return n, err
}

// Read decodes BGZF blocks from the underlying stream into p. It returns
// io.EOF once the BGZF EOF block has been consumed; an io.ErrUnexpectedEOF or
// ErrTruncated is returned for streams that end prematurely.
func (br *Reader) Read(p []byte) (int, error) {
	if br.streamDone {
		return 0, io.EOF
	}
	total := 0
	for total < len(p) {
		if br.off >= len(br.block) {
			if br.sawEOFBlock {
				br.streamDone = true
				if total > 0 {
					return total, nil
				}
				return 0, io.EOF
			}
			if err := br.nextBlock(); err != nil {
				if err == io.EOF {
					// Stream ended without an EOF block.
					br.streamDone = true
					if total > 0 {
						return total, ErrTruncated
					}
					return 0, ErrTruncated
				}
				return total, err
			}
			continue
		}
		n := copy(p[total:], br.block[br.off:])
		br.off += n
		total += n
	}
	return total, nil
}

// nextBlock decodes one gzip member into br.block.
func (br *Reader) nextBlock() error {
	// Record where this block starts in the compressed stream before any of
	// its bytes are consumed.
	br.blockCoff = br.nextBlockCoff
	hdr, err := readBlockHeader(br.counted)
	if err != nil {
		return err
	}

	// The compressed deflate payload is between the parsed header and the
	// trailing 8-byte CRC32+ISIZE footer.
	deflatedLen := hdr.compressedSize - hdr.headerLen - 8
	if deflatedLen < 0 {
		return fmt.Errorf("bgzf: invalid block layout (deflate length %d)", deflatedLen)
	}

	if cap(br.deflated) < int(deflatedLen) {
		br.deflated = make([]byte, deflatedLen)
	} else {
		br.deflated = br.deflated[:deflatedLen]
	}
	deflated := br.deflated
	if _, err := io.ReadFull(br.counted, deflated); err != nil {
		return ioErrUnexpected(err)
	}

	var footer [8]byte
	if _, err := io.ReadFull(br.counted, footer[:]); err != nil {
		return ioErrUnexpected(err)
	}
	// Now this entire block has been consumed from the underlying reader;
	// remember where the next block will start.
	br.nextBlockCoff = br.blockCoff + int64(hdr.compressedSize)
	wantCRC := binary.LittleEndian.Uint32(footer[0:4])
	wantISIZE := binary.LittleEndian.Uint32(footer[4:8])

	if wantISIZE == 0 && deflatedLen == 2 {
		// Empty block — this is the BGZF EOF marker. The 2-byte deflate
		// stream is a single empty fixed-Huffman final block. Recognise it
		// without re-decoding to avoid edge cases.
		br.block = br.block[:0]
		br.off = 0
		br.sawEOFBlock = true
		return nil
	}

	if br.fr == nil {
		br.fr = kflate.NewReader(bytes.NewReader(deflated))
	} else {
		if rs, ok := br.fr.(kflate.Resetter); ok {
			if err := rs.Reset(bytes.NewReader(deflated), nil); err != nil {
				return err
			}
		} else {
			br.fr = kflate.NewReader(bytes.NewReader(deflated))
		}
	}

	br.decBuf.Reset()
	br.decBuf.Grow(int(wantISIZE))
	if _, err := io.Copy(&br.decBuf, br.fr); err != nil {
		return err
	}
	decoded := br.decBuf.Bytes()
	if uint32(len(decoded)) != wantISIZE {
		return ErrISIZE
	}
	if crc32.ChecksumIEEE(decoded) != wantCRC {
		return ErrChecksum
	}
	br.block = decoded
	br.off = 0
	return nil
}

// Close releases internal resources. It does not close the underlying reader.
func (br *Reader) Close() error {
	if br.fr != nil {
		err := br.fr.Close()
		br.fr = nil
		return err
	}
	return nil
}

// DecompressedRemainder returns the bytes of the current block that have been
// decoded but not yet consumed through Read — i.e. br.block[br.off:]. After the
// caller has read a stream prefix that ended mid-block (e.g. a BAM header that
// shares its BGZF block with the first records), these are the leftover
// decompressed payload bytes of that block. The returned slice aliases the
// Reader's internal buffer and is only valid until the next Read/nextBlock, so
// callers must consume or copy it immediately. It is the Go analogue of
// htslib's `in->uncompressed_block + in->block_offset` tail that samtools cat
// re-emits before switching to a verbatim compressed-block copy.
func (br *Reader) DecompressedRemainder() []byte {
	if br.off >= len(br.block) {
		return nil
	}
	return br.block[br.off:]
}

// RawRemaining returns an io.Reader over the still-unconsumed *compressed* bytes
// of the underlying stream: every BGZF block that begins after the one currently
// decoded into the Reader. Because nextBlock reads whole blocks and never reads
// ahead, the underlying stream is positioned exactly at a block boundary, so the
// bytes read from here are complete BGZF blocks (including the input's trailing
// EOF block) that can be copied verbatim to another BGZF stream without
// inflating them. It is the Go analogue of htslib's bgzf_raw_read loop in
// bam_cat. Callers must first drain DecompressedRemainder (those bytes belong to
// the block *before* this boundary); mixing Read and RawRemaining reads is not
// supported.
func (br *Reader) RawRemaining() io.Reader {
	return br.counted
}

// CopyRawTrimEOF copies compressed bytes from src to dst, dropping the final
// 28-byte BGZF EOF block. It mirrors htslib bam_cat's ebuf/es hold-back: when
// concatenating BGZF streams verbatim, every input's terminating EOF marker
// must be removed so the joined output carries exactly one EOF block (which the
// caller appends after the last input). CopyRawTrimEOF streams through a rolling
// 28-byte hold-back buffer, so it never buffers the whole input. It returns the
// number of bytes written to dst.
//
// It returns ErrTruncated when src holds fewer than 28 bytes (no EOF block, so
// the input was truncated). If the trailing 28 bytes are not the canonical EOF
// block, they are preserved (written through) rather than dropped — matching
// htslib's "unexpected block structure" behaviour of emitting the bytes with a
// warning, so no data is silently lost.
func CopyRawTrimEOF(dst io.Writer, src io.Reader) (int64, error) {
	const eof = 28 // the BGZF EOF block is always exactly 28 bytes.
	buf := make([]byte, 64*1024)
	hold := make([]byte, 0, eof+len(buf))
	var written int64
	sawAny := false
	for {
		n, rerr := src.Read(buf)
		if n > 0 {
			sawAny = true
			hold = append(hold, buf[:n]...)
			if len(hold) > eof {
				emit := len(hold) - eof
				m, werr := dst.Write(hold[:emit])
				written += int64(m)
				if werr != nil {
					return written, werr
				}
				// Retain only the trailing eof bytes for the next round.
				hold = append(hold[:0], hold[emit:]...)
			}
		}
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			return written, rerr
		}
	}
	if !sawAny || len(hold) < eof {
		return written, ErrTruncated
	}
	// hold now holds exactly the trailing EOF block; drop it — unless it is not
	// the canonical marker, in which case preserve the bytes (matches htslib's
	// warn-and-emit path so unexpected trailing data is not lost).
	if !bytes.Equal(hold, EOFBlock) {
		m, werr := dst.Write(hold)
		written += int64(m)
		if werr != nil {
			return written, werr
		}
	}
	return written, nil
}

// VirtualOffset returns the BGZF virtual offset of the next byte that Read
// will deliver. The high 48 bits are the compressed-stream offset of the
// block that owns the next byte; the low 16 bits are the uncompressed
// in-block byte position.
//
// When the current block has been fully consumed (off == len(block)) the
// returned offset points at byte 0 of the next block, which is the
// canonical form for "between records" markers in BAI/TBI/CSI.
func (br *Reader) VirtualOffset() uint64 {
	coff := br.blockCoff
	uoff := br.off
	if uoff >= len(br.block) {
		// The next byte will come from the next block at uoff 0.
		coff = br.nextBlockCoff
		uoff = 0
	}
	return uint64(coff)<<16 | uint64(uoff)&0xFFFF
}

// blockHeader holds the parsed gzip+BC header bytes.
type blockHeader struct {
	// compressedSize is total block size on disk (BSIZE+1).
	compressedSize int
	// headerLen is the number of bytes already consumed (gzip fixed header
	// + extra field + optional FNAME/FCOMMENT/FHCRC).
	headerLen int
}

// readBlockHeader reads one BGZF block header from r, validating that the
// gzip magic, deflate method, FEXTRA flag, and BC subfield are present. It
// returns the parsed BSIZE-derived total block length and the number of bytes
// consumed from r.
func readBlockHeader(r io.Reader) (blockHeader, error) {
	var fixed [12]byte
	if _, err := io.ReadFull(r, fixed[:]); err != nil {
		return blockHeader{}, err
	}
	if fixed[0] != 0x1f || fixed[1] != 0x8b {
		return blockHeader{}, ErrBadMagic
	}
	if fixed[2] != 8 {
		return blockHeader{}, fmt.Errorf("bgzf: unsupported compression method %d", fixed[2])
	}
	flg := fixed[3]
	if flg&flagFEXTRA == 0 {
		return blockHeader{}, ErrNoExtra
	}
	xlen := int(binary.LittleEndian.Uint16(fixed[10:12]))
	if xlen < 6 {
		return blockHeader{}, ErrNoBCSubfield
	}

	extra := make([]byte, xlen)
	if _, err := io.ReadFull(r, extra); err != nil {
		return blockHeader{}, ioErrUnexpected(err)
	}

	bsize, ok := findBCSubfield(extra)
	if !ok {
		return blockHeader{}, ErrNoBCSubfield
	}
	blockLen := int(bsize) + 1
	consumed := 12 + xlen

	// Skip FNAME, FCOMMENT, FHCRC if present. These bytes count against the
	// header length; the deflate payload is the remainder of the block.
	if flg&flagFNAME != 0 {
		n, err := skipCString(r)
		if err != nil {
			return blockHeader{}, ioErrUnexpected(err)
		}
		consumed += n
	}
	if flg&flagFCOMMENT != 0 {
		n, err := skipCString(r)
		if err != nil {
			return blockHeader{}, ioErrUnexpected(err)
		}
		consumed += n
	}
	if flg&flagFHCRC != 0 {
		var crc [2]byte
		if _, err := io.ReadFull(r, crc[:]); err != nil {
			return blockHeader{}, ioErrUnexpected(err)
		}
		consumed += 2
	}

	if consumed > blockLen {
		return blockHeader{}, ErrBadBSIZE
	}
	return blockHeader{compressedSize: blockLen, headerLen: consumed}, nil
}

// findBCSubfield scans the gzip extra field for the BC subfield (SI = 'B','C',
// SLEN = 2) and returns its BSIZE value. The extra field is a sequence of
// (SI1, SI2, SLEN-le16, payload[SLEN]) records.
func findBCSubfield(extra []byte) (uint16, bool) {
	for len(extra) >= 4 {
		si1 := extra[0]
		si2 := extra[1]
		slen := int(binary.LittleEndian.Uint16(extra[2:4]))
		if 4+slen > len(extra) {
			return 0, false
		}
		if si1 == 'B' && si2 == 'C' && slen == 2 {
			return binary.LittleEndian.Uint16(extra[4:6]), true
		}
		extra = extra[4+slen:]
	}
	return 0, false
}

// skipCString reads bytes from r until (and including) a terminating NUL and
// returns the number of bytes consumed.
func skipCString(r io.Reader) (int, error) {
	var b [1]byte
	n := 0
	for {
		if _, err := io.ReadFull(r, b[:]); err != nil {
			return n, err
		}
		n++
		if b[0] == 0 {
			return n, nil
		}
	}
}

func ioErrUnexpected(err error) error {
	if err == io.EOF {
		return io.ErrUnexpectedEOF
	}
	return err
}
