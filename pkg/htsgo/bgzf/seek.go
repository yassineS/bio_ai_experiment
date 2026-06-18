package bgzf

import (
	"bytes"
	"encoding/binary"
	"errors"
	kflate "github.com/klauspost/compress/flate"
	"hash/crc32"
	"io"
	"sort"
)

// ErrSeekPastEnd indicates that a requested uncompressed offset is beyond the
// end of the decompressed stream.
var ErrSeekPastEnd = errors.New("bgzf: seek past end of uncompressed stream")

// SeekReader provides random access into a BGZF stream using a .gzi block
// index, decompressing only the blocks that overlap the requested region.
// It mirrors htslib's bgzf_useek + the `bgzip -b N -s M` workflow: rather than
// decompressing the whole file, it binary-searches the block index for the
// block that owns a given uncompressed offset, seeks the underlying reader to
// that block, and inflates from there.
//
// A SeekReader is not safe for concurrent use; create one per goroutine (the
// underlying io.ReaderAt may be shared since ReadAt is concurrency-safe).
type SeekReader struct {
	r io.ReaderAt

	// blocks is the in-memory block index, always including the implicit
	// leading (0,0) entry as blocks[0], sorted by UncompressedOffset. Each
	// entry's CompressedOffset/UncompressedOffset are populated; the size
	// fields are derived lazily and may be zero for a .gzi-loaded index.
	blocks []BlockOffset

	// curBlock holds the decompressed payload of the block currently loaded;
	// curUAddr is that block's uncompressed start offset and curCAddr is its
	// compressed start offset. curCAddr is tracked explicitly (rather than
	// re-derived from the index) so advanceBlock works even when the current
	// block is not named by the index — i.e. one reached via the sparse-index
	// walk-forward fallback.
	curBlock []byte
	curUAddr int64
	curCAddr int64
	loaded   bool

	// pos is the current uncompressed read position in the virtual stream.
	pos int64
}

// NewSeekReader builds a SeekReader over r using the block index produced by
// ReadGZI (which omits the implicit leading (0,0) entry — NewSeekReader adds
// it back). The index entries must be sorted by ascending CompressedOffset, as
// htslib and WriteGZI both emit them; NewSeekReader re-sorts defensively by
// UncompressedOffset so the binary search is correct regardless of input order.
//
// r must remain valid for the lifetime of the SeekReader.
func NewSeekReader(r io.ReaderAt, index []BlockOffset) *SeekReader {
	blocks := make([]BlockOffset, 0, len(index)+1)
	// Always prepend the implicit first block at (0,0) unless the caller
	// already supplied it.
	if len(index) == 0 || index[0].CompressedOffset != 0 || index[0].UncompressedOffset != 0 {
		blocks = append(blocks, BlockOffset{CompressedOffset: 0, UncompressedOffset: 0})
	}
	blocks = append(blocks, index...)
	sort.Slice(blocks, func(i, j int) bool {
		return blocks[i].UncompressedOffset < blocks[j].UncompressedOffset
	})
	return &SeekReader{r: r, blocks: blocks, pos: 0}
}

// blockIndexFor returns the index into sr.blocks of the block that owns the
// uncompressed offset uoffset — i.e. the block with the greatest
// UncompressedOffset that is <= uoffset. This is the Go equivalent of the
// binary search in htslib's bgzf_useek.
func (sr *SeekReader) blockIndexFor(uoffset int64) int {
	// sort.Search finds the first block whose UncompressedOffset > uoffset;
	// the owning block is the one before it.
	i := sort.Search(len(sr.blocks), func(i int) bool {
		return sr.blocks[i].UncompressedOffset > uoffset
	})
	if i == 0 {
		return 0
	}
	return i - 1
}

// SeekTo positions the reader at the given uncompressed-stream offset. Only
// absolute (SEEK_SET) semantics are supported. It decompresses the owning
// block if it is not already loaded but does not eagerly read past it. After
// SeekTo, Read delivers bytes starting at uoffset. (The method is named SeekTo
// rather than Seek so it does not collide with the io.Seeker signature, which
// this type intentionally does not implement.)
//
// Seeking exactly to the end of the stream is allowed and leaves Read at EOF.
func (sr *SeekReader) SeekTo(uoffset int64) error {
	if uoffset < 0 {
		return errors.New("bgzf: negative seek offset")
	}
	idx := sr.blockIndexFor(uoffset)
	target := sr.blocks[idx]
	// Load the nearest indexed block at or before the target if a different
	// one is currently loaded.
	if !sr.loaded || sr.curUAddr != target.UncompressedOffset {
		if err := sr.loadBlock(target.CompressedOffset, target.UncompressedOffset); err != nil {
			return err
		}
	}
	// The index may be sparse (it need not name the block that actually owns
	// uoffset), so walk forward block-by-block until the loaded block contains
	// the offset. The owning block is the one where uoffset lies in
	// [curUAddr, curUAddr+len(curBlock)); seeking exactly to the end of the
	// final block is allowed and leaves Read at EOF.
	for uoffset >= sr.curUAddr+int64(len(sr.curBlock)) {
		err := sr.advanceBlock()
		if err == io.EOF {
			// No more blocks. Only the exact end-of-stream offset is valid.
			if uoffset == sr.curUAddr+int64(len(sr.curBlock)) {
				sr.pos = uoffset
				return nil
			}
			return ErrSeekPastEnd
		}
		if err != nil {
			return err
		}
	}
	sr.pos = uoffset
	return nil
}

// loadBlock seeks the underlying reader to the compressed offset caddr and
// decompresses exactly one BGZF block, storing the payload in sr.curBlock and
// recording its uncompressed start uaddr.
func (sr *SeekReader) loadBlock(caddr, uaddr int64) error {
	hdr, err := readBlockHeaderAt(sr.r, caddr)
	if err != nil {
		return err
	}
	deflatedLen := hdr.compressedSize - hdr.headerLen - 8
	if deflatedLen < 0 {
		return errors.New("bgzf: invalid block layout in seek")
	}
	// The deflate payload starts at caddr+headerLen; the footer follows.
	frame := make([]byte, deflatedLen+8)
	// io.ReaderAt may return io.EOF together with a full read, so only treat
	// the read as failed when fewer than len(frame) bytes came back.
	if n, err := sr.r.ReadAt(frame, caddr+int64(hdr.headerLen)); err != nil && n < len(frame) {
		return ioErrUnexpected(err)
	}
	deflated := frame[:deflatedLen]
	footer := frame[deflatedLen:]
	wantCRC := binary.LittleEndian.Uint32(footer[0:4])
	wantISIZE := binary.LittleEndian.Uint32(footer[4:8])

	if wantISIZE == 0 && deflatedLen == 2 {
		// EOF marker block — empty payload.
		sr.curBlock = sr.curBlock[:0]
		sr.curUAddr = uaddr
		sr.curCAddr = caddr
		sr.loaded = true
		return nil
	}

	fr := kflate.NewReader(bytes.NewReader(deflated))
	defer fr.Close()
	out := make([]byte, 0, wantISIZE)
	buf := bytes.NewBuffer(out)
	if _, err := io.Copy(buf, fr); err != nil {
		return err
	}
	decoded := buf.Bytes()
	if uint32(len(decoded)) != wantISIZE {
		return ErrISIZE
	}
	if crc32.ChecksumIEEE(decoded) != wantCRC {
		return ErrChecksum
	}
	sr.curBlock = decoded
	sr.curUAddr = uaddr
	sr.curCAddr = caddr
	sr.loaded = true
	return nil
}

// Read fills p with decompressed bytes starting at the current position,
// advancing across block boundaries as needed. It returns io.EOF once the
// stream is exhausted. The block index need not enumerate every block: when
// the next block's offset is not in the index, Read walks the compressed
// stream forward from the last known block to find it.
func (sr *SeekReader) Read(p []byte) (int, error) {
	if len(p) == 0 {
		return 0, nil
	}
	if !sr.loaded {
		if err := sr.SeekTo(sr.pos); err != nil {
			return 0, err
		}
	}
	total := 0
	for total < len(p) {
		inBlock := sr.pos - sr.curUAddr
		if inBlock < int64(len(sr.curBlock)) {
			n := copy(p[total:], sr.curBlock[inBlock:])
			total += n
			sr.pos += int64(n)
			continue
		}
		// Current block exhausted — advance to the next block.
		if err := sr.advanceBlock(); err != nil {
			if err == io.EOF {
				if total > 0 {
					return total, nil
				}
				return 0, io.EOF
			}
			return total, err
		}
	}
	return total, nil
}

// advanceBlock loads the block immediately following the currently-loaded one.
// It prefers the index entry whose uncompressed offset equals the end of the
// current block; if the index is sparse (does not name that block), it falls
// back to reading the next block header directly from the compressed stream at
// the byte just past the current block.
func (sr *SeekReader) advanceBlock() error {
	nextUAddr := sr.curUAddr + int64(len(sr.curBlock))
	// If an index entry exactly names the next uncompressed offset, use it.
	// (blockIndexFor finds the entry owning curUAddr; the one after it is the
	// candidate successor.)
	idx := sr.blockIndexFor(sr.curUAddr)
	if idx+1 < len(sr.blocks) && sr.blocks[idx+1].UncompressedOffset == nextUAddr {
		nb := sr.blocks[idx+1]
		return sr.loadBlock(nb.CompressedOffset, nb.UncompressedOffset)
	}
	// Otherwise, derive the next block's compressed offset by re-reading the
	// current block's header to learn its on-disk size. We must use the current
	// block's true compressed offset (tracked in curCAddr), not an index lookup
	// — the current block may itself be one that the (sparse) index does not
	// name, in which case blockIndexFor(curUAddr) would point at an earlier
	// block and we would never make progress.
	curCAddr := sr.curCAddr
	hdr, err := readBlockHeaderAt(sr.r, curCAddr)
	if err != nil {
		return err
	}
	nextCAddr := curCAddr + int64(hdr.compressedSize)
	// Detect the EOF block (empty payload) to surface io.EOF cleanly.
	probe, perr := readBlockHeaderAt(sr.r, nextCAddr)
	if perr == io.EOF {
		return io.EOF
	}
	if perr != nil {
		return perr
	}
	if probe.compressedSize == len(EOFBlock) {
		// Could be the EOF marker; load it and check for empty payload.
		if err := sr.loadBlock(nextCAddr, nextUAddr); err != nil {
			return err
		}
		if len(sr.curBlock) == 0 {
			return io.EOF
		}
		return nil
	}
	return sr.loadBlock(nextCAddr, nextUAddr)
}

// ReadRegion is the direct equivalent of `bgzip -b uoffset -s size`: it returns
// exactly size decompressed bytes starting at uncompressed offset uoffset
// (fewer only if the stream ends first). A size <= 0 reads to end of stream.
func (sr *SeekReader) ReadRegion(uoffset, size int64) ([]byte, error) {
	if err := sr.SeekTo(uoffset); err != nil {
		return nil, err
	}
	if size <= 0 {
		return io.ReadAll(sr)
	}
	out := make([]byte, size)
	n, err := io.ReadFull(sr, out)
	if err == io.EOF || err == io.ErrUnexpectedEOF {
		// Stream ended before `size` bytes — return what we got.
		return out[:n], nil
	}
	if err != nil {
		return nil, err
	}
	return out[:n], nil
}

// Offset returns the current uncompressed-stream read position.
func (sr *SeekReader) Offset() int64 { return sr.pos }

// readBlockHeaderAt reads and parses a BGZF block header at the given
// compressed offset using ReadAt, returning the same blockHeader that the
// streaming reader produces. It reads a generous fixed window (enough for the
// fixed header plus any plausible extra/optional fields) and parses from it.
func readBlockHeaderAt(r io.ReaderAt, caddr int64) (blockHeader, error) {
	// The BGZF block header is at most 18 bytes for the canonical layout, but
	// gzip permits FNAME/FCOMMENT of arbitrary length. We first read the
	// fixed 12 bytes to learn XLEN, then read exactly what we need.
	var fixed [12]byte
	if _, err := r.ReadAt(fixed[:], caddr); err != nil {
		if err == io.EOF {
			return blockHeader{}, io.EOF
		}
		return blockHeader{}, ioErrUnexpected(err)
	}
	if fixed[0] != 0x1f || fixed[1] != 0x8b {
		return blockHeader{}, ErrBadMagic
	}
	if fixed[2] != 8 {
		return blockHeader{}, ErrBadMagic
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
	if _, err := r.ReadAt(extra, caddr+12); err != nil {
		return blockHeader{}, ioErrUnexpected(err)
	}
	bsize, ok := findBCSubfield(extra)
	if !ok {
		return blockHeader{}, ErrNoBCSubfield
	}
	blockLen := int(bsize) + 1
	consumed := 12 + xlen
	// Skip FNAME/FCOMMENT/FHCRC if present by scanning from disk.
	pos := caddr + int64(consumed)
	if flg&flagFNAME != 0 {
		n, err := skipCStringAt(r, pos)
		if err != nil {
			return blockHeader{}, ioErrUnexpected(err)
		}
		consumed += n
		pos += int64(n)
	}
	if flg&flagFCOMMENT != 0 {
		n, err := skipCStringAt(r, pos)
		if err != nil {
			return blockHeader{}, ioErrUnexpected(err)
		}
		consumed += n
		pos += int64(n)
	}
	if flg&flagFHCRC != 0 {
		consumed += 2
	}
	if consumed > blockLen {
		return blockHeader{}, ErrBadBSIZE
	}
	return blockHeader{compressedSize: blockLen, headerLen: consumed}, nil
}

// skipCStringAt counts the bytes of a NUL-terminated string starting at pos in
// r (including the terminator), reading one byte at a time via ReadAt.
func skipCStringAt(r io.ReaderAt, pos int64) (int, error) {
	var b [1]byte
	n := 0
	for {
		if _, err := r.ReadAt(b[:], pos+int64(n)); err != nil {
			return n, err
		}
		n++
		if b[0] == 0 {
			return n, nil
		}
	}
}
