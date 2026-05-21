package cram

import (
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"fmt"
	"hash/crc32"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram/codec"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// Version selects the CRAM format version a RecordWriter emits. Both
// versions share the v3 container, slice and record layout and the v3
// CRC32 fields; they differ only in the per-block compression codecs the
// writer is allowed to use.
type Version int

const (
	// VersionV30 is CRAM v3.0: the writer compresses each block with raw
	// or gzip. It is the default and the format every existing caller
	// gets.
	VersionV30 Version = iota
	// VersionV31 is CRAM v3.1: in addition to raw and gzip the writer may
	// compress a block with the rANS 4x16 codec (block method 5), the
	// distinguishing capability of v3.1.
	VersionV31
)

// fileDefinition returns the on-disk file-definition version for v. Both
// CRAM v3.0 and v3.1 carry major version 3; only the minor version
// differs.
func (v Version) fileDefinition() FileDefinition {
	switch v {
	case VersionV31:
		return FileDefinition{Major: 3, Minor: 1}
	default:
		return FileDefinition{Major: 3, Minor: 0}
	}
}

// String returns the version as a "major.minor" string.
func (v Version) String() string {
	return v.fileDefinition().VersionString()
}

// defaultRecordsPerSlice caps how many records the writer packs into one
// slice (and, since the writer emits one slice per container, one
// container). A new container is started once this many records have
// been buffered. The value is a balance between container-header
// overhead and per-slice memory; it has no effect on correctness.
const defaultRecordsPerSlice = 10000

// gzipBlockThreshold is the smallest external-block payload the writer
// bothers to gzip-compress. Below it the deflate header and trailer cost
// more than they save, so the block is written raw (method 0).
const gzipBlockThreshold = 32

// content-id constants. Every external data series the writer emits is
// keyed by a stable small integer; the compression header's encoding map
// points each two-character series at its block by this id. The values
// are an implementation detail — any distinct positive integers would
// do — but keeping them fixed makes the on-disk layout predictable and
// the writer easy to read against the reader.
const (
	cidSAMHeader = 0  // the SAM file-header block.
	cidBF        = 1  // bit flags (SAM FLAG).
	cidCF        = 2  // CRAM per-record flags.
	cidRI        = 3  // reference id (multi-reference slices).
	cidRL        = 4  // read length.
	cidAP        = 5  // alignment position.
	cidRG        = 6  // read group (always -1; RG travels as a tag).
	cidRN        = 7  // read names.
	cidMF        = 8  // mate flags.
	cidNS        = 9  // mate reference id.
	cidNP        = 10 // mate alignment position.
	cidTS        = 11 // template size.
	cidTL        = 12 // tag-line (tag-combination index).
	cidMQ        = 13 // mapping quality.
	cidFN        = 14 // read-feature count.
	cidFC        = 15 // read-feature codes.
	cidFP        = 16 // read-feature positions.
	cidBB        = 17 // base-stretch bytes (feature 'b').
	cidBBLen     = 18 // base-stretch lengths.
	cidBA        = 19 // single bases (unmapped reads).
	cidQS        = 20 // quality scores.
	cidIN        = 21 // inserted bases (feature 'I').
	cidINLen     = 22 // inserted-base lengths.
	cidSC        = 23 // soft-clipped bases (feature 'S').
	cidSCLen     = 24 // soft-clipped lengths.
	cidDL        = 25 // deletion lengths (feature 'D').
	cidRS        = 26 // reference-skip lengths (feature 'N').
	cidPD        = 27 // padding lengths (feature 'P').
	cidHC        = 28 // hard-clip lengths (feature 'H').
	cidTagBase   = 64 // first content id handed out to auxiliary tag series.
)

// tagContentIDs returns the two content ids an auxiliary tag's
// BYTE_ARRAY_LEN series uses: a length block and a value block. The
// i-th distinct tag gets the pair (cidTagBase+2i, cidTagBase+2i+1).
func tagContentIDs(i int) (lenID, valID int32) {
	base := cidTagBase + int32(i)*2
	return base, base + 1
}

// RecordWriter encodes a stream of sam.Record values as a CRAM v3.0
// file. It mirrors RecordReader: construct one with NewRecordWriter,
// feed records with Write, and finalise the file with Close.
//
// The writer produces reference-free CRAM — a mapped read's bases are
// stored literally rather than diffed against an external reference — so
// the file is fully self-contained and decodes without a reference
// FASTA. Records are buffered and flushed one slice (and one container)
// at a time; Close flushes the final partial slice and appends the CRAM
// EOF marker.
//
// A RecordWriter is not safe for concurrent use.
type RecordWriter struct {
	w      io.Writer
	closer io.Closer
	header *sam.Header

	// version is the CRAM format this writer emits. It is a per-writer
	// field — not package-level state — so two writers targeting
	// different versions can run concurrently. It governs both the
	// file-definition minor version and the per-block codec set.
	version Version

	// refIndex maps a reference name to its zero-based @SQ position, so a
	// record's RName / RNext can be turned into the integer ids the CRAM
	// data series store.
	refIndex map[string]int32

	// buf holds the records of the slice currently being assembled; it is
	// flushed to a container once it reaches recordsPerSlice entries.
	buf []*sam.Record
	// recordsPerSlice is the slice-size cap; defaultRecordsPerSlice
	// unless overridden for testing.
	recordsPerSlice int
	// recordCounter is the running total of records written to previous
	// containers, stored in each container and slice header.
	recordCounter int64

	// wroteHeader records whether the file-definition and SAM-header
	// container have been emitted yet.
	wroteHeader bool
	// closed records whether Close has run, so a double Close is a no-op.
	closed bool
	// err latches the first write error; every later call returns it.
	err error
}

// NewRecordWriter returns a RecordWriter that encodes records to w as a
// CRAM v3.0 file. The SAM header is written immediately as the first
// container's single block, so h must be complete before the call. It
// returns an error only if the initial header write to w fails.
//
// To target CRAM v3.1 instead, use NewRecordWriterVersion.
func NewRecordWriter(w io.Writer, h *sam.Header) (*RecordWriter, error) {
	return NewRecordWriterVersion(w, h, VersionV30)
}

// NewRecordWriterVersion returns a RecordWriter that encodes records to w
// as a CRAM file of the requested version. VersionV30 produces a v3.0
// file (raw/gzip block compression); VersionV31 produces a v3.1 file,
// which may additionally compress blocks with the rANS 4x16 codec. The
// SAM header is written immediately, so h must be complete before the
// call. It returns an error only if the initial header write to w fails.
func NewRecordWriterVersion(w io.Writer, h *sam.Header, version Version) (*RecordWriter, error) {
	if h == nil {
		h = &sam.Header{}
	}
	rw := &RecordWriter{
		w:               w,
		header:          h,
		version:         version,
		refIndex:        make(map[string]int32, len(h.Refs)),
		recordsPerSlice: defaultRecordsPerSlice,
	}
	for i, ref := range h.Refs {
		rw.refIndex[ref.Name] = int32(i)
	}
	if err := rw.writeFileHeader(); err != nil {
		return nil, err
	}
	return rw, nil
}

// CreateCRAM creates the named file and returns a RecordWriter over it
// targeting CRAM v3.0. The caller must call Close, which flushes the
// final container and releases the file handle.
func CreateCRAM(path string, h *sam.Header) (*RecordWriter, error) {
	return CreateCRAMVersion(path, h, VersionV30)
}

// CreateCRAMVersion creates the named file and returns a RecordWriter
// over it targeting the requested CRAM version. The caller must call
// Close, which flushes the final container and releases the file handle.
func CreateCRAMVersion(path string, h *sam.Header, version Version) (*RecordWriter, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	rw, err := NewRecordWriterVersion(f, h, version)
	if err != nil {
		f.Close()
		os.Remove(path)
		return nil, err
	}
	rw.closer = f
	return rw, nil
}

// WriteCRAM encodes header and every record in records to w as a
// complete CRAM v3.0 file, including the trailing EOF marker. It is a
// convenience wrapper over NewRecordWriter / Write / Close.
func WriteCRAM(w io.Writer, h *sam.Header, records []*sam.Record) error {
	return WriteCRAMVersion(w, h, records, VersionV30)
}

// WriteCRAMVersion encodes header and every record in records to w as a
// complete CRAM file of the requested version, including the trailing
// EOF marker. It is a convenience wrapper over NewRecordWriterVersion /
// Write / Close.
func WriteCRAMVersion(w io.Writer, h *sam.Header, records []*sam.Record, version Version) error {
	rw, err := NewRecordWriterVersion(w, h, version)
	if err != nil {
		return err
	}
	for i, rec := range records {
		if err := rw.Write(rec); err != nil {
			return fmt.Errorf("cram: writing record %d: %w", i, err)
		}
	}
	return rw.Close()
}

// Write buffers one record for encoding. Records accumulate until a
// slice is full, at which point a container is assembled and flushed to
// the underlying writer. A record shape the simple writer cannot encode
// (see the package documentation) is rejected here with a clear error
// and nothing is written.
func (rw *RecordWriter) Write(rec *sam.Record) error {
	if rw.err != nil {
		return rw.err
	}
	if rw.closed {
		return fmt.Errorf("cram: Write after Close")
	}
	if rec == nil {
		return fmt.Errorf("cram: cannot write a nil record")
	}
	if err := rw.checkRecord(rec); err != nil {
		return err
	}
	rw.buf = append(rw.buf, rec)
	if len(rw.buf) >= rw.recordsPerSlice {
		return rw.flushContainer()
	}
	return nil
}

// Close flushes the final partial slice, appends the CRAM v3 EOF marker,
// and releases the file handle when the writer was created by
// CreateCRAM. It is safe to call Close more than once.
func (rw *RecordWriter) Close() error {
	if rw.closed {
		return rw.err
	}
	rw.closed = true
	if rw.err == nil && len(rw.buf) > 0 {
		rw.flushContainer()
	}
	if rw.err == nil {
		if _, err := rw.w.Write(eofMarkerV3); err != nil {
			rw.err = fmt.Errorf("cram: writing EOF marker: %w", err)
		}
	}
	if rw.closer != nil {
		if cerr := rw.closer.Close(); cerr != nil && rw.err == nil {
			rw.err = cerr
		}
	}
	return rw.err
}

// checkRecord rejects a record the simple reference-free writer cannot
// faithfully encode. The writer aims for correctness over coverage: a
// shape it cannot round-trip is an explicit error, never silent data
// loss.
func (rw *RecordWriter) checkRecord(rec *sam.Record) error {
	if rec.RName != "" && rec.RName != "*" {
		if _, ok := rw.refIndex[rec.RName]; !ok {
			return fmt.Errorf("cram: record %q references %q, absent from the SAM header @SQ lines",
				rec.QName, rec.RName)
		}
	}
	if rec.RNext != "" && rec.RNext != "*" && rec.RNext != "=" {
		if _, ok := rw.refIndex[rec.RNext]; !ok {
			return fmt.Errorf("cram: record %q mate references %q, absent from the SAM header @SQ lines",
				rec.QName, rec.RNext)
		}
	}
	if rec.Flag&sam.FlagUnmapped == 0 {
		// A mapped record's CIGAR must consume exactly the SEQ length, or
		// the read-feature reconstruction cannot recover the sequence.
		if len(rec.Cigar) > 0 && rec.Seq != "" && rec.Seq != "*" {
			if rec.Cigar.QueryLength() != len(rec.Seq) {
				return fmt.Errorf("cram: record %q CIGAR query length %d does not match SEQ length %d",
					rec.QName, rec.Cigar.QueryLength(), len(rec.Seq))
			}
		}
	}
	for _, a := range rec.Aux {
		if _, err := encodeTagValue(a); err != nil {
			return fmt.Errorf("cram: record %q tag %q: %w", rec.QName, a.Tag, err)
		}
	}
	return nil
}

// writeFileHeader emits the 26-byte CRAM file definition and the
// file-header container that carries the SAM text header as its single
// block.
func (rw *RecordWriter) writeFileHeader() error {
	if rw.wroteHeader {
		return nil
	}
	fd := rw.version.fileDefinition()
	var def [fileDefSize]byte
	copy(def[0:4], fileDefMagic[:])
	def[4] = fd.Major
	def[5] = fd.Minor
	// FileID is left NUL — the conventional originating-file-name slot is
	// optional and a reader trims trailing NULs.
	if _, err := rw.w.Write(def[:]); err != nil {
		return fmt.Errorf("cram: writing file definition: %w", err)
	}

	// The SAM header block payload is a 4-byte little-endian text length
	// followed by the header text, matching what readSAMHeader expects.
	text := []byte(rw.header.Text())
	payload := make([]byte, 4+len(text))
	binary.LittleEndian.PutUint32(payload[:4], uint32(len(text)))
	copy(payload[4:], text)
	headerBlock := encodeBlock(rw.version, ContentFileHeader, cidSAMHeader, payload)

	// The file-header container holds exactly the header block. Its
	// reference fields are the unmapped/no-data sentinels.
	body := headerBlock
	hdr := containerHeaderBytes(containerFields{
		length:        int32(len(body)),
		refSeqID:      0,
		startPos:      0,
		alignmentSpan: 0,
		numRecords:    0,
		recordCounter: 0,
		numBases:      0,
		numBlocks:     1,
		landmarks:     nil,
	})
	if _, err := rw.w.Write(hdr); err != nil {
		return fmt.Errorf("cram: writing file-header container header: %w", err)
	}
	if _, err := rw.w.Write(body); err != nil {
		return fmt.Errorf("cram: writing file-header container body: %w", err)
	}
	rw.wroteHeader = true
	return nil
}

// flushContainer encodes every buffered record into one data container
// (a compression-header block plus a single slice) and writes it to the
// underlying writer. The buffer is cleared and the record counter
// advanced on success.
func (rw *RecordWriter) flushContainer() error {
	if rw.err != nil {
		return rw.err
	}
	if len(rw.buf) == 0 {
		return nil
	}
	container, err := encodeContainer(rw.version, rw.buf, rw.refIndex, rw.recordCounter)
	if err != nil {
		rw.err = err
		return err
	}
	if _, err := rw.w.Write(container); err != nil {
		rw.err = fmt.Errorf("cram: writing container: %w", err)
		return rw.err
	}
	rw.recordCounter += int64(len(rw.buf))
	rw.buf = rw.buf[:0]
	return nil
}

// encodeBlock assembles a complete on-disk CRAM v3 block from a content
// type, content id and uncompressed payload. The payload is compressed
// with whichever method chooseBlockCompression picks for the writer's
// version — never larger than raw — and the trailing IEEE CRC32 over the
// whole block is appended.
func encodeBlock(version Version, ct BlockContentType, contentID int32, payload []byte) []byte {
	method, stored := chooseBlockCompression(version, payload)
	var b []byte
	b = append(b, byte(method), byte(ct))
	b = appendITF8(b, contentID)
	b = appendITF8(b, int32(len(stored)))
	b = appendITF8(b, int32(len(payload)))
	b = append(b, stored...)
	crc := crc32.Checksum(b, crc32.IEEETable)
	var crcBuf [4]byte
	binary.LittleEndian.PutUint32(crcBuf[:], crc)
	return append(b, crcBuf[:]...)
}

// chooseBlockCompression picks the block compression method for a
// payload and returns it together with the bytes to store. The candidate
// set depends on the CRAM version: v3.0 considers raw (method 0) and
// gzip (method 1); v3.1 additionally considers rANS 4x16 (method 5), its
// distinguishing codec. The smallest candidate wins and raw is always in
// the running, so the stored payload is never larger than the input —
// the block stays decodable even if a codec misbehaves.
func chooseBlockCompression(version Version, payload []byte) (CompressionMethod, []byte) {
	method := CompRaw
	stored := payload
	// Below the threshold the per-codec framing overhead outweighs any
	// saving, so a tiny block is always stored raw.
	if len(payload) < gzipBlockThreshold {
		return method, stored
	}
	if gz := gzipCompress(payload); len(gz) < len(stored) {
		method, stored = CompGzip, gz
	}
	if version == VersionV31 {
		// rANS 4x16 is a v3.1-only codec; never offer it for v3.0. Order 0
		// is used — correctness and round-trip matter here, not ratio, and
		// the order-0 model round-trips every input including the empty
		// one.
		if r, err := codec.RANS4x16Encode(payload, 0); err == nil && len(r) < len(stored) {
			method, stored = CompRANS4x16, r
		}
	}
	return method, stored
}

// gzipCompress returns the gzip (RFC 1952) compression of in. It is the
// encode-side counterpart of block.go's gunzip.
func gzipCompress(in []byte) []byte {
	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	if _, err := zw.Write(in); err != nil {
		// gzip.Writer over a bytes.Buffer cannot fail; fall back to raw.
		return in
	}
	if err := zw.Close(); err != nil {
		return in
	}
	return buf.Bytes()
}

// containerFields collects the variable-length fields of a CRAM
// container header so containerHeaderBytes can serialise them.
type containerFields struct {
	length        int32
	refSeqID      int32
	startPos      int32
	alignmentSpan int32
	numRecords    int32
	recordCounter int64
	numBases      int64
	numBlocks     int32
	landmarks     []int32
}

// containerHeaderBytes serialises a CRAM v3 container header: the fixed
// 4-byte little-endian length, the ITF-8/LTF-8 fields, the landmark
// array, and the trailing IEEE CRC32 over all of it. It is the writer-
// side inverse of parseContainerHeader.
func containerHeaderBytes(f containerFields) []byte {
	var lenBuf [4]byte
	binary.LittleEndian.PutUint32(lenBuf[:], uint32(f.length))
	b := append([]byte(nil), lenBuf[:]...)
	b = appendITF8(b, f.refSeqID)
	b = appendITF8(b, f.startPos)
	b = appendITF8(b, f.alignmentSpan)
	b = appendITF8(b, f.numRecords)
	b = appendLTF8(b, f.recordCounter)
	b = appendLTF8(b, f.numBases)
	b = appendITF8(b, f.numBlocks)
	b = appendITF8(b, int32(len(f.landmarks)))
	for _, lm := range f.landmarks {
		b = appendITF8(b, lm)
	}
	crc := crc32.Checksum(b, crc32.IEEETable)
	var crcBuf [4]byte
	binary.LittleEndian.PutUint32(crcBuf[:], crc)
	return append(b, crcBuf[:]...)
}
