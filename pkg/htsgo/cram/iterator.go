package cram

import (
	"encoding/binary"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// RecordReader walks a CRAM stream record by record, reconstructing each
// CRAM alignment record into a sam.Record. It builds on the structural
// Reader and the C4a data-series layer: it parses every data container,
// decodes each slice's data series through the per-record traversal, and
// yields the reconstructed records in file order.
//
// A RecordReader is created with NewRecordReader, which also parses the
// embedded SAM header (available via Header). It is not safe for
// concurrent use; Read advances shared state.
type RecordReader struct {
	rd     *Reader
	header *sam.Header

	// refNames and readGroups are the @SQ names and @RG IDs the CRAM
	// data series index into, derived once from the header.
	refNames   []string
	readGroups []string

	// pending holds records decoded from the current slice that have not
	// yet been returned by Read; a slice is decoded in one shot so that
	// its interleaved data series are read in a single consistent pass.
	pending []*sam.Record
	// next is the index of the next pending record to return.
	next int
	// done is set once the stream is exhausted.
	done bool
	// recordIndex is the running record number, used to synthesise read
	// names when the file did not preserve them.
	recordIndex int64
	// needsReference is set once any decoded record reached a base an
	// external reference would supply.
	needsReference bool
}

// NewRecordReader reads the CRAM file definition and the embedded SAM
// header from r and returns a RecordReader positioned before the first
// alignment record. It returns an error if r is not a CRAM stream or the
// embedded header cannot be parsed.
func NewRecordReader(r io.Reader) (*RecordReader, error) {
	rd, err := NewReader(r)
	if err != nil {
		return nil, err
	}
	rr := &RecordReader{rd: rd}
	if err := rr.readSAMHeader(); err != nil {
		return nil, err
	}
	return rr, nil
}

// OpenRecords opens the named CRAM file and returns a RecordReader over
// it. The caller must call Close to release the file handle.
func OpenRecords(path string) (*RecordReader, error) {
	rd, err := Open(path)
	if err != nil {
		return nil, err
	}
	rr := &RecordReader{rd: rd}
	if err := rr.readSAMHeader(); err != nil {
		rd.Close()
		return nil, err
	}
	return rr, nil
}

// Close releases the underlying CRAM Reader's file handle, if any.
func (rr *RecordReader) Close() error { return rr.rd.Close() }

// Header returns the SAM header parsed from the CRAM file's first
// container. The header is available immediately after NewRecordReader.
func (rr *RecordReader) Header() *sam.Header { return rr.header }

// readSAMHeader reads the CRAM file's first container — the file-header
// container — and parses the SAM text header it carries. The header
// block payload is a 4-byte little-endian text length followed by that
// many bytes of SAM header text.
func (rr *RecordReader) readSAMHeader() error {
	c, err := rr.rd.Next()
	if err != nil {
		return wrapf(err, "reading the CRAM file-header container")
	}
	if len(c.Blocks) == 0 {
		return errFormat("the CRAM file-header container has no blocks")
	}
	first := &c.Blocks[0]
	if first.ContentType != ContentFileHeader {
		return errFormat("the first CRAM block is %s, not a SAM file header", first.ContentType)
	}
	payload, err := first.Decompress()
	if err != nil {
		return wrapf(err, "decompressing the SAM header block")
	}
	if len(payload) < 4 {
		return errFormat("SAM header block is %d bytes, too short for a 4-byte length prefix", len(payload))
	}
	textLen := binary.LittleEndian.Uint32(payload[:4])
	if int(textLen) > len(payload)-4 {
		return errFormat("SAM header declares %d text bytes but only %d follow the prefix", textLen, len(payload)-4)
	}
	h, err := sam.ParseHeaderText(string(payload[4 : 4+int(textLen)]))
	if err != nil {
		return wrapf(err, "parsing the embedded SAM header")
	}
	rr.header = h
	rr.refNames = make([]string, len(h.Refs))
	for i, ref := range h.Refs {
		rr.refNames[i] = ref.Name
	}
	rr.readGroups = make([]string, len(h.ReadGroups))
	for i, rg := range h.ReadGroups {
		rr.readGroups[i] = rg.ID
	}
	return nil
}

// Read returns the next reconstructed alignment record, or io.EOF when
// the stream is exhausted. Records are decoded a slice at a time so an
// error mid-slice is reported on the Read that first reaches it.
func (rr *RecordReader) Read() (*sam.Record, error) {
	for {
		if rr.next < len(rr.pending) {
			rec := rr.pending[rr.next]
			rr.next++
			return rec, nil
		}
		if rr.done {
			return nil, io.EOF
		}
		if err := rr.fillNextSlice(); err != nil {
			return nil, err
		}
	}
}

// ReadAll reads and returns every remaining record. It is a convenience
// wrapper over repeated Read calls; it returns whatever records it
// decoded alongside the first error encountered before io.EOF.
func (rr *RecordReader) ReadAll() ([]*sam.Record, error) {
	var out []*sam.Record
	for {
		rec, err := rr.Read()
		if err == io.EOF {
			return out, nil
		}
		if err != nil {
			return out, err
		}
		out = append(out, rec)
	}
}

// fillNextSlice advances to the next data container's next slice and
// decodes all of its records into the pending buffer. It skips
// non-data containers (the file-header container is already consumed)
// and sets done at end of stream.
func (rr *RecordReader) fillNextSlice() error {
	rr.pending = rr.pending[:0]
	rr.next = 0
	for {
		c, err := rr.rd.Next()
		if err == io.EOF {
			rr.done = true
			return nil
		}
		if err != nil {
			rr.done = true
			return err
		}
		if len(c.Blocks) == 0 || c.Blocks[0].ContentType != ContentCompressionHeader {
			continue // a non-data container; keep looking.
		}
		dc, err := ParseDataContainer(c)
		if err != nil {
			rr.done = true
			return wrapf(err, "container %d", c.Index)
		}
		for si, sl := range dc.Slices {
			if err := rr.decodeSlice(dc.Compression, sl, c.Index, si); err != nil {
				rr.done = true
				return err
			}
		}
		if len(rr.pending) > 0 {
			return nil
		}
		// A container whose slices held no records (rare) — keep reading.
	}
}

// decodeSlice decodes every record of one slice into the pending buffer.
func (rr *RecordReader) decodeSlice(h *CompressionHeader, sl *Slice, containerIdx, sliceIdx int) error {
	if sl.Header == nil {
		return errFormat("container %d slice %d has no header", containerIdx, sliceIdx)
	}
	if sl.Header.NumRecords < 0 {
		return errFormat("container %d slice %d declares a negative record count %d",
			containerIdx, sliceIdx, sl.Header.NumRecords)
	}
	src, err := sl.NewSource()
	if err != nil {
		return wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	dec, err := newRecordDecoder(h, sl.Header, src, rr.refNames, rr.readGroups)
	if err != nil {
		return wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	recs, err := dec.decodeSliceRecords(sl.Header.NumRecords)
	if err != nil {
		return wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	if dec.needsReference {
		rr.needsReference = true
	}
	rr.pending = append(rr.pending, recs...)
	rr.recordIndex += int64(len(recs))
	return nil
}

// WriteSAM decodes the whole CRAM stream and writes it as text SAM to w:
// the embedded SAM header followed by every reconstructed record. It is
// the convenience entry point the decode-to-SAM oracle exercises.
func (rr *RecordReader) WriteSAM(w io.Writer) error {
	sw := sam.NewSAMWriter(w)
	if err := sw.WriteHeader(rr.header); err != nil {
		return err
	}
	for {
		rec, err := rr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if err := sw.Write(rec); err != nil {
			return err
		}
	}
	return sw.Close()
}
