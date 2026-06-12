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
	// needsReference is set once any decoded record reached a base an
	// external reference would supply.
	needsReference bool
	// refResolver supplies external reference bases when one was set via
	// SetReference / SetRefCache; nil means decode in the C4b fallback
	// mode where reference-derived bases are filled with 'N'.
	refResolver *referenceResolver
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

// Close releases the underlying CRAM Reader's file handle, if any. If a
// reference FASTA was opened via SetReferenceFASTA, its file handle is
// released too.
func (rr *RecordReader) Close() error {
	if rr.refResolver != nil && rr.refResolver.fasta != nil {
		rr.refResolver.fasta.Close()
	}
	return rr.rd.Close()
}

// SetReference makes rr reconstruct reference-backed mapped reads from
// the supplied ReferenceSource instead of filling reference-derived
// bases with 'N'. Each slice's reference span is fetched and its MD5
// verified against the slice header; an MD5 mismatch fails the decode.
// SetReference must be called before the first Read.
//
// A *FASTAReference is recognised so its file handle is released by
// Close; any other ReferenceSource is used through its name-addressed
// Fetch method. To use the htslib REF_CACHE, call SetRefCache.
func (rr *RecordReader) SetReference(src ReferenceSource) {
	if rr.refResolver == nil {
		rr.refResolver = &referenceResolver{}
	}
	if f, ok := src.(*FASTAReference); ok {
		rr.refResolver.fasta = f
		rr.refResolver.custom = nil
		return
	}
	rr.refResolver.custom = src
}

// SetReferenceFASTA opens the named FASTA file as the decode reference
// and attaches it to rr. The FASTA's file handle is released by Close.
// It is a convenience wrapper over OpenFASTAReference + SetReference.
func (rr *RecordReader) SetReferenceFASTA(path string) error {
	f, err := OpenFASTAReference(path)
	if err != nil {
		return err
	}
	rr.SetReference(f)
	return nil
}

// SetRefCache attaches the htslib local reference cache rooted at dir
// (the REF_CACHE directory) as the decode reference, looked up by the
// MD5 each slice header records. SetRefCache and SetReference can both
// be set: an explicit FASTA is tried first, the cache second.
func (rr *RecordReader) SetRefCache(dir string) {
	if rr.refResolver == nil {
		rr.refResolver = &referenceResolver{}
	}
	rr.refResolver.cache = OpenRefCache(dir)
}

// UseRefCacheFromEnv attaches the REF_CACHE directory as a reference
// source when the REF_CACHE environment variable is set. It reports
// whether a cache was attached.
func (rr *RecordReader) UseRefCacheFromEnv() bool {
	c, ok := RefCacheFromEnv()
	if !ok {
		return false
	}
	if rr.refResolver == nil {
		rr.refResolver = &referenceResolver{}
	}
	rr.refResolver.cache = c
	return true
}

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
		if err := rr.decodeContainerInto(c, &rr.pending); err != nil {
			rr.done = true
			return err
		}
		if len(rr.pending) > 0 {
			return nil
		}
		// A container whose slices held no records (rare) — keep reading.
	}
}

// decodeContainerInto parses one structural data container and appends
// every reconstructed record of its slices to dst, in file order. It is
// the offset-agnostic core shared by the sequential fillNextSlice path
// and the seek-based RegionReader: a container at any byte offset is
// self-contained (it carries its own compression-header block), so this
// method needs only the per-file context (refNames, readGroups,
// refResolver) that RecordReader already gathered from offset 0.
func (rr *RecordReader) decodeContainerInto(c *Container, dst *[]*sam.Record) error {
	dc, err := ParseDataContainer(c)
	if err != nil {
		return wrapf(err, "container %d", c.Index)
	}
	for si, sl := range dc.Slices {
		recs, err := rr.decodeSlice(dc.Compression, sl, c.Index, si)
		if err != nil {
			return err
		}
		*dst = append(*dst, recs...)
	}
	return nil
}

// decodeSlice decodes and returns every record of one slice. The caller
// appends the returned records to its pending/result buffer; factoring
// the decode to return its records (rather than appending to a fixed
// field) lets both the sequential reader and the seek-based RegionReader
// reuse it.
func (rr *RecordReader) decodeSlice(h *CompressionHeader, sl *Slice, containerIdx, sliceIdx int) ([]*sam.Record, error) {
	if sl.Header == nil {
		return nil, errFormat("container %d slice %d has no header", containerIdx, sliceIdx)
	}
	if sl.Header.NumRecords < 0 {
		return nil, errFormat("container %d slice %d declares a negative record count %d",
			containerIdx, sliceIdx, sl.Header.NumRecords)
	}
	src, err := sl.NewSource()
	if err != nil {
		return nil, wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	refBases, refStart, err := rr.resolveSliceReference(sl)
	if err != nil {
		return nil, wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	dec, err := newRecordDecoder(h, sl.Header, src, rr.refNames, rr.readGroups, refBases, refStart)
	if err != nil {
		return nil, wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	recs, err := dec.decodeSliceRecords(sl.Header.NumRecords)
	if err != nil {
		return nil, wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	if dec.needsReference {
		rr.needsReference = true
	}
	return recs, nil
}

// resolveSliceReference resolves the reference span a slice covers. It
// returns the span bytes and the 1-based coordinate of the span's first
// base.
//
// An embedded reference (the slice's own per-span reference block,
// written by samtools' embed_ref mode) takes priority: it is
// self-contained, needs no external FASTA/REF_CACHE source, and — like
// htslib — is trusted verbatim without an MD5 cross-check. Only when no
// embedded reference is present does it consult the attached external
// sources, MD5-verifying the span.
//
// It resolves a single-reference slice (RefSeqID >= 0). An
// unmapped-reads slice (RefSeqID == -1) and a multi-reference slice
// (RefSeqID == -2) need no slice-level span — the former has no
// reference bases and the latter resolves its references per record
// against the contig table, both falling back to the C4b 'N' fill — so
// they return a nil span. A nil span with no source is the C4b path.
func (rr *RecordReader) resolveSliceReference(sl *Slice) ([]byte, int32, error) {
	sh := sl.Header
	if sh.RefSeqID < 0 {
		return nil, 0, nil
	}
	// An embedded reference is the slice's own copy of its reference
	// span. It is the most direct source and is honoured whether or not
	// an external reference is also configured.
	if sl.HasEmbeddedReference() {
		bases, err := sl.EmbeddedReference()
		if err != nil {
			return nil, 0, err
		}
		// The embedded block begins at AlignmentStart; trim any trailing
		// bytes past the slice span so refStart+span indexing matches the
		// span exactly (htslib indexes [ref_start, ref_end]).
		if sh.AlignmentSpan >= 0 && int(sh.AlignmentSpan) <= len(bases) {
			bases = bases[:sh.AlignmentSpan]
		}
		return bases, sh.AlignmentStart, nil
	}
	if !rr.refResolver.hasSource() {
		return nil, 0, nil
	}
	contig, err := rr.refNameByID(sh.RefSeqID)
	if err != nil {
		return nil, 0, err
	}
	bases, err := rr.refResolver.sliceReference(sh, contig, rr.contigMD5(sh.RefSeqID))
	if err != nil {
		return nil, 0, err
	}
	return bases, sh.AlignmentStart, nil
}

// contigMD5 returns the hex M5 tag of the @SQ entry for a reference id,
// or "" when the header carries no M5 for that reference. The M5 tag is
// the contig's whole-sequence MD5 — the key htslib's REF_CACHE uses.
func (rr *RecordReader) contigMD5(id int32) string {
	if rr.header == nil || id < 0 || int(id) >= len(rr.header.Refs) {
		return ""
	}
	for _, f := range rr.header.Refs[id].Extra {
		if f.Tag == "M5" {
			return f.Value
		}
	}
	return ""
}

// refNameByID resolves a reference id to its SAM @SQ name. It is the
// iterator-level counterpart of recordDecoder.refName.
func (rr *RecordReader) refNameByID(id int32) (string, error) {
	if id < 0 || int(id) >= len(rr.refNames) {
		return "", errFormat("reference id %d has no @SQ entry (%d known)", id, len(rr.refNames))
	}
	return rr.refNames[id], nil
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
