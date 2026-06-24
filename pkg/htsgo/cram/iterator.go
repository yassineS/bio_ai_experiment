package cram

import (
	"encoding/binary"
	"io"
	"sync"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/mdnm"
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

	// namePrefix is the basename of the opened file, used to synthesise
	// read names for records dropped by a lossy-names CRAM. It is set from
	// the underlying Reader (empty for a bare io.Reader with no path).
	namePrefix string

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

	// threads is the requested worker count for parallel CRAM decode, set
	// by SetThreads. A value below 2 (the default) keeps the reader on its
	// single-threaded path. par is the lazily-started parallel decode driver
	// when threads >= 2; it is nil on the single-threaded path. refMu guards
	// the shared referenceResolver memo and stateful reference handles while
	// the parallel workers resolve slice reference spans concurrently.
	threads int
	par     *parallelDriver
	refMu   sync.Mutex
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
	rr := &RecordReader{rd: rd, namePrefix: rd.namePrefix}
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
	rr := &RecordReader{rd: rd, namePrefix: rd.namePrefix}
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
	if rr.par != nil {
		// Unblock the feeder/workers/collector if the consumer abandoned the
		// stream before EOF so no decode goroutine is left running.
		rr.par.stop()
	}
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

// UseRefPathFromEnv attaches the network REF_PATH URL-fetch source named by the
// REF_PATH environment variable as a last-resort reference, consulted by MD5
// after an explicit FASTA and REF_CACHE both miss. It mirrors htslib's REF_PATH
// mechanism but is opt-in: it activates only when REF_PATH is set (so an
// offline decode never silently reaches out to the network). It reports whether
// a network source was attached.
func (rr *RecordReader) UseRefPathFromEnv() bool {
	p, ok := RefPathFromEnv()
	if !ok {
		return false
	}
	if rr.refResolver == nil {
		rr.refResolver = &referenceResolver{}
	}
	rr.refResolver.refpath = p
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
		if rr.threads >= 2 {
			if err := rr.fillNextSliceParallel(); err != nil {
				return nil, err
			}
		} else if err := rr.fillNextSlice(); err != nil {
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
	refBases, refStart, err := rr.resolveSliceReference(sl, h.Preservation.ReferenceRequired)
	if err != nil {
		return nil, wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	dec, err := newRecordDecoder(h, sl.Header, src, rr.refNames, rr.readGroups, refBases, refStart)
	if err != nil {
		return nil, wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	dec.namePrefix = rr.namePrefix
	recs, err := dec.decodeSliceRecords(sl.Header.NumRecords)
	if err != nil {
		return nil, wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	if dec.needsReference {
		rr.needsReference = true
	}
	// MD/NM are regenerated only against an EXTERNAL reference (a FASTA or
	// REF_CACHE attached via SetReference / SetRefCache), matching upstream
	// `samtools view -T ref file.cram`. An embedded reference is deliberately
	// excluded: upstream `samtools view file.cram` of an embed_ref CRAM emits
	// no MD/NM (it only regenerates them when given an external reference),
	// and an embed_ref=2 reduced reference does not even carry the full base
	// content the walk would need. So regeneration is suppressed for an
	// embedded-reference slice even though its bases reconstruct the SEQ.
	if !sl.HasEmbeddedReference() {
		regenerateMDNM(recs, refBases, refStart)
	}
	return recs, nil
}

// regenerateMDNM appends the reference-derived MD:Z and NM:i aux tags to
// every mapped record that lacks them, matching upstream `samtools view -T
// ref file.cram`. It is gated exactly as htslib's cram_decode_seq is (see
// reference_code/htslib/cram/cram_decode.c:1118-1119): a reference must be
// available for the slice (refBases non-nil), the record must be mapped
// (placed on a reference), and the record must not already carry the tag.
// Records with no reference span, unmapped records, and records that already
// carry the tag are left untouched.
//
// The MD and NM tags are inserted in the exact position htslib's cram_to_bam
// writes them: htslib appends MD then NM to the per-record aux block (after
// every dictionary tag), then appends the data-series RG:Z tag *after* that
// aux block (reference_code/htslib/cram/cram_decode.c:3183-3197). Our
// reconstruction appends that data-series RG as the final aux tag, so when
// the last existing tag is RG the new MD/NM are spliced in just before it;
// otherwise (no trailing data-series RG — including the case where RG was
// emitted in its tag-dictionary position) they go at the very end. Within the
// pair MD always precedes NM.
//
// refBases is the slice's reference span and refStart is the 1-based
// reference coordinate of refBases[0], so the 0-based offset passed to
// mdnm.Compute is refStart-1. A nil refBases (no external/embedded reference,
// or an unmapped/multi-reference slice) means "no reference available" and
// suppresses regeneration entirely.
func regenerateMDNM(recs []*sam.Record, refBases []byte, refStart int32) {
	if refBases == nil {
		return
	}
	refOffset := int(refStart) - 1
	for _, rec := range recs {
		if rec.IsUnmapped() || rec.RName == "" || rec.RName == "*" || len(rec.Cigar) == 0 {
			continue
		}
		_, hasMD := rec.GetAux("MD")
		_, hasNM := rec.GetAux("NM")
		if hasMD && hasNM {
			continue
		}
		md, nm := mdnm.Compute(rec, refBases, refOffset)
		var add []sam.Aux
		if !hasMD {
			add = append(add, sam.Aux{Tag: "MD", Type: 'Z', Value: md})
		}
		if !hasNM {
			add = append(add, sam.Aux{Tag: "NM", Type: 'i', Value: int64(nm)})
		}
		rec.Aux = insertBeforeTrailingRG(rec.Aux, add)
		rec.InvalidateAuxIndex()
	}
}

// insertBeforeTrailingRG returns aux with add spliced in immediately before a
// trailing RG tag, or appended at the end when the last tag is not RG. This
// mirrors htslib writing the data-series RG:Z tag after the rest of the aux
// block: a record whose final tag is the data-series RG gets MD/NM placed
// just ahead of it, exactly as `samtools view -T ref` emits them.
func insertBeforeTrailingRG(aux, add []sam.Aux) []sam.Aux {
	if len(add) == 0 {
		return aux
	}
	if n := len(aux); n > 0 && aux[n-1].Tag == "RG" {
		out := make([]sam.Aux, 0, n+len(add))
		out = append(out, aux[:n-1]...)
		out = append(out, add...)
		out = append(out, aux[n-1])
		return out
	}
	return append(aux, add...)
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
//
// refRequired is the container compression header's RR (reference-required)
// preservation-map entry. When it is false the records were encoded
// reference-free — their bases are carried verbatim, no implicit match runs
// need the reference — so the external FASTA/REF_CACHE is deliberately NOT
// consulted, mirroring htslib's cram_decode.c, which gates the whole
// cram_get_ref load on !comp_hdr->no_ref. This is what lets a CRAM whose
// contig is absent from the -T reference (so it was encoded reference-free)
// decode without a "contig not in index" error: an embedded reference, when
// present, is still honoured because it is self-contained.
func (rr *RecordReader) resolveSliceReference(sl *Slice, refRequired bool) ([]byte, int32, error) {
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
	// The records are reference-free (RR=0): their bases are stored verbatim,
	// so the external reference is not consulted — and must not be, since the
	// contig may be absent from it (a CRAM encoded with a -T reference whose
	// contigs do not match). htslib gates the equivalent cram_get_ref load on
	// !no_ref for exactly this case.
	if !refRequired {
		return nil, 0, nil
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
