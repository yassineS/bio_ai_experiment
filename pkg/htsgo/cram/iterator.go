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
	// Only ONE slice's records are held at a time — refilling overwrites the
	// buffer (pending[:0]) so the prior slice's records are released, keeping
	// the live working set to a single slice rather than a whole container.
	pending []*sam.Record
	// next is the index of the next pending record to return.
	next int
	// done is set once the stream is exhausted.
	done bool

	// curDC is the data container currently being streamed slice-by-slice;
	// nil before the first refill and once the stream is exhausted.
	curDC *DataContainer
	// curDCIndex is curDC's structural container index, kept for error context.
	curDCIndex int
	// curSliceIdx is the index of the next slice to decode within curDC.
	curSliceIdx int
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

// fillNextSlice decodes exactly ONE slice into the pending buffer, advancing
// through the current data container's slices across successive calls and
// moving on to the next container when the current one is exhausted. Decoding a
// single slice per refill — rather than a whole container — keeps the live
// working set bounded to one slice of fat sam.Records: the next refill
// overwrites pending (pending[:0]), releasing the prior slice's records.
//
// It skips non-data containers (the file-header container is already consumed)
// and sets done at end of stream. Record order is byte-identical to the
// whole-container path: containers in file order, slices in index order
// (0,1,2,…), records within a slice in order — exactly as decodeContainerInto
// produced. Empty slices (a slice that decodes to zero records, rare) are
// skipped transparently: the loop simply advances to the next slice/container
// until it has records or reaches EOF.
func (rr *RecordReader) fillNextSlice() error {
	rr.pending = rr.pending[:0]
	rr.next = 0
	for {
		// If a container is mid-stream, emit its next slice.
		if rr.curDC != nil && rr.curSliceIdx < len(rr.curDC.Slices) {
			si := rr.curSliceIdx
			recs, err := rr.decodeSlice(rr.curDC.Compression, rr.curDC.Slices[si], rr.curDCIndex, si)
			if err != nil {
				rr.done = true
				return err
			}
			rr.curSliceIdx++
			rr.pending = append(rr.pending, recs...)
			if len(rr.pending) > 0 {
				return nil
			}
			continue // an empty slice — advance to the next one.
		}
		// The current container is exhausted (or there is none): read the next.
		c, err := rr.rd.Next()
		if err == io.EOF {
			rr.done = true
			rr.curDC = nil
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
		rr.curDC = dc
		rr.curDCIndex = c.Index
		rr.curSliceIdx = 0
		// Loop again: the first branch now decodes slice 0.
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
	recs, suppress, err := dec.decodeSliceRecords(sl.Header.NumRecords)
	if err != nil {
		return nil, wrapf(err, "container %d slice %d", containerIdx, sliceIdx)
	}
	if dec.needsReference {
		rr.needsReference = true
	}
	// MD/NM are regenerated whenever the slice HAS a reference — external (a
	// FASTA or REF_CACHE attached via SetReference / SetRefCache) OR embedded
	// (the reference bytes stored in the slice itself) — and skipped only for a
	// reference-free (no_ref / RR=0) slice. This mirrors htslib exactly, which
	// gates regeneration on `s->ref != NULL` (cram_decode.c:1118): its embedded
	// branch sets `s->ref` for ANY embedded block (a real-reference embed OR a
	// consensus / embed_ref=1 embed), so upstream regenerates MD/NM from the
	// embedded bases in BOTH cases — only a no_ref slice is left bare.
	// resolveSliceReference returns the embedded bases for an embedded slice and
	// nil for a no_ref slice, and regenerateMDNM is a no-op when refBases is nil,
	// so calling it unconditionally reproduces htslib's gate precisely.
	//
	// suppress carries the per-record "cF" (CRAM-flags) bits an embed_ref=2
	// writer stores: when the source record had no MD (bit 1) or no NM (bit 2),
	// htslib leaves it bare even though the embedded reference is available
	// (cram_decode.c:2117-2122). regenerateMDNM honours those bits so a
	// reduced-reference CRAM whose reads carried no MD/NM stays byte-identical
	// to `samtools view`.
	regenerateMDNM(recs, suppress, refBases, refStart)
	return recs, nil
}

// mdnmSuppress carries the per-record MD/NM regeneration suppression an
// embed_ref=2 (reduced/consensus reference) CRAM encodes in its "cF" aux
// tag: md is set when the source record had no MD tag and nm when it had
// no NM (cram_encode.c:3417). regenerateMDNM skips the corresponding tag
// for such a record, reproducing htslib's per-record gate exactly.
type mdnmSuppress struct {
	md bool
	nm bool
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
//
// suppress, when non-nil, is indexed parallel to recs and carries the
// per-record "cF" suppression an embed_ref=2 (reduced/consensus) CRAM stores:
// suppress[i].md / .nm true means the source record had no MD / no NM tag, so
// that record is left bare even though the reference is available. htslib does
// the same by forcing has_MD/has_NM at cram_decode.c:2117-2122. A nil suppress
// (the common case, no cF bits) regenerates every eligible record.
func regenerateMDNM(recs []*sam.Record, suppress []mdnmSuppress, refBases []byte, refStart int32) {
	if refBases == nil {
		return
	}
	refOffset := int(refStart) - 1
	for i, rec := range recs {
		if rec.IsUnmapped() || rec.RName == "" || rec.RName == "*" || len(rec.Cigar) == 0 {
			continue
		}
		var sup mdnmSuppress
		if i < len(suppress) {
			sup = suppress[i]
		}
		// Direct linear scan over the (small, <=~12-tag) aux list rather than
		// rec.GetAux, which builds and discards an aux-index map per record;
		// this path runs for every mapped record of every slice (~10k/slice),
		// so the map churn is pure GC pressure under the decode memory cap.
		hasMD, hasNM := false, false
		for j := range rec.Aux {
			switch rec.Aux[j].Tag {
			case "MD":
				hasMD = true
			case "NM":
				hasNM = true
			}
		}
		// A cF "no MD"/"no NM" bit is equivalent to the tag already being
		// present for the regeneration gate: htslib forces has_MD/has_NM so
		// the `!has_MD`/`!has_NM` condition is false and nothing is emitted.
		wantMD := !hasMD && !sup.md
		wantNM := !hasNM && !sup.nm
		if !wantMD && !wantNM {
			continue
		}
		md, nm := mdnm.Compute(rec, refBases, refOffset)
		var add []sam.Aux
		if wantMD {
			add = append(add, sam.Aux{Tag: "MD", Type: 'Z', Value: md})
		}
		if wantNM {
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
	n := len(aux)
	// Splice point: just before a trailing data-series RG, else at the end.
	insertPos := n
	if n > 0 && aux[n-1].Tag == "RG" {
		insertPos = n - 1
	}
	need := n + len(add)
	// Reuse the existing backing array's spare capacity when present
	// (decodeTags over-allocates and mergeAux rarely fills it), splicing in
	// place with a single memmove instead of allocating a fresh aux slice per
	// record — the standalone reallocation here was the dominant resident
	// allocation of CRAM decode. `copy` is memmove-safe for overlapping
	// ranges, so shifting the [insertPos:n] tail right is correct. The emitted
	// elements and their order are byte-identical (… MD NM [RG]).
	if cap(aux) >= need {
		aux = aux[:need]
		copy(aux[insertPos+len(add):], aux[insertPos:n])
		copy(aux[insertPos:], add)
		return aux
	}
	// Insufficient capacity: one exact-sized allocation (still cheaper than the
	// previous make-plus-three-appends, and only when the slice can't grow).
	out := make([]sam.Aux, need)
	copy(out, aux[:insertPos])
	copy(out[insertPos:], add)
	copy(out[insertPos+len(add):], aux[insertPos:n])
	return out
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
	// !no_ref, and likewise does NOT regenerate MD/NM for a no_ref slice even
	// when -T is given (decode_md = s->ref && ..., and s->ref stays NULL here).
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
