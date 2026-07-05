package cram

import (
	"fmt"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// errFormat is the family of errors a CRAM record decode returns for a
// malformed or unsupported file. A decode error is always one of these,
// never a panic; the fuzz target depends on that guarantee.
var errFormat = func(format string, args ...interface{}) error {
	return fmt.Errorf("cram: "+format, args...)
}

// wrapf wraps err with a fmt.Errorf-style context prefix.
func wrapf(err error, format string, args ...interface{}) error {
	return fmt.Errorf("cram: "+fmt.Sprintf(format, args...)+": %w", err)
}

// CRAM per-record "CRAM flags" (CF) bits. The CF data series carries
// these alongside the SAM bit flags; they describe how the record was
// encoded rather than its alignment properties.
const (
	// cfQualityPreserved is set when the QS data series carries this
	// record's quality scores.
	cfQualityPreserved = 0x1
	// cfDetached is set when the record stores its own mate information
	// (NS, NP, TS) rather than pointing at a downstream mate.
	cfDetached = 0x2
	// cfHasMateDownstream is set when the record's mate is a later record
	// in the same slice and the NF series holds the records-to-skip count.
	cfHasMateDownstream = 0x4
	// cfNoSeq is set when the record carried no SEQ ("*"). The read
	// features and RL still describe the CIGAR (so a mapped no-SEQ record
	// round-trips its alignment), but the reconstructed bases and
	// qualities are discarded on decode: SEQ and QUAL render as "*". This
	// matches htslib's CRAM_FLAG_NO_SEQ (cram_structs.h: 1<<3), which the
	// encoder sets for any record with zero query bases and the decoder
	// honours by resetting the read length to zero.
	cfNoSeq = 0x8
)

// CRAM per-record "mate flags" (MF) bits, used only by a detached
// record to reconstruct the mate-related SAM flag bits.
const (
	mfMateReverse  = 0x1 // the mate is on the reverse strand.
	mfMateUnmapped = 0x2 // the mate is unmapped.
)

// recordDecoder threads the per-slice decode state through a slice's
// records: the compression header that names every encoding, the
// SeriesSource that holds the decompressed blocks and cursors, and the
// reference metadata a record needs to resolve its coordinates.
type recordDecoder struct {
	h          *CompressionHeader
	src        *SeriesSource
	slice      *SliceHeader
	refNames   []string // SAM @SQ names indexed by reference id.
	readGroups []string // SAM @RG IDs indexed by the CRAM RG series value.

	// prevAlignmentStart carries the running alignment start for the
	// AP-delta encoding: when the preservation map's AP entry is true,
	// each record's AP series value is a delta from the previous
	// record's start.
	prevAlignmentStart int32

	// tagDict is the parsed tag-combination dictionary: tagDict[i] is the
	// ordered list of three-byte tag keys a record with TL value i carries.
	tagDict [][]tagKey

	// needsReference is set when reconstructing a record reached a base
	// that an external reference would supply. A reference-free CRAM file
	// never sets it; a reference-backed file does, and the affected bases
	// are filled with 'N'. The iterator surfaces it via NeedsReference.
	needsReference bool

	// readLenLimit bounds a single record's read length. It is derived
	// from the slice's total decompressed block size so a malformed file
	// declaring an astronomically long read cannot trigger a huge
	// allocation; a genuine read is always far smaller than the data
	// blocks that encode it.
	readLenLimit int32

	// refBases holds the slice's reference span, upper-cased, indexed so
	// that refBases[i] is reference position SliceHeader.AlignmentStart+i.
	// It is nil when no external reference was supplied, in which case
	// reference-derived bases are filled with 'N' (the C4b fallback).
	refBases []byte
	// refStart is the 1-based reference coordinate of refBases[0].
	refStart int32
	// substMatrix decodes a substitution feature's BS code into a read
	// base relative to the reference base at that position.
	substMatrix substMatrix

	// qsEnc and baEnc are the resolved QS (quality score) and BA (base) data
	// series encodings, looked up once from the compression header when the
	// decoder is built rather than re-resolved through the data-series map on
	// every quality base and every unmapped base. They are the per-slice cache
	// the per-base hot loops (readQualityBlockInto, decodeUnmapped) read; both
	// are nil when the series is absent, in which case those loops fall back to
	// the unchanged per-value byteSeries path. They point into the container's
	// shared CompressionHeader and are only ever read after construction, so a
	// worker decoding its own slice never mutates them.
	qsEnc *Encoding
	baEnc *Encoding

	// namePrefix is the basename of the opened file. When a lossy-names
	// (read-names-not-preserved) CRAM drops a record's name, the decoder
	// synthesises "<namePrefix>:<record_number>", matching htslib's
	// cram_to_bam name generation. Empty when the file was opened from a
	// bare io.Reader with no path (htslib likewise has no prefix then).
	namePrefix string

	// featScratch is a reused backing array for decodeFeatures: each mapped
	// record's feature list is decoded into it and fully consumed by
	// reconstructMapped before the next record's decode reuses it. This
	// removes the dominant per-record allocation of CRAM->BAM decode (a
	// fresh []readFeature per mapped record). It is safe because a single
	// recordDecoder decodes its slice's records sequentially and the
	// returned slice is never retained past the same record's reconstruction
	// (the emitted sam.Record copies SEQ/QUAL/CIGAR out, never the features).
	featScratch []readFeature

	// rawAuxBAMSink mirrors the RecordReader flag of the same name: when set,
	// decodeTags builds each record's aux block directly as raw on-disk BAM aux
	// bytes (rec.RawAux) and leaves rec.Aux nil, the memory-lean CRAM→BAM view
	// passthrough. Default false keeps the eager []sam.Aux path unchanged.
	rawAuxBAMSink bool

	// arena is the per-slice record-field allocation arena, used only on the
	// rawAuxBAMSink (view -b -T) path. It collapses the per-record make() of a
	// record's packed SEQ, QUAL, CIGAR and raw aux into a handful of growable
	// backing buffers shared by every record of the slice: each record's fields
	// are 3-index sub-slices of those buffers, so the ~N per-record allocations
	// of a slice become ~O(log N) amortised arena growths. It is nil on the eager
	// path, where the per-record make() is kept verbatim.
	//
	// SAFETY: on the single-threaded path the arena is the RecordReader's
	// persistent one, reused across slices (reset, not reallocated, per slice);
	// on the parallel path each worker uses its own per-slice arena. Either way
	// the arena is reset at the start of a slice's decode, and a slice's records
	// are fully served+written before the next reset (single-threaded) or the
	// worker's slice is collected before its arena is reused — so no record's
	// reused field backing is read after a reset. The 3-index sub-slices cap each
	// field at its exact length, so a later append on any field (the MD/NM raw
	// splice, or a downstream consumer) reallocates rather than writing into a
	// neighbouring record's arena region.
	arena *recordArena
}

// recordArena holds the growable backing buffers the rawAuxBAMSink decode path
// sub-slices for each record's fields, replacing the per-record make(). Its
// buffers are rewound to empty (reset) at the start of each slice and the SAME
// backing is reused for the next slice, so steady-state decode allocates almost
// nothing — the lever that lowers the Go runtime's retained pages. On the
// single-threaded path one arena (on the RecordReader) serves the whole decode;
// on the parallel path each worker has its own (a per-worker arena, reset per
// slice the worker decodes).
type recordArena struct {
	qual  []byte        // packed run of every record's QUAL scores.
	seq   []byte        // packed run of every record's 4-bit SEQ nibble block.
	cigar []sam.CigarOp // run of every record's CIGAR ops.
	aux   []byte        // run of every record's raw on-disk BAM aux block.
	// recSlab and drSlab are per-slice slabs the decoder hands records out of,
	// so a slice's nRecords sam.Record / decodedRecord structs are one slab
	// allocation each instead of one make() per record. recIdx is the next free
	// slot. They are pre-sized to nRecords in decodeSliceRecords; nextRecord
	// grows them (rare — the count is known) if a decode overruns the estimate.
	recSlab []sam.Record
	drSlab  []decodedRecord
	recIdx  int
	// seqScratch and coveredScratch are reused (per record, not per slice)
	// reconstruction temporaries that never escape into a record: seqScratch is
	// the ASCII SEQ a mapped record is reconstructed into before being packed
	// into the seq arena, and coveredScratch tracks which read positions a
	// feature has filled. Reusing them removes the two remaining per-record
	// temporaries of the mapped reconstruction.
	seqScratch     []byte
	coveredScratch []bool
	// auxScratch is the reused []sam.Aux header slice decodeTags assembles a
	// record's tags into on the gated path. It never escapes: buildRawAux copies
	// every field's bytes into the aux byte arena before the next record's decode
	// reuses it, so a single reused header slice per decoder suffices in place of
	// the per-record make([]sam.Aux). It is the header backing only; the tag
	// value bytes are serialised out by buildRawAux.
	auxScratch []sam.Aux
}

// reset rewinds every per-slice run to empty so the SAME backing buffers are
// reused for the next slice instead of reallocating fresh ones. It is the
// cross-slice reuse the steady-state allocation rate (and thus the runtime's
// retained pages) depends on. It is called at the start of each slice's decode
// on the single-threaded path, where a slice's records are fully served and
// written before the next slice's decode begins — so no record's reused backing
// is read after this rewind (see fillNextSlice's safety note).
func (a *recordArena) reset() {
	a.qual = a.qual[:0]
	a.seq = a.seq[:0]
	a.cigar = a.cigar[:0]
	a.aux = a.aux[:0]
	a.recIdx = 0
}

// reserve sizes the per-slice record slabs to at least n (the slice's
// authoritative record count) and clears the slots that will be handed out, so
// nextRecord returns zero-valued structs from a reused backing array. The slab
// is grown (and re-pointed) only when n exceeds its capacity; within a slice it
// never reallocates, so the &slab[i] pointers handed out stay valid for that
// slice. Records of a PRIOR slice that aliased the same slab must already be
// dead (served+written) — guaranteed on the single-threaded view path.
func (a *recordArena) reserve(n int) {
	if cap(a.recSlab) < n {
		a.recSlab = make([]sam.Record, n)
		a.drSlab = make([]decodedRecord, n)
	} else {
		a.recSlab = a.recSlab[:n]
		a.drSlab = a.drSlab[:n]
		for i := range a.recSlab {
			a.recSlab[i] = sam.Record{}
			a.drSlab[i] = decodedRecord{}
		}
	}
	a.recIdx = 0
}

// nextRecord hands out the next zero-valued (sam.Record, decodedRecord) pair
// from the per-slice slabs. The slabs are pre-sized to the slice's record count
// by reserve, so this never reallocates them and the returned pointers stay
// valid for the slab's lifetime — which equals the records' (a slice's records
// reference the same slab), so retaining records past the slice boundary stays
// correct. A decode that somehow runs past the reserved count (it cannot: the
// loop runs exactly nRecords times) falls back to a standalone allocation rather
// than growing the slab, keeping the already-handed-out pointers valid.
func (a *recordArena) nextRecord() (*sam.Record, *decodedRecord) {
	if a.recIdx >= len(a.recSlab) {
		return &sam.Record{}, &decodedRecord{}
	}
	i := a.recIdx
	a.recIdx++
	return &a.recSlab[i], &a.drSlab[i]
}

// qualFor returns a length-n QUAL sub-slice backed by the arena, prefilled with
// the 0xff "no quality" sentinel. The returned slice is capped at exactly n
// (a 3-index slice) so a later in-place edit (reverseQual) stays within the
// record's own region and any append reallocates instead of corrupting a
// neighbouring record's QUAL.
func (a *recordArena) qualFor(n int) []byte {
	start := len(a.qual)
	a.qual = growBytes(a.qual, n)
	q := a.qual[start : start+n : start+n]
	for i := range q {
		q[i] = 0xff
	}
	return q
}

// packSeq packs n ASCII bases from src into the arena's SEQ run and returns the
// packed nibble block (capped at its exact length so a later append reallocates).
func (a *recordArena) packSeq(src []byte, n int) []byte {
	packed := (n + 1) / 2
	start := len(a.seq)
	a.seq = growBytes(a.seq, packed)
	dst := a.seq[start : start+packed : start+packed]
	return sam.PackSeqBytesInto(dst[:0], src, n)
}

// growBytes extends b by n bytes and returns it, reusing b's spare capacity when
// present and otherwise growing the backing (the standard append doubling). It
// avoids the temporary make() that `append(b, make([]byte, n)...)` allocates per
// call, which would defeat the arena's purpose of near-zero steady-state
// allocation. The n new bytes are not zeroed; callers overwrite them.
func growBytes(b []byte, n int) []byte {
	need := len(b) + n
	if cap(b) >= need {
		return b[:need]
	}
	nb := make([]byte, need, need*2)
	copy(nb, b)
	return nb
}

// reconScratch returns reusable length-n ASCII-SEQ and coverage scratch buffers
// for one record's reconstruction. They are overwritten by the next record and
// never escape into a sam.Record (storeSeq packs SEQ into the SEQ arena), so a
// single backing array per decoder suffices. The coverage buffer is zeroed
// before return so each record starts from a clean "nothing covered" state, as a
// fresh make([]bool) would.
func (a *recordArena) reconScratch(n int) (seq []byte, covered []bool) {
	if cap(a.seqScratch) < n {
		a.seqScratch = make([]byte, n)
	}
	seq = a.seqScratch[:n]
	if cap(a.coveredScratch) < n {
		a.coveredScratch = make([]bool, n)
	}
	covered = a.coveredScratch[:n]
	for i := range covered {
		covered[i] = false
	}
	return seq, covered
}

// newRecordDecoder builds a recordDecoder for one slice. It parses the
// tag dictionary up front so per-record tag reconstruction is a lookup,
// and builds the base-substitution matrix from the preservation map.
// refBases, when non-nil, is the slice's resolved and MD5-verified
// reference span; refStart is the 1-based coordinate of its first base.
func newRecordDecoder(h *CompressionHeader, sh *SliceHeader, src *SeriesSource, refNames, readGroups []string, refBases []byte, refStart int32) (*recordDecoder, error) {
	td, err := parseTagDictionary(h.Preservation.TagDictionary)
	if err != nil {
		return nil, err
	}
	return &recordDecoder{
		h:                  h,
		src:                src,
		slice:              sh,
		refNames:           refNames,
		readGroups:         readGroups,
		prevAlignmentStart: sh.AlignmentStart,
		tagDict:            td,
		readLenLimit:       src.s.totalBytes(),
		refBases:           refBases,
		refStart:           refStart,
		substMatrix:        newSubstMatrix(h.Preservation.SubstitutionMatrix),
		qsEnc:              h.Encoding("QS"),
		baEnc:              h.Encoding("BA"),
	}, nil
}

// parseTagDictionary parses the preservation map's TD entry: a run of
// NUL-separated lists, each list being a concatenation of three-byte tag
// keys. The list at index i is the set of tags a record whose TL series
// value is i carries.
func parseTagDictionary(td []byte) ([][]tagKey, error) {
	if len(td) == 0 {
		return nil, nil
	}
	var out [][]tagKey
	var cur []tagKey
	i := 0
	for i < len(td) {
		if td[i] == 0 {
			out = append(out, cur)
			cur = nil
			i++
			continue
		}
		if i+3 > len(td) {
			return nil, errFormat("tag dictionary entry %d: truncated three-byte tag key", len(out))
		}
		cur = append(cur, tagKey{td[i], td[i+1], td[i+2]})
		i += 3
	}
	// A dictionary that does not end on a NUL still yields its final list.
	if cur != nil {
		out = append(out, cur)
	}
	return out, nil
}

// intSeries decodes one integer value of the named two-character data
// series from the slice's source.
func (rd *recordDecoder) intSeries(key string) (int32, error) {
	enc := rd.h.Encoding(key)
	if enc == nil {
		return 0, errFormat("data series %q is required but absent from the encoding map", key)
	}
	v, err := enc.decodeInt(rd.src.s)
	if err != nil {
		return 0, wrapf(err, "data series %q", key)
	}
	return v, nil
}

// byteSeries decodes one single-byte value of the named data series.
func (rd *recordDecoder) byteSeries(key string) (byte, error) {
	enc := rd.h.Encoding(key)
	if enc == nil {
		return 0, errFormat("data series %q is required but absent from the encoding map", key)
	}
	switch enc.ID {
	case EncodingExternal:
		c, err := rd.src.s.cursor(enc.ExternalID)
		if err != nil {
			return 0, wrapf(err, "data series %q", key)
		}
		b, err := c.readByte()
		if err != nil {
			return 0, wrapf(err, "data series %q", key)
		}
		return b, nil
	case EncodingHuffman, EncodingBeta:
		v, err := enc.decodeInt(rd.src.s)
		if err != nil {
			return 0, wrapf(err, "data series %q", key)
		}
		return byte(v), nil
	case EncodingConstByte, EncodingConstInt:
		// CRAM v4: a constant byte-valued series (no block read).
		return byte(enc.ConstValue), nil
	default:
		return 0, errFormat("data series %q has %s encoding, which is not byte-valued", key, enc.ID)
	}
}

// byteArraySeries decodes one variable-length byte-array value of the
// named data series.
func (rd *recordDecoder) byteArraySeries(key string) ([]byte, error) {
	enc := rd.h.Encoding(key)
	if enc == nil {
		return nil, errFormat("data series %q is required but absent from the encoding map", key)
	}
	b, err := enc.decodeByteArray(rd.src.s)
	if err != nil {
		return nil, wrapf(err, "data series %q", key)
	}
	return b, nil
}

// byteArraySeriesBorrow is byteArraySeries returning a NO-COPY sub-slice of
// the block buffer (see Encoding.decodeByteArrayBorrow). The returned bytes
// alias the decompressed block and MUST be copied out (e.g. string()) before
// the slice's source is released; use it only for immediately-consumed values
// such as a read name that is copied into rec.QName.
func (rd *recordDecoder) byteArraySeriesBorrow(key string) ([]byte, error) {
	enc := rd.h.Encoding(key)
	if enc == nil {
		return nil, errFormat("data series %q is required but absent from the encoding map", key)
	}
	b, err := enc.decodeByteArrayBorrow(rd.src.s)
	if err != nil {
		return nil, wrapf(err, "data series %q", key)
	}
	return b, nil
}

// optIntSeries decodes one integer value of the named series, or returns
// zero (and ok=false) when the series is absent from the encoding map.
func (rd *recordDecoder) optIntSeries(key string) (int32, bool, error) {
	enc := rd.h.Encoding(key)
	if enc == nil {
		return 0, false, nil
	}
	v, err := enc.decodeInt(rd.src.s)
	if err != nil {
		return 0, false, wrapf(err, "data series %q", key)
	}
	return v, true, nil
}

// refName resolves a reference id to its SAM name. -1 is the unmapped
// sentinel and maps to "" (the SAM "*"). An out-of-range id is an error.
func (rd *recordDecoder) refName(id int32) (string, error) {
	if id < 0 {
		return "", nil
	}
	if int(id) >= len(rd.refNames) {
		return "", errFormat("reference id %d has no @SQ entry (%d known)", id, len(rd.refNames))
	}
	return rd.refNames[id], nil
}

// decodedRecord pairs a reconstructed sam.Record with the per-record
// CRAM bookkeeping a slice's second pass needs: the next-fragment
// distance to a downstream mate. It is the unit decodeRecord produces;
// resolveMates then walks a slice's decodedRecords to fill in the
// mate-related SAM fields for non-detached pairs.
type decodedRecord struct {
	rec *sam.Record
	// mateDownstream is true when the CRAM flags marked this record as
	// having its mate stored later in the same slice.
	mateDownstream bool
	// nextFragment is the records-to-skip distance to the downstream
	// mate, valid only when mateDownstream is true.
	nextFragment int32
	// nameStored is true when an explicit read name was decoded for this
	// record (from the RN series, whether at the normal read-name position
	// or, for a lossy-names detached record, inside the mate block). When
	// it is false the name was dropped by a lossy-names CRAM and must be
	// reconstructed after the slice is decoded — copied from a mate if one
	// carries a name, otherwise synthesised as "<prefix>:<record_number>",
	// mirroring htslib's cram_to_bam name generation.
	nameStored bool
	// mateIndex is the in-slice index of this record's mate when the pair
	// is stored within the slice (set by resolveMates for both directions
	// of a downstream-mate pair), or -1 when unknown. It is used to pick
	// the record number htslib's cram_to_bam uses when synthesising a
	// dropped name: the mate's index when the mate is an earlier record.
	mateIndex int

	// suppressMD and suppressNM mark that MD / NM must NOT be regenerated
	// for this record even though a reference is available, reproducing
	// htslib's per-record regeneration gate. They are computed by
	// mdnmSuppressFor from the version-specific signal: for CRAM < 4 the
	// per-record "cF" aux-tag bits an embed_ref=2 writer stores
	// (cram_encode.c:3417, consumed at cram_decode.c:2117-2122); for CRAM
	// >= 4 the absence of an MD*/NM* placeholder in the record's tag
	// dictionary (regeneration is requested only when the placeholder is
	// present, has_MD < 0 at cram_decode.c:2046).
	suppressMD bool
	suppressNM bool
}

// decodeRecord decodes one CRAM alignment record into a decodedRecord.
// The traversal follows the CRAM v3.0 record layout: the bit flags (BF)
// and CRAM flags (CF), the reference id for a multi-reference slice, the
// read length (RL), the alignment position (AP), the read group (RG),
// the read name, the mate information for a detached record, the tag
// values, and finally either the read-feature list (mapped) or the raw
// bases (unmapped).
func (rd *recordDecoder) decodeRecord(index int) (*decodedRecord, error) {
	var rec *sam.Record
	var dr *decodedRecord
	if rd.arena != nil {
		// Hand the record/decodedRecord pair out of the per-slice slabs instead
		// of a per-record make(); both come back zero-valued.
		rec, dr = rd.arena.nextRecord()
		dr.rec = rec
		dr.mateIndex = -1
	} else {
		rec = &sam.Record{}
		dr = &decodedRecord{rec: rec, mateIndex: -1}
	}

	bf, err := rd.intSeries("BF")
	if err != nil {
		return nil, wrapf(err, "record %d bit flags", index)
	}
	rec.Flag = uint16(bf)

	cf, err := rd.intSeries("CF")
	if err != nil {
		return nil, wrapf(err, "record %d CRAM flags", index)
	}

	// The reference id is per-record only for a multi-reference slice
	// (ref id -2); a single-reference slice reuses the slice header's id.
	refID := rd.slice.RefSeqID
	if rd.slice.RefSeqID == -2 {
		v, ok, ierr := rd.optIntSeries("RI")
		if ierr != nil {
			return nil, wrapf(ierr, "record %d reference id", index)
		}
		if !ok {
			return nil, errFormat("record %d: multi-reference slice needs an RI data series", index)
		}
		refID = v
	}
	rec.RName, err = rd.refName(refID)
	if err != nil {
		return nil, wrapf(err, "record %d", index)
	}

	readLen, err := rd.intSeries("RL")
	if err != nil {
		return nil, wrapf(err, "record %d read length", index)
	}
	if readLen < 0 {
		return nil, errFormat("record %d declares a negative read length %d", index, readLen)
	}
	// A read longer than every data block of the slice combined cannot
	// be real: the bases, qualities and features that make up the read
	// are all encoded in those blocks. Rejecting it here keeps a
	// malformed RL value from triggering a huge SEQ/QUAL allocation.
	if readLen > rd.readLenLimit {
		return nil, errFormat("record %d declares read length %d, larger than the slice's %d data bytes",
			index, readLen, rd.readLenLimit)
	}

	ap, err := rd.intSeries("AP")
	if err != nil {
		return nil, wrapf(err, "record %d alignment position", index)
	}
	if rd.h.Preservation.APDelta {
		rd.prevAlignmentStart += ap
		rec.Pos = int64(rd.prevAlignmentStart)
	} else {
		rec.Pos = int64(ap)
	}

	rgValue, err := rd.intSeries("RG")
	if err != nil {
		return nil, wrapf(err, "record %d read group", index)
	}

	if err := rd.decodeReadName(dr, cf, index); err != nil {
		return nil, err
	}
	if err := rd.decodeMate(dr, cf, index); err != nil {
		return nil, err
	}

	tl, err := rd.intSeries("TL")
	if err != nil {
		return nil, wrapf(err, "record %d tag line", index)
	}
	tags, rgEmitted, suppress, err := rd.decodeTags(tl, rgValue, index)
	if err != nil {
		return nil, err
	}
	dr.suppressMD = suppress.md
	dr.suppressNM = suppress.nm

	// The reconstructed sequence/quality/CIGAR depends on whether the
	// record is mapped: a mapped record stores read features, an
	// unmapped one stores its bases verbatim.
	if rec.Flag&sam.FlagUnmapped == 0 {
		mq, merr := rd.intSeries("MQ")
		if merr != nil {
			return nil, wrapf(merr, "record %d mapping quality", index)
		}
		rec.MapQ = uint8(mq)
		if err := rd.decodeMapped(rec, cf, readLen, index); err != nil {
			return nil, err
		}
	} else {
		if err := rd.decodeUnmapped(rec, cf, readLen, index); err != nil {
			return nil, err
		}
	}

	// When the QO preservation flag is off (quality scores stored in
	// reference/original orientation rather than read orientation), a
	// reverse-strand read's quality string was written reversed relative
	// to the emitted SEQ, so flip it back. This mirrors htslib's
	// cram_decode.c `!qs_seq_orient` handling. CRAM v3.1+ and v4 writers
	// commonly set QO=0; older files default qs_seq_orient on (QO true),
	// which makes this a no-op.
	if !rd.h.Preservation.QualityScoreSeqOrient && rec.Flag&sam.FlagReverse != 0 {
		reverseQual(rec.Qual)
	}

	// A no-SEQ record (CRAM_FLAG_NO_SEQ) carried no query bases: its RL and
	// read features exist only to describe the CIGAR, so the reconstructed
	// SEQ (reference- or 'N'-filled) and QUAL are discarded and the record
	// renders SEQ/QUAL as "*". htslib resets cr->len to 0 here for the same
	// effect; the CIGAR computed above is kept intact.
	if cf&cfNoSeq != 0 {
		// Discard any packed SEQ the passthrough stored so both paths converge
		// on the identical "*" record before the BAM writer sees it (eager sets
		// rec.Seq = "*", which the writer encodes as l_seq=1, packed 0xf0).
		rec.RawSeq = nil
		rec.SeqLen = 0
		rec.Seq = "*"
		rec.Qual = nil
	}

	// The auxiliary tags are assembled from the dictionary-stored tags
	// and the read-group data series. For CRAM v2/v3 the RG tag is not in
	// the dictionary, so it is spliced in here at the position samtools
	// emits it (immediately before the program tag, or last when none).
	// For CRAM v4 the dictionary already carries an "RG*" placeholder that
	// decodeTags expanded in place, so rgEmitted is set and no second copy
	// is added.
	var finalAux []sam.Aux
	if rgEmitted {
		finalAux = tags
	} else {
		finalAux = mergeAux(tags, rd.readGroupTag(rgValue))
	}
	if rd.rawAuxBAMSink {
		// Memory-lean view passthrough: serialise the assembled aux fields
		// straight to a raw on-disk BAM aux byte block and leave rec.Aux nil.
		// The bytes are byte-for-byte what encodeRecord would write for the
		// equivalent rec.Aux (both go through the same per-field serialiser), so
		// a record decoded here is identical to one decoded eagerly and
		// re-encoded. A trailing data-series RG is located later by a raw-byte
		// walk in regenerateMDNM so MD/NM splice in just before it.
		raw, rerr := rd.arena.buildRawAux(finalAux)
		if rerr != nil {
			return nil, wrapf(rerr, "record %d aux", index)
		}
		rec.RawAux = raw
	} else {
		rec.Aux = finalAux
	}
	return dr, nil
}

// storeSeq attaches the reconstructed SEQ bytes to rec. On the memory-lean
// view passthrough (rawAuxBAMSink) it packs seq directly into the on-disk BAM
// 4-bit nibble layout and stores it in rec.RawSeq + rec.SeqLen, leaving rec.Seq
// "" so the heavier SEQ string is never built — SEQ is the dominant fat
// per-record field, so this is where the decode's peak RSS is saved. The packed
// bytes are byte-for-byte what the BAM writer's PackSeqInto would emit for the
// equivalent rec.Seq, so a record decoded here writes identically to one decoded
// eagerly. Off the passthrough (the default) it stores the eager SEQ string,
// exactly as before. resolveMates touches only coordinate fields, so packing SEQ
// during decodeRecord is safe; mdnm reads the packed bases through sam.SeqBaseAt.
func (rd *recordDecoder) storeSeq(rec *sam.Record, seq []byte) {
	if rd.rawAuxBAMSink {
		// Pack into the per-slice SEQ arena (capped sub-slice) instead of a fresh
		// per-record make(); byte-identical to the standalone pack.
		rec.RawSeq = rd.arena.packSeq(seq, len(seq))
		rec.SeqLen = len(seq)
		rec.Seq = ""
		return
	}
	rec.Seq = string(seq)
}

// buildRawAux serialises an assembled aux list into a raw on-disk BAM aux byte
// block carved from the per-slice aux arena. The bytes are byte-for-byte what
// the BAM writer's encodeRecord emits for the same aux list — both use the
// shared per-field serialiser (sam.AppendBAMAux / encodeBAMAux) — so a record
// carrying this RawAux writes identically to one carrying the equivalent
// rec.Aux. The returned slice is capped at its exact length (a 3-index slice)
// so the later MD/NM raw splice, which appends, reallocates rather than writing
// into a neighbouring record's arena region.
func (a *recordArena) buildRawAux(aux []sam.Aux) ([]byte, error) {
	if len(aux) == 0 {
		return nil, nil
	}
	start := len(a.aux)
	for i := range aux {
		var err error
		// Append directly onto a.aux so its length keeps growing across records;
		// the record's block is the [start:end] window. (Appending into a
		// reslice that starts at `start` and reassigning would drop the prefix —
		// the window indices must stay relative to the full arena buffer.)
		a.aux, err = sam.AppendBAMAux(a.aux, aux[i])
		if err != nil {
			return nil, err
		}
	}
	end := len(a.aux)
	return a.aux[start:end:end], nil
}

// mergeAux appends the read-group aux tag, reconstructed from the RG
// data series, to the dictionary tag list. It mirrors htslib's
// cram_decode.c: when RG appears in the tag dictionary it is emitted in
// its dictionary position (so the data-series RG is suppressed here),
// otherwise the data-series RG tag is appended after every dictionary
// tag. A nil rg adds nothing.
func mergeAux(dictTags []sam.Aux, rg *sam.Aux) []sam.Aux {
	if rg == nil {
		return dictTags
	}
	for i := range dictTags {
		if dictTags[i].Tag == "RG" {
			// The dictionary already carries RG; htslib emits it in place
			// and does not also append the data-series value.
			return dictTags
		}
	}
	return append(dictTags, *rg)
}

// reverseQual reverses a quality slice in place. It is applied to a
// reverse-strand read when the container's QO preservation flag is off,
// so the emitted QUAL is in read orientation. An all-0xff "no quality"
// slice is unchanged by the reversal.
func reverseQual(q []byte) {
	for i, j := 0, len(q)-1; i < j; i, j = i+1, j-1 {
		q[i], q[j] = q[j], q[i]
	}
}

// decodeSliceRecords decodes every record of one slice and resolves the
// downstream-mate links between them. CRAM stores a read pair compactly
// when the two mates fall in the same slice: the upstream mate carries a
// next-fragment distance instead of explicit mate coordinates, and the
// mate-related SAM fields are reconstructed here by pairing each such
// record with the record that distance away.
func (rd *recordDecoder) decodeSliceRecords(nRecords int32) ([]*sam.Record, []mdnmSuppress, error) {
	if nRecords < 0 {
		return nil, nil, errFormat("slice declares a negative record count %d", nRecords)
	}
	// Fail fast on a count grossly larger than even the decompressed
	// series data, before the loop allocates anything.
	if total := rd.src.s.totalBytes(); nRecords > total {
		return nil, nil, errFormat("slice declares %d records but holds only %d bytes of series data",
			nRecords, total)
	}
	if rd.arena != nil {
		// Rewind the reused byte runs to empty, then (re)size the record slabs to
		// the authoritative count. Cross-slice reuse: the same backing buffers
		// serve every slice, so steady-state decode allocates almost nothing.
		rd.arena.reset()
		rd.arena.reserve(int(nRecords))
	}
	decoded := make([]*decodedRecord, 0, nRecords)
	prev := rd.src.s.consumed()
	for i := int32(0); i < nRecords; i++ {
		dr, err := rd.decodeRecord(int(i))
		if err != nil {
			return nil, nil, err
		}
		decoded = append(decoded, dr)
		// Every record must consume series input. If one did not, the
		// declared count has outrun the data — stop rather than loop to
		// nRecords emitting identical zero-byte records (a crafted
		// header could otherwise drive a multi-billion-iteration loop).
		c := rd.src.s.consumed()
		if c == prev && i+1 < nRecords {
			return nil, nil, errFormat("slice declares %d records but the series data is exhausted after record %d",
				nRecords, i)
		}
		prev = c
	}
	if err := resolveMates(decoded); err != nil {
		return nil, nil, err
	}
	rd.reconstructDroppedNames(decoded)
	out := make([]*sam.Record, len(decoded))
	suppress := make([]mdnmSuppress, len(decoded))
	for i, dr := range decoded {
		out[i] = dr.rec
		suppress[i] = mdnmSuppress{md: dr.suppressMD, nm: dr.suppressNM}
	}
	return out, suppress, nil
}

// reconstructDroppedNames assigns a QNAME to every record of a lossy-names
// slice whose name was dropped at encode time (nameStored == false). It
// runs after resolveMates, so a downstream-mate pair has already had the
// upstream record's name copied to its mate (linkMates). For any record
// still without a name, it synthesises "<prefix>:<record_number>", exactly
// as htslib's cram_to_bam does: the record number is the slice's running
// record counter plus the in-slice index (or the mate's index when the
// mate is an earlier record in the same slice), plus one. With names
// preserved every record's name was stored, so this is a no-op.
func (rd *recordDecoder) reconstructDroppedNames(decoded []*decodedRecord) {
	if rd.h.Preservation.ReadNamesIncluded {
		return
	}
	for i, dr := range decoded {
		if dr.nameStored || dr.rec.QName != "" {
			// A name was stored, or resolveMates copied one from the mate.
			continue
		}
		// htslib synthesises from the mate's index when the mate is an
		// earlier record in this slice, otherwise from this record's index.
		// Both members of a within-slice pair therefore resolve to the same
		// number — the upstream record's — so a dropped pair shares a name.
		recNo := int64(i)
		if dr.mateIndex >= 0 && dr.mateIndex < i {
			recNo = int64(dr.mateIndex)
		}
		num := rd.slice.RecordCounter + recNo + 1
		if rd.namePrefix != "" {
			dr.rec.QName = fmt.Sprintf("%s:%d", rd.namePrefix, num)
		} else {
			dr.rec.QName = fmt.Sprintf("%d", num)
		}
	}
}

// NeedsReference reports whether decoding has so far reached a record
// whose sequence an external reference would supply. A reference-free
// CRAM file (the common case for unaligned data and the decode-to-SAM
// oracle) never sets it; for a reference-backed file the bases an
// external reference would provide are filled with 'N' and this returns
// true so a caller knows the sequences are incomplete.
func (rr *RecordReader) NeedsReference() bool { return rr.needsReference }

// resolveMates fills in the RNEXT, PNEXT and TLEN fields and the
// mate-related SAM flag bits for every record whose mate is stored
// downstream in the same slice. The upstream record's next-fragment
// distance names the mate; the two records' coordinates then yield the
// template length (TLEN) with the sign convention SAM mandates — the
// leftmost mate gets the positive value.
func resolveMates(decoded []*decodedRecord) error {
	for i, dr := range decoded {
		if !dr.mateDownstream {
			continue
		}
		mateIdx := i + int(dr.nextFragment) + 1
		if mateIdx < 0 || mateIdx >= len(decoded) {
			return errFormat("record %d names a downstream mate at index %d, outside the slice (%d records)",
				i, mateIdx, len(decoded))
		}
		// Record the link in both directions so a lossy-names slice can
		// synthesise a dropped name from the mate's index, as htslib does
		// (it sets the downstream record's mate_line back to the upstream).
		dr.mateIndex = mateIdx
		decoded[mateIdx].mateIndex = i
		linkMates(dr.rec, decoded[mateIdx].rec)
	}
	return nil
}

// linkMates cross-fills the mate fields of an upstream record and its
// downstream mate. Each record's RNEXT/PNEXT points at the other, the
// mate-reverse and mate-unmapped flag bits are copied from the mate's
// own flags, and the template length is computed from the pair's span.
func linkMates(up, down *sam.Record) {
	// CRAM v4 deduplicates a within-slice pair's read name: the upstream
	// mate stores an empty name (a bare stop byte) and the decoder copies
	// it from the downstream mate (htslib cram_decode.c, the
	// `cr->name_len == 0` mate-copy branch). A pair always shares a QNAME,
	// so copying either way that is non-empty is safe; older CRAM versions
	// store both names and never hit this.
	if up.QName == "" && down.QName != "" {
		up.QName = down.QName
	} else if down.QName == "" && up.QName != "" {
		down.QName = up.QName
	}

	setMateFields(up, down)
	setMateFields(down, up)

	// TLEN spans from the leftmost mapped base of the pair to the rightmost.
	// The magnitude is max(up.end, down.end) - min(up.pos, down.pos) + 1: when
	// the lower-POS mate soft-clips or indels further right than the higher-POS
	// mate, that lower mate — not the higher one — sets the right end, so taking
	// the right end from only the higher-POS mate would understate the span.
	//
	// The sign/tie-break must reproduce htslib's cram_decode.c mate
	// cross-reference pass EXACTLY so a decoded within-slice pair matches
	// upstream byte-for-byte; this is the decode-side mirror of writeencode.go
	// computeTLenOverrides (keep the two in sync). `up` is the earlier (lower
	// in-slice index) record — htslib's `p` — and `down` its later mate (`cr`).
	if up.Flag&sam.FlagUnmapped == 0 && down.Flag&sam.FlagUnmapped == 0 && up.RName == down.RName {
		aleft := up.Pos
		if down.Pos < aleft {
			aleft = down.Pos
		}
		aright := up.EndPosition()
		if down.EndPosition() > aright {
			aright = down.EndPosition()
		}
		span := aright - aleft + 1
		switch {
		case up.Pos == aleft && (up.EndPosition() < aright || up.Pos != down.Pos):
			// up is leftmost and is either the strictly-lower start or not the
			// sole rightmost mate: it takes the positive span.
			up.TLen, down.TLen = span, -span
		case up.Pos == down.Pos && up.EndPosition() == down.EndPosition():
			// Full overlap (equal start AND equal end): resolve via READ1 so the
			// sign is independent of record order.
			if up.Flag&sam.FlagRead1 != 0 {
				up.TLen, down.TLen = span, -span
			} else {
				up.TLen, down.TLen = -span, span
			}
		default:
			// down is leftmost, or up is the sole rightmost mate at an equal
			// start: down takes the positive span.
			up.TLen, down.TLen = -span, span
		}
	}
}

// setMateFields points rec's RNEXT and PNEXT at its mate, and copies the
// mate's strand and mapped state into rec's mate-related flag bits.
func setMateFields(rec, mate *sam.Record) {
	if mate.RName != "" && mate.RName == rec.RName {
		rec.RNext = "="
	} else {
		rec.RNext = mate.RName
	}
	rec.PNext = mate.Pos
	if mate.Flag&sam.FlagReverse != 0 {
		rec.Flag |= sam.FlagMateReverse
	}
	if mate.Flag&sam.FlagUnmapped != 0 {
		rec.Flag |= sam.FlagMateUnmapped
	}
}

// decodeReadName reconstructs the record's QNAME. When the preservation
// map's RN flag is set, the name is read from the RN data series at this
// (the normal read-name) position in the record layout.
//
// When read names are NOT preserved (a lossy-names CRAM), htslib reads no
// RN value here: a detached record's name is read later, inside the mate
// block (decodeMate), and a non-detached record's name was dropped and is
// reconstructed once the whole slice is decoded (reconstructDroppedNames).
// This mirrors cram_decode.c, which only reads the RN series at the
// read-name position when comp_hdr->read_names_included is set.
func (rd *recordDecoder) decodeReadName(dr *decodedRecord, cf int32, index int) error {
	if rd.h.Preservation.ReadNamesIncluded {
		// The name is copied straight into QName below, so borrow it from the
		// block buffer instead of allocating a per-record copy.
		name, err := rd.byteArraySeriesBorrow("RN")
		if err != nil {
			return wrapf(err, "record %d read name", index)
		}
		dr.rec.QName = string(name)
		dr.nameStored = true
		return nil
	}
	// Names not preserved: defer. decodeMate fills a detached record's name
	// from the mate block; reconstructDroppedNames handles the rest.
	return nil
}

// optByteArraySeries decodes one variable-length byte-array value of the
// named data series, or returns ok=false when the series is absent from
// the encoding map.
func (rd *recordDecoder) optByteArraySeries(key string) ([]byte, bool, error) {
	enc := rd.h.Encoding(key)
	if enc == nil {
		return nil, false, nil
	}
	b, err := enc.decodeByteArray(rd.src.s)
	if err != nil {
		return nil, false, wrapf(err, "data series %q", key)
	}
	return b, true, nil
}

// optByteArraySeriesBorrow is optByteArraySeries returning a NO-COPY sub-slice
// of the block buffer (see byteArraySeriesBorrow's retention contract): the
// bytes must be copied out before the slice's source is released.
func (rd *recordDecoder) optByteArraySeriesBorrow(key string) ([]byte, bool, error) {
	enc := rd.h.Encoding(key)
	if enc == nil {
		return nil, false, nil
	}
	b, err := enc.decodeByteArrayBorrow(rd.src.s)
	if err != nil {
		return nil, false, wrapf(err, "data series %q", key)
	}
	return b, true, nil
}

// decodeMate reconstructs the mate-related fields (RNEXT, PNEXT, TLEN
// and the mate SAM flag bits). A detached record (CF bit 0x2) stores its
// own mate information; a non-detached record's mate is the next record
// downstream and NF holds the distance to it — that resolution is left
// to resolveMates, so here only the NF distance is captured.
func (rd *recordDecoder) decodeMate(dr *decodedRecord, cf int32, index int) error {
	rec := dr.rec
	switch {
	case cf&cfDetached != 0:
		mf, err := rd.intSeries("MF")
		if err != nil {
			return wrapf(err, "record %d mate flags", index)
		}
		if mf&mfMateReverse != 0 {
			rec.Flag |= sam.FlagMateReverse
		}
		if mf&mfMateUnmapped != 0 {
			rec.Flag |= sam.FlagMateUnmapped
		}
		// A detached record whose read names are not preserved carries its
		// name inside the mate block — after the MF series, before NS —
		// rather than at the normal read-name position. htslib's
		// cram_decode.c reads the RN series here when
		// !read_names_included, so a lossy-names CRAM keeps the real names
		// of its detached records (only fully-reconstructable non-detached
		// duplicates are dropped). When the RN series is present in the
		// encoding map, read it; when it is absent the name was dropped and
		// reconstructDroppedNames synthesises it after the slice decodes.
		if !rd.h.Preservation.ReadNamesIncluded {
			name, ok, nerr := rd.optByteArraySeriesBorrow("RN")
			if nerr != nil {
				return wrapf(nerr, "record %d detached read name", index)
			}
			if ok {
				rec.QName = string(name)
				dr.nameStored = true
			}
		}
		nsID, err := rd.intSeries("NS")
		if err != nil {
			return wrapf(err, "record %d mate reference id", index)
		}
		np, err := rd.intSeries("NP")
		if err != nil {
			return wrapf(err, "record %d mate alignment position", index)
		}
		ts, err := rd.intSeries("TS")
		if err != nil {
			return wrapf(err, "record %d template size", index)
		}
		mateName, nerr := rd.refName(nsID)
		if nerr != nil {
			return wrapf(nerr, "record %d mate", index)
		}
		// RNEXT is "=" when the mate maps to the same reference, matching
		// the SAM convention samtools emits.
		if mateName != "" && mateName == rec.RName {
			rec.RNext = "="
		} else {
			rec.RNext = mateName
		}
		rec.PNext = int64(np)
		rec.TLen = int64(ts)
	case cf&cfHasMateDownstream != 0:
		// The mate is a later record in the slice; NF gives the distance.
		// Its fields are filled in by resolveMates once the whole slice
		// is decoded.
		nf, err := rd.intSeries("NF")
		if err != nil {
			return wrapf(err, "record %d next-fragment distance", index)
		}
		if nf < 0 {
			return errFormat("record %d declares a negative next-fragment distance %d", index, nf)
		}
		dr.mateDownstream = true
		dr.nextFragment = nf
	}
	return nil
}

// decodeTags reconstructs the record's auxiliary tags from its tag-line
// (TL) value: TL indexes the tag dictionary, whose entry lists the tags
// the record carries, and each tag's value is then read from its
// per-tag data series.
//
// rgValue is the record's RG data-series value, threaded in for the CRAM
// v4 case where the read-group tag appears in the tag dictionary as a
// placeholder "RG*" (type byte '*'): the value is then emitted in its
// dictionary position rather than appended afterwards. The returned
// rgEmitted flag tells the caller to suppress the separate RG append in
// that case. CRAM v2/v3 leave RG out of the dictionary, so rgEmitted is
// false and the caller appends it as before.
func (rd *recordDecoder) decodeTags(tl int32, rgValue int32, index int) (aux []sam.Aux, rgEmitted bool, suppress mdnmSuppress, err error) {
	if tl < 0 {
		return nil, false, suppress, errFormat("record %d declares a negative tag-line index %d", index, tl)
	}
	if int(tl) >= len(rd.tagDict) {
		// A TL of 0 with an empty dictionary means "no tags".
		if tl == 0 {
			return nil, false, rd.mdnmSuppressFor(cramFlagsTag{}, false, false), nil
		}
		return nil, false, suppress, errFormat("record %d tag-line index %d has no tag-dictionary entry (%d entries)",
			index, tl, len(rd.tagDict))
	}
	keys := rd.tagDict[tl]
	// Over-allocate by one: mergeAux may append a synthesised RG aux to this
	// slice, and sizing the capacity to len(keys)+1 lets that append reuse
	// the backing array instead of reallocating and copying every tag once
	// per record (a measurable slice of the CRAM-decode allocation churn).
	var out []sam.Aux
	if rd.rawAuxBAMSink {
		// Reuse the per-arena header scratch: on the rawAuxBAMSink path the tags are
		// serialised into the aux byte arena (buildRawAux) before the next record's
		// decode reuses this header slice, so a fresh make() per record is not
		// needed. Keep the larger backing across records (and slices). This reuse is
		// gated on rawAuxBAMSink (NOT merely arena != nil) because on the eager arena
		// path `out` becomes rec.Aux verbatim — reusing one backing slice there would
		// have every record's aux alias, and be overwritten by, the next record's.
		if cap(rd.arena.auxScratch) < len(keys)+1 {
			rd.arena.auxScratch = make([]sam.Aux, 0, len(keys)+1)
		}
		out = rd.arena.auxScratch[:0]
	} else {
		// Eager path: `out` is retained as rec.Aux, so it gets its own backing.
		// Size it for the dictionary tags plus the three appends that follow
		// downstream without a reallocation — the data-series RG (mergeAux) and the
		// regenerated MD and NM (insertBeforeTrailingRG) — so those splices reuse
		// this backing in place instead of allocating a second, larger aux slice per
		// record (the standalone realloc was the single largest resident allocation
		// of eager CRAM decode). The emitted tags and their order are unchanged.
		out = make([]sam.Aux, 0, len(keys)+3)
	}
	// For CRAM v4 the presence of an MD*/NM* placeholder in this record's
	// tag-line dictionary is what authorises regeneration; for v<4 the cF
	// aux-tag bits authorise suppression. Both are folded into the final
	// mdnmSuppress by mdnmSuppressFor below.
	var cf cramFlagsTag
	var seenMDStar, seenNMStar bool
	for _, key := range keys {
		// CRAM v4 stores certain tags in the dictionary as placeholders
		// whose type byte is '*' (0x2a), to be filled in by the decoder
		// rather than read from a per-tag block (cram_decode.c, the
		// `TN[2] == '*'` branch). RG* is the read group; NM*/MD* are
		// MD/NM regeneration placeholders.
		if rd.h.major >= 4 && key[2] == '*' {
			if key[0] == 'R' && key[1] == 'G' {
				if rg := rd.readGroupTag(rgValue); rg != nil {
					out = append(out, *rg)
				}
				rgEmitted = true
				continue
			}
			// NM*/MD* carry no per-tag block: they are regenerated from the
			// reference if one is attached (see iterator.go regenerateMDNM),
			// so nothing is read or emitted here. Their *presence* is what a
			// v4 file uses to request regeneration (cram_decode.c:2046 sets
			// has_MD = -pos for the placeholder, and the v4 gate regenerates
			// only when has_MD < 0); record it so a record WITHOUT the
			// placeholder is left bare, matching upstream.
			if key[0] == 'M' && key[1] == 'D' {
				seenMDStar = true
			}
			if key[0] == 'N' && key[1] == 'M' {
				seenNMStar = true
			}
			continue
		}
		enc := rd.h.Tags[key]
		if enc == nil {
			return nil, false, suppress, errFormat("record %d tag %q is not in the tag-encoding map", index, key.String())
		}
		// raw is consumed immediately — isCRAMFlagsTag reads it and decodeTagValue
		// copies it into aux.Value (string() for Z/H/A, a numeric parse otherwise,
		// a fresh []interface{} for B) — and never retained or mutated, so borrow it
		// from the block buffer instead of allocating a per-tag copy.
		raw, derr := enc.decodeByteArrayBorrow(rd.src.s)
		if derr != nil {
			return nil, false, suppress, wrapf(derr, "record %d tag %q value", index, key.String())
		}
		// "cF" (CRAM flags, a single byte typed 'C') is an internal tag
		// htslib writes and then strips on read — it is not part of the
		// SAM record. Its value bytes must still be drained from the data
		// series (done above), but the tag is never emitted. Its bits gate
		// per-record MD/NM regeneration: bit 1 = "source had no MD", bit 2 =
		// "no NM" (cram_encode.c:3417, set only by an embed_ref=2 writer).
		// htslib turns these into has_MD=1/has_NM=1 to suppress regeneration
		// for that record (cram_decode.c:2117-2122); we capture them here so
		// mdnmSuppressFor can reproduce the same per-record suppression.
		// (htslib cram_decode.c "Remove cF tag".)
		if isCRAMFlagsTag(key, raw) {
			cf = parseCRAMFlagsTag(raw)
			continue
		}
		val, perr := decodeTagValue(key, raw)
		if perr != nil {
			return nil, false, suppress, wrapf(perr, "record %d tag %q", index, key.String())
		}
		out = append(out, val)
	}
	if rd.rawAuxBAMSink {
		rd.arena.auxScratch = out // keep any growth for the next record.
	}
	return out, rgEmitted, rd.mdnmSuppressFor(cf, seenMDStar, seenNMStar), nil
}

// mdnmSuppressFor folds the version-specific MD/NM regeneration policy into
// a single per-record mdnmSuppress.
//
// CRAM < 4: regeneration is the default (the gate `do_md` is true because
// `samtools view` sets decode_md=-1, so `do_md = decode_md != 0` is true).
// The per-record cF bits SUPPRESS it for a record whose source lacked the
// tag — htslib forces has_MD/has_NM, making `!has_MD` false
// (cram_decode.c:2117-2122). So suppress = the cF bits.
//
// CRAM >= 4: regeneration is OFF by default (`do_md = decode_md > 0` is
// false for the -1 "auto" setting), and is requested ONLY when an MD*/NM*
// placeholder sits in the record's tag dictionary (has_MD < 0 at
// cram_decode.c:2046). So suppress = NOT the placeholder presence.
func (rd *recordDecoder) mdnmSuppressFor(cf cramFlagsTag, seenMDStar, seenNMStar bool) mdnmSuppress {
	if rd.h.major >= 4 {
		return mdnmSuppress{md: !seenMDStar, nm: !seenNMStar}
	}
	return mdnmSuppress{md: cf.noMD, nm: cf.noNM}
}

// isCRAMFlagsTag reports whether a tag key/value is the internal "cF"
// CRAM-flags tag (a single byte typed 'C'), which htslib drains but does
// not surface as a SAM auxiliary field. It mirrors the
// `TN[-3]=='c' && TN[-2]=='F' && TN[-1]=='C' && out_sz == 1` guard in
// htslib's cram_decode.c.
func isCRAMFlagsTag(key tagKey, raw []byte) bool {
	return key[0] == 'c' && key[1] == 'F' && key[2] == 'C' && len(raw) == 1
}

// cramFlagsTag holds the decoded per-record "cF" (CRAM-flags) aux-tag
// suppression bits an embed_ref=2 writer stores. noMD mirrors cF bit 1
// ("the source record had no MD") and noNM mirrors cF bit 2 ("no NM"),
// exactly as cram_encode.c:3417 sets them and cram_decode.c:2117-2122
// consumes them. A zero cramFlagsTag (the common case: no cF tag, or a
// cF tag with no suppression bits) means "do not suppress".
type cramFlagsTag struct {
	noMD bool
	noNM bool
}

// parseCRAMFlagsTag decodes the single cF value byte into its MD/NM
// suppression bits. raw is the one-byte payload isCRAMFlagsTag matched;
// bit 1 (0x1) is "no MD" and bit 2 (0x2) is "no NM".
func parseCRAMFlagsTag(raw []byte) cramFlagsTag {
	if len(raw) != 1 {
		return cramFlagsTag{}
	}
	return cramFlagsTag{
		noMD: raw[0]&0x1 != 0,
		noNM: raw[0]&0x2 != 0,
	}
}

// readGroupTag turns a CRAM RG data-series value into the SAM "RG:Z:"
// aux tag, or nil when the record has no read group (value -1) or the
// index has no @RG header entry. An out-of-range index is dropped rather
// than treated as fatal: the alignment fields stay valid and only this
// one optional tag is affected.
func (rd *recordDecoder) readGroupTag(rgValue int32) *sam.Aux {
	if rgValue < 0 || int(rgValue) >= len(rd.readGroups) {
		return nil
	}
	return &sam.Aux{
		Tag:   "RG",
		Type:  'Z',
		Value: rd.readGroups[rgValue],
	}
}

// decodeMapped reconstructs the SEQ, QUAL and CIGAR of a mapped record
// from its read-feature list. In a reference-free file every base is
// carried by a feature, so the reconstruction needs no external
// sequence; a reference-backed file would copy reference bases for the
// runs the features do not cover.
func (rd *recordDecoder) decodeMapped(rec *sam.Record, cf int32, readLen int32, index int) error {
	fn, err := rd.intSeries("FN")
	if err != nil {
		return wrapf(err, "record %d feature count", index)
	}
	feats, err := rd.decodeFeatures(fn)
	if err != nil {
		return wrapf(err, "record %d", index)
	}
	seq, qual, cigar, err := rd.reconstructMapped(feats, readLen, int32(rec.Pos))
	if err != nil {
		return wrapf(err, "record %d", index)
	}
	rd.storeSeq(rec, seq)
	rec.Cigar = cigar

	// A mapped record's quality is read as one block of readLen scores
	// only when the CRAM flags say quality was preserved separately;
	// otherwise the per-base Q features (already gathered above) carry it.
	if cf&cfQualityPreserved != 0 {
		// On the arena path reconstructMapped already carved a readLen QUAL
		// buffer from the per-slice arena; the preserved block overwrites it in
		// place so no second arena buffer is consumed. Off the arena path qual is
		// a fresh make() reused the same way. Either way the bytes equal a fresh
		// readQualityBlock(readLen).
		q, qerr := rd.readQualityBlockInto(qual, readLen)
		if qerr != nil {
			return wrapf(qerr, "record %d quality", index)
		}
		rec.Qual = q
	} else {
		rec.Qual = qual
	}
	return nil
}

// decodeUnmapped reconstructs an unmapped record: its bases come
// verbatim from the BA data series and its qualities from QS. An
// unmapped record has no CIGAR and a position of zero.
func (rd *recordDecoder) decodeUnmapped(rec *sam.Record, cf int32, readLen int32, index int) error {
	// On the gated path the bases are read into the reused per-record scratch and
	// packed into the SEQ arena by storeSeq; off it they go into a fresh make().
	var seq []byte
	if rd.arena != nil {
		seq, _ = rd.arena.reconScratch(int(readLen))
	} else {
		seq = make([]byte, readLen)
	}
	// Fast path mirrors readQualityBlockInto: the BA (base) series is normally an
	// EXTERNAL byte block, so resolve its encoding and cursor once and read every
	// base directly from the cursor rather than through two map lookups per base.
	// The bytes and error wrapping match the per-base byteSeries("BA") loop it
	// replaces (byteSeries wraps a cursor/read error as `data series "BA"`, then
	// decodeUnmapped wraps that as `record N base i`).
	if enc := rd.baEnc; enc != nil && enc.ID == EncodingExternal {
		c, cerr := rd.src.s.cursor(enc.ExternalID)
		if cerr != nil {
			return wrapf(wrapf(cerr, "data series %q", "BA"), "record %d base %d", index, 0)
		}
		for i := int32(0); i < readLen; i++ {
			b, rerr := c.readByte()
			if rerr != nil {
				return wrapf(wrapf(rerr, "data series %q", "BA"), "record %d base %d", index, i)
			}
			seq[i] = b
		}
	} else {
		for i := int32(0); i < readLen; i++ {
			b, err := rd.byteSeries("BA")
			if err != nil {
				return wrapf(err, "record %d base %d", index, i)
			}
			seq[i] = b
		}
	}
	rd.storeSeq(rec, seq)
	if cf&cfQualityPreserved != 0 {
		q, err := rd.readQualityBlock(readLen)
		if err != nil {
			return wrapf(err, "record %d quality", index)
		}
		rec.Qual = q
	}
	return nil
}

// readQualityBlock reads readLen quality scores from the QS data series into a
// QUAL buffer: a fresh per-record make() off the gated path, or a per-slice
// arena buffer (capped at readLen) on it. The bytes are identical either way.
func (rd *recordDecoder) readQualityBlock(readLen int32) ([]byte, error) {
	var q []byte
	if rd.arena != nil {
		q = rd.arena.qualFor(int(readLen))
	} else {
		q = make([]byte, readLen)
	}
	return rd.readQualityBlockInto(q, readLen)
}

// readQualityBlockInto fills dst (which must be at least readLen long) with
// readLen quality scores from the QS data series and returns dst[:readLen]. It
// lets a caller reuse an already-carved QUAL buffer (e.g. the mapped
// reconstruction's arena qual) for the preserved-quality block read.
func (rd *recordDecoder) readQualityBlockInto(dst []byte, readLen int32) ([]byte, error) {
	q := dst[:readLen]
	// Fast path: the QS series is almost always an EXTERNAL byte block. Resolve
	// its encoding (memoised on the decoder) and its block cursor ONCE, then read
	// every quality base straight from the cursor — rather than re-resolving both
	// the encoding and the cursor through two map lookups on EVERY base, which
	// made this the single hottest loop of CRAM decode. The cursor is the same
	// cached, stateful cursor byteSeries("QS") would fetch on each iteration, so
	// it advances identically and the bytes (and the error wrapping) are exactly
	// what the per-base path produced.
	if enc := rd.qsEnc; enc != nil && enc.ID == EncodingExternal {
		c, err := rd.src.s.cursor(enc.ExternalID)
		if err != nil {
			return nil, wrapf(err, "data series %q", "QS")
		}
		for i := int32(0); i < readLen; i++ {
			b, rerr := c.readByte()
			if rerr != nil {
				return nil, wrapf(rerr, "data series %q", "QS")
			}
			q[i] = b
		}
		return q, nil
	}
	// Any other QS encoding (HUFFMAN/BETA/CONST, or an absent series) takes the
	// unchanged per-base path, byte-for-byte as before.
	for i := int32(0); i < readLen; i++ {
		b, err := rd.byteSeries("QS")
		if err != nil {
			return nil, err
		}
		q[i] = b
	}
	return q, nil
}
