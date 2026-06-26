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
}

// decodeRecord decodes one CRAM alignment record into a decodedRecord.
// The traversal follows the CRAM v3.0 record layout: the bit flags (BF)
// and CRAM flags (CF), the reference id for a multi-reference slice, the
// read length (RL), the alignment position (AP), the read group (RG),
// the read name, the mate information for a detached record, the tag
// values, and finally either the read-feature list (mapped) or the raw
// bases (unmapped).
func (rd *recordDecoder) decodeRecord(index int) (*decodedRecord, error) {
	rec := &sam.Record{}
	dr := &decodedRecord{rec: rec, mateIndex: -1}

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
	tags, rgEmitted, err := rd.decodeTags(tl, rgValue, index)
	if err != nil {
		return nil, err
	}

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
	if rgEmitted {
		rec.Aux = tags
	} else {
		rec.Aux = mergeAux(tags, rd.readGroupTag(rgValue))
	}
	return dr, nil
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
func (rd *recordDecoder) decodeSliceRecords(nRecords int32) ([]*sam.Record, error) {
	if nRecords < 0 {
		return nil, errFormat("slice declares a negative record count %d", nRecords)
	}
	// Fail fast on a count grossly larger than even the decompressed
	// series data, before the loop allocates anything.
	if total := rd.src.s.totalBytes(); nRecords > total {
		return nil, errFormat("slice declares %d records but holds only %d bytes of series data",
			nRecords, total)
	}
	decoded := make([]*decodedRecord, 0)
	prev := rd.src.s.consumed()
	for i := int32(0); i < nRecords; i++ {
		dr, err := rd.decodeRecord(int(i))
		if err != nil {
			return nil, err
		}
		decoded = append(decoded, dr)
		// Every record must consume series input. If one did not, the
		// declared count has outrun the data — stop rather than loop to
		// nRecords emitting identical zero-byte records (a crafted
		// header could otherwise drive a multi-billion-iteration loop).
		c := rd.src.s.consumed()
		if c == prev && i+1 < nRecords {
			return nil, errFormat("slice declares %d records but the series data is exhausted after record %d",
				nRecords, i)
		}
		prev = c
	}
	if err := resolveMates(decoded); err != nil {
		return nil, err
	}
	rd.reconstructDroppedNames(decoded)
	out := make([]*sam.Record, len(decoded))
	for i, dr := range decoded {
		out[i] = dr.rec
	}
	return out, nil
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

	// TLEN spans from the leftmost mapped base of the pair to the
	// rightmost; the upstream (smaller-coordinate) record carries the
	// positive value and the downstream record its negation.
	if up.Flag&sam.FlagUnmapped == 0 && down.Flag&sam.FlagUnmapped == 0 && up.RName == down.RName {
		left := up.Pos
		right := down.EndPosition()
		if down.Pos < left {
			left = down.Pos
			right = up.EndPosition()
		}
		span := right - left + 1
		if up.Pos <= down.Pos {
			up.TLen = span
			down.TLen = -span
		} else {
			up.TLen = -span
			down.TLen = span
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
		name, err := rd.byteArraySeries("RN")
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
			name, ok, nerr := rd.optByteArraySeries("RN")
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
func (rd *recordDecoder) decodeTags(tl int32, rgValue int32, index int) (aux []sam.Aux, rgEmitted bool, err error) {
	if tl < 0 {
		return nil, false, errFormat("record %d declares a negative tag-line index %d", index, tl)
	}
	if int(tl) >= len(rd.tagDict) {
		// A TL of 0 with an empty dictionary means "no tags".
		if tl == 0 {
			return nil, false, nil
		}
		return nil, false, errFormat("record %d tag-line index %d has no tag-dictionary entry (%d entries)",
			index, tl, len(rd.tagDict))
	}
	keys := rd.tagDict[tl]
	out := make([]sam.Aux, 0, len(keys))
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
			// NM*/MD* (and any other '*' placeholder) carry no per-tag
			// block: they are regenerated from the reference if one is
			// attached (see iterator.go regenerateMDNM), so nothing is
			// read or emitted here. This matches htslib inserting an empty
			// placeholder it later overwrites.
			continue
		}
		enc := rd.h.Tags[key]
		if enc == nil {
			return nil, false, errFormat("record %d tag %q is not in the tag-encoding map", index, key.String())
		}
		raw, derr := enc.decodeByteArray(rd.src.s)
		if derr != nil {
			return nil, false, wrapf(derr, "record %d tag %q value", index, key.String())
		}
		// "cF" (CRAM flags, a single byte typed 'C') is an internal tag
		// htslib writes and then strips on read — it is not part of the
		// SAM record. Its value bytes must still be drained from the data
		// series (done above), but the tag is never emitted. The bits
		// suppress auto-regeneration of MD/NM; this decoder does not
		// auto-generate either, so the suppression is implicit and only
		// the tag-drop matters here. (htslib cram_decode.c "Remove cF tag".)
		if isCRAMFlagsTag(key, raw) {
			continue
		}
		val, perr := decodeTagValue(key, raw)
		if perr != nil {
			return nil, false, wrapf(perr, "record %d tag %q", index, key.String())
		}
		out = append(out, val)
	}
	return out, rgEmitted, nil
}

// isCRAMFlagsTag reports whether a tag key/value is the internal "cF"
// CRAM-flags tag (a single byte typed 'C'), which htslib drains but does
// not surface as a SAM auxiliary field. It mirrors the
// `TN[-3]=='c' && TN[-2]=='F' && TN[-1]=='C' && out_sz == 1` guard in
// htslib's cram_decode.c.
func isCRAMFlagsTag(key tagKey, raw []byte) bool {
	return key[0] == 'c' && key[1] == 'F' && key[2] == 'C' && len(raw) == 1
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
	rec.Seq = string(seq)
	rec.Cigar = cigar

	// A mapped record's quality is read as one block of readLen scores
	// only when the CRAM flags say quality was preserved separately;
	// otherwise the per-base Q features (already gathered above) carry it.
	if cf&cfQualityPreserved != 0 {
		q, qerr := rd.readQualityBlock(readLen)
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
	seq := make([]byte, readLen)
	for i := int32(0); i < readLen; i++ {
		b, err := rd.byteSeries("BA")
		if err != nil {
			return wrapf(err, "record %d base %d", index, i)
		}
		seq[i] = b
	}
	rec.Seq = string(seq)
	if cf&cfQualityPreserved != 0 {
		q, err := rd.readQualityBlock(readLen)
		if err != nil {
			return wrapf(err, "record %d quality", index)
		}
		rec.Qual = q
	}
	return nil
}

// readQualityBlock reads readLen quality scores from the QS data series.
func (rd *recordDecoder) readQualityBlock(readLen int32) ([]byte, error) {
	q := make([]byte, readLen)
	for i := int32(0); i < readLen; i++ {
		b, err := rd.byteSeries("QS")
		if err != nil {
			return nil, err
		}
		q[i] = b
	}
	return q, nil
}
