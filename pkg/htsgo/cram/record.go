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
}

// newRecordDecoder builds a recordDecoder for one slice. It parses the
// tag dictionary up front so per-record tag reconstruction is a lookup.
func newRecordDecoder(h *CompressionHeader, sh *SliceHeader, src *SeriesSource, refNames, readGroups []string) (*recordDecoder, error) {
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
	dr := &decodedRecord{rec: rec}

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
		rec.Pos = rd.prevAlignmentStart
	} else {
		rec.Pos = ap
	}

	rgValue, err := rd.intSeries("RG")
	if err != nil {
		return nil, wrapf(err, "record %d read group", index)
	}

	if err := rd.decodeReadName(rec, cf, index); err != nil {
		return nil, err
	}
	if err := rd.decodeMate(dr, cf, index); err != nil {
		return nil, err
	}

	tl, err := rd.intSeries("TL")
	if err != nil {
		return nil, wrapf(err, "record %d tag line", index)
	}
	tags, err := rd.decodeTags(tl, index)
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

	// The auxiliary tags are assembled from the dictionary-stored tags
	// and the read-group data series. samtools emits the read-group tag
	// immediately before the program tag (or last when no program tag is
	// present), so the RG tag is spliced in at that position.
	rec.Aux = mergeAux(tags, rd.readGroupTag(rgValue))
	return dr, nil
}

// mergeAux interleaves the read-group aux tag into the dictionary tag
// list. samtools writes the RG:Z tag immediately before the program
// (PG) tag; when the record carries no PG tag, RG is appended last. A
// nil rg adds nothing.
func mergeAux(dictTags []sam.Aux, rg *sam.Aux) []sam.Aux {
	if rg == nil {
		return dictTags
	}
	out := make([]sam.Aux, 0, len(dictTags)+1)
	inserted := false
	for _, a := range dictTags {
		if !inserted && a.Tag == "PG" {
			out = append(out, *rg)
			inserted = true
		}
		out = append(out, a)
	}
	if !inserted {
		out = append(out, *rg)
	}
	return out
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
	decoded := make([]*decodedRecord, 0, nRecords)
	for i := int32(0); i < nRecords; i++ {
		dr, err := rd.decodeRecord(int(i))
		if err != nil {
			return nil, err
		}
		decoded = append(decoded, dr)
	}
	if err := resolveMates(decoded); err != nil {
		return nil, err
	}
	out := make([]*sam.Record, len(decoded))
	for i, dr := range decoded {
		out[i] = dr.rec
	}
	return out, nil
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
		linkMates(dr.rec, decoded[mateIdx].rec)
	}
	return nil
}

// linkMates cross-fills the mate fields of an upstream record and its
// downstream mate. Each record's RNEXT/PNEXT points at the other, the
// mate-reverse and mate-unmapped flag bits are copied from the mate's
// own flags, and the template length is computed from the pair's span.
func linkMates(up, down *sam.Record) {
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
// map's RN flag is set the name comes from the RN data series; otherwise
// the writer dropped the name and the decoder synthesises a numeric one
// from the slice's record counter, matching htslib's convention.
func (rd *recordDecoder) decodeReadName(rec *sam.Record, cf int32, index int) error {
	if rd.h.Preservation.ReadNamesIncluded {
		name, err := rd.byteArraySeries("RN")
		if err != nil {
			return wrapf(err, "record %d read name", index)
		}
		rec.QName = string(name)
		return nil
	}
	// Names were not preserved: htslib generates a name from the running
	// record number (slice record counter plus the in-slice index).
	rec.QName = fmt.Sprintf("%d", rd.slice.RecordCounter+int64(index))
	return nil
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
		// A detached record may still drop its read name; htslib reads RN
		// again here when names are not preserved. Names preserved is the
		// common case and is handled in decodeReadName, so nothing extra
		// is read for that path.
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
		rec.PNext = np
		rec.TLen = ts
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
func (rd *recordDecoder) decodeTags(tl int32, index int) ([]sam.Aux, error) {
	if tl < 0 {
		return nil, errFormat("record %d declares a negative tag-line index %d", index, tl)
	}
	if int(tl) >= len(rd.tagDict) {
		// A TL of 0 with an empty dictionary means "no tags".
		if tl == 0 {
			return nil, nil
		}
		return nil, errFormat("record %d tag-line index %d has no tag-dictionary entry (%d entries)",
			index, tl, len(rd.tagDict))
	}
	keys := rd.tagDict[tl]
	out := make([]sam.Aux, 0, len(keys))
	for _, key := range keys {
		enc := rd.h.Tags[key]
		if enc == nil {
			return nil, errFormat("record %d tag %q is not in the tag-encoding map", index, key.String())
		}
		raw, err := enc.decodeByteArray(rd.src.s)
		if err != nil {
			return nil, wrapf(err, "record %d tag %q value", index, key.String())
		}
		aux, perr := decodeTagValue(key, raw)
		if perr != nil {
			return nil, wrapf(perr, "record %d tag %q", index, key.String())
		}
		out = append(out, aux)
	}
	return out, nil
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
	seq, qual, cigar, err := rd.reconstructMapped(feats, readLen)
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
