package cram

import (
	"fmt"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// seriesBuffers accumulates the encoded bytes of every external data
// series of one slice. Each field is the verbatim payload of the
// external block that series is written to; an empty field means the
// series carries no values and its block is omitted.
//
// The integer series store one ITF-8 integer per value; the byte series
// store one raw byte per value; the *Len fields hold the ITF-8 lengths
// of the byte-array series whose values they precede. The layout is
// exactly what decode.go's EXTERNAL and BYTE_ARRAY_LEN decoders consume.
type seriesBuffers struct {
	bf, cf, ri, rl, ap, rg []byte // per-record integer series.
	rn                     []byte // read names, NUL-terminated (BYTE_ARRAY_STOP).
	mf, ns, np, ts         []byte // detached-mate integer series.
	tl, mq, fn, fp         []byte // tag-line, mapping quality, feature count/positions.
	fc                     []byte // read-feature codes (one byte each).
	bb, bbLen              []byte // base-stretch values and lengths.
	in, inLen              []byte // inserted-base values and lengths.
	sc, scLen              []byte // soft-clip values and lengths.
	dl, rs, pd, hc         []byte // deletion/skip/pad/hard-clip lengths.
	ba                     []byte // single bases (unmapped reads).
	bs                     []byte // base-substitution codes (feature 'X', reference-based encoding).
	qs                     []byte // quality scores.

	// tagLens and tagVals hold each auxiliary tag's BYTE_ARRAY_LEN series:
	// tagLens[key] is the run of ITF-8 value lengths and tagVals[key] is
	// the concatenated value bytes. A length-prefixed layout — rather than
	// a delimiter — is used because a fixed-width binary tag value can
	// contain any byte, so no stop byte would be unambiguous.
	tagLens map[tagKey][]byte
	tagVals map[tagKey][]byte
}

// newSeriesBuffers returns an empty seriesBuffers ready for appends.
func newSeriesBuffers() *seriesBuffers {
	return &seriesBuffers{
		tagLens: make(map[tagKey][]byte),
		tagVals: make(map[tagKey][]byte),
	}
}

// encodeContainer encodes a batch of records as one complete CRAM data
// container: a container header, a compression-header block, a slice-
// header block and the slice's external data blocks. version selects the
// per-block codec set (see chooseBlockCompression). binning selects the
// lossy quality-binning scheme applied to each record's QUAL (BinningNone
// leaves quality untouched). recordCounter is the running record total of
// all earlier containers.
func encodeContainer(version Version, binning QualityBinning, records []*sam.Record, refIndex map[string]int32, reference map[string][]byte, recordCounter int64) ([]byte, error) {
	if len(records) == 0 {
		return nil, fmt.Errorf("cram: cannot encode an empty container")
	}

	// Determine the slice's reference scope. A slice is single-reference
	// when every mapped record shares one reference id; otherwise it is a
	// multi-reference slice (-2) and each record carries its own RI.
	refID, multiRef := sliceRefScope(records, refIndex)

	enc := &recordEncoder{
		version:   version,
		intw:      newIntWriter(version),
		refIndex:  refIndex,
		multiRef:  multiRef,
		reference: reference,
		binning:   binning,
		buffers:   newSeriesBuffers(),
	}
	if err := enc.encodeAll(records); err != nil {
		return nil, err
	}

	// The tag dictionary lists, per distinct tag-combination, the ordered
	// tag keys; each record's TL series value indexes it.
	tagDict := enc.tagDictionary()

	compHeader := encodeCompressionHeader(version, multiRef, tagDict, enc.tagKeysSorted())
	compBlock := encodeBlock(version, ContentCompressionHeader, 0, compHeader)

	// Assemble the slice: a slice-header block followed by the slice's data
	// blocks. The simple writer never uses the CORE bitstream, but htslib's
	// decoder demands a slice's first data block (s->block[0]) be of type
	// CORE for EVERY CRAM version — without it `samtools view` fails with
	// "Failure to decode slice". So every slice leads with an empty CORE
	// block. It is counted in num_blocks but not listed among the
	// block-content ids, matching cram_encode.c.
	extBlocks, contentIDs := enc.buffers.blocks(version, enc.tagKeysSorted())
	dataBlocks := [][]byte{encodeBlock(version, ContentCoreData, 0, nil)}
	dataBlocks = append(dataBlocks, extBlocks...)

	startPos, span := sliceSpan(records)
	sliceHeader := encodeSliceHeader(version, sliceHeaderFields{
		refSeqID:       refID,
		alignmentStart: startPos,
		alignmentSpan:  span,
		numRecords:     int32(len(records)),
		recordCounter:  recordCounter,
		numBlocks:      int32(len(dataBlocks)),
		contentIDs:     contentIDs,
	})
	sliceHeaderBlock := encodeBlock(version, ContentMappedSlice, 0, sliceHeader)

	// The container body is the compression-header block, the slice-
	// header block, then the data blocks. The slice landmark is the byte
	// offset of the slice-header block from the start of the body.
	var body []byte
	body = append(body, compBlock...)
	landmark := int32(len(body))
	body = append(body, sliceHeaderBlock...)
	for _, db := range dataBlocks {
		body = append(body, db...)
	}

	numBlocks := int32(1 /*comp header*/ + 1 /*slice header*/ + len(dataBlocks))
	hdr := containerHeaderBytes(version, containerFields{
		length:        int32(len(body)),
		refSeqID:      refID,
		startPos:      startPos,
		alignmentSpan: span,
		numRecords:    int32(len(records)),
		recordCounter: recordCounter,
		numBases:      enc.numBases,
		numBlocks:     numBlocks,
		landmarks:     []int32{landmark},
	})

	out := make([]byte, 0, len(hdr)+len(body))
	out = append(out, hdr...)
	out = append(out, body...)
	return out, nil
}

// sliceRefScope determines the reference id a slice header records and
// whether the slice is multi-reference. A slice is single-reference when
// every record maps to the same reference (or all are unmapped);
// otherwise the slice id is -2 and each record stores its own RI.
func sliceRefScope(records []*sam.Record, refIndex map[string]int32) (refID int32, multiRef bool) {
	seen := int32(-3) // sentinel: no record inspected yet.
	for _, rec := range records {
		id := recordRefID(rec, refIndex)
		if seen == -3 {
			seen = id
			continue
		}
		if id != seen {
			return -2, true
		}
	}
	if seen == -3 {
		return -1, false
	}
	return seen, false
}

// recordRefID returns the reference id of a record: the zero-based @SQ
// index of its RName, or -1 when the record is unmapped or its RName is
// the SAM "*".
func recordRefID(rec *sam.Record, refIndex map[string]int32) int32 {
	if rec.RName == "" || rec.RName == "*" {
		return -1
	}
	if id, ok := refIndex[rec.RName]; ok {
		return id
	}
	return -1
}

// sliceSpan computes the alignment start and reference span a slice
// covers from its mapped records. An all-unmapped slice spans nothing
// and starts at position 0.
func sliceSpan(records []*sam.Record) (start, span int32) {
	minPos := int32(0)
	maxEnd := int32(0)
	any := false
	for _, rec := range records {
		if rec.Flag&sam.FlagUnmapped != 0 || rec.Pos <= 0 {
			continue
		}
		end := rec.EndPosition()
		if !any {
			minPos, maxEnd, any = rec.Pos, end, true
			continue
		}
		if rec.Pos < minPos {
			minPos = rec.Pos
		}
		if end > maxEnd {
			maxEnd = end
		}
	}
	if !any {
		return 0, 0
	}
	return minPos, maxEnd - minPos + 1
}

// sliceHeaderFields collects the values a slice header carries so
// encodeSliceHeader can serialise them.
type sliceHeaderFields struct {
	refSeqID       int32
	alignmentStart int32
	alignmentSpan  int32
	numRecords     int32
	recordCounter  int64
	// numBlocks is the count of the slice's data blocks INCLUDING the CORE
	// block; contentIDs lists only the EXTERNAL blocks' content ids. The two
	// differ by one for a CRAM v4 slice, which carries an (empty) CORE block
	// that htslib requires as s->block[0] but does not list among the
	// block-content ids. For v3 the writer emits no CORE block, so
	// numBlocks == len(contentIDs).
	numBlocks  int32
	contentIDs []int32
}

// encodeSliceHeader serialises a CRAM slice-header block payload: the
// reference scope, the record count and counter, the data-block count
// and content-id list, the embedded-reference id (-1, none) and a
// 16-byte zero reference MD5. It is the writer-side inverse of
// parseSliceHeader.
//
// For CRAM v2/v3 every field is ITF-8 / LTF-8. For CRAM v4.0 the ref_seq_id
// is a signed (zig-zag) uint7 varint, the alignment start and span are
// 64-bit uint7 varints, the record counter is a 64-bit uint7 varint, and
// the remaining counts and ids are unsigned uint7 varints — matching
// cram_decode.c's major>=4 slice-header read.
func encodeSliceHeader(version Version, f sliceHeaderFields) []byte {
	iw := newIntWriter(version)
	var b []byte
	// ref_seq_id is signed (-1 unmapped, -2 multi-reference).
	b = iw.s32(b, f.refSeqID)
	if version.usesUint7() {
		// The alignment start and span widen to 64-bit in v4.
		b = iw.u64(b, int64(f.alignmentStart))
		b = iw.u64(b, int64(f.alignmentSpan))
	} else {
		b = iw.u32(b, f.alignmentStart)
		b = iw.u32(b, f.alignmentSpan)
	}
	b = iw.u32(b, f.numRecords)
	b = iw.u64(b, f.recordCounter)
	// The slice header carries the block count (num_blocks, which includes
	// the CORE block) and, separately, the length of the block-content-id
	// array (num_content_ids, EXTERNAL blocks only). They coincide for v3
	// (no CORE block) and differ by one for v4.
	b = iw.u32(b, f.numBlocks)
	b = iw.u32(b, int32(len(f.contentIDs)))
	for _, id := range f.contentIDs {
		b = iw.u32(b, id)
	}
	// Embedded-reference block id (htslib ref_base_id), read as an UNSIGNED
	// varint. The "no embedded reference" sentinel differs by version: CRAM
	// v2/v3 use -1 (which ITF-8 round-trips as 0xffffffff), whereas CRAM v4
	// uses 0 — htslib writes ref_base_id with varint_put32 (unsigned) and
	// chooses 0 for the no-embedded-ref case (cram_encode.c, major>=4). A
	// signed -1 here would zig-zag to 1 and htslib would read it as a real
	// block id, derailing the decode.
	if version.usesUint7() {
		b = iw.u32(b, 0)
	} else {
		b = iw.u32(b, -1)
	}
	// A reference-free slice records an all-zero reference MD5; the
	// decoder only checks the MD5 when a reference source is attached.
	b = append(b, make([]byte, 16)...)
	return b
}

// putU appends an unsigned integer data-series value through the encoder's
// version-aware framing: ITF-8 for v2/v3, an unsigned uint7 varint for v4.
// It is used for every non-negative integer data series, whose v4 encoding
// is VARINT_UNSIGNED.
func (e *recordEncoder) putU(dst []byte, v int32) []byte { return e.intw.u32(dst, v) }

// putS appends a signed integer data-series value: ITF-8 (whose
// sign-extension round-trips a small negative) for v2/v3, a signed
// (zig-zag) uint7 varint for v4. It is used for the series that can carry a
// negative value (RI, RG, NS, TS), whose v4 encoding is VARINT_SIGNED.
func (e *recordEncoder) putS(dst []byte, v int32) []byte { return e.intw.s32(dst, v) }

// blocks turns the populated series buffers into on-disk external data
// blocks and the parallel list of their content ids. Only non-empty
// series produce a block; an empty series is omitted, which the reader
// tolerates (a series with no block carries no values for the slice).
// version selects the per-block codec set (see chooseBlockCompression).
// tagKeys lists the auxiliary tag keys in the order their content ids
// were assigned.
func (sb *seriesBuffers) blocks(version Version, tagKeys []tagKey) (data [][]byte, ids []int32) {
	add := func(id int32, payload []byte) {
		if len(payload) == 0 {
			return
		}
		data = append(data, encodeBlock(version, ContentExternal, id, payload))
		ids = append(ids, id)
	}
	add(cidBF, sb.bf)
	add(cidCF, sb.cf)
	add(cidRI, sb.ri)
	add(cidRL, sb.rl)
	add(cidAP, sb.ap)
	add(cidRG, sb.rg)
	add(cidRN, sb.rn)
	add(cidMF, sb.mf)
	add(cidNS, sb.ns)
	add(cidNP, sb.np)
	add(cidTS, sb.ts)
	add(cidTL, sb.tl)
	add(cidMQ, sb.mq)
	add(cidFN, sb.fn)
	add(cidFC, sb.fc)
	add(cidFP, sb.fp)
	add(cidBB, sb.bb)
	add(cidBBLen, sb.bbLen)
	add(cidBS, sb.bs)
	add(cidBA, sb.ba)
	add(cidQS, sb.qs)
	add(cidIN, sb.in)
	add(cidINLen, sb.inLen)
	add(cidSC, sb.sc)
	add(cidSCLen, sb.scLen)
	add(cidDL, sb.dl)
	add(cidRS, sb.rs)
	add(cidPD, sb.pd)
	add(cidHC, sb.hc)
	for i, key := range tagKeys {
		lenID, valID := tagContentIDs(i)
		add(lenID, sb.tagLens[key])
		add(valID, sb.tagVals[key])
	}
	return data, ids
}

// recordEncoder walks a batch of records and fills the per-series byte
// buffers. It also accumulates the set of distinct tag combinations (the
// tag dictionary) and the running base count.
type recordEncoder struct {
	// version is the CRAM format being written. It selects the integer
	// framing of every data series: ITF-8 for v2/v3, uint7 varints for v4.
	version Version
	// intw is the version-aware integer serialiser used for every integer
	// data series and length field.
	intw intWriter

	refIndex map[string]int32
	multiRef bool
	// reference maps a contig name to its full reference bases. When a
	// mapped record's contig is present here, encodeFeatures diffs the read
	// against the reference and emits a substitution feature only at each
	// mismatch (the matched bases are reconstructed from the reference on
	// decode), exactly as upstream CRAM does. When nil, or the contig is
	// absent, the writer falls back to the self-contained reference-free
	// encoding (every M/=/X run carried literally in a base-stretch feature).
	reference map[string][]byte
	// binning is the lossy quality-binning scheme applied to each
	// record's QUAL before it is appended to the QS series. BinningNone
	// leaves quality untouched.
	binning QualityBinning
	buffers *seriesBuffers

	// numBases is the running total of read bases, stored in the
	// container header.
	numBases int64

	// tagCombos lists each distinct ordered tag-key combination in first-
	// seen order; comboIndex maps a combination's string form to its
	// index. A record's TL value is its combination's index.
	tagCombos  [][]tagKey
	comboIndex map[string]int32
	// tagKeys is the set of every distinct tag key across all records,
	// each of which gets its own tag-encoding-map entry and block.
	tagKeySet map[tagKey]bool
}

// encodeAll encodes every record of the batch into the series buffers.
func (e *recordEncoder) encodeAll(records []*sam.Record) error {
	e.comboIndex = make(map[string]int32)
	e.tagKeySet = make(map[tagKey]bool)
	for i, rec := range records {
		if err := e.encodeRecord(rec); err != nil {
			return fmt.Errorf("cram: encoding record %d (%q): %w", i, rec.QName, err)
		}
	}
	return nil
}

// tagKeysSorted returns every distinct tag key seen across the batch in
// a stable sorted order, so the tag-encoding map and the tag blocks have
// a deterministic layout.
func (e *recordEncoder) tagKeysSorted() []tagKey {
	keys := make([]tagKey, 0, len(e.tagKeySet))
	for k := range e.tagKeySet {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool {
		a, b := keys[i], keys[j]
		for x := 0; x < 3; x++ {
			if a[x] != b[x] {
				return a[x] < b[x]
			}
		}
		return false
	})
	return keys
}

// tagDictionary returns the tag-combination dictionary in TD wire form:
// each combination's three-byte tag keys concatenated, every combination
// terminated by a NUL. It mirrors what parseTagDictionary expects.
func (e *recordEncoder) tagDictionary() []byte {
	var td []byte
	for _, combo := range e.tagCombos {
		for _, k := range combo {
			td = append(td, k[0], k[1], k[2])
		}
		td = append(td, 0)
	}
	if td == nil {
		// A batch with no tags at all still needs one empty entry so a
		// TL value of 0 has a dictionary slot.
		td = []byte{0}
	}
	return td
}

// encodeRecord encodes one record into every relevant data series.
func (e *recordEncoder) encodeRecord(rec *sam.Record) error {
	b := e.buffers
	mapped := rec.Flag&sam.FlagUnmapped == 0

	// CRAM per-record flags: quality is always preserved in its own
	// series, and every record is encoded detached so it carries its own
	// mate fields (the simple writer never uses the downstream-mate
	// optimisation).
	cf := int32(cfQualityPreserved | cfDetached)

	b.bf = e.putU(b.bf, int32(rec.Flag))
	b.cf = e.putU(b.cf, cf)

	refID := recordRefID(rec, e.refIndex)
	if e.multiRef {
		// RI can be -1 (unmapped within a multi-reference slice), so it is a
		// signed series (VARINT_SIGNED in v4).
		b.ri = e.putS(b.ri, refID)
	}

	readLen := len(rec.Seq)
	if rec.Seq == "*" {
		readLen = 0
	}
	b.rl = e.putU(b.rl, int32(readLen))
	e.numBases += int64(readLen)

	// Alignment position is stored absolute (the preservation map's AP
	// entry is false), so no running delta is needed.
	b.ap = e.putU(b.ap, rec.Pos)

	// The read group always travels as an ordinary auxiliary tag, so the
	// RG data series is the no-read-group sentinel -1 for every record. -1
	// makes RG a signed series (VARINT_SIGNED in v4).
	b.rg = e.putS(b.rg, -1)

	// Read names are preserved; an empty / "*" QNAME is stored as the
	// empty string, which round-trips to "".
	name := rec.QName
	if name == "*" {
		name = ""
	}
	b.rn = append(b.rn, name...)
	b.rn = append(b.rn, 0)

	if err := e.encodeMate(rec); err != nil {
		return err
	}

	tl := e.tagLine(rec)
	b.tl = e.putU(b.tl, tl)
	if err := e.encodeTags(rec); err != nil {
		return err
	}

	// Quality is normalised to a fixed-width slice, then mapped through
	// the lossy binning table. normaliseQuality allocates a fresh slice
	// and BinQuality is the identity for BinningNone, so the caller's
	// rec.Qual is never modified and the default writer stays lossless.
	quality := e.binning.BinQuality(normaliseQuality(rec.Qual, readLen))
	if mapped {
		b.mq = e.putU(b.mq, int32(rec.MapQ))
		if err := e.encodeFeatures(rec, readLen); err != nil {
			return err
		}
	} else {
		// An unmapped record stores its bases verbatim in BA.
		b.ba = append(b.ba, seqBytes(rec.Seq)...)
	}
	// Quality is written for every record (cfQualityPreserved is always
	// set): readLen bytes into the QS series.
	b.qs = append(b.qs, quality...)
	return nil
}

// encodeMate writes the detached-mate data series for a record. Every
// record is detached, so MF, NS, NP and TS are always emitted: the mate
// flags packed from the SAM flag bits, the mate reference id, the mate
// position and the template length.
func (e *recordEncoder) encodeMate(rec *sam.Record) error {
	b := e.buffers
	var mf int32
	if rec.Flag&sam.FlagMateReverse != 0 {
		mf |= mfMateReverse
	}
	if rec.Flag&sam.FlagMateUnmapped != 0 {
		mf |= mfMateUnmapped
	}
	b.mf = e.putU(b.mf, mf)

	// The mate reference id: RNEXT "=" means the record's own reference,
	// "*"/"" means none (-1), any other name indexes the @SQ lines.
	var nsID int32
	switch rec.RNext {
	case "", "*":
		nsID = -1
	case "=":
		nsID = recordRefID(rec, e.refIndex)
	default:
		if id, ok := e.refIndex[rec.RNext]; ok {
			nsID = id
		} else {
			nsID = -1
		}
	}
	// NS can be -1 and TS (the template length) can be negative, so both are
	// signed series (VARINT_SIGNED in v4); NP (mate position) is non-negative.
	b.ns = e.putS(b.ns, nsID)
	b.np = e.putU(b.np, rec.PNext)
	b.ts = e.putS(b.ts, rec.TLen)
	return nil
}

// tagLine returns the TL data-series value for a record: the index of
// the record's tag combination in the dictionary, registering a new
// combination the first time it is seen.
func (e *recordEncoder) tagLine(rec *sam.Record) int32 {
	keys := recordTagKeys(rec)
	for _, k := range keys {
		e.tagKeySet[k] = true
	}
	id := comboID(keys)
	if idx, ok := e.comboIndex[id]; ok {
		return idx
	}
	idx := int32(len(e.tagCombos))
	e.tagCombos = append(e.tagCombos, keys)
	e.comboIndex[id] = idx
	return idx
}

// recordTagKeys returns the ordered three-byte tag keys of a record's
// auxiliary fields. The order is the record's own Aux order, which the
// reader reproduces, so tag order round-trips.
func recordTagKeys(rec *sam.Record) []tagKey {
	keys := make([]tagKey, 0, len(rec.Aux))
	for _, a := range rec.Aux {
		keys = append(keys, auxTagKey(a))
	}
	return keys
}

// auxTagKey builds the three-byte CRAM tag key of an auxiliary field:
// the two tag-name bytes and the one-byte CRAM value-type letter.
func auxTagKey(a sam.Aux) tagKey {
	var name [2]byte
	if len(a.Tag) >= 2 {
		name[0], name[1] = a.Tag[0], a.Tag[1]
	}
	return tagKey{name[0], name[1], cramTagType(a)}
}

// comboID returns a string that uniquely identifies an ordered tag-key
// combination, for use as a map key.
func comboID(keys []tagKey) string {
	buf := make([]byte, 0, len(keys)*3)
	for _, k := range keys {
		buf = append(buf, k[0], k[1], k[2])
	}
	return string(buf)
}

// encodeTags writes every auxiliary tag value of a record into its
// per-tag BYTE_ARRAY_LEN series: the ITF-8 value length into the tag's
// length buffer and the value bytes into its value buffer.
func (e *recordEncoder) encodeTags(rec *sam.Record) error {
	for _, a := range rec.Aux {
		key := auxTagKey(a)
		raw, err := encodeTagValue(a)
		if err != nil {
			return fmt.Errorf("tag %q: %w", a.Tag, err)
		}
		e.buffers.tagLens[key] = e.putU(e.buffers.tagLens[key], int32(len(raw)))
		e.buffers.tagVals[key] = append(e.buffers.tagVals[key], raw...)
	}
	return nil
}

// normaliseQuality returns a readLen-byte quality slice for a record. A
// record with no quality ("*", an empty or all-0xff Qual) is stored as
// readLen bytes of 0xff, the SAM "no quality" sentinel, so it round-
// trips back to the same absent state.
func normaliseQuality(qual []byte, readLen int) []byte {
	out := make([]byte, readLen)
	for i := range out {
		out[i] = 0xff
	}
	if len(qual) == readLen {
		copy(out, qual)
	}
	return out
}

// seqBytes returns the bytes of a SEQ string, treating the SAM "*"
// no-sequence marker as empty.
func seqBytes(seq string) []byte {
	if seq == "*" {
		return nil
	}
	return []byte(seq)
}
