package cram

// appendEncoding serialises one CRAM encoding into dst and returns the
// extended slice. Each encoding is a codec id, a parameter-byte count, and
// the codec-specific parameters, all framed through iw (ITF-8 for v2/v3,
// uint7 varints for v4). It is the writer-side inverse of parseEncoding.
// The simple writer emits EXTERNAL, BYTE_ARRAY_STOP and BYTE_ARRAY_LEN for
// every version, and additionally VARINT_UNSIGNED / VARINT_SIGNED for the
// integer series of a CRAM v4 file.
func appendEncoding(iw intWriter, dst []byte, enc encSpec) []byte {
	params := enc.params(iw)
	dst = iw.u32(dst, int32(enc.id))
	dst = iw.u32(dst, int32(len(params)))
	return append(dst, params...)
}

// encSpec is a minimal description of one CRAM encoding the writer emits:
// an EXTERNAL block reference, a VARINT integer codec, a BYTE_ARRAY_STOP
// stream, or a BYTE_ARRAY_LEN pair of sub-encodings. It is deliberately
// narrower than the full Encoding struct — the writer never produces the
// CORE-bitstream codecs.
type encSpec struct {
	id EncodingID
	// externalID is the content id for an EXTERNAL, BYTE_ARRAY_STOP or
	// VARINT_UNSIGNED / VARINT_SIGNED encoding.
	externalID int32
	// varintOffset is the constant offset of a VARINT_UNSIGNED /
	// VARINT_SIGNED encoding. The writer always emits 0 — it stores raw
	// values — but the field exists so the serialisation mirrors the
	// decoder's two-parameter read.
	varintOffset int64
	// stopByte is the delimiter of a BYTE_ARRAY_STOP encoding.
	stopByte byte
	// lenEnc and valEnc are the sub-encodings of a BYTE_ARRAY_LEN; both
	// nil for the leaf encodings.
	lenEnc *encSpec
	valEnc *encSpec
}

// extSpec returns an EXTERNAL encoding drawing from the given content id.
func extSpec(id int32) encSpec {
	return encSpec{id: EncodingExternal, externalID: id}
}

// varintUnsignedSpec returns a CRAM v4 VARINT_UNSIGNED encoding over the
// given content id with a zero offset: an unsigned uint7 varint stream in
// its own external block. It is the v4 encoding for every non-negative
// integer data series.
func varintUnsignedSpec(id int32) encSpec {
	return encSpec{id: EncodingVarintUnsigned, externalID: id}
}

// varintSignedSpec returns a CRAM v4 VARINT_SIGNED encoding over the given
// content id with a zero offset: a zig-zag signed uint7 varint stream. It
// is the v4 encoding for an integer data series that can carry a negative
// value (RI, RG, NS, TS).
func varintSignedSpec(id int32) encSpec {
	return encSpec{id: EncodingVarintSigned, externalID: id}
}

// intSeriesSpec returns the encoding for an integer data series stored in
// the external block id. For CRAM v2/v3 it is EXTERNAL (the integers are
// ITF-8 in the block); for CRAM v4 it is VARINT_UNSIGNED or, when signed is
// true, VARINT_SIGNED (the integers are uint7 varints, EXTERNAL being
// byte-only in v4). This is the single place the writer maps an integer
// series onto the version-correct codec.
func intSeriesSpec(version Version, id int32, signed bool) encSpec {
	if !version.usesUint7() {
		return extSpec(id)
	}
	if signed {
		return varintSignedSpec(id)
	}
	return varintUnsignedSpec(id)
}

// byteArrayLenSpec returns a BYTE_ARRAY_LEN encoding whose length and value
// sub-encodings draw from the given content ids. The value sub-encoding is
// always EXTERNAL (raw bytes); the length sub-encoding is EXTERNAL for
// v2/v3 and VARINT_UNSIGNED for v4 (where the length integers are uint7
// varints), matching what the decoder's LenEnc.decodeInt expects.
func byteArrayLenSpec(version Version, lenID, valID int32) encSpec {
	l := intSeriesSpec(version, lenID, false)
	v := extSpec(valID)
	return encSpec{id: EncodingByteArrayLen, lenEnc: &l, valEnc: &v}
}

// byteArrayStopSpec returns a BYTE_ARRAY_STOP encoding: values are read
// from an external block up to a stop byte.
func byteArrayStopSpec(stop byte, id int32) encSpec {
	return encSpec{id: EncodingByteArrayStop, externalID: id, stopByte: stop}
}

// params serialises the codec-specific parameter bytes of the encoding,
// framing every integer through iw (ITF-8 for v2/v3, uint7 for v4).
func (e encSpec) params(iw intWriter) []byte {
	switch e.id {
	case EncodingExternal:
		return iw.u32(nil, e.externalID)
	case EncodingVarintUnsigned, EncodingVarintSigned:
		// CRAM v4: a content id then a signed 64-bit offset
		// (cram_codecs.c cram_varint_decode_init).
		out := iw.u32(nil, e.externalID)
		return iw.s64(out, e.varintOffset)
	case EncodingByteArrayStop:
		out := []byte{e.stopByte}
		return iw.u32(out, e.externalID)
	case EncodingByteArrayLen:
		var out []byte
		out = appendEncoding(iw, out, *e.lenEnc)
		out = appendEncoding(iw, out, *e.valEnc)
		return out
	default:
		return nil
	}
}

// encodeCompressionHeader serialises a CRAM compression-header block
// payload: the preservation map, the data-series encoding map and the
// tag-encoding map. It is the writer-side inverse of
// parseCompressionHeader. version selects the integer framing (ITF-8 for
// v2/v3, uint7 for v4) and the integer data-series codec (EXTERNAL vs
// VARINT).
//
// multiRef reports whether the slice is multi-reference (and so needs an
// RI data series). tagDict is the TD preservation-map entry; tagKeys
// lists every distinct auxiliary tag, in the order their content ids
// were assigned.
func encodeCompressionHeader(version Version, multiRef bool, tagDict []byte, tagKeys []tagKey) []byte {
	iw := newIntWriter(version)
	var out []byte
	out = append(out, encodePreservationMap(iw, tagDict)...)
	out = append(out, encodeDataSeriesMap(version, iw, multiRef, tagKeys)...)
	out = append(out, encodeTagEncodingMap(version, iw, tagKeys)...)
	return out
}

// encodePreservationMap serialises the preservation map. The writer
// preserves read names (RN true), stores alignment positions absolutely
// (AP false) and produces reference-free CRAM (RR false); it emits the
// default substitution matrix (SM) and the tag-combination dictionary
// (TD). The TD byte count and the map frame are integer-framed through iw.
func encodePreservationMap(iw intWriter, tagDict []byte) []byte {
	var entries []byte
	count := int32(0)

	entries = append(entries, 'R', 'N', 1) // read names included.
	count++
	entries = append(entries, 'A', 'P', 0) // alignment positions absolute.
	count++
	entries = append(entries, 'R', 'R', 0) // reference not required.
	count++
	// SM: the default substitution matrix (codes map straight to
	// candidates), five bytes of 0x1B.
	entries = append(entries, 'S', 'M', 0x1b, 0x1b, 0x1b, 0x1b, 0x1b)
	count++
	// TD: a varint byte count then the dictionary bytes.
	entries = append(entries, 'T', 'D')
	entries = iw.u32(entries, int32(len(tagDict)))
	entries = append(entries, tagDict...)
	count++

	return frameMap(iw, entries, count)
}

// dataSeriesEntry pairs a two-character data-series key with its
// encoding for the data-series encoding map.
type dataSeriesEntry struct {
	key dataSeriesKey
	enc encSpec
}

// encodeDataSeriesMap serialises the data-series encoding map: one entry
// per data series the writer's records use. Every series is declared even
// if a particular slice leaves it empty — an absent block is harmless, but
// a missing encoding-map entry is not.
//
// The integer series use intSeriesSpec, which is EXTERNAL for v2/v3 and the
// matching VARINT codec for v4 (signed for the series that can be negative —
// RG, NS, TS, RI — unsigned otherwise). The byte series (RN, BB/IN/SC
// values, BA, QS, FC) stay EXTERNAL / BYTE_ARRAY_* in every version.
func encodeDataSeriesMap(version Version, iw intWriter, multiRef bool, tagKeys []tagKey) []byte {
	entries := []dataSeriesEntry{
		{dataSeriesKey{'B', 'F'}, intSeriesSpec(version, cidBF, false)},
		{dataSeriesKey{'C', 'F'}, intSeriesSpec(version, cidCF, false)},
		{dataSeriesKey{'R', 'L'}, intSeriesSpec(version, cidRL, false)},
		{dataSeriesKey{'A', 'P'}, intSeriesSpec(version, cidAP, false)},
		{dataSeriesKey{'R', 'G'}, intSeriesSpec(version, cidRG, true)},
		{dataSeriesKey{'R', 'N'}, byteArrayStopSpec(0, cidRN)},
		{dataSeriesKey{'M', 'F'}, intSeriesSpec(version, cidMF, false)},
		{dataSeriesKey{'N', 'S'}, intSeriesSpec(version, cidNS, true)},
		{dataSeriesKey{'N', 'P'}, intSeriesSpec(version, cidNP, false)},
		{dataSeriesKey{'T', 'S'}, intSeriesSpec(version, cidTS, true)},
		{dataSeriesKey{'N', 'F'}, intSeriesSpec(version, cidMF, false)}, // unused; declared for completeness.
		{dataSeriesKey{'T', 'L'}, intSeriesSpec(version, cidTL, false)},
		{dataSeriesKey{'M', 'Q'}, intSeriesSpec(version, cidMQ, false)},
		{dataSeriesKey{'F', 'N'}, intSeriesSpec(version, cidFN, false)},
		{dataSeriesKey{'F', 'C'}, extSpec(cidFC)},
		{dataSeriesKey{'F', 'P'}, intSeriesSpec(version, cidFP, false)},
		{dataSeriesKey{'B', 'B'}, byteArrayLenSpec(version, cidBBLen, cidBB)},
		{dataSeriesKey{'I', 'N'}, byteArrayLenSpec(version, cidINLen, cidIN)},
		{dataSeriesKey{'S', 'C'}, byteArrayLenSpec(version, cidSCLen, cidSC)},
		{dataSeriesKey{'D', 'L'}, intSeriesSpec(version, cidDL, false)},
		{dataSeriesKey{'R', 'S'}, intSeriesSpec(version, cidRS, false)},
		{dataSeriesKey{'P', 'D'}, intSeriesSpec(version, cidPD, false)},
		{dataSeriesKey{'H', 'C'}, intSeriesSpec(version, cidHC, false)},
		{dataSeriesKey{'B', 'A'}, extSpec(cidBA)},
		{dataSeriesKey{'Q', 'S'}, extSpec(cidQS)},
	}
	if multiRef {
		entries = append(entries, dataSeriesEntry{dataSeriesKey{'R', 'I'}, intSeriesSpec(version, cidRI, true)})
	}

	var body []byte
	for _, e := range entries {
		body = append(body, e.key[0], e.key[1])
		body = appendEncoding(iw, body, e.enc)
	}
	return frameMap(iw, body, int32(len(entries)))
}

// encodeTagEncodingMap serialises the tag-encoding map: one entry per
// distinct auxiliary tag key, each keyed by the integer whose three low
// bytes are the tag key and each carrying a BYTE_ARRAY_LEN encoding over a
// length block and a value block. The content ids assigned to a tag are
// tagContentIDs of its index in tagKeys. The key and length integers are
// framed through iw (ITF-8 for v2/v3, uint7 for v4); the length
// sub-encoding is VARINT_UNSIGNED for v4 so the decoder reads the uint7
// length the writer emits.
func encodeTagEncodingMap(version Version, iw intWriter, tagKeys []tagKey) []byte {
	var body []byte
	for i, key := range tagKeys {
		id := int32(uint32(key[0])<<16 | uint32(key[1])<<8 | uint32(key[2]))
		body = iw.u32(body, id)
		lenID, valID := tagContentIDs(i)
		body = appendEncoding(iw, body, byteArrayLenSpec(version, lenID, valID))
	}
	return frameMap(iw, body, int32(len(tagKeys)))
}

// frameMap prefixes a CRAM compression-header map's entry bytes with the
// byte size and entry count the parser reads via mapFrame, both framed
// through iw (ITF-8 for v2/v3, uint7 for v4). The size field measures from
// the entry count onward.
func frameMap(iw intWriter, entries []byte, count int32) []byte {
	var counted []byte
	counted = iw.u32(counted, count)
	counted = append(counted, entries...)
	var out []byte
	out = iw.u32(out, int32(len(counted)))
	return append(out, counted...)
}
