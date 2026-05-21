package cram

// appendEncoding serialises one CRAM encoding into dst and returns the
// extended slice. Each encoding is an ITF-8 codec id, an ITF-8
// parameter-byte count, and the codec-specific parameters; it is the
// writer-side inverse of parseEncoding. The simple writer only ever
// emits EXTERNAL and BYTE_ARRAY_LEN encodings.
func appendEncoding(dst []byte, enc encSpec) []byte {
	params := enc.params()
	dst = appendITF8(dst, int32(enc.id))
	dst = appendITF8(dst, int32(len(params)))
	return append(dst, params...)
}

// encSpec is a minimal description of one CRAM encoding the writer
// emits: an EXTERNAL block reference or a BYTE_ARRAY_LEN pair of
// sub-encodings. It is deliberately narrower than the full Encoding
// struct — the writer never produces the CORE-bitstream codecs.
type encSpec struct {
	id EncodingID
	// externalID is the content id for an EXTERNAL or BYTE_ARRAY_STOP
	// encoding.
	externalID int32
	// stopByte is the delimiter of a BYTE_ARRAY_STOP encoding.
	stopByte byte
	// lenEnc and valEnc are the sub-encodings of a BYTE_ARRAY_LEN; both
	// nil for an EXTERNAL encoding.
	lenEnc *encSpec
	valEnc *encSpec
}

// extSpec returns an EXTERNAL encoding drawing from the given content id.
func extSpec(id int32) encSpec {
	return encSpec{id: EncodingExternal, externalID: id}
}

// byteArrayLenSpec returns a BYTE_ARRAY_LEN encoding whose length and
// value sub-encodings are EXTERNAL blocks.
func byteArrayLenSpec(lenID, valID int32) encSpec {
	l := extSpec(lenID)
	v := extSpec(valID)
	return encSpec{id: EncodingByteArrayLen, lenEnc: &l, valEnc: &v}
}

// byteArrayStopSpec returns a BYTE_ARRAY_STOP encoding: values are read
// from an external block up to a stop byte.
func byteArrayStopSpec(stop byte, id int32) encSpec {
	return encSpec{id: EncodingByteArrayStop, externalID: id, stopByte: stop}
}

// params serialises the codec-specific parameter bytes of the encoding.
func (e encSpec) params() []byte {
	switch e.id {
	case EncodingExternal:
		return appendITF8(nil, e.externalID)
	case EncodingByteArrayStop:
		out := []byte{e.stopByte}
		return appendITF8(out, e.externalID)
	case EncodingByteArrayLen:
		var out []byte
		out = appendEncoding(out, *e.lenEnc)
		out = appendEncoding(out, *e.valEnc)
		return out
	default:
		return nil
	}
}

// encodeCompressionHeader serialises a CRAM compression-header block
// payload: the preservation map, the data-series encoding map and the
// tag-encoding map. It is the writer-side inverse of
// parseCompressionHeader.
//
// multiRef reports whether the slice is multi-reference (and so needs an
// RI data series). tagDict is the TD preservation-map entry; tagKeys
// lists every distinct auxiliary tag, in the order their content ids
// were assigned.
func encodeCompressionHeader(multiRef bool, tagDict []byte, tagKeys []tagKey) []byte {
	var out []byte
	out = append(out, encodePreservationMap(tagDict)...)
	out = append(out, encodeDataSeriesMap(multiRef, tagKeys)...)
	out = append(out, encodeTagEncodingMap(tagKeys)...)
	return out
}

// encodePreservationMap serialises the preservation map. The writer
// preserves read names (RN true), stores alignment positions absolutely
// (AP false) and produces reference-free CRAM (RR false); it emits the
// default substitution matrix (SM) and the tag-combination dictionary
// (TD).
func encodePreservationMap(tagDict []byte) []byte {
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
	// TD: an ITF-8 byte count then the dictionary bytes.
	entries = append(entries, 'T', 'D')
	entries = appendITF8(entries, int32(len(tagDict)))
	entries = append(entries, tagDict...)
	count++

	return frameMap(entries, count)
}

// dataSeriesEntry pairs a two-character data-series key with its
// encoding for the data-series encoding map.
type dataSeriesEntry struct {
	key dataSeriesKey
	enc encSpec
}

// encodeDataSeriesMap serialises the data-series encoding map: one entry
// per data series the writer's records use, each an EXTERNAL or
// BYTE_ARRAY_LEN encoding pointing at the series' block. Every series is
// declared even if a particular slice leaves it empty — an absent block
// is harmless, but a missing encoding-map entry is not.
func encodeDataSeriesMap(multiRef bool, tagKeys []tagKey) []byte {
	entries := []dataSeriesEntry{
		{dataSeriesKey{'B', 'F'}, extSpec(cidBF)},
		{dataSeriesKey{'C', 'F'}, extSpec(cidCF)},
		{dataSeriesKey{'R', 'L'}, extSpec(cidRL)},
		{dataSeriesKey{'A', 'P'}, extSpec(cidAP)},
		{dataSeriesKey{'R', 'G'}, extSpec(cidRG)},
		{dataSeriesKey{'R', 'N'}, byteArrayStopSpec(0, cidRN)},
		{dataSeriesKey{'M', 'F'}, extSpec(cidMF)},
		{dataSeriesKey{'N', 'S'}, extSpec(cidNS)},
		{dataSeriesKey{'N', 'P'}, extSpec(cidNP)},
		{dataSeriesKey{'T', 'S'}, extSpec(cidTS)},
		{dataSeriesKey{'N', 'F'}, extSpec(cidMF)}, // unused; declared for completeness.
		{dataSeriesKey{'T', 'L'}, extSpec(cidTL)},
		{dataSeriesKey{'M', 'Q'}, extSpec(cidMQ)},
		{dataSeriesKey{'F', 'N'}, extSpec(cidFN)},
		{dataSeriesKey{'F', 'C'}, extSpec(cidFC)},
		{dataSeriesKey{'F', 'P'}, extSpec(cidFP)},
		{dataSeriesKey{'B', 'B'}, byteArrayLenSpec(cidBBLen, cidBB)},
		{dataSeriesKey{'I', 'N'}, byteArrayLenSpec(cidINLen, cidIN)},
		{dataSeriesKey{'S', 'C'}, byteArrayLenSpec(cidSCLen, cidSC)},
		{dataSeriesKey{'D', 'L'}, extSpec(cidDL)},
		{dataSeriesKey{'R', 'S'}, extSpec(cidRS)},
		{dataSeriesKey{'P', 'D'}, extSpec(cidPD)},
		{dataSeriesKey{'H', 'C'}, extSpec(cidHC)},
		{dataSeriesKey{'B', 'A'}, extSpec(cidBA)},
		{dataSeriesKey{'Q', 'S'}, extSpec(cidQS)},
	}
	if multiRef {
		entries = append(entries, dataSeriesEntry{dataSeriesKey{'R', 'I'}, extSpec(cidRI)})
	}

	var body []byte
	for _, e := range entries {
		body = append(body, e.key[0], e.key[1])
		body = appendEncoding(body, e.enc)
	}
	return frameMap(body, int32(len(entries)))
}

// encodeTagEncodingMap serialises the tag-encoding map: one entry per
// distinct auxiliary tag key, each keyed by the ITF-8 integer whose
// three low bytes are the tag key and each carrying a BYTE_ARRAY_LEN
// encoding over a length block and a value block. The content ids
// assigned to a tag are tagContentIDs of its index in tagKeys.
func encodeTagEncodingMap(tagKeys []tagKey) []byte {
	var body []byte
	for i, key := range tagKeys {
		id := int32(uint32(key[0])<<16 | uint32(key[1])<<8 | uint32(key[2]))
		body = appendITF8(body, id)
		lenID, valID := tagContentIDs(i)
		body = appendEncoding(body, byteArrayLenSpec(lenID, valID))
	}
	return frameMap(body, int32(len(tagKeys)))
}

// frameMap prefixes a CRAM compression-header map's entry bytes with the
// ITF-8 byte size and ITF-8 entry count the parser reads via mapFrame.
// The size field measures from the entry count onward.
func frameMap(entries []byte, count int32) []byte {
	var counted []byte
	counted = appendITF8(counted, count)
	counted = append(counted, entries...)
	var out []byte
	out = appendITF8(out, int32(len(counted)))
	return append(out, counted...)
}
