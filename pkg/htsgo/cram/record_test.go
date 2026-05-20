package cram

import (
	"bytes"
	"testing"
)

// extEnc returns an EXTERNAL encoding drawing from the given content id.
func extEnc(id int32) *Encoding { return &Encoding{ID: EncodingExternal, ExternalID: id} }

// stopEnc returns a BYTE_ARRAY_STOP encoding with a NUL delimiter.
func stopEnc(id int32) *Encoding {
	return &Encoding{ID: EncodingByteArrayStop, StopByte: 0, ExternalID: id}
}

// TestRecordSeriesHelpers exercises the per-series read helpers of the
// record decoder over hand-built external blocks.
func TestRecordSeriesHelpers(t *testing.T) {
	h := refFreeHeader()
	h.DataSeries[dataSeriesKey{'B', 'F'}] = extEnc(1)
	h.DataSeries[dataSeriesKey{'B', 'A'}] = extEnc(2)
	h.DataSeries[dataSeriesKey{'R', 'N'}] = stopEnc(3)

	var bf bytes.Buffer
	bf.Write(encITF8(99))
	src := newTestSource(nil, map[int32][]byte{
		1: bf.Bytes(),
		2: {'A'},
		3: append([]byte("name"), 0),
	})
	rd := &recordDecoder{h: h, src: &SeriesSource{s: src}, slice: &SliceHeader{}}

	if v, err := rd.intSeries("BF"); err != nil || v != 99 {
		t.Errorf("intSeries(BF) = %d, %v; want 99", v, err)
	}
	if b, err := rd.byteSeries("BA"); err != nil || b != 'A' {
		t.Errorf("byteSeries(BA) = %c, %v; want A", b, err)
	}
	if ba, err := rd.byteArraySeries("RN"); err != nil || string(ba) != "name" {
		t.Errorf("byteArraySeries(RN) = %q, %v; want name", ba, err)
	}

	// A series absent from the encoding map is an error, not a panic.
	if _, err := rd.intSeries("ZZ"); err == nil {
		t.Error("intSeries on an unknown series should error")
	}
	if _, err := rd.byteSeries("ZZ"); err == nil {
		t.Error("byteSeries on an unknown series should error")
	}
	if _, err := rd.byteArraySeries("ZZ"); err == nil {
		t.Error("byteArraySeries on an unknown series should error")
	}
}

// TestRefName checks reference-id resolution including the unmapped
// sentinel and an out-of-range id.
func TestRefName(t *testing.T) {
	rd := &recordDecoder{refNames: []string{"chr1", "chr2"}}
	if n, err := rd.refName(0); err != nil || n != "chr1" {
		t.Errorf("refName(0) = %q, %v; want chr1", n, err)
	}
	if n, err := rd.refName(-1); err != nil || n != "" {
		t.Errorf("refName(-1) = %q, %v; want \"\"", n, err)
	}
	if _, err := rd.refName(99); err == nil {
		t.Error("refName on an out-of-range id should error")
	}
}

// TestDecodeFeaturesPayloads drives decodeFeatures over a synthetic FC /
// FP / payload layout covering several feature codes.
func TestDecodeFeaturesPayloads(t *testing.T) {
	h := refFreeHeader()
	h.DataSeries[dataSeriesKey{'F', 'C'}] = extEnc(10) // feature codes
	h.DataSeries[dataSeriesKey{'F', 'P'}] = extEnc(11) // feature positions
	h.DataSeries[dataSeriesKey{'B', 'B'}] = &Encoding{
		ID: EncodingByteArrayLen, LenEnc: extEnc(12), ValEnc: extEnc(13),
	}
	h.DataSeries[dataSeriesKey{'I', 'N'}] = stopEnc(14)
	h.DataSeries[dataSeriesKey{'D', 'L'}] = extEnc(15)
	h.DataSeries[dataSeriesKey{'B', 'A'}] = extEnc(16)
	h.DataSeries[dataSeriesKey{'Q', 'S'}] = extEnc(17)

	// Three features: a base stretch "AC", an insertion "GG", a deletion.
	var fc, fp, bbLen, bbVal, inBlk, dl bytes.Buffer
	fc.WriteString("bID")
	fp.Write(encITF8(1)) // b at pos 1
	fp.Write(encITF8(2)) // I at pos 3 (delta 2)
	fp.Write(encITF8(0)) // D at pos 3 (delta 0)
	bbLen.Write(encITF8(2))
	bbVal.WriteString("AC")
	inBlk.WriteString("GG")
	inBlk.WriteByte(0)
	dl.Write(encITF8(4))

	src := newTestSource(nil, map[int32][]byte{
		10: fc.Bytes(), 11: fp.Bytes(),
		12: bbLen.Bytes(), 13: bbVal.Bytes(),
		14: inBlk.Bytes(), 15: dl.Bytes(),
	})
	rd := &recordDecoder{h: h, src: &SeriesSource{s: src}, slice: &SliceHeader{}}
	feats, err := rd.decodeFeatures(3)
	if err != nil {
		t.Fatalf("decodeFeatures: %v", err)
	}
	if len(feats) != 3 {
		t.Fatalf("decoded %d features, want 3", len(feats))
	}
	if feats[0].code != featBases || string(feats[0].bases) != "AC" {
		t.Errorf("feature 0 = %+v, want b/AC", feats[0])
	}
	if feats[1].code != featInsertion || string(feats[1].bases) != "GG" {
		t.Errorf("feature 1 = %+v, want I/GG", feats[1])
	}
	if feats[2].code != featDeletion || feats[2].length != 4 {
		t.Errorf("feature 2 = %+v, want D/length 4", feats[2])
	}
	// The positions accumulate from the FP deltas: 1, 3, 3.
	if feats[0].pos != 1 || feats[1].pos != 3 || feats[2].pos != 3 {
		t.Errorf("feature positions = %d/%d/%d, want 1/3/3",
			feats[0].pos, feats[1].pos, feats[2].pos)
	}

	// A negative feature count is rejected without a panic.
	if _, err := rd.decodeFeatures(-1); err == nil {
		t.Error("decodeFeatures(-1) should error")
	}
}

// TestDecodeFeatureUnknownCode checks that an unrecognised feature code
// surfaces as an error.
func TestDecodeFeatureUnknownCode(t *testing.T) {
	rd := &recordDecoder{h: refFreeHeader(), src: &SeriesSource{s: newTestSource(nil, nil)}}
	f := readFeature{code: 'Z'}
	if err := rd.decodeFeaturePayload(&f); err == nil {
		t.Error("an unknown feature code should error")
	}
}

// TestDecodeFeaturePayloadAllCodes drives decodeFeaturePayload over the
// single-payload feature codes (B, i, q, Q, X, S, N, P, H) so each
// branch of the dispatch is exercised.
func TestDecodeFeaturePayloadAllCodes(t *testing.T) {
	h := refFreeHeader()
	h.DataSeries[dataSeriesKey{'B', 'A'}] = extEnc(1)
	h.DataSeries[dataSeriesKey{'Q', 'S'}] = extEnc(2)
	h.DataSeries[dataSeriesKey{'Q', 'Q'}] = &Encoding{
		ID: EncodingByteArrayLen, LenEnc: extEnc(3), ValEnc: extEnc(4),
	}
	h.DataSeries[dataSeriesKey{'B', 'S'}] = extEnc(5)
	h.DataSeries[dataSeriesKey{'S', 'C'}] = stopEnc(6)
	h.DataSeries[dataSeriesKey{'R', 'S'}] = extEnc(7)
	h.DataSeries[dataSeriesKey{'P', 'D'}] = extEnc(8)
	h.DataSeries[dataSeriesKey{'H', 'C'}] = extEnc(9)

	var qqLen, qqVal, rs, pd, hc bytes.Buffer
	qqLen.Write(encITF8(2))
	qqVal.Write([]byte{7, 8})
	rs.Write(encITF8(12))
	pd.Write(encITF8(3))
	hc.Write(encITF8(5))
	src := newTestSource(nil, map[int32][]byte{
		1: {'G', 'C'},    // BA: two single-base values.
		2: {40, 41},      // QS: two quality scores.
		3: qqLen.Bytes(), // QQ length.
		4: qqVal.Bytes(), // QQ values.
		5: {2},           // BS: one substitution code.
		6: append([]byte("SOFT"), 0),
		7: rs.Bytes(),
		8: pd.Bytes(),
		9: hc.Bytes(),
	})
	rd := &recordDecoder{h: h, src: &SeriesSource{s: src}, slice: &SliceHeader{}}

	check := func(code byte, verify func(readFeature) bool) {
		f := readFeature{code: code}
		if err := rd.decodeFeaturePayload(&f); err != nil {
			t.Errorf("feature %c: %v", code, err)
			return
		}
		if !verify(f) {
			t.Errorf("feature %c: payload not as expected: %+v", code, f)
		}
	}
	check(featBase, func(f readFeature) bool { return f.base == 'G' && f.quality == 40 })
	check(featInsertBase, func(f readFeature) bool { return f.base == 'C' })
	check(featQualityScore, func(f readFeature) bool { return f.quality == 41 })
	check(featScores, func(f readFeature) bool { return bytes.Equal(f.bases, []byte{7, 8}) })
	check(featSubst, func(f readFeature) bool { return f.substCode == 2 })
	check(featSoftClip, func(f readFeature) bool { return string(f.bases) == "SOFT" })
	check(featRefSkip, func(f readFeature) bool { return f.length == 12 })
	check(featPadding, func(f readFeature) bool { return f.length == 3 })
	check(featHardClip, func(f readFeature) bool { return f.length == 5 })
}

// TestDecodeFeaturesErrorWraps drives the error paths of decodeFeatures:
// an exhausted feature-code block, an exhausted feature-position block,
// and an exhausted payload block. Each must surface as an error, never
// a panic.
func TestDecodeFeaturesErrorWraps(t *testing.T) {
	// Exhausted FC block: zero bytes for a one-feature record.
	h := refFreeHeader()
	h.DataSeries[dataSeriesKey{'F', 'C'}] = extEnc(1)
	h.DataSeries[dataSeriesKey{'F', 'P'}] = extEnc(2)
	rd := &recordDecoder{
		h:   h,
		src: &SeriesSource{s: newTestSource(nil, map[int32][]byte{1: {}, 2: {}})},
	}
	if _, err := rd.decodeFeatures(1); err == nil {
		t.Error("an exhausted feature-code block should error")
	}

	// FC present but FP exhausted.
	h2 := refFreeHeader()
	h2.DataSeries[dataSeriesKey{'F', 'C'}] = extEnc(1)
	h2.DataSeries[dataSeriesKey{'F', 'P'}] = extEnc(2)
	rd2 := &recordDecoder{
		h:   h2,
		src: &SeriesSource{s: newTestSource(nil, map[int32][]byte{1: {featBases}, 2: {}})},
	}
	if _, err := rd2.decodeFeatures(1); err == nil {
		t.Error("an exhausted feature-position block should error")
	}

	// FC and FP present but the BB payload block is exhausted.
	h3 := refFreeHeader()
	h3.DataSeries[dataSeriesKey{'F', 'C'}] = extEnc(1)
	h3.DataSeries[dataSeriesKey{'F', 'P'}] = extEnc(2)
	h3.DataSeries[dataSeriesKey{'B', 'B'}] = &Encoding{
		ID: EncodingByteArrayLen, LenEnc: extEnc(3), ValEnc: extEnc(4),
	}
	rd3 := &recordDecoder{
		h: h3,
		src: &SeriesSource{s: newTestSource(nil, map[int32][]byte{
			1: {featBases}, 2: encITF8(1), 3: {}, 4: {},
		})},
	}
	if _, err := rd3.decodeFeatures(1); err == nil {
		t.Error("an exhausted feature-payload block should error")
	}

	// A feature payload whose data series is absent from the map errors.
	h4 := refFreeHeader()
	h4.DataSeries[dataSeriesKey{'F', 'C'}] = extEnc(1)
	h4.DataSeries[dataSeriesKey{'F', 'P'}] = extEnc(2)
	rd4 := &recordDecoder{
		h: h4,
		src: &SeriesSource{s: newTestSource(nil, map[int32][]byte{
			1: {featInsertion}, 2: encITF8(1),
		})},
	}
	if _, err := rd4.decodeFeatures(1); err == nil {
		t.Error("a feature whose payload series is absent should error")
	}
}

// TestDecodeRecordMultiRef checks the per-record reference-id path used
// by a multi-reference slice (slice ref id -2).
func TestDecodeRecordMultiRef(t *testing.T) {
	h := refFreeHeader()
	h.Preservation.ReadNamesIncluded = true
	h.Preservation.APDelta = false
	h.Preservation.TagDictionary = []byte{0}
	h.DataSeries[dataSeriesKey{'B', 'F'}] = huffConst(int32(0))
	h.DataSeries[dataSeriesKey{'C', 'F'}] = huffConst(0)
	h.DataSeries[dataSeriesKey{'R', 'I'}] = extEnc(1)
	h.DataSeries[dataSeriesKey{'R', 'L'}] = huffConst(3)
	h.DataSeries[dataSeriesKey{'A', 'P'}] = huffConst(5)
	h.DataSeries[dataSeriesKey{'R', 'G'}] = huffConst(-1)
	h.DataSeries[dataSeriesKey{'R', 'N'}] = stopEnc(2)
	h.DataSeries[dataSeriesKey{'T', 'L'}] = huffConst(0)
	h.DataSeries[dataSeriesKey{'M', 'Q'}] = huffConst(0)
	h.DataSeries[dataSeriesKey{'F', 'N'}] = huffConst(1)
	h.DataSeries[dataSeriesKey{'F', 'C'}] = huffConst(int32(featBases))
	h.DataSeries[dataSeriesKey{'F', 'P'}] = huffConst(1)
	h.DataSeries[dataSeriesKey{'B', 'B'}] = &Encoding{
		ID: EncodingByteArrayLen, LenEnc: huffConst(3), ValEnc: extEnc(3),
	}
	var ri, rn bytes.Buffer
	ri.Write(encITF8(1)) // reference id 1.
	rn.WriteString("m")
	rn.WriteByte(0)
	src := newTestSource(nil, map[int32][]byte{
		1: ri.Bytes(), 2: rn.Bytes(), 3: []byte("ACG"),
	})
	rd := &recordDecoder{
		h: h, src: &SeriesSource{s: src},
		slice:        &SliceHeader{RefSeqID: -2},
		refNames:     []string{"chr1", "chr2"},
		readLenLimit: 1 << 20,
	}
	dr, err := rd.decodeRecord(0)
	if err != nil {
		t.Fatalf("decodeRecord: %v", err)
	}
	if dr.rec.RName != "chr2" {
		t.Errorf("multi-ref record RName = %q, want chr2", dr.rec.RName)
	}
}
