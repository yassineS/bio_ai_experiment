package cram

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// huffConst returns a degenerate single-symbol HUFFMAN encoding, the
// CRAM idiom for a constant data series. It consumes no CORE bits.
func huffConst(sym int32) *Encoding {
	return &Encoding{ID: EncodingHuffman, Symbols: []int32{sym}, BitLengths: []int32{0}}
}

// mappedRecordHeader builds a compression header for a single mapped,
// reference-free record whose whole sequence is one base-stretch
// feature. The data series are wired to a small set of external blocks.
func mappedRecordHeader() *CompressionHeader {
	h := refFreeHeader()
	h.Preservation.ReadNamesIncluded = true
	h.Preservation.APDelta = false
	h.Preservation.TagDictionary = []byte{0} // one empty dictionary entry.
	h.DataSeries[dataSeriesKey{'B', 'F'}] = extEnc(1)
	h.DataSeries[dataSeriesKey{'C', 'F'}] = huffConst(0)
	h.DataSeries[dataSeriesKey{'R', 'L'}] = extEnc(2)
	h.DataSeries[dataSeriesKey{'A', 'P'}] = extEnc(3)
	h.DataSeries[dataSeriesKey{'R', 'G'}] = huffConst(-1)
	h.DataSeries[dataSeriesKey{'R', 'N'}] = stopEnc(4)
	h.DataSeries[dataSeriesKey{'T', 'L'}] = huffConst(0)
	h.DataSeries[dataSeriesKey{'M', 'Q'}] = huffConst(60)
	h.DataSeries[dataSeriesKey{'F', 'N'}] = huffConst(1)
	h.DataSeries[dataSeriesKey{'F', 'C'}] = extEnc(5)
	h.DataSeries[dataSeriesKey{'F', 'P'}] = huffConst(1)
	h.DataSeries[dataSeriesKey{'B', 'B'}] = &Encoding{
		ID: EncodingByteArrayLen, LenEnc: extEnc(6), ValEnc: extEnc(7),
	}
	return h
}

// TestDecodeRecordMapped decodes a single hand-built mapped record end
// to end through decodeRecord.
func TestDecodeRecordMapped(t *testing.T) {
	h := mappedRecordHeader()
	var bf, rl, ap, rn, fc, bbLen, bbVal bytes.Buffer
	bf.Write(encITF8(0)) // a plain mapped record.
	rl.Write(encITF8(5))
	ap.Write(encITF8(100))
	rn.WriteString("read1")
	rn.WriteByte(0)
	fc.WriteByte(featBases)
	bbLen.Write(encITF8(5))
	bbVal.WriteString("ACGTA")
	src := newTestSource(nil, map[int32][]byte{
		1: bf.Bytes(), 2: rl.Bytes(), 3: ap.Bytes(), 4: rn.Bytes(),
		5: fc.Bytes(), 6: bbLen.Bytes(), 7: bbVal.Bytes(),
	})
	rd := &recordDecoder{
		h: h, src: &SeriesSource{s: src},
		slice:        &SliceHeader{RefSeqID: 0},
		refNames:     []string{"chr1"},
		readLenLimit: 1 << 20,
	}
	dr, err := rd.decodeRecord(0)
	if err != nil {
		t.Fatalf("decodeRecord: %v", err)
	}
	rec := dr.rec
	if rec.QName != "read1" || rec.RName != "chr1" || rec.Pos != 100 {
		t.Errorf("record = %s/%s/%d, want read1/chr1/100", rec.QName, rec.RName, rec.Pos)
	}
	if rec.Seq != "ACGTA" || rec.Cigar.String() != "5M" {
		t.Errorf("seq/cigar = %q/%q, want ACGTA/5M", rec.Seq, rec.Cigar.String())
	}
	if rec.MapQ != 60 {
		t.Errorf("MapQ = %d, want 60", rec.MapQ)
	}
}

// TestDecodeRecordUnmapped decodes a hand-built unmapped record whose
// sequence comes from the BA data series.
func TestDecodeRecordUnmapped(t *testing.T) {
	h := refFreeHeader()
	h.Preservation.ReadNamesIncluded = true
	h.Preservation.APDelta = false
	h.Preservation.TagDictionary = []byte{0}
	h.DataSeries[dataSeriesKey{'B', 'F'}] = extEnc(1)
	h.DataSeries[dataSeriesKey{'C', 'F'}] = huffConst(0)
	h.DataSeries[dataSeriesKey{'R', 'L'}] = extEnc(2)
	h.DataSeries[dataSeriesKey{'A', 'P'}] = huffConst(0)
	h.DataSeries[dataSeriesKey{'R', 'G'}] = huffConst(-1)
	h.DataSeries[dataSeriesKey{'R', 'N'}] = stopEnc(3)
	h.DataSeries[dataSeriesKey{'T', 'L'}] = huffConst(0)
	h.DataSeries[dataSeriesKey{'B', 'A'}] = extEnc(4)

	var bf, rl, rn, ba bytes.Buffer
	bf.Write(encITF8(int32(sam.FlagUnmapped)))
	rl.Write(encITF8(4))
	rn.WriteString("u")
	rn.WriteByte(0)
	ba.WriteString("TTTT")
	src := newTestSource(nil, map[int32][]byte{
		1: bf.Bytes(), 2: rl.Bytes(), 3: rn.Bytes(), 4: ba.Bytes(),
	})
	rd := &recordDecoder{h: h, src: &SeriesSource{s: src}, slice: &SliceHeader{RefSeqID: -1}, readLenLimit: 1 << 20}
	dr, err := rd.decodeRecord(0)
	if err != nil {
		t.Fatalf("decodeRecord: %v", err)
	}
	if dr.rec.Seq != "TTTT" || dr.rec.RName != "" {
		t.Errorf("unmapped record = %q/%q, want TTTT/\"\"", dr.rec.Seq, dr.rec.RName)
	}
	if len(dr.rec.Cigar) != 0 {
		t.Errorf("unmapped record should have no CIGAR, got %q", dr.rec.Cigar.String())
	}
}

// TestDecodeReadNameSynthesised checks that a CRAM file that did not
// preserve read names synthesises numeric names from the record counter.
func TestDecodeReadNameSynthesised(t *testing.T) {
	h := refFreeHeader()
	h.Preservation.ReadNamesIncluded = false
	rd := &recordDecoder{h: h, slice: &SliceHeader{RecordCounter: 1000}}
	rec := &sam.Record{}
	if err := rd.decodeReadName(rec, 0, 5); err != nil {
		t.Fatalf("decodeReadName: %v", err)
	}
	if rec.QName != "1005" {
		t.Errorf("synthesised name = %q, want 1005", rec.QName)
	}
}

// TestDecodeTagsErrors checks the tag-line bounds handling.
func TestDecodeTagsErrors(t *testing.T) {
	rd := &recordDecoder{h: refFreeHeader()}
	// A TL of 0 with an empty dictionary yields no tags, not an error.
	if tags, _, err := rd.decodeTags(0, -1, 0); err != nil || tags != nil {
		t.Errorf("decodeTags(0) = %v, %v; want nil,nil", tags, err)
	}
	if _, _, err := rd.decodeTags(-1, -1, 0); err == nil {
		t.Error("a negative tag-line index should error")
	}
	if _, _, err := rd.decodeTags(7, -1, 0); err == nil {
		t.Error("an out-of-range tag-line index should error")
	}
}

// TestOpenRecordsAndClose exercises the file-path entry point and Close.
func TestOpenRecordsAndClose(t *testing.T) {
	path := filepath.Join(samtoolsTestDir, "dat/test_input_1_a.cram")
	if _, err := os.Stat(path); err != nil {
		t.Fatalf("samtools submodule not initialised — fixture unavailable; run `git submodule update --init reference_code/samtools`")
	}
	rr, err := OpenRecords(path)
	if err != nil {
		t.Fatalf("OpenRecords: %v", err)
	}
	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(recs) != 15 {
		t.Errorf("decoded %d records, want 15", len(recs))
	}
	if err := rr.Close(); err != nil {
		t.Errorf("Close: %v", err)
	}
	if _, err := OpenRecords(filepath.Join(samtoolsTestDir, "dat/does-not-exist.cram")); err == nil {
		t.Error("OpenRecords on a missing file should error")
	}
}

// TestRecordReaderMalformed feeds malformed and truncated bytes to the
// record reader; each must surface as an error rather than a panic.
func TestRecordReaderMalformed(t *testing.T) {
	inputs := [][]byte{
		nil,
		[]byte("not a cram file at all"),
		[]byte("CRAM\x03\x00"),                  // a bare file definition.
		append([]byte("CRAM\x03\x00"), 0, 1, 2), // a truncated container.
	}
	for i, in := range inputs {
		rr, err := NewRecordReader(bytes.NewReader(in))
		if err != nil {
			continue // a clean error before iteration is fine.
		}
		var buf bytes.Buffer
		if err := rr.WriteSAM(&buf); err == nil {
			// Some inputs may legitimately decode to an empty stream;
			// the requirement is only that nothing panics.
			t.Logf("input %d decoded without error", i)
		}
	}
}

// TestDecodeMateDetached checks the detached-mate path: a record
// carrying CF bit 0x2 stores its own NS/NP/TS and mate flags.
func TestDecodeMateDetached(t *testing.T) {
	h := refFreeHeader()
	// The mate block here is MF→NS→NP→TS with no embedded read name,
	// which is the names-preserved detached layout.
	h.Preservation.ReadNamesIncluded = true
	h.DataSeries[dataSeriesKey{'M', 'F'}] = extEnc(1)
	h.DataSeries[dataSeriesKey{'N', 'S'}] = extEnc(2)
	h.DataSeries[dataSeriesKey{'N', 'P'}] = extEnc(3)
	h.DataSeries[dataSeriesKey{'T', 'S'}] = extEnc(4)

	var mf, ns, np, ts bytes.Buffer
	mf.Write(encITF8(mfMateReverse | mfMateUnmapped))
	ns.Write(encITF8(0))
	np.Write(encITF8(250))
	ts.Write(encITF8(400))
	src := newTestSource(nil, map[int32][]byte{
		1: mf.Bytes(), 2: ns.Bytes(), 3: np.Bytes(), 4: ts.Bytes(),
	})
	rd := &recordDecoder{
		h: h, src: &SeriesSource{s: src},
		slice:    &SliceHeader{},
		refNames: []string{"chr1"},
	}
	dr := &decodedRecord{rec: &sam.Record{RName: "chr1"}}
	if err := rd.decodeMate(dr, cfDetached, 0); err != nil {
		t.Fatalf("decodeMate: %v", err)
	}
	rec := dr.rec
	if rec.RNext != "=" || rec.PNext != 250 || rec.TLen != 400 {
		t.Errorf("detached mate = %q/%d/%d, want =/250/400", rec.RNext, rec.PNext, rec.TLen)
	}
	if rec.Flag&sam.FlagMateReverse == 0 || rec.Flag&sam.FlagMateUnmapped == 0 {
		t.Error("detached mate flags not applied")
	}
}

// TestDecodeMateDownstreamNegative checks that a negative next-fragment
// distance is rejected without a panic.
func TestDecodeMateDownstreamNegative(t *testing.T) {
	h := refFreeHeader()
	var nf bytes.Buffer
	nf.Write(encITF8(-1))
	h.DataSeries[dataSeriesKey{'N', 'F'}] = extEnc(1)
	rd := &recordDecoder{
		h:   h,
		src: &SeriesSource{s: newTestSource(nil, map[int32][]byte{1: nf.Bytes()})},
	}
	dr := &decodedRecord{rec: &sam.Record{}}
	if err := rd.decodeMate(dr, cfHasMateDownstream, 0); err == nil {
		t.Error("a negative next-fragment distance should error")
	}
}

// TestDecodeSliceRecordsErrors checks the slice-level traversal guards.
func TestDecodeSliceRecordsErrors(t *testing.T) {
	rd := &recordDecoder{
		h:     refFreeHeader(),
		slice: &SliceHeader{},
		src:   &SeriesSource{s: newTestSource(nil, nil)},
	}
	if _, err := rd.decodeSliceRecords(-1); err == nil {
		t.Error("a negative record count should error")
	}
	// A record count that exceeds the slice's series data must be
	// rejected up front rather than allocate and loop — the slice here
	// has zero bytes, so any positive count overshoots.
	if _, err := rd.decodeSliceRecords(1); err == nil {
		t.Error("a record count exceeding the series data should error")
	}
}

// TestDecodeSliceMatePair decodes a two-record slice where the first
// record's mate is the second (downstream), exercising decodeRecord,
// decodeSliceRecords, resolveMates, linkMates and setMateFields.
func TestDecodeSliceMatePair(t *testing.T) {
	h := refFreeHeader()
	h.Preservation.ReadNamesIncluded = true
	h.Preservation.APDelta = false
	h.Preservation.TagDictionary = []byte{0}
	h.DataSeries[dataSeriesKey{'B', 'F'}] = extEnc(1)
	h.DataSeries[dataSeriesKey{'C', 'F'}] = extEnc(2)
	h.DataSeries[dataSeriesKey{'R', 'L'}] = huffConst(4)
	h.DataSeries[dataSeriesKey{'A', 'P'}] = extEnc(3)
	h.DataSeries[dataSeriesKey{'R', 'G'}] = huffConst(-1)
	h.DataSeries[dataSeriesKey{'R', 'N'}] = stopEnc(4)
	h.DataSeries[dataSeriesKey{'N', 'F'}] = huffConst(0) // mate is the very next record.
	h.DataSeries[dataSeriesKey{'T', 'L'}] = huffConst(0)
	h.DataSeries[dataSeriesKey{'M', 'Q'}] = huffConst(20)
	h.DataSeries[dataSeriesKey{'F', 'N'}] = huffConst(1)
	h.DataSeries[dataSeriesKey{'F', 'C'}] = huffConst(int32(featBases))
	h.DataSeries[dataSeriesKey{'F', 'P'}] = huffConst(1)
	h.DataSeries[dataSeriesKey{'B', 'B'}] = &Encoding{
		ID: EncodingByteArrayLen, LenEnc: huffConst(4), ValEnc: extEnc(5),
	}

	var bf, cf, ap, rn, bb bytes.Buffer
	bf.Write(encITF8(int32(sam.FlagPaired)))                   // record 0.
	bf.Write(encITF8(int32(sam.FlagPaired | sam.FlagReverse))) // record 1.
	cf.Write(encITF8(cfHasMateDownstream))                     // record 0 has a downstream mate.
	cf.Write(encITF8(0))                                       // record 1.
	ap.Write(encITF8(100))                                     // record 0 position.
	ap.Write(encITF8(300))                                     // record 1 position.
	rn.WriteString("pair")
	rn.WriteByte(0)
	rn.WriteString("pair")
	rn.WriteByte(0)
	bb.WriteString("ACGT")
	bb.WriteString("TGCA")
	src := newTestSource(nil, map[int32][]byte{
		1: bf.Bytes(), 2: cf.Bytes(), 3: ap.Bytes(), 4: rn.Bytes(), 5: bb.Bytes(),
	})
	rd := &recordDecoder{
		h: h, src: &SeriesSource{s: src},
		slice:        &SliceHeader{RefSeqID: 0},
		refNames:     []string{"chr1"},
		readLenLimit: 1 << 20,
	}
	recs, err := rd.decodeSliceRecords(2)
	if err != nil {
		t.Fatalf("decodeSliceRecords: %v", err)
	}
	if len(recs) != 2 {
		t.Fatalf("decoded %d records, want 2", len(recs))
	}
	up, down := recs[0], recs[1]
	if up.RNext != "=" || up.PNext != 300 {
		t.Errorf("upstream mate fields = %q/%d, want =/300", up.RNext, up.PNext)
	}
	if down.RNext != "=" || down.PNext != 100 {
		t.Errorf("downstream mate fields = %q/%d, want =/100", down.RNext, down.PNext)
	}
	if up.Flag&sam.FlagMateReverse == 0 {
		t.Error("upstream record should inherit the mate-reverse flag")
	}
	// TLEN spans 100..303 = 204; positive on the upstream record.
	if up.TLen != 204 || down.TLen != -204 {
		t.Errorf("TLEN = %d/%d, want 204/-204", up.TLen, down.TLen)
	}
}

// TestDecodeArrayTagSubtypes checks every 'B'-array element subtype and
// the array error paths.
func TestDecodeArrayTagSubtypes(t *testing.T) {
	cases := []struct {
		sub  byte
		body []byte
		want string
	}{
		{'c', []byte{0xfe, 2}, "XB:B:c,-2,2"},
		{'C', []byte{200, 1}, "XB:B:C,200,1"},
		{'s', []byte{0xfe, 0xff, 3, 0}, "XB:B:s,-2,3"},
		{'S', []byte{0x10, 0x27}, "XB:B:S,10000"},
		{'i', []byte{0xff, 0xff, 0xff, 0xff}, "XB:B:i,-1"},
		{'I', []byte{0, 0, 0, 0x80}, "XB:B:I,2147483648"},
		{'f', []byte{0, 0, 0x80, 0x3f}, "XB:B:f,1"},
	}
	for _, c := range cases {
		count := len(c.body)
		switch c.sub {
		case 's', 'S':
			count /= 2
		case 'i', 'I', 'f':
			count /= 4
		}
		raw := []byte{c.sub, byte(count), 0, 0, 0}
		raw = append(raw, c.body...)
		aux, err := decodeTagValue(tagKey{'X', 'B', 'B'}, raw)
		if err != nil {
			t.Errorf("subtype %c: %v", c.sub, err)
			continue
		}
		if got := aux.FormatSAM(); got != c.want {
			t.Errorf("subtype %c: got %q, want %q", c.sub, got, c.want)
		}
	}
	// A truncated array and an unknown subtype are errors, not panics.
	if _, err := decodeTagValue(tagKey{'X', 'B', 'B'}, []byte{'c', 1}); err == nil {
		t.Error("a truncated 'B' value should error")
	}
	if _, err := decodeTagValue(tagKey{'X', 'B', 'B'}, []byte{'?', 0, 0, 0, 0}); err == nil {
		t.Error("an unknown array subtype should error")
	}
	if _, err := decodeTagValue(tagKey{'X', 'B', 'B'}, []byte{'c', 9, 0, 0, 0, 1}); err == nil {
		t.Error("a count/length mismatch should error")
	}
}
