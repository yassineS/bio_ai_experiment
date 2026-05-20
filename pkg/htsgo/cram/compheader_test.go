package cram

import (
	"bytes"
	"testing"
)

// buildMap frames a CRAM compression-header map: the ITF-8 byte size of
// (entry count + body), then the ITF-8 entry count, then the body.
func buildMap(count int32, body []byte) []byte {
	var inner bytes.Buffer
	inner.Write(encITF8(count))
	inner.Write(body)
	var out bytes.Buffer
	out.Write(encITF8(int32(inner.Len())))
	out.Write(inner.Bytes())
	return out.Bytes()
}

// buildCompressionHeader assembles a minimal compression-header payload
// from the three maps' bodies and entry counts.
func buildCompressionHeader(presCount int32, presBody []byte,
	dsCount int32, dsBody []byte, tagCount int32, tagBody []byte) []byte {
	var b bytes.Buffer
	b.Write(buildMap(presCount, presBody))
	b.Write(buildMap(dsCount, dsBody))
	b.Write(buildMap(tagCount, tagBody))
	return b.Bytes()
}

// TestParseCompressionHeaderBasic builds and parses a compression header
// with one preservation entry, one data series and one tag.
func TestParseCompressionHeaderBasic(t *testing.T) {
	pres := []byte{'R', 'N', 0} // RN = false
	var ds bytes.Buffer
	ds.WriteByte('B')
	ds.WriteByte('F')
	ds.Write(encEncoding(EncodingExternal, encITF8(11)))
	var tag bytes.Buffer
	// Tag key 'X','A','Z' packed into an ITF-8 integer.
	tag.Write(encITF8(int32('X')<<16 | int32('A')<<8 | int32('Z')))
	tag.Write(encEncoding(EncodingByteArrayStop, append([]byte{0}, encITF8(20)...)))

	payload := buildCompressionHeader(1, pres, 1, ds.Bytes(), 1, tag.Bytes())
	h, err := parseCompressionHeader(payload)
	if err != nil {
		t.Fatalf("parseCompressionHeader: %v", err)
	}
	if h.Preservation.ReadNamesIncluded {
		t.Errorf("RN should be false")
	}
	if !h.Preservation.hasRN {
		t.Errorf("RN should be marked present")
	}
	bf := h.Encoding("BF")
	if bf == nil || bf.ID != EncodingExternal || bf.ExternalID != 11 {
		t.Errorf("BF encoding wrong: %+v", bf)
	}
	tk := tagKey{'X', 'A', 'Z'}
	if te := h.Tags[tk]; te == nil || te.ID != EncodingByteArrayStop {
		t.Errorf("tag XAZ encoding wrong: %+v", te)
	}
}

// TestParseCompressionHeaderDefaults checks the preservation-map
// booleans take their CRAM defaults when no entry is written.
func TestParseCompressionHeaderDefaults(t *testing.T) {
	payload := buildCompressionHeader(0, nil, 0, nil, 0, nil)
	h, err := parseCompressionHeader(payload)
	if err != nil {
		t.Fatalf("parseCompressionHeader: %v", err)
	}
	pm := h.Preservation
	if !pm.ReadNamesIncluded || !pm.APDelta || !pm.ReferenceRequired {
		t.Errorf("preservation defaults wrong: %+v", pm)
	}
	if pm.hasRN || pm.hasAP || pm.hasRR {
		t.Errorf("no preservation entry should be marked present")
	}
}

// TestParseCompressionHeaderSMandTD checks the SM (substitution matrix)
// and TD (tag dictionary) raw values are retained.
func TestParseCompressionHeaderSMandTD(t *testing.T) {
	var pres bytes.Buffer
	pres.WriteString("SM")
	pres.Write([]byte{1, 2, 3, 4, 5})
	pres.WriteString("TD")
	td := []byte("RGZ\x00MDZ\x00")
	pres.Write(encITF8(int32(len(td))))
	pres.Write(td)
	pres.WriteString("AP")
	pres.WriteByte(1)

	payload := buildCompressionHeader(3, pres.Bytes(), 0, nil, 0, nil)
	h, err := parseCompressionHeader(payload)
	if err != nil {
		t.Fatalf("parseCompressionHeader: %v", err)
	}
	if !bytes.Equal(h.Preservation.SubstitutionMatrix, []byte{1, 2, 3, 4, 5}) {
		t.Errorf("SM = %v, want [1 2 3 4 5]", h.Preservation.SubstitutionMatrix)
	}
	if !bytes.Equal(h.Preservation.TagDictionary, td) {
		t.Errorf("TD = %q, want %q", h.Preservation.TagDictionary, td)
	}
	if !h.Preservation.APDelta {
		t.Errorf("AP should be true")
	}
}

// TestParseCompressionHeaderUnknownKey rejects an unknown preservation
// key, whose value length is not self-describing.
func TestParseCompressionHeaderUnknownKey(t *testing.T) {
	pres := []byte{'Z', 'Z', 0}
	payload := buildCompressionHeader(1, pres, 0, nil, 0, nil)
	if _, err := parseCompressionHeader(payload); err == nil {
		t.Errorf("unknown preservation key should be rejected")
	}
}

// TestParseCompressionHeaderTruncated checks the parser errors — never
// panics — on input truncated at every byte offset.
func TestParseCompressionHeaderTruncated(t *testing.T) {
	var ds bytes.Buffer
	ds.WriteString("BF")
	ds.Write(encEncoding(EncodingExternal, encITF8(11)))
	full := buildCompressionHeader(1, []byte{'R', 'N', 1}, 1, ds.Bytes(), 0, nil)
	for n := 0; n < len(full); n++ {
		func() {
			defer func() {
				if r := recover(); r != nil {
					t.Errorf("parseCompressionHeader panicked at truncation %d: %v", n, r)
				}
			}()
			if _, err := parseCompressionHeader(full[:n]); err == nil && n < len(full) {
				t.Errorf("expected error for compression header truncated to %d bytes", n)
			}
		}()
	}
}

// TestParseCompressionHeaderBadMapSize rejects a map whose declared size
// overruns the payload.
func TestParseCompressionHeaderBadMapSize(t *testing.T) {
	// A preservation map claiming a huge size.
	payload := append(encITF8(9999), encITF8(0)...)
	if _, err := parseCompressionHeader(payload); err == nil {
		t.Errorf("oversized map should be rejected")
	}
}

// TestCompressionHeaderEncodingHelper checks the Encoding accessor's key
// validation.
func TestCompressionHeaderEncodingHelper(t *testing.T) {
	h := &CompressionHeader{DataSeries: map[dataSeriesKey]*Encoding{
		{'B', 'F'}: {ID: EncodingExternal},
	}}
	if h.Encoding("BF") == nil {
		t.Errorf("Encoding(BF) should be found")
	}
	if h.Encoding("XYZ") != nil {
		t.Errorf("Encoding with a 3-char key should return nil")
	}
	if h.Encoding("ZZ") != nil {
		t.Errorf("Encoding for an absent series should return nil")
	}
}
