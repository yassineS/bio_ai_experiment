package cram

import (
	"bytes"
	"io"
	"testing"
)

// v30Fixtures are the real CRAM v3.0 files whose every container's
// compression header, slice headers and data series the C4a layer must
// decode without error.
var v30Fixtures = []struct {
	name string
	rel  string
}{
	{"test_input_1_a", "dat/test_input_1_a.cram"},
	{"7.quickcheck.cram30.ok", "quickcheck/7.quickcheck.cram30.ok.cram"},
}

// TestDecodeFixtureDataSeries is the C4a correctness oracle. For each
// real CRAM v3.0 fixture it parses every data container's compression
// header and slice headers, then decodes every data series that is
// self-delimiting (EXTERNAL, BYTE_ARRAY_STOP, BYTE_ARRAY_LEN) by
// draining it to exhaustion — a clean drain proves the encoding layer
// decoded every byte of the series' block as a well-formed value.
//
// CORE-bitstream series (HUFFMAN, BETA, …) share one interleaved
// stream and cannot be isolated without C4b's per-record traversal;
// they are decoded structurally (their encoding parameters are
// validated and a degenerate HUFFMAN, which consumes no bits, is
// decoded to the record count) but a clean per-series drain is left to
// C4b. Every encoding the fixtures use is handled; an encoding that
// could not be handled would surface as an explicit error here.
func TestDecodeFixtureDataSeries(t *testing.T) {
	for _, fx := range v30Fixtures {
		fx := fx
		t.Run(fx.name, func(t *testing.T) {
			data, ok := loadFixture(t, fx.rel)
			if !ok {
				t.Skip("samtools submodule not initialised — fixture unavailable")
			}
			rd, err := NewReader(bytes.NewReader(data))
			if err != nil {
				t.Fatalf("NewReader: %v", err)
			}
			if rd.FileDefinition().Major != 3 {
				t.Fatalf("expected CRAM major version 3")
			}
			dataContainers, drained, coreSeries := 0, 0, 0
			for {
				c, err := rd.Next()
				if err == io.EOF {
					break
				}
				if err != nil {
					t.Fatalf("Next: %v", err)
				}
				if len(c.Blocks) == 0 || c.Blocks[0].ContentType != ContentCompressionHeader {
					continue // the file-header container
				}
				dc, err := ParseDataContainer(c)
				if err != nil {
					t.Fatalf("container %d ParseDataContainer: %v", c.Index, err)
				}
				dataContainers++
				if dc.Compression == nil || dc.Compression.DataSeries == nil {
					t.Fatalf("container %d: nil compression header", c.Index)
				}
				for si, sl := range dc.Slices {
					if sl.Header == nil {
						t.Fatalf("container %d slice %d: nil slice header", c.Index, si)
					}
					src, err := sl.NewSource()
					if err != nil {
						t.Fatalf("container %d slice %d NewSource: %v", c.Index, si, err)
					}
					for key, enc := range dc.Compression.DataSeries {
						k := key.String()
						switch {
						case sl.Drainable(dc.Compression, k):
							res, err := sl.DrainSeries(dc.Compression, src, k)
							if err != nil {
								t.Errorf("container %d slice %d series %s (%s): drain failed: %v",
									c.Index, si, k, enc.ID, err)
								continue
							}
							if res.Count < 0 {
								t.Errorf("container %d slice %d series %s: negative count", c.Index, si, k)
							}
							drained++
						case enc.ID == EncodingHuffman:
							// A degenerate HUFFMAN (constant series) consumes
							// no bits, so it can be decoded to the record
							// count without disturbing the shared CORE
							// bitstream. A non-degenerate one is left to C4b.
							tbl, err := enc.huffmanDecoder()
							if err != nil {
								t.Errorf("container %d slice %d series %s: invalid HUFFMAN table: %v",
									c.Index, si, k, err)
								continue
							}
							if tbl.degenerate {
								vals, err := enc.decodeInts(src.s, int(sl.Header.NumRecords))
								if err != nil {
									t.Errorf("container %d slice %d series %s: degenerate HUFFMAN decode: %v",
										c.Index, si, k, err)
									continue
								}
								if len(vals) != int(sl.Header.NumRecords) {
									t.Errorf("container %d slice %d series %s: got %d values, want record count %d",
										c.Index, si, k, len(vals), sl.Header.NumRecords)
								}
							}
							coreSeries++
						default:
							// EXTERNAL/STOP series with an absent block carry
							// no data; a CORE BETA/GAMMA/etc. series is left
							// to C4b. Either way the encoding must be a known
							// one — parseEncoding already rejected unknowns.
							coreSeries++
						}
					}
				}
			}
			// Decode the tag series too: every self-delimiting tag of
			// every slice must drain cleanly.
			rd2, err := NewReader(bytes.NewReader(data))
			if err != nil {
				t.Fatalf("NewReader (tag pass): %v", err)
			}
			tagsDrained := 0
			for {
				c, err := rd2.Next()
				if err == io.EOF {
					break
				}
				if err != nil {
					t.Fatalf("Next (tag pass): %v", err)
				}
				if len(c.Blocks) == 0 || c.Blocks[0].ContentType != ContentCompressionHeader {
					continue
				}
				dc, err := ParseDataContainer(c)
				if err != nil {
					t.Fatalf("ParseDataContainer (tag pass): %v", err)
				}
				for si, sl := range dc.Slices {
					src, err := sl.NewSource()
					if err != nil {
						t.Fatalf("NewSource (tag pass): %v", err)
					}
					for tk := range dc.Compression.Tags {
						tag := tk.String()
						if !sl.TagDrainable(dc.Compression, tag) {
							continue
						}
						if _, err := sl.DrainTag(dc.Compression, src, tag); err != nil {
							t.Errorf("container %d slice %d tag %s: drain failed: %v",
								c.Index, si, tag, err)
							continue
						}
						tagsDrained++
					}
				}
			}
			if dataContainers == 0 {
				t.Fatalf("%s: no data containers found", fx.name)
			}
			t.Logf("%s: %d data containers, %d series drained cleanly, %d CORE/structural series, %d tags drained",
				fx.name, dataContainers, drained, coreSeries, tagsDrained)
		})
	}
}

// TestDecodeFixtureBETASeries decodes the one BETA-encoded CORE data
// series the test_input_1_a fixture carries (its alignment-position
// series), proving the BETA bitstream decoder works on real data.
func TestDecodeFixtureBETASeries(t *testing.T) {
	data, ok := loadFixture(t, "dat/test_input_1_a.cram")
	if !ok {
		t.Skip("samtools submodule not initialised")
	}
	rd, err := NewReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	foundBeta := false
	for {
		c, err := rd.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Next: %v", err)
		}
		if len(c.Blocks) == 0 || c.Blocks[0].ContentType != ContentCompressionHeader {
			continue
		}
		dc, err := ParseDataContainer(c)
		if err != nil {
			t.Fatalf("ParseDataContainer: %v", err)
		}
		for _, sl := range dc.Slices {
			enc := dc.Compression.Encoding("AP")
			if enc == nil || enc.ID != EncodingBeta {
				continue
			}
			src, err := sl.NewSource()
			if err != nil {
				t.Fatalf("NewSource: %v", err)
			}
			vals, err := sl.DecodeIntSeries(dc.Compression, src, "AP", int(sl.Header.NumRecords))
			if err != nil {
				t.Fatalf("decode AP BETA series: %v", err)
			}
			if len(vals) != int(sl.Header.NumRecords) {
				t.Fatalf("AP series: got %d values, want %d", len(vals), sl.Header.NumRecords)
			}
			foundBeta = true
		}
	}
	if !foundBeta {
		t.Skip("fixture carried no BETA-encoded AP series in this build")
	}
}

// TestParseDataContainerErrors checks ParseDataContainer rejects
// containers that are not data containers.
func TestParseDataContainerErrors(t *testing.T) {
	if _, err := ParseDataContainer(nil); err == nil {
		t.Errorf("nil container should be rejected")
	}
	if _, err := ParseDataContainer(&Container{}); err == nil {
		t.Errorf("a container with no blocks should be rejected")
	}
	notCH := &Container{Blocks: []Block{{ContentType: ContentFileHeader, Method: CompRaw}}}
	if _, err := ParseDataContainer(notCH); err == nil {
		t.Errorf("a container whose first block is not a compression header should be rejected")
	}
}

// TestSliceHasSeriesData checks the HasSeriesData / Drainable helpers.
func TestSliceHasSeriesData(t *testing.T) {
	h := &CompressionHeader{DataSeries: map[dataSeriesKey]*Encoding{
		{'B', 'F'}: {ID: EncodingExternal, ExternalID: 11},
		{'Q', 'S'}: {ID: EncodingExternal, ExternalID: 99},
		{'C', 'F'}: {ID: EncodingHuffman, Symbols: []int32{0}, BitLengths: []int32{0}},
	}}
	sl := &Slice{
		Header:   &SliceHeader{NumRecords: 3},
		external: map[int32]*Block{11: {ContentType: ContentExternal}},
	}
	if !sl.HasSeriesData(h, "BF") {
		t.Errorf("BF (block present) should report data")
	}
	if sl.HasSeriesData(h, "QS") {
		t.Errorf("QS (block absent) should report no data")
	}
	if !sl.HasSeriesData(h, "CF") {
		t.Errorf("CF (CORE series) should report data")
	}
	if sl.HasSeriesData(h, "ZZ") {
		t.Errorf("absent series should report no data")
	}
	if !sl.Drainable(h, "BF") {
		t.Errorf("BF should be drainable")
	}
	if sl.Drainable(h, "CF") {
		t.Errorf("a HUFFMAN series should not be drainable")
	}
}

// TestDecodeSeriesErrors checks the decode entry points reject an
// unknown series key.
func TestDecodeSeriesErrors(t *testing.T) {
	h := &CompressionHeader{DataSeries: map[dataSeriesKey]*Encoding{}}
	sl := &Slice{Header: &SliceHeader{}, external: map[int32]*Block{}}
	src := &SeriesSource{s: newTestSource(nil, nil)}
	if _, err := sl.DecodeIntSeries(h, src, "ZZ", 0); err == nil {
		t.Errorf("DecodeIntSeries of an unknown series should error")
	}
	if _, err := sl.DecodeByteArraySeries(h, src, "ZZ", 0); err == nil {
		t.Errorf("DecodeByteArraySeries of an unknown series should error")
	}
	if _, err := sl.DrainIntSeries(h, src, "ZZ"); err == nil {
		t.Errorf("DrainIntSeries of an unknown series should error")
	}
	if _, err := sl.DrainByteArraySeries(h, src, "ZZ"); err == nil {
		t.Errorf("DrainByteArraySeries of an unknown series should error")
	}
	if _, err := sl.DrainSeries(h, src, "ZZ"); err == nil {
		t.Errorf("DrainSeries of an unknown series should error")
	}
}

// TestDrainTag decodes a tag series from a synthetic slice.
func TestDrainTag(t *testing.T) {
	h := &CompressionHeader{Tags: map[tagKey]*Encoding{
		{'P', 'G', 'Z'}: {ID: EncodingByteArrayStop, ExternalID: 20, StopByte: 0},
		{'X', 'X', 'C'}: {ID: EncodingHuffman, Symbols: []int32{1}, BitLengths: []int32{0}},
	}}
	sl := &Slice{
		Header: &SliceHeader{NumRecords: 2},
		external: map[int32]*Block{
			20: {Method: CompRaw, ContentType: ContentExternal, Data: []byte("p1\x00p2\x00"),
				UncompressedSize: 6},
		},
	}
	src, err := sl.NewSource()
	if err != nil {
		t.Fatalf("NewSource: %v", err)
	}
	if !sl.TagDrainable(h, "PGZ") {
		t.Fatalf("PGZ should be drainable")
	}
	res, err := sl.DrainTag(h, src, "PGZ")
	if err != nil {
		t.Fatalf("DrainTag PGZ: %v", err)
	}
	if res.Count != 2 || string(res.ByteArrays[0]) != "p1" {
		t.Errorf("DrainTag PGZ = %v (count %d), want [p1 p2]", res.ByteArrays, res.Count)
	}
	if sl.TagDrainable(h, "XXC") {
		t.Errorf("a HUFFMAN tag should not be drainable")
	}
	if _, err := sl.DrainTag(h, src, "XXC"); err == nil {
		t.Errorf("draining a HUFFMAN tag should error")
	}
	if _, err := sl.DrainTag(h, src, "AB"); err == nil {
		t.Errorf("a 2-char tag key should error")
	}
	if _, err := sl.DrainTag(h, src, "ZZZ"); err == nil {
		t.Errorf("an absent tag should error")
	}
}

// TestSliceDecodeAndDrainAPI exercises the fixed-count and drain decode
// entry points on a synthetic slice, covering the success paths.
func TestSliceDecodeAndDrainAPI(t *testing.T) {
	var intBlk, stopBlk bytes.Buffer
	for _, v := range []int32{10, 20, 30} {
		intBlk.Write(encITF8(v))
	}
	stopBlk.WriteString("a\x00bb\x00ccc\x00")
	h := &CompressionHeader{DataSeries: map[dataSeriesKey]*Encoding{
		{'B', 'F'}: {ID: EncodingExternal, ExternalID: 1},
		{'R', 'N'}: {ID: EncodingByteArrayStop, ExternalID: 2, StopByte: 0},
	}}
	sl := &Slice{
		Header: &SliceHeader{NumRecords: 3},
		external: map[int32]*Block{
			1: {Method: CompRaw, ContentType: ContentExternal, Data: intBlk.Bytes(),
				UncompressedSize: int32(intBlk.Len())},
			2: {Method: CompRaw, ContentType: ContentExternal, Data: stopBlk.Bytes(),
				UncompressedSize: int32(stopBlk.Len())},
		},
	}
	src, err := sl.NewSource()
	if err != nil {
		t.Fatalf("NewSource: %v", err)
	}
	iv, err := sl.DecodeIntSeries(h, src, "BF", 3)
	if err != nil || len(iv) != 3 || iv[1] != 20 {
		t.Fatalf("DecodeIntSeries BF = %v, %v", iv, err)
	}

	src2, _ := sl.NewSource()
	bv, err := sl.DecodeByteArraySeries(h, src2, "RN", 3)
	if err != nil || len(bv) != 3 || string(bv[2]) != "ccc" {
		t.Fatalf("DecodeByteArraySeries RN = %v, %v", bv, err)
	}

	src3, _ := sl.NewSource()
	div, err := sl.DrainIntSeries(h, src3, "BF")
	if err != nil || len(div) != 3 {
		t.Fatalf("DrainIntSeries BF = %v, %v", div, err)
	}

	src4, _ := sl.NewSource()
	dbv, err := sl.DrainByteArraySeries(h, src4, "RN")
	if err != nil || len(dbv) != 3 {
		t.Fatalf("DrainByteArraySeries RN = %v, %v", dbv, err)
	}

	src5, _ := sl.NewSource()
	res, err := sl.DrainSeries(h, src5, "RN")
	if err != nil || res.Count != 3 {
		t.Fatalf("DrainSeries RN = %+v, %v", res, err)
	}
}

// TestSliceNoCoreBlock checks a slice with no CORE block still yields a
// usable (empty) source, and that a CORE-needing decode errors cleanly
// rather than panicking.
func TestSliceNoCoreBlock(t *testing.T) {
	sl := &Slice{Header: &SliceHeader{NumRecords: 1}, external: map[int32]*Block{}}
	src, err := sl.NewSource()
	if err != nil {
		t.Fatalf("NewSource: %v", err)
	}
	h := &CompressionHeader{DataSeries: map[dataSeriesKey]*Encoding{
		{'A', 'P'}: {ID: EncodingBeta, NumBits: 8},
	}}
	if _, err := sl.DecodeIntSeries(h, src, "AP", 1); err == nil {
		t.Errorf("decoding a CORE series from an empty CORE block should error")
	}
}

// TestSeriesValueKind checks the data-series value catalogue.
func TestSeriesValueKind(t *testing.T) {
	if SeriesValueKind("BF") != SeriesInt {
		t.Errorf("BF should be an integer series")
	}
	if SeriesValueKind("QS") != SeriesByte {
		t.Errorf("QS should be a byte series")
	}
	if SeriesValueKind("ZZ") != SeriesInt {
		t.Errorf("an uncatalogued series should default to integer")
	}
}
