package cram

import (
	"bytes"
	"compress/gzip"
	"io"
	"testing"
)

// FuzzCompressionHeader runs the compression-header parser over
// arbitrary input. The parser must never panic; any malformed map,
// encoding or overrun must surface as a returned error. When the parse
// succeeds, the resulting maps are walked to ensure no later operation
// (encoding-id stringer, Huffman table build) panics either.
func FuzzCompressionHeader(f *testing.F) {
	// Seed with a hand-built minimal header and the real compression
	// headers extracted from the v3.0 fixtures when available.
	var ds bytes.Buffer
	ds.WriteString("BF")
	ds.Write(encEncoding(EncodingExternal, encITF8(11)))
	f.Add(buildCompressionHeader(1, []byte{'R', 'N', 1}, 1, ds.Bytes(), 0, nil))
	f.Add(buildCompressionHeader(0, nil, 0, nil, 0, nil))
	f.Add([]byte{})
	for _, fx := range v30Fixtures {
		if data, ok := readFixtureNoT(fx.rel); ok {
			for _, p := range compressionHeaderPayloads(data) {
				f.Add(p)
			}
		}
	}
	f.Fuzz(func(t *testing.T, data []byte) {
		h, err := parseCompressionHeader(data)
		if err != nil {
			return
		}
		// A successful parse must yield walkable maps.
		for k, enc := range h.DataSeries {
			_ = k.String()
			_ = enc.idString()
			if enc.ID == EncodingHuffman {
				_, _ = enc.huffmanDecoder()
			}
		}
		for k, enc := range h.Tags {
			_ = k.String()
			_ = enc.idString()
		}
		_ = h.Encoding("BF")
	})
}

// FuzzSliceHeader runs the slice-header parser over arbitrary input. It
// must never panic; truncated or corrupt fields must surface as errors.
func FuzzSliceHeader(f *testing.F) {
	f.Add(buildSliceHeader(SliceHeader{
		RefSeqID: 1, NumBlocks: 2, BlockContentIDs: []int32{1, 2}, EmbeddedRefID: -1,
	}))
	f.Add([]byte{})
	f.Fuzz(func(t *testing.T, data []byte) {
		// Exercise both the v2 (ITF-8 counter) and v3+ (LTF-8 counter)
		// branches; neither may panic on arbitrary input.
		for _, major := range []uint8{2, 3} {
			sh, err := parseSliceHeader(data, major)
			if err != nil {
				continue
			}
			// A successful parse must produce a self-consistent header.
			if int(sh.NumBlocks) < 0 {
				t.Fatalf("parsed a negative block count %d", sh.NumBlocks)
			}
		}
	})
}

// FuzzRecordReader runs the full CRAM record decode — file definition,
// containers, compression headers, slices, the per-record traversal,
// read-feature decode and SAM emission — over arbitrary input. A
// malformed stream must surface as a returned error, never a panic; the
// decode-to-SAM path must be panic-free on every input. A successful
// decode is additionally re-emitted as SAM to exercise the writer.
func FuzzRecordReader(f *testing.F) {
	f.Add([]byte{})
	f.Add([]byte("CRAM"))
	for _, fx := range v30Fixtures {
		if data, ok := readFixtureNoT(fx.rel); ok {
			f.Add(data)
		}
	}
	f.Fuzz(func(t *testing.T, data []byte) {
		rr, err := NewRecordReader(bytes.NewReader(data))
		if err != nil {
			return
		}
		var buf bytes.Buffer
		// WriteSAM drives Read, the per-record traversal and the SAM
		// writer; any malformed structure must come back as an error.
		_ = rr.WriteSAM(&buf)
		_ = rr.NeedsReference()
	})
}

// FuzzTagValue runs the tag-value decoder over arbitrary bytes for every
// SAM value type. It must never panic; a truncated or malformed value
// must surface as a returned error.
func FuzzTagValue(f *testing.F) {
	f.Add([]byte{3, 0, 0, 0})
	f.Add([]byte("abc\x00"))
	f.Add([]byte{'S', 2, 0, 0, 0, 1, 2, 3, 4})
	f.Add([]byte{})
	f.Fuzz(func(t *testing.T, raw []byte) {
		for _, typ := range []byte{'A', 'c', 'C', 's', 'S', 'i', 'I', 'f', 'Z', 'H', 'B'} {
			aux, err := decodeTagValue(tagKey{'X', 'X', typ}, raw)
			if err == nil {
				// A clean decode must format without panicking.
				_ = aux.FormatSAM()
			}
		}
	})
}

// FuzzReadCRAI runs the .crai index parser over arbitrary input. A
// .crai is a gzip-compressed TSV; the parser must never panic — every
// malformed gzip stream, non-integer field, wrong field count or
// out-of-range value must surface as a returned error. When the parse
// succeeds the index is additionally queried to exercise the overlap
// arithmetic on the decoded entries.
func FuzzReadCRAI(f *testing.F) {
	// Seed with a valid gzip-compressed .crai body and some raw,
	// definitely-not-gzip inputs.
	seed := func(text string) []byte {
		var buf bytes.Buffer
		gw := gzip.NewWriter(&buf)
		gw.Write([]byte(text))
		gw.Close()
		return buf.Bytes()
	}
	f.Add(seed("0\t100\t50\t1024\t12\t300\n1\t10\t90\t4096\t0\t500\n"))
	f.Add(seed("-1\t0\t0\t8192\t0\t128\n"))
	f.Add(seed(""))
	f.Add([]byte{})
	f.Add([]byte("not gzip at all"))
	f.Add([]byte{0x1f, 0x8b}) // a truncated gzip magic.
	f.Fuzz(func(t *testing.T, data []byte) {
		idx, err := ReadCRAI(bytes.NewReader(data))
		if err != nil {
			return
		}
		// A successful parse must yield queryable entries; the query must
		// not panic on any decoded entry, including extreme coordinates.
		for _, e := range idx.Entries {
			if e.AlignmentSpan < 0 {
				t.Fatalf("parsed a negative alignment span %d", e.AlignmentSpan)
			}
			if e.ContainerOffset < 0 || e.SliceOffset < 0 || e.SliceSize < 0 {
				t.Fatalf("parsed a negative offset/size in %+v", e)
			}
		}
		_ = idx.Query(0, 0, 0)
		_ = idx.Query(-1, 0, 1<<40)
		_ = idx.Query(1<<20, -5, 5)
	})
}

// compressionHeaderPayloads extracts every compression-header block
// payload from a CRAM file, for use as fuzz corpus seeds.
func compressionHeaderPayloads(data []byte) [][]byte {
	var out [][]byte
	rd, err := NewReader(bytes.NewReader(data))
	if err != nil {
		return out
	}
	for {
		c, err := rd.Next()
		if err == io.EOF || err != nil {
			return out
		}
		for bi := range c.Blocks {
			if c.Blocks[bi].ContentType != ContentCompressionHeader {
				continue
			}
			if p, derr := c.Blocks[bi].Decompress(); derr == nil {
				out = append(out, p)
			}
		}
	}
}
