package cram

import (
	"bytes"
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
		sh, err := parseSliceHeader(data)
		if err != nil {
			return
		}
		// A successful parse must produce a self-consistent header.
		if int(sh.NumBlocks) < 0 {
			t.Fatalf("parsed a negative block count %d", sh.NumBlocks)
		}
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
