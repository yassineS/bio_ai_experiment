package cram

import (
	"bytes"
	"compress/gzip"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram/codec"
)

// v3Def is a CRAM v3.0 file definition value for tests that need a
// FileDefinition without a real file.
var v3Def = FileDefinition{Major: 3, Minor: 0}

// gzipBytes gzip-compresses b, for building synthetic gzip blocks.
func gzipBytes(t *testing.T, b []byte) []byte {
	t.Helper()
	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	if _, err := zw.Write(b); err != nil {
		t.Fatalf("gzip write: %v", err)
	}
	if err := zw.Close(); err != nil {
		t.Fatalf("gzip close: %v", err)
	}
	return buf.Bytes()
}

// TestCompressionMethodString exercises the CompressionMethod.String
// mapping including the unknown fallthrough.
func TestCompressionMethodString(t *testing.T) {
	cases := map[CompressionMethod]string{
		CompRaw:      "raw",
		CompGzip:     "gzip",
		CompBzip2:    "bzip2",
		CompLZMA:     "lzma",
		CompRANS4x8:  "rans4x8",
		CompRANS4x16: "rans4x16",
		CompArith:    "arith",
		CompFQZComp:  "fqzcomp",
		CompNameTok:  "name-tokeniser",
		200:          "unknown(200)",
	}
	for m, want := range cases {
		if got := m.String(); got != want {
			t.Errorf("CompressionMethod(%d).String() = %q, want %q", byte(m), got, want)
		}
	}
}

// TestBlockContentTypeString exercises the BlockContentType.String
// mapping including the unknown fallthrough.
func TestBlockContentTypeString(t *testing.T) {
	cases := map[BlockContentType]string{
		ContentFileHeader:        "file-header",
		ContentCompressionHeader: "compression-header",
		ContentMappedSlice:       "slice-header",
		ContentReserved:          "reserved",
		ContentExternal:          "external",
		ContentCoreData:          "core-data",
		99:                       "unknown(99)",
	}
	for ct, want := range cases {
		if got := ct.String(); got != want {
			t.Errorf("BlockContentType(%d).String() = %q, want %q", byte(ct), got, want)
		}
	}
}

// TestDecompressRaw checks a raw block returns its data unchanged.
func TestDecompressRaw(t *testing.T) {
	payload := []byte("hello cram")
	b := Block{Method: CompRaw, Data: payload, UncompressedSize: int32(len(payload))}
	out, err := b.Decompress()
	if err != nil {
		t.Fatalf("Decompress: %v", err)
	}
	if !bytes.Equal(out, payload) {
		t.Errorf("raw block: got %q, want %q", out, payload)
	}
}

// TestDecompressGzip checks a gzip block round-trips.
func TestDecompressGzip(t *testing.T) {
	payload := bytes.Repeat([]byte("ACGT"), 100)
	comp := gzipBytes(t, payload)
	b := Block{Method: CompGzip, Data: comp, UncompressedSize: int32(len(payload))}
	out, err := b.Decompress()
	if err != nil {
		t.Fatalf("Decompress: %v", err)
	}
	if !bytes.Equal(out, payload) {
		t.Errorf("gzip block did not round-trip")
	}
}

// TestDecompressSizeMismatch checks that a block whose decompressed
// length disagrees with its declared UncompressedSize is rejected.
func TestDecompressSizeMismatch(t *testing.T) {
	payload := []byte("0123456789")
	b := Block{Method: CompRaw, Data: payload, UncompressedSize: 999}
	if _, err := b.Decompress(); err == nil {
		t.Errorf("expected a size-mismatch error")
	}
}

// TestDecompressUnsupportedMethods checks every out-of-scope compression
// method returns a clear unsupported-method error rather than panicking.
func TestDecompressUnsupportedMethods(t *testing.T) {
	for _, m := range []CompressionMethod{250} {
		b := Block{Method: m, Data: []byte{1, 2, 3}}
		_, err := b.Decompress()
		if err == nil {
			t.Errorf("method %d: expected unsupported-method error", byte(m))
			continue
		}
		if !bytes.Contains([]byte(err.Error()), []byte("unsupported compression method")) {
			t.Errorf("method %d: error %q should mention unsupported method", byte(m), err)
		}
	}
}

// TestDecompressCorruptGzip checks a malformed gzip payload produces an
// error, not a panic.
func TestDecompressCorruptGzip(t *testing.T) {
	b := Block{Method: CompGzip, Data: []byte("not gzip at all"), UncompressedSize: 10}
	if _, err := b.Decompress(); err == nil {
		t.Errorf("expected error decompressing corrupt gzip")
	}
}

// TestSupportedMethod pins the SupportedMethod predicate.
func TestSupportedMethod(t *testing.T) {
	supported := []CompressionMethod{CompRaw, CompGzip, CompBzip2, CompLZMA, CompRANS4x8, CompRANS4x16, CompArith, CompFQZComp, CompNameTok}
	for _, m := range supported {
		if !(&Block{Method: m}).SupportedMethod() {
			t.Errorf("method %s should be supported", m)
		}
	}
	for _, m := range []CompressionMethod{250} {
		if (&Block{Method: m}).SupportedMethod() {
			t.Errorf("method %d should not be supported", byte(m))
		}
	}
}

// buildBlock assembles a CRAM v3 block on the wire: header, data and the
// trailing CRC32.
func buildBlock(method CompressionMethod, ct BlockContentType, cid int32, data []byte, uncompressed int32) []byte {
	var buf bytes.Buffer
	buf.WriteByte(byte(method))
	buf.WriteByte(byte(ct))
	buf.Write(encITF8(cid))
	buf.Write(encITF8(int32(len(data))))
	buf.Write(encITF8(uncompressed))
	buf.Write(data)
	return appendCRC(buf.Bytes())
}

// TestReadBlockRoundTrip builds a v3 block and reads it back, checking
// the CRC32 validates.
func TestReadBlockRoundTrip(t *testing.T) {
	data := []byte("payload-bytes")
	wire := buildBlock(CompRaw, ContentExternal, 42, data, int32(len(data)))
	b, err := readBlock(bytes.NewReader(wire), v3Def, 1<<20)
	if err != nil {
		t.Fatalf("readBlock: %v", err)
	}
	if b.Method != CompRaw || b.ContentType != ContentExternal || b.ContentID != 42 {
		t.Errorf("block header fields wrong: %+v", b)
	}
	if !bytes.Equal(b.Data, data) {
		t.Errorf("block data mismatch")
	}
}

// TestReadBlockCRCMismatch flips a CRC byte and checks readBlock rejects
// it.
func TestReadBlockCRCMismatch(t *testing.T) {
	wire := buildBlock(CompRaw, ContentExternal, 1, []byte("abc"), 3)
	wire[len(wire)-1] ^= 0xff
	if _, err := readBlock(bytes.NewReader(wire), v3Def, 1<<20); err == nil {
		t.Errorf("expected a CRC32 mismatch error")
	}
}

// TestReadBlockTruncated checks readBlock errors — never panics — on
// inputs truncated at each stage of the block.
func TestReadBlockTruncated(t *testing.T) {
	full := buildBlock(CompRaw, ContentExternal, 1, []byte("abcdef"), 6)
	for n := 0; n < len(full); n++ {
		if _, err := readBlock(bytes.NewReader(full[:n]), v3Def, 1<<20); err == nil {
			t.Errorf("expected error for block truncated to %d bytes", n)
		}
	}
}

// TestDecompressRANS4x8 checks the rANS 4x8 dispatch path round-trips
// through a block, using the codec's own encoder to build the payload.
func TestDecompressRANS4x8(t *testing.T) {
	payload := bytes.Repeat([]byte("ACGTACGTNN"), 200)
	comp, err := codec.RANS4x8Encode(payload, 0)
	if err != nil {
		t.Fatalf("RANS4x8Encode: %v", err)
	}
	b := Block{Method: CompRANS4x8, Data: comp, UncompressedSize: int32(len(payload))}
	out, err := b.Decompress()
	if err != nil {
		t.Fatalf("Decompress rANS4x8: %v", err)
	}
	if !bytes.Equal(out, payload) {
		t.Errorf("rANS4x8 block did not round-trip")
	}
}

// TestDecompressRANS4x16 checks the rANS 4x16 dispatch path round-trips
// through a block.
func TestDecompressRANS4x16(t *testing.T) {
	payload := bytes.Repeat([]byte("the quick brown fox "), 300)
	comp, err := codec.RANS4x16Encode(payload, 0)
	if err != nil {
		t.Fatalf("RANS4x16Encode: %v", err)
	}
	b := Block{Method: CompRANS4x16, Data: comp, UncompressedSize: int32(len(payload))}
	out, err := b.Decompress()
	if err != nil {
		t.Fatalf("Decompress rANS4x16: %v", err)
	}
	if !bytes.Equal(out, payload) {
		t.Errorf("rANS4x16 block did not round-trip")
	}
}

// TestDecompressLZMA checks the LZMA (method 3) dispatch path round-trips
// through a block, using the codec's own .xz encoder to build the payload.
func TestDecompressLZMA(t *testing.T) {
	payload := bytes.Repeat([]byte("GATTACA reference sequence "), 400)
	comp, err := codec.LZMAEncode(payload)
	if err != nil {
		t.Fatalf("LZMAEncode: %v", err)
	}
	b := Block{Method: CompLZMA, Data: comp, UncompressedSize: int32(len(payload))}
	out, err := b.Decompress()
	if err != nil {
		t.Fatalf("Decompress LZMA: %v", err)
	}
	if !bytes.Equal(out, payload) {
		t.Errorf("LZMA block did not round-trip")
	}
}

// TestDecompressFQZComp checks the fqzcomp (method 7) dispatch path
// round-trips through a block, using the codec's own encoder to build
// the payload. The fqzcomp decoded length is verified against the
// block's declared UncompressedSize like every other method.
func TestDecompressFQZComp(t *testing.T) {
	payload := bytes.Repeat([]byte{30, 30, 35, 2, 30, 35, 35, 2}, 400)
	for strat := 0; strat <= 3; strat++ {
		comp, err := codec.FQZCompEncode(payload, strat, nil)
		if err != nil {
			t.Fatalf("FQZCompEncode strat %d: %v", strat, err)
		}
		b := Block{Method: CompFQZComp, Data: comp, UncompressedSize: int32(len(payload))}
		out, err := b.Decompress()
		if err != nil {
			t.Fatalf("Decompress fqzcomp strat %d: %v", strat, err)
		}
		if !bytes.Equal(out, payload) {
			t.Errorf("fqzcomp block did not round-trip at strat %d", strat)
		}
	}
}

// TestDecompressFQZCompSizeMismatch checks an fqzcomp block whose
// decompressed length disagrees with its declared UncompressedSize is
// rejected.
func TestDecompressFQZCompSizeMismatch(t *testing.T) {
	comp, err := codec.FQZCompEncode([]byte{10, 11, 12, 13}, 0, nil)
	if err != nil {
		t.Fatalf("FQZCompEncode: %v", err)
	}
	b := Block{Method: CompFQZComp, Data: comp, UncompressedSize: 99}
	if _, err := b.Decompress(); err == nil {
		t.Errorf("expected size-mismatch error for fqzcomp block")
	}
}

// TestDecompressLZMASizeMismatch checks an LZMA block whose decompressed
// length disagrees with its declared UncompressedSize is rejected.
func TestDecompressLZMASizeMismatch(t *testing.T) {
	comp, err := codec.LZMAEncode([]byte("twelve bytes"))
	if err != nil {
		t.Fatalf("LZMAEncode: %v", err)
	}
	b := Block{Method: CompLZMA, Data: comp, UncompressedSize: 99}
	if _, err := b.Decompress(); err == nil {
		t.Errorf("expected size-mismatch error for LZMA block")
	}
}

// TestDecompressCorruptLZMA checks a malformed LZMA payload produces an
// error, not a panic.
func TestDecompressCorruptLZMA(t *testing.T) {
	b := Block{Method: CompLZMA, Data: []byte("not an xz stream"), UncompressedSize: 10}
	if _, err := b.Decompress(); err == nil {
		t.Errorf("expected error decompressing corrupt LZMA block")
	}
}

// TestReadBlockV2NoCRC checks a CRAM v2 block (no trailing CRC32) parses
// when the FileDefinition reports a v2 version.
func TestReadBlockV2NoCRC(t *testing.T) {
	v2Def := FileDefinition{Major: 2, Minor: 1}
	var buf bytes.Buffer
	buf.WriteByte(byte(CompRaw))
	buf.WriteByte(byte(ContentExternal))
	buf.Write(encITF8(3))
	data := []byte("v2data")
	buf.Write(encITF8(int32(len(data))))
	buf.Write(encITF8(int32(len(data))))
	buf.Write(data)
	b, err := readBlock(bytes.NewReader(buf.Bytes()), v2Def, 1<<20)
	if err != nil {
		t.Fatalf("readBlock v2: %v", err)
	}
	if !bytes.Equal(b.Data, data) || b.CRC != 0 {
		t.Errorf("v2 block parsed wrong: %+v", b)
	}
}

// TestReadBlockSizeExceedsContainer checks a block declaring a
// compressed size larger than the bytes left in the container is
// rejected before any large allocation.
func TestReadBlockSizeExceedsContainer(t *testing.T) {
	wire := buildBlock(CompRaw, ContentExternal, 1, []byte("abcdef"), 6)
	if _, err := readBlock(bytes.NewReader(wire), v3Def, 3); err == nil {
		t.Errorf("expected error when compressed size exceeds container budget")
	}
}

// TestReadBlockNegativeSize checks a block declaring a negative size is
// rejected rather than triggering a huge or panicking allocation.
func TestReadBlockNegativeSize(t *testing.T) {
	var buf bytes.Buffer
	buf.WriteByte(byte(CompRaw))
	buf.WriteByte(byte(ContentExternal))
	buf.Write(encITF8(0))
	// -1 as ITF-8.
	buf.Write([]byte{0xff, 0xff, 0xff, 0xff, 0x0f})
	buf.Write(encITF8(0))
	if _, err := readBlock(bytes.NewReader(buf.Bytes()), v3Def, 1<<20); err == nil {
		t.Errorf("expected error for negative compressed size")
	}
}

// TestDecompressNameTok checks the name-tokeniser (method 8) dispatch
// path: a block built with the codec's encoder decompresses through
// Block.Decompress to the NUL-joined read names.
func TestDecompressNameTok(t *testing.T) {
	names := []byte("HS25_09827:2:2215:4133:22216#49\n" +
		"HS25_09827:2:1212:15822:94146#49\n" +
		"HS25_09827:2:1209:9304:17097#49\n")
	want := bytes.ReplaceAll(names, []byte{'\n'}, []byte{0})
	for _, lvl := range []int{1, 9, 11, 19} {
		comp, err := codec.NameTokEncode(names, lvl)
		if err != nil {
			t.Fatalf("NameTokEncode L%d: %v", lvl, err)
		}
		b := Block{Method: CompNameTok, Data: comp, UncompressedSize: int32(len(want))}
		out, err := b.Decompress()
		if err != nil {
			t.Fatalf("Decompress name-tokeniser L%d: %v", lvl, err)
		}
		if !bytes.Equal(out, want) {
			t.Errorf("name-tokeniser block L%d did not round-trip", lvl)
		}
	}
}
