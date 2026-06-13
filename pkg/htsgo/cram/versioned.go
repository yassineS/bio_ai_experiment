package cram

// intReader is a version-aware decoder for the variable-length integers
// that pervade a CRAM container: the container header, block headers,
// slice header, compression-header maps and per-encoding parameters all
// read a stream of such integers. For CRAM v2.x/v3.x it dispatches to the
// ITF-8 / LTF-8 readers (itf8.go); for v4.0 it dispatches to the uint7
// varint readers (varint.go). Threading one intReader through the parsers
// keeps the v3 byte-for-byte path intact while letting v4 reinterpret
// every field, exactly as htslib's per-version varint_vec function table
// does (cram_io.c cram_init_varint).
type intReader struct {
	// uint7 selects the v4 uint7 path when true, the v3 ITF-8/LTF-8 path
	// when false.
	uint7 bool
	// majorVersion is the CRAM major version the reader was built for. It
	// is recorded onto each parsed Encoding so the decode side can pick
	// the matching integer format for external-block reads.
	majorVersion uint8
}

// newIntReader returns an intReader for the given CRAM major version.
func newIntReader(major uint8) intReader {
	return intReader{uint7: major >= 4, majorVersion: major}
}

// major returns the CRAM major version this reader decodes for.
func (r intReader) major() uint8 { return r.majorVersion }

// u32 decodes an unsigned 32-bit field at off within p, the version-aware
// analogue of itf8At. It returns the value and the number of bytes
// consumed. The v3 ITF-8 form already round-trips small negatives (it
// sign-extends), so a v3 caller reading a value it treats as signed (a
// ref id) reuses this; v4 callers that need zig-zag signedness call s32.
func (r intReader) u32(p []byte, off int) (int32, int, error) {
	if r.uint7 {
		return uint7At32(p, off)
	}
	return itf8At(p, off)
}

// s32 decodes a signed 32-bit field at off within p. For v4 it applies the
// uint7 zig-zag transform; for v3 it falls back to the ITF-8 reader, whose
// sign-extension of the top nibble already yields the correct signed value
// (this matches htslib mapping varint_get32s to safe_itf8_get for v3).
func (r intReader) s32(p []byte, off int) (int32, int, error) {
	if r.uint7 {
		return sint7At32(p, off)
	}
	return itf8At(p, off)
}

// u64 decodes an unsigned 64-bit field at off within p, the version-aware
// analogue of ltf8At.
func (r intReader) u64(p []byte, off int) (int64, int, error) {
	if r.uint7 {
		return uint7At64(p, off)
	}
	return ltf8At(p, off)
}

// s64 decodes a signed 64-bit field at off within p. For v4 it applies the
// uint7 zig-zag transform; for v3 it uses the LTF-8 reader.
func (r intReader) s64(p []byte, off int) (int64, int, error) {
	if r.uint7 {
		return sint7At64(p, off)
	}
	return ltf8At(p, off)
}

// intWriter is the encode-side mirror of intReader: a version-aware
// serialiser for the variable-length integers that frame a CRAM container.
// The v3 writer (uint7 false) appends ITF-8 / LTF-8 exactly as before; the
// v4.0 writer (uint7 true) appends uint7 LEB128 varints, with signed fields
// zig-zag coded. Threading one intWriter through the container, block,
// slice-header, compression-header and encoding serialisers keeps the v3
// on-disk bytes byte-for-byte unchanged while emitting valid v4 framing,
// the inverse of how intReader threads the decode side.
type intWriter struct {
	// uint7 selects the v4 uint7 put path when true, the v3 ITF-8/LTF-8
	// path when false.
	uint7 bool
}

// newIntWriter returns an intWriter for the given CRAM writer version.
func newIntWriter(v Version) intWriter { return intWriter{uint7: v.usesUint7()} }

// u32 appends an unsigned 32-bit field to dst. v4 uses the uint7 varint;
// v3 uses ITF-8 (whose sign-extension already round-trips a small negative
// a v3 caller treats as signed, e.g. a ref id, so v3 callers reuse this).
func (w intWriter) u32(dst []byte, v int32) []byte {
	if w.uint7 {
		return appendUint7(dst, uint64(uint32(v)))
	}
	return appendITF8(dst, v)
}

// s32 appends a signed 32-bit field to dst. v4 zig-zag encodes it as a
// uint7 varint; v3 falls back to ITF-8, whose sign-extension already
// yields the correct signed value (matching the decode-side intReader.s32).
func (w intWriter) s32(dst []byte, v int32) []byte {
	if w.uint7 {
		return appendSint7(dst, int64(v))
	}
	return appendITF8(dst, v)
}

// u64 appends an unsigned 64-bit field to dst (uint7 for v4, LTF-8 for v3).
func (w intWriter) u64(dst []byte, v int64) []byte {
	if w.uint7 {
		return appendUint7(dst, uint64(v))
	}
	return appendLTF8(dst, v)
}

// s64 appends a signed 64-bit field to dst (zig-zag uint7 for v4, LTF-8 for
// v3).
func (w intWriter) s64(dst []byte, v int64) []byte {
	if w.uint7 {
		return appendSint7(dst, v)
	}
	return appendLTF8(dst, v)
}
