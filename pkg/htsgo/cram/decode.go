package cram

import (
	"fmt"
)

// byteCursor is a forward-only reader over an external block's
// uncompressed payload. CRAM external data series are laid out
// back-to-back in their block, so a series decoder pulls values from a
// cursor that tracks its position; variable-length integers and raw byte
// runs are both read through it. The integer format is ITF-8 for CRAM
// v2/v3 and uint7 varints for v4, selected by the cursor's intReader.
type byteCursor struct {
	data []byte
	pos  int
	r    intReader
}

// remaining reports how many unread bytes the cursor still has.
func (c *byteCursor) remaining() int { return len(c.data) - c.pos }

// readByte reads and returns the next byte. It returns an error when the
// cursor is exhausted.
func (c *byteCursor) readByte() (byte, error) {
	if c.pos >= len(c.data) {
		return 0, fmt.Errorf("cram: external block exhausted")
	}
	b := c.data[c.pos]
	c.pos++
	return b, nil
}

// readN reads and returns the next n bytes as a fresh slice. It returns
// an error if fewer than n bytes remain.
func (c *byteCursor) readN(n int) ([]byte, error) {
	if n < 0 {
		return nil, fmt.Errorf("cram: external block read of negative length %d", n)
	}
	if c.pos+n > len(c.data) {
		return nil, fmt.Errorf("cram: external block exhausted: need %d bytes, %d remain", n, c.remaining())
	}
	out := append([]byte(nil), c.data[c.pos:c.pos+n]...)
	c.pos += n
	return out, nil
}

// readInt reads one version-aware unsigned 32-bit integer from the
// cursor (ITF-8 for v2/v3, uint7 for v4), advancing past it. It is the
// EXTERNAL-integer read used by the v2/v3 integer data series; v4 stores
// integers through the VARINT codecs instead (see readVarintInt /
// readVarintLong).
func (c *byteCursor) readInt() (int32, error) {
	v, n, err := c.r.u32(c.data, c.pos)
	if err != nil {
		return 0, err
	}
	c.pos += n
	return v, nil
}

// readVarintInt reads one v4 uint7 32-bit value (unsigned or, when signed
// is true, zig-zag signed) from the cursor, advancing past it. It backs
// the VARINT_UNSIGNED / VARINT_SIGNED integer codecs, whose offset the
// caller adds.
func (c *byteCursor) readVarintInt(signed bool) (int32, error) {
	var v int32
	var n int
	var err error
	if signed {
		v, n, err = sint7At32(c.data, c.pos)
	} else {
		v, n, err = uint7At32(c.data, c.pos)
	}
	if err != nil {
		return 0, err
	}
	c.pos += n
	return v, nil
}

// exhausted reports whether the cursor has consumed every byte.
func (c *byteCursor) exhausted() bool { return c.pos >= len(c.data) }

// seriesSource bundles the inputs a single data series decode needs: the
// shared CORE-block bitstream of the slice and a way to obtain the
// byteCursor for an external block keyed by content id. Each external
// block has exactly one cursor per decode pass so successive series
// reading the same block continue where the previous one stopped.
type seriesSource struct {
	core *bitReader
	// external lazily creates and caches one byteCursor per content id.
	external map[int32]*byteCursor
	// blocks maps a content id to the external block's decompressed
	// payload.
	blocks map[int32][]byte
	// reader is the version-aware integer reader handed to every cursor so
	// EXTERNAL integer values are read as ITF-8 (v2/v3) or uint7 (v4).
	reader intReader
}

// totalBytes returns the combined size of the slice's CORE bitstream and
// every external block. It bounds a per-record allocation: a single
// record's read length, base run or quality run cannot legitimately
// exceed the bytes of the blocks that encode the whole slice, so a
// declared length larger than this is corrupt. The result is clamped to
// the int32 range so it can bound an int32 length directly.
func (s *seriesSource) totalBytes() int32 {
	total := int64(len(s.core.data))
	for _, b := range s.blocks {
		total += int64(len(b))
	}
	if total > int64(^uint32(0)>>1) {
		return int32(^uint32(0) >> 1)
	}
	return int32(total)
}

// consumed returns a monotonically increasing measure of how much
// series input has been read: the CORE bitstream's bit position plus
// every external cursor's byte position. A per-record or per-feature
// decode loop compares it across iterations — an iteration that
// advances it by zero means the declared count has outrun the data, so
// the loop can stop instead of producing unbounded zero-byte items.
func (s *seriesSource) consumed() int64 {
	total := int64(s.core.pos)*8 - int64(s.core.nb)
	for _, c := range s.external {
		total += int64(c.pos) * 8
	}
	return total
}

// hasBlock reports whether an external block with the given content id
// is present in the slice. A data series whose encoding names an absent
// block carries no values in this slice (CRAM omits a series' block
// when the series contributes nothing).
func (s *seriesSource) hasBlock(id int32) bool {
	_, ok := s.blocks[id]
	return ok
}

// cursor returns the byteCursor for the external block with the given
// content id, creating it from the decompressed block payload on first
// use. It returns an error when no external block carries that id.
func (s *seriesSource) cursor(id int32) (*byteCursor, error) {
	if c, ok := s.external[id]; ok {
		return c, nil
	}
	data, ok := s.blocks[id]
	if !ok {
		return nil, fmt.Errorf("cram: no external block with content id %d", id)
	}
	c := &byteCursor{data: data, r: s.reader}
	s.external[id] = c
	return c, nil
}

// decodeInts decodes n integer values of a data series through enc,
// pulling from the CORE bitstream and/or external blocks as the
// encoding dictates. It is used for the data series whose values are
// integers (the majority of CRAM series).
func (enc *Encoding) decodeInts(s *seriesSource, n int) ([]int32, error) {
	if n < 0 {
		return nil, fmt.Errorf("cram: cannot decode a negative count (%d) of values", n)
	}
	out := make([]int32, 0, n)
	for i := 0; i < n; i++ {
		v, err := enc.decodeInt(s)
		if err != nil {
			return out, fmt.Errorf("cram: %s value %d: %w", enc.idString(), i, err)
		}
		out = append(out, v)
	}
	return out, nil
}

// drainInts decodes integer values through enc until the external block
// the encoding draws from is exhausted, returning every value. It is
// defined only for the EXTERNAL encoding, whose values come from a
// self-delimiting external block and so can be drained without record
// counts. A bitstream encoding (HUFFMAN, BETA, …) has no per-series
// boundary in the shared CORE block, so draining one is not meaningful
// and returns an error.
func (enc *Encoding) drainInts(s *seriesSource) ([]int32, error) {
	if enc == nil || enc.ID == EncodingNull {
		return nil, nil
	}
	switch enc.ID {
	case EncodingExternal:
		c, err := s.cursor(enc.ExternalID)
		if err != nil {
			return nil, err
		}
		var out []int32
		for !c.exhausted() {
			v, err := c.readInt()
			if err != nil {
				return out, fmt.Errorf("cram: EXTERNAL value %d: %w", len(out), err)
			}
			out = append(out, v)
		}
		return out, nil
	case EncodingVarintUnsigned, EncodingVarintSigned:
		// CRAM v4: a varint series is a self-delimiting uint7 stream in its
		// external block, so it drains like EXTERNAL with the offset applied.
		c, err := s.cursor(enc.ExternalID)
		if err != nil {
			return nil, err
		}
		var out []int32
		for !c.exhausted() {
			v, err := c.readVarintInt(enc.ID == EncodingVarintSigned)
			if err != nil {
				return out, fmt.Errorf("cram: VARINT value %d: %w", len(out), err)
			}
			out = append(out, v+int32(enc.VarintOffset))
		}
		return out, nil
	default:
		return nil, fmt.Errorf("cram: %s is not drainable; only EXTERNAL and VARINT series have a self-delimiting block", enc.ID)
	}
}

// drainByteArrayLen decodes every value of a BYTE_ARRAY_LEN series whose
// length sub-encoding is an EXTERNAL block. It handles the two layouts
// CRAM permits:
//
//   - separate blocks: the length and values sub-encodings draw from
//     different external blocks, so all the lengths sit at the front of
//     the length block and the values are read from the values block;
//   - a shared block: the length and values sub-encodings name the same
//     content id, so the length and the bytes of each value are
//     interleaved — read one length, then that many value bytes, repeat
//     until the block is exhausted.
//
// A BYTE_ARRAY_LEN whose length encoding lives in the CORE bitstream
// cannot be drained — its value count has no self-delimiting boundary —
// and returns an error.
func (enc *Encoding) drainByteArrayLen(s *seriesSource) ([][]byte, error) {
	if enc == nil || enc.ID != EncodingByteArrayLen {
		return nil, fmt.Errorf("cram: drainByteArrayLen called on a non-BYTE_ARRAY_LEN encoding")
	}
	if enc.LenEnc == nil || enc.LenEnc.ID != EncodingExternal {
		return nil, fmt.Errorf("cram: BYTE_ARRAY_LEN length sub-encoding is not EXTERNAL; cannot drain without a record count")
	}
	if enc.ValEnc != nil && enc.ValEnc.ID == EncodingExternal &&
		enc.ValEnc.ExternalID == enc.LenEnc.ExternalID {
		// Shared block: length and value bytes are interleaved.
		c, err := s.cursor(enc.LenEnc.ExternalID)
		if err != nil {
			return nil, err
		}
		var out [][]byte
		for !c.exhausted() {
			ln, lerr := c.readInt()
			if lerr != nil {
				return out, fmt.Errorf("cram: BYTE_ARRAY_LEN value %d length: %w", len(out), lerr)
			}
			if ln < 0 {
				return out, fmt.Errorf("cram: BYTE_ARRAY_LEN value %d has negative length %d", len(out), ln)
			}
			v, verr := c.readN(int(ln))
			if verr != nil {
				return out, fmt.Errorf("cram: BYTE_ARRAY_LEN value %d (%d bytes): %w", len(out), ln, verr)
			}
			out = append(out, v)
		}
		return out, nil
	}
	// Separate blocks: drain all the lengths, then read the values.
	lengths, err := enc.LenEnc.drainInts(s)
	if err != nil {
		return nil, fmt.Errorf("cram: BYTE_ARRAY_LEN lengths: %w", err)
	}
	out := make([][]byte, 0, len(lengths))
	for i, ln := range lengths {
		if ln < 0 {
			return out, fmt.Errorf("cram: BYTE_ARRAY_LEN value %d has negative length %d", i, ln)
		}
		v, err := enc.ValEnc.decodeRawBytes(s, int(ln))
		if err != nil {
			return out, fmt.Errorf("cram: BYTE_ARRAY_LEN value %d (%d bytes): %w", i, ln, err)
		}
		out = append(out, v)
	}
	return out, nil
}

// drainRawBytes returns every byte of an EXTERNAL series' block as
// individual single-byte values. It is the drain for a byte-valued
// EXTERNAL series (a quality-score or base series), whose block stores
// one verbatim byte per value rather than ITF-8 integers.
func (enc *Encoding) drainRawBytes(s *seriesSource) ([]byte, error) {
	if enc == nil || enc.ID == EncodingNull {
		return nil, nil
	}
	switch enc.ID {
	case EncodingExternal:
		c, err := s.cursor(enc.ExternalID)
		if err != nil {
			return nil, err
		}
		return c.readN(c.remaining())
	case EncodingXPack, EncodingXRLE:
		// CRAM v4 XPACK / XRLE expand a whole external block at once, so the
		// full expanded output is the drained byte series.
		out, err := enc.transformBytes(s)
		if err != nil {
			return nil, err
		}
		// Return the bytes not yet served positionally, draining the rest.
		tb := enc.transform
		rest := append([]byte(nil), out[tb.pos:]...)
		tb.pos = len(out)
		return rest, nil
	default:
		return nil, fmt.Errorf("cram: %s is not a raw-byte EXTERNAL or transform series", enc.ID)
	}
}

// drainByteArrays decodes byte-array values through enc until the
// external block the encoding draws from is exhausted. It is defined for
// BYTE_ARRAY_STOP, whose values are each delimited by a stop byte.
// BYTE_ARRAY_LEN cannot be drained — its length sub-encoding may live in
// the CORE bitstream — so it returns an error.
func (enc *Encoding) drainByteArrays(s *seriesSource) ([][]byte, error) {
	if enc == nil || enc.ID == EncodingNull {
		return nil, nil
	}
	if enc.ID != EncodingByteArrayStop {
		return nil, fmt.Errorf("cram: %s is not drainable; only BYTE_ARRAY_STOP has self-delimiting values", enc.ID)
	}
	c, err := s.cursor(enc.ExternalID)
	if err != nil {
		return nil, err
	}
	var out [][]byte
	for !c.exhausted() {
		var v []byte
		for {
			b, rerr := c.readByte()
			if rerr != nil {
				return out, fmt.Errorf("cram: BYTE_ARRAY_STOP array %d ran off the end before its stop byte", len(out))
			}
			if b == enc.StopByte {
				break
			}
			v = append(v, b)
		}
		out = append(out, v)
	}
	return out, nil
}

// idString returns the encoding's identifier name, treating a nil
// encoding as NULL for error messages.
func (enc *Encoding) idString() string {
	if enc == nil {
		return EncodingNull.String()
	}
	return enc.ID.String()
}

// decodeInt decodes a single integer value through the encoding.
func (enc *Encoding) decodeInt(s *seriesSource) (int32, error) {
	if enc == nil || enc.ID == EncodingNull {
		return 0, fmt.Errorf("cram: NULL encoding yields no integer value")
	}
	switch enc.ID {
	case EncodingExternal:
		c, err := s.cursor(enc.ExternalID)
		if err != nil {
			return 0, err
		}
		return c.readInt()
	case EncodingVarintUnsigned, EncodingVarintSigned:
		// CRAM v4: integers are stored as uint7 varints in an external
		// block, with a constant offset added. VARINT_SIGNED applies the
		// zig-zag transform before adding the offset.
		c, err := s.cursor(enc.ExternalID)
		if err != nil {
			return 0, err
		}
		v, err := c.readVarintInt(enc.ID == EncodingVarintSigned)
		if err != nil {
			return 0, err
		}
		return v + int32(enc.VarintOffset), nil
	case EncodingConstByte, EncodingConstInt:
		// CRAM v4: every value is the same constant; no block is read.
		return int32(enc.ConstValue), nil
	case EncodingXDelta:
		// CRAM v4 XDELTA value-at-a-time integer path
		// (cram_xdelta_decode_int): decode one delta through the sub-codec,
		// un-zig-zag it, and prefix-sum onto the running previous value.
		raw, err := enc.SubEnc.decodeInt(s)
		if err != nil {
			return 0, fmt.Errorf("cram: XDELTA delta: %w", err)
		}
		enc.deltaLast += int32(unzigzag32(uint32(raw)))
		return enc.deltaLast, nil
	case EncodingXPack, EncodingXRLE:
		// CRAM v4 XPACK / XRLE wrap a byte-valued series; an integer read
		// yields the next expanded byte (cram_x*_decode_char serves bytes
		// from the expanded block).
		b, err := enc.readTransform(s, 1)
		if err != nil {
			return 0, err
		}
		return int32(b[0]), nil
	case EncodingHuffman:
		t, err := enc.huffmanDecoder()
		if err != nil {
			return 0, err
		}
		return t.decode(s.core)
	case EncodingBeta:
		v, err := s.core.readBits(uint(enc.NumBits))
		if err != nil {
			return 0, err
		}
		return int32(v) - enc.Offset, nil
	case EncodingGamma:
		v, err := readGamma(s.core)
		if err != nil {
			return 0, err
		}
		return v - enc.Offset, nil
	case EncodingSubexp:
		v, err := readSubexp(s.core, int(enc.K))
		if err != nil {
			return 0, err
		}
		return v - enc.Offset, nil
	case EncodingGolomb:
		v, err := readGolomb(s.core, int(enc.M))
		if err != nil {
			return 0, err
		}
		return v - enc.Offset, nil
	case EncodingGolombRice:
		// GOLOMB_RICE stores log2(M); the divisor is 1 << K.
		v, err := readGolomb(s.core, 1<<uint(enc.K))
		if err != nil {
			return 0, err
		}
		return v - enc.Offset, nil
	case EncodingByteArrayLen, EncodingByteArrayStop:
		return 0, fmt.Errorf("cram: %s is a byte-array encoding, not an integer encoding", enc.ID)
	default:
		return 0, fmt.Errorf("cram: %s cannot decode an integer value", enc.ID)
	}
}

// huffmanDecoder lazily builds and caches the canonical Huffman decoder
// for a HUFFMAN encoding.
func (enc *Encoding) huffmanDecoder() (*huffmanTable, error) {
	if enc.huffman != nil {
		return enc.huffman, nil
	}
	t, err := newHuffmanTable(enc.Symbols, enc.BitLengths)
	if err != nil {
		return nil, err
	}
	enc.huffman = t
	return t, nil
}

// decodeByteArrays decodes n byte-array values of a data series through
// enc. It is used for the data series whose values are variable-length
// byte strings (read names, quality strings, tag values …). Only the
// BYTE_ARRAY_LEN and BYTE_ARRAY_STOP encodings — and EXTERNAL, read as
// raw bytes — produce byte arrays.
func (enc *Encoding) decodeByteArrays(s *seriesSource, n int) ([][]byte, error) {
	if n < 0 {
		return nil, fmt.Errorf("cram: cannot decode a negative count (%d) of byte arrays", n)
	}
	out := make([][]byte, 0, n)
	for i := 0; i < n; i++ {
		v, err := enc.decodeByteArray(s)
		if err != nil {
			return out, fmt.Errorf("cram: %s byte array %d: %w", enc.idString(), i, err)
		}
		out = append(out, v)
	}
	return out, nil
}

// decodeByteArray decodes a single variable-length byte array through
// the encoding.
func (enc *Encoding) decodeByteArray(s *seriesSource) ([]byte, error) {
	if enc == nil || enc.ID == EncodingNull {
		return nil, fmt.Errorf("cram: NULL encoding yields no byte array")
	}
	switch enc.ID {
	case EncodingByteArrayStop:
		c, err := s.cursor(enc.ExternalID)
		if err != nil {
			return nil, err
		}
		var out []byte
		for {
			b, err := c.readByte()
			if err != nil {
				return nil, fmt.Errorf("cram: BYTE_ARRAY_STOP ran off the end before its stop byte %#02x", enc.StopByte)
			}
			if b == enc.StopByte {
				return out, nil
			}
			out = append(out, b)
		}
	case EncodingByteArrayLen:
		// The length sub-encoding gives the array's length; the values
		// sub-encoding then yields that many bytes.
		ln, err := enc.LenEnc.decodeInt(s)
		if err != nil {
			return nil, fmt.Errorf("cram: BYTE_ARRAY_LEN length: %w", err)
		}
		if ln < 0 {
			return nil, fmt.Errorf("cram: BYTE_ARRAY_LEN decoded a negative length %d", ln)
		}
		return enc.ValEnc.decodeRawBytes(s, int(ln))
	case EncodingExternal:
		// An EXTERNAL series carrying raw bytes (rather than ITF-8 ints)
		// reads a single byte per value when used directly; callers that
		// want a fixed-length run use decodeRawBytes.
		c, err := s.cursor(enc.ExternalID)
		if err != nil {
			return nil, err
		}
		b, err := c.readByte()
		if err != nil {
			return nil, err
		}
		return []byte{b}, nil
	case EncodingXPack, EncodingXRLE:
		// A transform-wrapped byte series read one value at a time yields a
		// single expanded byte.
		return enc.readTransform(s, 1)
	case EncodingXDelta:
		v, err := enc.decodeInt(s)
		if err != nil {
			return nil, err
		}
		return []byte{byte(v)}, nil
	default:
		return nil, fmt.Errorf("cram: %s cannot decode a byte array", enc.ID)
	}
}

// decodeRawBytes reads exactly n raw bytes through the encoding. It is
// the values half of a BYTE_ARRAY_LEN decode: the values sub-encoding
// is typically an EXTERNAL byte block, but a HUFFMAN-of-bytes or BETA
// is also legal, so each byte is decoded through the sub-encoding.
func (enc *Encoding) decodeRawBytes(s *seriesSource, n int) ([]byte, error) {
	if enc == nil || enc.ID == EncodingNull {
		if n == 0 {
			return nil, nil
		}
		return nil, fmt.Errorf("cram: NULL encoding cannot supply %d bytes", n)
	}
	switch enc.ID {
	case EncodingExternal:
		c, err := s.cursor(enc.ExternalID)
		if err != nil {
			return nil, err
		}
		return c.readN(n)
	case EncodingByteArrayStop:
		// A stop-delimited values stream still yields fixed-length runs
		// when the caller knows the length; read n bytes verbatim.
		c, err := s.cursor(enc.ExternalID)
		if err != nil {
			return nil, err
		}
		return c.readN(n)
	case EncodingXPack, EncodingXRLE:
		// CRAM v4 XPACK / XRLE serve their expanded bytes positionally.
		return enc.readTransform(s, n)
	case EncodingXDelta:
		// CRAM v4 XDELTA over a byte series reconstructs each byte through
		// the value-at-a-time prefix-sum path (cram_xdelta_decode_int); the
		// 2-byte block path is reached via decodeByteArrayBlock instead.
		out := make([]byte, n)
		for i := 0; i < n; i++ {
			v, err := enc.decodeInt(s)
			if err != nil {
				return nil, err
			}
			out[i] = byte(v)
		}
		return out, nil
	case EncodingHuffman, EncodingBeta:
		// Byte-valued sub-encoding: each decoded int is a single byte.
		// A value outside 0-255 is truncated; callers use this only for
		// series the spec defines as byte-valued.
		out := make([]byte, n)
		for i := 0; i < n; i++ {
			v, err := enc.decodeInt(s)
			if err != nil {
				return nil, err
			}
			out[i] = byte(v)
		}
		return out, nil
	default:
		return nil, fmt.Errorf("cram: %s cannot supply raw bytes", enc.ID)
	}
}
