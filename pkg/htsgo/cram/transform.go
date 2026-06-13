package cram

import (
	"fmt"
)

// transformBlock is the fully expanded byte output of a CRAM v4 transform
// codec (XPACK, XRLE, or the XDELTA block path). htslib expands a
// transform once into slice->block_by_id[512+codec_id] and then serves the
// decoded bytes positionally as records are decoded; this mirrors that with
// a cached buffer and a read cursor.
type transformBlock struct {
	data []byte
	pos  int
}

// transformBytes lazily expands enc's transform output (caching it on the
// encoding) and returns its byte buffer. It is the equivalent of htslib's
// cram_*_decode_expand_char / get_block: the transform's source is decoded
// through its sub-codec(s), then the transform is reversed wholesale.
func (enc *Encoding) transformBytes(s *seriesSource) ([]byte, error) {
	if enc.transform != nil {
		return enc.transform.data, nil
	}
	var (
		out []byte
		err error
	)
	switch enc.ID {
	case EncodingXPack:
		out, err = enc.expandXPack(s)
	case EncodingXRLE:
		out, err = enc.expandXRLE(s)
	case EncodingXDelta:
		out, err = enc.expandXDelta(s)
	default:
		return nil, fmt.Errorf("cram: %s is not a byte-producing transform codec", enc.ID)
	}
	if err != nil {
		return nil, err
	}
	enc.transform = &transformBlock{data: out}
	return out, nil
}

// readTransform serves the next n bytes of enc's expanded transform output,
// advancing the transform's read cursor. It is how a byte-valued series that
// is wrapped by a transform codec pulls a fixed-length run.
func (enc *Encoding) readTransform(s *seriesSource, n int) ([]byte, error) {
	if n < 0 {
		return nil, fmt.Errorf("cram: %s transform read of negative length %d", enc.ID, n)
	}
	if _, err := enc.transformBytes(s); err != nil {
		return nil, err
	}
	tb := enc.transform
	if tb.pos+n > len(tb.data) {
		return nil, fmt.Errorf("cram: %s transform exhausted: need %d bytes, %d remain",
			enc.ID, n, len(tb.data)-tb.pos)
	}
	out := append([]byte(nil), tb.data[tb.pos:tb.pos+n]...)
	tb.pos += n
	return out, nil
}

// subCodecBytes decodes one sub-codec of a transform to its complete byte
// output (htslib's sub_codec->get_block). A transform's sub-codec reads a
// self-delimiting source — an EXTERNAL block (the byte stream verbatim), a
// BYTE_ARRAY_STOP block, or a constant — so the whole block can be produced
// at once. Nested transforms are expanded recursively.
func subCodecBytes(sub *Encoding, s *seriesSource) ([]byte, error) {
	if sub == nil || sub.ID == EncodingNull {
		return nil, nil
	}
	switch sub.ID {
	case EncodingExternal:
		c, err := s.cursor(sub.ExternalID)
		if err != nil {
			return nil, err
		}
		return c.readN(c.remaining())
	case EncodingByteArrayStop:
		c, err := s.cursor(sub.ExternalID)
		if err != nil {
			return nil, err
		}
		return c.readN(c.remaining())
	case EncodingConstByte:
		// A constant-byte sub-codec has no source block; htslib never pairs
		// one with a transform, but treat it as an empty stream rather than
		// erroring so a degenerate header decodes predictably.
		return nil, nil
	case EncodingXPack, EncodingXRLE, EncodingXDelta:
		return sub.transformBytes(s)
	default:
		return nil, fmt.Errorf("cram: %s sub-codec of a transform does not produce a self-delimiting byte block", sub.ID)
	}
}

// expandXPack reverses an XPACK transform (cram_codecs.c
// cram_xpack_decode_expand_char + htscodecs hts_unpack): the sub-codec
// supplies packed bytes, each holding 8/nbits symbols stored
// least-significant-first, and each symbol indexes the reverse map to its
// real byte value. nbits == 0 (a single-symbol alphabet) means every output
// byte is map[0].
func (enc *Encoding) expandXPack(s *seriesSource) ([]byte, error) {
	if len(enc.PackMap) == 0 {
		return nil, fmt.Errorf("cram: XPACK encoding has an empty map")
	}
	src, err := subCodecBytes(enc.SubEnc, s)
	if err != nil {
		return nil, fmt.Errorf("cram: XPACK sub-codec: %w", err)
	}
	if enc.PackBits == 0 {
		// Degenerate alphabet: nval <= 1, so every value is map[0]. There is
		// one packed byte per output byte (hts_pack stores them 1:1).
		out := make([]byte, len(src))
		for i := range out {
			out[i] = enc.PackMap[0]
		}
		return out, nil
	}
	nbits := int(enc.PackBits)
	per := 8 / nbits // symbols packed into each source byte
	mask := byte((1 << uint(nbits)) - 1)
	out := make([]byte, 0, len(src)*per)
	for _, b := range src {
		for k := 0; k < per; k++ {
			idx := (b >> uint(k*nbits)) & mask
			if int(idx) >= len(enc.PackMap) {
				return nil, fmt.Errorf("cram: XPACK symbol %d is outside the %d-entry map", idx, len(enc.PackMap))
			}
			out = append(out, enc.PackMap[idx])
		}
	}
	return out, nil
}

// expandXRLE reverses an XRLE transform (cram_codecs.c
// cram_xrle_decode_expand_char + htscodecs hts_rle_decode): the literal
// sub-codec supplies a literal byte stream, the length sub-codec a varint
// stream prefixed by the total output size; a literal whose byte is one of
// the RLE symbols consumes one run length L and expands to L+1 copies, and
// any other literal stands for a single byte.
func (enc *Encoding) expandXRLE(s *seriesSource) ([]byte, error) {
	lit, err := subCodecBytes(enc.LitSubEnc, s)
	if err != nil {
		return nil, fmt.Errorf("cram: XRLE literal sub-codec: %w", err)
	}
	lenBytes, err := subCodecBytes(enc.LenSubEnc, s)
	if err != nil {
		return nil, fmt.Errorf("cram: XRLE length sub-codec: %w", err)
	}
	// The length stream starts with a uint7 total output size, then one
	// uint7 run length per RLE-symbol literal.
	outSize, n, err := uint7At64(lenBytes, 0)
	if err != nil {
		return nil, fmt.Errorf("cram: XRLE output size: %w", err)
	}
	if outSize < 0 {
		return nil, fmt.Errorf("cram: XRLE declares negative output size %d", outSize)
	}
	runPos := n
	isRLE := make([]bool, 256)
	for _, sym := range enc.RLESyms {
		isRLE[sym] = true
	}
	out := make([]byte, 0, outSize)
	for _, b := range lit {
		if int64(len(out)) >= outSize {
			return nil, fmt.Errorf("cram: XRLE expanded past its declared output size %d", outSize)
		}
		if isRLE[b] {
			rlen, rn, rerr := uint7At32(lenBytes, runPos)
			if rerr != nil {
				return nil, fmt.Errorf("cram: XRLE run length: %w", rerr)
			}
			runPos += rn
			if rlen < 0 {
				return nil, fmt.Errorf("cram: XRLE negative run length %d", rlen)
			}
			total := int64(rlen) + 1
			if int64(len(out))+total > outSize {
				return nil, fmt.Errorf("cram: XRLE run overruns its declared output size %d", outSize)
			}
			for i := int64(0); i < total; i++ {
				out = append(out, b)
			}
		} else {
			out = append(out, b)
		}
	}
	if int64(len(out)) != outSize {
		return nil, fmt.Errorf("cram: XRLE produced %d bytes, header declares %d", len(out), outSize)
	}
	return out, nil
}

// expandXDelta reverses the XDELTA block path (cram_codecs.c
// cram_xdelta_decode_block): the sub-codec supplies a stream of zig-zag
// uint7 deltas of word_size bytes each, which are prefix-summed (modulo the
// word width) and emitted little-endian. htslib only implements word_size 2
// for the block path; the value-at-a-time int path (word size treated as a
// 32-bit word) is handled in decodeInt, so this covers the byte-producing
// case only.
func (enc *Encoding) expandXDelta(s *seriesSource) ([]byte, error) {
	src, err := subCodecBytes(enc.SubEnc, s)
	if err != nil {
		return nil, fmt.Errorf("cram: XDELTA sub-codec: %w", err)
	}
	w := int(enc.DeltaWordSize)
	if w != 2 {
		return nil, fmt.Errorf("cram: XDELTA block path supports word size 2 only (got %d)", w)
	}
	var out []byte
	var last int64
	pos := 0
	for pos < len(src) {
		v, n, derr := uint7At32(src, pos)
		if derr != nil {
			return nil, fmt.Errorf("cram: XDELTA delta: %w", derr)
		}
		pos += n
		d := int64(int16(unzigzag32(uint32(v))))
		last = int64(int16(last + d))
		out = append(out, byte(last), byte(last>>8))
	}
	return out, nil
}

// unzigzag32 reverses the 32-bit zig-zag map (htscodecs unzigzag32): even
// codes are non-negative, odd codes negative.
func unzigzag32(x uint32) int32 {
	return int32(x>>1) ^ -int32(x&1)
}
