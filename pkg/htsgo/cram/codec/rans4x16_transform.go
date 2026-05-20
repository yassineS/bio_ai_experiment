package codec

import "fmt"

// rANS 4x16 transform layer — the C2.2 slice. Ported from htscodecs
// rANS_static4x16pr.c (rans_compress_to_4x16 / rans_uncompress_to_4x16),
// pack.c (hts_pack / hts_unpack_meta / hts_unpack), rle.c (hts_rle_encode /
// hts_rle_decode / rle_find_syms) and utils.h (unstripe). The on-wire
// format is byte-identical to htscodecs — see rans4x16_test.go for the
// compliance vectors.
//
// The transforms wrap the order-0/order-1 rANS core. Their format-byte
// bits are:
//
//	X_PACK   0x80  bit-pack 1/2/4/8 symbols per byte (alphabet ≤ 16)
//	X_RLE    0x40  run-length encode, runs and literals in separate streams
//	X_STRIPE 0x08  transpose the input into N interleaved streams
//	X_NOSZ   0x10  the raw size is supplied by the caller, not stored
//
// On encode the pipeline is PACK then RLE then rANS; on decode it is
// reversed (rANS, un-RLE, un-PACK). STRIPE is handled separately, before
// the format byte is even consulted: it splits the input, recursively
// compresses each stripe with X_NOSZ and interleaves the results.
//
// X_32 (0x04, 32-way SIMD unrolling) is a distinct on-wire format and is
// still rejected; see RANS4x16Decode.

// --- public dispatch helpers -------------------------------------------------

// transformOrderRANS4x16 reports whether order requests any of the
// implemented transforms (PACK/RLE/STRIPE) and, if so, returns the
// rANS-4x16 stream produced by the transform pipeline. When ok is false
// the caller falls back to the plain order-0/order-1 path.
func transformOrderRANS4x16(in []byte, order int) (out []byte, ok bool, err error) {
	if order&(x4x16Pack|x4x16RLE|x4x16Stripe) == 0 {
		return nil, false, nil
	}
	out, err = compressToRANS4x16(in, order)
	return out, true, err
}

// --- STRIPE ------------------------------------------------------------------

// compressToRANS4x16 ports rans_compress_to_4x16's transform-aware path:
// it handles STRIPE (which recurses), then PACK and RLE around the rANS
// core. order carries the X_* transform bits OR'd with the rANS order
// (0 or 1). When NOSZ is set the raw-size varint is omitted.
func compressToRANS4x16(in []byte, order int) ([]byte, error) {
	// STRIPE on inputs of 20 bytes or fewer is disabled by htscodecs:
	// the per-stripe rANS overhead dwarfs any gain.
	if order&x4x16Stripe != 0 && len(in) <= 20 {
		order &^= x4x16Stripe
	}

	if order&x4x16Stripe != 0 {
		return compressStripeRANS4x16(in, order)
	}

	noSize := order&x4x16NoSz != 0
	doPack := order&x4x16Pack != 0
	doRLE := order&x4x16RLE != 0
	formatByte := byte(order)

	out := []byte{formatByte}
	if !noSize {
		out = varPutU32(out, uint32(len(in)))
	}

	ransOrder := order & 1
	data := in

	// PACK: bit-pack the input and emit the pack meta block plus a
	// varint of the packed length. hts_pack fails (returns ok=false)
	// when the alphabet exceeds 16 symbols.
	if doPack && len(data) != 0 {
		packed, meta, packedOK := htsPack(data)
		if !packedOK {
			formatByte &^= x4x16Pack
			out[0] = formatByte
			doPack = false
		} else {
			out = append(out, meta...)
			out = varPutU32(out, uint32(len(packed)))
			data = packed
		}
	} else if doPack {
		formatByte &^= x4x16Pack
		out[0] = formatByte
		doPack = false
	}

	// RLE: split runs from literals. The run-length meta block is
	// itself O0-rANS-compressed, or stored raw when that is smaller.
	if doRLE && len(data) != 0 {
		lits, run, syms := htsRLEEncode(data)
		meta := make([]byte, 0, len(run)+len(syms)+1)
		meta = append(meta, byte(len(syms)))
		meta = append(meta, syms...)
		meta = append(meta, run...)
		rmetaLen := len(meta)

		// htscodecs only keeps RLE when it shrinks the data enough to
		// be worth the speed hit (rle_len + rmeta_len < 0.99*in_size).
		if len(lits)+rmetaLen >= int(0.99*float64(len(data))) {
			formatByte &^= x4x16RLE
			out[0] = formatByte
			doRLE = false
		} else {
			cMeta := compressO0RANS4x16(meta)
			if cMeta != nil && len(cMeta) < rmetaLen {
				out = varPutU32(out, uint32(rmetaLen*2))
				out = varPutU32(out, uint32(len(lits)))
				out = varPutU32(out, uint32(len(cMeta)))
				out = append(out, cMeta...)
			} else {
				out = varPutU32(out, uint32(rmetaLen*2+1))
				out = varPutU32(out, uint32(len(lits)))
				out = append(out, meta...)
			}
			data = lits
		}
	} else if doRLE {
		formatByte &^= x4x16RLE
		out[0] = formatByte
		doRLE = false
	}

	// htscodecs drops order 1 to order 0 below 8 bytes of rANS input.
	if ransOrder != 0 && len(data) < 8 {
		ransOrder = 0
		formatByte &^= 1
		out[0] = formatByte
	}

	var payload []byte
	if ransOrder != 0 {
		payload = compressO1RANS4x16(data, 0)
	} else {
		payload = compressO0RANS4x16(data)
	}

	// X_CAT fallback when rANS would not shrink the (possibly already
	// transformed) data.
	if len(payload) >= len(data) {
		formatByte = (formatByte &^ 3) | x4x16Cat
		if noSize {
			formatByte |= x4x16NoSz
		}
		out[0] = formatByte
		out = append(out, data...)
		return out, nil
	}

	out = append(out, payload...)
	return out, nil
}

// compressStripeRANS4x16 ports the STRIPE branch of
// rans_compress_to_4x16: it transposes the input into N interleaved
// streams (N defaults to 4) and recursively compresses each, trying the
// methods {order-1, RLE, PACK, order-0} and keeping the smallest.
func compressStripeRANS4x16(in []byte, order int) ([]byte, error) {
	n := (order >> 8) & 0xff
	if n == 0 {
		n = 4
	}
	if n > len(in) {
		n = len(in)
	}
	if n < 1 {
		// An empty input cannot be striped; fall back to plain order-0.
		return frameRANS4x16(in, compressO0RANS4x16(in), 0x00), nil
	}

	inSize := len(in)
	partLen := make([]int, n)
	idx := make([]int, n)
	for i := 0; i < n; i++ {
		partLen[i] = inSize/n + boolToInt(inSize%n > i)
		if i > 0 {
			idx[i] = idx[i-1] + partLen[i-1]
		}
	}

	transposed := make([]byte, inSize)
	x := 0
	for i := 0; i+n <= inSize-1; i, x = i+n, x+1 {
		for j := 0; j < n; j++ {
			transposed[idx[j]+x] = in[i+j]
		}
	}
	// Trailing partial group: htscodecs runs i from the loop above to
	// in_size in steps of N, copying j while i+j < in_size.
	for i := x * n; i < inSize; i, x = i+n, x+1 {
		for j := 0; i+j < inSize; j++ {
			transposed[idx[j]+x] = in[i+j]
		}
	}

	out := []byte{byte(order &^ x4x16NoSz)}
	out = varPutU32(out, uint32(inSize))
	out = append(out, byte(n))

	// methods: {order-1, RLE, order-0+PACK, order-0}. Each substream is
	// compressed with X_NOSZ; the smallest wins.
	methods := []int{1, x4x16RLE, x4x16Pack, 0}
	var streams []byte
	for i := 0; i < n; i++ {
		var best []byte
		for _, m := range methods {
			if order&m != m {
				continue
			}
			enc, err := compressToRANS4x16(transposed[idx[i]:idx[i]+partLen[i]], m|x4x16NoSz)
			if err != nil {
				return nil, err
			}
			if best == nil || len(enc) < len(best) {
				best = enc
			}
		}
		if best == nil {
			return nil, fmt.Errorf("rans4x16: STRIPE stream %d has no viable method", i)
		}
		out = varPutU32(out, uint32(len(best)))
		streams = append(streams, best...)
	}
	out = append(out, streams...)
	return out, nil
}

// boolToInt returns 1 for true and 0 for false.
func boolToInt(b bool) int {
	if b {
		return 1
	}
	return 0
}

// uncompressStripeRANS4x16 ports the STRIPE branch of
// rans_uncompress_to_4x16: it reads the per-stripe compressed lengths,
// recursively decompresses each stripe with X_NOSZ and interleaves the
// results. The raw size is read from the stream (STRIPE is always
// top-level, never itself NOSZ).
func uncompressStripeRANS4x16(in []byte) ([]byte, error) {
	cMeta := 1
	ulen, cp, ok := varGetU32(in, cMeta)
	if !ok {
		return nil, fmt.Errorf("rans4x16: truncated STRIPE raw-size varint")
	}
	cMeta = cp
	if cMeta >= len(in) {
		return nil, fmt.Errorf("rans4x16: STRIPE stream truncated before stripe count")
	}
	if ulen > maxRANSRawSize {
		return nil, fmt.Errorf("rans4x16: STRIPE raw size %d exceeds the %d-byte safety ceiling",
			ulen, maxRANSRawSize)
	}
	n := int(in[cMeta])
	cMeta++
	if n < 1 {
		return nil, fmt.Errorf("rans4x16: STRIPE needs at least one stream")
	}

	ulenN := make([]int, n)
	idxN := make([]int, n)
	clenN := make([]int, n)
	clenTot := 0
	for i := 0; i < n; i++ {
		ulenN[i] = int(ulen)/n + boolToInt(int(ulen)%n > i)
		if i > 0 {
			idxN[i] = idxN[i-1] + ulenN[i-1]
		}
		clen, c, ok := varGetU32(in, cMeta)
		if !ok {
			return nil, fmt.Errorf("rans4x16: truncated STRIPE stream length %d", i)
		}
		cMeta = c
		clenN[i] = int(clen)
		clenTot += int(clen)
		if cMeta > len(in) || clenN[i] > len(in) || clenN[i] < 1 {
			return nil, fmt.Errorf("rans4x16: STRIPE stream %d length %d invalid", i, clenN[i])
		}
	}
	if cMeta+clenTot > len(in) {
		return nil, fmt.Errorf("rans4x16: STRIPE streams overrun the payload")
	}

	outN := make([]byte, ulen)
	for i := 0; i < n; i++ {
		stripe, err := decodeRANS4x16WithSize(in[cMeta:cMeta+clenN[i]], uint32(ulenN[i]))
		if err != nil {
			return nil, fmt.Errorf("rans4x16: STRIPE stream %d: %w", i, err)
		}
		if len(stripe) != ulenN[i] {
			return nil, fmt.Errorf("rans4x16: STRIPE stream %d decoded %d bytes, want %d",
				i, len(stripe), ulenN[i])
		}
		copy(outN[idxN[i]:], stripe)
		cMeta += clenN[i]
	}

	// unstripe: interleave the N streams back into the original order
	// (utils.h). The tuned C cases are just an unrolled round-robin.
	out := make([]byte, ulen)
	idx := make([]int, n)
	copy(idx, idxN)
	for j := 0; j < int(ulen); j++ {
		k := j % n
		out[j] = outN[idx[k]]
		idx[k]++
	}
	return out, nil
}

// --- decode pipeline ---------------------------------------------------------

// decodeRANS4x16WithSize decodes a rANS 4x16 stream whose raw size is
// supplied by the caller (the X_NOSZ case used by STRIPE sub-streams)
// rather than read from a varint. A stream that does store its own size
// ignores expectSize and uses the stored value.
func decodeRANS4x16WithSize(in []byte, expectSize uint32) ([]byte, error) {
	if len(in) == 0 {
		return nil, fmt.Errorf("rans4x16: empty input")
	}
	format := in[0]
	if format&x4x16X32 != 0 {
		return nil, fmt.Errorf("rans4x16: format byte 0x%02x uses the X_32 transform, "+
			"which is not implemented", format)
	}
	if format&x4x16Stripe != 0 {
		return uncompressStripeRANS4x16(in)
	}
	return uncompressTransformRANS4x16(in, expectSize)
}

// uncompressTransformRANS4x16 ports the non-STRIPE body of
// rans_uncompress_to_4x16: it parses the format byte, optional PACK and
// RLE meta blocks, runs the rANS core, then reverses the transforms
// (un-RLE, un-PACK). noSize streams take their raw size from expectSize.
func uncompressTransformRANS4x16(in []byte, expectSize uint32) ([]byte, error) {
	format := in[0]
	doPack := format&x4x16Pack != 0
	doRLE := format&x4x16RLE != 0
	doCat := format&x4x16Cat != 0
	noSize := format&x4x16NoSz != 0
	ransOrder := int(format & 1)
	cp := 1

	var osz uint32
	if !noSize {
		v, c, ok := varGetU32(in, cp)
		if !ok {
			return nil, fmt.Errorf("rans4x16: truncated raw-size varint")
		}
		osz, cp = v, c
	} else {
		osz = expectSize
	}
	if osz > maxRANSRawSize {
		return nil, fmt.Errorf("rans4x16: declared raw size %d exceeds the %d-byte safety ceiling",
			osz, maxRANSRawSize)
	}

	// PACK meta: the symbol map plus a varint of the packed length.
	var packMap []byte
	var packNSym int
	unpackedSz := osz
	if doPack {
		m, ns, used, ok := htsUnpackMeta(in[cp:])
		if !ok {
			return nil, fmt.Errorf("rans4x16: malformed PACK meta block")
		}
		packMap = m
		packNSym = ns
		cp += used
		v, c, ok := varGetU32(in, cp)
		if !ok {
			return nil, fmt.Errorf("rans4x16: truncated PACK length varint")
		}
		cp = c
		// v is the rANS input length when packing applies; it is also
		// the raw size for the un-RLE stage below.
		osz = v
	}

	// RLE meta: the run-length meta block, raw or O0-rANS-compressed.
	var rleMeta []byte
	var rleLitLen uint32
	if doRLE {
		uMetaSize, c, ok := varGetU32(in, cp)
		if !ok {
			return nil, fmt.Errorf("rans4x16: truncated RLE meta-size varint")
		}
		rleLen, c2, ok := varGetU32(in, c)
		if !ok {
			return nil, fmt.Errorf("rans4x16: truncated RLE literal-length varint")
		}
		rleLitLen = rleLen
		if uMetaSize&1 != 0 {
			// Raw (uncompressed) RLE meta block.
			metaLen := uMetaSize / 2
			if int(metaLen) > len(in)-c2 {
				metaLen = uint32(len(in) - c2)
			}
			rleMeta = in[c2 : c2+int(metaLen)]
			cp = c2 + int(metaLen)
		} else {
			cMetaSize, c3, ok := varGetU32(in, c2)
			if !ok {
				return nil, fmt.Errorf("rans4x16: truncated RLE compressed-meta-size varint")
			}
			uMeta := uMetaSize / 2
			if c3+int(cMetaSize) > len(in) {
				return nil, fmt.Errorf("rans4x16: RLE meta block overruns the payload")
			}
			decoded, err := uncompressO0RANS4x16(in[c3:c3+int(cMetaSize)], uMeta)
			if err != nil {
				return nil, fmt.Errorf("rans4x16: RLE meta decompress: %w", err)
			}
			rleMeta = decoded
			cp = c3 + int(cMetaSize)
		}
	}

	// rANS core: in[cp:] -> tmp1. When RLE applies, the rANS output is
	// the literal stream of length rleLitLen.
	tmp1Size := osz
	if doRLE {
		tmp1Size = rleLitLen
	}
	var tmp1 []byte
	if len(in)-cp > 0 {
		if doCat {
			if int(tmp1Size) > len(in)-cp {
				return nil, fmt.Errorf("rans4x16: X_CAT payload %d bytes, stream holds %d",
					tmp1Size, len(in)-cp)
			}
			tmp1 = make([]byte, tmp1Size)
			copy(tmp1, in[cp:cp+int(tmp1Size)])
		} else {
			var err error
			if ransOrder != 0 {
				tmp1, err = uncompressO1RANS4x16(in[cp:], tmp1Size)
			} else {
				tmp1, err = uncompressO0RANS4x16(in[cp:], tmp1Size)
			}
			if err != nil {
				return nil, err
			}
		}
	} else {
		tmp1 = nil
	}

	// un-RLE: tmp1 -> tmp2.
	data := tmp1
	if doRLE {
		out, err := htsRLEDecode(tmp1, rleMeta)
		if err != nil {
			return nil, err
		}
		data = out
	}

	// un-PACK: tmp2 -> tmp3.
	if doPack {
		if packNSym == 1 {
			unpackedSz = uint32(len(data))
		}
		out, err := htsUnpack(data, packNSym, packMap, unpackedSz)
		if err != nil {
			return nil, err
		}
		data = out
	}
	return data, nil
}

// --- PACK --------------------------------------------------------------------

// htsPack ports hts_pack: it bit-packs data into 1/2/4/8 symbols per
// byte when the alphabet is 16 symbols or fewer. It returns the packed
// bytes, the meta block (symbol count followed by the symbol map) and
// ok=false when the alphabet is too large to pack.
func htsPack(data []byte) (packed, meta []byte, ok bool) {
	var present [256]bool
	for _, b := range data {
		present[b] = true
	}
	var code [256]int
	n := 0
	meta = make([]byte, 1, 17)
	for i := 0; i < 256; i++ {
		if present[i] {
			code[i] = n
			n++
			meta = append(meta, byte(i))
		}
	}
	meta[0] = byte(n) // 256 wraps to 0, matching hts_pack.
	if n > 16 {
		return nil, nil, false
	}

	var valPerByte int
	switch {
	case n > 4:
		valPerByte = 2
	case n > 2:
		valPerByte = 4
	case n > 1:
		valPerByte = 8
	default:
		valPerByte = 0 // constant input: 0 packed bytes.
	}

	length := len(data)
	out := make([]byte, 0, length+1)
	switch valPerByte {
	case 2:
		i := 0
		for ; i < length&^1; i += 2 {
			out = append(out, byte(code[data[i]])|byte(code[data[i+1]])<<4)
		}
		if length-i == 1 {
			out = append(out, byte(code[data[i]]))
		}
	case 4:
		i := 0
		for ; i < length&^3; i += 4 {
			out = append(out, byte(code[data[i]])|byte(code[data[i+1]])<<2|
				byte(code[data[i+2]])<<4|byte(code[data[i+3]])<<6)
		}
		if s := length - i; s > 0 {
			var b byte
			for x := 0; x < s; x++ {
				b |= byte(code[data[i+x]]) << (2 * x)
			}
			out = append(out, b)
		}
	case 8:
		i := 0
		for ; i < length&^7; i += 8 {
			var b byte
			for x := 0; x < 8; x++ {
				b |= byte(code[data[i+x]]) << x
			}
			out = append(out, b)
		}
		if s := length - i; s > 0 {
			var b byte
			for x := 0; x < s; x++ {
				b |= byte(code[data[i+x]]) << x
			}
			out = append(out, b)
		}
	case 0:
		// Constant input packs to zero bytes.
	}
	return out, meta, true
}

// htsUnpackMeta ports hts_unpack_meta: it reads the PACK meta block (a
// symbol count followed by the symbol map). It returns the 16-entry
// symbol map (zero-padded past the declared symbols, matching the C
// caller's uint8_t map[16]={0} so an out-of-range packed nibble decodes
// to a defined symbol), the number of symbols per byte (0/1/2/4/8), the
// bytes consumed and ok=false on a malformed block.
func htsUnpackMeta(data []byte) (m []byte, nsym, used int, ok bool) {
	if len(data) == 0 {
		return nil, 0, 0, false
	}
	n := int(data[0])
	if n == 0 {
		n = 256
	}
	switch {
	case n <= 1:
		nsym = 0
	case n <= 2:
		nsym = 8
	case n <= 4:
		nsym = 4
	case n <= 16:
		nsym = 2
	default:
		// More than 16 symbols: no packing, one byte of meta consumed.
		return nil, 1, 1, true
	}
	if len(data) <= 1 {
		return nil, 0, 0, false
	}
	m = make([]byte, 16)
	j := 1
	for c := 0; c < n && j < len(data); c, j = c+1, j+1 {
		m[c] = data[j]
	}
	if j-1 < n {
		return nil, 0, 0, false
	}
	return m, nsym, j, true
}

// htsUnpack ports hts_unpack: it expands a bit-packed stream back to the
// original bytes using the symbol map. unpackedSz is the expected
// output length.
func htsUnpack(data []byte, nsym int, m []byte, unpackedSz uint32) ([]byte, error) {
	if nsym == 1 {
		out := make([]byte, len(data))
		copy(out, data)
		return out, nil
	}
	out := make([]byte, unpackedSz)
	outLen := int(unpackedSz)
	switch nsym {
	case 8:
		if (outLen+7)/8 > len(data) {
			return nil, fmt.Errorf("rans4x16: PACK(8) stream too short")
		}
		olen := outLen &^ 7
		j := 0
		for i := 0; i < olen; i += 8 {
			c := data[j]
			j++
			for x := 0; x < 8; x++ {
				out[i+x] = m[c>>x&1]
			}
		}
		if outLen != olen {
			c := data[j]
			for i := olen; i < outLen; i++ {
				out[i] = m[c&1]
				c >>= 1
			}
		}
	case 4:
		if (outLen+3)/4 > len(data) {
			return nil, fmt.Errorf("rans4x16: PACK(4) stream too short")
		}
		olen := outLen &^ 3
		j := 0
		for i := 0; i < olen; i += 4 {
			c := data[j]
			j++
			out[i] = m[c&3]
			out[i+1] = m[c>>2&3]
			out[i+2] = m[c>>4&3]
			out[i+3] = m[c>>6&3]
		}
		if outLen != olen {
			c := data[j]
			for i := olen; i < outLen; i++ {
				out[i] = m[c&3]
				c >>= 2
			}
		}
	case 2:
		if (outLen+1)/2 > len(data) {
			return nil, fmt.Errorf("rans4x16: PACK(2) stream too short")
		}
		olen := outLen &^ 1
		j := 0
		for i := 0; i < olen; i += 2 {
			c := data[j]
			j++
			out[i] = m[c&15]
			out[i+1] = m[c>>4&15]
		}
		if outLen != olen {
			c := data[j]
			out[olen] = m[c&15]
		}
	case 0:
		// Constant input: every byte is the single mapped symbol.
		for i := range out {
			out[i] = m[0]
		}
	default:
		return nil, fmt.Errorf("rans4x16: invalid PACK symbols-per-byte %d", nsym)
	}
	return out, nil
}

// --- RLE ---------------------------------------------------------------------

// rleFindSyms ports rle_find_syms: it picks the symbols worth run-length
// encoding. A symbol scores +1 for each repeat and -1 for each non-run
// occurrence; a positive score means RLE pays off.
func rleFindSyms(data []byte) []byte {
	var saved [256]int64
	last := -1
	for _, b := range data {
		if int(b) == last {
			saved[b]++
		} else {
			saved[b]--
			last = int(b)
		}
	}
	var syms []byte
	for i := 0; i < 256; i++ {
		if saved[i] > 0 {
			syms = append(syms, byte(i))
		}
	}
	return syms
}

// htsRLEEncode ports hts_rle_encode: it splits data into a literal
// stream and a run-length stream. Every byte is copied to the literals;
// a byte chosen for RLE is additionally followed in the run stream by a
// varint of (run length − 1). The chosen symbols are returned too.
func htsRLEEncode(data []byte) (lits, run, syms []byte) {
	syms = rleFindSyms(data)
	var isRLE [256]bool
	for _, s := range syms {
		isRLE[s] = true
	}
	for i := 0; i < len(data); i++ {
		lits = append(lits, data[i])
		if isRLE[data[i]] {
			start := i
			last := data[i]
			for i < len(data) && data[i] == last {
				i++
			}
			i--
			run = varPutU32(run, uint32(i-start))
		}
	}
	return lits, run, syms
}

// htsRLEDecode ports hts_rle_decode: it expands the literal stream using
// the run-length stream. meta is [nsyms][syms...][runs...] — the same
// layout the encoder builds. A run-length symbol consumes one run-length
// varint and is emitted (run+1) times.
func htsRLEDecode(lit, meta []byte) ([]byte, error) {
	if len(meta) == 0 {
		return nil, fmt.Errorf("rans4x16: empty RLE meta block")
	}
	nsyms := int(meta[0])
	if nsyms == 0 {
		nsyms = 256
	}
	if len(meta) < 1+nsyms {
		return nil, fmt.Errorf("rans4x16: RLE meta block too short for %d symbols", nsyms)
	}
	var isRLE [256]bool
	for _, s := range meta[1 : 1+nsyms] {
		isRLE[s] = true
	}
	run := meta[1+nsyms:]
	rp := 0

	out := make([]byte, 0, len(lit)*2)
	for _, b := range lit {
		if isRLE[b] {
			rlen, c, ok := varGetU32(run, rp)
			if !ok {
				return nil, fmt.Errorf("rans4x16: truncated RLE run-length varint")
			}
			rp = c
			for k := uint32(0); k <= rlen; k++ {
				out = append(out, b)
			}
			if int64(len(out)) > maxRANSRawSize {
				return nil, fmt.Errorf("rans4x16: RLE expansion exceeds the safety ceiling")
			}
		} else {
			out = append(out, b)
		}
	}
	return out, nil
}
