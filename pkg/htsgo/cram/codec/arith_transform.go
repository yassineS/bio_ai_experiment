package codec

import (
	"bytes"
	"compress/bzip2"
	"fmt"
	"io"
)

// x4x16Ext is the X_EXT transform bit: the payload is compressed with an
// external codec (bzip2 in htscodecs) rather than the range coder.
const x4x16Ext = 0x04

// arith_dynamic transform layer — ported from htscodecs arith_dynamic.c
// (arith_compress_to / arith_uncompress_to). arith_dynamic reuses the rANS
// 4x16 format byte and transform bits verbatim:
//
//	X_PACK   0x80  bit-pack 1/2/4/8 symbols per byte (alphabet ≤ 16)
//	X_RLE    0x40  run-length encode (runs coded with the 258-symbol model)
//	X_CAT    0x20  store the payload verbatim (tiny/incompressible inputs)
//	X_NOSZ   0x10  raw size is supplied by the caller, not stored
//	X_STRIPE 0x08  transpose the input into N interleaved streams
//	X_ORDER  0x03  order (0 or 1)
//
// The PACK helpers (htsPack / htsUnpackMeta / htsUnpack) and the varint
// helpers are shared with rANS 4x16 — only the entropy core (the adaptive
// range coder, arith.go) differs. STRIPE here recurses through
// arithUncompressTo, mirroring arith_uncompress_to.

// ArithDecode decompresses a complete arith_dynamic stream (format byte
// included) and returns the raw bytes. Order-0, order-1, the X_CAT
// store-uncompressed form and the PACK/RLE/STRIPE transforms are all
// supported. It is the CRAM compression-method-6 entry point.
func ArithDecode(in []byte) ([]byte, error) {
	return arithUncompressTo(in, 0, 0)
}

// ArithEncode compresses in with the given order byte and returns a
// complete arith_dynamic stream. order carries the X_* transform bits
// OR'd with the entropy order (0 or 1); it mirrors the order argument of
// htscodecs' arith_compress. The output is byte-identical to the
// reference encoder.
func ArithEncode(in []byte, order int) ([]byte, error) {
	return arithCompressTo(in, order)
}

// --- decode ------------------------------------------------------------------

// arithUncompressTo ports arith_uncompress_to. expectSize supplies the raw
// size for X_NOSZ streams (STRIPE sub-streams); top-level callers pass 0.
// depth bounds STRIPE recursion.
func arithUncompressTo(in []byte, expectSize uint32, depth int) ([]byte, error) {
	if len(in) == 0 {
		return nil, fmt.Errorf("arith: empty input")
	}
	if in[0]&x4x16Stripe != 0 {
		return arithUncompressStripe(in, depth)
	}

	format := in[0]
	doPack := format&x4x16Pack != 0
	doRLE := format&x4x16RLE != 0
	doCat := format&x4x16Cat != 0
	doExt := format&x4x16Ext != 0
	noSize := format&x4x16NoSz != 0
	order := int(format & 3)
	cp := 1

	var osz uint32
	if !noSize {
		v, c, ok := varGetU32(in, cp)
		if !ok {
			return nil, fmt.Errorf("arith: truncated raw-size varint")
		}
		osz, cp = v, c
	} else {
		osz = expectSize
	}
	if osz > maxRANSRawSize {
		return nil, fmt.Errorf("arith: declared raw size %d exceeds the %d-byte safety ceiling",
			osz, maxRANSRawSize)
	}

	// PACK meta: the symbol map plus a varint of the packed length.
	var packMap []byte
	var packNSym int
	unpackedSz := osz
	if doPack {
		m, ns, used, ok := htsUnpackMeta(in[cp:])
		if !ok {
			return nil, fmt.Errorf("arith: malformed PACK meta block")
		}
		packMap = m
		packNSym = ns
		cp += used
		v, c, ok := varGetU32(in, cp)
		if !ok {
			return nil, fmt.Errorf("arith: truncated PACK length varint")
		}
		cp = c
		if v > maxRANSRawSize {
			return nil, fmt.Errorf("arith: PACK length %d exceeds the %d-byte safety ceiling",
				v, maxRANSRawSize)
		}
		osz = v
	}

	// Entropy core: in[cp:] -> data. osz is the post-RLE / packed length.
	var data []byte
	if len(in)-cp > 0 {
		payload := in[cp:]
		switch {
		case doCat:
			if int(osz) > len(payload) {
				return nil, fmt.Errorf("arith: X_CAT payload %d bytes, stream holds %d",
					osz, len(payload))
			}
			data = make([]byte, osz)
			copy(data, payload[:osz])
		case doExt:
			// X_EXT: the payload is bzip2 (htscodecs uses libbz2).
			d, err := io.ReadAll(bzip2.NewReader(bytes.NewReader(payload)))
			if err != nil {
				return nil, fmt.Errorf("arith: X_EXT bzip2 decode: %w", err)
			}
			if uint32(len(d)) != osz {
				return nil, fmt.Errorf("arith: X_EXT decoded %d bytes, expected %d", len(d), osz)
			}
			data = d
		case doRLE && order == 1:
			d, err := arithUncompressO1RLE(payload, osz)
			if err != nil {
				return nil, err
			}
			data = d
		case doRLE:
			d, err := arithUncompressO0RLE(payload, osz)
			if err != nil {
				return nil, err
			}
			data = d
		case order == 1:
			d, err := arithUncompressO1(payload, osz)
			if err != nil {
				return nil, err
			}
			data = d
		default:
			d, err := arithUncompressO0(payload, osz)
			if err != nil {
				return nil, err
			}
			data = d
		}
	} else {
		data = nil
	}

	// un-PACK: data -> unpacked.
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

// arithUncompressStripe ports the STRIPE branch of arith_uncompress_to: it
// reads per-stripe compressed lengths, recursively decompresses each
// stripe with X_NOSZ and interleaves the results.
func arithUncompressStripe(in []byte, depth int) ([]byte, error) {
	if depth > maxStripeDepth {
		return nil, fmt.Errorf("arith: STRIPE nested deeper than %d levels", maxStripeDepth)
	}
	cMeta := 1
	ulen, cp, ok := varGetU32(in, cMeta)
	if !ok {
		return nil, fmt.Errorf("arith: truncated STRIPE raw-size varint")
	}
	cMeta = cp
	if cMeta >= len(in) {
		return nil, fmt.Errorf("arith: STRIPE stream truncated before stripe count")
	}
	if ulen > maxRANSRawSize {
		return nil, fmt.Errorf("arith: STRIPE raw size %d exceeds the %d-byte safety ceiling",
			ulen, maxRANSRawSize)
	}
	n := int(in[cMeta])
	cMeta++
	if n < 1 {
		return nil, fmt.Errorf("arith: STRIPE needs at least one stream")
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
			return nil, fmt.Errorf("arith: truncated STRIPE stream length %d", i)
		}
		cMeta = c
		clenN[i] = int(clen)
		clenTot += int(clen)
		if cMeta > len(in) || clenN[i] > len(in) || clenN[i] < 1 {
			return nil, fmt.Errorf("arith: STRIPE stream %d length %d invalid", i, clenN[i])
		}
	}
	if cMeta+clenTot > len(in) {
		return nil, fmt.Errorf("arith: STRIPE streams overrun the payload")
	}

	outN := make([]byte, ulen)
	for i := 0; i < n; i++ {
		stripe, err := arithUncompressTo(in[cMeta:cMeta+clenN[i]], uint32(ulenN[i]), depth+1)
		if err != nil {
			return nil, fmt.Errorf("arith: STRIPE stream %d: %w", i, err)
		}
		if len(stripe) != ulenN[i] {
			return nil, fmt.Errorf("arith: STRIPE stream %d decoded %d bytes, want %d",
				i, len(stripe), ulenN[i])
		}
		copy(outN[idxN[i]:], stripe)
		cMeta += clenN[i]
	}

	// unstripe: interleave the N streams back into the original order.
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

// --- encode ------------------------------------------------------------------

// arithCompressTo ports arith_compress_to: it handles STRIPE (which
// recurses), then PACK and RLE around the adaptive range-coder core. order
// carries the X_* transform bits OR'd with the entropy order (0 or 1).
func arithCompressTo(in []byte, order int) ([]byte, error) {
	if len(in) <= 20 {
		order &^= x4x16Stripe
	}

	if order&x4x16Cat != 0 {
		out := []byte{x4x16Cat}
		out = varPutU32(out, uint32(len(in)))
		out = append(out, in...)
		return out, nil
	}

	if order&x4x16Ext != 0 {
		// X_EXT compresses with libbz2 in htscodecs. Go's standard
		// library ships a bzip2 *decoder* only, and adding a third-party
		// bzip2 encoder is outside the sanctioned dependency set
		// (CLAUDE.md). Decode of X_EXT streams is fully supported; encode
		// is the documented gap.
		return nil, fmt.Errorf("arith: X_EXT (bzip2) encode is unsupported — " +
			"Go has no standard-library bzip2 encoder; decode is supported")
	}

	if order&x4x16Stripe != 0 {
		return arithCompressStripe(in, order)
	}

	doPack := order&x4x16Pack != 0
	doRLE := order&x4x16RLE != 0
	noSize := order&x4x16NoSz != 0
	formatByte := byte(order)

	out := []byte{formatByte}
	if !noSize {
		out = varPutU32(out, uint32(len(in)))
	}

	entOrder := order & 3
	data := in

	// PACK: bit-pack the input and emit the pack meta block plus a varint
	// of the packed length. hts_pack declines when the alphabet > 16.
	if doPack && len(data) != 0 {
		packed, meta, ok := htsPack(data)
		if !ok {
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

	if doRLE && len(data) == 0 {
		formatByte &^= x4x16RLE
		out[0] = formatByte
		doRLE = false
	}

	// htscodecs drops order 1 to order 0 when the input is below 8 bytes.
	if entOrder != 0 && len(data) < 8 {
		entOrder = 0
		formatByte &^= 3
		out[0] = formatByte
	}

	var payload []byte
	switch {
	case doRLE && entOrder == 1:
		payload = arithCompressO1RLE(data)
	case doRLE:
		payload = arithCompressO0RLE(data)
	case entOrder == 1:
		payload = arithCompressO1(data)
	default:
		payload = arithCompressO0(data)
	}

	// X_CAT fallback when the range coder would not shrink the (possibly
	// already transformed) data. htscodecs keeps PACK but clears the
	// order/RLE bits and sets X_CAT.
	if len(payload) >= len(data) {
		formatByte = (formatByte &^ 3) &^ x4x16RLE
		formatByte |= x4x16Cat
		out[0] = formatByte
		out = append(out, data...)
		return out, nil
	}

	out = append(out, payload...)
	return out, nil
}

// arithCompressStripe ports the STRIPE branch of arith_compress_to: it
// transposes the input into N interleaved streams (N defaults to 4) and
// recursively compresses each, trying the per-stripe methods from the C
// reference's m[][] table and keeping the smallest.
func arithCompressStripe(in []byte, order int) ([]byte, error) {
	n := (order >> 8) & 0xff
	if n == 0 {
		n = 4
	}
	if n > len(in) {
		n = len(in)
	}
	if n < 1 {
		// Unreachable in practice: arithCompressTo clears X_STRIPE for
		// inputs of 20 bytes or fewer, so len(in) > 20 and n >= 1 here.
		// Kept as defensive cover against a future caller bypassing that.
		return nil, fmt.Errorf("arith: STRIPE needs a non-empty input")
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
	for i := x * n; i < inSize; i, x = i+n, x+1 {
		for j := 0; i+j < inSize; j++ {
			transposed[idx[j]+x] = in[i+j]
		}
	}

	out := []byte{byte(order &^ x4x16NoSz)}
	out = varPutU32(out, uint32(inSize))
	out = append(out, byte(n))

	// Per-stripe method table from arith_compress_to: stripe 0 tries
	// {1,64,0}, stripes 1+ try {1,0} / {1,128}. Methods with bit 1 set
	// (order-1) are skipped when the requested order is 0.
	mTable := [][]int{
		{1, 64, 0},
		{1, 0},
		{1, 128},
		{1, 128},
	}
	var streams []byte
	for i := 0; i < n; i++ {
		// arith_compress_to caps the stripe index at 3 (MIN(i,3)) when
		// indexing the per-stripe method table.
		mi := i
		if mi > 3 {
			mi = 3
		}
		methods := mTable[mi]
		var best []byte
		for _, m := range methods {
			if order&3 == 0 && m&1 != 0 {
				continue
			}
			enc, err := arithCompressTo(transposed[idx[i]:idx[i]+partLen[i]], m|x4x16NoSz)
			if err != nil {
				return nil, err
			}
			if best == nil || len(enc) < len(best) {
				best = enc
			}
		}
		if best == nil {
			return nil, fmt.Errorf("arith: STRIPE stream %d has no viable method", i)
		}
		out = varPutU32(out, uint32(len(best)))
		streams = append(streams, best...)
	}
	out = append(out, streams...)
	return out, nil
}
