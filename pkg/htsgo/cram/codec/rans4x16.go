package codec

import "fmt"

// rANS 4x16 — the entropy coder used by CRAM v3.1 blocks.
//
// Ported from samtools/htscodecs rANS_static4x16pr.c, rANS_word.h,
// rANS_static16_int.h and varint.h. The on-wire format is byte-identical
// to htscodecs so the htscodecs test corpus
// (reference_code/htscodecs/tests/dat/r4x16/) is the compliance oracle —
// see rans4x16_test.go.
//
// Order 0 and order 1 are both implemented (C2 + C2.1). The format-byte
// transform layer (X_PACK, X_RLE, X_STRIPE, X_32) is deferred to a later
// slice (C2.2) — see docs/CRAM_ROADMAP.md. A stream that asks for any of
// those is rejected with a clear error rather than mis-decoded. The
// order-1 model lives in rans4x16_o1.go.
//
// Stream layout (order-0):
//
//	[formatByte:1][rawSize:varint][freq table][rANS bytes]
//
// formatByte is 0x00 for a plain order-0 stream, 0x01 for order-1, or
// 0x20 (X_CAT) for the store-uncompressed fallback the encoder picks
// when rANS would expand the input. rawSize is a big-endian base-128
// varint. Unlike rANS 4x8
// the renormalisation is 16-bit ("word"): the encoder emits two bytes at
// a time and the decoder refills two bytes at a time.
//
// Frequency table: the alphabet is delta-RLE encoded (a symbol byte,
// optionally followed by a run length when it abuts the previous
// symbol), terminated by a 0 byte; the frequencies then follow as
// big-endian varints, one per present symbol. The stored frequencies
// sum to a power of two ≤ 4096 (round2(rawSize) capped at TOTFREQ); both
// encoder and decoder scale them up to TOTFREQ by a left shift.

const (
	// rans4x16ByteL is the lower bound of the 4x16 normalisation
	// interval. It is 1<<15, not 4x8's 1<<23: the 4x16 coder renorms a
	// 16-bit word at a time, so the interval is 16 bits wide.
	rans4x16ByteL = 1 << 15
)

// transform-layer bits of the format byte (see rANS_static16_int.h).
const (
	x4x16Pack   = 0x80
	x4x16RLE    = 0x40
	x4x16Cat    = 0x20
	x4x16NoSz   = 0x10
	x4x16Stripe = 0x08
	x4x16X32    = 0x04
	// x4x16Unsupported collects every format-byte bit C2 cannot handle.
	// The low bit (order-1) is checked separately.
	x4x16Unsupported = x4x16Pack | x4x16RLE | x4x16NoSz | x4x16Stripe | x4x16X32
)

// RANS4x16Decode decompresses a complete rANS 4x16 stream (format byte
// included) and returns the raw bytes. Order-0, order-1 and the X_CAT
// store-uncompressed form are supported; a stream using a transform
// (PACK/RLE/STRIPE/X32/NOSZ) is rejected with an error.
func RANS4x16Decode(in []byte) ([]byte, error) {
	if len(in) == 0 {
		return nil, fmt.Errorf("rans4x16: empty input")
	}
	format := in[0]
	if format&x4x16Unsupported != 0 {
		return nil, fmt.Errorf("rans4x16: format byte 0x%02x uses an unsupported transform "+
			"(PACK/RLE/STRIPE/X32/NOSZ); only order-0 and order-1 are implemented", format)
	}

	rawSize, cp, ok := varGetU32(in, 1)
	if !ok {
		return nil, fmt.Errorf("rans4x16: truncated raw-size varint")
	}
	if rawSize > maxRANSRawSize {
		return nil, fmt.Errorf("rans4x16: declared raw size %d exceeds the %d-byte safety ceiling",
			rawSize, maxRANSRawSize)
	}

	if format&x4x16Cat != 0 {
		if cp+int(rawSize) > len(in) {
			return nil, fmt.Errorf("rans4x16: X_CAT payload %d bytes, stream holds %d",
				rawSize, len(in)-cp)
		}
		out := make([]byte, rawSize)
		copy(out, in[cp:cp+int(rawSize)])
		return out, nil
	}

	if format&1 != 0 {
		return uncompressO1RANS4x16(in[cp:], rawSize)
	}
	return uncompressO0RANS4x16(in[cp:], rawSize)
}

// RANS4x16Encode compresses in with the rANS 4x16 codec (order 0 or 1)
// and returns a complete stream. The output is byte-identical to
// htscodecs' rans_compress_to_4x16, including its X_CAT fallback when
// rANS would not shrink the input and its downgrade of order 1 to
// order 0 for inputs below 8 bytes.
func RANS4x16Encode(in []byte, order int) ([]byte, error) {
	switch order {
	case 0:
		return frameRANS4x16(in, compressO0RANS4x16(in), 0x00), nil
	case 1:
		// htscodecs downgrades order 1 to order 0 below 8 bytes: the
		// four-way split has too little data to model a context.
		if len(in) < 8 {
			return frameRANS4x16(in, compressO0RANS4x16(in), 0x00), nil
		}
		return frameRANS4x16(in, compressO1RANS4x16(in, 0), 0x01), nil
	default:
		return nil, fmt.Errorf("rans4x16: order %d is not implemented (want 0 or 1)", order)
	}
}

// frameRANS4x16 wraps a rANS payload in the on-wire framing: a format
// byte and the big-endian raw-size varint. When the payload did not
// shrink the input it falls back to X_CAT (store verbatim), matching
// htscodecs' rans_compress_to_4x16.
func frameRANS4x16(in, payload []byte, format byte) []byte {
	if len(payload) >= len(in) {
		out := []byte{x4x16Cat}
		out = varPutU32(out, uint32(len(in)))
		return append(out, in...)
	}
	out := []byte{format}
	out = varPutU32(out, uint32(len(in)))
	return append(out, payload...)
}

// --- order-0 decode ----------------------------------------------------------

// uncompressO0RANS4x16 implements rans_uncompress_O0_4x16: in is the
// payload after the format byte and raw-size varint, rawSize is the
// declared decompressed length. The maxRANSRawSize ceiling is enforced
// here, not only in RANS4x16Decode, because the order-1 decoder reaches
// this function recursively with an attacker-controlled size when the
// frequency header is itself rANS-compressed.
func uncompressO0RANS4x16(in []byte, rawSize uint32) ([]byte, error) {
	if rawSize > maxRANSRawSize {
		return nil, fmt.Errorf("rans4x16: declared raw size %d exceeds the %d-byte safety ceiling",
			rawSize, maxRANSRawSize)
	}
	if len(in) < 16 {
		return nil, fmt.Errorf("rans4x16: order-0 payload %d bytes, need ≥16 for four states", len(in))
	}

	F, fsum, cp, err := decodeFreqRANS4x16(in)
	if err != nil {
		return nil, err
	}
	normaliseFreqShiftRANS4x16(&F, fsum, ransTotFreq)

	// Reverse lookup, indexed by the low TF_SHIFT bits of the state.
	var ssym [ransTotFreq]byte
	var sfreq [ransTotFreq]uint32
	var sbase [ransTotFreq]uint32
	x := uint32(0)
	for j := 0; j < 256; j++ {
		if F[j] == 0 {
			continue
		}
		if F[j] > ransTotFreq-x {
			return nil, fmt.Errorf("rans4x16: cumulative frequency overflow at symbol %d", j)
		}
		for y := uint32(0); y < F[j]; y++ {
			ssym[x+y] = byte(j)
			sfreq[x+y] = F[j]
			sbase[x+y] = y
		}
		x += F[j]
	}
	if x != ransTotFreq {
		return nil, fmt.Errorf("rans4x16: frequencies sum to %d, want %d", x, ransTotFreq)
	}

	var r [4]uint32
	for k := 0; k < 4; k++ {
		if cp+4 > len(in) {
			return nil, fmt.Errorf("rans4x16: truncated rANS state")
		}
		r[k] = uint32(in[cp]) | uint32(in[cp+1])<<8 | uint32(in[cp+2])<<16 | uint32(in[cp+3])<<24
		cp += 4
		if r[k] < rans4x16ByteL {
			return nil, fmt.Errorf("rans4x16: rANS state %d below normalisation bound", r[k])
		}
	}

	out := make([]byte, rawSize)
	mask := uint32(ransTotFreq - 1)
	for i := 0; i < int(rawSize); i++ {
		k := i & 3
		m := r[k] & mask
		out[i] = ssym[m]
		r[k] = sfreq[m]*(r[k]>>ransTFShift) + sbase[m]
		// 16-bit renorm: refill one word when the state drops below L.
		// Past end-of-input the state is left untouched, matching
		// htscodecs' RansDecRenormSafe — a valid stream never reaches
		// here, and a truncated one decodes to garbage either way.
		if r[k] < rans4x16ByteL && cp+1 < len(in) {
			w := uint32(in[cp]) | uint32(in[cp+1])<<8
			cp += 2
			r[k] = (r[k] << 16) | w
		}
	}
	return out, nil
}

// decodeFreqRANS4x16 implements decode_freq: it parses the delta-RLE
// alphabet then a big-endian varint per present symbol. It returns the
// frequency table, their sum, and the number of bytes consumed.
func decodeFreqRANS4x16(in []byte) (F [256]uint32, fsum uint32, n int, err error) {
	cp, err := decodeAlphabetRANS4x16(in, &F)
	if err != nil {
		return F, 0, 0, err
	}
	tot := uint32(0)
	for j := 0; j < 256; j++ {
		if F[j] == 0 {
			continue
		}
		var v uint32
		var ok bool
		v, cp, ok = varGetU32(in, cp)
		if !ok {
			return F, 0, 0, fmt.Errorf("rans4x16: truncated frequency varint")
		}
		F[j] = v
		tot += v
	}
	return F, tot, cp, nil
}

// decodeAlphabetRANS4x16 implements decode_alphabet: it marks F[sym]=1
// for every present symbol. Symbols are delta-RLE encoded — when a
// symbol byte is immediately followed by its successor value, that
// successor is read as a run start and the next byte as the run length.
// The table is terminated by a 0 byte.
func decodeAlphabetRANS4x16(in []byte, F *[256]uint32) (int, error) {
	if len(in) == 0 {
		return 0, fmt.Errorf("rans4x16: empty frequency table")
	}
	cp := 0
	rle := 0
	j := int(in[cp])
	cp++
	for {
		F[j] = 1
		if cp >= len(in) {
			return 0, fmt.Errorf("rans4x16: truncated alphabet")
		}
		if rle == 0 && j+1 == int(in[cp]) {
			j = int(in[cp])
			cp++
			if cp >= len(in) {
				return 0, fmt.Errorf("rans4x16: truncated alphabet run length")
			}
			rle = int(in[cp])
			cp++
		} else if rle != 0 {
			rle--
			j++
			if j > 255 {
				return 0, fmt.Errorf("rans4x16: alphabet run overflowed symbol 255")
			}
		} else {
			j = int(in[cp])
			cp++
		}
		if j == 0 {
			break
		}
	}
	return cp, nil
}

// --- order-0 encode ----------------------------------------------------------

// compressO0RANS4x16 implements rans_compress_O0_4x16. It returns the
// freq table + rANS bytes (no format byte / size — those are framing,
// added by RANS4x16Encode). An empty input yields an empty payload.
func compressO0RANS4x16(in []byte) []byte {
	if len(in) == 0 {
		return nil
	}

	var F [256]uint32
	for _, b := range in {
		F[b]++
	}

	// Normalise to a power-of-two ceiling, store that compact table,
	// then re-normalise the working copy up to TOTFREQ for coding.
	fsum := uint32(len(in))
	maxVal := round2u32(fsum)
	if maxVal > ransTotFreq {
		maxVal = ransTotFreq
	}
	normaliseFreqRANS4x16(&F, fsum, maxVal)
	tab := encodeFreqRANS4x16(&F)
	normaliseFreqRANS4x16(&F, maxVal, ransTotFreq)

	cum := make([]uint32, 256)
	c := uint32(0)
	for j := 0; j < 256; j++ {
		cum[j] = c
		c += F[j]
	}

	// rANS encodes in reverse, writing bytes back-to-front. Encode
	// order: the 0-3 tail bytes (states i-1..0), then the 4-way main
	// loop (states 3..0 per group, groups last-first), then the four
	// 4-byte state flushes (3..0). revBuf reverses it so the forward
	// result is [flushes][main renorm bytes][tail renorm bytes].
	rev := newRevBuf(len(in) + 64)
	var r [4]uint32
	for k := range r {
		r[k] = rans4x16ByteL
	}
	n := len(in)
	i := n & 3
	if i == 3 {
		s := in[n-1]
		r[2] = ransEncPutRANS4x16(r[2], rev, cum[s], F[s], ransTFShift)
	}
	if i >= 2 {
		s := in[n-(i-1)]
		r[1] = ransEncPutRANS4x16(r[1], rev, cum[s], F[s], ransTFShift)
	}
	if i >= 1 {
		s := in[n-i]
		r[0] = ransEncPutRANS4x16(r[0], rev, cum[s], F[s], ransTFShift)
	}
	for b := n &^ 3; b > 0; b -= 4 {
		s3, s2, s1, s0 := in[b-1], in[b-2], in[b-3], in[b-4]
		r[3] = ransEncPutRANS4x16(r[3], rev, cum[s3], F[s3], ransTFShift)
		r[2] = ransEncPutRANS4x16(r[2], rev, cum[s2], F[s2], ransTFShift)
		r[1] = ransEncPutRANS4x16(r[1], rev, cum[s1], F[s1], ransTFShift)
		r[0] = ransEncPutRANS4x16(r[0], rev, cum[s0], F[s0], ransTFShift)
	}
	for k := 3; k >= 0; k-- {
		rev.writeState(r[k])
	}

	return append(tab, rev.bytes()...)
}

// ransEncPutRANS4x16 encodes one symbol into the state, emitting a
// 16-bit renorm word into rev when the state would overflow the
// symbol's interval. shift is the frequency-table precision (12 for
// order-0, 10 or 12 for order-1). Mirrors RansEncRenorm + RansEncPut
// from rANS_word.h. The plain (divide-based) form is byte-identical to
// the reciprocal-based RansEncPutSymbol the reference encoder uses.
func ransEncPutRANS4x16(x uint32, rev *revBuf, start, freq uint32, shift uint) uint32 {
	xMax := ((uint32(rans4x16ByteL) >> shift) << 16) * freq
	if x >= xMax {
		// writeByte prepends, so emit the high byte first to land
		// [low, high] in forward order.
		rev.writeByte(byte(x >> 8))
		rev.writeByte(byte(x))
		x >>= 16
	}
	return ((x / freq) << shift) + (x % freq) + start
}

// encodeAlphabetRANS4x16 implements encode_alphabet: it writes the set
// of present symbols (F[j] != 0) as a delta-RLE list — a symbol byte,
// optionally followed by a run length when it abuts the previous
// symbol — terminated by a 0 byte.
func encodeAlphabetRANS4x16(F *[256]uint32) []byte {
	var cp []byte
	rle := 0
	for j := 0; j < 256; j++ {
		if F[j] == 0 {
			continue
		}
		if rle != 0 {
			rle--
		} else {
			cp = append(cp, byte(j))
			if j > 0 && F[j-1] != 0 {
				run := j + 1
				for run < 256 && F[run] != 0 {
					run++
				}
				run -= j + 1
				rle = run
				cp = append(cp, byte(run))
			}
		}
	}
	return append(cp, 0)
}

// encodeFreqRANS4x16 implements encode_freq: the delta-RLE alphabet
// (terminated by 0) followed by one big-endian varint per present
// symbol.
func encodeFreqRANS4x16(F *[256]uint32) []byte {
	cp := encodeAlphabetRANS4x16(F)
	for j := 0; j < 256; j++ {
		if F[j] != 0 {
			cp = varPutU32(cp, F[j])
		}
	}
	return cp
}

// --- frequency normalisation -------------------------------------------------

// round2u32 rounds v up to the next power of two (round2 in
// rANS_static16_int.h).
func round2u32(v uint32) uint32 {
	v--
	v |= v >> 1
	v |= v >> 2
	v |= v >> 4
	v |= v >> 8
	v |= v >> 16
	v++
	return v
}

// normaliseFreqRANS4x16 rescales F in place so the present frequencies
// sum exactly to tot, given they currently sum to size. It ports
// normalise_freq: an integer reciprocal scale, a one-shot "normalise
// harder" retry that compounds on the already-scaled array, and a
// final spread-the-deficit pass. Called twice by the encoder — first
// to the stored power-of-two ceiling, then up to TOTFREQ.
//
// htscodecs' normalise_freq returns -1 when the largest bucket would
// go non-positive; that path is unreachable here because every caller
// passes tot >= the count of present symbols, so the spread-the-deficit
// pass can always leave each bucket >= 1. The signal is therefore not
// propagated.
func normaliseFreqRANS4x16(F *[256]uint32, size, tot uint32) {
	if size == 0 {
		return
	}
	loop := 0
	for {
		tr := (uint64(tot)<<31)/uint64(size) + (uint64(1)<<30)/uint64(size)
		var newSize uint32
		m, M := 0, 0
		for j := 0; j < 256; j++ {
			if F[j] == 0 {
				continue
			}
			if uint32(m) < F[j] {
				m = int(F[j])
				M = j
			}
			F[j] = uint32((uint64(F[j]) * tr) >> 31)
			if F[j] == 0 {
				F[j] = 1
			}
			newSize += F[j]
		}
		size = newSize
		adjust := int(tot) - int(newSize)
		if adjust > 0 {
			F[M] += uint32(adjust)
		} else if adjust < 0 {
			if int(F[M]) > -adjust && (loop == 1 || int(F[M])/2 >= -adjust) {
				F[M] = uint32(int(F[M]) + adjust)
			} else if loop < 1 {
				loop++
				continue
			} else {
				adjust += int(F[M]) - 1
				F[M] = 1
				for j := 0; adjust != 0 && j < 256; j++ {
					if F[j] < 2 {
						continue
					}
					var mm int
					if int(F[j]) > -adjust {
						mm = adjust
					} else {
						mm = 1 - int(F[j])
					}
					F[j] = uint32(int(F[j]) + mm)
					adjust -= mm
				}
			}
		}
		return
	}
}

// normaliseFreqShiftRANS4x16 ports normalise_freq_shift: when the
// stored frequencies already sum to a power of two, scaling them up to
// maxTot is a single left shift.
func normaliseFreqShiftRANS4x16(F *[256]uint32, size, maxTot uint32) {
	if size == 0 || size == maxTot {
		return
	}
	shift := 0
	for size < maxTot {
		size *= 2
		shift++
	}
	for i := 0; i < 256; i++ {
		F[i] <<= shift
	}
}

// --- big-endian base-128 varints ---------------------------------------------

// varPutU32 appends i to cp as a big-endian base-128 varint (var_put_u32
// with BIG_END, the htscodecs default).
func varPutU32(cp []byte, i uint32) []byte {
	switch {
	case i < 1<<7:
		return append(cp, byte(i))
	case i < 1<<14:
		return append(cp, byte(i>>7)|128, byte(i&0x7f))
	case i < 1<<21:
		return append(cp, byte(i>>14)|128, byte((i>>7)&0x7f)|128, byte(i&0x7f))
	case i < 1<<28:
		return append(cp, byte(i>>21)|128, byte((i>>14)&0x7f)|128,
			byte((i>>7)&0x7f)|128, byte(i&0x7f))
	default:
		return append(cp, byte(i>>28)|128, byte((i>>21)&0x7f)|128,
			byte((i>>14)&0x7f)|128, byte((i>>7)&0x7f)|128, byte(i&0x7f))
	}
}

// varGetU32 reads a big-endian base-128 varint from in at offset cp. It
// returns the value, the offset past it, and false if the input is
// truncated. At most five bytes are consumed (a uint32 needs no more).
func varGetU32(in []byte, cp int) (uint32, int, bool) {
	var j uint32
	for n := 0; n < 5; n++ {
		if cp >= len(in) {
			return 0, cp, false
		}
		c := in[cp]
		cp++
		j = (j << 7) | uint32(c&0x7f)
		if c&0x80 == 0 {
			return j, cp, true
		}
	}
	return j, cp, true
}
