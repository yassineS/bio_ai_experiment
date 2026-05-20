package codec

import "fmt"

// rANS Nx16 (32-way) — the C2.3 slice. Ported from samtools/htscodecs
// rANS_static32x16pr.c (rans_compress_O0_32x16 / rans_uncompress_O0_32x16,
// rans_compress_O1_32x16 / rans_uncompress_O1_32x16, NX==32). The on-wire
// format is byte-identical to htscodecs — see rans4x16_test.go for the
// compliance vectors (the .4 / .5 suffixes).
//
// The 32-way coder is a distinct on-wire layout selected by the X_32
// format bit (0x04): 32 rANS states are initialised, interleaved and
// flushed instead of the 4x16 coder's 4. Everything else is shared with
// the 4x16 coder — the same delta-RLE frequency-table format, the same
// normalise_freq / normalise_freq_shift, the same 12-bit table precision
// and the same 16-bit ("word") renormalisation. Only the scalar
// (non-SIMD) NX-way path is ported; htscodecs' AVX2/AVX512/NEON variants
// are byte-identical to it, just faster.
//
// Order-0 layout differs from the 4x16 form only in the state count: the
// payload is [freq table][32 little-endian state words emitted last,
// state 0 frontmost]. Order-1 likewise interleaves 32 states; the
// frequency table is built by the shared encodeFreq1RANS4x16 with the
// Nway argument set to 32 (it changes the n/NX stride and the count of
// context-0 phantom symbols).

// ransNX is the 32-way state count (NX in rANS_static32x16pr.c).
const ransNX = 32

// --- framing -----------------------------------------------------------------

// frameRANS4x16X32 wraps a 32-way rANS payload in the on-wire framing: a
// format byte (with the X_32 bit set) and the big-endian raw-size
// varint. As with frameRANS4x16 it falls back to X_CAT when the payload
// did not shrink the input.
func frameRANS4x16X32(in, payload []byte, format byte) []byte {
	if len(payload) >= len(in) {
		out := []byte{x4x16Cat}
		out = varPutU32(out, uint32(len(in)))
		return append(out, in...)
	}
	out := []byte{format | x4x16X32}
	out = varPutU32(out, uint32(len(in)))
	return append(out, payload...)
}

// encodeRANS4x16X32 produces a complete 32-way rANS stream for the plain
// (no PACK/RLE/STRIPE) order-0 / order-1 paths. order is 0 or 1.
//
// htscodecs forces order-1 down to order-0 when the input is shorter
// than NX bytes: rans_compress_O1_32x16 returns NULL for in_size < NX
// and the caller retries at order 0.
func encodeRANS4x16X32(in []byte, order int) []byte {
	if order&1 != 0 && len(in) >= ransNX {
		return frameRANS4x16X32(in, compressO1RANS4x16X32(in, 0), 0x01)
	}
	return frameRANS4x16X32(in, compressO0RANS4x16X32(in), 0x00)
}

// --- order-0 encode ----------------------------------------------------------

// compressO0RANS4x16X32 implements rans_compress_O0_32x16. It returns the
// freq table + rANS bytes (no format byte / size — those are framing).
// An empty input yields an empty payload.
func compressO0RANS4x16X32(in []byte) []byte {
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

	// rANS encodes in reverse. The C encoder writes (ptr decrementing):
	// first the 0..NX-1 tail symbols for states i-1..0, then the NX-way
	// main loop (states NX-1..0 per group, groups last-first), then the
	// 32 four-byte state flushes (NX-1..0). revBuf prepends every write,
	// so writing in that same order yields the forward result
	// [state0..state31][main renorm bytes][tail renorm bytes].
	rev := newRevBuf(len(in) + 256)
	var r [ransNX]uint32
	for k := range r {
		r[k] = rans4x16ByteL
	}
	n := len(in)
	i := n & (ransNX - 1)
	// Tail: state z (z from i-1 down to 0) encodes in[n-i+z].
	for z := i - 1; z >= 0; z-- {
		s := in[n-i+z]
		r[z] = ransEncPutRANS4x16(r[z], rev, cum[s], F[s], ransTFShift)
	}
	// Main loop: i runs high to low in steps of NX; state z encodes
	// in[i-NX+z].
	for b := n &^ (ransNX - 1); b > 0; b -= ransNX {
		for z := ransNX - 1; z >= 0; z-- {
			s := in[b-ransNX+z]
			r[z] = ransEncPutRANS4x16(r[z], rev, cum[s], F[s], ransTFShift)
		}
	}
	for z := ransNX - 1; z >= 0; z-- {
		rev.writeState(r[z])
	}

	return append(tab, rev.bytes()...)
}

// --- order-0 decode ----------------------------------------------------------

// uncompressO0RANS4x16X32 implements rans_uncompress_O0_32x16: in is the
// payload after the format byte and raw-size varint, rawSize is the
// declared decompressed length.
func uncompressO0RANS4x16X32(in []byte, rawSize uint32) ([]byte, error) {
	if rawSize > maxRANSRawSize {
		return nil, fmt.Errorf("rans4x16: declared raw size %d exceeds the %d-byte safety ceiling",
			rawSize, maxRANSRawSize)
	}
	if len(in) < 16 {
		return nil, fmt.Errorf("rans4x16: 32-way order-0 payload %d bytes, need ≥16", len(in))
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

	if len(in)-cp < ransNX*4 {
		return nil, fmt.Errorf("rans4x16: 32-way order-0 stream truncated before the %d rANS states", ransNX)
	}
	var r [ransNX]uint32
	for k := 0; k < ransNX; k++ {
		r[k] = uint32(in[cp]) | uint32(in[cp+1])<<8 | uint32(in[cp+2])<<16 | uint32(in[cp+3])<<24
		cp += 4
		if r[k] < rans4x16ByteL {
			return nil, fmt.Errorf("rans4x16: 32-way rANS state %d below normalisation bound", r[k])
		}
	}

	out := make([]byte, rawSize)
	mask := uint32(ransTotFreq - 1)
	outEnd := int(rawSize) &^ (ransNX - 1)
	for i := 0; i < outEnd; i += ransNX {
		for z := 0; z < ransNX; z++ {
			m := r[z] & mask
			out[i+z] = ssym[m]
			r[z] = sfreq[m]*(r[z]>>ransTFShift) + sbase[m]
			// 16-bit renorm: refill one word when the state drops
			// below L. Past end-of-input the state is left untouched,
			// matching htscodecs' RansDecRenormSafe.
			if r[z] < rans4x16ByteL && cp+1 < len(in) {
				w := uint32(in[cp]) | uint32(in[cp+1])<<8
				cp += 2
				r[z] = (r[z] << 16) | w
			}
		}
	}
	// Remainder: states 0..(rawSize&31)-1, no further renorm (the C
	// reads `out[out_end+z] = s3[R[z]&mask]` for the tail).
	for z := 0; z < int(rawSize)&(ransNX-1); z++ {
		out[outEnd+z] = ssym[r[z]&mask]
	}
	return out, nil
}

// --- order-1 encode ----------------------------------------------------------

// compressO1RANS4x16X32 implements rans_compress_O1_32x16. It returns the
// order-1 payload (shift byte + frequency table + rANS bytes); the
// framing format byte and raw-size varint are added by the caller.
//
// forceShift is 0 for the normal auto-tuned precision; a test may pass
// 10 or 12 to pin the table precision. Production callers pass 0.
//
// The caller must ensure len(in) >= NX; rans_compress_O1_32x16 returns
// NULL otherwise (handled in encodeRANS4x16X32 and compressToRANS4x16X32).
func compressO1RANS4x16X32(in []byte, forceShift int) []byte {
	n := len(in)
	header, shift, cum, freq := encodeFreq1RANS4x16(in, forceShift, ransNX)
	sh := uint(shift)

	rev := newRevBuf(n + 256)
	var r [ransNX]uint32
	for k := range r {
		r[k] = rans4x16ByteL
	}

	// Per-state cursors. State z scans the input region [z*isz4, ...)
	// and starts two before its quarter-of-32 boundary; the C sets
	// iN[z] = (z+1)*isz4-2 and lN[z] = in[iN[z]+1].
	isz4 := n / ransNX
	var iN [ransNX]int
	var lN [ransNX]byte
	for z := 0; z < ransNX; z++ {
		iN[z] = (z+1)*isz4 - 2
		lN[z] = in[iN[z]+1]
	}

	// Remainder: state NX-1 absorbs everything past NX*isz4.
	z := ransNX - 1
	lN[z] = in[n-1]
	for iN[z] = n - 2; iN[z] > ransNX*isz4-2; iN[z]-- {
		c := in[iN[z]]
		r[z] = ransEncPutRANS4x16(r[z], rev, cum[c][lN[z]], freq[c][lN[z]], sh)
		lN[z] = c
	}

	// Main loop: every state walks its cursor down to 0. The C drives
	// it by `i32[0] >= in`, i.e. while state 0's cursor is still inside
	// the input.
	for iN[0] >= 0 {
		for z := ransNX - 1; z >= 0; z-- {
			c := in[iN[z]]
			l := lN[z]
			r[z] = ransEncPutRANS4x16(r[z], rev, cum[c][l], freq[c][l], sh)
			lN[z] = c
			iN[z]--
		}
	}

	// The first symbol of each interleaved part decodes with context 0.
	for z := ransNX - 1; z >= 0; z-- {
		r[z] = ransEncPutRANS4x16(r[z], rev, cum[0][lN[z]], freq[0][lN[z]], sh)
	}
	for z := ransNX - 1; z >= 0; z-- {
		rev.writeState(r[z])
	}
	return append(header, rev.bytes()...)
}

// --- order-1 decode ----------------------------------------------------------

// uncompressO1RANS4x16X32 implements rans_uncompress_O1_32x16: in is the
// payload after the framing format byte and raw-size varint, rawSize is
// the declared decompressed length.
func uncompressO1RANS4x16X32(in []byte, rawSize uint32) ([]byte, error) {
	if len(in) < ransNX*4 {
		return nil, fmt.Errorf("rans4x16: 32-way order-1 payload %d bytes, need ≥%d", len(in), ransNX*4)
	}
	if rawSize > maxRANSRawSize {
		return nil, fmt.Errorf("rans4x16: declared raw size %d exceeds the %d-byte safety ceiling",
			rawSize, maxRANSRawSize)
	}
	shift := int(in[0] >> 4)
	if shift != tfShiftO1 && shift != tfShiftO1Fast {
		return nil, fmt.Errorf("rans4x16: 32-way order-1 table precision %d invalid (want 10 or 12)", shift)
	}
	sh := uint(shift)
	cp := 1

	// hdr is the buffer the frequency table is read from: the payload
	// itself, or a decompressed copy when the table was rANS-O0
	// compressed. tabEnd, when set, is where the rANS bytes resume.
	hdr := in
	hcp := cp
	tabEnd := -1
	if in[0]&1 != 0 {
		uFreqSz, c, ok := varGetU32(in, cp)
		if !ok {
			return nil, fmt.Errorf("rans4x16: truncated 32-way order-1 uncompressed-table size")
		}
		cFreqSz, c2, ok := varGetU32(in, c)
		if !ok {
			return nil, fmt.Errorf("rans4x16: truncated 32-way order-1 compressed-table size")
		}
		cp = c2
		if int(cFreqSz) > len(in)-cp {
			return nil, fmt.Errorf("rans4x16: 32-way order-1 compressed table %d bytes exceeds payload", cFreqSz)
		}
		tabEnd = cp + int(cFreqSz)
		decoded, err := uncompressO0RANS4x16(in[cp:cp+int(cFreqSz)], uFreqSz)
		if err != nil {
			return nil, fmt.Errorf("rans4x16: 32-way order-1 table decompress: %w", err)
		}
		hdr = decoded
		hcp = 0
	}

	// Order-0 alphabet of contexts.
	var F0 [256]uint32
	consumed, err := decodeAlphabetRANS4x16(hdr[hcp:], &F0)
	if err != nil {
		return nil, err
	}
	hcp += consumed
	if hcp >= len(hdr) {
		return nil, fmt.Errorf("rans4x16: 32-way order-1 frequency table truncated after the alphabet")
	}

	// Per-context reverse-lookup tables.
	symTab := new([256][1 << tfShiftO1]byte)
	freqTab := new([256][256]uint32)
	baseTab := new([256][256]uint32)
	for i := 0; i < 256; i++ {
		if F0[i] == 0 {
			continue
		}
		F, total, nb, err := decodeFreqDRANS4x16(hdr, hcp, &F0)
		if err != nil {
			return nil, err
		}
		if nb == 0 {
			return nil, fmt.Errorf("rans4x16: 32-way order-1 frequency table truncated at context %d", i)
		}
		hcp += nb
		if total == 0 {
			continue
		}
		normaliseFreqShiftRANS4x16(&F, total, 1<<sh)
		x := uint32(0)
		for j := 0; j < 256; j++ {
			if F[j] == 0 {
				continue
			}
			if F[j] > (1<<sh)-x {
				return nil, fmt.Errorf("rans4x16: 32-way order-1 context %d cumulative frequency overflow", i)
			}
			for y := uint32(0); y < F[j]; y++ {
				symTab[i][x+y] = byte(j)
			}
			freqTab[i][j] = F[j]
			baseTab[i][j] = x
			x += F[j]
		}
		if x != 1<<sh {
			return nil, fmt.Errorf("rans4x16: 32-way order-1 context %d frequencies sum to %d, want %d",
				i, x, 1<<sh)
		}
	}

	// The rANS bytes resume either where the table parse ended (plain
	// table) or past the compressed table (compressed table).
	dcp := hcp
	if tabEnd >= 0 {
		dcp = tabEnd
	}
	if dcp+ransNX*4 > len(in) {
		return nil, fmt.Errorf("rans4x16: 32-way order-1 stream truncated before the %d rANS states", ransNX)
	}
	var r [ransNX]uint32
	for k := 0; k < ransNX; k++ {
		r[k] = uint32(in[dcp]) | uint32(in[dcp+1])<<8 | uint32(in[dcp+2])<<16 | uint32(in[dcp+3])<<24
		dcp += 4
		if r[k] < rans4x16ByteL {
			return nil, fmt.Errorf("rans4x16: 32-way order-1 rANS state %d below normalisation bound", r[k])
		}
	}

	out := make([]byte, rawSize)
	mask := uint32((1 << sh) - 1)
	isz4 := int(rawSize) / ransNX
	var l [ransNX]int
	var i4 [ransNX]int
	for z := 0; z < ransNX; z++ {
		i4[z] = z * isz4
	}
	decode := func(k int) {
		m := r[k] & mask
		c := symTab[l[k]][m]
		out[i4[k]] = c
		r[k] = freqTab[l[k]][c]*(r[k]>>sh) + m - baseTab[l[k]][c]
		l[k] = int(c)
		if r[k] < rans4x16ByteL && dcp+1 < len(in) {
			w := uint32(in[dcp]) | uint32(in[dcp+1])<<8
			dcp += 2
			r[k] = (r[k] << 16) | w
		}
	}
	for i4[0] < isz4 {
		for z := 0; z < ransNX; z++ {
			decode(z)
			i4[z]++
		}
	}
	// Remainder: state NX-1 absorbs the tail.
	for ; i4[ransNX-1] < int(rawSize); i4[ransNX-1]++ {
		decode(ransNX - 1)
	}
	return out, nil
}
