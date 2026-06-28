package codec

import (
	"fmt"
	"math"
	"sync"
)

// rANS 4x16 order-1 (the C2.1 slice). Ported from htscodecs
// rANS_static4x16pr.c (rans_compress_O1_4x16 / rans_uncompress_O1_4x16,
// rans_compute_shift) and rANS_static16_int.h (encode_freq1,
// encode_freq_d, decode_freq_d). The on-wire format is byte-identical to
// htscodecs — see rans4x16_test.go for the compliance vectors.
//
// Order-1 models each byte against its predecessor (context 0 for the
// first byte of each of the four interleaved quarters). The frequency
// table is therefore 256 per-context tables; htscodecs auto-tunes the
// table precision between 10 and 12 bits (rans_compute_shift) and may
// rANS-O0-compress the whole table when it exceeds 1000 bytes.
//
// Order-1 payload layout (the bytes after the framing format byte and
// raw-size varint):
//
//	[shift<<4 | hdrCompressed:1][frequency table][rANS bytes]
//
// shift is 10 or 12. When hdrCompressed is set the frequency table is
// itself an order-0 rANS stream, preceded by two varints (uncompressed
// and compressed table sizes).

const (
	tfShiftO1     = 12
	tfShiftO1Fast = 10
	totFreqO1     = 1 << tfShiftO1     // 4096
	totFreqO1Fast = 1 << tfShiftO1Fast // 1024

	// maxO1FreqTableSize bounds the uncompressed size of an order-1
	// frequency table. A legitimate table is at most ~128 KiB (256
	// contexts × at most 256 two-byte entries, plus the alphabet), so
	// the rANS-compressed-header path must reject a larger declared
	// size before it drives a decode allocation: the general
	// maxRANSRawSize ceiling (1 GiB) is far too loose for a frequency
	// table and lets a tiny crafted stream demand a ~1 GiB buffer.
	maxO1FreqTableSize = 1 << 20
)

// --- order-1 encode ----------------------------------------------------------

// compressO1RANS4x16 implements rans_compress_O1_4x16. It returns the
// order-1 payload (shift byte + frequency table + rANS bytes); the
// framing format byte and raw-size varint are added by the caller.
//
// forceShift is 0 for the normal auto-tuned precision; a test may pass
// 10 or 12 to pin the table precision and exercise the matching decode
// path directly. Production callers always pass 0.
func compressO1RANS4x16(in []byte, forceShift int) []byte {
	n := len(in)
	header, shift, cum, freq := encodeFreq1RANS4x16(in, forceShift, 4)
	sh := uint(shift)

	rev := newRevBuf(n + 64)
	var r [4]uint32
	for k := range r {
		r[k] = rans4x16ByteL
	}
	isz4 := n >> 2
	i0, i1, i2, i3 := isz4-2, 2*isz4-2, 3*isz4-2, 4*isz4-2
	l0, l1, l2 := in[i0+1], in[i1+1], in[i2+1]
	// Quarter 3 absorbs the remainder, so its last symbol is in[n-1].
	l3 := in[n-1]
	for x := n - 2; x > 4*isz4-2; x-- {
		c3 := in[x]
		r[3] = ransEncPutRANS4x16(r[3], rev, cum[c3][l3], freq[c3][l3], sh)
		l3 = c3
	}
	for ; i0 >= 0; i0, i1, i2, i3 = i0-1, i1-1, i2-1, i3-1 {
		c0, c1, c2, c3 := in[i0], in[i1], in[i2], in[i3]
		r[3] = ransEncPutRANS4x16(r[3], rev, cum[c3][l3], freq[c3][l3], sh)
		r[2] = ransEncPutRANS4x16(r[2], rev, cum[c2][l2], freq[c2][l2], sh)
		r[1] = ransEncPutRANS4x16(r[1], rev, cum[c1][l1], freq[c1][l1], sh)
		r[0] = ransEncPutRANS4x16(r[0], rev, cum[c0][l0], freq[c0][l0], sh)
		l0, l1, l2, l3 = c0, c1, c2, c3
	}
	// The first symbol of each quarter decodes with context 0.
	r[3] = ransEncPutRANS4x16(r[3], rev, cum[0][l3], freq[0][l3], sh)
	r[2] = ransEncPutRANS4x16(r[2], rev, cum[0][l2], freq[0][l2], sh)
	r[1] = ransEncPutRANS4x16(r[1], rev, cum[0][l1], freq[0][l1], sh)
	r[0] = ransEncPutRANS4x16(r[0], rev, cum[0][l0], freq[0][l0], sh)
	for k := 3; k >= 0; k-- {
		rev.writeState(r[k])
	}
	return append(header, rev.bytes()...)
}

// encodeFreq1RANS4x16 implements encode_freq1: it builds the order-1
// histogram, picks the table precision, encodes the frequency table
// (optionally rANS-O0-compressed) and returns it along with the chosen
// shift and the per-context cumulative-start / frequency tables the
// rANS coder needs. forceShift, when non-zero, overrides the auto-tuned
// precision (see compressO1RANS4x16).
//
// nway is the rANS state count — 4 for the 4x16 coder, 32 for the 32x16
// coder. It mirrors the C encode_freq1's Nway argument: it sets the
// quarter/thirty-second stride (isz4 = n/nway) and the number of
// context-0 phantom counts the first symbol of each interleaved part
// contributes.
//
// Deviation note. The per-context total T[i] is the exact frequency-row
// sum — the number of bytes that follow byte value i. Current htscodecs
// instead derives T[i] from the shared utils.h hist1_4, which adds one
// extra count for the input's final byte. That extra count perturbs the
// normalise_freq rounding for the final byte's context.
//
// This is a verified, deliberate divergence from current htscodecs
// *source*, chosen to match the vendored compliance *vectors*: the
// r4x16/q8.129 and q8.193 vectors are byte-exact reproducible only with
// T[i] == sum(F[i]); adding the final-byte count makes the encoder
// diverge from them. The vectors are static fixtures and htscodecs'
// own test harness never checks encode-equals-vector, so they can — and
// here demonstrably do — predate the current hist1_4. The vectors are
// this project's compliance oracle, so the encoder matches them.
//
// Consequence: for an input whose final byte value never recurs as a
// context, our encoder's output differs by a few bytes from current
// htscodecs. Both streams are valid and mutually decodable (our decoder
// reads genuine htscodecs output and vice versa); only encoder
// byte-equality with the latest C source is affected. The order-0
// context alphabet (present) is tracked separately so the final byte is
// always encodable even when T[final] is 0.
func encodeFreq1RANS4x16(in []byte, forceShift, nway int) (header []byte, shift int, cum, freq *[256][256]uint32) {
	n := len(in)
	isz4 := n / nway

	// Order-1 histogram: F[prev][cur]; T[i] is the count of i as a
	// context (the number of bytes that follow an i). htscodecs'
	// encode_freq1 normalises each context's row to T[i], so T must be
	// exactly the row sum — see the deviation note below.
	F := new([256][256]uint32)
	var T [256]uint32
	prev := 0
	for _, b := range in {
		F[prev][b]++
		prev = int(b)
	}
	for i := 0; i < 256; i++ {
		var tt uint32
		for j := 0; j < 256; j++ {
			tt += F[i][j]
		}
		T[i] += tt
	}
	// Phantom counts for the first symbol of interleaved parts 1..nway-1:
	// each decodes with context 0.
	for z := 1; z < nway; z++ {
		F[0][in[z*isz4]]++
	}
	T[0] += uint32(nway - 1)

	// present marks every symbol that occurs anywhere in the input; it
	// is the order-0 alphabet of contexts (htscodecs' present8). A
	// symbol that appears only as the input's final byte is a valid
	// symbol but never a context, so T[final] can be 0 while
	// present[final] must still be set. Context 0 is always present —
	// T[0] += 3 above guarantees it — so it needs no explicit mark.
	present := T
	present[in[n-1]] = 1

	// Order-0 alphabet of contexts, preceded by the (placeholder) shift
	// byte.
	body := []byte{0}
	body = append(body, encodeAlphabetRANS4x16(&present)...)

	var S [256]uint32
	shift = ransComputeShift(&present, F, &T, &S)
	if forceShift != 0 {
		shift = forceShift
	}
	sh := uint(shift)

	cum = new([256][256]uint32)
	freq = new([256][256]uint32)
	for i := 0; i < 256; i++ {
		if present[i] == 0 {
			continue
		}
		maxVal := S[i]
		if shift == tfShiftO1Fast && maxVal > totFreqO1Fast {
			maxVal = totFreqO1Fast
		}
		Fi := F[i] // copy of this context's row
		normaliseFreqRANS4x16(&Fi, T[i], maxVal)
		T[i] = maxVal
		body = encodeFreqDRANS4x16(body, &present, &Fi)
		normaliseFreqShiftRANS4x16(&Fi, T[i], 1<<sh)
		T[i] = 1 << sh
		c := uint32(0)
		for j := 0; j < 256; j++ {
			cum[i][j] = c
			freq[i][j] = Fi[j]
			c += Fi[j]
		}
	}
	body[0] = byte(shift << 4)

	// htscodecs rANS-O0-compresses the frequency table when it exceeds
	// 1000 bytes and the compressed form (plus its two size varints) is
	// genuinely smaller.
	if len(body) > 1000 {
		uFreq := body[1:]
		cFreq := compressO0RANS4x16(uFreq)
		if cFreq != nil && len(cFreq)+6 < len(body) {
			nb := []byte{byte(shift<<4) | 1}
			nb = varPutU32(nb, uint32(len(uFreq)))
			nb = varPutU32(nb, uint32(len(cFreq)))
			nb = append(nb, cFreq...)
			body = nb
		}
	}
	return body, shift, cum, freq
}

// ransComputeShift implements rans_compute_shift: it estimates the
// compressed entropy at 10-bit and 12-bit table precision and picks the
// smaller. S[i] receives the per-context power-of-two scale the stored
// frequencies are normalised to. F0 marks which contexts are present
// and T holds their totals (htscodecs passes the same array for both).
func ransComputeShift(F0 *[256]uint32, F *[256][256]uint32, T *[256]uint32, S *[256]uint32) int {
	e10, e12 := 0.0, 0.0
	maxTot := uint32(0)
	for i := 0; i < 256; i++ {
		if F0[i] == 0 {
			continue
		}
		maxVal := round2u32(T[i])
		sm10, sm12 := 0, 0
		for j := 0; j < 256; j++ {
			if F[i][j] != 0 && maxVal/F[i][j] > totFreqO1Fast {
				sm10++
			}
			if F[i][j] != 0 && maxVal/F[i][j] > totFreqO1 {
				sm12++
			}
		}
		l10 := math.Log(float64(totFreqO1Fast + sm10))
		l12 := math.Log(float64(totFreqO1 + sm12))
		tSlow := float64(totFreqO1) / float64(T[i])
		tFast := float64(totFreqO1Fast) / float64(T[i])
		ns := 0
		for j := 0; j < 256; j++ {
			if F[i][j] != 0 {
				ns++
				e10 -= float64(F[i][j]) * (fastLog(math.Max(float64(F[i][j])*tFast, 1)) - l10)
				e12 -= float64(F[i][j]) * (fastLog(math.Max(float64(F[i][j])*tSlow, 1)) - l12)
				e10 += 1.3
				e12 += 4.7
			}
		}
		if ns < 64 && maxVal > 128 {
			maxVal /= 2
		}
		if maxVal > 1024 {
			maxVal /= 2
		}
		if maxVal > totFreqO1 {
			maxVal = totFreqO1
		}
		S[i] = maxVal
		if maxTot < maxVal {
			maxTot = maxVal
		}
	}
	if e10/e12 < 1.01 || maxTot <= totFreqO1Fast {
		return tfShiftO1Fast
	}
	return tfShiftO1
}

// fastLog is htscodecs' fast approximate logarithm (utils.h): a bit
// reinterpretation of the IEEE-754 double. It is replicated exactly so
// the entropy estimate — and therefore the chosen shift — matches the
// reference encoder bit-for-bit.
func fastLog(a float64) float64 {
	x := int64(math.Float64bits(a))
	return float64(x-4606921278410026770) * 1.539095918623324e-16
}

// encodeFreqDRANS4x16 implements encode_freq_d: it writes one context's
// frequency row, restricted to the symbols that are valid contexts
// (F0[j] != 0). Absent symbols are run-length encoded as a 0 byte
// followed by run-length-minus-one.
func encodeFreqDRANS4x16(cp []byte, F0, F *[256]uint32) []byte {
	dz := 0
	for j := 0; j < 256; j++ {
		if F0[j] == 0 {
			continue
		}
		if F[j] != 0 {
			if dz != 0 {
				// Collapse the dz zero bytes just written into
				// [0][dz-1].
				cp = cp[:len(cp)-(dz-1)]
				cp = append(cp, byte(dz-1))
			}
			dz = 0
			cp = varPutU32(cp, F[j])
		} else {
			dz++
			cp = append(cp, 0)
		}
	}
	if dz != 0 {
		cp = cp[:len(cp)-(dz-1)]
		cp = append(cp, byte(dz-1))
	}
	return cp
}

// --- order-1 decode ----------------------------------------------------------

// ransO1Tables holds the per-context reverse-lookup tables an order-1
// decode builds (the symbol map plus the per-symbol frequency and
// cumulative-base rows). At the 12-bit table precision these are ~1.5 MiB
// in total — allocating them fresh on every order-1 block dominated the
// decoder's GC churn (issue #59).
//
// They are pooled and reused across decodes WITHOUT being zeroed on Get.
// This mirrors htscodecs, which keeps the equivalent tables in a
// per-thread scratch buffer and reuses them un-cleared: every present
// context's entire [0, 1<<shift) span is fully repopulated by the
// table-build loop before any decode reads it (the `x != 1<<sh`
// validation proves the span is filled), and the decoder only ever
// indexes a context whose symbols were emitted — i.e. cells written this
// call. Stale data from a previous decode is therefore never read, so
// the reuse is byte-exact. Do NOT add a clear/memset here: #54 showed the
// clear cost cancels the allocation saving exactly, defeating the point.
type ransO1Tables struct {
	sym  *[256][1 << tfShiftO1]byte
	freq *[256][256]uint32
	base *[256][256]uint32
}

// ransO1TablePool reuses the order-1 decode scratch tables across calls.
// sync.Pool is concurrency-safe, so the parallel (multi-thread) decode
// path draws one independent table set per in-flight decode. Shared by
// both the 4x16 and 32x16 order-1 decoders — their tables are identically
// shaped.
var ransO1TablePool = sync.Pool{New: func() any {
	return &ransO1Tables{
		sym:  new([256][1 << tfShiftO1]byte),
		freq: new([256][256]uint32),
		base: new([256][256]uint32),
	}
}}

// uncompressO1RANS4x16 implements rans_uncompress_O1_4x16: in is the
// payload after the framing format byte and raw-size varint, rawSize is
// the declared decompressed length.
func uncompressO1RANS4x16(in []byte, rawSize uint32) ([]byte, error) {
	if len(in) < 16 {
		return nil, fmt.Errorf("rans4x16: order-1 payload %d bytes, need ≥16 for four states", len(in))
	}
	shift := int(in[0] >> 4)
	if shift != tfShiftO1 && shift != tfShiftO1Fast {
		return nil, fmt.Errorf("rans4x16: order-1 table precision %d invalid (want 10 or 12)", shift)
	}
	sh := uint(shift)
	cp := 1

	// hdr is the buffer the frequency table is read from: the payload
	// itself, or a decompressed copy when the table was rANS-O0
	// compressed. tabEnd, when set, is where the rANS bytes resume in
	// the original payload.
	hdr := in
	hcp := cp
	tabEnd := -1
	if in[0]&1 != 0 {
		uFreqSz, c, ok := varGetU32(in, cp)
		if !ok {
			return nil, fmt.Errorf("rans4x16: truncated order-1 uncompressed-table size")
		}
		if uFreqSz > maxO1FreqTableSize {
			return nil, fmt.Errorf("rans4x16: order-1 frequency table size %d exceeds the %d-byte limit",
				uFreqSz, maxO1FreqTableSize)
		}
		cFreqSz, c2, ok := varGetU32(in, c)
		if !ok {
			return nil, fmt.Errorf("rans4x16: truncated order-1 compressed-table size")
		}
		cp = c2
		if int(cFreqSz) > len(in)-cp {
			return nil, fmt.Errorf("rans4x16: order-1 compressed table %d bytes exceeds payload", cFreqSz)
		}
		tabEnd = cp + int(cFreqSz)
		decoded, err := uncompressO0RANS4x16(in[cp:cp+int(cFreqSz)], uFreqSz)
		if err != nil {
			return nil, fmt.Errorf("rans4x16: order-1 table decompress: %w", err)
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
		return nil, fmt.Errorf("rans4x16: order-1 frequency table truncated after the alphabet")
	}

	// Per-context reverse-lookup tables, drawn from a pool and reused
	// without zeroing (see ransO1Tables).
	t := ransO1TablePool.Get().(*ransO1Tables)
	defer ransO1TablePool.Put(t)
	symTab := t.sym
	freqTab := t.freq
	baseTab := t.base
	for i := 0; i < 256; i++ {
		if F0[i] == 0 {
			continue
		}
		F, total, n, err := decodeFreqDRANS4x16(hdr, hcp, &F0)
		if err != nil {
			return nil, err
		}
		if n == 0 {
			return nil, fmt.Errorf("rans4x16: order-1 frequency table truncated at context %d", i)
		}
		hcp += n
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
				return nil, fmt.Errorf("rans4x16: order-1 context %d cumulative frequency overflow", i)
			}
			for y := uint32(0); y < F[j]; y++ {
				symTab[i][x+y] = byte(j)
			}
			freqTab[i][j] = F[j]
			baseTab[i][j] = x
			x += F[j]
		}
		if x != 1<<sh {
			return nil, fmt.Errorf("rans4x16: order-1 context %d frequencies sum to %d, want %d",
				i, x, 1<<sh)
		}
	}

	// The rANS bytes resume either where the table parse ended (plain
	// table) or past the compressed table (compressed table).
	dcp := hcp
	if tabEnd >= 0 {
		dcp = tabEnd
	}
	if dcp+16 > len(in) {
		return nil, fmt.Errorf("rans4x16: order-1 stream truncated before the four rANS states")
	}
	var r [4]uint32
	for k := 0; k < 4; k++ {
		r[k] = uint32(in[dcp]) | uint32(in[dcp+1])<<8 | uint32(in[dcp+2])<<16 | uint32(in[dcp+3])<<24
		dcp += 4
		if r[k] < rans4x16ByteL {
			return nil, fmt.Errorf("rans4x16: order-1 rANS state %d below normalisation bound", r[k])
		}
	}

	out := make([]byte, rawSize)
	mask := uint32((1 << sh) - 1)
	isz4 := int(rawSize) >> 2
	var l [4]int
	i4 := [4]int{0, isz4, 2 * isz4, 3 * isz4}
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
	for ; i4[0] < isz4; i4[0], i4[1], i4[2], i4[3] = i4[0]+1, i4[1]+1, i4[2]+1, i4[3]+1 {
		decode(0)
		decode(1)
		decode(2)
		decode(3)
	}
	for ; i4[3] < int(rawSize); i4[3]++ {
		decode(3)
	}
	return out, nil
}

// decodeFreqDRANS4x16 implements decode_freq_d: it reads one context's
// frequency row from hdr starting at cp, restricted to the symbols
// F0[j] != 0. It returns the row, its total, and the bytes consumed.
func decodeFreqDRANS4x16(hdr []byte, cp int, F0 *[256]uint32) (F [256]uint32, total uint32, n int, err error) {
	op := cp
	dz := 0
	for j := 0; j < 256 && cp < len(hdr); j++ {
		if F0[j] == 0 {
			continue
		}
		var f uint32
		if dz != 0 {
			dz--
		} else {
			var ok bool
			f, cp, ok = varGetU32(hdr, cp)
			if !ok {
				return F, 0, 0, fmt.Errorf("rans4x16: truncated order-1 frequency varint")
			}
			if f == 0 {
				if cp >= len(hdr) {
					return F, 0, 0, fmt.Errorf("rans4x16: truncated order-1 zero-run length")
				}
				dz = int(hdr[cp])
				cp++
			}
		}
		F[j] = f
		total += f
	}
	return F, total, cp - op, nil
}
