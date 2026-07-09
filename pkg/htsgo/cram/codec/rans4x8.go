package codec

import (
	"encoding/binary"
	"fmt"
	"sync"
)

// rANS 4x8 — the entropy coder used by CRAM v3.0 blocks.
//
// Ported from samtools/htscodecs rANS_static.c + rANS_byte.h. The
// on-wire format is byte-identical to htscodecs so the htscodecs test
// corpus (reference_code/htscodecs/tests/dat/r4x8/) is the compliance
// oracle — see rans4x8_test.go.
//
// Stream layout:
//
//	[order:1][compressedSize:uint32 LE][rawSize:uint32 LE][freq table][rANS bytes]
//
// `order` is 0 or 1; `compressedSize` counts the bytes after the
// 9-byte header (freq table + rANS bytes); `rawSize` is the
// decompressed length.

const (
	ransTFShift  = 12 // frequency-table precision: 12 bits.
	ransTotFreq  = 1 << ransTFShift
	ransByteL    = 1 << 23 // lower bound of the normalisation interval.
	ransHeaderSz = 9

	// maxRANSRawSize is a defensive ceiling on the decompressed size a
	// stream header may declare. rANS can legitimately expand highly
	// redundant data by a large factor (a 30-byte stream of one symbol
	// decodes to megabytes), so the ceiling can't be derived from the
	// compressed size — it has to be absolute. 1 GiB is orders of
	// magnitude above any real CRAM block (slices are typically ≤10 MB)
	// while still rejecting a malicious header that would otherwise
	// drive an unbounded allocation. htscodecs guards with INT_MAX;
	// this is the stricter, OOM-safe equivalent.
	maxRANSRawSize = 1 << 30
)

// RANS4x8Decode decompresses a complete rANS 4x8 stream (header
// included) and returns the raw bytes.
func RANS4x8Decode(in []byte) ([]byte, error) {
	if len(in) < ransHeaderSz {
		return nil, fmt.Errorf("rans4x8: input shorter than the 9-byte header (%d bytes)", len(in))
	}
	order := in[0]
	compSize := binary.LittleEndian.Uint32(in[1:5])
	rawSize := binary.LittleEndian.Uint32(in[5:9])
	if int(compSize) != len(in)-ransHeaderSz {
		return nil, fmt.Errorf("rans4x8: compressed-size header %d != actual payload %d",
			compSize, len(in)-ransHeaderSz)
	}
	if rawSize > maxRANSRawSize {
		return nil, fmt.Errorf("rans4x8: declared raw size %d exceeds the %d-byte safety ceiling",
			rawSize, maxRANSRawSize)
	}
	switch order {
	case 0:
		return ransDecodeO0(in, rawSize)
	case 1:
		return ransDecodeO1(in, rawSize)
	default:
		return nil, fmt.Errorf("rans4x8: unknown order byte %d (want 0 or 1)", order)
	}
}

// RANS4x8Encode compresses in with the given order (0 or 1) and
// returns a complete rANS 4x8 stream.
func RANS4x8Encode(in []byte, order int) ([]byte, error) {
	switch order {
	case 0:
		return ransEncodeO0(in), nil
	case 1:
		// Order-1 needs ≥4 bytes for the quarter split; htscodecs
		// falls back to order-0 below that threshold.
		if len(in) < 4 {
			return ransEncodeO0(in), nil
		}
		return ransEncodeO1(in), nil
	default:
		return nil, fmt.Errorf("rans4x8: unknown order %d (want 0 or 1)", order)
	}
}

// --- decode -----------------------------------------------------------------

// ransDecodeO0 implements rans_uncompress_O0.
func ransDecodeO0(in []byte, rawSize uint32) ([]byte, error) {
	cp := ransHeaderSz // cursor into in.

	// Reverse-lookup tables, indexed by the low TF_SHIFT bits of the
	// rANS state: symbol, frequency, and within-symbol base offset.
	var ssym [ransTotFreq + 32]byte
	var sfreq [ransTotFreq + 32]uint32
	var sbase [ransTotFreq + 32]uint32

	// Frequency-table parse. Symbols are delta-RLE encoded; freqs are
	// 12-bit var-length (high bit of the first byte = "another byte
	// follows"). The table is terminated by symbol 0.
	x := uint32(0)
	rle := 0
	if cp >= len(in) {
		return nil, fmt.Errorf("rans4x8: truncated frequency table")
	}
	j := int(in[cp])
	cp++
	for {
		if cp >= len(in) {
			return nil, fmt.Errorf("rans4x8: truncated frequency table")
		}
		f := int(in[cp])
		cp++
		if f >= 128 {
			if cp >= len(in) {
				return nil, fmt.Errorf("rans4x8: truncated frequency table")
			}
			f = ((f & 127) << 8) | int(in[cp])
			cp++
		}
		c := x
		if x+uint32(f) > ransTotFreq {
			return nil, fmt.Errorf("rans4x8: order-0 cumulative frequency overflow")
		}
		for y := uint32(0); y < uint32(f); y++ {
			ssym[y+c] = byte(j)
			sfreq[y+c] = uint32(f)
			sbase[y+c] = y
		}
		x += uint32(f)

		if cp >= len(in) {
			return nil, fmt.Errorf("rans4x8: truncated frequency table")
		}
		if rle == 0 && j+1 == int(in[cp]) {
			j = int(in[cp])
			cp++
			if cp >= len(in) {
				return nil, fmt.Errorf("rans4x8: truncated frequency table")
			}
			rle = int(in[cp])
			cp++
		} else if rle != 0 {
			rle--
			j++
			if j > 255 {
				return nil, fmt.Errorf("rans4x8: order-0 frequency-table symbol overflow")
			}
		} else {
			j = int(in[cp])
			cp++
		}
		if j == 0 {
			break
		}
	}
	if x < ransTotFreq-1 || x > ransTotFreq {
		return nil, fmt.Errorf("rans4x8: order-0 frequencies sum to %d, want %d", x, ransTotFreq)
	}
	if x != ransTotFreq {
		// htscodecs historically fills 4095 not 4096; pad the last slot.
		ssym[x] = ssym[x-1]
		sfreq[x] = sfreq[x-1]
		sbase[x] = sbase[x-1] + 1
	}

	var r [4]uint32
	for k := 0; k < 4; k++ {
		var err error
		r[k], cp, err = ransDecInit(in, cp)
		if err != nil {
			return nil, err
		}
	}

	out := make([]byte, rawSize)
	mask := uint32(ransTotFreq - 1)
	outEnd := int(rawSize) &^ 3
	for i := 0; i < outEnd; i += 4 {
		for k := 0; k < 4; k++ {
			m := r[k] & mask
			r[k] = sfreq[m]*(r[k]>>ransTFShift) + sbase[m]
			out[i+k] = ssym[m]
			r[k], cp = ransDecRenorm(r[k], in, cp)
		}
	}
	// Tail: 0-3 leftover bytes.
	for k := 0; k < int(rawSize)&3; k++ {
		out[outEnd+k] = ssym[r[k]&mask]
	}
	return out, nil
}

// ransDecodeO1 implements rans_uncompress_O1. The htscodecs `map`
// symbol-remapping is a cache-locality optimisation only; we index the
// per-context tables by the raw byte value directly, which is
// functionally identical.
func ransDecodeO1(in []byte, rawSize uint32) ([]byte, error) {
	cp := ransHeaderSz

	// D[ctx] is the TF_SHIFT-wide reverse lookup (cumfreq → symbol);
	// freq[ctx][sym] / start[ctx][sym] advance the state. D is
	// 256×4096 = 1 MiB — large for a local, but Go's growable stacks
	// handle it and keeping it on the stack avoids a per-call heap
	// allocation on the decode hot path.
	var D [256][ransTotFreq]byte
	var freq [256][256]uint32
	var start [256][256]uint32

	if cp >= len(in) {
		return nil, fmt.Errorf("rans4x8: truncated order-1 frequency table")
	}
	rleI := 0
	i := int(in[cp])
	cp++
	for {
		x := uint32(0)
		rleJ := 0
		if cp >= len(in) {
			return nil, fmt.Errorf("rans4x8: truncated order-1 frequency table")
		}
		j := int(in[cp])
		cp++
		for {
			if cp >= len(in) {
				return nil, fmt.Errorf("rans4x8: truncated order-1 frequency table")
			}
			f := int(in[cp])
			cp++
			if f >= 128 {
				if cp >= len(in) {
					return nil, fmt.Errorf("rans4x8: truncated order-1 frequency table")
				}
				f = ((f & 127) << 8) | int(in[cp])
				cp++
			}
			c := x
			if f == 0 {
				// A zero freq means the whole range — a
				// single-symbol context.
				f = ransTotFreq
			}
			start[i][j] = c
			freq[i][j] = uint32(f)
			if x+uint32(f) > ransTotFreq {
				return nil, fmt.Errorf("rans4x8: order-1 cumulative frequency overflow")
			}
			for y := uint32(0); y < uint32(f); y++ {
				D[i][c+y] = byte(j)
			}
			x += uint32(f)

			if cp >= len(in) {
				return nil, fmt.Errorf("rans4x8: truncated order-1 frequency table")
			}
			if rleJ == 0 && j+1 == int(in[cp]) {
				j = int(in[cp])
				cp++
				if cp >= len(in) {
					return nil, fmt.Errorf("rans4x8: truncated order-1 frequency table")
				}
				rleJ = int(in[cp])
				cp++
			} else if rleJ != 0 {
				rleJ--
				j++
				if j > 255 {
					return nil, fmt.Errorf("rans4x8: order-1 inner symbol overflow")
				}
			} else {
				j = int(in[cp])
				cp++
			}
			if j == 0 {
				break
			}
		}
		if x < ransTotFreq-1 || x > ransTotFreq {
			return nil, fmt.Errorf("rans4x8: order-1 context %d frequencies sum to %d", i, x)
		}
		if x < ransTotFreq {
			D[i][x] = D[i][x-1]
		}

		if cp >= len(in) {
			return nil, fmt.Errorf("rans4x8: truncated order-1 frequency table")
		}
		if rleI == 0 && i+1 == int(in[cp]) {
			i = int(in[cp])
			cp++
			if cp >= len(in) {
				return nil, fmt.Errorf("rans4x8: truncated order-1 frequency table")
			}
			rleI = int(in[cp])
			cp++
		} else if rleI != 0 {
			rleI--
			i++
			if i > 255 {
				return nil, fmt.Errorf("rans4x8: order-1 context symbol overflow")
			}
		} else {
			i = int(in[cp])
			cp++
		}
		if i == 0 {
			break
		}
	}

	var r [4]uint32
	for k := 0; k < 4; k++ {
		var err error
		r[k], cp, err = ransDecInit(in, cp)
		if err != nil {
			return nil, err
		}
	}

	out := make([]byte, rawSize)
	mask := uint32(ransTotFreq - 1)
	isz4 := int(rawSize) >> 2
	// Each of the 4 states decodes one quarter; ctx[k] is the previous
	// symbol in quarter k (0 to start).
	var ctx [4]uint32
	idx := [4]int{0, isz4, 2 * isz4, 3 * isz4}

	for ; idx[0] < isz4; idx[0], idx[1], idx[2], idx[3] = idx[0]+1, idx[1]+1, idx[2]+1, idx[3]+1 {
		for k := 0; k < 4; k++ {
			cc := D[ctx[k]][r[k]&mask]
			out[idx[k]] = cc
			m := r[k] & mask
			r[k] = freq[ctx[k]][cc]*(r[k]>>ransTFShift) + m - start[ctx[k]][cc]
			r[k], cp = ransDecRenorm(r[k], in, cp)
			ctx[k] = uint32(cc)
		}
	}
	// Remainder bytes belong to quarter 3 / state 3.
	for ; idx[3] < int(rawSize); idx[3]++ {
		cc := D[ctx[3]][r[3]&mask]
		out[idx[3]] = cc
		m := r[3] & mask
		r[3] = freq[ctx[3]][cc]*(r[3]>>ransTFShift) + m - start[ctx[3]][cc]
		r[3], cp = ransDecRenorm(r[3], in, cp)
		ctx[3] = uint32(cc)
	}
	return out, nil
}

// ransDecInit reads a 4-byte little-endian rANS state and validates it
// is above the normalisation lower bound.
func ransDecInit(in []byte, cp int) (uint32, int, error) {
	if cp+4 > len(in) {
		return 0, cp, fmt.Errorf("rans4x8: truncated rANS state")
	}
	x := uint32(in[cp]) | uint32(in[cp+1])<<8 | uint32(in[cp+2])<<16 | uint32(in[cp+3])<<24
	if x < ransByteL {
		return 0, cp, fmt.Errorf("rans4x8: rANS state %d below normalisation bound", x)
	}
	return x, cp + 4, nil
}

// ransDecRenorm refills the rANS state to at least ransByteL by
// shifting in input bytes. Past end-of-input it shifts in zeros, which
// matches htscodecs' RansDecRenormSafe tail behaviour.
func ransDecRenorm(x uint32, in []byte, cp int) (uint32, int) {
	for x < ransByteL {
		var b uint32
		if cp < len(in) {
			b = uint32(in[cp])
		}
		x = (x << 8) | b
		cp++
	}
	return x, cp
}

// --- encode -----------------------------------------------------------------

// normaliseFreqsO0 rescales raw symbol counts so they sum to
// ransTotFreq-1 (the "historically fill 4095" convention), mirroring
// rans_compress_O0. It works on a mutable copy of raw and — critically
// — the "normalise_harder" retry compounds on the already-scaled
// values exactly as the C does (F[j] is mutated in place there, and a
// goto re-runs the loop over the mutated array). m / M track the
// largest count and its symbol.
func normaliseFreqsO0(raw []int) []int {
	total := 0
	for _, v := range raw {
		total += v
	}
	f := append([]int(nil), raw...)
	var tr uint64
	if total != 0 {
		tr = (uint64(ransTotFreq)<<31)/uint64(total) + (uint64(1)<<30)/uint64(total)
	}
	for {
		fsum, m, M := 0, 0, 0
		for j := range f {
			if f[j] == 0 {
				continue
			}
			if m < f[j] {
				m, M = f[j], j
			}
			f[j] = int((uint64(f[j]) * tr) >> 31)
			if f[j] == 0 {
				f[j] = 1
			}
			fsum += f[j]
		}
		fsum++
		if fsum < ransTotFreq {
			f[M] += ransTotFreq - fsum
			return f
		}
		if fsum-ransTotFreq > f[M]/2 {
			tr = 2104533975 // ≈ *0.98 — the "normalise_harder" retry.
			continue
		}
		f[M] -= fsum - ransTotFreq
		return f
	}
}

// normaliseFreqsO1Into is the per-context normaliser for rans_compress_O1,
// writing into a caller-supplied dst (which must be at least len(raw)) instead
// of allocating a fresh row. Unlike the O0 path it scales with a floating-point
// factor p (and the retry sets p=0.98 and compounds), and its overshoot test
// uses `>=`. Callers only invoke it for contexts with a non-zero total. It lets
// the order-1 encoder reuse one scratch row across all 256 contexts. dst is
// fully overwritten from raw before any use, so its prior contents do not
// matter.
func normaliseFreqsO1Into(raw, dst []int) []int {
	total := 0
	for _, v := range raw {
		total += v
	}
	f := dst[:len(raw)]
	copy(f, raw)
	if total == 0 {
		return f
	}
	p := float64(ransTotFreq) / float64(total)
	for {
		t2, m, M := 0, 0, 0
		for j := range f {
			if f[j] == 0 {
				continue
			}
			if m < f[j] {
				m, M = f[j], j
			}
			f[j] = int(float64(f[j]) * p)
			if f[j] == 0 {
				f[j] = 1
			}
			t2 += f[j]
		}
		t2++
		if t2 < ransTotFreq {
			f[M] += ransTotFreq - t2
			return f
		}
		if t2-ransTotFreq >= f[M]/2 {
			p = 0.98
			continue
		}
		f[M] -= t2 - ransTotFreq
		return f
	}
}

// ransEncodeO0 implements rans_compress_O0.
func ransEncodeO0(in []byte) []byte {
	raw := make([]int, 256)
	for _, b := range in {
		raw[b]++
	}
	f := normaliseFreqsO0(raw)

	// Cumulative frequency starts per symbol.
	cum := make([]uint32, 256)
	c := uint32(0)
	for j := 0; j < 256; j++ {
		cum[j] = c
		c += uint32(f[j])
	}

	// Frequency table (delta-RLE on symbols, var-length 12-bit freqs).
	var tab []byte
	rle := 0
	for j := 0; j < 256; j++ {
		if f[j] == 0 {
			continue
		}
		if rle != 0 {
			rle--
		} else {
			tab = append(tab, byte(j))
			if j > 0 && f[j-1] != 0 {
				run := j + 1
				for run < 256 && f[run] != 0 {
					run++
				}
				run -= j + 1
				rle = run
				tab = append(tab, byte(run))
			}
		}
		if f[j] < 128 {
			tab = append(tab, byte(f[j]))
		} else {
			tab = append(tab, byte(128|(f[j]>>8)), byte(f[j]&0xff))
		}
	}
	tab = append(tab, 0)

	// rANS-encode in reverse, writing bytes backwards. The decoder
	// reads states 0..3 from the front then renorms forward, so the
	// encode call order must be: tail bytes (states rem-1..0), then
	// the main 4-way loop (states 3..0 per group, groups last-first),
	// then the four state flushes (3..0). revBuf writes back-to-front,
	// so the last thing written lands at the start of the body.
	rev := newRevBuf(len(in) + 64)
	var r [4]uint32
	for k := range r {
		r[k] = ransByteL
	}
	n := len(in)
	rem := n & 3
	// Tail: the last `rem` bytes go to states rem-1..0.
	for k := rem - 1; k >= 0; k-- {
		sym := in[(n&^3)+k]
		r[k] = ransEncPut(r[k], rev, cum[sym], uint32(f[sym]))
	}
	for i := n &^ 3; i > 0; i -= 4 {
		for k := 3; k >= 0; k-- {
			sym := in[i-4+k]
			r[k] = ransEncPut(r[k], rev, cum[sym], uint32(f[sym]))
		}
	}
	for k := 3; k >= 0; k-- {
		rev.writeState(r[k])
	}
	body := rev.bytes()

	return assembleRans(0, tab, body, len(in))
}

// ransO1Scratch holds the large per-block tables the order-1 encoder needs.
// The order-1 model is 256 contexts × 256 symbols, so a fresh encode would
// otherwise allocate ~1.5 MiB of frequency/cumulative tables per block; on a
// real CRAM the writer compresses thousands of blocks, so those allocations
// dominate the encoder's GC churn. Pooling one scratch per encode reuses that
// memory across blocks. Every field is zero-reset before use — the freq/cum
// tables are additive (built by += over histogram counts), so any stale entry
// left from a previous block would corrupt the frequency model and break the
// round-trip, hence reset() clears every table in full.
type ransO1Scratch struct {
	rawF [256][256]int
	f    [256][256]int
	rawT [256]int
	cum  [256][256]uint32
	// norm is a per-context scratch for normaliseFreqsO1Into so it does not
	// allocate a fresh row on every context.
	norm [256]int
}

// reset zeroes every table so a pooled scratch carries no stale counts from a
// previous block. The zeroing is critical: rawF/rawT accumulate histogram
// counts and cum/f are rebuilt from them, so a single stale non-zero entry
// would silently corrupt the frequency model.
func (s *ransO1Scratch) reset() {
	s.rawF = [256][256]int{}
	s.f = [256][256]int{}
	s.rawT = [256]int{}
	s.cum = [256][256]uint32{}
	// norm is fully overwritten by normaliseFreqsO1Into before every read,
	// so it needs no reset.
}

// ransO1ScratchPool recycles order-1 encoder scratch across blocks.
var ransO1ScratchPool = sync.Pool{New: func() any { return new(ransO1Scratch) }}

// ransEncodeO1 implements rans_compress_O1: a per-previous-byte context
// model over four interleaved quarters.
func ransEncodeO1(in []byte) []byte {
	n := len(in)
	sc := ransO1ScratchPool.Get().(*ransO1Scratch)
	sc.reset()
	defer ransO1ScratchPool.Put(sc)

	// Per-context raw counts.
	rawF := &sc.rawF
	rawT := &sc.rawT
	// Context 0 is the previous byte; the first byte of the stream has
	// context 0.
	prev := 0
	for _, b := range in {
		rawF[prev][b]++
		rawT[prev]++
		prev = int(b)
	}
	// htscodecs adds 3 phantom context-0 counts for the first symbol of
	// quarters 1,2,3 (they decode with starting context 0).
	isz4 := n >> 2
	for q := 1; q <= 3; q++ {
		s := in[q*isz4]
		rawF[0][s]++
		rawT[0]++
	}

	f := &sc.f
	cum := &sc.cum
	var tab []byte
	rleI := 0
	for i := 0; i < 256; i++ {
		if rawT[i] == 0 {
			continue
		}
		row := normaliseFreqsO1Into(rawF[i][:], sc.norm[:])
		copy(f[i][:], row)
		c := uint32(0)
		for j := 0; j < 256; j++ {
			cum[i][j] = c
			c += uint32(f[i][j])
		}
		// Context-symbol delta-RLE.
		if rleI != 0 {
			rleI--
		} else {
			tab = append(tab, byte(i))
			if i > 0 && rawT[i-1] != 0 {
				run := i + 1
				for run < 256 && rawT[run] != 0 {
					run++
				}
				run -= i + 1
				rleI = run
				tab = append(tab, byte(run))
			}
		}
		// Inner table for context i.
		rleJ := 0
		for j := 0; j < 256; j++ {
			if f[i][j] == 0 {
				continue
			}
			if rleJ != 0 {
				rleJ--
			} else {
				tab = append(tab, byte(j))
				if j > 0 && f[i][j-1] != 0 {
					run := j + 1
					for run < 256 && f[i][run] != 0 {
						run++
					}
					run -= j + 1
					rleJ = run
					tab = append(tab, byte(run))
				}
			}
			if f[i][j] < 128 {
				tab = append(tab, byte(f[i][j]))
			} else {
				tab = append(tab, byte(128|(f[i][j]>>8)), byte(f[i][j]&0xff))
			}
		}
		tab = append(tab, 0)
	}
	tab = append(tab, 0)

	// Encode in reverse over the 4 quarters.
	rev := newRevBuf(n + 64)
	var r [4]uint32
	for k := range r {
		r[k] = ransByteL
	}
	i0, i1, i2, i3 := 1*isz4-2, 2*isz4-2, 3*isz4-2, 4*isz4-2
	l0, l1, l2, l3 := in[i0+1], in[i1+1], in[i2+1], in[i3+1]
	// Remainder past 4*isz4 belongs to state 3.
	l3 = in[n-1]
	for x := n - 2; x > 4*isz4-2; x-- {
		c3 := in[x]
		r[3] = ransEncPut(r[3], rev, cum[c3][l3], uint32(f[c3][l3]))
		l3 = c3
	}
	for ; i0 >= 0; i0, i1, i2, i3 = i0-1, i1-1, i2-1, i3-1 {
		c0, c1, c2, c3 := in[i0], in[i1], in[i2], in[i3]
		r[3] = ransEncPut(r[3], rev, cum[c3][l3], uint32(f[c3][l3]))
		r[2] = ransEncPut(r[2], rev, cum[c2][l2], uint32(f[c2][l2]))
		r[1] = ransEncPut(r[1], rev, cum[c1][l1], uint32(f[c1][l1]))
		r[0] = ransEncPut(r[0], rev, cum[c0][l0], uint32(f[c0][l0]))
		l0, l1, l2, l3 = c0, c1, c2, c3
	}
	// Final symbols use context 0.
	r[3] = ransEncPut(r[3], rev, cum[0][l3], uint32(f[0][l3]))
	r[2] = ransEncPut(r[2], rev, cum[0][l2], uint32(f[0][l2]))
	r[1] = ransEncPut(r[1], rev, cum[0][l1], uint32(f[0][l1]))
	r[0] = ransEncPut(r[0], rev, cum[0][l0], uint32(f[0][l0]))
	for k := 3; k >= 0; k-- {
		rev.writeState(r[k])
	}
	body := rev.bytes()

	return assembleRans(1, tab, body, n)
}

// assembleRans builds the final 9-byte header + freq table + body.
func assembleRans(order byte, tab, body []byte, rawSize int) []byte {
	compSize := len(tab) + len(body)
	out := make([]byte, ransHeaderSz+compSize)
	out[0] = order
	binary.LittleEndian.PutUint32(out[1:5], uint32(compSize))
	binary.LittleEndian.PutUint32(out[5:9], uint32(rawSize))
	copy(out[9:], tab)
	copy(out[9+len(tab):], body)
	return out
}

// ransEncPut encodes one symbol (cumulative start, frequency) into the
// state, renormalising into rev. Mirrors rANS_byte.h's RansEncRenorm +
// RansEncPut.
func ransEncPut(x uint32, rev *revBuf, start, freq uint32) uint32 {
	// Renormalise: emit low bytes until x fits the symbol's interval.
	xMax := ((uint32(ransByteL) >> ransTFShift) << 8) * freq
	for x >= xMax {
		rev.writeByte(byte(x & 0xff))
		x >>= 8
	}
	return ((x / freq) << ransTFShift) + (x % freq) + start
}

// revBuf accumulates bytes written back-to-front (rANS emits its output
// in reverse). bytes() returns them in forward order.
type revBuf struct {
	buf []byte
	pos int // next write index, decremented as we go.
}

func newRevBuf(capHint int) *revBuf {
	if capHint < 64 {
		capHint = 64
	}
	b := make([]byte, capHint)
	return &revBuf{buf: b, pos: capHint}
}

func (r *revBuf) ensure(n int) {
	if r.pos >= n {
		return
	}
	grow := len(r.buf)
	if grow < n {
		grow = n
	}
	nb := make([]byte, len(r.buf)+grow)
	copy(nb[len(nb)-len(r.buf):], r.buf)
	r.pos += len(nb) - len(r.buf)
	r.buf = nb
}

func (r *revBuf) writeByte(b byte) {
	r.ensure(1)
	r.pos--
	r.buf[r.pos] = b
}

// writeState writes a 4-byte rANS state little-endian, back-to-front so
// the forward-order result matches RansEncFlush.
func (r *revBuf) writeState(x uint32) {
	r.ensure(4)
	r.pos -= 4
	r.buf[r.pos+0] = byte(x)
	r.buf[r.pos+1] = byte(x >> 8)
	r.buf[r.pos+2] = byte(x >> 16)
	r.buf[r.pos+3] = byte(x >> 24)
}

func (r *revBuf) bytes() []byte {
	return r.buf[r.pos:]
}
