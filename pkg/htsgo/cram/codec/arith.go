package codec

import "fmt"

// arith_dynamic — the adaptive range coder used by CRAM v3.1 blocks
// (compression method 6).
//
// Ported from samtools/htscodecs c_range_coder.h (the carry-less byte-wise
// range coder), c_simple_model.h (the adaptive frequency model) and
// arith_dynamic.c (arith_compress_to / arith_uncompress_to). The on-wire
// format is byte-identical to htscodecs, so the htscodecs test corpus
// (reference_code/htscodecs/tests/dat/arith/) serves as the compliance
// oracle — see arith_test.go.
//
// Unlike rANS, which writes a static frequency table up front, arith_dynamic
// adapts its model as it codes: every symbol bumps the coded symbol's
// frequency and periodically the model is renormalised. Encoder and decoder
// run the identical model update so they stay in lockstep.
//
// The format byte and the PACK/RLE/STRIPE/CAT/NOSZ transform layer are
// shared with rANS 4x16 — arith_dynamic reuses the same X_* bit definitions
// (x4x16Pack, x4x16RLE, …), the same big-endian varints (varPutU32 /
// varGetU32) and the same hts_pack / unstripe helpers. Only the entropy
// core differs.
//
// Stream layout (order-0 / order-1, non-transform):
//
//	[formatByte:1][rawSize:varint][maxSym:1][range-coder bytes]
//
// maxSym is one past the largest byte value present in the input (0 means
// 256). The decoder seeds its adaptive model with maxSym symbols.

// --- range coder -------------------------------------------------------------

// rcTop is the renormalisation threshold of the range coder: whenever the
// range drops below 1<<24 a byte is shifted out (encode) or in (decode).
const rcTop = 1 << 24

// rcThres is 255*TOP, the carry threshold in RC_ShiftLow.
const rcThres = 255 * uint32(rcTop)

// rangeEncoder is the carry-less byte-wise range encoder from
// c_range_coder.h. It buffers a pending top byte (cache) plus a count of
// consecutive 0xFF bytes so a late carry can ripple through them.
type rangeEncoder struct {
	low   uint32
	rng   uint32
	ffNum uint32 // count of pending 0xFF bytes
	cache uint32 // top byte of low, awaiting flush
	carry uint32 // 1 when the buffered bytes must be incremented
	out   []byte
}

// newRangeEncoder returns a range encoder ready for RC_StartEncode.
func newRangeEncoder() *rangeEncoder {
	return &rangeEncoder{rng: 0xFFFFFFFF}
}

// shiftLow ports RC_ShiftLow: it emits the buffered top byte (plus any
// carry), flushes pending 0xFF bytes, and shifts low up by eight bits.
// RC_StartEncode leaves cache=0, so — exactly like the C — the first byte
// emitted is always 0; the matching RC_StartDecode reads it back as the
// first of five priming bytes, so we keep it on the wire verbatim.
func (e *rangeEncoder) shiftLow() {
	if e.low < rcThres || e.carry != 0 {
		e.out = append(e.out, byte(e.cache+e.carry))
		for e.ffNum != 0 {
			e.out = append(e.out, byte(e.carry-1))
			e.ffNum--
		}
		e.cache = e.low >> 24
		e.carry = 0
	} else {
		e.ffNum++
	}
	e.low <<= 8
}

// encode ports RC_Encode: it narrows the coding interval to the sub-range
// [cumFreq, cumFreq+freq) of totFreq, renormalising as needed.
func (e *rangeEncoder) encode(cumFreq, freq, totFreq uint32) {
	tmp := e.low
	e.rng /= totFreq
	e.low += cumFreq * e.rng
	e.rng *= freq
	if e.low < tmp {
		e.carry++
	}
	for e.rng < rcTop {
		e.rng <<= 8
		e.shiftLow()
	}
}

// finish ports RC_FinishEncode: five shiftLow calls flush the remaining
// state.
func (e *rangeEncoder) finish() []byte {
	for i := 0; i < 5; i++ {
		e.shiftLow()
	}
	return e.out
}

// rangeDecoder is the byte-wise range decoder from c_range_coder.h.
type rangeDecoder struct {
	code uint32
	rng  uint32
	in   []byte
	pos  int
	err  bool
}

// newRangeDecoder ports RC_StartDecode: it primes the decoder with the
// first five stream bytes (the first of which the encoder always emits as
// 0). A stream shorter than five bytes leaves the decoder flagged so the
// first symbol decode reports a corrupt stream.
func newRangeDecoder(in []byte) *rangeDecoder {
	d := &rangeDecoder{rng: 0xFFFFFFFF, in: in}
	if len(in) < 5 {
		d.err = true
		d.pos = len(in)
		return d
	}
	for i := 0; i < 5; i++ {
		d.code = (d.code << 8) | uint32(d.in[d.pos])
		d.pos++
	}
	return d
}

// getFreq ports RC_GetFreq: it returns the cumulative frequency the decoder
// is currently positioned at, dividing the range by totFreq as a side
// effect (exactly as the C does).
func (d *rangeDecoder) getFreq(totFreq uint32) uint32 {
	if totFreq != 0 && d.rng >= totFreq {
		d.rng /= totFreq
		return d.code / d.rng
	}
	return 0
}

// decode ports RC_Decode: it consumes the sub-range [cumFreq, cumFreq+freq)
// and renormalises, refilling code one byte at a time.
func (d *rangeDecoder) decode(cumFreq, freq uint32) {
	d.code -= cumFreq * d.rng
	d.rng *= freq
	for d.rng < rcTop {
		if d.pos >= len(d.in) {
			d.err = true
			return
		}
		d.code = (d.code << 8) + uint32(d.in[d.pos])
		d.pos++
		d.rng <<= 8
	}
}

// --- adaptive frequency model ------------------------------------------------

const (
	// arithMaxFreq is MAX_FREQ from c_simple_model.h: the total frequency
	// is renormalised once it would exceed this.
	arithMaxFreq = (1 << 16) - 17
	// arithStep is STEP: the increment applied to a symbol's frequency
	// each time it is coded.
	arithStep = 16
)

// symFreq is a (symbol, frequency) pair, mirroring SymFreqs.
type symFreq struct {
	freq uint16
	sym  uint16
}

// arithModel is the adaptive frequency model from c_simple_model.h. The
// symbol list F is kept approximately sorted by frequency: each coded
// symbol bumps its frequency and, if it now outranks its predecessor,
// swaps one place up the list (a single bubble-sort step). nsym is the
// alphabet size (256 for byte models, 258 for the RLE run models).
type arithModel struct {
	totFreq uint32
	nsym    int
	// f has nsym+1 entries; f[nsym] is the terminator with freq 0.
	f []symFreq
}

// newArithModel ports SIMPLE_MODEL(NSYM,_init): symbols 0..maxSym-1 start
// with frequency 1, the rest with 0.
func newArithModel(nsym, maxSym int) *arithModel {
	m := &arithModel{nsym: nsym, f: make([]symFreq, nsym+1)}
	i := 0
	for ; i < maxSym && i < nsym; i++ {
		m.f[i] = symFreq{sym: uint16(i), freq: 1}
	}
	for ; i < nsym; i++ {
		m.f[i] = symFreq{sym: uint16(i), freq: 0}
	}
	m.totFreq = uint32(maxSym)
	m.f[nsym] = symFreq{sym: 0, freq: 0} // terminates the normalize loop.
	return m
}

// normalize ports SIMPLE_MODEL(NSYM,_normalize): every non-zero frequency
// is halved (rounding up) and the total recomputed.
func (m *arithModel) normalize() {
	m.totFreq = 0
	for i := 0; i < len(m.f) && m.f[i].freq != 0; i++ {
		m.f[i].freq -= m.f[i].freq >> 1
		m.totFreq += uint32(m.f[i].freq)
	}
}

// encodeSymbol ports SIMPLE_MODEL(NSYM,_encodeSymbol): it range-encodes
// sym, updates the model and keeps the list approximately sorted.
func (m *arithModel) encodeSymbol(rc *rangeEncoder, sym uint16) {
	var acc uint32
	i := 0
	for m.f[i].sym != sym {
		acc += uint32(m.f[i].freq)
		i++
	}
	rc.encode(acc, uint32(m.f[i].freq), m.totFreq)
	m.f[i].freq += arithStep
	m.totFreq += arithStep
	if m.totFreq > arithMaxFreq {
		m.normalize()
	}
	if i > 0 && m.f[i].freq > m.f[i-1].freq {
		m.f[i], m.f[i-1] = m.f[i-1], m.f[i]
	}
}

// decodeSymbol ports SIMPLE_MODEL(NSYM,_decodeSymbol): it range-decodes one
// symbol, updates the model and keeps the list approximately sorted. ok is
// false on a corrupt stream (a frequency outside the model).
func (m *arithModel) decodeSymbol(rc *rangeDecoder) (sym uint16, ok bool) {
	freq := rc.getFreq(m.totFreq)
	if freq > arithMaxFreq {
		return 0, false
	}
	var acc uint32
	i := 0
	for {
		acc += uint32(m.f[i].freq)
		if acc > freq {
			break
		}
		i++
		if i > m.nsym {
			return 0, false
		}
	}
	acc -= uint32(m.f[i].freq)
	rc.decode(acc, uint32(m.f[i].freq))
	m.f[i].freq += arithStep
	m.totFreq += arithStep
	if m.totFreq > arithMaxFreq {
		m.normalize()
	}
	if i > 0 && m.f[i].freq > m.f[i-1].freq {
		t := m.f[i]
		m.f[i], m.f[i-1] = m.f[i-1], m.f[i]
		return t.sym, true
	}
	return m.f[i].sym, true
}

// --- order-0 / order-1 cores -------------------------------------------------

// arithCompressO0 ports arith_compress_O0: one adaptive byte model.
func arithCompressO0(in []byte) []byte {
	m := byte(0)
	for _, b := range in {
		if b > m {
			m = b
		}
	}
	mv := int(m) + 1
	out := []byte{byte(mv)}

	model := newArithModel(256, mv)
	rc := newRangeEncoder()
	for _, b := range in {
		model.encodeSymbol(rc, uint16(b))
	}
	return append(out, rc.finish()...)
}

// arithUncompressO0 ports arith_uncompress_O0.
func arithUncompressO0(in []byte, outSize uint32) ([]byte, error) {
	if len(in) < 1 {
		return nil, fmt.Errorf("arith: order-0 stream missing max-symbol byte")
	}
	mv := int(in[0])
	if mv == 0 {
		mv = 256
	}
	model := newArithModel(256, mv)
	rc := newRangeDecoder(in[1:])
	out := make([]byte, outSize)
	for i := range out {
		sym, ok := model.decodeSymbol(rc)
		if !ok || rc.err {
			return nil, fmt.Errorf("arith: corrupt order-0 range-coded stream at byte %d", i)
		}
		out[i] = byte(sym)
	}
	return out, nil
}

// arithCompressO1 ports arith_compress_O1: 256 adaptive byte models, one
// per preceding byte.
func arithCompressO1(in []byte) []byte {
	m := byte(0)
	for _, b := range in {
		if b > m {
			m = b
		}
	}
	mv := int(m) + 1
	out := []byte{byte(mv)}

	models := make([]*arithModel, 256)
	for i := range models {
		models[i] = newArithModel(256, mv)
	}
	rc := newRangeEncoder()
	last := byte(0)
	for _, b := range in {
		models[last].encodeSymbol(rc, uint16(b))
		last = b
	}
	return append(out, rc.finish()...)
}

// arithUncompressO1 ports arith_uncompress_O1.
func arithUncompressO1(in []byte, outSize uint32) ([]byte, error) {
	if len(in) < 1 {
		return nil, fmt.Errorf("arith: order-1 stream missing max-symbol byte")
	}
	mv := int(in[0])
	if mv == 0 {
		mv = 256
	}
	models := make([]*arithModel, 256)
	for i := range models {
		models[i] = newArithModel(256, mv)
	}
	rc := newRangeDecoder(in[1:])
	out := make([]byte, outSize)
	last := byte(0)
	for i := range out {
		sym, ok := models[last].decodeSymbol(rc)
		if !ok || rc.err {
			return nil, fmt.Errorf("arith: corrupt order-1 range-coded stream at byte %d", i)
		}
		out[i] = byte(sym)
		last = out[i]
	}
	return out, nil
}

// --- order-0 / order-1 RLE cores ---------------------------------------------

// arithRLENSym is NSYM for the run-length models (258 in arith_dynamic.c).
const arithRLENSym = 258

// arithRLEMaxRun is MAX_RUN: runs are coded in chunks of at most
// MAX_RUN-1, with a 0 terminator after a maximal chunk.
const arithRLEMaxRun = 4

// arithCompressO0RLE ports arith_compress_O0_RLE.
func arithCompressO0RLE(in []byte) []byte {
	m := byte(0)
	for _, b := range in {
		if b > m {
			m = b
		}
	}
	mv := int(m) + 1
	out := []byte{byte(mv)}

	byteModel := newArithModel(256, mv)
	runModel := make([]*arithModel, arithRLENSym)
	for i := range runModel {
		runModel[i] = newArithModel(arithRLENSym, arithRLEMaxRun)
	}
	rc := newRangeEncoder()
	for i := 0; i < len(in); {
		byteModel.encodeSymbol(rc, uint16(in[i]))
		run := 0
		last := in[i]
		i++
		for i < len(in) && in[i] == last {
			run++
			i++
		}
		rctx := int(last)
		for {
			c := run
			if c >= arithRLEMaxRun {
				c = arithRLEMaxRun - 1
			}
			runModel[rctx].encodeSymbol(rc, uint16(c))
			run -= c
			if rctx == int(last) {
				rctx = 256
			} else if rctx < arithRLENSym-1 {
				rctx++
			}
			if c == arithRLEMaxRun-1 && run == 0 {
				runModel[rctx].encodeSymbol(rc, 0)
			}
			if run == 0 {
				break
			}
		}
	}
	return append(out, rc.finish()...)
}

// arithUncompressO0RLE ports arith_uncompress_O0_RLE.
func arithUncompressO0RLE(in []byte, outSize uint32) ([]byte, error) {
	if len(in) < 1 {
		return nil, fmt.Errorf("arith: order-0 RLE stream missing max-symbol byte")
	}
	mv := int(in[0])
	if mv == 0 {
		mv = 256
	}
	byteModel := newArithModel(256, mv)
	runModel := make([]*arithModel, arithRLENSym)
	for i := range runModel {
		runModel[i] = newArithModel(arithRLENSym, arithRLEMaxRun)
	}
	rc := newRangeDecoder(in[1:])
	out := make([]byte, outSize)
	for i := 0; i < len(out); i++ {
		sym, ok := byteModel.decodeSymbol(rc)
		if !ok || rc.err {
			return nil, fmt.Errorf("arith: corrupt order-0 RLE stream at byte %d", i)
		}
		last := byte(sym)
		out[i] = last
		run, r, rctx := 0, 0, int(out[i])
		for {
			rs, ok := runModel[rctx].decodeSymbol(rc)
			if !ok || rc.err {
				return nil, fmt.Errorf("arith: corrupt order-0 RLE run at byte %d", i)
			}
			r = int(rs)
			if rctx == int(last) {
				rctx = 256
			} else if rctx < arithRLENSym-1 {
				rctx++
			}
			run += r
			if !(r == arithRLEMaxRun-1 && run < len(out)) {
				break
			}
		}
		for run > 0 && i+1 < len(out) {
			i++
			out[i] = last
			run--
		}
	}
	return out, nil
}

// arithCompressO1RLE ports arith_compress_O1_RLE.
func arithCompressO1RLE(in []byte) []byte {
	m := byte(0)
	for _, b := range in {
		if b > m {
			m = b
		}
	}
	mv := int(m) + 1
	out := []byte{byte(mv)}

	byteModel := make([]*arithModel, 256)
	for i := range byteModel {
		byteModel[i] = newArithModel(256, mv)
	}
	runModel := make([]*arithModel, arithRLENSym)
	for i := range runModel {
		runModel[i] = newArithModel(arithRLENSym, arithRLEMaxRun)
	}
	rc := newRangeEncoder()
	last := byte(0)
	for i := 0; i < len(in); {
		byteModel[last].encodeSymbol(rc, uint16(in[i]))
		run := 0
		last = in[i]
		i++
		for i < len(in) && in[i] == last {
			run++
			i++
		}
		rctx := int(last)
		for {
			c := run
			if c >= arithRLEMaxRun {
				c = arithRLEMaxRun - 1
			}
			runModel[rctx].encodeSymbol(rc, uint16(c))
			run -= c
			if rctx == int(last) {
				rctx = 256
			} else if rctx < arithRLENSym-1 {
				rctx++
			}
			if c == arithRLEMaxRun-1 && run == 0 {
				runModel[rctx].encodeSymbol(rc, 0)
			}
			if run == 0 {
				break
			}
		}
	}
	return append(out, rc.finish()...)
}

// arithUncompressO1RLE ports arith_uncompress_O1_RLE.
func arithUncompressO1RLE(in []byte, outSize uint32) ([]byte, error) {
	if len(in) < 1 {
		return nil, fmt.Errorf("arith: order-1 RLE stream missing max-symbol byte")
	}
	mv := int(in[0])
	if mv == 0 {
		mv = 256
	}
	byteModel := make([]*arithModel, 256)
	for i := range byteModel {
		byteModel[i] = newArithModel(256, mv)
	}
	runModel := make([]*arithModel, arithRLENSym)
	for i := range runModel {
		runModel[i] = newArithModel(arithRLENSym, arithRLEMaxRun)
	}
	rc := newRangeDecoder(in[1:])
	out := make([]byte, outSize)
	last := byte(0)
	for i := 0; i < len(out); i++ {
		sym, ok := byteModel[last].decodeSymbol(rc)
		if !ok || rc.err {
			return nil, fmt.Errorf("arith: corrupt order-1 RLE stream at byte %d", i)
		}
		out[i] = byte(sym)
		last = out[i]
		run, r, rctx := 0, 0, int(last)
		for {
			rs, ok := runModel[rctx].decodeSymbol(rc)
			if !ok || rc.err {
				return nil, fmt.Errorf("arith: corrupt order-1 RLE run at byte %d", i)
			}
			r = int(rs)
			if rctx == int(last) {
				rctx = 256
			} else if rctx < arithRLENSym-1 {
				rctx++
			}
			run += r
			if !(r == arithRLEMaxRun-1 && run < len(out)) {
				break
			}
		}
		for run > 0 && i+1 < len(out) {
			i++
			out[i] = last
			run--
		}
	}
	return out, nil
}
