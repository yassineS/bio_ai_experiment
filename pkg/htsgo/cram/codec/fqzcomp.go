package codec

import "fmt"

// fqzcomp — the htscodecs quality-score codec, CRAM block compression
// method 7.
//
// Ported from samtools/htscodecs fqzcomp_qual.c (fqz_compress /
// fqz_decompress and the parameter/model machinery). fqzcomp is a
// per-read quality-score compressor: a parameterised context model —
// built from the quality value, the position in the read, a
// running quality-change delta, the read selector, etc. — drives the
// adaptive range coder. The on-wire format is byte-identical to
// htscodecs, so the htscodecs test corpus
// (reference_code/htscodecs/tests/dat/fqzcomp/) is the compliance
// oracle; see fqzcomp_test.go.
//
// The range coder (rangeEncoder / rangeDecoder) and the adaptive
// frequency model (arithModel) are reused verbatim from arith.go — the
// fqzcomp "SIMPLE_MODEL" is the same c_simple_model.h type the
// arith_dynamic codec already ports. fqzcomp only adds the
// context-selection layer (fqzParam / fqzGParams) on top.
//
// The fqzcomp stream is self-describing: its header stores the
// uncompressed size, the parameter block(s) and — via the per-read
// length and selector models — every read boundary the decoder needs.
// A CRAM method-7 block therefore decodes standalone, without any
// external slice metadata.
//
// Stream layout:
//
//	[uncompressedSize:varint][parameter block(s)][range-coder bytes]

// fqzVers is FQZ_VERS: the only fqzcomp format version this codec
// accepts.
const fqzVers = 5

// fqzCtxBits / fqzCtxSize mirror CTX_BITS / CTX_SIZE: the context space
// addressed by the per-context quality models.
const (
	fqzCtxBits = 16
	fqzCtxSize = 1 << fqzCtxBits
)

// fqzMaxSym is QMAX: the alphabet size of the quality models.
const fqzMaxSym = 256

// fqzMaxDecodeRatio / fqzMinDecodeSlack bound the decoded size a
// declared header may claim relative to the compressed input, so a
// corrupt or hostile stream cannot trigger an unbounded allocation. The
// ratio is far above any realistic fqzcomp compression factor on
// quality data; the slack covers tiny inputs whose header overhead
// dominates.
const (
	fqzMaxDecodeRatio = 256
	fqzMinDecodeSlack = 1 << 16
)

// Global parameter flags (GFLAG_*).
const (
	fqzGFlagMultiParam = 1
	fqzGFlagHaveStab   = 2
	fqzGFlagDoRev      = 4
)

// Per-parameter-block flags (PFLAG_*).
const (
	fqzPFlagDoDedup  = 2
	fqzPFlagDoLen    = 4
	fqzPFlagDoSel    = 8
	fqzPFlagHaveQMap = 16
	fqzPFlagHavePTab = 32
	fqzPFlagHaveDTab = 64
	fqzPFlagHaveQTab = 128
)

// Record flag bits, mirroring the BAM flags fqzcomp reuses.
const (
	fqzFReverse = 16
	fqzFRead2   = 128
)

// fqzParam is one fqzcomp parameter block — the context-model
// configuration for a group of reads (fqz_param in fqzcomp_qual.c).
type fqzParam struct {
	context uint16 // starting context value

	pflags                                      int
	doSel, doDedup, storeQMap, fixedLen         bool
	useQTab, useDTab, usePTab                   bool
	qbits, qloc, pbits, ploc, dbits, dloc, sloc int
	pshift, qshift, dshift                      int
	maxSym, nsym, maxSel                        int
	qmask                                       int
	doR2, doQA                                  int

	qmap [256]int
	qtab [256]uint32
	ptab [1024]uint32
	dtab [256]uint32
}

// fqzGParams is the global fqzcomp parameter set: a collection of
// parameter blocks plus stream-wide metadata (fqz_gparams).
type fqzGParams struct {
	vers   int
	gflags int
	nparam int
	maxSel int
	maxSym int
	stab   [256]uint32
	p      []fqzParam
}

// fqzState is the per-read decoding cursor (fqz_state).
type fqzState struct {
	qctx     uint32 // quality sub-context
	p        int    // position: bytes remaining in the read
	delta    uint32 // running quality-change total
	prevq    uint32 // previous quality value
	s        int    // selector
	firstLen bool
	lastLen  int
	rec      int
	ctx      uint32
}

// fqzModel bundles every adaptive model the codec uses (fqz_model). The
// per-context quality models dominate: CTX_SIZE of them.
type fqzModel struct {
	qual    []*arithModel
	length  [4]*arithModel
	revcomp *arithModel
	sel     *arithModel
	dup     *arithModel
}

// newFQZModel builds the model set for the given global parameters,
// mirroring fqz_create_models.
func newFQZModel(gp *fqzGParams) *fqzModel {
	m := &fqzModel{}
	// The CTX_SIZE quality models share one contiguous backing array,
	// mirroring the single htscodecs_tls_alloc the C codec uses.
	m.qual = newArithModelBatch(fqzCtxSize, fqzMaxSym, gp.maxSym+1)
	for i := range m.length {
		m.length[i] = newArithModel(256, 256)
	}
	m.revcomp = newArithModel(2, 2)
	m.dup = newArithModel(2, 2)
	if gp.maxSel > 0 {
		m.sel = newArithModel(256, gp.maxSel+1)
	}
	return m
}

// --- varint -----------------------------------------------------------------

// fqzGetU32 reads one big-endian 7-bit varint (var_get_u32). It returns
// the value and the number of bytes consumed, or ok=false if the buffer
// is exhausted mid-value.
func fqzGetU32(in []byte) (val uint32, n int, ok bool) {
	var j uint32
	for n < 5 && n < len(in) {
		c := in[n]
		j = (j << 7) | uint32(c&0x7f)
		n++
		if c&0x80 == 0 {
			return j, n, true
		}
	}
	if n == 0 {
		return 0, 0, false
	}
	// Either ran out of buffer or hit the 5-byte cap; treat a clean cap
	// (last byte had no continuation bit) as success.
	if n <= len(in) && in[n-1]&0x80 == 0 {
		return j, n, true
	}
	return 0, n, false
}

// fqzPutU32 appends val to out as a big-endian 7-bit varint
// (var_put_u32).
func fqzPutU32(out []byte, val uint32) []byte {
	switch {
	case val < 1<<7:
		return append(out, byte(val))
	case val < 1<<14:
		return append(out, byte((val>>7)&0x7f|128), byte(val&0x7f))
	case val < 1<<21:
		return append(out, byte((val>>14)&0x7f|128), byte((val>>7)&0x7f|128),
			byte(val&0x7f))
	case val < 1<<28:
		return append(out, byte((val>>21)&0x7f|128), byte((val>>14)&0x7f|128),
			byte((val>>7)&0x7f|128), byte(val&0x7f))
	default:
		return append(out, byte((val>>28)&0x7f|128), byte((val>>21)&0x7f|128),
			byte((val>>14)&0x7f|128), byte((val>>7)&0x7f|128), byte(val&0x7f))
	}
}

// --- run-length tables -------------------------------------------------------

// fqzReadArray decodes a doubly run-length-encoded table (read_array).
// It fills array[0:size] and returns the number of input bytes consumed,
// or ok=false on a malformed stream.
func fqzReadArray(in []byte, array []uint32, size int) (used int, ok bool) {
	if size > 1024 {
		size = 1024
	}
	var r [1024]int
	i, j, z := 0, 0, 0
	last := -1
	// Remove level one of run-length encoding.
	for z < size && i < len(in) {
		run := int(in[i])
		r[j] = run
		j++
		z += run
		if run == last {
			i++
			if i >= len(in) {
				return 0, false
			}
			copies := int(in[i])
			z += run * copies
			for copies > 0 && z <= size && j < 1024 {
				r[j] = run
				j++
				copies--
			}
		}
		if j >= 1024 {
			return 0, false
		}
		last = run
		i++
	}
	nb := i

	// Expand the inner level of run-length encoding.
	rMax := j
	i, j, z = 0, 0, 0
	for j < size {
		runLen := 0
		var runPart int
		if z >= rMax {
			return 0, false
		}
		for {
			runPart = r[z]
			z++
			runLen += runPart
			if runPart != 255 || z >= rMax {
				break
			}
		}
		if runPart == 255 {
			return 0, false
		}
		for runLen > 0 && j < size {
			runLen--
			array[j] = uint32(i)
			j++
		}
		i++
	}
	return nb, true
}

// fqzStoreArray encodes array[0:size] with the same doubly run-length
// scheme (store_array) and appends it to out.
func fqzStoreArray(out []byte, array []uint32, size int) []byte {
	var tmp [2048]byte
	i, j, k := 0, 0, 0
	for i < size {
		runStart := i
		for i < size && int(array[i]) == j {
			i++
		}
		runLen := i - runStart
		for {
			r := runLen
			if r > 255 {
				r = 255
			}
			tmp[k] = byte(r)
			k++
			runLen -= r
			if r != 255 {
				break
			}
		}
		j++
	}
	for i < size {
		tmp[k] = 0
		k++
		j++
	}

	// RLE on the output.
	enc := make([]byte, 0, k+8)
	last := -1
	for j = 0; j < k; {
		b := tmp[j]
		j++
		enc = append(enc, b)
		if int(b) == last {
			n := j
			for j < k && tmp[j] == byte(last) {
				j++
			}
			enc = append(enc, byte(j-n))
		} else {
			last = int(b)
		}
	}
	return append(out, enc...)
}

// --- parameter (de)serialisation --------------------------------------------

// fqzReadParameters1 parses one parameter block (fqz_read_parameters1).
func fqzReadParameters1(pm *fqzParam, in []byte) (used int, ok bool) {
	if len(in) < 7 {
		return 0, false
	}
	idx := 0
	pm.context = uint16(in[idx]) | uint16(in[idx+1])<<8
	idx += 2

	pm.pflags = int(in[idx])
	idx++
	pm.useQTab = pm.pflags&fqzPFlagHaveQTab != 0
	pm.useDTab = pm.pflags&fqzPFlagHaveDTab != 0
	pm.usePTab = pm.pflags&fqzPFlagHavePTab != 0
	pm.doSel = pm.pflags&fqzPFlagDoSel != 0
	pm.fixedLen = pm.pflags&fqzPFlagDoLen != 0
	pm.doDedup = pm.pflags&fqzPFlagDoDedup != 0
	pm.storeQMap = pm.pflags&fqzPFlagHaveQMap != 0
	pm.maxSym = int(in[idx])
	idx++

	pm.qbits = int(in[idx]) >> 4
	pm.qmask = (1 << uint(pm.qbits)) - 1
	pm.qshift = int(in[idx]) & 15
	idx++
	pm.qloc = int(in[idx]) >> 4
	pm.sloc = int(in[idx]) & 15
	idx++
	pm.ploc = int(in[idx]) >> 4
	pm.dloc = int(in[idx]) & 15
	idx++

	if pm.storeQMap {
		for i := range pm.qmap {
			pm.qmap[i] = intMax
		}
		if idx+pm.maxSym > len(in) {
			return 0, false
		}
		for i := 0; i < pm.maxSym; i++ {
			pm.qmap[i] = int(in[idx])
			idx++
		}
	} else {
		for i := range pm.qmap {
			pm.qmap[i] = i
		}
	}

	if pm.qbits != 0 {
		if pm.useQTab {
			n, good := fqzReadArray(in[idx:], pm.qtab[:], 256)
			if !good {
				return 0, false
			}
			idx += n
		} else {
			for i := range pm.qtab {
				pm.qtab[i] = uint32(i)
			}
		}
	}

	if pm.usePTab {
		n, good := fqzReadArray(in[idx:], pm.ptab[:], 1024)
		if !good {
			return 0, false
		}
		idx += n
	} else {
		for i := range pm.ptab {
			pm.ptab[i] = 0
		}
	}

	if pm.useDTab {
		n, good := fqzReadArray(in[idx:], pm.dtab[:], 256)
		if !good {
			return 0, false
		}
		idx += n
	} else {
		for i := range pm.dtab {
			pm.dtab[i] = 0
		}
	}

	return idx, true
}

// intMax is INT_MAX, used as the "no mapping" sentinel in qmap.
const intMax = int(^uint32(0) >> 1)

// fqzReadParameters parses the global parameter set (fqz_read_parameters)
// and returns the number of header bytes consumed.
func fqzReadParameters(gp *fqzGParams, in []byte) (used int, ok bool) {
	if len(in) < 10 {
		return 0, false
	}
	idx := 0
	gp.vers = int(in[idx])
	idx++
	if gp.vers != fqzVers {
		return 0, false
	}
	gp.gflags = int(in[idx])
	idx++

	if gp.gflags&fqzGFlagMultiParam != 0 {
		gp.nparam = int(in[idx])
		idx++
	} else {
		gp.nparam = 1
	}
	if gp.nparam <= 0 {
		return 0, false
	}
	if gp.nparam > 1 {
		gp.maxSel = gp.nparam
	} else {
		gp.maxSel = 0
	}

	if gp.gflags&fqzGFlagHaveStab != 0 {
		if idx >= len(in) {
			return 0, false
		}
		gp.maxSel = int(in[idx])
		idx++
		n, good := fqzReadArray(in[idx:], gp.stab[:], 256)
		if !good {
			return 0, false
		}
		idx += n
	} else {
		for i := 0; i < gp.nparam && i < 256; i++ {
			gp.stab[i] = uint32(i)
		}
		for i := gp.nparam; i < 256; i++ {
			gp.stab[i] = uint32(gp.nparam - 1)
		}
	}

	gp.p = make([]fqzParam, gp.nparam)
	gp.maxSym = 0
	for i := 0; i < gp.nparam; i++ {
		n, good := fqzReadParameters1(&gp.p[i], in[idx:])
		if !good {
			return 0, false
		}
		if gp.p[i].doSel && gp.maxSel == 0 {
			return 0, false // inconsistent
		}
		idx += n
		if gp.maxSym < gp.p[i].maxSym {
			gp.maxSym = gp.p[i].maxSym
		}
	}
	return idx, true
}

// fqzStoreParameters1 serialises one parameter block
// (fqz_store_parameters1).
func fqzStoreParameters1(out []byte, pm *fqzParam) []byte {
	out = append(out,
		byte(pm.context), byte(pm.context>>8),
		byte(pm.pflags), byte(pm.maxSym),
		byte(pm.qbits<<4|pm.qshift),
		byte(pm.qloc<<4|pm.sloc),
		byte(pm.ploc<<4|pm.dloc))

	if pm.storeQMap {
		for i := 0; i < 256; i++ {
			if pm.qmap[i] != intMax {
				out = append(out, byte(i))
			}
		}
	}
	if pm.qbits != 0 && pm.useQTab {
		out = fqzStoreArray(out, pm.qtab[:], 256)
	}
	if pm.pbits != 0 && pm.usePTab {
		out = fqzStoreArray(out, pm.ptab[:], 1024)
	}
	if pm.dbits != 0 && pm.useDTab {
		out = fqzStoreArray(out, pm.dtab[:], 256)
	}
	return out
}

// fqzStoreParameters serialises the global parameter set
// (fqz_store_parameters).
func fqzStoreParameters(out []byte, gp *fqzGParams) []byte {
	out = append(out, byte(gp.vers), byte(gp.gflags))
	if gp.gflags&fqzGFlagMultiParam != 0 {
		out = append(out, byte(gp.nparam))
	}
	if gp.gflags&fqzGFlagHaveStab != 0 {
		out = append(out, byte(gp.maxSel))
		out = fqzStoreArray(out, gp.stab[:], 256)
	}
	for i := range gp.p {
		out = fqzStoreParameters1(out, &gp.p[i])
	}
	return out
}

// --- context update ----------------------------------------------------------

// fqzUpdateCtx ports fqz_update_ctx: it folds the just-coded quality q
// into the context state and returns the next quality-model context.
// ptab and dtab must already be pre-shifted by ploc / dloc (as the C
// does once up front).
func fqzUpdateCtx(pm *fqzParam, st *fqzState, q int) uint32 {
	var last uint32
	st.qctx = (st.qctx << uint(pm.qshift)) + pm.qtab[q]
	last += (st.qctx & uint32(pm.qmask)) << uint(pm.qloc)

	pIdx := st.p
	if pIdx > 1023 {
		pIdx = 1023
	}
	last += pm.ptab[pIdx]
	dIdx := st.delta
	if dIdx > 255 {
		dIdx = 255
	}
	last += pm.dtab[dIdx]
	last += uint32(st.s) << uint(pm.sloc)

	if st.prevq != uint32(q) {
		st.delta++
	}
	st.prevq = uint32(q)
	st.p--

	return last & (fqzCtxSize - 1)
}

// --- decode ------------------------------------------------------------------

// FQZCompDecode decodes an fqzcomp-compressed quality-score block (CRAM
// block compression method 7) back to the raw concatenated quality
// buffer. The fqzcomp stream is self-describing, so no external slice
// metadata is required.
func FQZCompDecode(in []byte) ([]byte, error) {
	out, _, err := fqzDecode(in)
	return out, err
}

// fqzDecode is the workhorse behind FQZCompDecode; it also returns the
// per-read length array (the lengths fqz_decompress fills out), which
// the compliance tests use to reconstruct the newline-delimited input.
func fqzDecode(in []byte) ([]byte, []int, error) {
	size, n, ok := fqzGetU32(in)
	if !ok {
		return nil, nil, fmt.Errorf("fqzcomp: truncated uncompressed-size header")
	}
	idx := n

	// Guard against a hostile uncompressed-size header forcing a huge
	// allocation. Every decoded byte must consume range-coder progress,
	// so the decoded length cannot meaningfully exceed a generous
	// multiple of the compressed input; fqzcomp's real-world ratio on
	// quality data is well under 100x. CRAM's own UncompressedSize check
	// in block.go re-validates the exact size afterwards.
	if uint64(size) > fqzMaxDecodeRatio*uint64(len(in))+fqzMinDecodeSlack {
		return nil, nil, fmt.Errorf("fqzcomp: declared uncompressed size %d implausible for %d input bytes",
			size, len(in))
	}

	var gp fqzGParams
	used, ok := fqzReadParameters(&gp, in[idx:])
	if !ok {
		return nil, nil, fmt.Errorf("fqzcomp: malformed parameter block")
	}
	idx += used

	// Pre-shift ptab / dtab so the hot loop needs no shift, exactly as
	// uncompress_block_fqz2f does.
	for i := range gp.p {
		pm := &gp.p[i]
		for j := range pm.ptab {
			pm.ptab[j] <<= uint(pm.ploc)
		}
		for j := range pm.dtab {
			pm.dtab[j] <<= uint(pm.dloc)
		}
	}

	model := newFQZModel(&gp)
	rc := newRangeDecoder(in[idx:])

	out := make([]byte, size)
	var lengths []int

	st := fqzState{firstLen: true}
	// uncompress_block_fqz2f keeps pm pointed at gp.p[0] for the whole
	// quality-decode loop; decompress_new_read selects its own block
	// internally for length/dedup, but the context model uses gp.p[0].
	pm := &gp.p[0]
	rev := 0
	type revRec struct{ rev, length int }
	var revRecs []revRec

	i := 0
	for i < int(size) {
		if st.p == 0 {
			r, err := fqzDecodeNewRead(&gp, &st, model, rc, out, int(size), &i,
				&rev, &lengths)
			if err != nil {
				return nil, nil, err
			}
			if gp.gflags&fqzGFlagDoRev != 0 {
				// st.lastLen is this read's length: fqzDecodeNewRead has
				// just set it, either from the length model or the
				// fixed-length reuse. It is the C decoder's per-read
				// len_a[rec].
				revRecs = append(revRecs, revRec{rev, st.lastLen})
			}
			if r == 1 {
				continue
			}
		}

		for {
			q, ok := model.qual[st.ctx].decodeSymbol(rc)
			if !ok || rc.err {
				return nil, nil, fmt.Errorf("fqzcomp: corrupt range-coded stream at byte %d", i)
			}
			st.ctx = fqzUpdateCtx(pm, &st, int(q))
			out[i] = byte(pm.qmap[q])
			i++
			if st.p == 0 || i >= int(size) {
				break
			}
		}
	}

	if gp.gflags&fqzGFlagDoRev != 0 {
		pos := 0
		for _, rr := range revRecs {
			if rr.rev != 0 {
				lo, hi := pos, pos+rr.length-1
				for lo < hi {
					out[lo], out[hi] = out[hi], out[lo]
					lo, hi = lo+1, hi-1
				}
			}
			pos += rr.length
			if pos >= int(size) {
				break
			}
		}
	}

	if rc.err {
		return nil, nil, fmt.Errorf("fqzcomp: corrupt range-coded stream")
	}
	return out, lengths, nil
}

// fqzDecodeNewRead ports decompress_new_read: it starts a new read,
// decoding the selector, length, reverse flag and dedup bit. It returns
// 1 if the read was a dedup copy (caller should continue), 0 otherwise.
//
// As in the C, the do_sel decision and the length/dedup models are
// driven by gp.p[0] (the caller's pm), while the per-read length block
// is selected via the stab/selector.
func fqzDecodeNewRead(gp *fqzGParams, st *fqzState, model *fqzModel,
	rc *rangeDecoder, uncomp []byte, outSize int, ip *int, rev *int,
	lengths *[]int) (int, error) {

	i := *ip
	if gp.p[0].doSel {
		s, ok := model.sel.decodeSymbol(rc)
		if !ok || rc.err {
			return 0, fmt.Errorf("fqzcomp: corrupt selector")
		}
		st.s = int(s)
	} else {
		st.s = 0
	}

	x := st.s
	if gp.gflags&fqzGFlagHaveStab != 0 {
		si := st.s
		if si > 255 {
			si = 255
		}
		x = int(gp.stab[si])
	}
	if x >= gp.nparam {
		return 0, fmt.Errorf("fqzcomp: selector %d out of range", x)
	}
	pm := &gp.p[x]

	length := st.lastLen
	if !pm.fixedLen || st.firstLen {
		l0, ok0 := model.length[0].decodeSymbol(rc)
		l1, ok1 := model.length[1].decodeSymbol(rc)
		l2, ok2 := model.length[2].decodeSymbol(rc)
		l3, ok3 := model.length[3].decodeSymbol(rc)
		if !ok0 || !ok1 || !ok2 || !ok3 || rc.err {
			return 0, fmt.Errorf("fqzcomp: corrupt read length")
		}
		length = int(l0) | int(l1)<<8 | int(l2)<<16 | int(uint32(l3)<<24)
		st.firstLen = false
		st.lastLen = length
	}
	if length > outSize-i || length <= 0 {
		return 0, fmt.Errorf("fqzcomp: invalid read length %d", length)
	}

	if lengths != nil {
		*lengths = append(*lengths, length)
	}

	if gp.gflags&fqzGFlagDoRev != 0 {
		r, ok := model.revcomp.decodeSymbol(rc)
		if !ok || rc.err {
			return 0, fmt.Errorf("fqzcomp: corrupt reverse flag")
		}
		*rev = int(r)
	}

	if pm.doDedup {
		d, ok := model.dup.decodeSymbol(rc)
		if !ok || rc.err {
			return 0, fmt.Errorf("fqzcomp: corrupt dedup flag")
		}
		if d != 0 {
			if length > i {
				return 0, fmt.Errorf("fqzcomp: dedup before any data")
			}
			copy(uncomp[i:i+length], uncomp[i-length:i])
			i += length
			st.p = 0
			st.rec++
			*ip = i
			return 1, nil
		}
	}

	st.rec++
	st.p = length
	st.delta = 0
	st.prevq = 0
	st.qctx = 0
	st.ctx = uint32(pm.context)

	*ip = i
	return 0, nil
}
