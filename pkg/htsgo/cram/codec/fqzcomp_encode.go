package codec

import (
	"fmt"
	"math"
)

// fqzcomp encoder — the compress side of CRAM block compression method
// 7, ported from samtools/htscodecs fqzcomp_qual.c (fqz_compress,
// compress_block_fqz2f, fqz_pick_parameters, fqz_qual_stats and the
// store_array / store_parameters machinery).
//
// fqzcomp's parameter selection is entirely deterministic: given the
// strategy level (the -s flag, 0..3), the input quality buffer and the
// per-read length array, fqz_pick_parameters chooses a fixed parameter
// block, with no randomness. Reproducing that selection bit-for-bit
// makes the encoder byte-exact against the htscodecs vectors.
//
// fqzSlice mirrors the C fqz_slice: the per-read metadata the encoder
// needs (length, and the read1/2 + selector flags). For the standalone
// compliance vectors flags is all-zero — the fqzcomp_qual test tool's
// `cut -f1` transform discards the read2/selector columns.

// fqzSlice is the per-read metadata an fqzcomp encode works from.
type fqzSlice struct {
	numRecords int
	length     []uint32
	flags      []uint32
}

// fqzStratOpts is strat_opts[][12]: the predefined strategy parameter
// templates. Index by strategy (0..3); index 4 is the "custom" slot.
var fqzStratOpts = [5][12]int{
	{10, 5, 4, -1, 2, 1, 0, 14, 10, 14, 0, -1}, // basic (level < 7)
	{8, 5, 7, 0, 0, 0, 0, 14, 8, 14, 1, -1},    // e.g. HiSeq 2000
	{12, 6, 2, 0, 2, 3, 0, 9, 12, 14, 0, 0},    // e.g. MiSeq
	{12, 6, 0, 0, 0, 0, 0, 12, 0, 0, 0, 0},     // e.g. IonTorrent
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},       // custom
}

// fqzNStrats is the number of strategy templates.
const fqzNStrats = 5

// fqzFastLog ports htscodecs' fast_log: a bit-trick approximation of the
// natural logarithm. fqz_qual_stats compares entropy estimates computed
// with this exact function, so a byte-exact encoder must reproduce it.
func fqzFastLog(a float64) float64 {
	x := int64(math.Float64bits(a))
	return float64(x-4606921278410026770) * 1.539095918623324e-16
}

// fqzQualStats ports fqz_qual_stats for one_param == -1: it scans the
// quality buffer, fills qhist, sets nsym / max_sym / do_dedup, and
// auto-tunes do_sel and the read1/2 split via entropy comparison.
func fqzQualStats(s *fqzSlice, in []byte, pm *fqzParam, qhist *[256]uint32) {
	const np = 32

	var qhistb, qhist1, qhist2 [np][256]uint32
	var t1, t2 [np]uint64

	avg := make([]uint32, 2560)
	dir := 0
	lastLen := 0
	doDedup := 0
	numRec := 0

	hasR2 := 0
	maxSel := 0
	for rec := 0; rec < s.numRecords; rec++ {
		numRec++
		if maxSel < int(s.flags[rec]>>16) {
			maxSel = int(s.flags[rec] >> 16)
		}
		if s.flags[rec]&fqzFRead2 != 0 {
			hasR2 = 1
		}
	}

	avgQual := make([]int, s.numRecords+1)

	rec := 0
	i, j := 0, 0
	for i < len(in) {
		if rec < s.numRecords {
			j = int(s.length[rec])
			if s.flags[rec]&fqzFRead2 != 0 {
				dir = 1
			} else {
				dir = 0
			}
			if i > 0 && j == lastLen &&
				bytesEqual(in[i-lastLen:i], in[i:i+j]) {
				doDedup++
			}
		} else {
			j = len(in) - i
			dir = 0
		}
		lastLen = j

		var tot uint32
		for i < len(in) && j > 0 {
			b := in[i]
			tot += uint32(b)
			qhist[b]++
			qhistb[j&(np-1)][b]++
			if dir != 0 {
				qhist2[j&(np-1)][b]++
				t2[j&(np-1)]++
			} else {
				qhist1[j&(np-1)][b]++
				t1[j&(np-1)]++
			}
			i++
			j--
		}
		if lastLen != 0 {
			tot = uint32(float64(tot)*10.0/float64(lastLen) + 0.5)
		} else {
			tot = 0
		}
		avgQual[rec] = int(tot)
		avg[minInt(2559, int(tot))]++
		rec++
	}
	if (doDedup+1) != 0 && (rec+1)/(doDedup+1) < 500 {
		pm.doDedup = true
	} else {
		pm.doDedup = false
	}

	lastLen = 0

	// Unique symbol count.
	pm.maxSym = 0
	pm.nsym = 0
	for k := 0; k < 256; k++ {
		if qhist[k] != 0 {
			pm.maxSym = k
			pm.nsym++
		}
	}

	// Auto-tune: does average quality help?
	if pm.doQA != 0 {
		var qf0, qf1, qf2 float64
		if pm.nsym > 8 {
			qf0, qf1, qf2 = 0.2, 0.5, 0.8
		} else {
			qf0, qf1, qf2 = 0.05, 0.22, 0.60
		}

		total := 0
		k := 0
		for k < 2560 {
			total += int(avg[k])
			if float64(total) > qf0*float64(numRec) {
				break
			}
			avg[k] = 0
			k++
		}
		for k < 2560 {
			total += int(avg[k])
			if float64(total) > qf1*float64(numRec) {
				break
			}
			avg[k] = 1
			k++
		}
		for k < 2560 {
			total += int(avg[k])
			if float64(total) > qf2*float64(numRec) {
				break
			}
			avg[k] = 2
			k++
		}
		for k < 2560 {
			avg[k] = 3
			k++
		}

		var qbin4 [4][np][256]int
		var qbin2 [2][np][256]int
		var qbin1 [np][256]int
		var qcnt4 [4][np]int
		var qcnt2 [4][np]int
		var qcnt1 [np]int

		i = 0
		rec = 0
		for i < len(in) {
			if (rec&7) != 0 && rec < s.numRecords {
				i += int(s.length[rec])
				rec++
				continue
			}
			if rec < s.numRecords {
				j = int(s.length[rec])
			} else {
				j = len(in) - i
			}
			lastLen = j

			tot := avgQual[rec]
			qb4 := int(avg[minInt(2559, tot)])
			qb2 := qb4 / 2

			for i < len(in) && j > 0 {
				x := j & (np - 1)
				b := in[i]
				qbin4[qb4][x][b]++
				qcnt4[qb4][x]++
				qbin2[qb2][x][b]++
				qcnt2[qb2][x]++
				qbin1[x][b]++
				qcnt1[x]++
				i++
				j--
			}
			rec++
		}

		var e1, e2, e4 float64
		for jj := 0; jj < np; jj++ {
			for ii := 0; ii < 256; ii++ {
				if qbin1[jj][ii] != 0 {
					e1 += float64(qbin1[jj][ii]) * fqzFastLog(float64(qbin1[jj][ii])/float64(qcnt1[jj]))
				}
				if qbin2[0][jj][ii] != 0 {
					e2 += float64(qbin2[0][jj][ii]) * fqzFastLog(float64(qbin2[0][jj][ii])/float64(qcnt2[0][jj]))
				}
				if qbin2[1][jj][ii] != 0 {
					e2 += float64(qbin2[1][jj][ii]) * fqzFastLog(float64(qbin2[1][jj][ii])/float64(qcnt2[1][jj]))
				}
				if qbin4[0][jj][ii] != 0 {
					e4 += float64(qbin4[0][jj][ii]) * fqzFastLog(float64(qbin4[0][jj][ii])/float64(qcnt4[0][jj]))
				}
				if qbin4[1][jj][ii] != 0 {
					e4 += float64(qbin4[1][jj][ii]) * fqzFastLog(float64(qbin4[1][jj][ii])/float64(qcnt4[1][jj]))
				}
				if qbin4[2][jj][ii] != 0 {
					e4 += float64(qbin4[2][jj][ii]) * fqzFastLog(float64(qbin4[2][jj][ii])/float64(qcnt4[2][jj]))
				}
				if qbin4[3][jj][ii] != 0 {
					e4 += float64(qbin4[3][jj][ii]) * fqzFastLog(float64(qbin4[3][jj][ii])/float64(qcnt4[3][jj]))
				}
			}
		}
		e1 /= -math.Log(2) / 8
		e2 /= -math.Log(2) / 8
		e4 /= -math.Log(2) / 8

		qm := 0.98
		if pm.doQA > 0 {
			qm = 1
		}
		if (pm.doQA == -1 || pm.doQA >= 4) &&
			e4+float64(s.numRecords)/4 < e2*qm+float64(s.numRecords)/8 &&
			e4+float64(s.numRecords)/4 < e1*qm {
			for k := 0; k < s.numRecords; k++ {
				s.flags[k] |= avg[minInt(2559, avgQual[k])] << 16
			}
			pm.doSel = true
			maxSel = 3
		} else if (pm.doQA == -1 || pm.doQA >= 2) && e2+float64(s.numRecords)/8 < e1*qm {
			for k := 0; k < s.numRecords; k++ {
				s.flags[k] |= (avg[minInt(2559, avgQual[k])] >> 1) << 16
			}
			pm.doSel = true
			maxSel = 1
		}

		if pm.doQA == -1 {
			if pm.pbits > 0 && pm.dbits > 0 {
				pm.sloc = pm.dloc - 1
				pm.pbits--
				pm.dbits--
				pm.dloc++
			} else if pm.dbits >= 2 {
				pm.sloc = pm.dloc
				pm.dbits -= 2
				pm.dloc += 2
			} else if pm.qbits >= 2 {
				pm.qbits -= 2
				pm.ploc -= 2
				pm.sloc = 16 - 2 - pm.doR2
				if pm.qbits == 6 && pm.qshift == 5 {
					pm.qbits--
				}
			}
			pm.doQA = 4
		}
	}

	// Auto-tune: does splitting READ1 / READ2 help?
	if hasR2 != 0 || pm.doR2 != 0 {
		var e1, e2 float64
		for jj := 0; jj < np; jj++ {
			if t1[jj] == 0 || t2[jj] == 0 {
				continue
			}
			for ii := 0; ii < 256; ii++ {
				if qhistb[jj][ii] == 0 {
					continue
				}
				e1 -= float64(qhistb[jj][ii]) * math.Log(float64(qhistb[jj][ii])/float64(t1[jj]+t2[jj]))
				if qhist1[jj][ii] != 0 {
					e2 -= float64(qhist1[jj][ii]) * math.Log(float64(qhist1[jj][ii])/float64(t1[jj]))
				}
				if qhist2[jj][ii] != 0 {
					e2 -= float64(qhist2[jj][ii]) * math.Log(float64(qhist2[jj][ii])/float64(t2[jj]))
				}
			}
		}
		e1 /= math.Log(2) * 8
		e2 /= math.Log(2) * 8

		qm := 0.95
		if pm.doR2 > 0 {
			qm = 1
		}
		if e2+float64(8+s.numRecords/8) < e1*qm {
			for r := 0; r < s.numRecords; r++ {
				sel := int(s.flags[r] >> 16)
				if s.flags[r]&fqzFRead2 != 0 {
					s.flags[r] = (s.flags[r] & 0xffff) | uint32((sel*2)+1)<<16
				} else {
					s.flags[r] = (s.flags[r] & 0xffff) | uint32((sel*2)+0)<<16
				}
				if maxSel < int(s.flags[r]>>16) {
					maxSel = int(s.flags[r] >> 16)
				}
			}
		}
	}

	if maxSel > 0 {
		pm.doSel = true
		pm.maxSel = maxSel
	}
}

// bytesEqual reports whether two byte slices have identical contents.
func bytesEqual(a, b []byte) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// minInt / maxInt are the obvious integer helpers (no generics: CI pins
// Go 1.21 but this keeps it explicit).
func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// fqzPickParameters ports fqz_pick_parameters: it builds the global
// parameter set for the given strategy and input. strat is the -s level
// (0..3); vers selects the CRAM generation (4 for v4.0, 3 for v3.1).
func fqzPickParameters(strat, vers int, s *fqzSlice, in []byte) (*fqzGParams, error) {
	dsqr := []uint32{
		0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	}
	var qhist [256]uint32

	if strat >= fqzNStrats {
		strat = fqzNStrats - 1
	}

	gp := &fqzGParams{}
	gp.vers = fqzVers
	gp.p = make([]fqzParam, 1)
	gp.nparam = 1
	gp.maxSel = 0

	if vers == 3 {
		gp.gflags |= fqzGFlagDoRev
	}

	pm := &gp.p[0]
	pm.qbits = fqzStratOpts[strat][0]
	pm.qshift = fqzStratOpts[strat][1]
	pm.pbits = fqzStratOpts[strat][2]
	pm.pshift = fqzStratOpts[strat][3]
	pm.dbits = fqzStratOpts[strat][4]
	pm.dshift = fqzStratOpts[strat][5]
	pm.qloc = fqzStratOpts[strat][6]
	pm.sloc = fqzStratOpts[strat][7]
	pm.ploc = fqzStratOpts[strat][8]
	pm.dloc = fqzStratOpts[strat][9]
	pm.doR2 = fqzStratOpts[strat][10]
	pm.doQA = fqzStratOpts[strat][11]

	// Validity-check input lengths and buffer size.
	tlen := 0
	for i := 0; i < s.numRecords; i++ {
		if tlen+int(s.length[i]) > len(in) {
			s.length[i] = uint32(len(in) - tlen)
		}
		tlen += int(s.length[i])
	}
	if s.numRecords > 0 && tlen < len(in) {
		s.length[s.numRecords-1] += uint32(len(in) - tlen)
	}

	fqzQualStats(s, in, pm, &qhist)

	pm.storeQMap = pm.nsym <= 8 && pm.nsym*2 < pm.maxSym

	// Fixed-length check.
	firstLen := uint32(0)
	if s.numRecords > 0 {
		firstLen = s.length[0]
	}
	idx := 1
	for ; idx < s.numRecords; idx++ {
		if s.length[idx] != firstLen {
			break
		}
	}
	pm.fixedLen = idx == s.numRecords
	pm.useQTab = false

	if strat < fqzNStrats-1 {
		if pm.pshift < 0 {
			l0 := 1.0
			if s.numRecords > 0 {
				l0 = float64(s.length[0])
			}
			pm.pshift = maxInt(0, int(math.Log(l0/float64(int(1)<<uint(pm.pbits)))/math.Log(2)+0.5))
		}

		if pm.nsym <= 4 {
			pm.qshift = 2
			if len(in) < 5000000 {
				pm.pbits = 2
				pm.pshift = 5
			}
		} else if pm.nsym <= 8 {
			pm.qbits = minInt(pm.qbits, 9)
			pm.qshift = 3
			if len(in) < 5000000 {
				pm.qbits = 6
			}
		}

		if len(in) < 300000 {
			pm.qbits = pm.qshift
			pm.dbits = 2
		}
	}

	for i := range dsqr {
		if dsqr[i] > uint32((1<<uint(pm.dbits))-1) {
			dsqr[i] = uint32((1 << uint(pm.dbits)) - 1)
		}
	}

	if pm.storeQMap {
		j := 0
		for i := 0; i < 256; i++ {
			if qhist[i] != 0 {
				pm.qmap[i] = j
				j++
			} else {
				pm.qmap[i] = intMax
			}
		}
		pm.maxSym = pm.nsym
	} else {
		pm.nsym = 255
		for i := 0; i < 256; i++ {
			pm.qmap[i] = i
		}
	}
	if gp.maxSym < pm.maxSym {
		gp.maxSym = pm.maxSym
	}

	if pm.qbits != 0 {
		for i := 0; i < 256; i++ {
			pm.qtab[i] = uint32(i)
		}
	}
	pm.qmask = (1 << uint(pm.qbits)) - 1

	if pm.pbits != 0 {
		for i := 0; i < 1024; i++ {
			pm.ptab[i] = uint32(minInt((1<<uint(pm.pbits))-1, i>>uint(pm.pshift)))
		}
	}

	if pm.dbits != 0 {
		for i := 0; i < 256; i++ {
			pm.dtab[i] = dsqr[minInt(len(dsqr)-1, i>>uint(pm.dshift))]
		}
	}

	pm.usePTab = pm.pbits > 0
	pm.useDTab = pm.dbits > 0

	pm.pflags = 0
	if pm.useQTab {
		pm.pflags |= fqzPFlagHaveQTab
	}
	if pm.useDTab {
		pm.pflags |= fqzPFlagHaveDTab
	}
	if pm.usePTab {
		pm.pflags |= fqzPFlagHavePTab
	}
	if pm.doSel {
		pm.pflags |= fqzPFlagDoSel
	}
	if pm.fixedLen {
		pm.pflags |= fqzPFlagDoLen
	}
	if pm.doDedup {
		pm.pflags |= fqzPFlagDoDedup
	}
	if pm.storeQMap {
		pm.pflags |= fqzPFlagHaveQMap
	}

	gp.maxSel = 0
	if pm.doSel {
		gp.maxSel = 1
		gp.gflags |= fqzGFlagHaveStab
	}

	if gp.maxSel != 0 && s.numRecords != 0 {
		mx := 0
		for i := 0; i < s.numRecords; i++ {
			if mx < int(s.flags[i]>>16) {
				mx = int(s.flags[i] >> 16)
			}
		}
		gp.maxSel = mx
	}

	return gp, nil
}

// fqzEncodeNewRead ports compress_new_read: it emits the per-read header
// (selector, length, reverse flag, dedup bit) and resets the state. It
// returns true if the read was emitted as a dedup of its predecessor.
func fqzEncodeNewRead(s *fqzSlice, st *fqzState, gp *fqzGParams,
	model *fqzModel, rc *rangeEncoder, in []byte, ip *int) bool {

	rec := st.rec
	i := *ip

	if gp.p[0].doSel || (gp.gflags&fqzGFlagMultiParam != 0) {
		if rec < s.numRecords {
			st.s = int(s.flags[rec] >> 16)
		} else {
			st.s = 0
		}
		model.sel.encodeSymbol(rc, uint16(st.s))
	} else {
		st.s = 0
	}

	x := st.s
	if gp.gflags&fqzGFlagHaveStab != 0 {
		x = int(gp.stab[st.s])
	}
	pm := &gp.p[x]

	length := int(s.length[rec])
	if !pm.fixedLen || st.firstLen {
		model.length[0].encodeSymbol(rc, uint16((length>>0)&0xff))
		model.length[1].encodeSymbol(rc, uint16((length>>8)&0xff))
		model.length[2].encodeSymbol(rc, uint16((length>>16)&0xff))
		model.length[3].encodeSymbol(rc, uint16((length>>24)&0xff))
		st.firstLen = false
	}

	if gp.gflags&fqzGFlagDoRev != 0 {
		if s.flags[rec]&fqzFReverse != 0 {
			model.revcomp.encodeSymbol(rc, 1)
		} else {
			model.revcomp.encodeSymbol(rc, 0)
		}
	}

	st.rec++
	st.p = length
	st.delta = 0
	st.qctx = 0
	st.prevq = 0
	st.ctx = uint32(pm.context)

	if pm.doDedup {
		if i != 0 && length == st.lastLen &&
			bytesEqual(in[i-st.lastLen:i], in[i:i+length]) {
			model.dup.encodeSymbol(rc, 1)
			i += length - 1
			st.p = 0
			*ip = i
			return true
		}
		model.dup.encodeSymbol(rc, 0)
		st.lastLen = length
	}

	*ip = i
	return false
}

// FQZCompEncode compresses a raw quality-score buffer with the fqzcomp
// codec (CRAM block compression method 7), using strategy level strat
// (0..3) and the per-read metadata in s. The output is byte-identical
// to htscodecs fqz_compress for the same inputs and strategy.
//
// When s is nil the buffer is treated as a single read.
func FQZCompEncode(in []byte, strat int, s *fqzSlice) ([]byte, error) {
	if len(in) > math.MaxInt32 {
		return nil, fmt.Errorf("fqzcomp: input too large")
	}
	if strat < 0 || strat > fqzMaxStrat {
		return nil, fmt.Errorf("fqzcomp: strategy %d out of range 0..%d", strat, fqzMaxStrat)
	}
	if s == nil {
		s = &fqzSlice{
			numRecords: 1,
			length:     []uint32{uint32(len(in))},
			flags:      []uint32{0},
		}
	}
	// Work on a copy of the length/flags: fqz_pick_parameters and
	// fqz_qual_stats mutate them (oversize clamping, selector bits).
	work := &fqzSlice{
		numRecords: s.numRecords,
		length:     append([]uint32(nil), s.length...),
		flags:      append([]uint32(nil), s.flags...),
	}

	gp, err := fqzPickParameters(strat, 4, work, in)
	if err != nil {
		return nil, err
	}

	out := fqzPutU32(nil, uint32(len(in)))
	out = fqzStoreParameters(out, gp)

	// Pre-shift ptab / dtab as compress_block_fqz2f does.
	for i := range gp.p {
		pm := &gp.p[i]
		for j := range pm.ptab {
			pm.ptab[j] <<= uint(pm.ploc)
		}
		for j := range pm.dtab {
			pm.dtab[j] <<= uint(pm.dloc)
		}
	}

	model := newFQZModel(gp)
	rc := newRangeEncoder()

	st := fqzState{firstLen: true}
	pm := &gp.p[0]

	for i := 0; i < len(in); {
		if st.p == 0 {
			if st.rec >= work.numRecords || work.length[st.rec] <= 0 {
				return nil, fmt.Errorf("fqzcomp: record count mismatch at byte %d", i)
			}
			if fqzEncodeNewRead(work, &st, gp, model, rc, in, &i) {
				i++
				continue
			}
		}
		for st.p > 0 {
			q := pm.qmap[in[i]]
			ctx := st.ctx
			st.ctx = fqzUpdateCtx(pm, &st, q)
			model.qual[ctx].encodeSymbol(rc, uint16(q))
			i++
		}
	}

	out = append(out, rc.finish()...)
	return out, nil
}

// fqzMaxStrat is FQZ_MAX_STRAT: the highest strategy level the encoder
// accepts.
const fqzMaxStrat = 3
