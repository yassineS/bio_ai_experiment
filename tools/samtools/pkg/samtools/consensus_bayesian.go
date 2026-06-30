package samtools

import (
	"math"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// This file ports the Gap5-derived bayesian consensus caller from upstream
// samtools (bam_consensus.c::calculate_consensus_gap5 / _gap5m and the
// nm_init / nm_local / poly_len localised-MAPQ machinery). It implements
// the four bayesian modes (BAYES_116, RECALL, PRECISE, MIXED); MODE_RECALL
// is the upstream default and the mode the bare `samtools consensus`
// invocation runs.

// bayesianMode enumerates the upstream MODE_* constants for the bayesian
// caller. MODE_SIMPLE (0) is handled by the frequency-counting path and is
// not represented here.
type bayesianMode int

const (
	// modeBayes116 reproduces samtools 1.16 (no separate indel parameter).
	modeBayes116 bayesianMode = 1
	// modeRecall is the Gap5 parameter set; the upstream default.
	modeRecall bayesianMode = 2
	// modePrecise favours precision (fewer false positives) over recall.
	modePrecise bayesianMode = 3
	// modeMixed runs recall and precise together and blends the calls.
	modeMixed bayesianMode = 4
)

// Upstream P_HET / P_INDEL / P_HET_SCALE / P_HOMOPOLY defaults.
const (
	defaultPHet     = 1e-3
	defaultPIndel   = 2e-4
	defaultHetScale = 1.0
	defaultPHomopol = 0.5
)

// bayesOptions carries the bayesian-specific knobs, mirroring the relevant
// fields of upstream's consensus_opts struct.
type bayesOptions struct {
	mode        bayesianMode
	useMQual    bool
	adjQual     bool
	nmAdjust    bool
	nmHalo      int
	scCost      int
	scaleMQual  float64
	lowMQual    int
	highMQual   int
	defaultQual int
	minQual     int // --min-BQ
	minDepth    int
	consCutoff  int
	ambig       bool
	pHet        float64
	pIndel      float64
	hetScale    float64
	homopolyFix float64 // 0 disables; otherwise the poly_adj multiplier
	homopolyRed float64 // homopoly-redux: poly_mul for RECALL mode
}

// consProbs holds the precomputed log-probability matrices for one
// parameter set, mirroring upstream's cons_probs struct.
type consProbs struct {
	lprior15 [15]float64
	pMM      [101]float64
	pxx      [101]float64
	pxM      [101]float64
	pox      [101]float64
	poM      [101]float64
	poo      [101]float64
	puu      [101]float64
	pum      [101]float64
	pmm      [101]float64
	polyMul  float64
}

// consResult holds the per-position bayesian result, mirroring upstream's
// consensus_t struct.
type consResult struct {
	call      int // A=0 C=1 G=2 T=3 *=4
	hetCall   int // 5x5 index: "ACGT*"[het/5] / "ACGT*"[het%5]
	hetLogOdd int
	phred     int
	depth     int
}

// Numerical constants from upstream bam_consensus.c.
const (
	tenLog2OverLog10 = 3.0103
	dblMin           = 2.2250738585072014e-308
)

// minEExp mirrors upstream's `DBL_MIN_EXP * log(2) + 1` guard used to
// detect underflow before exp().
var minEExp = float64(-1021)*math.Ln2 + 1

// consensusInit builds a consProbs matrix from the parameter set, mirroring
// upstream consensus_init. qcalS/qcalU/qcalO are the substitution/undercall/
// overcall quality calibration maps (the FLAT identity table by default).
func consensusInit(pHet, pIndel, hetScale, polyMul float64,
	qcalS, qcalU, qcalO *[101]int, mode bayesianMode) *consProbs {

	cp := &consProbs{polyMul: polyMul}

	// Priors: 25-entry ACGT* x ACGT* matrix, then folded to 15 combos.
	var prior [25]float64
	for i := 0; i < 25; i++ {
		prior[i] = pHet / 6
	}
	prior[0], prior[6], prior[12], prior[18], prior[24] = 1, 1, 1, 1, 1
	for i := 4; i < 24; i += 5 {
		prior[i] = pIndel / 6
	}
	for i := 20; i < 24; i++ {
		prior[i] = pIndel / 6
	}
	idx15 := [15]int{0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 14, 18, 19, 24}
	for j, k := range idx15 {
		cp.lprior15[j] = math.Log(prior[k])
	}

	for i := 1; i < 101; i++ {
		prob := 1 - math.Pow(10, -float64(qcalS[i])/10.0)
		cp.pMM[i] = math.Log(prob)
		cp.pxx[i] = math.Log((1 - prob) / 3)
		cp.pxM[i] = math.Log((math.Exp(cp.pMM[i]) + math.Exp(cp.pxx[i])) / 2)
		cp.pxM[i] += math.Log(hetScale)

		if mode == modeBayes116 {
			cp.pmm[i] = cp.pMM[i]
			cp.poM[i] = cp.pxM[i]
			cp.pum[i] = cp.pxM[i]
			cp.pox[i] = cp.pxx[i]
			cp.poo[i] = cp.pxx[i]
			cp.puu[i] = cp.pxx[i]
			continue
		}

		prob = 1 - math.Pow(10, -float64(qcalO[i])/10.0)
		cp.poo[i] = math.Log((1 - prob) / 3)
		if cp.poo[i] > cp.pMM[i]-.5 {
			cp.poo[i] = cp.pMM[i] - .5
		}
		cp.pox[i] = math.Log((math.Exp(cp.poo[i]) + math.Exp(cp.pxx[i])) / 2)
		cp.poM[i] = math.Log((math.Exp(cp.poo[i]) + math.Exp(cp.pMM[i])) / 2)
		if cp.poM[i] > cp.pxM[i]+.5 {
			cp.poM[i] = cp.pxM[i] + .5
		}

		prob = 1 - math.Pow(10, -float64(qcalU[i])/10.0)
		cp.pmm[i] = math.Log(prob)
		cp.puu[i] = math.Log((1 - prob) / 3)
		if cp.puu[i] > cp.pMM[i]-.5 {
			cp.puu[i] = cp.pMM[i] - .5
		}
		cp.pum[i] = math.Log((math.Exp(cp.puu[i]) + math.Exp(cp.pmm[i])) / 2)
	}

	cp.pMM[0] = cp.pMM[1]
	cp.pxx[0] = cp.pxx[1]
	cp.pxM[0] = cp.pxM[1]
	cp.pmm[0] = cp.pmm[1]
	cp.poo[0] = cp.poo[1]
	cp.pox[0] = cp.pox[1]
	cp.poM[0] = cp.poM[1]
	cp.puu[0] = cp.puu[1]
	cp.pum[0] = cp.pum[1]

	return cp
}

// flatQCal is the identity (FLAT) quality calibration table; index i maps
// to min(i,99). Upstream's static_qcal[QCAL_FLAT] is the identity 0..99.
var flatQCal = func() [101]int {
	var t [101]int
	for i := 0; i < 101; i++ {
		if i < 100 {
			t[i] = i
		} else {
			t[i] = 99
		}
	}
	return t
}()

// bayesProbSet bundles the recall and precise parameter sets needed by
// calculate_consensus_gap5m. cpPrecise is nil unless the mode needs it.
type bayesProbSet struct {
	recall  *consProbs
	precise *consProbs
}

// buildBayesProbSet constructs the consProbs matrices for a given mode,
// mirroring the consensus_init calls in upstream main_consensus.
func buildBayesProbSet(o bayesOptions) bayesProbSet {
	q := &flatQCal
	var ps bayesProbSet
	switch o.mode {
	case modePrecise:
		ps.precise = consensusInit(o.pHet, o.pIndel, 0.3*o.hetScale,
			o.homopolyRed, q, q, q, modePrecise)
	case modeMixed:
		ps.precise = consensusInit(math.Pow(o.pHet, 0.7), math.Pow(o.pIndel, 0.7),
			0.3*o.hetScale, o.homopolyRed, q, q, q, modePrecise)
	}
	recallPoly := 0.01
	if o.mode == modeRecall {
		recallPoly = o.homopolyRed
	}
	ps.recall = consensusInit(o.pHet, o.pIndel, o.hetScale, recallPoly,
		q, q, q, modeRecall)
	return ps
}

// bayesRead holds the per-read precomputed state needed by the bayesian
// caller: the localised NM array and (homopoly-fixed) qualities.
type bayesRead struct {
	// localNM[i] packs the SNP-score adjustment in bits 0..23 and a
	// poly-X length in bits 24+, exactly like upstream's local_nm[].
	localNM []int32
	// qual is the per-base quality, post homopoly-fix if enabled.
	qual []byte
}

// nmInit precomputes the localised-NM array for one read, porting upstream
// nm_init. It returns nil when use-MQ is disabled or the read has no SEQ.
func nmInit(rec *sam.Record, o bayesOptions) *bayesRead {
	if !o.useMQual {
		return nil
	}
	qlen := len(rec.Seq)
	if qlen <= 0 {
		return nil
	}
	localNM := make([]int32, qlen)
	// Copy qualities so homopoly-fix doesn't mutate the shared record.
	qual := make([]byte, qlen)
	copy(qual, rec.Qual)
	for len(qual) < qlen {
		qual = append(qual, 0)
	}

	seqi := func(i int) int { return baseToSeqi(upper(rec.Seq[i])) }

	polyAdj := 1.0
	if o.homopolyFix != 0 {
		polyAdj = o.homopolyFix
	}

	if o.adjQual {
		const qhalo = 8
		const qhalop = 2
		qmin := int(qual[0])
		qminp := int(qual[0])
		base := seqi(0)
		polyl, polyr := 0, 0

		i := 1
		for ; i < qlen; i++ {
			if seqi(i) != base {
				break
			}
			if i < qhalop && qminp > int(qual[i]) {
				qminp = int(qual[i])
			}
		}

		i = 0
		for ; i < qlen && i < qhalo; i++ {
			if qmin > int(qual[i]) {
				qmin = int(qual[i])
			}
		}

		for ; i < qlen-qhalo; i++ {
			if o.homopolyFix != 0 && seqi(i) != base {
				polyl = i
				base = seqi(i)
				qminp = int(qual[i])
				j := i + 1
				for ; j < qlen; j++ {
					if seqi(j) != base {
						break
					}
					if i < qhalop && qminp > int(qual[j]) {
						qminp = int(qual[j])
					}
				}
				polyr = j - 1
			} else {
				polyr = polyl
			}
			pl := polyr - polyl

			var t int
			if o.mode == modeBayes116 {
				t = (int(qual[i]) + 5*qmin) / 4
			} else {
				t = int(qual[i])/3 + int(float64(qminp-pl*2)*polyAdj)
			}
			if t < int(qual[i]) {
				localNM[i] += int32(int(qual[i]) - t)
			}

			qminp = int(qual[i])
			lo := polyl
			if i-qhalop > lo {
				lo = i - qhalop
			}
			hi := polyr
			if i+qhalop < hi {
				hi = i + qhalop
			}
			for k := lo; k <= hi; k++ {
				if k >= 0 && k < qlen && qminp > int(qual[k]) {
					qminp = int(qual[k])
				}
			}

			if qmin > int(qual[i+qhalo]) {
				qmin = int(qual[i+qhalo])
			} else if qmin <= int(qual[i-qhalo]) {
				qmin = 99
				for j := i - qhalo + 1; j <= i+qhalo; j++ {
					if qmin > int(qual[j]) {
						qmin = int(qual[j])
					}
				}
			}
		}
		for ; i < qlen; i++ {
			var t int
			if o.mode == modeBayes116 {
				t = (int(qual[i]) + 5*qmin) / 4
			} else {
				t = int(qual[i])/3 + int(float64(qminp)*polyAdj)
			}
			if t < int(qual[i]) {
				localNM[i] += int32(int(qual[i]) - t)
			}
		}
	}

	if o.homopolyFix != 0 {
		homopolyQualFix(rec.Seq, qual)
	}

	// Poly-X length recorded in the top byte of local_nm.
	for i := 0; i < qlen; {
		base := seqi(i)
		j := i + 1
		for ; j < qlen; j++ {
			if seqi(j) != base {
				break
			}
		}
		poly := j - i - 1
		if poly > 100 {
			poly = 100
		}
		for k := i; k < j; k++ {
			if k >= 0 && k < qlen {
				cur := int(localNM[k] >> 24)
				if poly > cur {
					cur = poly
				}
				localNM[k] = int32(cur<<24) | (localNM[k] & ((1 << 24) - 1))
			}
		}
		i = j
	}

	// Soft-clip cost.
	halo := o.nmHalo
	cig := rec.Cigar
	ncig := len(cig)
	if ncig > 0 {
		firstSC := cig[0].Op() == sam.CigarSoftClip ||
			(cig[0].Op() == sam.CigarHardClip && ncig > 1 && cig[1].Op() == sam.CigarSoftClip)
		if firstSC {
			i := 0
			for ; i < halo && i < qlen; i++ {
				localNM[i] += int32(o.scCost)
			}
			for ; i < halo*2 && i < qlen; i++ {
				localNM[i] += int32(o.scCost >> 1)
			}
		}
		lastSC := cig[ncig-1].Op() == sam.CigarSoftClip ||
			(cig[ncig-1].Op() == sam.CigarHardClip && ncig > 1 && cig[ncig-2].Op() == sam.CigarSoftClip)
		if lastSC {
			i := qlen - 1
			for ; i >= qlen-halo && i >= 0; i-- {
				localNM[i] += int32(o.scCost)
			}
			for ; i >= qlen-halo*2 && i >= 0; i-- {
				localNM[i] += int32(o.scCost >> 1)
			}
		}
	}

	// MD-tag substitution costs.
	md, ok := rec.GetAux("MD")
	if ok {
		if s, isStr := md.String(); isStr {
			applyMDCosts(localNM, s, halo, qlen)
		}
	}

	return &bayesRead{localNM: localNM, qual: qual}
}

// applyMDCosts walks the MD tag and bumps localNM around each substitution,
// porting the MD loop of upstream nm_init.
func applyMDCosts(localNM []int32, md string, halo, qlen int) {
	pos := 0
	i := 0
	for i < len(md) {
		c := md[i]
		if c >= '0' && c <= '9' {
			n := 0
			for i < len(md) && md[i] >= '0' && md[i] <= '9' {
				n = n*10 + int(md[i]-'0')
				i++
			}
			pos += n
			continue
		}
		if c == '^' {
			i++
			for i < len(md) && !(md[i] >= '0' && md[i] <= '9') {
				i++
			}
			continue
		}
		// Substitution.
		k := pos - halo*2
		if k < 0 {
			k = 0
		}
		for ; k < pos-halo && k < qlen; k++ {
			localNM[k] += 5
		}
		for ; k < pos+halo && k < qlen; k++ {
			localNM[k] += 10
		}
		for ; k < pos+halo*2 && k < qlen; k++ {
			localNM[k] += 5
		}
		// Upstream's MD loop advances only the md-string cursor after a
		// substitution (md++); it does NOT advance pos. Consecutive
		// substitutions therefore share the same query offset, so the
		// halo bands stack at one centre rather than sliding one base per
		// substitution. Incrementing pos here over-counted the +5 outer
		// band by one base per extra substitution, inflating localNM at a
		// read's mismatch-dense left edge and depressing the adjusted MAPQ
		// just enough to mask a coverage-island's first callable base.
		i++
	}
}

// homopolyQualFix redistributes qualities within homopolymer runs, porting
// upstream homopoly_qual_fix. seq is the read sequence; qual is mutated.
func homopolyQualFix(seq string, qual []byte) {
	for i := 0; i < len(seq); i++ {
		s := i
		base := upper(seq[i])
		for i+1 < len(seq) && upper(seq[i+1]) == base {
			i++
		}
		if s == i {
			continue
		}
		for j, k := s, i; j < k; j, k = j+1, k-1 {
			e := math.Pow(10, -float64(qual[j])/10.0) + math.Pow(10, -float64(qual[k])/10.0)
			v := -fastLog2(e/2)*3.0104 + .49
			if v < 0 {
				v = 0
			}
			if v > 255 {
				v = 255
			}
			qual[j] = byte(v)
			qual[k] = byte(v)
		}
	}
}

// nmLocal returns the localised NM figure within +/- halo of seqOffset,
// porting upstream nm_local. seqOffset is the 0-based query index.
func nmLocal(br *bayesRead, seqOffset int) float64 {
	if br == nil {
		return 0
	}
	nm := br.localNM
	if len(nm) == 0 {
		return 0
	}
	if seqOffset < 0 {
		return float64(nm[0] & ((1 << 24) - 1))
	}
	if seqOffset >= len(nm) {
		return float64(nm[len(nm)-1] & ((1 << 24) - 1))
	}
	return float64(nm[seqOffset]&((1<<24)-1)) / 10.0
}

// polyLen returns the poly-X length recorded for the given query offset,
// porting upstream poly_len.
func polyLen(br *bayesRead, seqOffset int) int {
	if br == nil {
		return 0
	}
	if seqOffset >= 0 && seqOffset < len(br.localNM) {
		return int(br.localNM[seqOffset] >> 24)
	}
	return 0
}

// fastLog2 is a fast approximate base-2 logarithm, ported from upstream's
// degree-3 Taylor implementation so phred conversions are byte-faithful.
func fastLog2(val float64) float64 {
	bits := math.Float64bits(val)
	e := int((bits>>52)&2047) - 1024
	bits &^= 2047 << 52
	bits += 1023 << 52
	m := math.Float64frombits(bits)
	return float64(e) + ((-1.0/3.0)*m+2)*m - 2.0/3.0
}

// phLog mirrors upstream's ph_log macro: -10*log10(x) computed via fastLog2.
func phLog(x float64) float64 {
	return -tenLog2OverLog10 * fastLog2(x)
}

// bayesPileupBase is one read's contribution to a column, pre-resolved for
// the bayesian caller.
type bayesPileupBase struct {
	base4    byte // sam nt16 code (or 16 for '*')
	qual     byte
	mapQ     uint8
	read     *bayesRead
	seqOff   int  // 0-based query offset of this base
	readPos0 int  // read's 0-based reference start
	refSkip  bool // true for N CIGAR positions (skipped)
}

// calculateConsensusGap5 ports upstream calculate_consensus_gap5: the core
// bayesian posterior caller. flags carries CONS_MQUAL when MQUAL is enabled.
func calculateConsensusGap5(bases []bayesPileupBase, useMQual bool, td int,
	o bayesOptions, cp *consProbs) consResult {

	// L: convert sam nt16 base to acgt*n order (0..5).
	// =ACM GRSV TWYH KDBN, then * bucket (16..31) -> 4.
	L := [32]int{
		5, 0, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	}

	var S [15]float64
	var counts [6]int
	depth := 0

	mapSing := [15]int{0, 5, 5, 5, 5, 1, 5, 5, 5, 2, 5, 5, 3, 5, 4}
	mapHet := [15]int{0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 14, 18, 19, 24}

	for _, p := range bases {
		if int(p.qual) < o.minQual {
			continue
		}
		if p.refSkip {
			continue
		}

		qual := int(p.qual)
		// default_qual substitution: upstream uses default when qual==255
		// or (qual==0 && first qual byte == 255). We treat a missing qual
		// (255) as the default; the all-0xff "*" QUAL case is normalised
		// to 0 quals upstream of here.
		if qual == 255 {
			qual = o.defaultQual
		}

		base := L[p.base4&31]

		if useMQual {
			mqual := float64(p.mapQ)
			if o.nmAdjust {
				// Upstream calls nm_local with seq_offset+1 — it indexes
				// the base just past the pileup cursor (bam_consensus.c:1382).
				// nmLocal clamps out-of-range offsets like nm_local's guard.
				mqual /= nmLocal(p.read, p.seqOff+1) + 1
				dt := td
				if dt > 30 {
					dt = 30
				}
				mqual *= 1 + 2*(0.5-float64(dt)/60.0)
			}
			mqual *= o.scaleMQual
			if mqual < float64(o.lowMQual) {
				mqual = float64(o.lowMQual)
			}
			if mqual > float64(o.highMQual) {
				mqual = float64(o.highMQual)
			}
			// Upstream bam_consensus.c:1402 indexes q2p[qual]; qual is a
			// uint8_t floored to >=1, and q2p[] is sized for the q range valid
			// SAM ever sees (Phred 0..93). C silently reads past its static
			// array for a non-conformant qual>100 (undefined behaviour, no
			// crash); Go panics on the [101] table instead, so clamp to 100 as
			// a crash-safety guard — byte-inert for valid SAM, but a malformed
			// BAM with qual>100 no longer aborts the process.
			if qual > 100 {
				qual = 100
			}
			P := q2pTab[qual]
			mi := int(mqual)
			if mi < 0 {
				mi = 0
			}
			if mi > 255 {
				mi = 255
			}
			M := mqualPow1m[mi]
			qual = int(phLog(P + .75*M - P*M))
		}

		// Upstream bam_consensus.c:1424-1425 floors qual to >=1 ("Quality 0
		// should never be permitted as it breaks the maths") and applies NO
		// upper clamp before cp->pMM[qual] (line 1426). Mirror that: floor
		// only.
		if qual < 1 {
			qual = 1
		}
		// Crash-safety guard (see the q2pTab note above): clamp qual to the
		// cp.p*[101] table bound. Byte-inert for valid SAM (Phred 0..93);
		// prevents an index-out-of-range panic on a non-conformant BAM.
		if qual > 100 {
			qual = 100
		}

		// Upstream poly_len is likewise called with seq_offset+1
		// (bam_consensus.c:1414); polyLen returns 0 when out of range.
		poly := float64(polyLen(p.read, p.seqOff+1))
		// Upstream bam_consensus.c:1423:
		//   int qual2 = MAX(1, qual-(poly-2)*cp->poly_mul);
		// The subexpression qual-(poly-2)*poly_mul is a double, MAX(1,double)
		// keeps the double, and assignment to int qual2 truncates toward
		// zero. There is NO upper clamp on qual2.
		t := float64(qual) - (poly-2)*cp.polyMul
		qual2 := 0
		if t < 1 {
			qual2 = 1
		} else {
			qual2 = int(t)
		}
		// Crash-safety guard: qual2 indexes the cp.p*[101] tables, and the
		// poly term can push it past 100 on a non-conformant BAM (qual>100).
		// Clamp to keep the index in range. Byte-inert for valid SAM input.
		if qual2 > 100 {
			qual2 = 100
		}

		// CHANGE #3: do NOT re-parenthesise the cp.p*[qual]-xx grouping below
		// (upstream bam_consensus.c:1426-1434). xx is subtracted from each of
		// the eight log-prob terms in exactly this order; the grouping is
		// already byte-faithful.
		xx := cp.pxx[qual]
		MM := cp.pMM[qual] - xx
		xM := cp.pxM[qual] - xx
		oo := cp.poo[qual2] - xx
		oM := cp.poM[qual2] - xx
		ox := cp.pox[qual2] - xx
		uu := cp.puu[qual2] - xx
		um := cp.pum[qual2] - xx
		mm := cp.pmm[qual2] - xx

		counts[base]++

		switch base {
		case 0: // A
			S[0] += MM
			S[1] += xM
			S[2] += xM
			S[3] += xM
			S[4] += oM
			S[8] += ox
			S[11] += ox
			S[13] += ox
			S[14] += oo
		case 1: // C
			S[1] += xM
			S[5] += MM
			S[6] += xM
			S[7] += xM
			S[8] += oM
			S[4] += ox
			S[11] += ox
			S[13] += ox
			S[14] += oo
		case 2: // G
			S[2] += xM
			S[6] += xM
			S[9] += MM
			S[10] += xM
			S[11] += oM
			S[4] += ox
			S[8] += ox
			S[13] += ox
			S[14] += oo
		case 3: // T
			S[3] += xM
			S[7] += xM
			S[10] += xM
			S[12] += MM
			S[13] += oM
			S[4] += ox
			S[8] += ox
			S[11] += ox
			S[14] += oo
		case 4: // *
			S[0] += uu
			S[1] += uu
			S[2] += uu
			S[3] += uu
			S[4] += um
			S[5] += uu
			S[6] += uu
			S[7] += uu
			S[8] += um
			S[9] += uu
			S[10] += uu
			S[11] += um
			S[12] += uu
			S[13] += um
			S[14] += mm
		case 5: // N
			S[0] += MM
			S[1] += MM
			S[2] += MM
			S[3] += MM
			S[4] += oM
			S[5] += MM
			S[6] += MM
			S[7] += MM
			S[8] += oM
			S[9] += MM
			S[10] += MM
			S[11] += oM
			S[12] += MM
			S[13] += oM
			S[14] += oo
		}
		depth++
	}

	shift := -math.MaxFloat64
	max := -math.MaxFloat64
	maxHet := -math.MaxFloat64
	call, hetCall := 0, 0

	pure := func(j int) bool {
		return j == 0 || j == 5 || j == 9 || j == 12 || j == 14
	}

	for j := 0; j < 15; j++ {
		S[j] += cp.lprior15[j]
		if shift < S[j] {
			shift = S[j]
		}
		if !pure(j) {
			if maxHet < S[j] {
				maxHet = S[j]
				hetCall = j
			}
			continue
		}
		if max < S[j] {
			max = S[j]
			call = j
		}
	}

	for j := 0; j < 15; j++ {
		S[j] -= shift
		e := fastExp(S[j])
		if S[j] > minEExp {
			S[j] = e
		} else {
			S[j] = dblMin
		}
	}

	var norm [15]float64
	tot1, tot2 := 0.0, 0.0
	for j := 0; j < 15; j++ {
		norm[j] += tot1
		norm[14-j] += tot2
		tot1 += S[j]
		tot2 += S[14-j]
	}

	var res consResult
	if depth == 0 || depth == counts[5] {
		res.call = 4
		return res
	}
	res.depth = depth

	if norm[call] == 0 {
		norm[call] = dblMin
	}
	var ph float64
	if S[call] == 1 && norm[call] < .01 {
		ph = phLog(norm[call]) + .5
	} else {
		ph = phLog(1-S[call]/(norm[call]+S[call])) + .5
	}
	res.call = mapSing[call]
	if ph < 0 {
		res.phred = 0
	} else {
		res.phred = int(ph)
	}

	if norm[hetCall] == 0 {
		norm[hetCall] = dblMin
	}
	ph = tenLog2OverLog10*(fastLog2(S[hetCall])-fastLog2(norm[hetCall])) + .5
	res.hetCall = mapHet[hetCall]
	res.hetLogOdd = int(ph)

	return res
}

// fastExp ports upstream's fast_exp: a quantized exp() lookup. For
// |y| <= 50 it returns exp(trunc(y*10)/10); otherwise it clamps y to
// [-500, 500] and returns exp(trunc(y)). The 0.1-resolution quantization
// is load-bearing for byte-parity with upstream.
func fastExp(y float64) float64 {
	if y >= -50 && y <= 50 {
		return math.Exp(float64(int(y*10)) / 10.0)
	}
	if y < -500 {
		y = -500
	}
	if y > 500 {
		y = 500
	}
	return math.Exp(float64(int(y)))
}

// calculateConsensusGap5m ports upstream calculate_consensus_gap5m: it runs
// the gap5 caller once for non-mixed modes, twice (and blends) for MIXED.
func calculateConsensusGap5m(bases []bayesPileupBase, useMQual bool, td int,
	o bayesOptions, ps bayesProbSet) consResult {

	if o.mode != modeMixed {
		cp := ps.recall
		if o.mode == modePrecise {
			cp = ps.precise
		}
		return calculateConsensusGap5(bases, useMQual, td, o, cp)
	}

	consP := calculateConsensusGap5(bases, useMQual, td, o, ps.precise)
	consR := calculateConsensusGap5(bases, useMQual, td, o, ps.recall)

	cons := consP
	switch {
	case consP.phred > 0 && consR.phred > 0 && consP.call == consR.call:
		cons.phred += minInt(20, consR.phred)
	case consP.hetLogOdd >= 0 && consR.hetLogOdd >= 0 && consP.hetCall == consR.hetCall:
		cons.hetLogOdd += minInt(20, consR.hetLogOdd)
	case consP.hetLogOdd >= 0:
		q2 := maxInt(consR.phred, consR.hetLogOdd)
		cons.hetLogOdd = maxInt(1, cons.hetLogOdd-q2/2)
	case consR.hetLogOdd >= 70:
		q1 := consP.phred
		q2 := consR.hetLogOdd
		cons = consR
		// Upstream bam_consensus.c:1847:
		//   cons->het_logodd = MIN(15, MAX((q2-q1*2)/2, 1+q2/(q1+1.0)));
		// (q2-q1*2)/2 is INT division (q1,q2 are int), but 1+q2/(q1+1.0)
		// uses DOUBLE division (q1+1.0 is a double). MAX is taken on the
		// doubles, then truncated to int by the assignment.
		cons.hetLogOdd = minInt(15, int(math.Max(float64((q2-q1*2)/2), 1+float64(q2)/(float64(q1)+1.0))))
	case consR.hetLogOdd >= 0:
		q1 := consP.phred
		q2 := consR.hetLogOdd
		cons = consR
		eq := 0
		if consP.hetCall == consR.hetCall {
			eq = 5
		}
		// CHANGE #5: do NOT re-parenthesise. Upstream bam_consensus.c:1853-1854:
		//   cons->het_logodd = MAX(1, q2 - 0.3*q1) + 5*(consP.het_call==consR.het_call);
		// MAX(1, q2-0.3*q1) is computed on the double q2-0.3*q1 then truncated
		// to int (MAX of an int 1 and a double promotes to double); the eq
		// term (0 or 5) is added afterwards. The grouping is already faithful.
		cons.hetLogOdd = maxInt(1, int(float64(q2)-0.3*float64(q1))) + eq
		cons.phred = 0
	case consR.hetLogOdd < 0:
		consR.phred = consR.phred / 2
		if consR.phred > consP.phred {
			cons = consR
		}
		cons.phred = maxInt(10, cons.phred)
	}
	return cons
}

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

// bayesHetTab maps a heterozygous call index to its IUPAC ambiguity code,
// mirroring upstream consensus_base's static "het" string.
var bayesHetTab = []byte("AMRWa" + "MCSYc" + "RSGKg" + "WYKTt" + "acgt*")

// bayesCallToBase converts a consResult to the output base/qual pair,
// porting upstream consensus_base's bayesian branch.
func bayesCallToBase(cons consResult, o bayesOptions) (byte, int) {
	var cb byte
	var cq int
	hetTab := bayesHetTab
	switch {
	case cons.depth < o.minDepth && cons.call != 4:
		cb = 'N'
		cq = 0
	case cons.hetLogOdd > 0 && o.ambig:
		cb = hetTab[cons.hetCall]
		cq = cons.hetLogOdd
	default:
		cb = "ACGT*"[cons.call]
		cq = cons.phred
	}
	if cq < o.consCutoff && cb != '*' &&
		cons.hetCall%5 != 4 && cons.hetCall/5 != 4 {
		cb = 'N'
		cq = 0
	}
	return cb, cq
}
