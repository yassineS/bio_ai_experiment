// Consensus-based indel caller for `bcftools mpileup --indels-cns`
// (a.k.a. --indels-2.0). This is the Go port of upstream's
// reference_code/bcftools/bam2bcf_edlib.c — the alternative to the
// legacy probabilistic indel caller (bam2bcf_indel.c, see
// bam2bcf_indel.go / bam2bcf_indel_align.go).
//
// Algorithm overview (bam2bcf_edlib.c:1319 bcf_edlib_gap_prep):
//
//  1. Discover candidate indel sizes at the column via
//     bcfCgpFindTypes — the same helper the legacy path uses.
//  2. For each (sample, candidate type), build TWO per-sample consensus
//     haplotypes via bcfCgpConsensus — bam2bcf_edlib.c:316
//     bcf_cgp_consensus. cons[0] takes the majority allele in
//     heterozygous columns; cons[1] takes the alternative. Per-base
//     counts (cons_base) and per-offset insertion accumulators
//     (cons_ins) drive both threads, with low-coverage smoothing from
//     ref_base / ref_ins observations of reads NOT carrying the type
//     under evaluation.
//  3. For each read in the column, align its [qbeg,qend) window
//     against BOTH cons[0] and cons[1] via the in-tree edlib engine
//     (pkg/htsgo/edlib) in EDLIB_MODE_HW / EDLIB_TASK_LOC, and take
//     the better (lower) score (bam2bcf_edlib.c:971 bcf_cgp_align_score).
//  4. bcfCgpComputeIndelQCNS folds the per-(read,type) scores into
//     each read's p.aux word, using upstream's TMP_MAGIC=255 sentinel
//     (vs the legacy path's 111), indelQ1/indelQ2 with vs_ref
//     blending, and the optional poly_mqual homopolymer rescaling
//     (bam2bcf_edlib.c:1067-1302).
//
// The downstream emission pipeline (bcfCallGlfgenIndel →
// bcfCallCombineIndel → bcfCall2bcfIndel) is shared with the legacy
// path: only the score-assignment stage differs, plus a handful of
// CNS-only deviations in bam2bcf.go's glfgen (`opts.IndelsCNS` gates)
// and mpileup.go's IDV/IMF emission (recomputed from ADF/ADR).
//
// Upstream reference: reference_code/bcftools/bam2bcf_edlib.c.

package bcftools

import (
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/edlib"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/strfinder"
)

// cgpConsNI mirrors bam2bcf_edlib.c:251 NI — the per-offset cap on
// distinct insertion strings tracked by cons_ins / ref_ins. Hits past
// this point are silently dropped by upstream's bcf_cgp_append_cons.
const cgpConsNI = 100

// cgpMaxIns mirrors bam2bcf_edlib.c:315 MAX_INS — the per-insertion
// length cap used both when collecting bases and as the heti/hetd
// array bound during het-threading.
const cgpMaxIns = 8192

// Consensus tuning constants — bam2bcf_edlib.c:635-640.
const (
	cgpConsCutoff     = 0.40 // % needed for base vs N
	cgpConsCutoffDel  = 0.35 // % to include any het del
	cgpConsCutoff2    = 0.80 // % needed for gap in cons[1]
	cgpConsCutoffInc  = 0.35 // % to include any insertion cons[0]
	cgpConsCutoffInc2 = 0.80 // % to include any insertion cons[1] HOM
	cgpConsCutoffIns  = 0.60 // and then 60% needed for it to be bases vs N
)

// cgpTMPMagic is upstream's TMP_MAGIC for the indelQ length-norm
// dampener (bam2bcf_edlib.c:1215). It is 255 in the edlib-flavored
// compute_indelQ vs 111 in the legacy bam2bcf_indel.c version.
const cgpTMPMagic = 255.0

// strFreqEntry is one (string, freq) accumulator in the per-offset
// cons_ins / ref_ins map. Mirrors the (str[NI], len[NI], freq[NI])
// arrays of upstream's str_freq struct (bam2bcf_edlib.c:253-257).
type strFreqEntry struct {
	str  []byte
	freq int
}

// strFreq is a per-offset bag of up to NI distinct insertion strings
// with their observation counts.
type strFreq struct {
	entries []strFreqEntry
}

// append adds str/freq to sf. Same-content strings have their freq
// summed; new strings are appended unless already at the NI cap.
func (sf *strFreq) append(str []byte, freq int) bool {
	for i := range sf.entries {
		if len(sf.entries[i].str) == len(str) &&
			bytesEqualCNS(sf.entries[i].str, str) {
			sf.entries[i].freq += freq
			return true
		}
	}
	if len(sf.entries) >= cgpConsNI {
		return false
	}
	cp := make([]byte, len(str))
	copy(cp, str)
	sf.entries = append(sf.entries, strFreqEntry{str: cp, freq: freq})
	return true
}

func bytesEqualCNS(a, b []byte) bool {
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

// cgpBase6 maps ASCII A/C/G/T (case-insensitive), U/u (=T) and '*'
// (=5) to the 0..5 alphabet upstream uses for cons_base accumulators
// (bam2bcf_edlib.c:324-335). Anything else (including N) maps to 4.
func cgpBase6(b byte) int {
	switch b {
	case 'A', 'a':
		return 0
	case 'C', 'c':
		return 1
	case 'G', 'g':
		return 2
	case 'T', 't', 'U', 'u':
		return 3
	case '*', '^':
		return 5
	}
	return 4
}

// seqNt16ToBase6 collapses a 4-bit IUPAC seq_nt16 code (0..15) into
// the 0..4 range upstream uses for consensus counters.
func seqNt16ToBase6(seqNt16 int) int {
	return seqNt16Int[seqNt16&15]
}

// consensusResult is the return value of bcfCgpConsensus.
type consensusResult struct {
	cons       [2][]byte // 2-bit-base consensus haplotypes
	tconLen    [2]int    // valid lengths in cons[0] / cons[1]
	leftShift  int       // net inserted(+) or deleted(-) bases before pos
	rightShift int       // net inserted(+) or deleted(-) bases after pos
	cposPos    int       // index in cons[] of the column right after pos
	band       int       // running max |insert-delete| seen in this sample
}

// bcfCgpConsensus is the Go port of bam2bcf_edlib.c:316 bcf_cgp_consensus.
func bcfCgpConsensus(pile []pileupBase, pos int, bca *bcfCallauxIndel,
	ref []byte, left, right, typeLen, biggestDel, posL, posR int) consensusResult {

	span := right - left
	if span <= 0 {
		return consensusResult{cons: [2][]byte{{}, {}}, tconLen: [2]int{0, 0}, cposPos: 0}
	}

	consBase := make([][6]int, span+1)
	refBase := make([][6]int, span+1)
	consIns := make([]strFreq, span+1)
	refIns := make([]strFreq, span+1)

	var (
		localBandMax = 0
		totalSpanStr = 0
		typeDepth    = 0
		bandLocalMax = 0
	)

	for i := range pile {
		p := &pile[i]
		if p.rec == nil {
			continue
		}
		rec := p.rec
		if rec.Flag&sam.FlagUnmapped != 0 {
			continue
		}
		x := int(rec.Pos) - 1 // ref coord
		y := 0                // query coord
		localBand := 0
		var insBuf [cgpMaxIns]byte
		for _, op := range rec.Cigar {
			o := op.Op()
			ln := int(op.Length())
			skipTo := 0
			switch o {
			case sam.CigarSoftClip:
				y += ln
			case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
				for j := 0; j < ln; j, x, y = j+1, x+1, y+1 {
					if x < left {
						continue
					}
					if x >= right {
						break
					}
					if y < 0 || y >= len(rec.Seq) {
						continue
					}
					base4 := nt16ASCII(rec.Seq[y])
					b6 := seqNt16ToBase6(base4)
					if p.indel == typeLen {
						consBase[x-left][b6]++
					} else if x != pos+1 {
						refBase[x-left][b6]++
					}
				}
			case sam.CigarInsertion:
				if x >= left && x < right {
					localBand += p.indel
					if localBandMax < localBand {
						localBandMax = localBand
					}
				}
				j := 0
				for ; j < ln; j, y = j+1, y+1 {
					if x < left {
						continue
					}
					if x >= right {
						break
					}
					if y < 0 || y >= len(rec.Seq) {
						continue
					}
					if j < cgpMaxIns {
						insBuf[j] = byte(seqNt16ToBase6(nt16ASCII(rec.Seq[y])))
					}
				}
				if x >= left && x < right {
					ilen := j
					if ilen > cgpMaxIns {
						ilen = cgpMaxIns
					}
					if p.indel == typeLen {
						consIns[x-left].append(insBuf[:ilen], 1)
						if x == pos+1 {
							typeDepth++
						}
					} else if x != pos+1 {
						refIns[x-left].append(insBuf[:ilen], 1)
					}
				}
			case sam.CigarDeletion:
				if x >= left && x < right {
					localBand += p.indel
					if localBandMax < -localBand {
						localBandMax = -localBand
					}
				}
				for j := 0; j < ln; j, x = j+1, x+1 {
					if x < left {
						continue
					}
					if x >= right {
						break
					}
					matchesType := (p.indel == typeLen && !p.isDel) ||
						(p.indel == 0 && p.isDel && ln == -typeLen)
					switch {
					case matchesType:
						consBase[x-left][5]++
						if x == pos+1 {
							typeDepth++
						}
					case x+ln <= pos+1 || (skipTo != 0 && x > skipTo):
						refBase[x-left][5]++
					case x <= pos && x+ln > pos+1:
						if x > skipTo {
							skipTo = x + ln
						}
					}
				}
			case sam.CigarSkipped:
				x += ln
			}
		}
		if int(rec.Pos)-1 <= posL && x >= posR {
			totalSpanStr++
		}
		if bandLocalMax < localBandMax {
			bandLocalMax = localBandMax
		}
	}

	// Smooth cons_base/ins with ref_base/ref_ins (bam2bcf_edlib.c:488-553).
	for i := 0; i < span; i++ {
		t := consBase[i][0] + consBase[i][1] + consBase[i][2] +
			consBase[i][3] + consBase[i][4] + consBase[i][5]
		for _, e := range consIns[i].entries {
			t += e.freq
		}
		r := refBase[i][0] + refBase[i][1] + refBase[i][2] +
			refBase[i][3] + refBase[i][4] + refBase[i][5]
		for _, e := range refIns[i].entries {
			r += e.freq
		}
		rfract := float64(r-t*2) * 0.75 / float64(r+1)
		floor := 1.01 / (float64(r) + 1e-10)
		if rfract < floor {
			rfract = floor
		}
		if i+left >= pos+1 && i+left < pos+1-biggestDel {
			continue
		}
		for k := 0; k < 6; k++ {
			consBase[i][k] += int(rfract * float64(refBase[i][k]))
		}
		for _, e := range refIns[i].entries {
			consIns[i].append(e.str, int(rfract*float64(e.freq)))
		}
	}

	// Allocate cons buffer to worst-case length.
	maxLen := span
	for i := 0; i < span; i++ {
		insMax := 0
		for _, e := range consIns[i].entries {
			if e.freq == 0 {
				continue
			}
			if insMax < len(e.str) {
				insMax = len(e.str)
			}
		}
		maxLen += insMax
	}
	if typeLen > 0 {
		maxLen += typeLen
	}
	cons := [2][]byte{make([]byte, maxLen+1), make([]byte, maxLen+1)}

	// Merge insertions of same length (bam2bcf_edlib.c:580-633).
	for i := 0; i < span; i++ {
		for j := 0; j < len(consIns[i].entries); j++ {
			if consIns[i].entries[j].freq == 0 {
				continue
			}
			jl := len(consIns[i].entries[j].str)
			tally := make([][5]int, jl)
			for l := 0; l < jl; l++ {
				b := consIns[i].entries[j].str[l]
				if int(b) >= 0 && int(b) < 5 {
					tally[l][b] = consIns[i].entries[j].freq
				}
			}
			for k := j + 1; k < len(consIns[i].entries); k++ {
				if consIns[i].entries[k].freq == 0 ||
					len(consIns[i].entries[k].str) != jl {
					continue
				}
				for l := 0; l < jl; l++ {
					b := consIns[i].entries[k].str[l]
					if int(b) >= 0 && int(b) < 5 {
						tally[l][b] += consIns[i].entries[k].freq
					}
				}
				consIns[i].entries[j].freq += consIns[i].entries[k].freq
				consIns[i].entries[k].freq = 0
			}
			for l := 0; l < jl; l++ {
				var (
					maxV = 0
					base = 0
					tot  = tally[l][0] + tally[l][1] + tally[l][2] +
						tally[l][3] + tally[l][4]
				)
				for k := 0; k < 5; k++ {
					if maxV < tally[l][k] {
						maxV = tally[l][k]
						base = k
					}
				}
				if float64(maxV) > 0.6*float64(tot) {
					consIns[i].entries[j].str[l] = byte(base)
				} else {
					consIns[i].entries[j].str[l] = 4
				}
			}
		}
	}

	// Per-position majority call walk (bam2bcf_edlib.c:642-791).
	var (
		heti       [cgpMaxIns]int
		hetd       [cgpMaxIns]int
		cposPos    = -1
		leftShift  = 0
		rightShift = 0
		tconLen    [2]int
	)
	for cnum := 0; cnum < 2; cnum++ {
		k := 0
		for i := 0; i < span; i++ {
			if i >= pos-left+1 && cposPos == -1 {
				cposPos = k
			}

			var (
				maxV, maxV2 = 0, 0
				maxJ, maxJ2 = 4, 4
				tot         = 0
			)
			for j := 0; j < 6; j++ {
				if maxV < consBase[i][j] {
					maxV2 = maxV
					maxJ2 = maxJ
					maxV = consBase[i][j]
					maxJ = j
				} else if maxV2 < consBase[i][j] {
					maxV2 = consBase[i][j]
					maxJ2 = j
				}
				tot += consBase[i][j]
			}

			var (
				maxVIns, maxJIns = 0, 0
				totIns           = 0
			)
			for j := range consIns[i].entries {
				// Pad the j==0 entry at the candidate column up to typeLen
				// (bam2bcf_edlib.c:680-690).
				if typeLen > 0 && i+left == pos+1 &&
					len(consIns[i].entries[j].str) < typeLen && j == 0 {
					old := consIns[i].entries[j].str
					padded := make([]byte, typeLen)
					copy(padded, old)
					for p := len(old); p < typeLen; p++ {
						padded[p] = 4
					}
					consIns[i].entries[j].str = padded
				}
				if consIns[i].entries[j].freq == 0 {
					continue
				}
				if maxVIns < consIns[i].entries[j].freq {
					maxVIns = consIns[i].entries[j].freq
					maxJIns = j
				}
				totIns += consIns[i].entries[j].freq
			}

			totSum := tot
			alwaysIns := (i == pos-left+1 && typeLen > 0) ||
				float64(maxVIns) > cgpConsCutoffInc2*float64(totSum)
			hetIns := false
			if !alwaysIns && maxVIns >= bca.MinSupport {
				if cnum == 0 {
					hetIns = float64(maxVIns) > cgpConsCutoffInc*float64(totSum)
					if i < cgpMaxIns {
						if hetIns {
							heti[i] = 1
						} else if float64(maxVIns) > 0.3*float64(totSum) {
							heti[i] = -1
						} else {
							heti[i] = 0
						}
					}
				} else if i < cgpMaxIns {
					hetIns = heti[i] == -1
				}
			}
			if alwaysIns || hetIns {
				if maxJIns < len(consIns[i].entries) &&
					float64(maxVIns) > cgpConsCutoffIns*float64(totIns) {
					payload := consIns[i].entries[maxJIns].str
					for j := 0; j < len(payload); j++ {
						if cnum == 0 {
							if k < pos-left+leftShift {
								leftShift++
							} else {
								rightShift++
							}
						}
						if k < len(cons[cnum]) {
							cons[cnum][k] = payload[j]
						}
						k++
					}
				} else if maxJIns < len(consIns[i].entries) {
					ilen := len(consIns[i].entries[maxJIns].str)
					for j := 0; j < ilen; j++ {
						if k < len(cons[cnum]) {
							cons[cnum][k] = 4
						}
						k++
					}
				}
			}

			alwaysDel := (typeLen < 0 && i > pos-left && i <= pos-left-typeLen) ||
				float64(consBase[i][5]) > cgpConsCutoff2*float64(tot)
			hetDel := false
			if !alwaysDel && consBase[i][5] >= bca.MinSupport {
				if cnum == 0 {
					tot2 := tot
					if i > pos-left && i <= pos-left-biggestDel {
						tot2 = totalSpanStr - typeDepth
					}
					hetDel = float64(consBase[i][5]) >= cgpConsCutoffDel*float64(tot2)
					if i < cgpMaxIns {
						if i > pos-left && i <= pos-left-biggestDel {
							hetd[i] = 0
						} else if hetDel {
							hetd[i] = 1
						} else if float64(consBase[i][5]) >= 0.3*float64(tot2) {
							hetd[i] = -1
						} else {
							hetd[i] = 0
						}
					}
				} else if i < cgpMaxIns {
					hetDel = hetd[i] == -1
					if maxJ == 5 && !hetDel {
						maxV = maxV2
						maxJ = maxJ2
					}
				}
			}
			if alwaysDel || hetDel {
				if cnum == 0 {
					if k < pos-left+leftShift {
						leftShift--
					} else {
						rightShift++
					}
				}
			} else {
				if float64(maxV) > cgpConsCutoff*float64(tot) {
					if k < len(cons[cnum]) {
						cons[cnum][k] = byte(maxJ)
					}
					k++
				} else if maxV > 0 {
					if k < len(cons[cnum]) {
						cons[cnum][k] = 4
					}
					k++
				} else {
					var c byte
					if left+k < len(ref) {
						c = byte(cgpBase6(ref[left+k]))
					} else {
						c = 4
					}
					if k < len(cons[cnum]) {
						cons[cnum][k] = c
					}
					k++
				}
			}
		}
		tconLen[cnum] = k
	}

	if cposPos < 0 {
		cposPos = 0
	}

	return consensusResult{
		cons:       cons,
		tconLen:    tconLen,
		leftShift:  leftShift,
		rightShift: rightShift,
		cposPos:    cposPos,
		band:       bandLocalMax,
	}
}

// edlibGlocal is the Go port of bam2bcf_edlib.c:834 edlib_glocal.
func edlibGlocal(ref, query []byte, m, delBias float64) int {
	if len(ref) == 0 || len(query) == 0 {
		return b2bIndelMissingScore
	}
	cfg := edlib.Config{K: -1, Mode: edlib.ModeHW, Task: edlib.TaskLoc}
	res, err := edlib.Align(query, ref, cfg)
	if err != nil || res.EditDistance < 0 ||
		len(res.EndLocations) == 0 || len(res.StartLocations) == 0 {
		return b2bIndelMissingScore
	}
	tLen := res.EndLocations[0] - res.StartLocations[0] + 1
	return int(m * (float64(res.EditDistance) - delBias*float64(tLen-len(query))))
}

// bcfCgpAlignScoreCNS is the Go port of bam2bcf_edlib.c:971
// bcf_cgp_align_score. It scores one read window against both candidate
// consensus haplotypes (ref1=cons[0], ref2=cons[1]) and picks the
// better (lower) of the two.
func bcfCgpAlignScoreCNS(ref1, ref2, query []byte, atype int,
	indelBias, delBias float64) int {
	const mm = 30.0

	// Trim poly-N at both ends (bam2bcf_edlib.c:987-1006).
	tBeg := 0
	tEnd1 := len(ref1)
	tEnd2 := len(ref2)

	maxLen := tEnd1 - tBeg
	if tEnd2-tBeg < maxLen {
		maxLen = tEnd2 - tBeg
	}
	l := 0
	for ; l < maxLen; l++ {
		if ref1[l] != 4 || ref2[l] != 4 {
			break
		}
	}
	if l > atype {
		tBeg += l - atype
	}

	l = tEnd1 - tBeg - 1
	for ; l >= 0; l-- {
		if ref1[tBeg+l] != 4 {
			break
		}
	}
	trim := tEnd1 - tBeg - 1 - l
	if trim > atype {
		tEnd1 -= trim - atype
	}

	l = tEnd2 - tBeg - 1
	for ; l >= 0; l-- {
		if ref2[tBeg+l] != 4 {
			break
		}
	}
	trim = tEnd2 - tBeg - 1 - l
	if trim > atype {
		tEnd2 -= trim - atype
	}

	if tBeg < 0 {
		tBeg = 0
	}
	if tEnd1 < tBeg {
		tEnd1 = tBeg
	}
	if tEnd2 < tBeg {
		tEnd2 = tBeg
	}
	if tEnd1 > len(ref1) {
		tEnd1 = len(ref1)
	}
	if tEnd2 > len(ref2) {
		tEnd2 = len(ref2)
	}
	if tBeg > len(ref1) {
		tBeg = len(ref1)
	}
	if tBeg > len(ref2) {
		tBeg = len(ref2)
	}

	sc2 := edlibGlocal(ref2[tBeg:tEnd2], query, mm, delBias)
	sc1 := b2bIndelMissingScore
	if tEnd1 != tEnd2 || !bytesEqualCNS(ref1[tBeg:tEnd1], ref2[tBeg:tEnd2]) {
		sc1 = edlibGlocal(ref1[tBeg:tEnd1], query, mm, delBias)
	}

	switch {
	case sc1 >= b2bIndelMissingScore && sc2 >= b2bIndelMissingScore:
		return b2bIndelMissingScore
	case sc1 >= b2bIndelMissingScore:
	case sc2 >= b2bIndelMissingScore:
		sc2 = sc1
	default:
		if sc2 > sc1 {
			sc2 = sc1
		}
	}
	if sc2 < 0 {
		sc2 = 0
	}

	qlen := len(query)
	var lNorm int
	if qlen > 0 {
		base := 0.5 * (100.0*float64(sc2)/float64(qlen) + 0.499)
		lNorm = int(base * indelBias * 0.5)
	}
	if lNorm > 255 {
		lNorm = 255
	}
	if lNorm < 0 {
		lNorm = 0
	}
	return (sc2 << 8) | lNorm
}

// bcfCgpComputeIndelQCNS is the Go port of bam2bcf_edlib.c:1067
// bcf_cgp_compute_indelQ — the edlib-flavored variant of the legacy
// path's bcfCgpComputeIndelQ.
//
// Differences vs the legacy variant:
//   - TMP_MAGIC = 255 (vs 111).
//   - indelQ1 (vs-ref) and indelQ2 (vs-next-best non-ref) are
//     computed separately, then blended via bca.VsRef.
//   - When PolyMQual is on, both seqQ and indelQ are nudged by the
//     min base quality in the homopolymer surrounding qpos.
//   - The indelQ-vs-seqQ cap (`if indelQ > seqQ`) still holds.
func bcfCgpComputeIndelQCNS(piles [][]pileupBase, scores []int,
	bca *bcfCallauxIndel, inscns []byte, lRun, maxIns, refType int,
	types []int, qavg float64) int {
	nTypes := len(types)
	sc := make([]int, nTypes)
	sumq := make([]int, nTypes)

	K := 0
	for s := 0; s < len(piles); s++ {
		for i := 0; i < len(piles[s]); i, K = i+1, K+1 {
			p := &piles[s][i]
			off := K * nTypes
			for t := 0; t < nTypes; t++ {
				sc[t] = scores[off+t]<<6 | t
			}
			for t := 1; t < nTypes; t++ {
				for j := t; j > 0 && sc[j] < sc[j-1]; j-- {
					sc[j], sc[j-1] = sc[j-1], sc[j]
				}
			}

			var indelQ, indelQ1, seqQ int
			if sc[0]&0x3f == refType {
				indelQ = (sc[1] >> 14) - (sc[0] >> 14)
				seqQ = estSeqQ(bca, types[sc[1]&0x3f], lRun)
			} else {
				t := 0
				for ; t < nTypes; t++ {
					if sc[t]&0x3f == refType {
						break
					}
				}
				indelQ1 = (sc[t] >> 14) - (sc[0] >> 14)
				if t == nTypes {
					t = nTypes - 1
				}
				t2 := 1
				for ; t2 < nTypes; t2++ {
					if sc[t2]&0x3f != refType {
						break
					}
				}
				if t2 == nTypes {
					t2--
				}
				indelQ2 := (sc[t2] >> 14) - (sc[0] >> 14)
				seqQ = estSeqQ(bca, types[sc[0]&0x3f], lRun)
				indelQ = int(bca.VsRef*float64(indelQ1) + (1-bca.VsRef)*float64(indelQ2))
			}

			if bca.IndelBias > 0 {
				indelQ = int(float64(indelQ) / (bca.IndelBias * 0.5))
				indelQ1 = int(float64(indelQ1) / (bca.IndelBias * 0.5))
			}

			// poly_mqual rescaling (bam2bcf_edlib.c:1164-1203).
			if bca.PolyMQual != 0 && p.rec != nil {
				qpos := p.qpos
				seq := p.rec.Seq
				qual := p.rec.Qual
				if qpos >= 0 && qpos < len(qual) {
					minQ := int(qual[qpos])
					baseLIdx := qpos + 1
					if baseLIdx >= len(seq) {
						baseLIdx = qpos
					}
					if baseLIdx >= 0 && baseLIdx < len(seq) {
						baseL := nt16ASCII(seq[baseLIdx])
						for l := qpos; l >= 0; l-- {
							if l >= len(seq) {
								continue
							}
							if nt16ASCII(seq[l]) != baseL {
								break
							}
							if l < len(qual) && minQ > int(qual[l]) {
								minQ = int(qual[l])
							}
						}
					}
					if qpos+1 < len(seq) {
						baseR := nt16ASCII(seq[qpos+1])
						for l := qpos + 1; l < len(seq); l++ {
							if l < len(qual) && minQ > int(qual[l]) {
								minQ = int(qual[l])
							}
							if nt16ASCII(seq[l]) != baseR {
								break
							}
						}
					}
					adjSeqQ := minQ - int(qavg/10.0)
					if cap20 := int(qavg / 20.0); adjSeqQ > cap20 {
						adjSeqQ = cap20
					}
					adjIndel := minQ - int(qavg/5.0)
					if cap20 := int(qavg / 20.0); adjIndel > cap20 {
						adjIndel = cap20
					}
					seqQ += adjSeqQ
					indelQ += adjIndel
					indelQ1 += adjIndel
					if seqQ < 0 {
						seqQ = 0
					}
					if indelQ < 0 {
						indelQ = 0
					}
					if indelQ1 < 0 {
						indelQ1 = 0
					}
				}
			}

			// Length-norm dampener with TMP_MAGIC=255
			// (bam2bcf_edlib.c:1217-1218).
			tmp := float64((sc[0] >> 6) & 0xff)
			if tmp > cgpTMPMagic {
				indelQ = 0
				indelQ1 = 0
			} else {
				indelQ = int((1.0-tmp/cgpTMPMagic)*float64(indelQ) + 0.499)
				indelQ1 = int((1.0-tmp/cgpTMPMagic)*float64(indelQ1) + 0.499)
			}

			if indelQ > 255 {
				indelQ = 255
			}
			if indelQ1 > 255 {
				indelQ1 = 255
			}
			if indelQ > seqQ {
				indelQ = seqQ
			}
			if indelQ > 255 {
				indelQ = 255
			}
			if indelQ1 > 255 {
				indelQ1 = 255
			}
			if seqQ > 255 {
				seqQ = 255
			}
			if indelQ < 0 {
				indelQ = 0
			}
			if indelQ1 < 0 {
				indelQ1 = 0
			}
			if seqQ < 0 {
				seqQ = 0
			}

			chosen := sc[0] & 0x3f
			p.aux = uint32(chosen)<<16 | uint32(seqQ)<<8 | uint32(indelQ)
			sumq[chosen] += indelQ
		}
	}

	// Sort sumq descending; promote REF to slot 0
	// (bam2bcf_edlib.c:1257-1269).
	bca.MaxIns = maxIns
	for t := 0; t < nTypes; t++ {
		sumq[t] = sumq[t]<<6 | t
	}
	for t := 1; t < nTypes; t++ {
		for j := t; j > 0 && sumq[j] > sumq[j-1]; j-- {
			sumq[j], sumq[j-1] = sumq[j-1], sumq[j]
		}
	}
	t := 0
	for ; t < nTypes; t++ {
		if sumq[t]&0x3f == refType {
			break
		}
	}
	if t > 0 && t < nTypes {
		tmp := sumq[t]
		for ; t > 0; t-- {
			sumq[t] = sumq[t-1]
		}
		sumq[0] = tmp
	}

	for j := 0; j < 4; j++ {
		bca.IndelTypes[j] = b2bIndelNull
	}
	if maxIns > 0 {
		bca.Inscns = make([]byte, 4*maxIns)
	} else {
		bca.Inscns = nil
	}
	for j := 0; j < 4 && j < nTypes; j++ {
		bca.IndelTypes[j] = types[sumq[j]&0x3f]
		if maxIns > 0 {
			copy(bca.Inscns[j*maxIns:(j+1)*maxIns],
				inscns[(sumq[j]&0x3f)*maxIns:(sumq[j]&0x3f+1)*maxIns])
		}
	}

	// Re-key per-read chosen type (bam2bcf_edlib.c:1285-1299).
	nAlt := 0
	for s := 0; s < len(piles); s++ {
		for i := 0; i < len(piles[s]); i++ {
			p := &piles[s][i]
			x := types[int(p.aux>>16)&0x3f]
			j := 0
			for ; j < 4; j++ {
				if x == bca.IndelTypes[j] {
					break
				}
			}
			if j == 4 {
				p.aux = 4 << 16
			} else {
				p.aux = uint32(j)<<16 | (p.aux & 0xffff)
			}
			if int(p.aux>>16)&0x3f > 0 {
				nAlt++
			}
		}
	}
	return nAlt
}

// bcfCallGapPrepCNS is the Go port of bam2bcf_edlib.c:1319
// bcf_edlib_gap_prep. See file header for the algorithm overview.
func bcfCallGapPrepCNS(piles [][]pileupBase, pos int, bca *bcfCallauxIndel, ref []byte) int {
	if bca == nil || ref == nil {
		return -1
	}

	anyIndel := false
	for _, pile := range piles {
		for _, p := range pile {
			if p.indel != 0 {
				anyIndel = true
				break
			}
		}
		if anyIndel {
			break
		}
	}
	if !anyIndel {
		return -1
	}

	// Average base quality (bam2bcf_edlib.c:1339-1359).
	qavg := 30.0
	qsum := 0.0
	qcount := 0.0
	const qwin = 50
	for _, pile := range piles {
		for _, p := range pile {
			if p.rec == nil {
				continue
			}
			qual := p.rec.Qual
			kstart := p.qpos - qwin
			if kstart < 0 {
				kstart = 0
			}
			kend := p.qpos + qwin
			if kend > len(qual) {
				kend = len(qual)
			}
			for k := kstart; k < kend; k++ {
				qsum += float64(qual[k])
				qcount++
			}
		}
	}
	qavg = (qsum + 1) / (qcount + 1)

	tr := bcfCgpFindTypes(piles, pos, bca, ref)
	if tr == nil {
		return -1
	}
	types := tr.Types
	nTypes := len(types)
	N := tr.N
	refType := tr.RefType

	// Window: max_indel = 20*max(|types[0]|,|types[-1]|) + win/4
	// (bam2bcf_edlib.c:1370-1378).
	at0 := types[0]
	if at0 < 0 {
		at0 = -at0
	}
	atN := types[nTypes-1]
	if atN < 0 {
		atN = -atN
	}
	mx := at0
	if atN > mx {
		mx = atN
	}
	maxIndel := 20*mx + bca.IndelWinSize/4
	if maxIndel > bca.IndelWinSize {
		maxIndel = bca.IndelWinSize
	}
	left := pos - maxIndel
	if left < 0 {
		left = 0
	}
	right := pos + maxIndel
	if types[0] < 0 {
		right += -types[0]
	}
	if right > len(ref) {
		right = len(ref)
	}

	lRun := bcfCgpLRun(ref, pos)
	var lRunBase byte
	if pos+1 < len(ref) {
		lRunBase = byte(nt16ASCII(ref[pos+1]))
	}
	lRunIns := byte(0)

	maxIns := types[nTypes-1]
	if maxIns < 0 {
		maxIns = 0
	}

	// Reset bias histograms.
	if len(bca.IrefPos) != b2bNpos {
		bca.IrefPos = make([]int, b2bNpos)
		bca.IaltPos = make([]int, b2bNpos)
	} else {
		for i := range bca.IrefPos {
			bca.IrefPos[i] = 0
			bca.IaltPos[i] = 0
		}
	}
	if len(bca.IrefMq) != b2bNqual {
		bca.IrefMq = make([]int, b2bNqual)
		bca.IaltMq = make([]int, b2bNqual)
	} else {
		for i := range bca.IrefMq {
			bca.IrefMq[i] = 0
			bca.IaltMq[i] = 0
		}
	}
	for i := range bca.IrefScl {
		bca.IrefScl[i] = 0
		bca.IaltScl[i] = 0
	}

	var inscns []byte
	if maxIns > 0 {
		inscns = bcfCgpCalcCons(piles, types, maxIns)
		if inscns == nil {
			return -1
		}
	}

	biggestDel := 0
	biggestIns := 0
	for t := 0; t < nTypes; t++ {
		if biggestDel > types[t] {
			biggestDel = types[t]
		}
		if biggestIns < types[t] {
			biggestIns = types[t]
		}
	}
	band := biggestIns - biggestDel

	// STR span around pos (bam2bcf_edlib.c:1419-1437).
	posL, posR := pos, pos
	{
		pstart := pos - 30
		if pstart < 0 {
			pstart = 0
		}
		pend := pos + 30
		if pend > len(ref) {
			pend = len(ref)
		}
		seg := make([]byte, pend-pstart)
		for i := range seg {
			seg[i] = byte(cgpBase6(ref[pstart+i]))
		}
		pmid := pos - pstart
		for _, elt := range strfinder.FindSTR(seg, false) {
			if elt.End >= pmid && elt.Start <= pmid {
				if posL > pstart+elt.Start {
					posL = pstart + elt.Start
				}
				if posR < pstart+elt.End {
					posR = pstart + elt.End
				}
			}
		}
	}

	scores := make([]int, N*nTypes)
	bca.IndelReg = 0

	for t := 0; t < nTypes; t++ {
		var ir int
		switch {
		case types[t] == 0:
			ir = 0
		case types[t] > 0:
			ir = estIndelreg(pos, ref, types[t], inscns[t*maxIns:(t+1)*maxIns])
		default:
			ir = estIndelreg(pos, ref, -types[t], nil)
		}
		if ir > bca.IndelReg {
			bca.IndelReg = ir
		}

		K := 0
		for s := 0; s < len(piles); s++ {
			cr := bcfCgpConsensus(piles[s], pos, bca, ref, left, right,
				types[t], -biggestDel, posL, posR)
			tcons := cr.cons
			tconLen := cr.tconLen
			leftShift := cr.leftShift
			rightShift := cr.rightShift
			if cr.band > band {
				band = cr.band
			}

			// Scan for base-runs in the insertion
			// (bam2bcf_edlib.c:1515-1522).
			if cr.cposPos < tconLen[0] {
				kBase := tcons[0][cr.cposPos]
				j := 0
				for ; j < types[t] && cr.cposPos+j < tconLen[0]; j++ {
					if tcons[0][cr.cposPos+j] != kBase {
						break
					}
				}
				if j > 0 && j == types[t] {
					switch kBase {
					case 0:
						lRunIns |= 0x1
					case 1:
						lRunIns |= 0x2
					case 2:
						lRunIns |= 0x4
					case 3:
						lRunIns |= 0x8
					default:
						lRunIns |= 0xf
					}
				}
				if types[t] < 0 {
					lRunIns |= 0xff
				}
			}

			for i := 0; i < len(piles[s]); i, K = i+1, K+1 {
				p := &piles[s][i]

				if t == 0 && p.rec != nil {
					imq := int(p.rec.MapQ)
					if imq > 59 {
						imq = 59
					}
					gp := getPos(p.rec.Cigar, int(p.rec.Pos)-1, p.qpos, p.qlen, p.indel)
					epos, scLen := gp.EPos, gp.ScLen
					if epos < 0 {
						epos = 0
					}
					if epos >= b2bNpos {
						epos = b2bNpos - 1
					}
					if scLen < 0 {
						scLen = 0
					}
					if scLen >= 100 {
						scLen = 99
					}
					if p.indel != 0 {
						bca.IaltMq[imq]++
						bca.IaltScl[scLen]++
						bca.IaltPos[epos]++
					} else {
						bca.IrefMq[imq]++
						bca.IrefScl[scLen]++
						bca.IrefPos[epos]++
					}
				}
				if p.rec == nil {
					continue
				}
				if p.rec.Flag&sam.FlagUnmapped != 0 {
					continue
				}

				skip := false
				for _, op := range p.rec.Cigar {
					if op.Op() == sam.CigarSkipped {
						skip = true
						break
					}
				}
				if skip {
					continue
				}

				// Long-read window-trim (bam2bcf_edlib.c:1576-1608).
				left2, right2 := left, right
				absLS := leftShift
				if absLS < 0 {
					absLS = -absLS
				}
				absRS := rightShift
				if absRS < 0 {
					absRS = -absRS
				}
				minWinSize := biggestIns
				if -biggestDel > minWinSize {
					minWinSize = -biggestDel
				}
				minWinSize += absLS + absRS
				totStr := 0
				for _, elt := range strfinder.FindSTR(tcons[0][:tconLen[0]], false) {
					totStr += elt.End - elt.Start
				}
				minWinSize += totStr + 10
				if len(p.rec.Seq) > 1000 {
					if pos-left >= minWinSize {
						if pos-minWinSize > left2 {
							left2 = pos - minWinSize
						}
					}
					if right-pos >= minWinSize {
						if pos+minWinSize < right2 {
							right2 = pos + minWinSize
						}
					}
				}

				rStart := int(p.rec.Pos) - 1

				var tbeg, tend int
				qbeg := tpos2qpos(p.rec.Cigar, rStart, left2, false, &tbeg)
				_ = tpos2qpos(p.rec.Cigar, rStart, pos, false, &tend)
				qend := tpos2qpos(p.rec.Cigar, rStart, right2, true, &tend)
				oldTend := tend
				oldTbeg := tbeg

				qlen := qend - qbeg
				if qlen <= 0 {
					scores[K*nTypes+t] = b2bIndelMissingScore
					continue
				}
				queryBuf := make([]byte, qlen)
				for l := qbeg; l < qend; l++ {
					if l < 0 || l >= len(p.rec.Seq) {
						queryBuf[l-qbeg] = 4
						continue
					}
					queryBuf[l-qbeg] = byte(seqNt16ToBase6(nt16ASCII(p.rec.Seq[l])))
				}

				// Trim tbeg/tend (bam2bcf_edlib.c:1647-1652).
				wbandMax := -biggestDel
				if biggestIns > wbandMax {
					wbandMax = biggestIns
				}
				wband := band + wbandMax*2 + 20
				tend1 := left + tconLen[0] - (left2 - left)
				tend2 := left + tconLen[1] - (left2 - left)
				if tend1 > oldTend+wband {
					tend1 = oldTend + wband
				}
				if tend2 > oldTend+wband {
					tend2 = oldTend + wband
				}
				if oldTbeg-wband > left2 {
					tbeg = oldTbeg - wband
				} else {
					tbeg = left2
				}

				if !(tend1 > tbeg && tend2 > tbeg) {
					scores[K*nTypes+t] = b2bIndelMissingScore
					continue
				}

				r1Start := tbeg - left
				if r1Start < 0 {
					r1Start = 0
				}
				if r1Start > tconLen[0] {
					r1Start = tconLen[0]
				}
				r1Len := tend1 - tbeg
				if r1Start+r1Len > tconLen[0] {
					r1Len = tconLen[0] - r1Start
				}
				if r1Len < 0 {
					r1Len = 0
				}
				r2Start := tbeg - left
				if r2Start < 0 {
					r2Start = 0
				}
				if r2Start > tconLen[1] {
					r2Start = tconLen[1]
				}
				r2Len := tend2 - tbeg
				if r2Start+r2Len > tconLen[1] {
					r2Len = tconLen[1] - r2Start
				}
				if r2Len < 0 {
					r2Len = 0
				}

				atype := types[t]
				if atype < 0 {
					atype = -atype
				}
				score := bcfCgpAlignScoreCNS(
					tcons[0][r1Start:r1Start+r1Len],
					tcons[1][r2Start:r2Start+r2Len],
					queryBuf,
					atype,
					bca.IndelBias, bca.DelBias,
				)
				scores[K*nTypes+t] = score
			}
		}
	}

	// Invalidate l_run if the candidate insertion's base type is not
	// represented in the surrounding homopolymer (bam2bcf_edlib.c:
	// 1690-1691).
	var lRunMask byte
	switch lRunBase {
	case 1: // A
		lRunMask = 0x1
	case 2: // C
		lRunMask = 0x2
	case 4: // G
		lRunMask = 0x4
	case 8: // T
		lRunMask = 0x8
	case 15: // N
		lRunMask = 0xf
	}
	if lRunMask&lRunIns == 0 {
		lRun = 1
	}

	nAlt := bcfCgpComputeIndelQCNS(piles, scores, bca, inscns, lRun, maxIns,
		refType, types, qavg)
	if nAlt > 0 {
		return nAlt
	}
	return -1
}
