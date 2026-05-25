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
//  2. For each (sample, candidate type), build a consensus haplotype
//     spanning [left, right): the reference window with the candidate
//     insertion bases inserted (for +N) or with the deleted span
//     removed (for -N). Upstream's bcf_cgp_consensus does a much
//     fancier per-base / per-insertion majority vote across all reads
//     in the column, producing two consensuses (cnum=0/1) to capture
//     hets. This port uses a simpler "reference + indel payload"
//     haplotype: it is correct in the homozygous case and gives
//     reasonable scores in the het case. The full consensus builder
//     is a follow-up slice.
//  3. For each read in the column, run edlib in EDLIB_MODE_HW
//     (infix) / EDLIB_TASK_DISTANCE against each candidate
//     haplotype. The edit distance (scaled into upstream's
//     `(sc << 8) | normalised` bit layout) becomes that read's score
//     for that type. See edlib_glocal at bam2bcf_edlib.c:834.
//  4. bcfCgpComputeIndelQ folds the per-(read,type) scores into each
//     read's p.aux word (chosen-type | seqQ | indelQ), exactly as the
//     legacy path does. The downstream emission pipeline
//     (bcfCallGlfgenIndel → bcfCallCombineIndel → bcfCall2bcfIndel)
//     is shared with the legacy path: only the score-assignment
//     stage swaps.
//
// Scoping for this slice: we land a functional end-to-end CNS path
// that uses the in-tree edlib engine (pkg/htsgo/edlib) and produces
// non-empty INDEL output. Upstream's full bam2bcf_edlib.c carries a
// significantly elaborated consensus builder
// (bcf_cgp_consensus with cons[0]/cons[1] heterozygous threading)
// and its own variant of compute_indelQ that drives indelQ1/indelQ2,
// vs_ref, poly_mqual and TMP_MAGIC=255. Those refinements are not
// in this slice — they shift residual byte-level scores at
// homopolymer columns. The dispatch wiring, edlib engagement and
// the p.aux/scoring contract are all in place so the next slice
// can land the full upstream consensus + computeIndelQ on top.

package bcftools

import (
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/edlib"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// edlibGlocal is the Go port of bam2bcf_edlib.c:834 edlib_glocal. It
// aligns `query` against `ref` in EDLIB_MODE_HW / EDLIB_TASK_DISTANCE
// and returns a score that combines the raw edit distance with a
// del-bias adjustment based on the aligned target span.
//
// Upstream's formula (bam2bcf_edlib.c:925-927):
//
//	t_len = endLocations[0] - startLocations[0] + 1
//	score = m * (editDistance - delBias * (t_len - lQuery))
//
// m is a constant 30 in upstream (bam2bcf_edlib.c:1015) — kept as a
// parameter here so future tuning is local. delBias is upstream's
// bca->del_bias (CLI --del-bias). On failure (alignment errored or
// returned no locations) the function returns the sentinel score
// b2bIndelMissingScore.
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
	score := int(m*(float64(res.EditDistance)-delBias*float64(tLen-len(query))) + 0.5)
	if score < 0 {
		score = 0
	}
	return score
}

// bcfCgpAlignScoreCNS is the Go port of bam2bcf_edlib.c:971
// bcf_cgp_align_score. It scores one read window against a candidate
// haplotype using edlibGlocal, then packs the result into upstream's
// `(sc << 8) | min(255, l)` bit layout.
//
// Upstream tries two candidate consensuses (the per-sample cons[0]
// and cons[1]) and picks the better of the two. This port hands a
// single haplotype (the simpler "ref + indel payload" consensus) so
// there is no second alignment to take a min over. The scoring math
// otherwise matches upstream byte-for-byte.
func bcfCgpAlignScoreCNS(ref, query []byte, qlen int, indelBias, delBias float64) int {
	const mm = 30.0
	sc := edlibGlocal(ref, query, mm, delBias)
	if sc >= b2bIndelMissingScore {
		return b2bIndelMissingScore
	}

	// Upstream bam2bcf_edlib.c:1053-1055:
	//   l = .5 * (100. * sc / (qend - qbeg) + .499);
	//   score = (sc<<8) | (int)MIN(255, l * bca->indel_bias * .5);
	var lNorm int
	if qlen > 0 {
		base := 0.5 * (100.0*float64(sc)/float64(qlen) + 0.499)
		lNorm = int(base * indelBias * 0.5)
	}
	if lNorm > 255 {
		lNorm = 255
	}
	if lNorm < 0 {
		lNorm = 0
	}
	return (sc << 8) | lNorm
}

// bcfCallGapPrepCNS is the Go port of bam2bcf_edlib.c:1319
// bcf_edlib_gap_prep — the consensus-based analogue of bcfCallGapPrep.
// Returns the number of reads chosen as non-REF (n_alt > 0), or -1
// when the site is not a candidate.
//
// It orchestrates the CNS path at one reference position pos:
//   - bails out (-1) when no read carries an indel at this column;
//   - discovers indel types via bcfCgpFindTypes (shared with legacy);
//   - builds the per-type insertion consensus (bcfCgpCalcCons, shared
//     with legacy — needed for bcfCall2bcfIndel's ALT strings);
//   - for each (read, type) pair, builds the candidate haplotype as
//     `ref[left..pos] + payload(type) + ref[pos..right)` (a simpler
//     analogue of upstream's bcf_cgp_consensus), extracts the read
//     window, and calls bcfCgpAlignScoreCNS;
//   - folds the per-read scores into per-allele indel qualities via
//     the shared bcfCgpComputeIndelQ.
//
// The function populates bca.IndelTypes, bca.Inscns, bca.MaxIns and
// each read's p.aux word — the same contract bcfCallGlfgenIndel /
// bcfCallCombineIndel / bcfCall2bcfIndel consume, so emission code
// is shared with the legacy path.
func bcfCallGapPrepCNS(piles [][]pileupBase, pos int, bca *bcfCallauxIndel, ref []byte) int {
	if bca == nil || ref == nil {
		return -1
	}

	// Cheap reject: no indels in any pile (bam2bcf_edlib.c:1330-1337).
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

	tr := bcfCgpFindTypes(piles, pos, bca, ref)
	if tr == nil {
		return -1
	}
	types := tr.Types
	nTypes := len(types)
	maxRdLen := tr.MaxRdLen
	N := tr.N
	refType := tr.RefType

	// Compute left/right window. Upstream (bam2bcf_edlib.c:1370-1375)
	// uses `max_indel = 20 * max(|types[0]|,|types[-1]|) + win/4`,
	// capped at indel_win_size. We use the full indel_win_size for
	// simplicity — slightly more reference context than upstream, but
	// edlib's HW mode tolerates flanking sequence freely.
	left := pos - bca.IndelWinSize
	if left < 0 {
		left = 0
	}
	right := pos + bca.IndelWinSize
	if types[0] < 0 {
		right -= types[0] // types[0] is negative → widens right
	}
	if right > len(ref) {
		right = len(ref)
	}

	lRun := bcfCgpLRun(ref, pos)

	maxIns := types[nTypes-1]
	if maxIns < 0 {
		maxIns = 0
	}

	// Reset the indel-flavored bias histograms (bam2bcf_edlib.c
	// inherits bam2bcf.h:119 layout; same as the legacy path).
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

	// Allocate score matrix. Layout matches the legacy path:
	// scores[K*nTypes + t] for read index K and type t.
	scores := make([]int, N*nTypes)

	// max_ref2 covers the longest haplotype (ref window + room for
	// either a max-length insertion or deletion).
	maxDelExt := -types[0]
	if maxDelExt < 0 {
		maxDelExt = 0
	}
	pad := maxIns
	if maxDelExt > pad {
		pad = maxDelExt
	}
	maxRef2 := (right - left) + 2 + 2*pad
	if maxRef2 < 1 {
		maxRef2 = 1
	}
	ref2 := make([]byte, maxRef2)
	queryBuf := make([]byte, (right-left)+maxRdLen+maxIns+2)
	bca.IndelReg = 0

	for t := 0; t < nTypes; t++ {
		// Update bca.IndelReg.
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
			// Build ref2 for this type. The haplotype is a 2-bit nt
			// sequence covering [left, right) with the indel applied
			// at pos. Layout:
			//
			//   ref[left .. pos] + payload(t) + ref[pos+1+|del| .. right)
			//
			// where payload(t) is the inserted bases for an insertion
			// of length t, or empty for a deletion (the deleted span
			// is just skipped).
			k := 0
			j := left
			for ; j <= pos && j < len(ref); j++ {
				ref2[k] = byte(seqNt16Int[baseToNt16(upperByte(ref[j]))])
				k++
			}
			if types[t] <= 0 {
				j += -types[t]
			} else {
				for l := 0; l < types[t]; l++ {
					ref2[k] = inscns[t*maxIns+l]
					k++
				}
			}
			for ; j < right && j < len(ref); j++ {
				ref2[k] = byte(seqNt16Int[baseToNt16(upperByte(ref[j]))])
				k++
			}
			tlen := k

			for i := 0; i < len(piles[s]); i, K = i+1, K+1 {
				p := &piles[s][i]

				// Per-read iref/ialt bias tallies (bam2bcf_edlib.c:
				// 1547-1559, mirrors bam2bcf_indel.c:826-849). Only on
				// the t==0 pass since they share storage.
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
					// Leave score at 0 (matches legacy / upstream).
					continue
				}
				if p.rec.Flag&sam.FlagUnmapped != 0 {
					continue
				}

				// Reject reads with CREF_SKIP (N) in their CIGAR
				// (bam2bcf_edlib.c:1569-1573).
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

				// Map left/right genomic coordinates to qbeg/qend
				// query coordinates (bam2bcf_edlib.c:1624-1629).
				rStart := int(p.rec.Pos) - 1
				var tbeg, tend int
				qbeg := tpos2qpos(p.rec.Cigar, rStart, left, false, &tbeg)
				_ = tpos2qpos(p.rec.Cigar, rStart, pos, false, &tend)
				qend := tpos2qpos(p.rec.Cigar, rStart, right, true, &tend)

				qlen := qend - qbeg
				if qlen <= 0 {
					scores[K*nTypes+t] = b2bIndelMissingScore
					continue
				}
				if qlen > len(queryBuf) {
					queryBuf = make([]byte, qlen)
				}
				// Write the query window in 2-bit nt codes
				// (bam2bcf_edlib.c:1635-1636).
				for l := qbeg; l < qend; l++ {
					if l < 0 || l >= len(p.rec.Seq) {
						queryBuf[l-qbeg] = 4
						continue
					}
					queryBuf[l-qbeg] = byte(seqNt16Int[baseToNt16(upperByte(p.rec.Seq[l]))])
				}

				// Slice ref2 to its valid prefix.
				refSlice := ref2[:tlen]
				score := bcfCgpAlignScoreCNS(refSlice, queryBuf[:qlen], qlen,
					bca.IndelBias, bca.DelBias)
				scores[K*nTypes+t] = score
			}
		}
	}

	nAlt := bcfCgpComputeIndelQ(piles, scores, bca, inscns, lRun, maxIns, refType, types)
	if nAlt > 0 {
		return nAlt
	}
	return -1
}
