// Per-haplotype alignment scoring and gap-prep orchestration for the
// `bcftools mpileup` indel caller (sub-slices 4c + 4d).
//
// This file ports the algorithmic core of bam2bcf_indel.c:
//
//   - bcfCgpRefSample   (~bam2bcf_indel.c:282 bcf_cgp_ref_sample) —
//     samples the reference window into a per-sample byte buffer of
//     4-bit IUPAC codes and masks REF positions where ≥70% of the reads
//     disagree (a poor-man's IUPAC code, marked as 15 = N).
//   - bcfCgpCalcCons    (~bam2bcf_indel.c:436 bcf_cgp_calc_cons) —
//     builds the per-type insertion consensus by majority rule across
//     all reads carrying that exact insertion length, in 2-bit nt codes
//     (0=A,1=C,2=G,3=T,4=N). An insertion whose consensus contains an
//     N is dropped (its type[t] is zeroed) so it stops being a
//     candidate downstream.
//   - bcfCgpAlignScore  (~bam2bcf_indel.c:497 bcf_cgp_align_score) —
//     aligns one read against a candidate haplotype using
//     baq.ProbalnGlocal, applies the indel-bias clamp to the
//     length-normalised score, then folds in the STR-finder fudge.
//     Reproduces upstream's `(sc << 8) | min(255, l)` bit-pattern
//     byte-for-byte.
//   - bcfCallGapPrep    (~bam2bcf_indel.c:698 bcf_call_gap_prep) — the
//     orchestrator: discovers indel types, samples a per-sample
//     reference, builds per-type insertion consensus, then walks every
//     read against every candidate haplotype to fill a score matrix
//     before handing off to bcfCgpComputeIndelQ.
//   - bcfCgpComputeIndelQ (~bam2bcf_indel.c:597 bcf_cgp_compute_indelQ)
//     folds the per-haplotype scores into per-read p.aux words and
//     populates bca.IndelTypes / bca.Inscns / bca.MaxIns.
//
// No emission: integration into bcfCallGlfgen / bcfCall2bcf and the
// per-site INDEL records lands in sub-slice 4e.
//
// Upstream reference: reference_code/bcftools/bam2bcf_indel.c.

package bcftools

import (
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/baq"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/strfinder"
)

// b2bIndelMissingScore is the sentinel value upstream stores when a
// read's alignment-score computation fails (probaln_glocal returns a
// negative value, bam2bcf_indel.c:540-541, 910-913) or when the read
// covers the candidate region entirely within a deletion so that
// tend <= tbeg (bam2bcf_indel.c:914-918). It matches the 0xffffff cap
// used by upstream. Note: upstream does NOT use this sentinel for
// unmapped reads, CREF_SKIP reads, or reads with p->b == NULL; those
// `continue` and leave the calloc'd score at 0.
const b2bIndelMissingScore = 0xffffff

// b2bIndelNull is upstream's B2B_INDEL_NULL sentinel for an absent indel
// type slot (bam2bcf.h:44, used at bam2bcf_indel.c:660). Upstream defines
// it as the literal value 10000 (an out-of-domain indel length); match
// upstream exactly.
const b2bIndelNull = 10000

// nt16ASCII maps an ASCII nucleotide byte (case-insensitive) to its
// 4-bit IUPAC code, mirroring htslib's seq_nt16_table for the subset of
// characters that appear in reference / read sequences. The seq_nt16_int
// mapping (seqNt16Int, in bam2bcf.go) then collapses these 4-bit codes
// to 2-bit nt indices.
func nt16ASCII(b byte) int {
	switch b {
	case 'A', 'a':
		return 1
	case 'C', 'c':
		return 2
	case 'M', 'm':
		return 3
	case 'G', 'g':
		return 4
	case 'R', 'r':
		return 5
	case 'S', 's':
		return 6
	case 'V', 'v':
		return 7
	case 'T', 't':
		return 8
	case 'W', 'w':
		return 9
	case 'Y', 'y':
		return 10
	case 'H', 'h':
		return 11
	case 'K', 'k':
		return 12
	case 'D', 'd':
		return 13
	case 'B', 'b':
		return 14
	}
	return 15
}

// bcfCgpRefSample is the Go port of bam2bcf_indel.c:282 bcf_cgp_ref_sample.
//
// It returns a per-sample byte buffer of length right-left (in 4-bit
// IUPAC codes) where positions at which the deepest ALT-supporting
// fraction reaches ≥70% are replaced with N (code 15) so the aligner
// neither rewards nor penalises that position. left and right are
// inclusive/exclusive in the same way upstream uses them: ref0
// occupies indices 0..right-left-1 of the returned buffer.
//
// piles is the per-sample pile (one slice per sample). ref is the full
// reference 0-indexed in ASCII.
func bcfCgpRefSample(piles [][]pileupBase, ref []byte, left, right int) [][]byte {
	n := len(piles)
	L := right - left
	if L <= 0 {
		out := make([][]byte, n)
		for s := range out {
			out[s] = []byte{}
		}
		return out
	}

	// Build the 4-bit ref0 once.
	ref0 := make([]byte, L)
	for i := 0; i < L; i++ {
		if i+left < 0 || i+left >= len(ref) {
			ref0[i] = 15
			continue
		}
		ref0[i] = byte(nt16ASCII(ref[i+left]))
	}

	out := make([][]byte, n)
	for s := 0; s < n; s++ {
		// Reset the cns counters per sample. Layout matches upstream:
		// low 16 bits hold the REF count, high 16 bits hold the ALT
		// count (bam2bcf_indel.c:349-351).
		cns := make([]uint32, L)

		for _, p := range piles[s] {
			if p.rec == nil {
				continue
			}
			rec := p.rec
			x := int(rec.Pos) - 1 // 0-based reference start
			y := 0
			for _, op := range rec.Cigar {
				o := op.Op()
				l := int(op.Length())
				switch o {
				case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
					if x+l >= left {
						j := 0
						if left-x > 0 {
							j = left - x
						}
						jEnd := l
						if right-x < l {
							jEnd = right - x
						}
						for ; j < jEnd; j++ {
							seqIdx := y + j
							if seqIdx < 0 || seqIdx >= len(rec.Seq) {
								continue
							}
							idx := x + j - left
							if idx < 0 || idx >= L {
								continue
							}
							b4 := byte(nt16ASCII(rec.Seq[seqIdx]))
							if b4 == ref0[idx] {
								cns[idx] += 1
							} else {
								cns[idx] += 1 << 16
							}
						}
					}
					x += l
					y += l
				case sam.CigarDeletion, sam.CigarSkipped:
					x += l
				case sam.CigarInsertion, sam.CigarSoftClip:
					y += l
				}
				if x > right {
					break
				}
			}
		}

		// Start the sample-specific reference as a copy of ref0.
		r := make([]byte, L)
		copy(r, ref0)

		// Find deepest and second-deepest ALT positions
		// (bam2bcf_indel.c:370-376). The comparison uses the high 16
		// bits as the ALT count.
		var max, max2 uint32
		maxI, max2I := -1, -1
		for i := 0; i < L; i++ {
			altCnt := cns[i] >> 16
			if altCnt >= max>>16 {
				max2 = max
				max2I = maxI
				max = cns[i]
				maxI = i
			} else if altCnt >= max2>>16 {
				max2 = cns[i]
				max2I = i
			}
		}
		// Mask the two deepest sites with N (code 15) only when ALT
		// >= 30% of the coverage there (i.e. REF/(REF+ALT) < 0.7).
		// Upstream's predicate (bam2bcf_indel.c:385-388):
		//   (double)(max&0xffff) / ((max&0xffff) + (max>>16)) >= 0.7
		//     ? max_i = -1 : /* keep max_i */
		// When total == 0 the division is NaN and `NaN >= 0.7` is false
		// in both C and Go, so max_i is NOT cleared and r[max_i] = 15
		// still runs. Follow upstream literally (no total>0 guard) so
		// the all-zero-coverage branch matches byte-for-byte.
		if maxI >= 0 {
			ref := max & 0xffff
			alt := max >> 16
			total := ref + alt
			if float64(ref)/float64(total) >= 0.7 {
				maxI = -1
			}
		}
		if max2I >= 0 {
			ref := max2 & 0xffff
			alt := max2 >> 16
			total := ref + alt
			if float64(ref)/float64(total) >= 0.7 {
				max2I = -1
			}
		}
		if maxI >= 0 {
			r[maxI] = 15
		}
		if max2I >= 0 {
			r[max2I] = 15
		}
		out[s] = r
	}
	return out
}

// bcfCgpCalcCons is the Go port of bam2bcf_indel.c:436 bcf_cgp_calc_cons.
//
// It builds a per-type insertion consensus by majority rule across all
// reads carrying that exact insertion length, returning a flat
// n_types * max_ins byte slice in 2-bit nt codes (0=A,1=C,2=G,3=T,4=N).
// Insertion types whose consensus contains an N are zeroed in the
// caller-owned types slice (matching upstream, which also rewrites
// types[t] = 0).
func bcfCgpCalcCons(piles [][]pileupBase, types []int, maxIns int) []byte {
	nTypes := len(types)
	if nTypes == 0 || maxIns <= 0 {
		return nil
	}
	// inscnsAux: counts of each base (0..4) at each (type, position)
	// slot. Indexed [(t*maxIns + j)*5 + c].
	inscnsAux := make([]int, 5*nTypes*maxIns)

	for t, tLen := range types {
		if tLen <= 0 {
			continue
		}
		for _, pile := range piles {
			for _, p := range pile {
				if p.indel != tLen || p.rec == nil {
					continue
				}
				seq := p.rec.Seq
				for k := 1; k <= tLen; k++ {
					if p.qpos+k < 0 || p.qpos+k >= len(seq) {
						continue
					}
					c := seqNt16Int[nt16ASCII(seq[p.qpos+k])]
					if c < 0 || c >= 5 {
						c = 4
					}
					inscnsAux[(t*maxIns+(k-1))*5+c]++
				}
			}
		}
	}

	inscns := make([]byte, nTypes*maxIns)
	for t := 0; t < nTypes; t++ {
		tLen := types[t]
		if tLen <= 0 {
			continue
		}
		for j := 0; j < tLen; j++ {
			max := 0
			maxK := -1
			off := (t*maxIns + j) * 5
			for k := 0; k < 5; k++ {
				if inscnsAux[off+k] > max {
					max = inscnsAux[off+k]
					maxK = k
				}
			}
			if maxK < 0 {
				inscns[t*maxIns+j] = 4
			} else {
				inscns[t*maxIns+j] = byte(maxK)
			}
			if maxK == 4 {
				// Discard insertions containing N's.
				types[t] = 0
				break
			}
		}
	}
	return inscns
}

// bcfCgpAlignScore is the Go port of bam2bcf_indel.c:497 bcf_cgp_align_score.
//
// It aligns one read's query window against a candidate haplotype using
// baq.ProbalnGlocal, then folds in the STR-finder fudge. The returned
// score word reproduces upstream's bit-pattern exactly:
//
//	score = (sc << 8) | min(255, l)
//
// where sc is the raw Phred-scaled probaln_glocal return and l is the
// length-normalised score after the indel-bias clamp. The STR fudge
// then rewrites the low byte with min(255, (score&0xff)*0.8 + iscore*2).
//
// On probaln failure (a negative return) the score is clamped to
// 0xffffff, matching upstream.
//
// ref2 holds the haplotype in 2-bit nt codes (0..4); query holds the
// read window likewise. qq holds the read window's clamped per-base
// qualities (already in [7,30]). tEnd-tBeg is the haplotype "real"
// reference span; "+ typeLen" widens it to accommodate the insertion
// payload, mirroring upstream's `tend - tbeg + type` length argument.
//
// indelBias is bca.IndelBias (a float multiplier on l); typeLen is
// abs(types[t]). longRead toggles the PacBio CCS parameter set
// (bam2bcf_indel.c:506-514). qpos is the read's query position
// (post-qbeg shift) of the candidate indel; rStart and rEnd are the
// read's reference span (0-based, inclusive) — used by the STR fudge to
// detect repeats that extend to the read edge.
func bcfCgpAlignScore(ref2 []byte, query []byte, qq []byte,
	typeLen int, longRead bool, indelBias float64,
	qpos, tBeg, tEnd, rStart, rEnd int) int {

	if typeLen < 0 {
		typeLen = -typeLen
	}

	par := baq.Par{D: 1e-4, E: 1e-2, BW: typeLen + 3}
	if longRead {
		par.D = 1e-3
		par.E = 1e-1
	}

	// Upstream (bam2bcf_indel.c:538-539) passes an explicit reference
	// length of `tend - tbeg + type` to probaln_glocal:
	//   sc = probaln_glocal(ref2 + tbeg - left, tend - tbeg + type,
	//                       query, qend - qbeg, qq, &apf, 0, 0);
	// Our ProbalnGlocal derives lRef from len(ref), and the orchestrator
	// hands us a slice that is typically longer than `tend - tbeg + type`
	// (it covers the full padded haplotype). Truncate to the upstream
	// length so sc matches byte-for-byte.
	tlen := (tEnd - tBeg) + typeLen
	if tlen > len(ref2) {
		tlen = len(ref2)
	}
	if tlen < 0 {
		tlen = 0
	}
	sc, err := baq.ProbalnGlocal(ref2[:tlen], query, qq, par, nil, nil)
	if err != nil || sc < 0 {
		return b2bIndelMissingScore
	}

	// Upstream's adjustment of indelQ via length-normalised score:
	//   l = (int)(100. * sc / (qend - qbeg) + .499) * bca->indel_bias
	//   score = sc<<8 | MIN(255, l)
	// In C this is `int l; l = (int)(...) * bca->indel_bias;` which
	// promotes l to double for the multiplication then truncates back
	// to int on assignment. Reproduce that exactly.
	qlen := len(query)
	var lNorm int
	if qlen > 0 {
		base := int(100.0*float64(sc)/float64(qlen) + 0.499)
		lNorm = int(float64(base) * indelBias)
	}
	if lNorm > 255 {
		lNorm = 255
	}
	if lNorm < 0 {
		lNorm = 0
	}
	score := (sc << 8) | lNorm

	// STR fudge (bam2bcf_indel.c:550-586). seg = ref2[tBeg-left ..],
	// but here we have already been handed the relevant slice — the
	// caller positions ref2 so that the segment of interest starts at
	// index 0. seg_len = (tEnd - tBeg) + typeLen.
	segLen := (tEnd - tBeg) + typeLen
	if segLen > len(ref2) {
		segLen = len(ref2)
	}
	if segLen < 0 {
		segLen = 0
	}
	reps := strfinder.FindSTR(ref2[:segLen], false)
	iscore := 0
	for _, r := range reps {
		// qpos here is the read's qpos relative to qbeg (i.e. the
		// candidate's query position into the window). Upstream's
		// condition compares against qpos in seg coordinates; this is
		// the same comparison since both `qpos` and `r.Start/End` are
		// measured into the same window. Upstream lines:
		//   if (elt->start <= qpos && elt->end >= qpos) { ... }
		//   if (elt->start+tbeg <= r_start ||
		//       elt->end+tbeg   >= r_end) iscore += 2*(elt->end-elt->start);
		if r.Start <= qpos && r.End >= qpos {
			repLen := r.RepLen
			if repLen <= 0 {
				repLen = 1
			}
			iscore += (r.End - r.Start) / repLen
			if r.Start+tBeg <= rStart || r.End+tBeg >= rEnd {
				iscore += 2 * (r.End - r.Start)
			}
		}
	}

	// score = (score & ~0xff) | MIN(255, (score&0xff)*.8 + iscore*2)
	l2 := int(float64(score&0xff)*0.8 + float64(iscore*2))
	if l2 > 255 {
		l2 = 255
	}
	if l2 < 0 {
		l2 = 0
	}
	score = (score & ^0xff) | l2
	return score
}

// bcfCgpComputeIndelQ is the Go port of bam2bcf_indel.c:597
// bcf_cgp_compute_indelQ.
//
// scores is laid out as N*nTypes ints (one per read per type) matching
// upstream's `score[K*n_types + t]`. piles is the per-sample pile.
// inscns is the per-type insertion consensus from bcfCgpCalcCons.
// l_run is bcfCgpLRun. ref_type is the index into types[] of the REF
// (0-length) entry. The function fills bca.IndelTypes, bca.Inscns,
// bca.MaxIns and writes each read's p.aux word; it returns the number
// of reads whose chosen type is non-REF (n_alt).
func bcfCgpComputeIndelQ(piles [][]pileupBase, scores []int,
	bca *bcfCallauxIndel, inscns []byte, lRun, maxIns, refType int,
	types []int) int {
	nTypes := len(types)
	sc := make([]int, nTypes)
	sumq := make([]int, nTypes)

	K := 0
	for s := 0; s < len(piles); s++ {
		for i := 0; i < len(piles[s]); i, K = i+1, K+1 {
			p := &piles[s][i]
			// sct = &score[K*n_types]
			off := K * nTypes
			for t := 0; t < nTypes; t++ {
				sc[t] = scores[off+t]<<6 | t
			}
			// Insertion sort ascending.
			for t := 1; t < nTypes; t++ {
				for j := t; j > 0 && sc[j] < sc[j-1]; j-- {
					sc[j], sc[j-1] = sc[j-1], sc[j]
				}
			}

			var indelQ, seqQ int
			if sc[0]&0x3f == refType {
				indelQ = (sc[1] >> 14) - (sc[0] >> 14)
				seqQ = estSeqQ(bca, types[sc[1]&0x3f], lRun)
			} else {
				// Find REF position.
				t := 0
				for ; t < nTypes; t++ {
					if sc[t]&0x3f == refType {
						break
					}
				}
				if t == nTypes {
					t = 0
				}
				indelQ = (sc[t] >> 14) - (sc[0] >> 14)
				seqQ = estSeqQ(bca, types[sc[0]&0x3f], lRun)
			}
			tmp := (sc[0] >> 6) & 0xff
			if tmp > 111 {
				indelQ = 0
			} else {
				indelQ = int((1.0-float64(tmp)/111.0)*float64(indelQ) + 0.499)
			}
			if indelQ > seqQ {
				indelQ = seqQ
			}
			if indelQ > 255 {
				indelQ = 255
			}
			if indelQ < 0 {
				indelQ = 0
			}
			if seqQ > 255 {
				seqQ = 255
			}
			if seqQ < 0 {
				seqQ = 0
			}
			chosen := sc[0] & 0x3f
			p.aux = uint32(chosen)<<16 | uint32(seqQ)<<8 | uint32(indelQ)
			minQ := indelQ
			if seqQ < minQ {
				minQ = seqQ
			}
			sumq[chosen] += minQ
		}
	}

	// Determine bca.IndelTypes and bca.Inscns: sort sumq descending
	// while encoding the original type index in the low 6 bits, then
	// move the REF type to position 0 (upstream bam2bcf_indel.c:653-659).
	bca.MaxIns = maxIns
	for t := 0; t < nTypes; t++ {
		sumq[t] = sumq[t]<<6 | t
	}
	// Insertion sort descending.
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

	// Update p->aux: re-key the chosen type from "index into types"
	// to "index into bca.IndelTypes". Reads whose chosen type is no
	// longer in the top-4 get aux = 4<<16 (so the SNP-only inspection
	// of (aux>>16) sees a sentinel and skips them, matching upstream).
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

// bcfCallGapPrep is the Go port of bam2bcf_indel.c:698 bcf_call_gap_prep.
//
// It orchestrates indel calling at one reference position pos:
//   - bails out (-1) when no read carries an indel at this column;
//   - discovers indel types via bcfCgpFindTypes;
//   - samples a per-sample reference (bcfCgpRefSample);
//   - builds the per-type insertion consensus (bcfCgpCalcCons);
//   - for each (read, type) pair, builds the haplotype, extracts the
//     read window and calls bcfCgpAlignScore;
//   - folds the per-read scores into per-allele indel qualities via
//     bcfCgpComputeIndelQ.
//
// It populates bca.IndelTypes, bca.Inscns, bca.MaxIns and each read's
// p.aux word, and returns the number of non-REF reads (n_alt), or -1
// when the site is not a candidate.
func bcfCallGapPrep(piles [][]pileupBase, pos int, bca *bcfCallauxIndel, ref []byte) int {
	if bca == nil || ref == nil {
		return -1
	}

	// Cheap reject: no indels in any pile.
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

	// Left/right window.
	left := pos - bca.IndelWinSize
	if left < 0 {
		left = 0
	}
	right := pos + bca.IndelWinSize
	if types[0] < 0 {
		right -= types[0] // types[0] is negative -> widens right
	}
	if right > len(ref) {
		right = len(ref)
	}
	// Trim right at a '\0'-equivalent end (we use len(ref) as the
	// natural end — upstream's `ref[i] == 0` check is moot here).

	refSample := bcfCgpRefSample(piles, ref, left, right)
	lRun := bcfCgpLRun(ref, pos)

	maxIns := types[nTypes-1]
	if maxIns < 0 {
		maxIns = 0
	}
	var inscns []byte
	if maxIns > 0 {
		inscns = bcfCgpCalcCons(piles, types, maxIns)
		if inscns == nil {
			return -1
		}
	}

	// Allocate score matrix and scratch buffers.
	scores := make([]int, N*nTypes)
	// max_ref2 covers the longest haplotype (ref window + room for
	// either a max-length insertion or deletion). +2 for the closed/open
	// boundaries.
	maxDelExt := -types[0]
	if maxDelExt < 0 {
		maxDelExt = 0
	}
	pad := maxIns
	if maxDelExt > pad {
		pad = maxDelExt
	}
	maxRef2 := (right - left) + 2 + 2*pad
	ref2 := make([]byte, maxRef2)
	queryBuf := make([]byte, (right-left)+maxRdLen+maxIns+2)
	bca.IndelReg = 0

	// Determine max deletion for the column.
	maxDeletion := 0
	for _, pile := range piles {
		for _, p := range pile {
			if -p.indel > maxDeletion {
				maxDeletion = -p.indel
			}
		}
	}
	_ = maxDeletion // (kept for future use; the score path doesn't need it)

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
			// Build ref2 for this (sample, type): left..pos from
			// ref_sample[s], then insertion payload (or skip deletion
			// bases), then pos..right.
			k := 0
			j := left
			for ; j <= pos && j-left < len(refSample[s]); j++ {
				code4 := refSample[s][j-left]
				ref2[k] = byte(seqNt16Int[code4])
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
			for ; j < right && j-left < len(refSample[s]); j++ {
				code4 := refSample[s][j-left]
				ref2[k] = byte(seqNt16Int[code4])
				k++
			}
			for ; k < maxRef2; k++ {
				ref2[k] = 4
			}
			rRight := right
			if right > j {
				rRight = j
			}

			for _, p := range piles[s] {
				if p.rec == nil {
					// Upstream (bam2bcf_indel.c) calloc's `score` and
					// `continue`s on p->b == NULL / BAM_FUNMAP / CREF_SKIP
					// (lines 854 and 858-862), leaving the slot at 0.
					// Writing 0xffffff here would rank these reads as
					// worst-possible in the sc[t] = score<<6 | t sort.
					// Leave the Go zero value (0) so we match upstream.
					K++
					continue
				}
				if p.rec.Flag&sam.FlagUnmapped != 0 {
					// bam2bcf_indel.c:854 — `continue`. Score stays 0.
					K++
					continue
				}

				// Reject reads with CREF_SKIP (N) in their CIGAR.
				skip := false
				for _, op := range p.rec.Cigar {
					if op.Op() == sam.CigarSkipped {
						skip = true
						break
					}
				}
				if skip {
					// bam2bcf_indel.c:858-862 — `continue`. Score stays 0.
					K++
					continue
				}

				// Determine alignment window. Long reads use a
				// half-window heuristic (bam2bcf_indel.c:866-875).
				left2, right2 := left, rRight
				if len(p.rec.Seq) > 1000 {
					if pos-left >= bca.IndelWinSize {
						left2 += bca.IndelWinSize / 2
					}
					if rRight-pos >= bca.IndelWinSize {
						right2 -= bca.IndelWinSize / 2
					}
				}

				rStart := int(p.rec.Pos) - 1
				rEnd := rStart - 1
				for _, op := range p.rec.Cigar {
					switch op.Op() {
					case sam.CigarMatch, sam.CigarDeletion,
						sam.CigarSkipped, sam.CigarEqual, sam.CigarMismatch:
						rEnd += int(op.Length())
					}
				}

				var tbeg, tend int
				qbeg := tpos2qpos(p.rec.Cigar, rStart, left2, false, &tbeg)
				qposLocal := tpos2qpos(p.rec.Cigar, rStart, pos, false, &tend) - qbeg
				qend := tpos2qpos(p.rec.Cigar, rStart, right2, true, &tend)

				if types[t] < 0 {
					l := -types[t]
					if tbeg-l > left {
						tbeg = tbeg - l
					} else {
						tbeg = left
					}
				}

				// Write the query window.
				qlen := qend - qbeg
				if qlen <= 0 || tend <= tbeg {
					scores[K*nTypes+t] = b2bIndelMissingScore
					K++
					continue
				}
				if qlen > len(queryBuf) {
					queryBuf = make([]byte, qlen)
				}
				for l := qbeg; l < qend; l++ {
					if l < 0 || l >= len(p.rec.Seq) {
						queryBuf[l-qbeg] = 4
						continue
					}
					queryBuf[l-qbeg] = byte(seqNt16Int[nt16ASCII(p.rec.Seq[l])])
				}

				// Build qq: clamp [7,30].
				qq := make([]byte, qlen)
				for l := qbeg; l < qend; l++ {
					var qv int
					if l < len(p.rec.Qual) {
						qv = int(p.rec.Qual[l])
					}
					if qv > 30 {
						qv = 30
					}
					if qv < 7 {
						qv = 7
					}
					qq[l-qbeg] = byte(qv)
				}

				longRead := len(p.rec.Seq) > 1000

				// ref2 slice positioned at tbeg-left.
				refStart := tbeg - left
				if refStart < 0 {
					refStart = 0
				}
				if refStart > len(ref2) {
					refStart = len(ref2)
				}
				score := bcfCgpAlignScore(
					ref2[refStart:],
					queryBuf[:qlen],
					qq,
					types[t], longRead, bca.IndelBias,
					qposLocal, tbeg, tend, rStart, rEnd,
				)
				scores[K*nTypes+t] = score
				K++
			}
		}
	}

	nAlt := bcfCgpComputeIndelQ(piles, scores, bca, inscns, lRun, maxIns, refType, types)
	if nAlt > 0 {
		return nAlt
	}
	return -1
}
