// Indel-calling state and static helpers for bcftools mpileup.
//
// This file is the foundation (sub-slices 4a + 4b) of the indel-calling
// port of bam2bcf_indel.c. It defines:
//
//   - bcfCallauxIndel: the indel-specific fields that upstream stores on
//     bcf_callaux_t. The struct is exposed at this slice so the CLI knobs
//     can already be parsed into a stable shape; sub-slices 4c+4d
//     populate the per-call scoring state.
//   - The static helpers est_seqQ, est_indelreg, bcf_cgp_l_run,
//     bcf_cgp_find_types, tpos2qpos and get_pos (upstream all `static`
//     in bam2bcf_indel.c). Each is faithful to upstream and unit-tested
//     in bam2bcf_indel_test.go.
//
// The big alignment-scoring core (bcf_cgp_ref_sample, bcf_cgp_calc_cons,
// bcf_cgp_align_score, bcf_call_gap_prep, bcf_cgp_compute_indelQ) is
// deferred to sub-slices 4c+4d, and the indel-aware glfgen / combine /
// 2bcf integration to sub-slice 4e.
//
// Upstream reference: reference_code/bcftools/bam2bcf_indel.c.

package bcftools

import (
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// b2bIndelTypesMax mirrors bam2bcf_indel.c:44 MAX_TYPES — the upper
// bound on the number of distinct indel lengths bcf_cgp_find_types is
// allowed to report. Beyond this the call is dropped.
const b2bIndelTypesMax = 64

// bcfCallauxIndel is the Go counterpart of the indel-specific subset of
// bam2bcf.h's bcf_callaux_t. Fields are named to match upstream as
// closely as Go conventions allow; comments cite upstream line ranges.
type bcfCallauxIndel struct {
	// OpenQ, ExtQ, TandemQ are the indel-quality model parameters
	// (bam2bcf.h:113). openQ is the cost of opening a gap, extQ the
	// per-base extension cost, tandemQ the tandem-repeat penalty.
	OpenQ   int
	ExtQ    int
	TandemQ int

	// MinSupport is the minimum number of indel-supporting reads
	// before a position is considered an indel candidate
	// (bam2bcf.h:114, used at bam2bcf_indel.c:203).
	MinSupport int
	// MaxSupport tracks the observed maximum (bam2bcf.h:114, filled
	// by bcf_cgp_find_types at line 207).
	MaxSupport int
	// MinFrac is the minimum supporting-read fraction (bam2bcf.h:115).
	MinFrac float64
	// MaxFrac tracks the observed maximum (bam2bcf.h:116).
	MaxFrac float64

	// IndelBias adjusts the indel score threshold (bam2bcf.h:133).
	// Lower => call more.
	IndelBias float64

	// IndelWinSize is the window of reference around the candidate
	// site that the indel aligner considers (bam2bcf.h:125).
	IndelWinSize int

	// IndelTypes holds the discovered indel lengths (bam2bcf.h:124).
	// Layout is sorted ascending — negative values are deletions,
	// positive values are insertions, and 0 is REF (always present).
	IndelTypes [4]int

	// MaxIns is the longest insertion length observed; the inscns
	// buffer below sizes against this (bam2bcf.h:127).
	MaxIns int
	// Inscns is the concatenated per-type insertion-consensus
	// sequence buffer (bam2bcf.h:129). Layout: n_types * max_ins bytes,
	// with type t's consensus at inscns[t*max_ins .. t*max_ins+ins_len].
	Inscns []byte

	// IndelReg is the est_indelreg result for the current candidate
	// (bam2bcf.h:127) — the length of the affected reference region.
	IndelReg int

	// PerSampleFlt is the indel filtering strategy (bam2bcf.h:117).
	PerSampleFlt int

	// AmbigReads / Edlib / Indels20 / PolyMQual / SeqQOffset / ReadLen
	// mirror the corresponding bam2bcf.h scalars. They are accepted at
	// the CLI but not yet driven.
	AmbigReads int
	Edlib      int
	Indels20   int
	PolyMQual  int
	SeqQOffset int
	ReadLen    int

	// IrefPos / IaltPos are the indel read-position bias tallies
	// (bam2bcf.h:119). Sized by bca->npos in upstream; we let the
	// indel-calling sub-slices size them at use-time.
	IrefPos []int
	IaltPos []int
	// IrefMq / IaltMq are the indel mapping-quality bias tallies
	// (bam2bcf.h:119), sized by bca->nqual.
	IrefMq []int
	IaltMq []int
	// IrefScl / IaltScl are the indel soft-clip length bias tallies
	// (bam2bcf.h:121). Upstream fixes the size at 100; we follow.
	IrefScl [100]int
	IaltScl [100]int

	// Nnm / Nm count NM observations on the ref and alt reads
	// (bam2bcf.h:137-138): index 0 is REF, index 1 is ALT.
	Nnm [2]uint32
	Nm  [2]float64

	// Chr is the current chromosome name, threaded through the indel
	// helpers (bam2bcf.h:140).
	Chr string

	// DelBias is upstream's del_bias (bam2bcf.h:134); >0 prefers
	// deletions, <0 prefers insertions.
	DelBias float64
	// VsRef is upstream's vs_ref (bam2bcf.h:135); 0 scores vs the
	// next-best allele, 1 vs the reference.
	VsRef float64
}

// newBcfCallauxIndel returns a fresh indel-aux struct populated with
// the upstream mpileup defaults (mpileup.c:1381-1383 and the indel
// getopt long table). The CLI plumbing in subcmds_mpileup.go fills
// MpileupOptions and validateMpileupOptions feeds those through here.
func newBcfCallauxIndel(opts MpileupOptions) *bcfCallauxIndel {
	openQ := opts.OpenProb
	if openQ == 0 {
		openQ = DefaultMpileupOpenProb
	}
	extQ := opts.ExtProb
	if extQ == 0 {
		extQ = DefaultMpileupExtProb
	}
	tandemQ := opts.TandemQual
	if tandemQ == 0 {
		tandemQ = DefaultMpileupTandemQual
	}
	minSupport := opts.MinIReads
	if minSupport <= 0 {
		minSupport = DefaultMpileupMinIReads
	}
	minFrac := opts.GapFrac
	if minFrac == 0 {
		minFrac = DefaultMpileupGapFrac
	}
	indelBias := opts.IndelBias
	if indelBias == 0 {
		indelBias = DefaultMpileupIndelBias
	}
	indelWin := opts.IndelSize
	if indelWin == 0 {
		indelWin = DefaultMpileupIndelSize
	}
	polyMQ := 0
	if opts.PolyMQual {
		polyMQ = 1
	}
	return &bcfCallauxIndel{
		OpenQ:        openQ,
		ExtQ:         extQ,
		TandemQ:      tandemQ,
		MinSupport:   minSupport,
		MinFrac:      minFrac,
		IndelBias:    indelBias,
		IndelWinSize: indelWin,
		DelBias:      opts.DelBias,
		ReadLen:      opts.MaxReadLen,
		SeqQOffset:   opts.SeqQOffset,
		PolyMQual:    polyMQ,
		VsRef:        opts.ScoreVsRef,
	}
}

// estSeqQ is the Go port of bam2bcf_indel.c:82 est_seqQ. It returns a
// phred-scaled quality estimate for an indel of relative length l on a
// homopolymer run of length lRun. l can be negative (deletion); only
// |l| matters. lRun is the length of the homopolymer at the candidate
// site, as returned by bcfCgpLRun.
func estSeqQ(bca *bcfCallauxIndel, l, lRun int) int {
	al := l
	if al < 0 {
		al = -al
	}
	q := bca.OpenQ + bca.ExtQ*(al-1)
	qh := 1000
	if lRun >= 3 {
		// Round-to-nearest mirroring upstream's +0.499 trick.
		qh = int(float64(bca.TandemQ)*float64(al)/float64(lRun) + 0.499)
	}
	if q < qh {
		return q
	}
	return qh
}

// estIndelreg is the Go port of bam2bcf_indel.c:90 est_indelreg. It
// walks the reference forward from pos+1 looking for the extent over
// which the candidate insertion (or, when ins4 is nil, a deletion of
// length l) is consistent with the local sequence, scoring +1 per
// matching base and -10 per mismatch and returning the offset of the
// rightmost max. ref is the 0-indexed reference sequence, with N as
// the wildcard sentinel.
func estIndelreg(pos int, ref []byte, l int, ins4 []byte) int {
	if l < 0 {
		l = -l
	}
	if l == 0 || pos < 0 || pos >= len(ref) {
		return 0
	}
	max := 0
	maxI := pos
	score := 0
	// "ACGTN"[(int)ins4[j%l]] in upstream. j tracks position within the
	// repeat unit.
	for i, j := pos+1, 0; i < len(ref); i, j = i+1, j+1 {
		var want byte
		if ins4 != nil {
			if l == 0 {
				return 0
			}
			c := ins4[j%l]
			if int(c) < 5 {
				want = "ACGTN"[c]
			} else {
				want = 'N'
			}
		} else {
			ri := pos + 1 + j%l
			if ri < 0 || ri >= len(ref) {
				break
			}
			want = upperByte(ref[ri])
		}
		got := upperByte(ref[i])
		if got != want {
			score -= 10
		} else {
			score++
		}
		if score < 0 {
			break
		}
		if max < score {
			max = score
			maxI = i
		}
	}
	return maxI - pos
}

// bcfCgpLRun is the Go port of bam2bcf_indel.c:415 bcf_cgp_l_run. It
// returns the length of the homopolymer surrounding pos+1: i.e. the
// number of consecutive reference bases equal to ref[pos+1], extended
// in both directions from pos+1. When ref[pos+1] is N (any) the answer
// is 1.
func bcfCgpLRun(ref []byte, pos int) int {
	if pos+1 < 0 || pos+1 >= len(ref) {
		return 1
	}
	c := baseToNt16(upperByte(ref[pos+1]))
	if c == 15 {
		return 1
	}
	// Extend right from pos+2 while bases match.
	i := pos + 2
	for i < len(ref) {
		if baseToNt16(upperByte(ref[i])) != c {
			break
		}
		i++
	}
	lRun := i
	// Extend left from pos.
	for i = pos; i >= 0; i-- {
		if baseToNt16(upperByte(ref[i])) != c {
			break
		}
	}
	lRun -= i + 1
	return lRun
}

// indelObs is one indel observation extracted from a column's
// pileupBase slice. It mirrors the (sample-id, indel-length) tuples
// upstream stores in the per-call aux array.
type indelObs struct {
	sample int
	indel  int
}

// indelTypesResult carries the outputs of bcfCgpFindTypes.
type indelTypesResult struct {
	// Types is the sorted, deduplicated list of indel lengths
	// observed at the candidate site, with 0 (REF) always included.
	Types []int
	// MaxRdLen is the longest read query length seen across the
	// per-sample piles, used to size the alignment window.
	MaxRdLen int
	// RefType is the index in Types of the 0-length REF entry.
	RefType int
	// N is the total number of reads in the per-sample piles.
	N int
}

// bcfCgpFindTypes is the Go port of bam2bcf_indel.c:165
// bcf_cgp_find_types. piles is a per-sample slice of pileup columns at
// the same reference position. ref is the reference sequence (0-based)
// and pos is the 0-based reference position. The function fills the
// per-call indel-aux state (max_support, max_frac) on bca and returns
// an indelTypesResult, or nil if the candidate site should be skipped
// (insufficient support, too many types, or too many N's in the ref).
func bcfCgpFindTypes(piles [][]pileupBase, pos int, bca *bcfCallauxIndel,
	ref []byte) *indelTypesResult {

	// Reset per-site tallies the function is documented to fill
	// (bam2bcf_indel.c:178).
	bca.MaxSupport = 0
	bca.MaxFrac = 0

	// Tally distinct indel lengths together with the read-pile total.
	// types[0] always carries the REF (0-length) sentinel.
	N := 0
	for _, p := range piles {
		N += len(p)
	}

	// Collect observed indel lengths and per-sample (na, nt) counts.
	aux := make([]int, 0, N+1)
	aux = append(aux, 0) // REF
	maxRdLen := 0
	nAlt := 0
	nTot := 0
	indelSupportOK := false
	for _, pile := range piles {
		na := 0
		nt := 0
		for _, p := range pile {
			nt++
			if p.indel != 0 {
				na++
				aux = append(aux, p.indel)
			}
			// Read length: upstream calls bam_cigar2qlen on the
			// record. When the rec back-pointer is set we can do the
			// same; otherwise fall back to qlen.
			rdLen := p.qlen
			if p.rec != nil {
				rdLen = 0
				for _, op := range p.rec.Cigar {
					switch op.Op() {
					case sam.CigarMatch, sam.CigarInsertion, sam.CigarSoftClip,
						sam.CigarEqual, sam.CigarMismatch:
						rdLen += int(op.Length())
					}
				}
			}
			if rdLen > maxRdLen {
				maxRdLen = rdLen
			}
		}
		// Per-sample support test (bam2bcf_indel.c:203).
		if nt > 0 {
			frac := float64(na) / float64(nt)
			if !indelSupportOK && na >= bca.MinSupport && frac >= bca.MinFrac {
				indelSupportOK = true
			}
			if na > bca.MaxSupport && frac > 0 {
				bca.MaxSupport = na
				bca.MaxFrac = frac
			}
		}
		nAlt += na
		nTot += nt
	}

	// Sort + dedup. Upstream uses a 32-bit MINUS_CONST trick to keep
	// negative deletions valid for an unsigned sort; integer sort here
	// is equivalent without the bias.
	sort.Ints(aux)
	nTypes := 0
	for i := range aux {
		if i == 0 || aux[i] != aux[i-1] {
			nTypes++
		}
	}

	// Site-level support test (bam2bcf_indel.c:219). When per_sample
	// filtering is OFF (PerSampleFlt == 0), upstream OVERRIDES the
	// per-sample-OK flag with the site-level test.
	if bca.PerSampleFlt == 0 {
		ok := true
		if nTot == 0 || float64(nAlt)/float64(nTot) < bca.MinFrac || nAlt < bca.MinSupport {
			ok = false
		}
		indelSupportOK = ok
	}
	if nTypes == 1 || !indelSupportOK {
		return nil
	}
	if nTypes >= b2bIndelTypesMax {
		return nil
	}

	// Reject sites with too many N's in the ref window
	// (bam2bcf_indel.c:241).
	iEnd := pos + 2*bca.IndelWinSize
	if 2*bca.IndelWinSize > maxRdLen {
		iEnd = pos + maxRdLen
	}
	if iEnd > len(ref) {
		iEnd = len(ref)
	}
	nN := 0
	i := pos
	for ; i < iEnd; i++ {
		if upperByte(ref[i]) == 'N' {
			nN++
		}
	}
	if nN*2 > (i - pos) {
		return nil
	}

	// Fill out types[] in sorted order.
	types := make([]int, 0, nTypes)
	for i := range aux {
		if i == 0 || aux[i] != aux[i-1] {
			types = append(types, aux[i])
		}
	}
	refType := 0
	for t, v := range types {
		if v == 0 {
			refType = t
			break
		}
	}
	return &indelTypesResult{
		Types:    types,
		MaxRdLen: maxRdLen,
		RefType:  refType,
		N:        N,
	}
}

// tpos2qpos is the Go port of bam2bcf_indel.c:51. It converts a
// reference position tpos into a query (sequence) position using the
// CIGAR string anchored at recPos (0-based reference start of the
// read). For deletions tpos may not be covered by the read; isLeft
// chooses whether the returned query position should be the deletion's
// left edge (true) or right edge (false). The clamped reference
// position used in the conversion is written to *tposOut.
func tpos2qpos(cigar sam.Cigar, recPos int, tpos int, isLeft bool, tposOut *int) int {
	x := recPos
	y := 0
	lastY := 0
	*tposOut = recPos
	for _, op := range cigar {
		o := op.Op()
		l := int(op.Length())
		switch o {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			if recPos > tpos {
				return y
			}
			if x+l > tpos {
				*tposOut = tpos
				return y + (tpos - x)
			}
			x += l
			y += l
			lastY = y
		case sam.CigarInsertion, sam.CigarSoftClip:
			y += l
		case sam.CigarDeletion, sam.CigarSkipped:
			if x+l > tpos {
				if isLeft {
					*tposOut = x
				} else {
					*tposOut = x + l
				}
				return y
			}
			x += l
		}
	}
	*tposOut = x
	return lastY
}

// getPosResult carries the outputs of getPos.
type getPosResult struct {
	ScLen int // soft-clip-length bias bin
	SLen  int // clipped sequence length (read length minus soft-clips)
	EPos  int // read-position bin (rescaled into bca->npos bins)
	End   int // -1 if no soft-clip, 0 if left clip closer, 1 if right
}

// getPos is the Go port of bam2bcf_indel.c:104 get_pos. It computes the
// soft-clip-length bias bin, the clipped sequence length, and the
// read-position bin for a pileup base. The caller supplies the per-base
// state via a pileupBase together with the bias-bin count bca.npos
// (b2bNpos). The returned ScLen is capped at 99 per upstream
// (bam2bcf_indel.c:149).
//
// indel matches bam_pileup1_t.indel: positive for an insertion of that
// length immediately after this column, negative for a deletion, 0
// otherwise. cigar/recPos are taken from the originating read.
func getPos(cigar sam.Cigar, recPos int, qpos, qlen, indel int) getPosResult {
	scLen := 0
	scDist := -1
	atLeft := true
	epos := qpos
	slen := qlen
	end := -1
	for _, op := range cigar {
		o := op.Op()
		l := int(op.Length())
		switch o {
		case sam.CigarSoftClip:
			slen -= l
			if atLeft {
				scLen += l
				epos -= scLen
				scDist = epos
				end = 0
			} else {
				rd := qlen - l - qpos
				if scDist < 0 || scDist > rd {
					scDist = rd
					scLen = l
					end = 1
				}
			}
		case sam.CigarHardClip:
			// no change to atLeft
		default:
			atLeft = false
		}
		// Unused-vars guard for recPos (kept for parity with upstream's
		// signature; the cigar walk does not need the reference start).
		_ = recPos
	}

	if indel > 0 && slen-(epos+indel) < epos {
		epos += indel - 1
	}

	// Rescale epos into bca->npos bins.
	eposBin := 0
	if slen+1 > 0 {
		eposBin = int(float64(epos) / float64(slen+1) * float64(b2bNpos))
	}

	scLenBin := 0
	if scLen > 0 {
		scLenBin = 15 * scLen / (scDist + 1)
		if scLenBin > 99 {
			scLenBin = 99
		}
	}
	return getPosResult{
		ScLen: scLenBin,
		SLen:  slen,
		EPos:  eposBin,
		End:   end,
	}
}
