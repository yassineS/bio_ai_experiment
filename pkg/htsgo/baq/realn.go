package baq

import (
	"math"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// Realn flag bits, mirroring htslib's enum htsRealnFlags (sam.h). They are
// combined into the flag argument of SamProbRealn.
const (
	// FlagApply caps the base qualities by the computed BAQ and records the
	// adjustment in a ZQ:Z: tag (htslib BAQ_APPLY, the calmd -A flag).
	FlagApply = 1
	// FlagExtend selects extended-BAQ mode, which can raise qualities above
	// the input values (htslib BAQ_EXTEND, the calmd -E flag).
	FlagExtend = 2
	// FlagRedo forces recomputation, discarding any pre-existing BQ tag
	// (htslib BAQ_REDO).
	FlagRedo = 4
)

// platformShift isolates the platform subfield of the realn flag (bits 3+).
const platformShift = 3

// baqIllumina is the platform subfield value for Illumina reads; values above
// it (or reads longer than 1000 bp) switch to long-read tuned parameters.
const baqIllumina = 1 << platformShift

// asciiToResidue maps an ASCII nucleotide byte to the 0/1/2/3/4 residue
// encoding probaln_glocal expects (A/C/G/T, 4 for anything else). It folds
// case and reproduces htslib's seq_nt16_int[seq_nt16_table[...]] composition.
var asciiToResidue = func() [256]byte {
	var t [256]byte
	for i := range t {
		t[i] = 4
	}
	t['A'], t['a'] = 0, 0
	t['C'], t['c'] = 1, 1
	t['G'], t['g'] = 2, 2
	t['T'], t['t'] = 3, 3
	return t
}()

// SamProbRealn is a faithful port of htslib's sam_prob_realn. It computes BAQ
// for rec against ref (the full 0-based reference contig sequence), and writes
// the result back into rec: a BQ:Z: aux tag, or — when FlagApply is set — the
// base qualities lowered by BAQ with a ZQ:Z: tag recording the delta.
//
// flag is a bitwise-OR of FlagApply / FlagExtend / FlagRedo plus an optional
// platform subfield. The return value mirrors upstream: 0 on success, -1 when
// the read needs no work (unmapped, empty, no qualities, no M/=/X op, or a
// reference skip), -3 when BAQ is already in the requested state, and -4 on
// failure.
func SamProbRealn(rec *sam.Record, ref []byte, flag int) int {
	applyBAQ := flag&FlagApply != 0
	extendBAQ := flag&FlagExtend != 0
	redoBAQ := flag&FlagRedo != 0
	system := flag & (0xff << platformShift)

	// lqseq is the query length (htslib core.l_qseq). Callers must pass a
	// record with a real SEQ — a record whose SEQ is "*" would mis-size
	// here; samtools calmd filters SEQ=="*" before invoking BAQ.
	lqseq := len(rec.Seq)
	refLen := len(ref)

	// d(I) e(M) band — Illumina defaults.
	conf := Par{D: 0.001, E: 0.1, BW: 10}
	if lqseq > 1000 || system > baqIllumina {
		// PacBio CCS-tuned parameters for long reads.
		conf.D = 1e-7
		conf.E = 1e-1
	}

	if rec.IsUnmapped() || lqseq == 0 || len(rec.Qual) == 0 {
		return -1
	}
	// A QUAL of "*" is stored as all 0xff; upstream tests qual[0].
	if rec.Qual[0] == 0xff {
		return -1
	}

	qual := rec.Qual

	// Locate any existing BQ / ZQ aux tags.
	bqIdx := auxIndexOf(rec, "BQ")
	zqIdx := auxIndexOf(rec, "ZQ")
	fixBQ := false
	if bqIdx >= 0 {
		if !redoBAQ && !validBAQTag(rec.Aux[bqIdx], lqseq) {
			fixBQ = true
		}
	}
	if zqIdx >= 0 {
		if !validBAQTag(rec.Aux[zqIdx], lqseq) {
			return -4
		}
	}
	if bqIdx >= 0 && redoBAQ {
		delAux(rec, "BQ")
		bqIdx = -1
		zqIdx = auxIndexOf(rec, "ZQ")
	}
	if bqIdx >= 0 && zqIdx >= 0 {
		// Remove the ZQ tag, keeping BQ.
		delAux(rec, "ZQ")
		zqIdx = -1
	}
	if zqIdx < 0 && fixBQ {
		// Invalid BQ tag with no ZQ: drop it and realign.
		delAux(rec, "BQ")
		bqIdx = -1
		fixBQ = false
	}

	if bqIdx >= 0 || zqIdx >= 0 {
		if (applyBAQ && zqIdx >= 0) || (!applyBAQ && bqIdx >= 0) {
			return -3 // already in the requested state
		}
		if bqIdx >= 0 && applyBAQ {
			// Convert BQ to ZQ: lower the qualities, retag.
			bq, _ := rec.Aux[bqIdx].String()
			for i := 0; i < lqseq; i++ {
				if int(qual[i])+64 < int(bq[i]) {
					qual[i] = 0
				} else {
					qual[i] = qual[i] - byte(int(bq[i])-64)
				}
			}
			rec.Aux[bqIdx].Tag = "ZQ"
			rec.InvalidateAuxIndex()
		} else if zqIdx >= 0 && !applyBAQ {
			// Convert ZQ to BQ: restore the qualities, retag.
			zq, _ := rec.Aux[zqIdx].String()
			for i := 0; i < lqseq; i++ {
				qual[i] += byte(int(zq[i]) - 64)
			}
			rec.Aux[zqIdx].Tag = "BQ"
			rec.InvalidateAuxIndex()
		}
		return 0
	}

	// Find the start and end of the alignment.
	x := int(rec.Pos) - 1
	if x < 0 {
		x = 0
	}
	y := 0
	yb, ye, xb, xe := -1, -1, -1, -1
	for _, op := range rec.Cigar {
		opCode := op.Op()
		l := int(op.Length())
		switch opCode {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			if yb < 0 {
				yb = y
			}
			if xb < 0 {
				xb = x
			}
			ye = y + l
			xe = x + l
			x += l
			y += l
		case sam.CigarSoftClip, sam.CigarInsertion:
			y += l
		case sam.CigarDeletion:
			x += l
		case sam.CigarSkipped:
			return -1 // do nothing if there is a reference skip
		}
	}
	if xb == -1 { // No matches in CIGAR.
		return -1
	}

	// Bandwidth and window start / end.
	bw := 7
	if d := absInt((xe - xb) - (ye - yb)); d > bw {
		bw = d + 3
	}
	conf.BW = bw

	xb -= yb + bw/2
	if xb < 0 {
		xb = 0
	}
	xe += lqseq - ye + bw/2
	if xe-xb-lqseq > bw {
		xb += (xe - xb - lqseq - bw) / 2
		xe -= (xe - xb - lqseq - bw) / 2
	}

	// Translate read and reference windows into 0/1/2/3/4 residues.
	tseq := make([]byte, lqseq)
	seq := rec.Seq
	for i := 0; i < lqseq; i++ {
		tseq[i] = asciiToResidue[seq[i]]
	}
	tref := make([]byte, 0, xe-xb)
	for i := xb; i < xe; i++ {
		if i >= refLen || ref[i] == 0 {
			xe = i
			break
		}
		tref = append(tref, asciiToResidue[ref[i]])
	}

	state := make([]int, lqseq)
	q := make([]byte, lqseq)
	if _, err := ProbalnGlocal(tref, tseq, qual, conf, state, q); err != nil {
		return -4
	}

	// bq starts as a copy of the base qualities (htslib's memcpy(bq, qual,
	// l_qseq)); soft-clip / insertion positions keep that value because the
	// CIGAR walks below only overwrite alignment-match positions.
	bq := make([]byte, lqseq)
	copy(bq, qual)
	if !extendBAQ {
		// bq[] is capped by base quality qual[].
		x, y = int(rec.Pos)-1, 0
		if x < 0 {
			x = 0
		}
		for _, op := range rec.Cigar {
			opCode := op.Op()
			l := int(op.Length())
			if l == 0 {
				continue
			}
			switch opCode {
			case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
				if l > lqseq-y {
					l = lqseq - y
				}
				for i := y; i < y+l; i++ {
					if state[i]&3 != 0 || state[i]>>2 != x-xb+(i-y) {
						bq[i] = 0
					} else if q[i] < qual[i] {
						bq[i] = q[i]
					} else {
						bq[i] = qual[i]
					}
				}
				x += l
				y += l
			case sam.CigarSoftClip, sam.CigarInsertion:
				if l > lqseq-y {
					l = lqseq - y
				}
				y += l
			case sam.CigarDeletion:
				x += l
			}
		}
		for i := 0; i < lqseq; i++ {
			bq[i] = qual[i] - bq[i] + 64 // finalize BQ
		}
	} else {
		// bq[] is BAQ that can be larger than qual[].
		left := make([]byte, lqseq)
		rght := make([]byte, lqseq)
		x, y = int(rec.Pos)-1, 0
		if x < 0 {
			x = 0
		}
		runLen := 0
		cigar := rec.Cigar
		for k := 0; k < len(cigar); k++ {
			opCode := cigar[k].Op()
			l := int(cigar[k].Length())

			// Concatenate consecutive alignment-match ops so that e.g.
			// 50M50M behaves identically to 100M.
			if opCode == sam.CigarMatch || opCode == sam.CigarEqual || opCode == sam.CigarMismatch {
				if k+1 < len(cigar) {
					nextOp := cigar[k+1].Op()
					if nextOp == sam.CigarMatch || nextOp == sam.CigarEqual || nextOp == sam.CigarMismatch {
						runLen += l
						continue
					}
				}
				l += runLen
				runLen = 0
			}

			if l == 0 {
				continue
			}
			switch opCode {
			case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
				if l > lqseq-y {
					l = lqseq - y
				}
				for i := y; i < y+l; i++ {
					if state[i]&3 != 0 || state[i]>>2 != x-xb+(i-y) {
						bq[i] = 0
					} else {
						bq[i] = q[i]
					}
				}
				left[y] = bq[y]
				for i := y + 1; i < y+l; i++ {
					if bq[i] > left[i-1] {
						left[i] = bq[i]
					} else {
						left[i] = left[i-1]
					}
				}
				rght[y+l-1] = bq[y+l-1]
				for i := y + l - 2; i >= y; i-- {
					if bq[i] > rght[i+1] {
						rght[i] = bq[i]
					} else {
						rght[i] = rght[i+1]
					}
				}
				for i := y; i < y+l; i++ {
					if left[i] < rght[i] {
						bq[i] = left[i]
					} else {
						bq[i] = rght[i]
					}
				}
				x += l
				y += l
			case sam.CigarSoftClip, sam.CigarInsertion:
				if l > lqseq-y {
					l = lqseq - y
				}
				y += l
			case sam.CigarDeletion:
				x += l
			}
		}
		for i := 0; i < lqseq; i++ {
			if qual[i] <= bq[i] {
				bq[i] = 64
			} else {
				bq[i] = 64 + qual[i] - bq[i] // finalize BQ
			}
		}
	}

	if applyBAQ {
		for i := 0; i < lqseq; i++ {
			qual[i] -= bq[i] - 64 // modify qual
		}
		rec.Aux = append(rec.Aux, sam.Aux{Tag: "ZQ", Type: 'Z', Value: string(bq)})
	} else {
		rec.Aux = append(rec.Aux, sam.Aux{Tag: "BQ", Type: 'Z', Value: string(bq)})
	}
	rec.InvalidateAuxIndex()
	return 0
}

// SamCapMapq is a faithful port of htslib's sam_cap_mapq. It returns a capped
// mapping quality derived from the count and quality of mismatches between rec
// and ref. A return of -1 means "no cap" (the read is clean enough that thres
// is not exceeded). thres < 0 selects the upstream default of 40.
func SamCapMapq(rec *sam.Record, ref []byte, thres int) int {
	if thres < 0 {
		thres = 40
	}
	seq := rec.Seq
	qual := rec.Qual
	refLen := len(ref)

	mm, q, length, clipL, clipQ := 0, 0, 0, 0, 0
	x := int(rec.Pos) - 1
	if x < 0 {
		x = 0
	}
	y := 0
	for _, op := range rec.Cigar {
		opCode := op.Op()
		l := int(op.Length())
		switch opCode {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			j := 0
			for ; j < l; j++ {
				z := y + j
				if x+j >= refLen || ref[x+j] == 0 || z >= len(seq) || z >= len(qual) {
					break
				}
				c1 := asciiToResidue[seq[z]] // 0-3, or 4 for ambiguous
				c2 := asciiToResidue[ref[x+j]]
				// Ambiguous residues (4) correspond to upstream's nt16
				// code 15; skip them.
				if c2 != 4 && c1 != 4 && qual[z] >= 13 {
					length++
					if seq[z] != '=' && c1 != c2 {
						mm++
						if qual[z] > 33 {
							q += 33
						} else {
							q += int(qual[z])
						}
					}
				}
			}
			if j < l {
				return capMapqResult(mm, q, length, clipQ, thres)
			}
			x += l
			y += l
			length += l
		case sam.CigarDeletion:
			j := 0
			for ; j < l; j++ {
				if x+j >= refLen || ref[x+j] == 0 {
					break
				}
			}
			if j < l {
				return capMapqResult(mm, q, length, clipQ, thres)
			}
			x += l
		case sam.CigarSoftClip:
			for j := 0; j < l; j++ {
				if y+j < len(qual) {
					clipQ += int(qual[y+j])
				}
			}
			clipL += l
			y += l
		case sam.CigarHardClip:
			clipQ += 13 * l
			clipL += l
		case sam.CigarInsertion:
			y += l
		case sam.CigarSkipped:
			x += l
		}
	}
	return capMapqResult(mm, q, length, clipQ, thres)
}

// capMapqResult evaluates the final sam_cap_mapq formula given the accumulated
// mismatch / clip statistics.
func capMapqResult(mm, q, length, clipQ, thres int) int {
	t := 1.0
	for i := 0; i < mm; i++ {
		t *= float64(length) / float64(i+1)
	}
	tf := float64(q) - 4.343*math.Log(t) + float64(clipQ)/5.0
	if tf > float64(thres) {
		return -1
	}
	if tf < 0 {
		tf = 0
	}
	tf = math.Sqrt((float64(thres)-tf)/float64(thres)) * float64(thres)
	return int(tf + 0.499)
}

// auxIndexOf returns the slice index of the aux field with the given tag, or
// -1 if absent.
func auxIndexOf(rec *sam.Record, tag string) int {
	for i := range rec.Aux {
		if rec.Aux[i].Tag == tag {
			return i
		}
	}
	return -1
}

// delAux removes the first aux field with the given tag, if present, and
// invalidates the record's lazy aux index.
func delAux(rec *sam.Record, tag string) {
	for i := range rec.Aux {
		if rec.Aux[i].Tag == tag {
			rec.Aux = append(rec.Aux[:i], rec.Aux[i+1:]...)
			rec.InvalidateAuxIndex()
			return
		}
	}
}

// validBAQTag reports whether a is a string aux of exactly lqseq characters,
// the sanity check htslib's realn_check_tag applies to BQ / ZQ tags.
func validBAQTag(a sam.Aux, lqseq int) bool {
	if a.Type != 'Z' {
		return false
	}
	s, ok := a.String()
	return ok && len(s) == lqseq
}

// absInt returns the absolute value of an int.
func absInt(v int) int {
	if v < 0 {
		return -v
	}
	return v
}
