package samtools

import "github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"

// Mate-pair overlap detection and removal, a faithful port of htslib's
// tweak_overlap_quality / overlap_push (sam.c). samtools mpileup applies it by
// DEFAULT (disabled with -x/--ignore-overlaps-removal): where the two mates of
// a read pair overlap on the reference, one mate's base qualities are zeroed so
// that base is dropped by the -Q min-base-quality filter, and — for agreeing
// bases — the kept mate's quality is raised to the (capped) SUM of the two, so
// the overlapping fragment is counted once with the combined confidence. Which
// mate is kept is chosen by a hash of the read name, exactly as htslib does, so
// the result is byte-identical to upstream.

// x31Hash is khash.h's __ac_X31_hash_string over the read name.
func x31Hash(s string) uint32 {
	if len(s) == 0 {
		return 0
	}
	h := uint32(s[0])
	for i := 1; i < len(s); i++ {
		h = (h << 5) - h + uint32(s[i])
	}
	return h
}

// wangHash is khash.h's __ac_Wang_hash, Thomas Wang's 32-bit integer hash.
func wangHash(key uint32) uint32 {
	key += ^(key << 15)
	key ^= key >> 10
	key += key << 3
	key ^= key >> 6
	key += ^(key << 11)
	key ^= key >> 16
	return key
}

// cigIter walks a record's CIGAR mapping reference offset to sequence index,
// a port of htslib's cigar_iref2iseq_set / cigar_iref2iseq_next. idx is the
// current CIGAR-op index (htslib's advancing cigar pointer), icig the position
// within that op, iseq the sequence index and iref the reference offset from
// the read's start.
type cigIter struct {
	cig  sam.Cigar
	idx  int
	icig int64
	iseq int64
	iref int64
}

// cigMatch is the success return (BAM_CMATCH == 0); negatives mean done (-1) or
// a malformed CIGAR (-2), matching the htslib helpers.
const cigMatch = 0

// set seeks to the first CIGAR match at or after reference offset target,
// resetting the per-record cursors (cigar_iref2iseq_set).
func (it *cigIter) set(target int64) int {
	pos := target
	if pos < 0 {
		return -1
	}
	it.icig, it.iseq, it.iref = 0, 0, 0
	for it.idx < len(it.cig) {
		op := it.cig[it.idx].Op()
		ncig := int64(it.cig[it.idx].Length())
		switch op {
		case sam.CigarSoftClip:
			it.idx++
			it.iseq += ncig
			it.icig = 0
		case sam.CigarHardClip, sam.CigarPadding:
			it.idx++
			it.icig = 0
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			pos -= ncig
			if pos < 0 {
				it.icig = ncig + pos
				it.iseq += it.icig
				it.iref += it.icig
				return cigMatch
			}
			it.idx++
			it.iseq += ncig
			it.icig = 0
			it.iref += ncig
		case sam.CigarInsertion:
			it.idx++
			it.iseq += ncig
			it.icig = 0
		case sam.CigarDeletion, sam.CigarSkipped:
			pos -= ncig
			if pos < 0 {
				pos = 0
			}
			it.idx++
			it.icig = 0
			it.iref += ncig
		default:
			return -2
		}
	}
	it.iseq = -1
	return -1
}

// next advances to the following CIGAR match base (cigar_iref2iseq_next).
func (it *cigIter) next() int {
	for it.idx < len(it.cig) {
		op := it.cig[it.idx].Op()
		ncig := int64(it.cig[it.idx].Length())
		switch op {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			if it.icig >= ncig-1 {
				it.icig = -1
				it.idx++
				continue
			}
			it.iseq++
			it.icig++
			it.iref++
			return cigMatch
		case sam.CigarDeletion, sam.CigarSkipped:
			it.idx++
			it.iref += ncig
			it.icig = -1
		case sam.CigarInsertion, sam.CigarSoftClip:
			it.idx++
			it.iseq += ncig
			it.icig = -1
		case sam.CigarHardClip, sam.CigarPadding:
			it.idx++
			it.icig = -1
		default:
			return -2
		}
	}
	it.iseq, it.iref = -1, -1
	return -1
}

// prevOpIsDel reports whether the CIGAR op just consumed (cigar[-1]) was a
// deletion, the htslib guard for the deletion catch-up branches.
func (it *cigIter) prevOpIsDel() bool {
	return it.idx > 0 && it.cig[it.idx-1].Op() == sam.CigarDeletion
}

// negRet maps a helper return to htslib's -(ret<-1): -1 for a malformed CIGAR
// (-2), else 0.
func negRet(r int) int {
	if r < -1 {
		return -1
	}
	return 0
}

// tweakOverlapQuality adjusts the base qualities of the two overlapping mates a
// (left) and b (right) in place, a port of htslib's tweak_overlap_quality.
func tweakOverlapQuality(a, b *sam.Record) int {
	aIt := &cigIter{cig: a.Cigar}
	bIt := &cigIter{cig: b.Cigar}
	aPos := a.Pos - 1
	bPos := b.Pos - 1
	iref := bPos

	aRet := aIt.set(iref - aPos)
	if aRet < 0 {
		return negRet(aRet)
	}
	bRet := bIt.set(iref - bPos)
	if bRet < 0 {
		return negRet(bRet)
	}

	// Pick which mate keeps quality (mul==1) and which is zeroed (mul==0).
	var amul, bmul int
	if wangHash(x31Hash(a.QName))&1 != 0 {
		amul, bmul = 1, 0
	} else {
		amul, bmul = 0, 1
	}

	aQual, bQual := a.Qual, b.Qual
	aSeq, bSeq := a.Seq, b.Seq

	for {
		for aRet >= 0 && aIt.iref >= 0 && aIt.iref < iref-aPos {
			aRet = aIt.next()
		}
		if aRet < 0 {
			return negRet(aRet)
		}
		for bRet >= 0 && bIt.iref >= 0 && bIt.iref < iref-bPos {
			bRet = bIt.next()
		}
		if bRet < 0 {
			return negRet(bRet)
		}

		if iref < aIt.iref+aPos {
			iref = aIt.iref + aPos
		}
		if iref < bIt.iref+bPos {
			iref = bIt.iref + bPos
		}
		iref++

		// A deletion in one mate moves it further along the reference; catch the
		// other up, amending its qualities with the mismatch (0.8x / zero) rule.
		if aIt.iref+aPos != bIt.iref+bPos {
			switch {
			case aIt.iref+aPos < bIt.iref+bPos && bIt.prevOpIsDel():
				for {
					mulQual(aQual, aIt.iseq, amul)
					aRet = aIt.next()
					if aRet < 0 {
						return negRet(aRet)
					}
					if aIt.iref+aPos >= bIt.iref+bPos {
						break
					}
				}
			case aIt.prevOpIsDel():
				for {
					mulQual(bQual, bIt.iseq, bmul)
					bRet = bIt.next()
					if bRet < 0 {
						return negRet(bRet)
					}
					if bIt.iref+bPos >= aIt.iref+aPos {
						break
					}
				}
			default:
				// ref-skip etc — unsupported here, as in htslib.
				continue
			}
		}

		if aIt.iseq > int64(len(aSeq)) || bIt.iseq > int64(len(bSeq)) {
			return -1
		}

		if aSeq[aIt.iseq] == bSeq[bIt.iseq] {
			// Confident: keep the (capped) sum of the two qualities on the keeper.
			qual := int(aQual[aIt.iseq]) + int(bQual[bIt.iseq])
			if qual > 200 {
				qual = 200
			}
			aQual[aIt.iseq] = byte(amul * qual)
			bQual[bIt.iseq] = byte(bmul * qual)
		} else {
			// Disagreement: keep the higher-quality base at 0.8x, zero the other.
			switch {
			case aQual[aIt.iseq] > bQual[bIt.iseq]:
				aQual[aIt.iseq] = byte(0.8 * float64(aQual[aIt.iseq]))
				bQual[bIt.iseq] = 0
			case aQual[aIt.iseq] < bQual[bIt.iseq]:
				bQual[bIt.iseq] = byte(0.8 * float64(bQual[bIt.iseq]))
				aQual[aIt.iseq] = 0
			default:
				aQual[aIt.iseq] = byte(float64(amul) * 0.8 * float64(aQual[aIt.iseq]))
				bQual[bIt.iseq] = byte(float64(bmul) * 0.8 * float64(bQual[bIt.iseq]))
			}
		}
	}
}

// mulQual applies the deletion catch-up rule: keep the quality at 0.8x when
// mul==1, else zero it (htslib's `q = mul ? q*0.8 : 0`).
func mulQual(qual []byte, i int64, mul int) {
	if mul != 0 {
		qual[i] = byte(0.8 * float64(qual[i]))
	} else {
		qual[i] = 0
	}
}

// applyOverlapRemoval runs overlapPush over a chrom's coordinate-sorted records
// (one input), pairing mates and de-weighting their overlap in place. It is the
// buffered walk's equivalent of the streaming path's per-read overlapPush.
func applyOverlapRemoval(recs []*sam.Record) {
	overlaps := make(map[string]*sam.Record)
	for _, rec := range recs {
		overlapPush(rec, overlaps)
	}
}

// overlapPush pairs a record with its already-seen mate and applies
// tweakOverlapQuality, a port of htslib's overlap_push. Records are presented in
// coordinate order; the earlier-arriving mate of each pair is held in overlaps
// until the later one is seen. Only proper pairs with a mapped, same-contig mate
// that can actually overlap are considered.
func overlapPush(rec *sam.Record, overlaps map[string]*sam.Record) {
	if rec.Flag&sam.FlagMateUnmapped != 0 || rec.Flag&sam.FlagProperPair == 0 {
		return
	}
	// Mate on a different contig cannot overlap.
	if rec.RNext != "=" && rec.RNext != rec.RName {
		return
	}
	mpos0 := rec.PNext - 1 // 0-based mate position; PNext==0 -> -1 (absent)
	lqseq := int64(len(rec.Seq))
	isize := rec.TLen
	if isize < 0 {
		isize = -isize
	}
	// No overlap possible for a wild CIGAR whose insert size already clears the
	// read and whose mate starts past this read's end.
	if isize >= 2*lqseq && mpos0 >= rec.EndPosition() {
		return
	}

	if a, ok := overlaps[rec.QName]; ok {
		tweakOverlapQuality(a, rec)
		delete(overlaps, rec.QName)
		return
	}
	// Only buffer the first-seen mate when its pair is still to arrive.
	if mpos0 >= rec.Pos-1 || (rec.Flag&sam.FlagPaired != 0 && mpos0 == -1) {
		overlaps[rec.QName] = rec
	}
}
