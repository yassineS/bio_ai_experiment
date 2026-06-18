// Full splice_csq with set_refalt for the bcftools csq haplotype engine.
//
// The haplotype engine needs the complete upstream splice behaviour:
// splice_build_hap, the synonymous-at-splice refinement, and the
// ref_beg / ref_end clipping that hap_init consumes. hapSplice ports
// csq.c's splice_t together with splice_csq / splice_csq_{mnp,complex,
// ins,del} run with set_refalt set.

package bcftools

import "github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"

// hapSplice is the Go analogue of csq.c's splice_t, including the
// kref/kalt trimmed-allele strings built by splice_build_hap.
type hapSplice struct {
	ht   *hapTranscript
	pos  int // VCF position, 0-based
	rlen int
	alen int
	ref  string
	alt  string

	// eng / recV / ial wire the splice arms back to csq_stage_splice so
	// splice/region/acceptor/donor consequences are staged into the
	// variant buffer exactly where upstream stages them. They are nil
	// for the set_refalt-only call inside hap_init, where staging is
	// driven separately.
	eng  *hapEngine
	recV *vcf.Variant
	ial  int

	checkAcceptor bool
	checkDonor    bool
	checkStart    bool
	checkStop     bool
	checkRegBeg   bool
	checkRegEnd   bool

	csq    uint32
	tbeg   int
	tend   int
	refBeg int // 0-based
	refEnd int // 0-based
	kref   string
	kalt   string
}

// refByte returns the padded-reference base at 0-based genomic
// coordinate pos (ht.ref[nRefPad + pos - beg0]).
func (s *hapSplice) refByte(pos int) byte {
	i := nRefPad + pos - s.ht.beg0
	if i < 0 || i >= len(s.ht.ref) {
		return 'N'
	}
	return s.ht.ref[i]
}

// refSlice returns length bytes of the padded reference starting at
// 0-based genomic coordinate pos.
func (s *hapSplice) refSlice(pos, length int) []byte {
	if length <= 0 {
		return nil
	}
	out := make([]byte, length)
	for i := 0; i < length; i++ {
		out[i] = s.refByte(pos + i)
	}
	return out
}

// spliceBuildHap ports splice_build_hap: it assembles the trimmed
// kref / kalt alleles aligned to the spliced reference around the
// splice region.
func (s *hapSplice) spliceBuildHap(beg, length int) {
	var rlen, alen, rbeg, abeg int
	if length < 0 {
		rlen, alen = -length, -length
		rbeg = beg - rlen + 1
		dlen := s.alen - s.rlen
		if dlen < 0 && beg < s.refEnd {
			dlen += s.refEnd - beg
		}
		abeg = rbeg + dlen
	} else {
		if beg < s.ht.beg0 {
			rbeg, abeg = s.ht.beg0, s.ht.beg0
			rlen, alen = 0, 0
		} else {
			rbeg, abeg = beg, beg
			rlen, alen = length, length
		}
	}

	var kref, kalt []byte

	// kref: part before vcf.ref, in vcf.ref, after vcf.ref.
	var roff int
	if rbeg < s.pos {
		kref = append(kref, s.refSlice(rbeg, s.pos-rbeg)...)
		roff = 0
	} else {
		roff = rbeg - s.pos
	}
	if roff < s.rlen && len(kref) < rlen {
		take := s.rlen - roff
		if take > rlen-len(kref) {
			take = rlen - len(kref)
		}
		kref = append(kref, s.ref[roff:roff+take]...)
	}
	end := s.pos + s.rlen
	if len(kref) < rlen {
		if end+rlen-len(kref)-1 > s.ht.end0 {
			rlen -= end + rlen - len(kref) - 1 - s.ht.end0
		}
		if len(kref) < rlen {
			kref = append(kref, s.refSlice(end, rlen-len(kref))...)
		}
	}

	// kalt.
	var aoff int
	if abeg < s.pos {
		kalt = append(kalt, s.refSlice(abeg, s.pos-abeg)...)
		aoff = 0
	} else {
		aoff = abeg - s.pos
	}
	if aoff < s.alen && len(kalt) < alen {
		take := s.alen - aoff
		if take > alen-len(kalt) {
			take = alen - len(kalt)
		}
		kalt = append(kalt, s.alt[aoff:aoff+take]...)
		aoff -= take
	}
	if aoff < 0 {
		aoff = 0
	} else {
		aoff--
	}
	end = s.pos + s.rlen
	if len(kalt) < alen {
		if end+alen+aoff-len(kalt)-1 > s.ht.end0 {
			alen -= end + alen + aoff - len(kalt) - 1 - s.ht.end0
		}
		if alen > 0 && alen > len(kalt) {
			kalt = append(kalt, s.refSlice(end+aoff, alen-len(kalt))...)
		}
	}

	s.kref = string(kref)
	s.kalt = string(kalt)
}

// hapSpliceCSQ ports splice_csq run with set_refalt: it trims the
// alleles then dispatches to the mnp / complex / ins / del arm.
func (s *hapSplice) run(exBeg, exEnd int) int {
	s.alen = len(s.alt)
	s.rlen = len(s.ref)
	rlen1, alen1 := s.rlen-1, s.alen-1
	s.tbeg, s.tend = 0, 0

	i := 0
	for i <= rlen1 && i <= alen1 {
		if s.ref[rlen1-i] != s.alt[alen1-i] {
			break
		}
		i++
	}
	s.tend = i
	rlen1 -= i
	alen1 -= i
	i = 0
	for i <= rlen1 && i <= alen1 {
		if s.ref[i] != s.alt[i] {
			break
		}
		i++
	}
	s.tbeg = i

	rtrim := s.rlen - s.tbeg - s.tend
	atrim := s.alen - s.tbeg - s.tend

	switch {
	case s.rlen == s.alen:
		return s.mnp(exBeg, exEnd)
	case rtrim > 1 && atrim > 1:
		s.csq |= s.complexBit()
		return s.mnp(exBeg, exEnd)
	case s.rlen < s.alen:
		return s.ins(exBeg, exEnd)
	default:
		return s.del(exBeg, exEnd)
	}
}

func (s *hapSplice) complexBit() uint32 {
	if s.rlen > s.alen {
		return csqTruncation
	}
	return csqElongation
}

// mnp ports splice_csq_mnp run with set_refalt.
func (s *hapSplice) mnp(exBeg, exEnd int) int {
	if s.tbeg+s.tend == s.rlen {
		return spliceVarRef
	}
	s.refBeg = s.pos + s.tbeg
	s.refEnd = s.pos + s.rlen - s.tend - 1

	ret := spliceInside
	if s.refBeg < exBeg {
		if s.checkRegBeg {
			if s.refEnd >= exBeg-nSpliceRegionIntron && s.refBeg < exBeg-nSpliceDonor {
				s.csq |= csqSpliceRegion
			}
			if s.refEnd >= exBeg-nSpliceDonor {
				if s.checkDonor && s.rev() {
					s.csq |= csqSpliceDonor
				}
				if s.checkAcceptor && !s.rev() {
					s.csq |= csqSpliceAcceptor
				}
			}
		}
		if s.refEnd >= exBeg {
			s.tbeg = s.refBeg - s.pos
			s.refBeg = exBeg
			ret = spliceOverlap
		}
	}
	if exEnd < s.refEnd {
		if s.checkRegEnd {
			if s.refBeg <= exEnd+nSpliceRegionIntron && s.refEnd > exEnd+nSpliceDonor {
				s.csq |= csqSpliceRegion
			}
			if s.refBeg <= exEnd+nSpliceDonor {
				if s.checkDonor && !s.rev() {
					s.csq |= csqSpliceDonor
				}
				if s.checkAcceptor && s.rev() {
					s.csq |= csqSpliceAcceptor
				}
			}
		}
		if s.refBeg <= exEnd {
			s.tend = s.rlen - (s.refEnd - s.pos + 1)
			s.refEnd = exEnd
			ret = spliceOverlap
		}
	}
	if s.refEnd < exBeg || s.refBeg > exEnd {
		s.stageSplice()
		return spliceOutside
	}
	if s.refBeg < exBeg+nSpliceRegionExon {
		if s.checkRegBeg {
			s.csq |= csqSpliceRegion
		}
		if !s.rev() {
			if s.checkStart {
				s.csq |= csqStartLost
			}
		} else if s.checkStop {
			s.csq |= csqStopLost
		}
	}
	if s.refEnd > exEnd-nSpliceRegionExon {
		if s.checkRegEnd {
			s.csq |= csqSpliceRegion
		}
		if s.rev() {
			if s.checkStart {
				s.csq |= csqStartLost
			}
		} else if s.checkStop {
			s.csq |= csqStopLost
		}
	}
	// set_refalt: build kref / kalt from the trimmed alleles.
	s.rlen -= s.tbeg + s.tend
	s.alen -= s.tbeg + s.tend
	s.kref = s.ref[s.tbeg : s.tbeg+s.rlen]
	s.kalt = s.alt[s.tbeg : s.tbeg+s.alen]
	s.stageSplice()
	return ret
}

// ins ports splice_csq_ins run with set_refalt.
func (s *hapSplice) ins(exBeg, exEnd int) int {
	if s.tbeg != 0 || s.ref[0] != s.alt[0] {
		s.refBeg = s.pos + s.tbeg - 1
		s.refEnd = s.pos + s.rlen - s.tend
	} else {
		if s.tend != 0 {
			s.tend--
		}
		s.refBeg = s.pos
		s.refEnd = s.pos + s.rlen - s.tend
	}

	if s.refBeg >= exEnd {
		if !s.checkRegEnd {
			return spliceOutside
		}
		s.spliceBuildHap(exEnd+1, nSpliceRegionIntron)
		ref, alt := s.kref, s.kalt
		if s.refBeg < exEnd+nSpliceRegionIntron && s.refEnd > exEnd+nSpliceDonor {
			s.csq |= csqSpliceRegion
			if eqPrefix(ref, alt, nSpliceRegionIntron) {
				s.csq |= csqSynonymous
			}
		}
		if s.refBeg < exEnd+nSpliceDonor {
			if s.checkDonor && !s.rev() {
				s.csq |= csqSpliceDonor
			}
			if s.checkAcceptor && s.rev() {
				s.csq |= csqSpliceAcceptor
			}
			if eqPrefix(ref, alt, nSpliceDonor) {
				s.csq |= csqSynonymous
			}
		}
		s.stageSplice()
		return spliceOutside
	}
	if s.refEnd < exBeg || (s.refEnd == exBeg && !s.checkRegBeg) {
		if !s.checkRegBeg {
			return spliceOutside
		}
		s.spliceBuildHap(exBeg-nSpliceRegionIntron, nSpliceRegionIntron)
		ref, alt := s.kref, s.kalt
		if s.refEnd > exBeg-nSpliceRegionIntron && s.refBeg < exBeg-nSpliceDonor {
			s.csq |= csqSpliceRegion
			if eqPrefix(ref, alt, nSpliceRegionIntron) {
				s.csq |= csqSynonymous
			}
		}
		if s.refEnd > exBeg-nSpliceDonor {
			if s.checkDonor && s.rev() {
				s.csq |= csqSpliceDonor
			}
			if s.checkAcceptor && !s.rev() {
				s.csq |= csqSpliceAcceptor
			}
			noff := nSpliceRegionIntron - nSpliceDonor
			if eqRange(ref, alt, noff, nSpliceDonor) {
				s.csq |= csqSynonymous
			}
		}
		s.stageSplice()
		return spliceOutside
	}
	if s.refBeg <= exBeg+2 {
		if s.checkRegBeg {
			s.csq |= csqSpliceRegion
		}
		if !s.rev() {
			if s.checkStart {
				s.csq |= csqStartLost
			}
		} else if s.checkStop {
			s.csq |= csqStopLost
		}
	}
	if s.refEnd > exEnd-2 {
		if s.checkRegEnd {
			s.csq |= csqSpliceRegion
		}
		if s.rev() {
			if s.checkStart {
				s.csq |= csqStartLost
			}
		} else if s.checkStop {
			s.csq |= csqStopLost
		}
	}
	// set_refalt block.
	if s.refBeg < s.pos {
		dlen := s.pos - s.refBeg
		s.tbeg += dlen
		if s.tbeg+s.tend == s.rlen {
			s.tend -= dlen
		}
		s.refBeg = s.pos
	}
	if s.refEnd == exBeg {
		s.tend--
	}
	s.spliceBuildHap(s.refBeg, s.alen-s.tend-s.tbeg+1)
	s.rlen -= s.tbeg + s.tend - 1
	if len(s.kref) > s.rlen {
		s.kref = s.kref[:s.rlen]
	}
	s.stageSplice()
	return spliceInside
}

// del ports splice_csq_del run with set_refalt.
func (s *hapSplice) del(exBeg, exEnd int) int {
	if s.checkStart && s.shiftedDelSynonymous(exBeg, exEnd) {
		s.csq |= csqStartRetained
		return spliceOverlap
	}

	s.refBeg = s.pos + s.tbeg - 1
	s.refEnd = s.pos + s.rlen - s.tend - 1

	if s.refBeg+1 < exBeg {
		if s.checkRegBeg {
			s.spliceBuildHap(exBeg-nSpliceRegionIntron, nSpliceRegionIntron)
			ref, alt := s.kref, s.kalt
			if s.refEnd >= exBeg-nSpliceRegionIntron && s.refBeg < exBeg-nSpliceDonor {
				s.csq |= csqSpliceRegion
				if eqPrefix(ref, alt, nSpliceRegionIntron) {
					s.csq |= csqSynonymous
				}
			}
			if s.refEnd >= exBeg-nSpliceDonor {
				if s.checkDonor && s.rev() {
					s.csq |= csqSpliceDonor
				}
				if s.checkAcceptor && !s.rev() {
					s.csq |= csqSpliceAcceptor
				}
				noff := nSpliceRegionIntron - nSpliceDonor
				if noff < len(s.kref) && noff < len(s.kalt) && eqRange(ref, alt, noff, nSpliceDonor) {
					s.csq |= csqSynonymous
				}
			}
		}
		if s.refEnd >= exBeg {
			s.tbeg = s.refBeg - s.pos + 1
			s.refBeg = exBeg - 1
			if s.tbeg+s.tend == s.alen {
				if s.tend == 0 {
					s.csq |= csqCodingSeq
					return spliceOverlap
				}
				s.tend--
			}
		}
	}
	if exEnd < s.refEnd {
		if s.checkRegEnd {
			s.spliceBuildHap(exEnd+1, nSpliceRegionIntron)
			ref, alt := s.kref, s.kalt
			if s.refBeg < exEnd+nSpliceRegionIntron && s.refEnd > exEnd+nSpliceDonor {
				s.csq |= csqSpliceRegion
				if eqPrefix(ref, alt, nSpliceRegionIntron) {
					s.csq |= csqSynonymous
				}
			}
			if s.refBeg < exEnd+nSpliceDonor {
				if s.checkDonor && !s.rev() {
					s.csq |= csqSpliceDonor
				}
				if s.checkAcceptor && s.rev() {
					s.csq |= csqSpliceAcceptor
				}
				if eqRange(ref, alt, nSpliceRegionIntron-nSpliceDonor, nSpliceDonor) {
					s.csq |= csqSynonymous
				}
			}
		}
		if s.refBeg < exEnd {
			s.tend = s.rlen - (s.refEnd - s.pos + 1)
			s.refEnd = exEnd
		}
	}
	if s.refEnd < exBeg || s.refBeg >= exEnd {
		s.stageSplice()
		return spliceOutside
	}
	if s.refBeg < exBeg+2 {
		if s.checkRegBeg {
			s.csq |= csqSpliceRegion
		}
		if !s.rev() {
			if s.checkStart {
				s.csq |= csqStartLost
			}
		} else if s.checkStop {
			s.csq |= csqStopLost
		}
	}
	if s.refEnd > exEnd-3 {
		if s.checkRegEnd {
			s.csq |= csqSpliceRegion
		}
		if s.rev() {
			if s.checkStart {
				s.csq |= csqStartLost
			}
		} else if s.checkStop {
			s.csq |= csqStopLost
		}
	}
	// set_refalt block.
	if s.tbeg > 0 {
		s.tbeg--
	}
	if s.rlen > s.tbeg+s.tend && s.alen > s.tbeg+s.tend {
		s.rlen -= s.tbeg + s.tend
		s.alen -= s.tbeg + s.tend
	}
	s.kref = s.ref[s.tbeg : s.tbeg+s.rlen]
	s.kalt = s.alt[s.tbeg : s.tbeg+s.alen]
	if (s.refBeg+1 < exBeg && s.refEnd >= exBeg) || (s.refBeg+1 < exEnd && s.refEnd >= exEnd) {
		refBeg := s.refBeg + len(s.kalt) - 1
		if refBeg < s.refEnd {
			if (s.refEnd-refBeg)%3 != 0 {
				s.csq |= csqFrameshift
			} else {
				s.csq |= csqInframeDel
			}
		}
		return spliceOverlap
	}
	s.stageSplice()
	return spliceInside
}

// shiftedDelSynonymous ports shifted_del_synonymous: it detects a
// deletion that, although it touches the start codon, leaves the
// translated start unchanged because the deleted bases are immediately
// re-supplied by identical downstream reference.
func (s *hapSplice) shiftedDelSynonymous(exBeg, exEnd int) bool {
	rev := s.rev()
	if rev && s.pos+s.rlen+2 <= exEnd {
		return false
	}
	if !rev && s.pos >= exBeg+3 {
		return false
	}
	refLen := len(s.ref)
	altLen := len(s.alt)
	if refLen <= altLen {
		return false
	}
	ndel := refLen - altLen
	if rev {
		vcfRefEnd := s.pos + refLen - 1
		trRefEnd := s.ht.end0 + nRefPad
		if vcfRefEnd+ndel > trRefEnd {
			return false
		}
		ptrVcf := s.ref[altLen:]
		ptrRef := s.refSlice(vcfRefEnd+1, len(ptrVcf))
		for i := 0; i < len(ptrVcf); i++ {
			if ptrVcf[i] != ptrRef[i] {
				return false
			}
		}
		return true
	}
	vcfBlockBeg := s.pos + refLen - 2*ndel
	if vcfBlockBeg < 0 {
		return false
	}
	if nRefPad+vcfBlockBeg < exBeg {
		return false
	}
	ptrVcf := s.ref[altLen:]
	ptrRef := s.refSlice(vcfBlockBeg, len(ptrVcf))
	for i := 0; i < len(ptrVcf); i++ {
		if ptrVcf[i] != ptrRef[i] {
			return false
		}
	}
	return true
}

func (s *hapSplice) rev() bool {
	return strandCode(s.ht.tr.Strand) == strandRev
}

// stageSplice ports csq_stage_splice: it pushes the accumulated splice
// consequence bits into the variant buffer. It is a no-op when the
// splice context is not wired to the engine (the hap_init set_refalt
// call) or when no splice bit was set.
func (s *hapSplice) stageSplice() {
	if s.eng == nil || s.recV == nil || s.csq == 0 {
		return
	}
	entry := &csqEntry{pos: s.pos}
	entry.typ.typ = s.csq
	entry.typ.biotype = s.ht.tr.Biotype
	entry.typ.strand = strandCode(s.ht.tr.Strand)
	entry.typ.trid = s.ht.tr.ID
	entry.typ.vcfIal = s.ial
	entry.typ.gene = s.ht.tr.Gene
	// csq_stage_splice -> csq_stage: deduplicate into INFO/BCSQ, then (VCF) set
	// the per-sample FORMAT/BCSQ bit for every haplotype carrying this ALT, or
	// (-O t) stage a text tuple. The shared pushSimple path does exactly this —
	// the splice consequence participates in FORMAT/BCSQ like any other, so the
	// VCF bitmask path must NOT be skipped.
	s.eng.pushSimple(entry, s.recV)
}

// eqPrefix reports whether the first n bytes of a and b are equal,
// mirroring strncmp(a,b,n)==0.
func eqPrefix(a, b string, n int) bool {
	if len(a) < n || len(b) < n {
		return false
	}
	return a[:n] == b[:n]
}

// eqRange reports whether a[off:off+n] equals b[off:off+n].
func eqRange(a, b string, off, n int) bool {
	if off < 0 || len(a) < off+n || len(b) < off+n {
		return false
	}
	return a[off:off+n] == b[off:off+n]
}
