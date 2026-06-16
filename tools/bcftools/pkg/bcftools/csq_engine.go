// Core haplotype-tree logic for the bcftools csq engine: hap_init,
// cds_translate, hap_finalize, hap_add_csq, csq_push and the kput_vcsq
// renderer. See csq_hap.go for the data model and csq_splice.go for the
// set_refalt splice machinery these functions consume.

package bcftools

import (
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// hapInit ports csq.c's hap_init: it applies one VCF allele to a CDS
// exon, builds the child haplotype node's spliced segment, and returns
// a status code (0 added, 1 overlapping variant, 2 silent discard).
func (e *hapEngine) hapInit(ht *hapTranscript, parent, child *hapNode, icds int, rec *hapRecord, ial int) int {
	tr := ht.tr
	cds := tr.CDSExons[icds]
	cdsBeg0 := cds.Start - 1
	cdsLen := cds.End - cds.Start + 1
	child.icds = icds
	child.vcfIal = ial

	s := &hapSplice{
		ht:            ht,
		pos:           rec.pos,
		ref:           rec.ref,
		alt:           rec.alt[ial],
		checkAcceptor: true,
		checkDonor:    true,
		eng:           e,
		recV:          rec.v,
		ial:           ial,
	}
	rev := strandCode(tr.Strand) == strandRev
	if !tr.Trim5 {
		if (!rev && icds == 0) || (rev && icds == len(tr.CDSExons)-1) {
			s.checkStart = true
		}
	}
	if !tr.Trim3 {
		if (!rev && icds == len(tr.CDSExons)-1) || (rev && icds == 0) {
			s.checkStop = true
		}
	}
	if s.checkStart {
		// Do not check starts in incomplete CDS (first codon not M).
		if !rev {
			c := ht.ref[nRefPad+cdsBeg0-ht.beg0:]
			if e.gencode.dna2stop(c[0], c[1], c[2]) != 'M' {
				s.checkStart = false
			}
		} else {
			base := nRefPad + cdsBeg0 - ht.beg0 + cdsLen - 3
			c := ht.ref[base:]
			if e.gencode.cdna2stop(c[0], c[1], c[2]) != 'M' {
				s.checkStart = false
			}
		}
	}
	if icds != 0 {
		s.checkRegBeg = true
	}
	if icds != len(tr.CDSExons)-1 {
		s.checkRegEnd = true
	}

	ret := s.run(cdsBeg0, cdsBeg0+cdsLen-1)

	if ret == spliceVarRef {
		return 2
	}
	if ret == spliceOutside || ret == spliceOverlap {
		if s.csq == 0 {
			return 2
		}
		child.seq = ""
		child.sbeg = 0
		child.rbeg = rec.pos
		child.rlen = 0
		child.dlen = 0
		child.variant = rec.ref + ">" + rec.alt[ial]
		child.typ = hapSSS
		child.csq = s.csq
		child.rec = rec.v
		return 0
	}

	if s.csq&csqSynonymous != 0 {
		s.csq &^= csqSynonymous
	}

	dbeg := 0
	kref := s.kref
	kalt := s.kalt
	refBeg := s.refBeg
	if refBeg < cdsBeg0 {
		dbeg = cdsBeg0 - refBeg
		kref = kref[:len(kref)-dbeg]
		refBeg = cdsBeg0
	}

	var sb strings.Builder
	if parent.typ == hapCDS {
		i := parent.icds
		if i != icds {
			length := tr.CDSExons[i].End - tr.CDSExons[i].Start + 1 - parent.rbeg - parent.rlen + (tr.CDSExons[i].Start - 1)
			if length > 0 {
				sb.Write(ht.ref[nRefPad+parent.rbeg+parent.rlen-ht.beg0 : nRefPad+parent.rbeg+parent.rlen-ht.beg0+length])
			}
		}
		for i++; i < icds; i++ {
			c := tr.CDSExons[i]
			clen := c.End - c.Start + 1
			sb.Write(ht.ref[nRefPad+(c.Start-1)-ht.beg0 : nRefPad+(c.Start-1)-ht.beg0+clen])
		}
		if parent.icds == icds {
			length := refBeg - parent.rbeg - parent.rlen
			if length < 0 {
				return 1 // overlapping variants
			}
			sb.Write(ht.ref[nRefPad+parent.rbeg+parent.rlen-ht.beg0 : nRefPad+parent.rbeg+parent.rlen-ht.beg0+length])
		} else {
			length := refBeg - cdsBeg0
			sb.Write(ht.ref[nRefPad+cdsBeg0-ht.beg0 : nRefPad+cdsBeg0-ht.beg0+length])
		}
	}
	sb.WriteString(kalt[dbeg:])

	child.seq = sb.String()
	child.sbeg = ht.cdsPos(icds) + (refBeg - cdsBeg0)
	child.rbeg = refBeg
	child.rlen = len(kref)
	child.typ = hapCDS
	child.prev = parent
	child.rec = rec.v
	child.csq = s.csq

	rlen := len(rec.ref)
	alen := len(rec.alt[ial])
	child.dlen = alen - rlen
	child.variant = rec.ref + ">" + rec.alt[ial]

	if child.rbeg+child.rlen > cdsBeg0+cdsLen {
		child.typ = hapSSS
		if child.csq == 0 {
			child.csq |= csqCodingSeq
		}
	}
	return 0
}

// translated holds the translation outputs of cdsTranslate.
type translated struct {
	aa   []byte // amino acids
	stop []byte // stop/start annotations ('*','M','-')
}

// cdsTranslate ports csq.c's cds_translate. ref is the spliced
// reference (with nRefPad on each end); seq is the slice of the spliced
// query transcript to translate; sbeg is its offset within the full
// query transcript; rbeg/rend bound seq within the reference; seqM is
// the total query transcript length. fill triggers frameshift fill to
// the transcript end. gc selects the genetic-code table used for codon
// translation.
func cdsTranslate(gc *gencode, ref, seq []byte, sbeg, rbeg, rend, strand int, fill int, seqM int) translated {
	var tseq, tstop []byte
	if len(seq) == 0 {
		return translated{aa: []byte{'?'}, stop: []byte{'?'}}
	}
	var tmp [3]byte

	if strand == strandFwd {
		npad := sbeg % 3
		i := 0
		for ; i < npad; i++ {
			tmp[i] = ref[rbeg+i-npad+nRefPad]
		}
		for ; i < 3 && i-npad < len(seq); i++ {
			tmp[i] = seq[i-npad]
		}
		length := len(seq) - i + npad
		var codon int // index into seq
		if i == 3 {
			tseq = append(tseq, gc.dna2aa(tmp[0], tmp[1], tmp[2]))
			tstop = append(tstop, gc.dna2stop(tmp[0], tmp[1], tmp[2]))
			codon = 3 - npad
			end := codon + length - 1 - (length % 3)
			for codon < end {
				tseq = append(tseq, gc.dna2aa(seq[codon], seq[codon+1], seq[codon+2]))
				tstop = append(tstop, gc.dna2stop(seq[codon], seq[codon+1], seq[codon+2]))
				codon += 3
			}
			i = 0
			for codon+i <= len(seq)-1 {
				tmp[i] = seq[codon+i]
				i++
			}
		}
		// right padding from the reference
		rcodon := rend + nRefPad
		if i > 0 {
			for ; i < 3; i++ {
				tmp[i] = ref[rcodon]
				rcodon++
			}
			tseq = append(tseq, gc.dna2aa(tmp[0], tmp[1], tmp[2]))
			tstop = append(tstop, gc.dna2stop(tmp[0], tmp[1], tmp[2]))
		}
		if fill != 0 {
			end := len(ref) - nRefPad
			for rcodon+3 <= end {
				tseq = append(tseq, gc.dna2aa(ref[rcodon], ref[rcodon+1], ref[rcodon+2]))
				tstop = append(tstop, gc.dna2stop(ref[rcodon], ref[rcodon+1], ref[rcodon+2]))
				rcodon += 3
			}
		}
	} else {
		npad := (seqM - (sbeg + len(seq))) % 3
		var i int
		if npad == 2 {
			tmp[1] = ref[rend+nRefPad]
			tmp[2] = ref[rend+nRefPad+1]
			i = 0
		} else if npad == 1 {
			tmp[2] = ref[rend+nRefPad]
			i = 1
		} else {
			i = 2
		}
		end := len(seq) // index after last seq byte
		for i >= 0 && end > 0 {
			end--
			tmp[i] = seq[end]
			i--
		}
		var codon int // index into seq, start of a codon
		if i == -1 {
			tseq = append(tseq, gc.cdna2aa(tmp[0], tmp[1], tmp[2]))
			tstop = append(tstop, gc.cdna2stop(tmp[0], tmp[1], tmp[2]))
			codon = end - 3
			for codon >= 0 {
				tseq = append(tseq, gc.cdna2aa(seq[codon], seq[codon+1], seq[codon+2]))
				tstop = append(tstop, gc.cdna2stop(seq[codon], seq[codon+1], seq[codon+2]))
				codon -= 3
			}
			// 0 - codon == leftover bases at the start of seq
			switch -codon {
			case 2:
				tmp[2] = seq[0]
				i = 1
			case 1:
				tmp[1] = seq[0]
				tmp[2] = seq[1]
				i = 0
			default:
				i = -1
			}
		}
		// left padding from the reference
		refEnd := nRefPad + rbeg
		if i >= 0 {
			for i >= 0 && refEnd > 0 {
				refEnd--
				tmp[i] = ref[refEnd]
				i--
			}
			tseq = append(tseq, gc.cdna2aa(tmp[0], tmp[1], tmp[2]))
			tstop = append(tstop, gc.cdna2stop(tmp[0], tmp[1], tmp[2]))
		}
		if fill != 0 {
			c := refEnd - 3
			for c >= nRefPad {
				tseq = append(tseq, gc.cdna2aa(ref[c], ref[c+1], ref[c+2]))
				tstop = append(tstop, gc.cdna2stop(ref[c], ref[c+1], ref[c+2]))
				c -= 3
			}
		}
	}
	return translated{aa: tseq, stop: tstop}
}

// hstack is one frame of the hap_finalize DFS stack.
type hstack struct {
	node   *hapNode
	ichild int
	dlen   int
	slen   int
}

// hapWalk carries the mutable hap_finalize / hap_add_csq state, the Go
// analogue of csq.c's hap_t.
type hapWalk struct {
	stack        []hstack
	ht           *hapTranscript
	sseq         []byte
	tseq         translated
	tref         translated
	sbeg         int
	upstreamStop bool
}

// hapFinalize ports csq.c's hap_finalize: a DFS of the haplotype tree
// that breaks each haplotype's spliced sequence into compound parts at
// codon boundaries and calls hapAddCsq for each part.
func (e *hapEngine) hapFinalize(ht *hapTranscript) {
	w := &hapWalk{ht: ht}

	w.stack = []hstack{{node: ht.root, ichild: -1}}
	istack := 0

	for istack >= 0 {
		st := &w.stack[istack]
		node := st.node
		for {
			st.ichild++
			if st.ichild >= len(node.child) {
				break
			}
			if node.child[st.ichild] != nil {
				break
			}
		}
		if st.ichild >= len(node.child) {
			istack--
			continue
		}
		node = node.child[st.ichild]
		prevSlen := w.stack[istack].slen
		prevDlen := w.stack[istack].dlen

		istack++
		if istack >= len(w.stack) {
			w.stack = append(w.stack, hstack{})
		}
		w.stack[istack] = hstack{node: node, ichild: -1}

		// Rebuild the running spliced sequence: truncate to the parent's
		// length, then append this node's CDS segment. Mirrors
		// hap->sseq.l = stack->slen; kputs(node->seq) in hap_finalize.
		w.sseq = w.sseq[:prevSlen]
		if node.typ == hapCDS {
			w.sseq = append(w.sseq, node.seq...)
		}
		w.stack[istack].slen = len(w.sseq)
		w.stack[istack].dlen = prevDlen + node.dlen

		if node.nend == 0 {
			continue
		}
		e.flushHaplotype(w, istack)
	}
}

// flushHaplotype emits the consequences for the leaf haplotype at stack
// depth istack, mirroring the per-strand body of hap_finalize.
func (e *hapEngine) flushHaplotype(w *hapWalk, istack int) {
	ht := w.ht
	tr := ht.tr
	sref := ht.sref
	srefM := len(sref) - 2*nRefPad

	sseqM := srefM + w.stack[istack].dlen
	w.upstreamStop = false

	rev := strandCode(tr.Strand) == strandRev
	w.sbeg = w.stack[1].node.sbeg

	if !rev {
		i := 0
		ibeg := -1
		dlen := 0
		indel := false
		for i < istack {
			i++
			dlen += w.stack[i].node.dlen
			if w.stack[i].node.dlen != 0 {
				indel = true
			}
			if i < istack {
				if dlen%3 != 0 {
					if ibeg == -1 {
						ibeg = i
					}
					continue
				}
				icur := w.node2sbeg(i)
				inext := w.node2sbeg(i + 1)
				if w.stack[i].node.dlen > 0 {
					icur += w.stack[i].node.dlen
				} else if w.stack[i].node.dlen < 0 {
					icur++
				}
				if icur/3 == inext/3 {
					if ibeg == -1 {
						ibeg = i
					}
					continue
				}
			}
			if ibeg < 0 {
				ibeg = i
			}
			ioff := w.node2soff(ibeg)
			icur := w.node2sbeg(ibeg)
			rbeg := w.node2rbeg(ibeg)
			rend := w.node2rend(i)
			fill := dlen % 3

			var altSeq []byte
			if len(w.sseq) != 0 {
				altSeq = w.sseq[ioff:w.stack[i].slen]
			} else {
				fill = 0
			}
			w.tseq = cdsTranslate(e.gencode, sref, altSeq, icur, rbeg, rend, strandFwd, fill, sseqM)
			refSeq := sref[nRefPad+rbeg : nRefPad+w.node2rend(i)]
			w.tref = cdsTranslate(e.gencode, sref, refSeq, rbeg, rbeg, rend, strandFwd, fill, srefM)

			e.hapAddCsq(w, w.stack[istack].node, 0, ibeg, i, dlen, indel)
			ibeg = -1
			dlen = 0
			indel = false
		}
	} else {
		i := istack + 1
		ibeg := -1
		dlen := 0
		indel := false
		for i > 1 {
			i--
			dlen += w.stack[i].node.dlen
			if w.stack[i].node.dlen != 0 {
				indel = true
			}
			if i > 1 {
				if dlen%3 != 0 {
					if ibeg == -1 {
						ibeg = i
					}
					continue
				}
				icur := sseqM - 1 - w.node2sbeg(i)
				inext := sseqM - 1 - w.node2sbeg(i-1)
				if w.stack[i].node.dlen > 0 {
					icur += w.stack[i].node.dlen - 1
				} else if w.stack[i].node.dlen < 0 {
					icur -= w.stack[i].node.dlen
				}
				if w.stack[i-1].node.dlen > 0 {
					inext -= w.stack[i-1].node.dlen
				}
				if icur/3 == inext/3 {
					if ibeg == -1 {
						ibeg = i
					}
					continue
				}
			}
			if ibeg < 0 {
				ibeg = i
			}
			ioff := w.node2soff(i)
			icur := w.node2sbeg(i)
			rbeg := w.node2rbeg(i)
			rend := w.node2rend(ibeg)
			fill := dlen % 3

			var altSeq []byte
			if len(w.sseq) != 0 {
				altSeq = w.sseq[ioff:w.stack[ibeg].slen]
			} else {
				fill = 0
			}
			w.tseq = cdsTranslate(e.gencode, sref, altSeq, icur, rbeg, rend, strandRev, fill, sseqM)
			refSeq := sref[nRefPad+rbeg : nRefPad+w.node2rend(ibeg)]
			w.tref = cdsTranslate(e.gencode, sref, refSeq, rbeg, rbeg, rend, strandRev, fill, srefM)

			e.hapAddCsq(w, w.stack[istack].node, sseqM, i, ibeg, dlen, indel)
			ibeg = -1
			dlen = 0
			indel = false
		}
	}
}

// node2* helpers mirror the macros in csq.c around hap_add_csq.
func (w *hapWalk) node2soff(i int) int {
	return w.stack[i].slen - (w.stack[i].node.rlen + w.stack[i].node.dlen)
}
func (w *hapWalk) node2sbeg(i int) int { return w.sbeg + w.node2soff(i) }
func (w *hapWalk) node2send(i int) int { return w.sbeg + w.stack[i].slen }
func (w *hapWalk) node2rbeg(i int) int { return w.stack[i].node.sbeg }
func (w *hapWalk) node2rend(i int) int { return w.stack[i].node.sbeg + w.stack[i].node.rlen }
func (w *hapWalk) node2rpos(i int) int { return w.stack[i].node.rbeg }

// hapAddCsq ports csq.c's hap_add_csq: given a compound variant group
// [ibeg,iend], it diffs the translated reference and haplotype
// sequences, determines the true consequence bits, builds the aa+dna
// variant string and stages the csq (plus any @pos back-references).
func (e *hapEngine) hapAddCsq(w *hapWalk, node *hapNode, tlen, ibeg, iend, dlen int, indel bool) {
	ht := w.ht
	tr := ht.tr
	rev := strandCode(tr.Strand) == strandRev
	refNode := ibeg
	if rev {
		refNode = iend
	}

	entry := &csqEntry{}
	entry.pos = w.stack[refNode].node.rbeg
	entry.typ.trid = tr.ID
	entry.typ.vcfIal = node.vcfIal
	entry.typ.gene = tr.Gene
	entry.typ.strand = strandCode(tr.Strand)
	entry.typ.biotype = tr.Biotype
	node.csqList = append(node.csqList, entry)

	var rmCsq uint32
	var typ uint32
	for i := ibeg; i <= iend; i++ {
		typ |= w.stack[i].node.csq & csqCompoundFull
	}
	if dlen == 0 && indel {
		typ |= csqInframeAlter
	}

	hasUpstreamStop := w.upstreamStop
	tref := &w.tref
	tseq := &w.tseq
	if w.stack[ibeg].node.typ != hapSSS {
		// truncate at the first ref stop
		for i := 0; i < len(tref.stop); i++ {
			if tref.stop[i] == '*' {
				tref.aa = tref.aa[:i+1]
				tref.stop = tref.stop[:i+1]
				break
			}
		}
		for i := 0; i < len(tseq.stop); i++ {
			if tseq.stop[i] == '*' {
				tseq.aa = tseq.aa[:i+1]
				tseq.stop = tseq.stop[:i+1]
				w.upstreamStop = true
				break
			}
		}
		if typ&csqStopLost != 0 {
			if tref.stop[len(tref.stop)-1] == '*' && tref.stop[len(tref.stop)-1] == tseq.stop[len(tseq.stop)-1] {
				rmCsq |= csqStopLost
				typ |= csqStopRetained
			} else if tref.stop[len(tref.stop)-1] != '*' {
				if tseq.stop[len(tseq.stop)-1] == '*' {
					rmCsq |= csqStopGained
					typ |= csqStopRetained
				} else {
					typ |= csqIncompleteCDS
				}
			}
		}
		if typ&csqStartLost != 0 {
			if tref.stop[len(tref.stop)-1] == 'M' && tref.stop[len(tref.stop)-1] == tseq.stop[len(tseq.stop)-1] {
				rmCsq |= csqStartLost
				typ |= csqStartRetained
			}
		}
		if dlen != 0 {
			if dlen%3 != 0 {
				typ |= csqFrameshift
			} else if dlen < 0 {
				typ |= csqInframeDel
			} else {
				typ |= csqInframeIns
			}
			if tref.stop[len(tref.stop)-1] != '*' && tseq.stop[len(tseq.stop)-1] == '*' {
				typ |= csqStopGained
			}
		} else {
			aaChange := false
			for i := 0; i < len(tref.aa) && i < len(tseq.aa); i++ {
				if tref.aa[i] == tseq.aa[i] {
					continue
				}
				aaChange = true
				if tref.stop[i] == '*' {
					typ |= csqStopLost
				} else if tseq.stop[i] == '*' {
					typ |= csqStopGained
				} else {
					typ |= csqMissense
				}
			}
			if !aaChange {
				typ |= csqSynonymous
			}
		}
	}
	suppressAA := false
	if ibeg != iend && typ&(csqInframeDel|csqInframeIns|csqInframeAlter) != 0 && tseq.stop[len(tseq.stop)-1] == '*' {
		rmCsq |= csqInframeDel | csqInframeIns | csqInframeAlter
		typ |= csqFrameshift | csqStopGained
	}
	if typ&csqFrameshift != 0 && typ&csqStartLost != 0 {
		rmCsq |= csqFrameshift
		suppressAA = true
		if ibeg == iend {
			w.stack[ibeg].node.typ = hapSSS
		}
	}
	if hasUpstreamStop {
		typ |= csqUpstreamStop
	}
	typ &^= rmCsq

	if w.stack[ibeg].node.typ == hapSSS {
		entry.typ.typ |= typ | (w.stack[ibeg].node.csq &^ rmCsq)
		entry.typ.ref = w.stack[ibeg].node.rec
		entry.typ.biotype = tr.Biotype
		e.csqPush(entry, w.stack[ibeg].node.rec)
		return
	}

	entry.typ.typ = typ

	var str strings.Builder
	if suppressAA {
		str.WriteString("||")
	} else {
		var aaRbeg, aaSbeg int
		if !rev {
			aaRbeg = w.node2rbeg(ibeg)/3 + 1
			aaSbeg = w.node2sbeg(ibeg)/3 + 1
		} else {
			aaRbeg = (ht.nsref-2*nRefPad-w.node2rend(iend))/3 + 1
			aaSbeg = (tlen-w.node2send(iend))/3 + 1
		}
		str.WriteByte('|')
		str.WriteString(strconv.Itoa(aaRbeg))
		e.kprintAAPrediction(&str, aaRbeg, tref.aa, tref.stop)
		if typ&csqSynonymous == 0 {
			str.WriteByte('>')
			str.WriteString(strconv.Itoa(aaSbeg))
			e.kprintAAPrediction(&str, aaSbeg, tseq.aa, tseq.stop)
		}
		str.WriteByte('|')
	}
	for i := ibeg; i <= iend; i++ {
		if i > ibeg {
			str.WriteByte('+')
		}
		str.WriteString(strconv.Itoa(w.node2rpos(i) + 1))
		str.WriteString(w.stack[i].node.variant)
	}
	entry.typ.vstr = str.String()
	entry.typ.hasVstr = true
	e.csqPush(entry, w.stack[refNode].node.rec)

	dnaStr := str.String()
	for i := ibeg; i <= iend; i++ {
		if w.stack[i].node.csq&^csqCompoundFull != 0 {
			tmp := &csqEntry{pos: w.stack[i].node.rbeg}
			tmp.typ.trid = tr.ID
			tmp.typ.gene = tr.Gene
			tmp.typ.strand = strandCode(tr.Strand)
			tmp.typ.typ = w.stack[i].node.csq &^ csqCompoundFull &^ rmCsq
			tmp.typ.biotype = tr.Biotype
			tmp.typ.vstr = dnaStr
			tmp.typ.hasVstr = true
			node.csqList = append(node.csqList, tmp)
			e.csqPush(tmp, w.stack[i].node.rec)
		}
		if i != refNode && (entry.typ.typ&csqCompoundFull != 0 || w.stack[i].node.csq&^csqCompoundFull == 0) {
			tmp := &csqEntry{pos: w.stack[i].node.rbeg}
			tmp.typ.trid = tr.ID
			tmp.typ.gene = tr.Gene
			tmp.typ.strand = strandCode(tr.Strand)
			tmp.typ.typ = csqPrintedUpstream | w.stack[i].node.csq
			tmp.typ.biotype = tr.Biotype
			tmp.typ.ref = w.stack[refNode].node.rec
			node.csqList = append(node.csqList, tmp)
			e.csqPush(tmp, w.stack[i].node.rec)
		}
	}
}

// csqPush ports csq.c's csq_push: it inserts a consequence into the
// matching buffered VCF record, applying the upstream deduplication and
// merge rules. It returns existed == true when the consequence merged
// into (or matched) an already-present entry, and false when a new entry
// was appended — mirroring upstream's non-zero/zero return. When the
// target buffer or record cannot be found, existed is reported true so
// callers treat it as a no-op skip (upstream's -1 path).
func (e *hapEngine) csqPush(csq *csqEntry, rec *vcf.Variant) (existed bool) {
	vb := e.pos2vbuf[csq.pos]
	if vb == nil {
		return true
	}
	if csq.typ.typ&csqInframeIns != 0 && csq.typ.typ&csqElongation != 0 {
		csq.typ.typ &^= csqInframeIns
	}
	if csq.typ.typ&csqInframeDel != 0 && csq.typ.typ&csqTruncation != 0 {
		csq.typ.typ &^= csqInframeDel
	}
	var vrec *vrecBuf
	for _, vr := range vb.vrec {
		if vr.rec == rec {
			vrec = vr
			break
		}
	}
	if vrec == nil {
		return true
	}
	if csq.typ.typ&csqSpliceRegion != 0 && csq.typ.typ&(csqSpliceDonor|csqSpliceAcceptor) != 0 {
		csq.typ.typ &^= csqSpliceRegion
	}

	const prnTscript = ^uint32(csqIntron | csqNonCoding)

	if csq.typ.typ&csqPrintedUpstream != 0 {
		for i := range vrec.vcsq {
			if csq.typ.typ&csqStartStop != 0 && vrec.vcsq[i].typ&csqStartStop != 0 {
				vrec.vcsq[i] = csq.typ
				csq.vrec, csq.idx = vrec, i
				return true
			}
			if vrec.vcsq[i].typ&csqPrintedUpstream == 0 {
				continue
			}
			if csq.typ.ref != vrec.vcsq[i].ref {
				continue
			}
			csq.vrec, csq.idx = vrec, i
			return true
		}
	} else if csq.typ.typ&csqCompoundFull != 0 {
		for i := range vrec.vcsq {
			vc := &vrec.vcsq[i]
			if csq.typ.trid != vc.trid && (csq.typ.typ|vc.typ)&prnTscript != 0 {
				continue
			}
			if csq.typ.biotype != vc.biotype {
				continue
			}
			if csq.typ.gene != vc.gene {
				continue
			}
			if csq.typ.vcfIal != vc.vcfIal {
				continue
			}
			if (csq.typ.typ&csqUpstreamStop)^(vc.typ&csqUpstreamStop) != 0 {
				continue
			}
			if csq.typ.hasVstr || vc.hasVstr {
				if !csq.typ.hasVstr || !vc.hasVstr {
					if csq.typ.typ&csqStartStop != 0 && vc.typ&csqStartStop != 0 {
						vc.typ |= csq.typ.typ
						if vc.typ&csqStopRetained != 0 {
							vc.typ &^= csqStopLost | csqSynonymous
						}
						if vc.typ&csqStartRetained != 0 {
							vc.typ &^= csqStartLost | csqSynonymous
						}
						if !vc.hasVstr {
							vc.vstr = csq.typ.vstr
							vc.hasVstr = csq.typ.hasVstr
						}
						csq.vrec, csq.idx = vrec, i
						return true
					}
					continue
				}
				if csq.typ.vstr != vc.vstr {
					continue
				}
			}
			vc.typ |= csq.typ.typ
			csq.vrec, csq.idx = vrec, i
			return true
		}
	} else {
		for i := range vrec.vcsq {
			vc := &vrec.vcsq[i]
			if csq.typ.trid != vc.trid && (csq.typ.typ|vc.typ)&prnTscript != 0 {
				continue
			}
			if csq.typ.biotype != vc.biotype {
				continue
			}
			if vc.typ&csqCompoundFull == 0 {
				vc.typ |= csq.typ.typ
				csq.vrec, csq.idx = vrec, i
				return true
			}
			if vc.typ == (vc.typ | csq.typ.typ) {
				csq.vrec, csq.idx = vrec, i
				return true
			}
		}
	}
	csq.vrec = vrec
	csq.idx = len(vrec.vcsq)
	vrec.vcsq = append(vrec.vcsq, csq.typ)
	return false
}

// kprintAAPrediction ports csq.c's kprint_aa_prediction: it appends an
// amino-acid prediction to str. With brief predictions disabled (or a
// prediction too short to abbreviate) the full amino-acid string is
// written. Otherwise only the first e.brief residues are written,
// followed by ".." and the 1-based residue index just past the
// prediction (dropping a trailing stop codon from the length, as
// upstream does).
func (e *hapEngine) kprintAAPrediction(str *strings.Builder, beg int, aa, stop []byte) {
	if e.brief == 0 || len(aa)-e.brief < 3 {
		str.Write(aa)
		return
	}
	length := len(aa)
	if length > 0 && stop[length-1] == '*' {
		length--
	}
	for i := 0; i < length && i < e.brief; i++ {
		str.WriteByte(aa[i])
	}
	str.WriteString("..")
	str.WriteString(strconv.Itoa(beg + length))
}

// kputVcsq ports csq.c's kput_vcsq: it renders one vcsq into the
// pipe-delimited BCSQ form.
func kputVcsq(c *vcsq) string {
	typ := c.typ

	if typ&csqIncompleteCDS != 0 && typ&^(csqStartStop|csqIncompleteCDS|csqUpstreamStop) != 0 {
		typ &^= csqStartStop | csqIncompleteCDS
	}
	if typ&csqStartStop != 0 && typ&csqMissense != 0 {
		typ &^= csqMissense
	}

	var sb strings.Builder
	if typ&csqPrintedUpstream != 0 && c.ref != nil {
		sb.WriteByte('@')
		sb.WriteString(strconv.Itoa(c.ref.Pos))
		return sb.String()
	}
	if typ&csqUpstreamStop != 0 {
		sb.WriteByte('*')
	}

	first := true
	for i := 1; i < len(csqStrings); i++ {
		if csqStrings[i] == "" || typ&(1<<uint(i)) == 0 {
			continue
		}
		if first {
			sb.WriteString(csqStrings[i])
			first = false
		} else {
			sb.WriteByte('&')
			sb.WriteString(csqStrings[i])
		}
	}

	sb.WriteByte('|')
	sb.WriteString(c.gene)

	sb.WriteByte('|')
	const prnTscript = ^uint32(csqIntron | csqNonCoding)
	if typ&prnTscript != 0 {
		sb.WriteString(c.trid)
	}

	sb.WriteByte('|')
	sb.WriteString(c.biotype)

	prnStrand := typ&csqCompoundFull != 0 &&
		typ&(csqSpliceAcceptor|csqSpliceDonor|csqSpliceRegion|csqElongation|csqTruncation) == 0
	if prnStrand || c.hasVstr {
		switch c.strand {
		case strandFwd:
			sb.WriteString("|+")
		case strandRev:
			sb.WriteString("|-")
		default:
			sb.WriteString("|.")
		}
	}
	if c.hasVstr {
		sb.WriteString(c.vstr)
	}
	return sb.String()
}
