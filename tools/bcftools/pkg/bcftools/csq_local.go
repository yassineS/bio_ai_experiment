// Per-record (non-haplotype-aware) consequence caller for bcftools csq,
// selected by -l/--local-csq. This ports csq.c's test_cds_local: each
// VCF record is annotated independently, without phasing variants onto
// per-sample haplotypes. The coding consequence for a record is derived
// purely from that record's own ref/alt against the spliced reference,
// so compound consequences spanning several records are not called
// jointly (unlike the default haplotype-aware test_cds path).
//
// The UTR / splice / intron / non-coding tests (test_utr / test_splice /
// test_tscript) are identical to the haplotype-aware path and are reused
// unchanged; only the CDS test differs.

package bcftools

import (
	"strconv"
	"strings"
)

// testCDSLocal ports csq.c's test_cds_local. For each coding transcript
// whose CDS overlaps the record, it applies every ALT allele
// independently via hapInit against a throwaway root, translates the
// altered vs reference spliced CDS, classifies the consequence and
// stages it. It returns true when at least one coding transcript
// overlapped the record (mirroring upstream's ret).
func (e *hapEngine) testCDSLocal(rec *hapRecord) bool {
	hit := false
	for _, t := range e.idx.ByChrom[rec.v.Chrom] {
		if !t.Coding {
			continue
		}
		for icds, cds := range t.CDSExons {
			if !overlapsPad(rec.pos, rec.rlen, cds.Start-1, cds.End-1) {
				continue
			}
			// Upstream sets ret=1 the moment a coding CDS overlaps, before
			// the transcript reference is built (csq.c:2966), so a record
			// over a coding region never falls through to test_tscript even
			// if the reference is unavailable. Mirror that ordering here.
			hit = true
			ht := e.getTranscriptForSplice(t)
			if ht == nil {
				continue
			}
			for ial := 1; ial < len(rec.alt); ial++ {
				if rec.alt[ial] == "" || rec.alt[ial][0] == '<' || rec.alt[ial] == "*" {
					continue
				}
				e.cdsLocalAllele(ht, icds, rec, ial)
			}
		}
	}
	return hit
}

// cdsLocalAllele applies one ALT allele to one CDS exon of a transcript
// and stages its consequence, mirroring the per-allele body of
// test_cds_local. It uses a fresh root/child node pair (no haplotype
// tree is retained) so the call is independent of any other record.
func (e *hapEngine) cdsLocalAllele(ht *hapTranscript, icds int, rec *hapRecord, ial int) {
	root := &hapNode{typ: hapRoot}
	node := &hapNode{}
	if e.hapInit(ht, root, node, icds, rec, ial) != 0 {
		return
	}

	tr := ht.tr
	rev := strandCode(tr.Strand) == strandRev

	entry := &csqEntry{pos: rec.pos}
	entry.typ.biotype = tr.Biotype
	entry.typ.strand = strandCode(tr.Strand)
	entry.typ.trid = tr.ID
	entry.typ.vcfIal = ial
	entry.typ.gene = tr.Gene

	csqType := node.csq

	// A start/stop/splice node carries its consequence directly, with no
	// codon prediction. Stage it as-is.
	if node.typ == hapSSS {
		entry.typ.typ = csqType
		e.pushSimple(entry, rec.v)
		return
	}

	sref := ht.sref
	srefM := ht.nsref - 2*nRefPad

	// Translate the altered spliced CDS (node.seq) and, separately, the
	// matching reference slice, both anchored at node.sbeg. Unlike the
	// haplotype path this is a single, self-contained translation per
	// record.
	alt := []byte(node.seq)
	alen := len(alt)
	fill := 0
	if node.dlen%3 != 0 && alen != 0 {
		fill = 1
	}
	strand := strandCode(tr.Strand)
	tseq := cdsTranslate(e.gencode, sref, alt, node.sbeg, node.sbeg, node.sbeg+node.rlen, strand, fill, srefM+node.dlen)

	refSlice := sref[nRefPad+node.sbeg : nRefPad+node.sbeg+node.rlen]
	tref := cdsTranslate(e.gencode, sref, refSlice, node.sbeg, node.sbeg, node.sbeg+node.rlen, strand, fill, srefM)

	// Truncate both translations at their first stop codon, mirroring the
	// truncation loops in test_cds_local (which clip tref/tseq once a '*'
	// stop annotation is seen).
	truncateAtStop(&tref)
	truncateAtStop(&tseq)

	if csqType&csqStopLost != 0 {
		if lastStop(tref) == '*' && lastStop(tref) == lastStop(tseq) {
			csqType &^= csqStopLost
			csqType |= csqStopRetained
		} else if lastStop(tref) != '*' {
			// CDS 3' incomplete: a change to a stop is stop_retained,
			// otherwise the CDS is flagged incomplete.
			if lastStop(tseq) == '*' {
				csqType &^= csqStopGained
				csqType |= csqStopRetained
			} else {
				csqType |= csqIncompleteCDS
			}
		}
	}
	if csqType&csqStartLost != 0 && stopAt(tref, 0) != 'M' {
		csqType &^= csqStartLost
	}

	if node.dlen != 0 {
		if node.dlen%3 != 0 {
			csqType |= csqFrameshift
		} else if node.dlen < 0 {
			csqType |= csqInframeDel
		} else {
			csqType |= csqInframeIns
		}
		if lastStop(tref) != '*' && lastStop(tseq) == '*' {
			csqType |= csqStopGained
		}
		if csqType&csqStartLost != 0 && csqType&csqFrameshift != 0 {
			// Prevent a spurious start_lost|...|1M>1? frameshift call: a
			// start-disrupting frameshift is reduced to a bare start_lost
			// SSS consequence (csq.c #1475227917-adjacent comment).
			csqType &^= csqFrameshift
			entry.typ.typ = csqType
			e.pushSimple(entry, rec.v)
			return
		}
	} else {
		aaChange := false
		for j := 0; j < len(tref.aa) && j < len(tseq.aa); j++ {
			if tref.aa[j] == tseq.aa[j] {
				continue
			}
			aaChange = true
			if stopAt(tref, j) == '*' {
				csqType |= csqStopLost
			} else if stopAt(tseq, j) == '*' {
				csqType |= csqStopGained
			} else {
				csqType |= csqMissense
			}
		}
		if !aaChange {
			csqType |= csqSynonymous
		}
	}

	if csqType&csqCompoundFull != 0 {
		// Build the amino-acid + dna change string. The aa positions are
		// derived from node.sbeg directly (per-record, no compound span).
		var aaRbeg, aaSbeg int
		if !rev {
			aaRbeg = node.sbeg/3 + 1
			aaSbeg = node.sbeg/3 + 1
		} else {
			aaRbeg = (ht.nsref-2*nRefPad-node.sbeg-node.rlen)/3 + 1
			aaSbeg = (ht.nsref-2*nRefPad+node.dlen-node.sbeg-alen)/3 + 1
		}

		var str strings.Builder
		str.WriteByte('|')
		str.WriteString(strconv.Itoa(aaRbeg))
		e.kprintAAPrediction(&str, aaRbeg, tref.aa, tref.stop)
		if csqType&csqSynonymous == 0 {
			str.WriteByte('>')
			str.WriteString(strconv.Itoa(aaSbeg))
			e.kprintAAPrediction(&str, aaSbeg, tseq.aa, tseq.stop)
		}
		str.WriteByte('|')
		str.WriteString(strconv.Itoa(rec.pos + 1))
		str.WriteString(node.variant)

		entry.typ.vstr = str.String()
		entry.typ.hasVstr = true
		entry.typ.typ = csqType & csqCompoundFull
		e.pushSimple(entry, rec.v)
	}

	if csqType&^csqCompoundFull != 0 {
		simple := &csqEntry{pos: rec.pos}
		simple.typ.biotype = tr.Biotype
		simple.typ.strand = strandCode(tr.Strand)
		simple.typ.trid = tr.ID
		simple.typ.vcfIal = ial
		simple.typ.gene = tr.Gene
		simple.typ.typ = csqType &^ csqCompoundFull
		e.pushSimple(simple, rec.v)
	}
}

// truncateAtStop clips a translated sequence at (and including) its
// first stop codon, mirroring the tref/tseq truncation in
// test_cds_local. With no stop the sequence is left unchanged.
func truncateAtStop(t *translated) {
	for j := 0; j < len(t.stop); j++ {
		if t.stop[j] == '*' {
			t.aa = t.aa[:j+1]
			t.stop = t.stop[:j+1]
			return
		}
	}
}

// lastStop returns the stop annotation of the final residue, or 0 when
// the translation is empty.
func lastStop(t translated) byte {
	if len(t.stop) == 0 {
		return 0
	}
	return t.stop[len(t.stop)-1]
}

// stopAt returns the stop annotation at index i, or 0 when i is out of
// range.
func stopAt(t translated, i int) byte {
	if i < 0 || i >= len(t.stop) {
		return 0
	}
	return t.stop[i]
}
