// Per-record consequence classifier for bcftools csq.
//
// This file ports the per-record (single-variant, non-haplotype-phased)
// consequence determination from upstream csq.c: the splice_csq family
// (splice_csq / splice_csq_ins / splice_csq_del / splice_csq_mnp /
// splice_csq_complex) and the test_cds / test_utr / test_splice /
// test_tscript orchestration, together with the SO-term precedence
// ordering used by kput_vcsq.
//
// SCOPE. This is the per-record slice of the multi-slice csq port (see
// docs/PARITY_ROADMAP.md "csq full-parity slicing plan"). It does NOT
// implement the haplotype engine (the hap_node_t tree, cds_translate,
// compound consequences, -p/--phase modes, -n/--ncsq). Output is
// per-record: one BCSQ entry per matching transcript. Haplotype-aware
// compound consequences and byte-for-byte golden parity arrive with the
// next (engine) slice.
//
// SO terms covered: synonymous, missense, stop_gained, stop_lost,
// start_lost, splice_acceptor, splice_donor, splice_region,
// 5_prime_utr, 3_prime_utr, intron, non_coding, inframe_deletion,
// inframe_insertion, frameshift, feature_elongation, feature_truncation,
// coding_sequence.

package bcftools

import (
	"fmt"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/gff"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Consequence-type bits, mirroring the CSQ_* macros in csq.c. The bit
// positions match upstream so the precedence walk in formatConsequence
// can iterate them in the same order kput_vcsq does.
const (
	csqSynonymous     = 1 << 1
	csqMissense       = 1 << 2
	csqStopLost       = 1 << 3
	csqStopGained     = 1 << 4
	csqInframeDel     = 1 << 5
	csqInframeIns     = 1 << 6
	csqFrameshift     = 1 << 7
	csqSpliceAcceptor = 1 << 8
	csqSpliceDonor    = 1 << 9
	csqStartLost      = 1 << 10
	csqSpliceRegion   = 1 << 11
	csqStopRetained   = 1 << 12
	csqUTR5           = 1 << 13
	csqUTR3           = 1 << 14
	csqNonCoding      = 1 << 15
	csqIntron         = 1 << 16
	csqInframeAlter   = 1 << 18
	csqCodingSeq      = 1 << 21
	csqElongation     = 1 << 22
	csqTruncation     = 1 << 23
	csqStartRetained  = 1 << 24
)

// csqStrings maps a consequence bit index to its SO-term string. The
// slot ordering is verbatim from csq.c's csq_strings[]; the precedence
// walk in formatConsequence relies on it.
var csqStrings = [...]string{
	1:  "synonymous",
	2:  "missense",
	3:  "stop_lost",
	4:  "stop_gained",
	5:  "inframe_deletion",
	6:  "inframe_insertion",
	7:  "frameshift",
	8:  "splice_acceptor",
	9:  "splice_donor",
	10: "start_lost",
	11: "splice_region",
	12: "stop_retained",
	13: "5_prime_utr",
	14: "3_prime_utr",
	15: "non_coding",
	16: "intron",
	18: "inframe_altering",
	21: "coding_sequence",
	22: "feature_elongation",
	23: "feature_truncation",
	24: "start_retained",
}

// splice region/donor window sizes, from gff.h.
const (
	nSpliceDonor        = 2 // 2bp at the intron edge -> donor/acceptor
	nSpliceRegionIntron = 8 // up to 8bp into the intron -> splice_region
	nSpliceRegionExon   = 3 // up to 3bp into the exon -> splice_region
)

// splice-csq return codes, mirroring the SPLICE_* macros.
const (
	spliceVarRef  = 0 // ref==alt, not a variant
	spliceOutside = 1 // csq set, does not overlap the coding region
	spliceInside  = 2 // overlaps the coding region, prediction needed
	spliceOverlap = 3 // indel overlaps the region boundary
)

// spliceCSQNonSplice is the set of consequence bits that this port's
// spliceCSQIns / spliceCSQDel set provisionally from raw allele length
// but which upstream's splice_csq_{ins,del} do NOT set at the splice
// layer (upstream recomputes frameshift vs inframe in the haplotype
// engine — see docs/PARITY_ROADMAP.md csq §). The test_splice arm of
// classifyTranscriptVariant masks these off so a CDS-internal indel is
// staged once by the CDS arm and not double-staged as a splice hit.
const spliceCSQNonSplice = csqFrameshift | csqInframeIns | csqInframeDel |
	csqElongation | csqTruncation

// spliceCtx carries the per-(variant,exon) splice state. It is the Go
// equivalent of csq.c's splice_t, trimmed to the fields the per-record
// classifier needs (no haplotype kref/kalt building — those are only
// used by the engine slice for synonymous-at-splice refinement).
//
// All coordinates are 1-based. Upstream csq.c works in 0-based
// coordinates, but every comparison in the splice arithmetic is
// translation-invariant (refBeg < exBeg, refEnd >= exBeg-N, ...), so
// shifting pos and exon coords together by +1 yields identical results
// while keeping consistency with the rest of this package's 1-based
// VCF/GFF handling.
type spliceCtx struct {
	tr   *CSQTranscript
	pos  int    // 1-based VCF position
	rlen int    // REF allele length
	alen int    // ALT allele length
	ref  string // REF allele
	alt  string // ALT allele

	checkAcceptor bool
	checkDonor    bool
	checkStart    bool
	checkStop     bool
	checkRegBeg   bool
	checkRegEnd   bool

	csq    uint32
	tbeg   int // bases trimmed from the start of ref/alt
	tend   int // bases trimmed from the end of ref/alt
	refBeg int // conservative-csq coordinates (1-based)
	refEnd int
}

// classifyCSQRecord returns the per-record BCSQ entries for a variant
// (one per matching transcript). It is the per-record analogue of the
// upstream test_cds + test_utr + test_splice + test_tscript dispatch.
func classifyCSQRecord(v *vcf.Variant, idx *CSQIndex) []string {
	if v == nil || idx == nil {
		return nil
	}
	ref := strings.ToUpper(v.Ref)
	if ref == "" {
		return nil
	}
	alts := make([]string, 0, len(v.Alt))
	for _, a := range v.Alt {
		a = strings.ToUpper(a)
		if a == "" || a == "." || a == "*" || strings.HasPrefix(a, "<") || a == ref {
			continue
		}
		alts = append(alts, a)
	}
	if len(alts) == 0 {
		return nil
	}
	refSeq := idx.Refs[v.Chrom]
	out := []string{}
	for _, t := range idx.ByChrom[v.Chrom] {
		// Off-by-one rlen extension matches upstream (accounts for
		// insertions that sit just past the feature boundary).
		if v.Pos+len(ref) <= t.Beg || v.Pos > t.End+1 {
			continue
		}
		for _, alt := range alts {
			for _, csqBits := range classifyTranscriptVariant(t, refSeq, v.Pos, ref, alt) {
				if csqBits == 0 {
					continue
				}
				out = append(out, formatConsequence(t, csqBits))
			}
		}
	}
	return out
}

// classifyTranscriptVariant computes the per-record consequence entries
// for one (transcript, allele) pair. Each returned bitmask is a SEPARATE
// staged consequence — the Go analogue of one upstream csq_stage call —
// so a single variant can yield several BCSQ entries (e.g. a UTR-exon
// terminal variant stages both a 5'/3'-UTR entry AND a splice_region
// entry). Upstream's per_variant dispatch runs test_cds, test_utr AND
// test_splice unconditionally (csq.c ~3733-3736: `hit += ...`), and only
// falls through to test_tscript when none of them produced a hit; the
// CDS / UTR / splice / transcript walk below mirrors that staging order.
func classifyTranscriptVariant(t *CSQTranscript, refSeq []byte, pos int, ref, alt string) []uint32 {
	var entries []uint32

	// 1. CDS (test_cds): a variant overlapping a CDS exon. For SNPs
	//    this is the codon-level missense/synonymous/stop/start
	//    classifier; for indels it is the splice_csq indel arm.
	for _, e := range t.CDSExons {
		if !overlaps(pos, len(ref), e.Start, e.End) {
			continue
		}
		if c, hit := classifyCDS(t, refSeq, pos, ref, alt, e); hit && c != 0 {
			entries = append(entries, c)
		}
	}

	// 2. UTR (test_utr): a variant overlapping a 5'/3' UTR region.
	//    Each matching UTR stages its own CSQ_UTR5 / CSQ_UTR3 entry,
	//    independent of any splice consequence (csq.c:3442-3454).
	for _, u := range t.UTRs {
		if !overlaps(pos, len(ref), u.Start, u.End) {
			continue
		}
		sc := newSpliceCtx(t, pos, ref, alt)
		ret := spliceCSQ(&sc, u.Start, u.End)
		if ret == spliceInside || ret == spliceOverlap {
			if u.Prime5 {
				entries = append(entries, csqUTR5)
			} else {
				entries = append(entries, csqUTR3)
			}
		}
	}

	// 3. Splice (test_splice): a variant near (but possibly outside) a
	//    full exon — splice_donor / splice_acceptor / splice_region.
	//    Upstream's idx_exon index pads each exon by
	//    N_SPLICE_REGION_INTRON on both sides (gff.c regidx_push), so a
	//    variant up to 8bp into the flanking intron still selects the
	//    exon. test_splice runs unconditionally; a splice hit is staged
	//    as its own entry alongside any CDS/UTR entry above
	//    (csq.c:3461-3488).
	//
	//    Only genuine splice bits are staged here. Upstream's
	//    splice_csq_{ins,del,mnp} never set a frame / elongation /
	//    truncation bit, so test_splice's `if (splice.csq) ret=1` is
	//    driven purely by splice/region/start/stop bits. This port's
	//    spliceCSQIns/Del DO provisionally set a per-record frame bit
	//    (see docs/PARITY_ROADMAP.md csq §), so we mask spliceCSQNonSplice
	//    off before deciding whether test_splice produced a hit — a
	//    CDS-internal indel must not be double-staged by the CDS arm and
	//    again here.
	for _, e := range t.Exons {
		if !overlaps(pos, len(ref)+1, e.Start-nSpliceRegionIntron, e.End+nSpliceRegionIntron) {
			continue
		}
		if !t.Coding {
			continue // splice sites only matter for coding transcripts
		}
		sc := newSpliceCtx(t, pos, ref, alt)
		sc.checkAcceptor, sc.checkDonor = true, true
		sc.checkRegBeg = t.Beg != e.Start
		sc.checkRegEnd = t.End != e.End
		spliceCSQ(&sc, e.Start, e.End)
		if spliceBits := sc.csq &^ spliceCSQNonSplice; spliceBits != 0 {
			entries = append(entries, spliceBits)
		}
	}

	// 4. Transcript fall-through (test_tscript): intron (coding) or
	//    non_coding. Upstream runs this ONLY when test_cds + test_utr +
	//    test_splice all returned no hit (csq.c:3736 `if (!hit)`).
	if len(entries) == 0 {
		sc := newSpliceCtx(t, pos, ref, alt)
		ret := spliceCSQ(&sc, t.Beg, t.End)
		if ret == spliceInside || ret == spliceOverlap {
			if t.Coding {
				entries = append(entries, csqIntron)
			} else {
				entries = append(entries, csqNonCoding)
			}
		}
	}
	return entries
}

// classifyCDS handles a variant that overlaps a CDS exon. The boolean
// return reports whether the variant was fully classified here (true)
// or whether further dispatch is still warranted (false).
func classifyCDS(t *CSQTranscript, refSeq []byte, pos int, ref, alt string, e CSQExon) (uint32, bool) {
	icds := cdsExonIndex(t, e)
	first := t.Strand != gff.StrandReverse && icds == 0 || t.Strand == gff.StrandReverse && icds == len(t.CDSExons)-1
	last := t.Strand != gff.StrandReverse && icds == len(t.CDSExons)-1 || t.Strand == gff.StrandReverse && icds == 0

	sc := newSpliceCtx(t, pos, ref, alt)
	sc.checkAcceptor, sc.checkDonor = true, true
	if first && !t.Trim5 {
		sc.checkStart = true
	}
	if last && !t.Trim3 {
		sc.checkStop = true
	}
	sc.checkRegBeg = icds != 0
	sc.checkRegEnd = icds != len(t.CDSExons)-1

	ret := spliceCSQ(&sc, e.Start, e.End)
	if ret == spliceVarRef {
		return 0, true
	}
	if ret == spliceOutside || ret == spliceOverlap {
		// Splice / start / stop at the CDS boundary; no codon-level
		// prediction. Mirrors hap_init's SPLICE_OUTSIDE/OVERLAP arm.
		return sc.csq, sc.csq != 0
	}

	// ret == spliceInside: the variant sits inside the coding region.
	if len(ref) == 1 && len(alt) == 1 {
		// SNP: codon-level missense / synonymous / stop / start.
		c := classifySNPCodon(t, refSeq, pos, ref[0], alt[0])
		return sc.csq | c, true
	}
	// Indel inside the CDS: splice_csq already decided frameshift vs
	// inframe insertion/deletion via the ref_end-ref_beg %% 3 test.
	return sc.csq, true
}

// classifySNPCodon computes the codon-level consequence for a SNP that
// sits inside a CDS. pos is the 1-based VCF/genomic position used
// directly by the codon machinery (cdsOffset). The bit result feeds
// the precedence-ordered formatter.
func classifySNPCodon(t *CSQTranscript, refSeq []byte, pos int, refBase, altBase byte) uint32 {
	codingOff, ok := cdsOffset(t, pos)
	if !ok || codingOff < 0 {
		return csqCodingSeq
	}
	codonIdx := codingOff / 3
	withinCodon := codingOff % 3
	codonStart := codonIdx * 3

	var codon [3]byte
	for i := 0; i < 3; i++ {
		g, ok := cdsToGenomic(t, codonStart+i)
		if !ok || g < 1 || g > len(refSeq) {
			return csqCodingSeq
		}
		codon[i] = upper(refSeq[g-1])
	}
	if t.Strand == gff.StrandReverse {
		codon = revcompCodon(codon)
		altBase = complementBase(altBase)
		withinCodon = 2 - withinCodon
	}
	mut := codon
	mut[withinCodon] = upper(altBase)

	refAA := translateCodon(codon)
	altAA := translateCodon(mut)

	switch {
	case refAA == '*' && altAA == '*':
		return csqStopRetained
	case refAA == '*' && altAA != '*':
		return csqStopLost
	case altAA == '*' && refAA != '*':
		return csqStopGained
	case codonIdx == 0 && refAA == 'M' && altAA != 'M':
		return csqStartLost
	case codonIdx == 0 && refAA == 'M' && altAA == 'M':
		return csqSynonymous
	case refAA != altAA:
		return csqMissense
	default:
		return csqSynonymous
	}
}

// newSpliceCtx initialises a spliceCtx and computes the trim counts,
// matching the leading work of csq.c's splice_csq.
func newSpliceCtx(t *CSQTranscript, pos int, ref, alt string) spliceCtx {
	sc := spliceCtx{
		tr:   t,
		pos:  pos,
		rlen: len(ref),
		alen: len(alt),
		ref:  ref,
		alt:  alt,
	}
	return sc
}

// spliceCSQ is the Go port of csq.c's splice_csq dispatcher: it trims
// the alleles, then routes to the mnp / complex / ins / del arm.
func spliceCSQ(sc *spliceCtx, exBeg, exEnd int) int {
	sc.alen = len(sc.alt)
	sc.rlen = len(sc.ref)

	rlen1, alen1 := sc.rlen-1, sc.alen-1
	sc.tbeg, sc.tend = 0, 0

	// Trim identical bases from the right, then from the left.
	i := 0
	for i <= rlen1 && i <= alen1 {
		if sc.ref[rlen1-i] != sc.alt[alen1-i] {
			break
		}
		i++
	}
	sc.tend = i
	rlen1 -= i
	alen1 -= i
	i = 0
	for i <= rlen1 && i <= alen1 {
		if sc.ref[i] != sc.alt[i] {
			break
		}
		i++
	}
	sc.tbeg = i

	rtrim := sc.rlen - sc.tbeg - sc.tend
	atrim := sc.alen - sc.tbeg - sc.tend

	switch {
	case sc.rlen == sc.alen:
		return spliceCSQMNP(sc, exBeg, exEnd)
	case rtrim > 1 && atrim > 1:
		return spliceCSQComplex(sc, exBeg, exEnd)
	case sc.rlen < sc.alen:
		return spliceCSQIns(sc, exBeg, exEnd)
	default:
		return spliceCSQDel(sc, exBeg, exEnd)
	}
}

// spliceCSQMNP ports splice_csq_mnp: SNPs and multi-nucleotide
// substitutions (REF and ALT the same length).
func spliceCSQMNP(sc *spliceCtx, exBeg, exEnd int) int {
	if sc.tbeg+sc.tend == sc.rlen {
		return spliceVarRef // not a real variant, eg ACGT>ACGT
	}
	sc.refBeg = sc.pos + sc.tbeg
	sc.refEnd = sc.pos + sc.rlen - sc.tend - 1

	ret := spliceInside
	if sc.refBeg < exBeg { // part before the exon
		if sc.checkRegBeg {
			if sc.refEnd >= exBeg-nSpliceRegionIntron && sc.refBeg < exBeg-nSpliceDonor {
				sc.csq |= csqSpliceRegion
			}
			if sc.refEnd >= exBeg-nSpliceDonor {
				if sc.checkDonor && sc.tr.Strand == gff.StrandReverse {
					sc.csq |= csqSpliceDonor
				}
				if sc.checkAcceptor && sc.tr.Strand != gff.StrandReverse {
					sc.csq |= csqSpliceAcceptor
				}
			}
		}
		if sc.refEnd >= exBeg {
			sc.refBeg = exBeg
			ret = spliceOverlap
		}
	}
	if exEnd < sc.refEnd { // part after the exon
		if sc.checkRegEnd {
			if sc.refBeg <= exEnd+nSpliceRegionIntron && sc.refEnd > exEnd+nSpliceDonor {
				sc.csq |= csqSpliceRegion
			}
			if sc.refBeg <= exEnd+nSpliceDonor {
				if sc.checkDonor && sc.tr.Strand != gff.StrandReverse {
					sc.csq |= csqSpliceDonor
				}
				if sc.checkAcceptor && sc.tr.Strand == gff.StrandReverse {
					sc.csq |= csqSpliceAcceptor
				}
			}
		}
		if sc.refBeg <= exEnd {
			sc.refEnd = exEnd
			ret = spliceOverlap
		}
	}
	if sc.refEnd < exBeg || sc.refBeg > exEnd {
		return spliceOutside
	}
	// Coordinate-invariance note: upstream csq.c is 0-based and tests
	// `ref_beg < ex_beg+3` / `ref_end > ex_end-3` for the first/last 3
	// exon bases (the splice_region exon window, N_SPLICE_REGION_EXON).
	// This port shifts pos and every exon coordinate by +1 into 1-based
	// space; since the literal window width 3 and the strict `<` / `>`
	// comparisons are unchanged, adding +1 to BOTH sides of each
	// inequality cancels — `refBeg < exBeg+3` selects the exact same
	// three bases as upstream. No off-by-one is introduced.
	if sc.refBeg < exBeg+nSpliceRegionExon {
		if sc.checkRegBeg {
			sc.csq |= csqSpliceRegion
		}
		if sc.tr.Strand != gff.StrandReverse {
			if sc.checkStart {
				sc.csq |= csqStartLost
			}
		} else if sc.checkStop {
			sc.csq |= csqStopLost
		}
	}
	if sc.refEnd > exEnd-nSpliceRegionExon {
		if sc.checkRegEnd {
			sc.csq |= csqSpliceRegion
		}
		if sc.tr.Strand == gff.StrandReverse {
			if sc.checkStart {
				sc.csq |= csqStartLost
			}
		} else if sc.checkStop {
			sc.csq |= csqStopLost
		}
	}
	return ret
}

// spliceCSQComplex ports splice_csq_complex: a complex substitution
// where both alleles have >1 trimmed base. Flagged elongation or
// truncation, then handled like an MNP.
func spliceCSQComplex(sc *spliceCtx, exBeg, exEnd int) int {
	if sc.rlen > sc.alen {
		sc.csq |= csqTruncation
	} else {
		sc.csq |= csqElongation
	}
	return spliceCSQMNP(sc, exBeg, exEnd)
}

// spliceCSQIns ports splice_csq_ins: insertions (ALT longer than REF).
func spliceCSQIns(sc *spliceCtx, exBeg, exEnd int) int {
	if sc.tbeg != 0 || sc.ref[0] != sc.alt[0] {
		sc.refBeg = sc.pos + sc.tbeg - 1
		sc.refEnd = sc.pos + sc.rlen - sc.tend
	} else {
		if sc.tend != 0 {
			sc.tend--
		}
		sc.refBeg = sc.pos
		sc.refEnd = sc.pos + sc.rlen - sc.tend
	}

	if sc.refBeg >= exEnd { // fully outside, beyond the exon
		if !sc.checkRegEnd {
			return spliceOutside
		}
		if sc.refBeg < exEnd+nSpliceRegionIntron && sc.refEnd > exEnd+nSpliceDonor {
			sc.csq |= csqSpliceRegion
		}
		if sc.refBeg < exEnd+nSpliceDonor {
			if sc.checkDonor && sc.tr.Strand != gff.StrandReverse {
				sc.csq |= csqSpliceDonor
			}
			if sc.checkAcceptor && sc.tr.Strand == gff.StrandReverse {
				sc.csq |= csqSpliceAcceptor
			}
		}
		return spliceOutside
	}
	if sc.refEnd < exBeg || (sc.refEnd == exBeg && !sc.checkRegBeg) { // fully outside, before the exon
		if !sc.checkRegBeg {
			return spliceOutside
		}
		if sc.refEnd > exBeg-nSpliceRegionIntron && sc.refBeg < exBeg-nSpliceDonor {
			sc.csq |= csqSpliceRegion
		}
		if sc.refEnd > exBeg-nSpliceDonor {
			if sc.checkDonor && sc.tr.Strand == gff.StrandReverse {
				sc.csq |= csqSpliceDonor
			}
			if sc.checkAcceptor && sc.tr.Strand != gff.StrandReverse {
				sc.csq |= csqSpliceAcceptor
			}
		}
		return spliceOutside
	}
	// Overlaps or sits inside the exon.
	if sc.refBeg <= exBeg+2 { // within the first 3bp
		if sc.checkRegBeg {
			sc.csq |= csqSpliceRegion
		}
		if sc.tr.Strand != gff.StrandReverse {
			if sc.checkStart {
				sc.csq |= csqStartLost
			}
		} else if sc.checkStop {
			sc.csq |= csqStopLost
		}
	}
	if sc.refEnd > exEnd-2 {
		if sc.checkRegEnd {
			sc.csq |= csqSpliceRegion
		}
		if sc.tr.Strand == gff.StrandReverse {
			if sc.checkStart {
				sc.csq |= csqStartLost
			}
		} else if sc.checkStop {
			sc.csq |= csqStopLost
		}
	}
	// An insertion inside the CDS is inframe when the inserted length
	// is a multiple of 3, frameshift otherwise (csq.c determines this
	// in the haplotype tree; for a per-record call the inserted-base
	// count is the deciding quantity).
	ins := sc.alen - sc.rlen
	if ins%3 != 0 {
		sc.csq |= csqFrameshift
	} else {
		sc.csq |= csqInframeIns
	}
	return spliceInside
}

// spliceCSQDel ports splice_csq_del: deletions (REF longer than ALT).
func spliceCSQDel(sc *spliceCtx, exBeg, exEnd int) int {
	sc.refBeg = sc.pos + sc.tbeg - 1           // 1bp before the deleted base
	sc.refEnd = sc.pos + sc.rlen - sc.tend - 1 // the last deleted base

	if sc.refBeg+1 < exBeg { // part before the exon
		if sc.checkRegBeg {
			if sc.refEnd >= exBeg-nSpliceRegionIntron && sc.refBeg < exBeg-nSpliceDonor {
				sc.csq |= csqSpliceRegion
			}
			if sc.refEnd >= exBeg-nSpliceDonor {
				if sc.checkDonor && sc.tr.Strand == gff.StrandReverse {
					sc.csq |= csqSpliceDonor
				}
				if sc.checkAcceptor && sc.tr.Strand != gff.StrandReverse {
					sc.csq |= csqSpliceAcceptor
				}
			}
		}
		if sc.refEnd >= exBeg {
			sc.refBeg = exBeg - 1
		}
	}
	if exEnd < sc.refEnd { // part after the exon
		if sc.checkRegEnd {
			if sc.refBeg < exEnd+nSpliceRegionIntron && sc.refEnd > exEnd+nSpliceDonor {
				sc.csq |= csqSpliceRegion
			}
			if sc.refBeg < exEnd+nSpliceDonor {
				if sc.checkDonor && sc.tr.Strand != gff.StrandReverse {
					sc.csq |= csqSpliceDonor
				}
				if sc.checkAcceptor && sc.tr.Strand == gff.StrandReverse {
					sc.csq |= csqSpliceAcceptor
				}
			}
		}
		if sc.refBeg < exEnd {
			sc.refEnd = exEnd
		}
	}
	if sc.refEnd < exBeg || sc.refBeg >= exEnd {
		return spliceOutside
	}
	if sc.refBeg < exBeg+2 { // ref_beg is off by -1
		if sc.checkRegBeg {
			sc.csq |= csqSpliceRegion
		}
		if sc.tr.Strand != gff.StrandReverse {
			if sc.checkStart {
				sc.csq |= csqStartLost
			}
		} else if sc.checkStop {
			sc.csq |= csqStopLost
		}
	}
	if sc.refEnd > exEnd-3 {
		if sc.checkRegEnd {
			sc.csq |= csqSpliceRegion
		}
		if sc.tr.Strand == gff.StrandReverse {
			if sc.checkStart {
				sc.csq |= csqStartLost
			}
		} else if sc.checkStop {
			sc.csq |= csqStopLost
		}
	}
	// A deletion overlapping the coding region is inframe when the
	// deleted span is a multiple of 3, frameshift otherwise. Upstream
	// computes (ref_end - ref_beg) %% 3 over the clipped coordinates.
	if sc.refBeg < sc.refEnd {
		if (sc.refEnd-sc.refBeg)%3 != 0 {
			sc.csq |= csqFrameshift
		} else {
			sc.csq |= csqInframeDel
		}
	} else {
		del := sc.rlen - sc.alen
		if del%3 != 0 {
			sc.csq |= csqFrameshift
		} else {
			sc.csq |= csqInframeDel
		}
	}
	return spliceInside
}

// overlaps reports whether a variant at 1-based pos spanning rlen bases
// touches the 1-based inclusive feature [beg, end].
func overlaps(pos, rlen, beg, end int) bool {
	vbeg := pos            // 1-based first ref base
	vend := pos + rlen - 1 // 1-based last ref base
	if rlen <= 0 {
		vend = vbeg
	}
	return vbeg <= end && vend >= beg
}

// cdsExonIndex returns the index of e within t.CDSExons (genomic
// order), or -1 if not found.
func cdsExonIndex(t *CSQTranscript, e CSQExon) int {
	for i, c := range t.CDSExons {
		if c.Start == e.Start && c.End == e.End {
			return i
		}
	}
	return -1
}

// csqCompound mirrors csq.c's CSQ_COMPOUND macro: the set of
// haplotype-aware consequence bits. Used by CSQ_PRN_STRAND below.
const csqCompound = csqSynonymous | csqMissense | csqStopLost | csqStopGained |
	csqInframeDel | csqInframeIns | csqFrameshift |
	csqStartLost | csqStopRetained | csqInframeAlter |
	csqStartRetained | csqElongation | csqTruncation

// formatConsequence renders one BCSQ entry from a consequence bitmask,
// applying the SO-term precedence ordering of csq.c's kput_vcsq: the
// FIRST set bit (lowest index) becomes the leading term, the rest are
// appended with `&`. It also replicates kput_vcsq's field gating — the
// transcript and strand fields are conditionally omitted (see below).
func formatConsequence(t *CSQTranscript, csq uint32) string {
	// Drop missense when a start/stop change is also present, matching
	// kput_vcsq's "remove missense from start/stops" (csq.c:2152).
	const startStop = csqStopLost | csqStopGained | csqStopRetained | csqStartLost | csqStartRetained
	if csq&startStop != 0 && csq&csqMissense != 0 {
		csq &^= csqMissense
	}

	var sb strings.Builder
	first := true
	for i := 1; i < len(csqStrings); i++ {
		if csqStrings[i] == "" || csq&(1<<uint(i)) == 0 {
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
	consequence := sb.String()
	if consequence == "" {
		consequence = "coding_sequence"
	}

	strand := "+"
	switch t.Strand {
	case gff.StrandReverse:
		strand = "-"
	case gff.StrandUnknown:
		strand = "."
	}

	// Field gating, mirroring kput_vcsq (csq.c:2176-2186):
	//
	//   consequence|gene|transcript|biotype[|strand]
	//
	//  - The transcript ID is printed only when CSQ_PRN_TSCRIPT holds,
	//    i.e. the consequence has some bit set OTHER than CSQ_INTRON /
	//    CSQ_NON_CODING. So a pure intron / non_coding record prints an
	//    EMPTY transcript field (e.g. "intron|GENE||protein_coding").
	//  - The "|±strand" field is printed only when CSQ_PRN_STRAND holds
	//    (a CSQ_COMPOUND bit is set AND no splice/elongation/truncation
	//    bit is set) or when a vstr amino-acid change string exists.
	//    The per-record slice has no vstr, so strand is gated solely on
	//    CSQ_PRN_STRAND. The engine slice adds the vstr suffix and will
	//    print strand whenever a vstr is present.
	transcript := ""
	if csq&^(csqIntron|csqNonCoding) != 0 {
		transcript = t.ID
	}
	prnStrand := csq&csqCompound != 0 &&
		csq&(csqSpliceAcceptor|csqSpliceDonor|csqSpliceRegion|csqElongation|csqTruncation) == 0

	entry := fmt.Sprintf("%s|%s|%s|%s", consequence, t.Gene, transcript, t.Biotype)
	if prnStrand {
		entry += "|" + strand
	}
	return entry
}
