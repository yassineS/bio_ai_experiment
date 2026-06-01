// Top-level processing loop for the bcftools csq haplotype engine: the
// vbuf / pos2vbuf variant buffer, csqStage, the test_cds / test_utr /
// test_splice / test_tscript dispatch and the transcript-boundary
// flushing. Ports csq.c's process / vbuf_push / vbuf_flush / hap_flush.

package bcftools

import (
	"fmt"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/gff"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// hapRecord wraps a VCF variant with the engine's working fields: a
// 0-based position and trimmed/uppercased ref/alt allele strings.
type hapRecord struct {
	v    *vcf.Variant
	pos  int      // 0-based
	rlen int      // len(REF)
	ref  string   // REF, uppercased
	alt  []string // index 0 unused; alt[i] is the i-th ALT
}

// newHapRecord builds a hapRecord from a VCF variant.
func newHapRecord(v *vcf.Variant) *hapRecord {
	ref := upperStr(v.Ref)
	alts := make([]string, len(v.Alt)+1)
	for i, a := range v.Alt {
		alts[i+1] = upperStr(a)
	}
	return &hapRecord{
		v:    v,
		pos:  v.Pos - 1,
		rlen: len(ref),
		ref:  ref,
		alt:  alts,
	}
}

// upperStr uppercases an ASCII string.
func upperStr(s string) string {
	b := []byte(s)
	for i := range b {
		b[i] = upper(b[i])
	}
	return string(b)
}

// process is the per-record entry point. It buffers the record, runs
// the consequence dispatch, then flushes any transcripts and buffered
// records that the new position has moved past.
func (e *hapEngine) process(v *vcf.Variant) error {
	rec := newHapRecord(v)

	// gVCF / no-alt records carry no consequence.
	callCsq := len(rec.alt) > 1
	if callCsq && len(rec.alt) == 2 && (rec.alt[1] == "" || rec.alt[1] == "*" || rec.alt[1] == ".") {
		callCsq = false
	}
	// -i / -e gate consequence calling, not record emission: a record
	// failing -i (or matching -e) is still written out, just without a
	// BCSQ annotation. Ports csq.c:3709-3713.
	if callCsq && e.includeF != nil && !e.includeF.Eval(v) {
		callCsq = false
	}
	if callCsq && e.excludeF != nil && e.excludeF.Eval(v) {
		callCsq = false
	}
	for i := 1; i < len(rec.alt); i++ {
		if rec.alt[i] == "" {
			rec.alt[i] = "."
		}
	}

	if !callCsq {
		e.vbufPush(rec)
		e.hapFlush(rec.pos - 1)
		e.vbufFlush(rec.pos - 1)
		return nil
	}

	if e.rid != "" && e.rid != v.Chrom {
		e.hapFlush(maxPos)
		e.vbufFlush(maxPos)
	}
	e.rid = v.Chrom
	vb := e.vbufPush(rec)

	if !hasSymbolicAlt(rec) {
		hit, err := e.testCDS(rec, vb)
		if err != nil {
			return err
		}
		hit = e.testUTR(rec) || hit
		hit = e.testSplice(rec) || hit
		if !hit {
			e.testTscript(rec)
		}
	}

	if rec.pos > 0 {
		e.hapFlush(rec.pos - 1)
		e.vbufFlush(rec.pos - 1)
	}
	return nil
}

// finish flushes all remaining transcripts and buffered records.
func (e *hapEngine) finish() {
	e.hapFlush(maxPos)
	e.vbufFlush(maxPos)
}

const maxPos = int(^uint(0) >> 1)

// hasSymbolicAlt reports whether the first ALT is a symbolic allele.
func hasSymbolicAlt(rec *hapRecord) bool {
	return len(rec.alt) > 1 && len(rec.alt[1]) > 0 && rec.alt[1][0] == '<'
}

// vbufPush buffers a VCF record, clustering records that share a
// position into one vbuf. Ports vbuf_push.
func (e *hapEngine) vbufPush(rec *hapRecord) *vbuf {
	var vb *vbuf
	if n := len(e.vcfBuf); n > 0 && e.vcfBuf[n-1].vrec[0].rec.Pos == rec.v.Pos {
		vb = e.vcfBuf[n-1]
	} else {
		vb = &vbuf{}
		e.vcfBuf = append(e.vcfBuf, vb)
	}
	vr := &vrecBuf{rec: rec.v}
	vb.vrec = append(vb.vrec, vr)
	e.pos2vbuf[rec.pos] = vb
	return vb
}

// hapFlush finalises transcripts whose end has been passed, running the
// haplotype tree and staging compound consequences. Ports hap_flush.
func (e *hapEngine) hapFlush(pos int) {
	var remaining []*hapTranscript
	// Sort active transcripts by end so the smallest end flushes first.
	sort.SliceStable(e.activeTr, func(i, j int) bool {
		return e.activeTr[i].end0 < e.activeTr[j].end0
	})
	for _, ht := range e.activeTr {
		if ht.end0 > pos {
			remaining = append(remaining, ht)
			continue
		}
		if ht.root != nil && len(ht.root.child) > 0 {
			e.hapFinalize(ht)
			if e.phase != phaseDropGT {
				for i, hdrIdx := range e.samples {
					for j := 0; j < 2; j++ {
						e.hapStageVCF(ht.hap[i*2+j], hdrIdx, j)
					}
				}
			}
		}
	}
	e.activeTr = remaining
}

// hapStageVCF walks a leaf hapNode's csqList and sets the FORMAT/BCSQ
// bits for the given sample/haplotype, mirroring upstream's
// hap_stage_vcf. The leaf-index encoding is hi = 2*ismpl + ihap, so we
// pass (ismpl, ihap) in directly. When icsq2 (= 2*csq.idx+ihap) hits
// the ncsq2 cap we break the loop, just like upstream.
func (e *hapEngine) hapStageVCF(node *hapNode, ismpl, ihap int) {
	if node == nil || len(node.csqList) == 0 || ismpl < 0 {
		return
	}
	for _, csq := range node.csqList {
		if csq.vrec == nil {
			continue
		}
		icsq2 := 2*csq.idx + ihap
		if icsq2 >= e.ncsq2 {
			break
		}
		e.setFmtBit(csq.vrec, ismpl, icsq2)
	}
}

// setFmtBit allocates vrec.fmtBM lazily and sets one (ismpl,icsq2) bit.
func (e *hapEngine) setFmtBit(vrec *vrecBuf, ismpl, icsq2 int) {
	if e.nfmtBcsq <= 0 || len(e.samples) == 0 {
		return
	}
	if vrec.fmtBM == nil {
		// Keyed by header sample index (ismpl), so the stride spans the
		// full header even when -s/-S restricts which samples are
		// processed. Mirrors upstream's calloc(hdr_nsmpl, ...).
		vrec.fmtBM = make([]uint32, e.nsmpl*e.nfmtBcsq)
	}
	ival, ibit := icsq2ToBit(icsq2)
	if 1+ival > vrec.nfmt {
		vrec.nfmt = 1 + ival
	}
	vrec.fmtBM[ismpl*e.nfmtBcsq+ival] |= 1 << uint(ibit)
}

// vbufFlush emits buffered VCF records whose annotations are complete
// (no active transcript still overlaps them). Ports vbuf_flush.
func (e *hapEngine) vbufFlush(pos int) {
	for len(e.vcfBuf) > 0 {
		vb := e.vcfBuf[0]
		if len(e.activeTr) > 0 {
			// Cannot emit a record while a transcript that may still
			// receive variants for it is active.
			if vb.keepUntil > pos {
				break
			}
		}
		e.vcfBuf = e.vcfBuf[1:]
		vbPos := -1
		if len(vb.vrec) > 0 {
			vbPos = vb.vrec[0].rec.Pos - 1
		}
		for _, vr := range vb.vrec {
			if len(vr.vcsq) > 0 {
				var parts []string
				for i := range vr.vcsq {
					parts = append(parts, kputVcsq(&vr.vcsq[i]))
				}
				setInfoField(vr.rec, e.opts.CustomTag, joinComma(parts))
				e.emitFmtBCSQ(vr)
			}
			e.out = append(e.out, vr.rec)
		}
		if vbPos >= 0 {
			delete(e.pos2vbuf, vbPos)
		}
	}
}

// emitFmtBCSQ serialises the per-record FORMAT/BCSQ bitmask into the
// variant's per-sample FORMAT slot. Mirrors upstream's
// bcf_update_format_int32 in vbuf_flush (csq.c:2833-2839): the per-
// sample stride is trimmed to vrec.nfmt (so common short cases emit a
// single int rather than the full nfmtBcsq width), and rows are
// memmoved down accordingly.
func (e *hapEngine) emitFmtBCSQ(vr *vrecBuf) {
	if e.phase == phaseDropGT || len(e.samples) == 0 || vr.fmtBM == nil || vr.nfmt == 0 {
		return
	}
	stride := vr.nfmt
	rec := vr.rec
	if rec == nil {
		return
	}
	// Build per-sample formatted values. Encoded as the raw int32 in
	// decimal, ',' joined when stride > 1, matching how upstream's
	// bcftools query reads FORMAT/BCSQ as a comma-separated int list.
	if len(rec.Samples) < len(e.hdr.Samples) {
		// Pad missing trailing samples (rare for csq inputs but cheap).
		for i := len(rec.Samples); i < len(e.hdr.Samples); i++ {
			rec.Samples = append(rec.Samples, vcf.Sample{Name: e.hdr.Samples[i], Data: map[string]string{}})
		}
	}
	tag := e.opts.CustomTag
	hasTag := false
	for _, f := range rec.Format {
		if f == tag {
			hasTag = true
			break
		}
	}
	if !hasTag {
		rec.Format = append(rec.Format, tag)
	}
	// FORMAT/BCSQ spans every header sample, not just the -s/-S subset:
	// unprocessed samples carry a 0 bitmask (mirrors upstream's
	// bcf_update_format_int32 over the full hdr_nsmpl-wide array).
	for hdrIdx := 0; hdrIdx < e.nsmpl; hdrIdx++ {
		if hdrIdx >= len(rec.Samples) {
			continue
		}
		if rec.Samples[hdrIdx].Data == nil {
			rec.Samples[hdrIdx].Data = map[string]string{}
		}
		var sb strings.Builder
		for k := 0; k < stride; k++ {
			if k > 0 {
				sb.WriteByte(',')
			}
			val := vr.fmtBM[hdrIdx*e.nfmtBcsq+k]
			sb.WriteString(strconv.FormatUint(uint64(val), 10))
		}
		rec.Samples[hdrIdx].Data[tag] = sb.String()
	}
}

// joinComma joins consequence strings with commas.
func joinComma(parts []string) string {
	out := ""
	for i, p := range parts {
		if i > 0 {
			out += ","
		}
		out += p
	}
	return out
}

// stageOneCsq stages a non-haplotype consequence (UTR / splice /
// intron / non-coding) against a record. Ports the csq_stage path for
// these simple consequences: csqPush deduplicates into vrec.vcsq and
// returns the canonical INFO/BCSQ idx; then for each sample whose GT
// carries this ALT, set the corresponding (ismpl,icsq2) bit in
// vrec.fmtBM (mirroring upstream's per-haplotype loop in csq_stage).
func (e *hapEngine) stageOneCsq(rec *hapRecord, typ uint32, t *CSQTranscript, ial int) {
	entry := &csqEntry{pos: rec.pos}
	entry.typ.typ = typ
	entry.typ.biotype = t.Biotype
	entry.typ.strand = strandCode(t.Strand)
	entry.typ.trid = t.ID
	entry.typ.vcfIal = ial
	entry.typ.gene = t.Gene
	e.csqPush(entry, rec.v)
	e.stageSimpleFmtBits(entry, rec.v)
}

// stageSimpleFmtBits sets FORMAT/BCSQ bits for a simple consequence
// (one not produced by the haplotype tree). Mirrors the per-sample
// loop in upstream's csq_stage (csq.c:3382-3419) for VCF output.
func (e *hapEngine) stageSimpleFmtBits(entry *csqEntry, rec *vcf.Variant) {
	if entry.vrec == nil || e.phase == phaseDropGT || len(e.samples) == 0 {
		return
	}
	for _, hdrIdx := range e.samples {
		if hdrIdx >= len(rec.Samples) {
			continue
		}
		gt := rec.Samples[hdrIdx].Data["GT"]
		alleles, _ := csqGTAlleles(gt)
		for j, ial := range alleles {
			if j >= 2 {
				break
			}
			if ial <= 0 || ial != entry.typ.vcfIal {
				continue
			}
			icsq2 := 2*entry.idx + j
			if icsq2 >= e.ncsq2 {
				break
			}
			e.setFmtBit(entry.vrec, hdrIdx, icsq2)
		}
	}
}

// testCDS applies coding variants to the haplotype trees of overlapping
// transcripts. Ports test_cds.
func (e *hapEngine) testCDS(rec *hapRecord, vb *vbuf) (bool, error) {
	hit := false
	for _, t := range e.idx.ByChrom[rec.v.Chrom] {
		if !t.Coding {
			continue
		}
		for icds, cds := range t.CDSExons {
			if !overlapsPad(rec.pos, rec.rlen, cds.Start-1, cds.End-1) {
				continue
			}
			ht := e.getTranscript(t)
			if ht == nil {
				continue
			}
			if ht.end0+1 > vb.keepUntil {
				vb.keepUntil = ht.end0 + 1
			}
			hit = true
			if e.phase == phaseDropGT {
				e.applyCDS(ht, icds, rec)
			} else if err := e.applyCDSSamples(ht, icds, rec); err != nil {
				return hit, err
			}
		}
	}
	return hit, nil
}

// applyCDS extends the haplotype tree for the no-GT (phase drop) case.
func (e *hapEngine) applyCDS(ht *hapTranscript, icds int, rec *hapRecord) {
	if rec.alt[1] == "" || rec.alt[1][0] == '<' || rec.alt[1] == "*" {
		return
	}
	parent := ht.hap[0]
	if parent == nil {
		parent = ht.root
	}
	child := &hapNode{}
	ret := e.hapInit(ht, parent, child, icds, rec, 1)
	if ret != 0 {
		return
	}
	if child.typ == hapSSS {
		entry := &csqEntry{pos: rec.pos}
		entry.typ.typ = child.csq
		entry.typ.biotype = ht.tr.Biotype
		entry.typ.strand = strandCode(ht.tr.Strand)
		entry.typ.trid = ht.tr.ID
		entry.typ.vcfIal = 1
		entry.typ.gene = ht.tr.Gene
		e.csqPush(entry, rec.v)
		return
	}
	parent.nend--
	parent.child = []*hapNode{child}
	ht.hap[0] = child
	child.nend = 1
}

// applyCDSSamples extends the haplotype tree across all samples /
// haplotypes for a record with genotypes. Ports the per-sample loop of
// test_cds. It returns an error for an unphased heterozygous genotype
// under -p require, mirroring upstream csq.c (~line 3274).
func (e *hapEngine) applyCDSSamples(ht *hapTranscript, icds int, rec *hapRecord) error {
	for si, hdrIdx := range e.samples {
		gt := ""
		if hdrIdx < len(rec.v.Samples) {
			gt = rec.v.Samples[hdrIdx].Data["GT"]
		}
		alleles, phased := csqGTAlleles(gt)
		if len(alleles) == 0 || alleles[0] < 0 {
			continue
		}
		ngts := len(alleles)
		if ngts != 1 && ngts != 2 {
			continue
		}
		if ngts > 1 && alleles[1] >= 0 && alleles[0] != alleles[1] {
			if e.phase == phaseMerge {
				if alleles[0] == 0 {
					alleles[0] = alleles[1]
				}
			}
			if !phased {
				switch e.phase {
				case phaseSkip:
					continue
				case phaseNonRef:
					if alleles[0] == 0 {
						alleles[0] = alleles[1]
					} else if alleles[1] == 0 {
						alleles[1] = alleles[0]
					}
				case phaseRequire:
					return fmt.Errorf("Unphased heterozygous genotype at %s:%d, sample %s. See the --phase option.",
						rec.v.Chrom, rec.pos+1, e.hdr.Samples[hdrIdx])
				}
			}
		}
		for ihap := 0; ihap < ngts; ihap++ {
			ial := alleles[ihap]
			if ial < 0 {
				continue
			}
			if ial == 0 {
				continue
			}
			if ial >= len(rec.alt) {
				continue
			}
			if rec.alt[ial] == "" || rec.alt[ial][0] == '<' || rec.alt[ial] == "*" {
				continue
			}
			hi := 2*si + ihap
			parent := ht.hap[hi]
			if parent == nil {
				parent = ht.root
			}
			if parent.curRec == rec.v && ial < len(parent.curChild) && parent.curChild[ial] >= 0 {
				ht.hap[hi] = parent.child[parent.curChild[ial]]
				ht.hap[hi].nend++
				parent.nend--
				continue
			}
			child := &hapNode{}
			ret := e.hapInit(ht, parent, child, icds, rec, ial)
			if ret != 0 {
				continue
			}
			if child.typ == hapSSS {
				entry := &csqEntry{pos: rec.pos}
				entry.typ.typ = child.csq
				entry.typ.biotype = ht.tr.Biotype
				entry.typ.strand = strandCode(ht.tr.Strand)
				entry.typ.trid = ht.tr.ID
				entry.typ.vcfIal = ial
				entry.typ.gene = ht.tr.Gene
				e.csqPush(entry, rec.v)
				e.stageSimpleFmtBits(entry, rec.v)
				continue
			}
			if parent.curRec != rec.v {
				parent.curChild = make([]int, len(rec.alt))
				for j := range parent.curChild {
					parent.curChild[j] = -1
				}
				parent.curRec = rec.v
			}
			j := len(parent.child)
			parent.child = append(parent.child, child)
			parent.curChild[ial] = j
			ht.hap[hi] = child
			child.nend++
			parent.nend--
		}
	}
	return nil
}

// testUTR stages 5'/3' UTR consequences. Ports test_utr.
func (e *hapEngine) testUTR(rec *hapRecord) bool {
	hit := false
	for _, t := range e.idx.ByChrom[rec.v.Chrom] {
		for _, u := range t.UTRs {
			if !overlapsPad(rec.pos, rec.rlen, u.Start-1, u.End-1) {
				continue
			}
			ht := e.getTranscriptForSplice(t)
			for ial := 1; ial < len(rec.alt); ial++ {
				if rec.alt[ial] == "" || rec.alt[ial][0] == '<' || rec.alt[ial] == "*" {
					continue
				}
				s := e.newHapSplice(ht, rec, ial)
				ret := s.run(u.Start-1, u.End-1)
				if ret != spliceInside && ret != spliceOverlap {
					continue
				}
				typ := uint32(csqUTR3)
				if u.Prime5 {
					typ = csqUTR5
				}
				e.stageOneCsq(rec, typ, t, ial)
				hit = true
			}
		}
	}
	return hit
}

// testSplice stages splice consequences for variants near exon edges.
// Ports test_splice.
func (e *hapEngine) testSplice(rec *hapRecord) bool {
	hit := false
	for _, t := range e.idx.ByChrom[rec.v.Chrom] {
		if !t.Coding {
			continue
		}
		for _, ex := range t.Exons {
			if !overlapsPad(rec.pos, rec.rlen, ex.Start-1-nSpliceRegionIntron, ex.End-1+nSpliceRegionIntron) {
				continue
			}
			ht := e.getTranscriptForSplice(t)
			for ial := 1; ial < len(rec.alt); ial++ {
				if rec.alt[ial] == "" || rec.alt[ial][0] == '<' || rec.alt[ial] == "*" {
					continue
				}
				s := e.newHapSplice(ht, rec, ial)
				s.checkAcceptor, s.checkDonor = true, true
				s.checkRegBeg = t.Beg != ex.Start
				s.checkRegEnd = t.End != ex.End
				s.run(ex.Start-1, ex.End-1)
				if s.csq != 0 {
					hit = true
				}
			}
		}
	}
	return hit
}

// testTscript stages intron / non_coding consequences. Ports
// test_tscript.
func (e *hapEngine) testTscript(rec *hapRecord) bool {
	hit := false
	for _, t := range e.idx.ByChrom[rec.v.Chrom] {
		if !overlapsPad(rec.pos, rec.rlen, t.Beg-1, t.End-1) {
			continue
		}
		ht := e.getTranscriptForSplice(t)
		for ial := 1; ial < len(rec.alt); ial++ {
			if rec.alt[ial] == "" || rec.alt[ial][0] == '<' || rec.alt[ial] == "*" {
				continue
			}
			s := e.newHapSplice(ht, rec, ial)
			ret := s.run(t.Beg-1, t.End-1)
			if ret != spliceInside && ret != spliceOverlap {
				continue
			}
			typ := uint32(csqNonCoding)
			if t.Coding {
				typ = csqIntron
			}
			e.stageOneCsq(rec, typ, t, ial)
			hit = true
		}
	}
	return hit
}

// getTranscriptForSplice returns a hapTranscript for the splice/utr/
// intron tests, building it (and its reference) on demand. Unlike the
// CDS path it does not register the transcript as active, since splice
// consequences are staged immediately and never need the tree.
func (e *hapEngine) getTranscriptForSplice(t *CSQTranscript) *hapTranscript {
	if ht, ok := e.hapTr[t.ID]; ok && ht != nil {
		return ht
	}
	refSeq := e.idx.Refs[t.Chrom]
	if len(refSeq) == 0 {
		return nil
	}
	ht := e.getTranscript(t)
	return ht
}

// newHapSplice builds a hapSplice context for the given allele, wired
// to the engine so its splice arms stage consequences via csq_push.
func (e *hapEngine) newHapSplice(ht *hapTranscript, rec *hapRecord, ial int) *hapSplice {
	return &hapSplice{
		ht:   ht,
		pos:  rec.pos,
		ref:  rec.ref,
		alt:  rec.alt[ial],
		rlen: rec.rlen,
		alen: len(rec.alt[ial]),
		eng:  e,
		recV: rec.v,
		ial:  ial,
	}
}

// overlapsPad reports whether a variant at 0-based pos spanning rlen
// bases (with the upstream off-by-one rlen extension for insertions)
// touches the 0-based inclusive feature [beg,end].
func overlapsPad(pos, rlen, beg, end int) bool {
	vbeg := pos
	vend := pos + rlen // deliberate +1 extension, as in regidx_overlap
	return vbeg <= end && vend >= beg
}

// ensure gff import is used.
var _ = gff.StrandForward
