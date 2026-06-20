package samtools

import (
	"bufio"
	"fmt"
	"math"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/errmod"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// upstreamPhaseRunner is the state struct that mirrors phase.c's
// `phaseg_t`. Configuration fields are populated from PhaseOptions in
// runUpstreamPhase; runtime state is updated by the streaming loop.
type upstreamPhaseRunner struct {
	// configuration
	k          int
	minBaseQ   uint8
	minVarLOD  int
	maxDepth   int
	fixChimera bool

	// Optional site filter (CLI -l/-e). When siteSet is non-nil a
	// position is in_set iff its (chrom, 1-based pos) appears in
	// the set; listExcl gates whether out-of-set positions are
	// dropped entirely (FLAG_LIST_EXCL).
	siteSet  map[string]map[int]struct{}
	listExcl bool

	// runtime state
	vposShift int // total het sites emitted on the current chromosome
}

// runUpstreamPhase drives the byte-faithful phase emit pipeline on a
// pre-grouped slice of records. The records must all be on a single
// reference, sorted by Pos. On return, all PS/FL/M/EV/// lines for
// this reference will have been written to bw in upstream order. The
// CC banner is emitted once for the whole stream, not per-reference;
// it is the caller's responsibility (see Phase).
//
// When bs is non-nil, after each block's phaseEmit completes (and
// once more at the end of the reference) the per-read dump_aln
// routing is invoked, dispatching each record into one of the three
// BAM outputs. This matches upstream phase.c's interleaving of
// phase()+dump_aln. rng is consumed for the evidence-less and is_flip
// branches in dump_aln. opts is used for opts.DropAmbiguous
// (FLAG_DROP_AMBI).
//
// Returns the number of het sites emitted on this reference.
func runUpstreamPhase(g *upstreamPhaseRunner, hash *fragKhash, recs []*sam.Record, rname string, isLastRef bool, bw *bufio.Writer, bs *bamSplitWriter, rng phaseRNG, opts PhaseOptions) (int, error) {
	pp := newPhaseStreamPileup(recs, rname)
	em := errmod.Init(1.0 - 0.83)
	bases := make([]uint16, 0, g.maxDepth)
	cns := make([]uint64, 0, 256)
	vpos := 0
	// g.vposShift is NOT reset here: upstream resets vpos_shift to 0 at a tid
	// change *before* flushing the previous chromosome's trailing buffer (so
	// that buffer's hets are renumbered from 0 and the new chromosome continues
	// from there). The reset is therefore applied just before this reference's
	// final flush below, not at its start — except for the last reference,
	// whose trailing buffer upstream flushes at end-of-stream with no reset.
	emitted := 0
	q := make([]float32, 16)
	// cursor over `recs` — head of the per-reference queue not yet
	// passed to dump_aln. Only used when bs is non-nil.
	cursor := 0

	for {
		col := pp.next()
		if col == nil {
			break
		}
		pos0 := int(col.pos) - 1
		n := len(col.plp)
		if n > g.maxDepth {
			continue
		}
		// Fill bases array (phase.c:738-752).
		bases = bases[:0]
		for i := 0; i < n; i++ {
			p := &col.plp[i]
			if p.isDel || p.isRefSkip {
				continue
			}
			if int(p.qpos) >= len(p.b.Qual) || int(p.qpos) >= len(p.b.Seq) {
				continue
			}
			packed, ok := admitPhaseBase(
				p.b.Seq[p.qpos], p.b.Qual[p.qpos], p.b.MapQ,
				p.b.Flag&sam.FlagReverse != 0, g.minBaseQ)
			if !ok {
				continue
			}
			bases = append(bases, packed)
		}
		if len(bases) == 0 {
			continue
		}
		for i := range q {
			q[i] = 0
		}
		em.Cal(bases, 4, q)
		c := gl2cns(q)
		// Site-list gating (phase.c:757-758). When the list is
		// non-nil, a position is "in_set" iff its (chrom, 1-based
		// pos) appears there. -e (listExcl) drops positions
		// outside the list entirely; otherwise the normal LOD
		// test is bypassed for in-set positions so the user can
		// force phasing at chosen sites.
		inSet := false
		if g.siteSet != nil {
			if perChrom, ok := g.siteSet[rname]; ok {
				if _, ok := perChrom[pos0+1]; ok {
					inSet = true
				}
			}
			if g.listExcl && !inSet {
				continue
			}
		}
		// Not a variant? (phase.c:758)
		if !isPhaseVariantColumn(c, g.minVarLOD, inSet) {
			continue
		}
		// Push variant to cns.
		if vpos < len(cns) {
			cns[vpos] = uint64(pos0)<<32 | uint64(c)
		} else {
			cns = append(cns, uint64(pos0)<<32|uint64(c))
		}
		// Walk reads to populate frag hash.
		dophase := true
		for i := 0; i < n; i++ {
			p := &col.plp[i]
			if p.isDel || p.isRefSkip {
				continue
			}
			if p.b.MapQ == 0 {
				continue
			}
			if int(p.qpos) >= len(p.b.Seq) {
				continue
			}
			b, ok := nt16ToInt(p.b.Seq[p.qpos])
			if !ok {
				continue
			}
			var ac int8
			cur := cns[vpos]
			if uint32(b) == uint32(cur&3) {
				ac = 1
			} else if uint32(b) == uint32(cur>>16&3) {
				ac = 2
			} else {
				ac = 0
			}
			key := x31HashString(p.b.QName)
			kk, isNew := hash.put(key)
			f := &hash.vals[kk]
			if !isNew {
				if vpos-int(f.vpos)+1 < maxVars {
					f.vlen = uint16(vpos - int(f.vpos) + 1)
					f.seq[f.vlen-1] = ac
					f.end = int32(p.b.EndPosition())
				}
				dophase = false
			} else {
				f.beg = int32(p.b.Pos) - 1
				f.end = int32(p.b.EndPosition())
				f.vpos = int32(vpos)
				f.vlen = 1
				f.seq[0] = ac
				f.single = 0
				f.phased = 0
				f.flip = 0
				f.ambig = 0
			}
		}
		if dophase {
			n2, err := phaseEmit(g, rname, vpos, cns, hash, bw)
			if err != nil {
				return emitted, err
			}
			emitted += n2
			// dump_aln (phase.c:484) — drain queue front by current
			// hash state, before updateVpos slides the indices. min_pos
			// is cns[vpos]>>32 when a frag with vpos>=block_vpos exists
			// (phase.c:411), else MaxInt32 to flush all queued reads.
			if bs != nil {
				minPos := int32(math.MaxInt32)
				if hashHasFragSpanning(vpos, hash) {
					minPos = int32(cns[vpos] >> 32)
				}
				cursor, err = dumpAln(recs, cursor, minPos, hash, bs, rng, opts.DropAmbiguous)
				if err != nil {
					return emitted, err
				}
			}
			updateVpos(vpos, hash)
			cns[0] = cns[vpos]
			vpos = 0
		}
		vpos++
	}
	// Final flush for any remaining vpos > 0. Upstream resets vpos_shift to 0
	// at the next chromosome's arrival, just before this flush — so the
	// trailing block is renumbered from 0 and the next reference continues from
	// there. The last reference has no following chromosome, so it keeps the
	// accumulated counter (upstream's end-of-stream phase() with no reset).
	if !isLastRef {
		g.vposShift = 0
	}
	if vpos > 0 {
		n2, err := phaseEmit(g, rname, vpos, cns, hash, bw)
		if err != nil {
			return emitted, err
		}
		emitted += n2
		if bs != nil {
			minPos := int32(math.MaxInt32)
			if hashHasFragSpanning(vpos, hash) {
				minPos = int32(cns[vpos] >> 32)
			}
			cursor, err = dumpAln(recs, cursor, minPos, hash, bs, rng, opts.DropAmbiguous)
			if err != nil {
				return emitted, err
			}
		}
	}
	// Final tail flush: any remaining queued reads beyond the last
	// block on this reference go via dump_aln with min_pos = INT_MAX,
	// matching upstream's after-loop drain (samtools phase.c:807-811
	// implicit: phase() is called at the end and dump_aln drains all).
	if bs != nil {
		var err error
		cursor, err = dumpAln(recs, cursor, int32(math.MaxInt32), hash, bs, rng, opts.DropAmbiguous)
		if err != nil {
			return emitted, err
		}
		_ = cursor
	}
	// Chromosome boundary: upstream deletes every fragment from the shared hash
	// (update_vpos(0x7fffffff), phase.c:732) once the previous chr's trailing
	// block has been flushed and its reads drained, leaving the table at its
	// current n_buckets full of tombstones. The next reference reuses that table.
	// The last reference has no following chr, so upstream never clears it.
	if !isLastRef {
		updateVpos(math.MaxInt32, hash)
	}
	return emitted, nil
}

// hashHasFragSpanning returns true iff any live frag in `hash` has
// vpos >= the given block-vpos. Mirrors the `i` flag computed by
// upstream cleanSeqs (phase.c:410). Since our cleanSeqs runs inside
// phaseEmit and discards its return value, the check is repeated here.
func hashHasFragSpanning(vpos int, hash *fragKhash) bool {
	for k := uint32(0); k < hash.end(); k++ {
		if !hash.exist(k) {
			continue
		}
		f := &hash.vals[k]
		if int(f.vpos) >= vpos {
			return true
		}
	}
	return false
}

// admitPhaseBase decides whether a single pileup base is admitted into
// the het-detection consensus, and if so returns its packed errmod
// observation word. It is the pure, byte-faithful port of the per-base
// body of phase.c's het-detection loop (phase.c:738-751), factored out
// so the admission rule can be unit-tested without an upstream binary.
//
// seqChar is the ASCII nucleotide at the column (A/C/G/T/N/...), baseQ
// the Phred base quality, mapQ the read's MAPQ, rev whether the read is
// on the reverse strand, and minBaseQ the -Q/--min-BQ threshold.
//
// A base is dropped (ok=false) when its quality is below minBaseQ or it
// is not one of A/C/G/T (seq_nt16_int code > 3 upstream — N and other
// IUPAC codes never reach errmod). The packed word is
// q<<5 | rev<<4 | base, where q = clamp(min(baseQ, mapQ), 4, 63) and
// base is the 0..3 code, exactly as phase.c builds bases[k].
func admitPhaseBase(seqChar, baseQ, mapQ byte, rev bool, minBaseQ uint8) (uint16, bool) {
	if baseQ < minBaseQ {
		return 0, false
	}
	b, ok := nt16ToInt(seqChar)
	if !ok || b > 3 {
		return 0, false
	}
	qq := int(baseQ)
	if int(mapQ) < qq {
		qq = int(mapQ)
	}
	if qq < 4 {
		qq = 4
	} else if qq > 63 {
		qq = 63
	}
	r := uint16(0)
	if rev {
		r = 1
	}
	return uint16(qq)<<5 | r<<4 | uint16(b), true
}

// phaseHetLOD extracts the heterozygous Phred-scaled LOD from a gl2cns
// result word. It mirrors phase.c's variant-column test expression
// `(c&0xffff)>>2` (phase.c:758). A column is a variant (het) site iff
// this LOD is >= the -q/min_varLOD threshold (or the site is in the
// forced -l list).
func phaseHetLOD(c uint32) int { return int((c & 0xffff) >> 2) }

// isPhaseVariantColumn reports whether a column with gl2cns result c is
// admitted as a variant (het) site given the -q minVarLOD threshold and
// whether the position is forced by the -l site list (inSet). This is
// the exact upstream gate at phase.c:758: forced sites always pass; all
// others must reach the LOD threshold.
func isPhaseVariantColumn(c uint32, minVarLOD int, inSet bool) bool {
	return inSet || phaseHetLOD(c) >= minVarLOD
}

// nt16ToInt converts a single ACGT/N character to a 0..4 code (A=0,
// C=1, G=2, T=3, N=4). Returns false for non-IUPAC characters.
func nt16ToInt(b byte) (int, bool) {
	switch b {
	case 'A', 'a':
		return 0, true
	case 'C', 'c':
		return 1, true
	case 'G', 'g':
		return 2, true
	case 'T', 't':
		return 3, true
	case 'N', 'n':
		return 4, true
	}
	return 0, false
}

// phaseEmit writes the PS / FL / M / EV lines for one block of vpos
// het positions. Mirrors the body of `phase()` in phase.c. Returns the
// number of het sites consumed (which equals vpos).
func phaseEmit(g *upstreamPhaseRunner, chr string, vpos int, cns []uint64, hash *fragKhash, bw *bufio.Writer) (int, error) {
	if vpos == 0 {
		return 0, nil
	}
	cleanSeqs(vpos, hash)
	if vpos == 1 {
		c0 := int(cns[0] >> 32)
		fmt.Fprintf(bw, "PS\t%s\t%d\t%d\n", chr, c0+1, c0+1)
		fmt.Fprintf(bw, "M0\t%s\t%d\t%d\t%c\t%c\t%d\t0\t0\t0\t0\n//\n",
			chr, c0+1, c0+1,
			"ACGTX"[cns[0]&3], "ACGTX"[cns[0]>>16&3], g.vposShift+1)
		for k := uint32(0); k < hash.end(); k++ {
			if !hash.exist(k) {
				continue
			}
			f := &hash.vals[k]
			if f.vpos != 0 {
				continue
			}
			f.flip = 0
			if f.seq[0] == 0 {
				f.phased = 0
			} else {
				f.phased = 1
				if f.seq[0] == 1 {
					f.phase = 0
				} else {
					f.phase = 1
				}
			}
		}
		g.vposShift++
		return 1, nil
	}

	first := int(cns[0] >> 32)
	last := int(cns[vpos-1] >> 32)
	fmt.Fprintf(bw, "PS\t%s\t%d\t%d\n", chr, first+1, last+1)

	siteMask := make([]int8, vpos)
	cntMat := countAll(g.k, vpos, hash)
	path := dynaprog(g.k, vpos, cntMat)
	pcnt := fragphase(vpos, path, hash, false)
	var nMasked int
	mask := genmask(vpos, pcnt, &nMasked)
	regMask := make([]uint64, nMasked)
	for i := 0; i < nMasked; i++ {
		regMask[i] = (cns[mask[i]>>32]>>32)<<32 | (cns[uint32(mask[i])] >> 32)
		for j := int(mask[i] >> 32); j <= int(int32(mask[i])); j++ {
			siteMask[j] = 1
		}
	}
	if g.fixChimera {
		pcnt = fragphase(vpos, path, hash, true)
	}

	for i := 0; i < nMasked; i++ {
		fmt.Fprintf(bw, "FL\t%s\t%d\t%d\n", chr, int(regMask[i]>>32)+1, int(uint32(regMask[i]))+1)
	}
	for i := 0; i < vpos; i++ {
		x := pcnt[i]
		var c [2]int8
		if (cns[i]&0xffff)>>2 == 0 {
			c[0] = 4
		} else {
			c[0] = int8(cns[i] & 3)
		}
		if (cns[i]>>16&0xffff)>>2 == 0 {
			c[1] = 4
		} else {
			c[1] = int8(cns[i] >> 16 & 3)
		}
		fmt.Fprintf(bw, "M%d\t%s\t%d\t%d\t%c\t%c\t%d\t%d\t%d\t%d\t%d\n",
			siteMask[i]+1, chr, first+1, int(cns[i]>>32)+1,
			"ACGTX"[c[path[i]]], "ACGTX"[c[1-path[i]]],
			i+g.vposShift+1,
			int(x&0xffff), int(x>>16&0xffff),
			int(x>>32&0xffff), int(x>>48&0xffff))
	}

	nSeqs := kSize(hash)
	seqs := make([]fragPtr, 0, nSeqs)
	for k := uint32(0); k < hash.end(); k++ {
		if !hash.exist(k) {
			continue
		}
		f := &hash.vals[k]
		if int(f.vpos) < vpos && f.single == 0 {
			seqs = append(seqs, fragPtr{f: f, bucket: k})
		}
	}
	ksortRseq(seqs)
	for _, sp := range seqs {
		f := sp.f
		fmt.Fprintf(bw, "EV\t0\t%s\t%d\t40\t%dM\t*\t0\t0\t",
			chr, int(f.vpos)+1+g.vposShift, int(f.vlen))
		for j := 0; j < int(f.vlen); j++ {
			c := cns[int(f.vpos)+j]
			if f.seq[j] == 0 {
				bw.WriteByte('N')
			} else if f.seq[j] == 1 {
				bw.WriteByte("ACGT"[c&3])
			} else {
				bw.WriteByte("ACGT"[c>>16&3])
			}
		}
		fmt.Fprintf(bw, "\t*\tYP:i:%d\tYF:i:%d\tYI:i:%d\tYO:i:%d\tYS:i:%d\n",
			f.phase, f.flip, f.in, f.out, f.beg+1)
	}
	bw.WriteString("//\n")
	g.vposShift += vpos
	return vpos, nil
}

// emitPhaseBanner writes the CC header banner upstream prints once at
// the start of phase output.
func emitPhaseBanner(bw *bufio.Writer) {
	bw.WriteString("CC\n")
	bw.WriteString("CC\tDescriptions:\nCC\n")
	bw.WriteString("CC\t  CC      comments\n")
	bw.WriteString("CC\t  PS      start of a phase set\n")
	bw.WriteString("CC\t  FL      filtered region\n")
	bw.WriteString("CC\t  M[012]  markers; 0 for singletons, 1 for phased and 2 for filtered\n")
	bw.WriteString("CC\t  EV      supporting reads; SAM format\n")
	bw.WriteString("CC\t  //      end of a phase set\nCC\n")
	bw.WriteString("CC\tFormats of PS, FL and M[012] lines (1-based coordinates):\nCC\n")
	bw.WriteString("CC\t  PS  chr  phaseSetStart  phaseSetEnd\n")
	bw.WriteString("CC\t  FL  chr  filterStart    filterEnd\n")
	bw.WriteString("CC\t  M?  chr  PS  pos  allele0  allele1  hetIndex  #supports0  #errors0  #supp1  #err1\n")
	bw.WriteString("CC\nCC\n")
}
