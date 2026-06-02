package samtools

import (
	"bufio"
	"fmt"

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
// Returns the number of het sites emitted on this reference.
func runUpstreamPhase(g *upstreamPhaseRunner, recs []*sam.Record, rname string, bw *bufio.Writer) (int, error) {
	pp := newPhaseStreamPileup(recs, rname)
	em := errmod.Init(1.0 - 0.83)
	bases := make([]uint16, 0, g.maxDepth)
	cns := make([]uint64, 0, 256)
	vpos := 0
	hash := newFragKhash()
	g.vposShift = 0
	emitted := 0
	q := make([]float32, 16)

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
			if int(p.qpos) >= len(p.b.Qual) {
				continue
			}
			baseQ := p.b.Qual[p.qpos]
			if baseQ < g.minBaseQ {
				continue
			}
			if int(p.qpos) >= len(p.b.Seq) {
				continue
			}
			b, ok := nt16ToInt(p.b.Seq[p.qpos])
			if !ok || b > 3 {
				continue
			}
			qq := int(baseQ)
			if int(p.b.MapQ) < qq {
				qq = int(p.b.MapQ)
			}
			if qq < 4 {
				qq = 4
			} else if qq > 63 {
				qq = 63
			}
			rev := uint16(0)
			if p.b.Flag&sam.FlagReverse != 0 {
				rev = 1
			}
			bases = append(bases, uint16(qq)<<5|rev<<4|uint16(b))
		}
		if len(bases) == 0 {
			continue
		}
		for i := range q {
			q[i] = 0
		}
		em.Cal(bases, 4, q)
		c := gl2cns(q)
		// Not a variant?
		if int((c&0xffff)>>2) < g.minVarLOD {
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
					f.end = p.b.EndPosition()
				}
				dophase = false
			} else {
				f.beg = p.b.Pos - 1
				f.end = p.b.EndPosition()
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
			updateVpos(vpos, hash)
			cns[0] = cns[vpos]
			vpos = 0
		}
		vpos++
	}
	// Final flush for any remaining vpos > 0.
	if vpos > 0 {
		n2, err := phaseEmit(g, rname, vpos, cns, hash, bw)
		if err != nil {
			return emitted, err
		}
		emitted += n2
	}
	return emitted, nil
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
