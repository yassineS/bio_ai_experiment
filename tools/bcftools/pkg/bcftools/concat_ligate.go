package bcftools

import (
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// This file implements `bcftools concat --ligate` (the phased-concat mode).
//
// The input files are phased chunks that overlap at their boundaries: the
// trailing records of chunk N cover the same sites as the leading records of
// chunk N+1. Naively concatenating them would duplicate the overlap and could
// leave the haplotypes of adjacent chunks phased inconsistently. Ligation
// streams the records in genomic order, votes (per diploid phased sample) over
// the overlap whether chunk N+1's haplotypes agree with chunk N as-is or
// swapped, flips the swapped samples in all subsequent records of the trailing
// chunk, emits each overlap site exactly once, and annotates the output with
// FORMAT/PS (phase set) and FORMAT/PQ (phase quality) fields.
//
// The algorithm mirrors vcfconcat.c's phased_concat / phased_flush /
// phase_update so the output matches upstream byte-for-byte (modulo the
// non-reproducible ##bcftools_concat provenance header line).

// The PS/PQ FORMAT header lines, matched verbatim to those bcftools appends in
// vcfconcat.c:init_data so the merged header is byte-identical to upstream's.
const (
	ligatePQHeader = `##FORMAT=<ID=PQ,Number=1,Type=Integer,Description="Phasing Quality (bigger is better)">`
	ligatePSHeader = `##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set">`
)

// ligateState carries the per-sample phase bookkeeping across chunk boundaries.
type ligateState struct {
	nsmpl     int
	minPQ     int
	swapPhase []int   // per-sample: 1 if the trailing chunk's haplotypes are flipped
	nmatch    []int   // per-sample overlap vote tallies (reset per boundary)
	nmism     []int   //
	phaseQual []int32 // per-sample phasing quality at the most recent boundary
	phaseSet  []int32 // per-sample current phase-set start (1-based POS)
	prevChrom string  // chromosome of the previously emitted records
}

// ligateConcat performs phased concatenation over the per-file variant groups.
// It returns the ligated record stream (with PS/PQ injected and swaps applied)
// in output order. Sample count comes from merged.
func ligateConcat(merged *vcf.Header, groups [][]*vcf.Variant, opts ConcatOptions) ([]*vcf.Variant, error) {
	st := &ligateState{
		nsmpl:     len(merged.Samples),
		minPQ:     opts.MinPQ,
		swapPhase: make([]int, len(merged.Samples)),
		nmatch:    make([]int, len(merged.Samples)),
		nmism:     make([]int, len(merged.Samples)),
		phaseQual: make([]int32, len(merged.Samples)),
		phaseSet:  make([]int32, len(merged.Samples)),
	}

	var out []*vcf.Variant
	if len(groups) == 0 {
		return out, nil
	}

	// Process chunks pairwise. carry holds the as-yet-unwritten tail of the
	// current logical stream: for the first file it is the whole file;
	// thereafter it is whatever survived the previous overlap merge (the
	// trailing non-overlap of the previous "next" chunk). For each adjacent
	// pair we split into (head-of-N outside overlap) + (overlap) +
	// (tail-of-(N+1) outside overlap), emit the head and the resolved overlap,
	// and carry the tail. This mirrors upstream's two-open-files window.
	st.startChromIfNeeded(groups[0])

	carry := groups[0]
	for fi := 1; fi < len(groups); fi++ {
		next := groups[fi]
		emitted, newCarry, err := st.ligatePair(carry, next, opts)
		if err != nil {
			return nil, err
		}
		out = append(out, emitted...)
		carry = newCarry
	}
	// Flush the final carry: no further overlap, write under current swap with
	// the running PS.
	out = append(out, st.flushNonOverlap(carry)...)
	return out, nil
}

// startChromIfNeeded initialises phase_set when the first record introduces a
// new chromosome (mirrors the prev_chr handling in phased_push).
func (st *ligateState) startChromIfNeeded(recs []*vcf.Variant) {
	if len(recs) == 0 {
		return
	}
	first := recs[0]
	if st.prevChrom == "" || st.prevChrom != first.Chrom {
		for i := 0; i < st.nsmpl; i++ {
			st.phaseSet[i] = int32(first.Pos)
		}
		st.prevChrom = first.Chrom
	}
}

// ligatePair ligates the current carry against the next chunk and returns the
// records safe to emit now plus the carry (the trailing non-overlap of next).
func (st *ligateState) ligatePair(cur, next []*vcf.Variant, opts ConcatOptions) (emit, carry []*vcf.Variant, err error) {
	curOverStart, nextOverEnd, ok := matchOverlap(cur, next)
	if !ok {
		// No overlap detected. Upstream errors unless --ligate-force/-warn.
		if !opts.LigateForce && !opts.LigateWarn {
			site := "(empty)"
			if len(next) > 0 {
				site = fmt.Sprintf("%s:%d", next[0].Chrom, next[0].Pos)
			}
			return nil, nil, fmt.Errorf("the --ligate option is intended for VCFs with perfect overlap, the site %s breaks the assumption", site)
		}
		// Force/warn: treat as a clean join with no swap re-evaluation.
		emit = st.flushNonOverlap(cur)
		st.startChromIfNeeded(next)
		return emit, next, nil
	}

	// Records of cur strictly before the overlap are written now under the
	// current swap state with the running PS.
	head := cur[:curOverStart]
	emit = append(emit, st.flushNonOverlap(head)...)

	// Build the overlap pair buffer.
	aOver := cur[curOverStart:] // from chunk N (file0)
	bOver := next[:nextOverEnd] // from chunk N+1 (file1)
	pairs := buildPairs(aOver, bOver)

	// Vote (phased_flush pass 1): tally nmatch/nmism over pairs present in both.
	st.reset()
	for _, p := range pairs {
		if p.a == nil || p.b == nil {
			continue
		}
		st.vote(p.a, p.b)
	}

	// First half of the overlap buffer: written under the PREVIOUS swap state.
	//
	// Upstream stores pairs in a flat slot array of length nbuf = 2*npairs and
	// splits it with `for(i=0; i<nbuf/2; i+=2)` then `for(; i<nbuf; i+=2)`.
	// Stepping i by 2 over [0, nbuf/2) visits ceil(npairs/2) pairs in the first
	// half, leaving floor(npairs/2) for the second. Reproduce that split point
	// exactly so the PQ boundary record matches upstream.
	npairs := len(pairs)
	half := (npairs + 1) / 2
	nswapPrev := st.nswap()
	for i := 0; i < half; i++ {
		p := pairs[i]
		rec := p.a
		fromA := true
		if rec == nil {
			rec = p.b
			fromA = false
		}
		if nswapPrev > 0 && fromA {
			st.phaseUpdate(rec)
		}
		st.applyPS(rec)
		emit = append(emit, rec)
	}

	// Decide the new swap state and phase quality (phased_flush mid section).
	st.decideSwap()

	// Second half: written under the NEW swap state; the first mask==3 record
	// carries PQ and triggers low-PQ phase-set breaks.
	pqPrinted := false
	nswapNew := st.nswap()
	for i := half; i < npairs; i++ {
		p := pairs[i]
		rec := p.b
		if rec == nil {
			rec = p.a
		}
		both := p.a != nil && p.b != nil
		if !pqPrinted && both {
			st.applyPQ(rec)
			pqPrinted = true
			for j := 0; j < st.nsmpl; j++ {
				if st.phaseQual[j] < int32(st.minPQ) {
					st.phaseSet[j] = int32(rec.Pos)
				}
			}
		}
		if nswapNew > 0 {
			st.phaseUpdate(rec)
		}
		st.applyPS(rec)
		emit = append(emit, rec)
	}

	carry = next[nextOverEnd:]
	st.startChromIfNeeded(carry)
	return emit, carry, nil
}

// flushNonOverlap emits records outside any overlap under the current swap
// state, applying PS to each (compact_PS is off by default so every record
// carries PS).
func (st *ligateState) flushNonOverlap(recs []*vcf.Variant) []*vcf.Variant {
	out := make([]*vcf.Variant, 0, len(recs))
	for _, r := range recs {
		// A record may begin a new chromosome; reset PS accordingly.
		if st.prevChrom != "" && r.Chrom != st.prevChrom {
			for i := 0; i < st.nsmpl; i++ {
				st.phaseSet[i] = int32(r.Pos)
			}
			st.prevChrom = r.Chrom
		}
		if st.nswap() > 0 {
			st.phaseUpdate(r)
		}
		st.applyPS(r)
		out = append(out, r)
	}
	return out
}

// reset zeroes the per-boundary vote tallies.
func (st *ligateState) reset() {
	for i := range st.nmatch {
		st.nmatch[i] = 0
		st.nmism[i] = 0
	}
}

// nswap reports how many samples are currently swapped.
func (st *ligateState) nswap() int {
	n := 0
	for _, s := range st.swapPhase {
		if s != 0 {
			n++
		}
	}
	return n
}

// vote tallies nmatch/nmism for one overlapping pair of records (phased_flush
// inner loop). Only diploid, phased, heterozygous, non-missing genotypes vote.
func (st *ligateState) vote(a, b *vcf.Variant) {
	for j := 0; j < st.nsmpl; j++ {
		ga, oka := diploidGT(a, j)
		gb, okb := diploidGT(b, j)
		if !oka || !okb {
			continue
		}
		if !ga.phased || !gb.phased {
			continue
		}
		if ga.missing() || gb.missing() {
			continue
		}
		if ga.a0 == ga.a1 || gb.a0 == gb.a1 {
			continue // homozygous on either side
		}
		if ga.a0 == gb.a0 && ga.a1 == gb.a1 {
			if st.swapPhase[j] != 0 {
				st.nmism[j]++
			} else {
				st.nmatch[j]++
			}
		}
		if ga.a0 == gb.a1 && ga.a1 == gb.a0 {
			if st.swapPhase[j] != 0 {
				st.nmatch[j]++
			} else {
				st.nmism[j]++
			}
		}
	}
}

// decideSwap finalises swap_phase and phase_qual for the boundary using the
// accumulated votes (phased_flush mid section).
func (st *ligateState) decideSwap() {
	for j := 0; j < st.nsmpl; j++ {
		if st.nmatch[j] >= st.nmism[j] {
			st.swapPhase[j] = 0
		} else {
			st.swapPhase[j] = 1
		}
		if st.nmatch[j] != 0 && st.nmism[j] != 0 {
			f := float64(st.nmatch[j]) / float64(st.nmatch[j]+st.nmism[j])
			st.phaseQual[j] = int32(99 * (0.7 + f*math.Log(f) + (1-f)*math.Log(1-f)) / 0.7)
		} else {
			st.phaseQual[j] = 99
		}
		st.nmatch[j] = 0
		st.nmism[j] = 0
	}
}

// phaseUpdate flips the haplotypes of swapped samples in rec, keeping them
// phased (phase_update in vcfconcat.c).
func (st *ligateState) phaseUpdate(rec *vcf.Variant) {
	if formatIndex(rec.Format, "GT") < 0 {
		return
	}
	for j := 0; j < st.nsmpl && j < len(rec.Samples); j++ {
		if st.swapPhase[j] == 0 {
			continue
		}
		g, ok := diploidGT(rec, j)
		if !ok {
			continue
		}
		if g.missing() || !g.phased {
			continue
		}
		// Flip and keep phased.
		rec.Samples[j].Data["GT"] = fmt.Sprintf("%s|%s", alleleStr(g.a1, g.a1miss), alleleStr(g.a0, g.a0miss))
	}
}

// applyPS writes the current per-sample phase set into rec's FORMAT/PS.
func (st *ligateState) applyPS(rec *vcf.Variant) {
	setFormatInt32(rec, "PS", st.phaseSet, st.nsmpl)
}

// applyPQ writes the current per-sample phase quality into rec's FORMAT/PQ.
func (st *ligateState) applyPQ(rec *vcf.Variant) {
	setFormatInt32(rec, "PQ", st.phaseQual, st.nsmpl)
}

// pair is one overlap site: a is the chunk-N record, b the chunk-(N+1) record;
// either may be nil if the site is present in only one chunk.
type pair struct {
	a, b *vcf.Variant
}

// buildPairs aligns the tail of chunk N (aOver) with the head of chunk N+1
// (bOver) by (CHROM, POS, REF, ALT), producing the merged overlap buffer in
// genomic order. Sites present in only one chunk get a nil counterpart.
func buildPairs(aOver, bOver []*vcf.Variant) []pair {
	var pairs []pair
	i, j := 0, 0
	for i < len(aOver) || j < len(bOver) {
		switch {
		case i >= len(aOver):
			pairs = append(pairs, pair{a: nil, b: bOver[j]})
			j++
		case j >= len(bOver):
			pairs = append(pairs, pair{a: aOver[i], b: nil})
			i++
		default:
			a, b := aOver[i], bOver[j]
			switch {
			case sameSite(a, b):
				pairs = append(pairs, pair{a: a, b: b})
				i++
				j++
			case siteLess(a, b):
				pairs = append(pairs, pair{a: a, b: nil})
				i++
			default:
				pairs = append(pairs, pair{a: nil, b: b})
				j++
			}
		}
	}
	return pairs
}

// matchOverlap finds the overlap between the tail of cur and the head of next.
// It returns the index in cur where the overlap begins and the index in next
// where the overlap ends, plus ok=false when no overlap exists. The overlap is
// the cur-suffix whose first site is on next[0]'s chromosome with POS >=
// next[0].Pos, matched against next's leading sites up to cur's last position.
func matchOverlap(cur, next []*vcf.Variant) (curStart, nextEnd int, ok bool) {
	if len(cur) == 0 || len(next) == 0 {
		return 0, 0, false
	}
	nb := next[0]
	curStart = len(cur)
	for i := 0; i < len(cur); i++ {
		if cur[i].Chrom == nb.Chrom && cur[i].Pos >= nb.Pos {
			curStart = i
			break
		}
	}
	if curStart >= len(cur) {
		return 0, 0, false // cur ends before next starts: no overlap
	}
	lastPos := cur[len(cur)-1].Pos
	lastChrom := cur[len(cur)-1].Chrom
	for j := 0; j < len(next); j++ {
		if next[j].Chrom == lastChrom && next[j].Pos <= lastPos {
			nextEnd = j + 1
		} else {
			break
		}
	}
	if nextEnd == 0 {
		return 0, 0, false
	}
	return curStart, nextEnd, true
}

// sameSite reports whether two records share CHROM, POS, REF and ALT.
func sameSite(a, b *vcf.Variant) bool {
	if a.Chrom != b.Chrom || a.Pos != b.Pos || a.Ref != b.Ref {
		return false
	}
	if len(a.Alt) != len(b.Alt) {
		return false
	}
	for i := range a.Alt {
		if a.Alt[i] != b.Alt[i] {
			return false
		}
	}
	return true
}

// siteLess orders two records by (POS, REF, ALT) within a shared chromosome.
func siteLess(a, b *vcf.Variant) bool {
	if a.Pos != b.Pos {
		return a.Pos < b.Pos
	}
	if a.Ref != b.Ref {
		return a.Ref < b.Ref
	}
	return strings.Join(a.Alt, ",") < strings.Join(b.Alt, ",")
}

// gt holds a parsed diploid genotype.
type gt struct {
	a0, a1         int
	a0miss, a1miss bool
	phased         bool
}

func (g gt) missing() bool { return g.a0miss || g.a1miss }

// diploidGT parses sample j's GT field as a diploid genotype. ok is false when
// the sample has no GT, is not diploid, or cannot be parsed.
func diploidGT(v *vcf.Variant, j int) (gt, bool) {
	if j >= len(v.Samples) {
		return gt{}, false
	}
	s, ok := v.Samples[j].Data["GT"]
	if !ok {
		return gt{}, false
	}
	phased := strings.IndexByte(s, '|') >= 0
	var sep byte = '/'
	if phased {
		sep = '|'
	}
	idx := strings.IndexByte(s, sep)
	if idx < 0 {
		return gt{}, false // haploid or malformed
	}
	left, right := s[:idx], s[idx+1:]
	// A diploid GT must not contain a further separator.
	if strings.IndexAny(right, "|/") >= 0 {
		return gt{}, false
	}
	a0, m0 := parseAllele(left)
	a1, m1 := parseAllele(right)
	return gt{a0: a0, a1: a1, a0miss: m0, a1miss: m1, phased: phased}, true
}

// alleleStr renders an allele index back to its string form.
func alleleStr(a int, miss bool) string {
	if miss {
		return "."
	}
	return strconv.Itoa(a)
}

// setFormatInt32 writes an Integer FORMAT field for every sample, appending the
// tag to FORMAT if absent. A missing per-sample value (INT32 min sentinel) is
// rendered as ".". Mirrors bcf_update_format_int32 semantics for our text
// model.
func setFormatInt32(v *vcf.Variant, tag string, vals []int32, nsmpl int) {
	if formatIndex(v.Format, tag) < 0 {
		v.Format = append(v.Format, tag)
	}
	for j := 0; j < nsmpl && j < len(v.Samples); j++ {
		if v.Samples[j].Data == nil {
			v.Samples[j].Data = map[string]string{}
		}
		var val int32
		if j < len(vals) {
			val = vals[j]
		}
		if val == int32(math.MinInt32) {
			v.Samples[j].Data[tag] = "."
		} else {
			v.Samples[j].Data[tag] = strconv.FormatInt(int64(val), 10)
		}
	}
}

// ensureLigateHeaders adds the PS and PQ FORMAT header lines to merged if they
// are not already present, matching the order and text upstream emits.
func ensureLigateHeaders(merged *vcf.Header) {
	hasPQ, hasPS := false, false
	for _, m := range merged.MetaInfo {
		if k, id := structuredID(m); k == "FORMAT" {
			switch id {
			case "PQ":
				hasPQ = true
			case "PS":
				hasPS = true
			}
		}
	}
	if !hasPQ {
		merged.MetaInfo = append(merged.MetaInfo, ligatePQHeader)
	}
	if !hasPS {
		merged.MetaInfo = append(merged.MetaInfo, ligatePSHeader)
	}
}
