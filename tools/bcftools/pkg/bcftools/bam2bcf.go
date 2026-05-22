// Per-site genotype-likelihood pipeline for `bcftools mpileup`, ported
// from bcftools' bam2bcf.c. This is slice 2 of the MAQ-model work: it
// wires the errmod error model (slice 1, errmod.go) into the three
// upstream functions that turn a pileup column into a BCF/VCF record:
//
//   - bcfCallGlfgen   (bam2bcf.c:250 bcf_call_glfgen) — one sample's
//     pileup column → packed base array → errmod → per-call QS / I16 /
//     genotype-likelihood annotations.
//   - bcfCallCombine  (bam2bcf.c:959 bcf_call_combine) — combine the
//     per-sample calls, order alleles by coverage-normalised QS, append
//     the `<*>` unseen allele, and build the multi-allelic PL grid.
//   - bcfCall2bcf     (bam2bcf.c:1200 bcf_call2bcf) — emit the record:
//     REF, ALT (incl. `<*>`), QUAL, INFO/DP/I16/QS, FORMAT/PL.
//
// SNP-only: indel branches of the upstream functions are skipped. BAQ
// realignment (slice 3) and the bias annotations VDB/SGB/RPBZ/... (slice
// 4) are deliberately out of scope; records are emitted without those
// INFO tags.

package bcftools

import (
	"math"
)

// b2bMaxAlleles is upstream's B2B_MAX_ALLELES: REF + up to 3 ALT + the
// `<*>` unseen allele, i.e. a 5x5 genotype matrix.
const b2bMaxAlleles = 5

// b2bDefMapQ mirrors bam2bcf.c's DEF_MAPQ: the mapping quality assigned
// to reads whose MAPQ is the "unavailable" sentinel 255.
const b2bDefMapQ = 20

// b2bCapDist mirrors bam2bcf.c's CAP_DIST: the cap on the tail-distance
// annotation (distance of a base from the nearer read end).
const b2bCapDist = 25

// b2bCapQ mirrors the bca->capQ value set unconditionally to 60 by
// bcf_call_init (bam2bcf.c:50). For the SNP path it is never raised; it
// is the upper bound on mapping quality applied in the glfgen cap
// sequence (bam2bcf.c:459-461).
const b2bCapQ = 60

// b2bSeqQ is the seqQ ceiling on q for the SNP path. Indels recompute
// seqQ from the indel model; for SNPs it stays at 99 (bam2bcf.c:459).
const b2bSeqQ = 99

// seqNt16Int maps a 4-bit IUPAC base code to the 2-bit base index used
// by the genotype matrix (0=A,1=C,2=G,3=T, 4=other/N). It mirrors
// htslib's seq_nt16_int table.
var seqNt16Int = [16]int{4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4}

// baseToNt16 maps an uppercase ASCII nucleotide to its 4-bit IUPAC code.
// Anything not A/C/G/T maps to 15 (N), which seqNt16Int turns into 4.
func baseToNt16(b byte) int {
	switch b {
	case 'A':
		return 1
	case 'C':
		return 2
	case 'G':
		return 4
	case 'T':
		return 8
	}
	return 15
}

// pileupBase is one read's contribution to one reference position, with
// every value the glfgen stage needs already resolved: the 2-bit base
// index, mapping quality, strand, and the raw query quality plus query
// position so glfgen can apply the delta_baseQ neighbour-quality cap.
type pileupBase struct {
	base4   int    // 2-bit base index 0..4
	rawQual uint8  // base quality at this query position (pre-cap)
	prevQ   int    // quality of the neighbour at qpos-1, or -1 if none
	nextQ   int    // quality of the neighbour at qpos+1, or -1 if none
	mapq    uint8  // read mapping quality (255 sentinel handled in glfgen)
	reverse bool   // read is on the reverse strand
	qpos    int    // 0-based query position of this base
	qlen    int    // read query length (for tail-distance annotation)
	qname   string // read name, for overlap detection
}

// bcfCallret is the per-sample result of bcfCallGlfgen, mirroring
// bcftools' bcf_callret1_t for the SNP path. The indel-only fields are
// omitted.
type bcfCallret struct {
	// p holds the 25-entry (5x5) phred-scaled genotype-likelihood matrix
	// produced by errmod; p[i*5+j] is genotype (i,j).
	p [25]float32
	// qs is the per-allele quality sum (FORMAT/INFO QS), indexed by
	// 2-bit base.
	qs [4]float64
	// anno is the 16-slot I16 accumulator. Slots are grouped in fours:
	// [0..3] strand-split depth (ref-fwd, ref-rev, alt-fwd, alt-rev),
	// [4..7] baseQ sum/sum-of-squares (ref, alt), [8..11] mapQ likewise,
	// [12..15] tail-distance likewise — matching bam2bcf.h's compute_I16.
	anno [16]float64
	// oriDepth is the read depth before --min-BQ filtering.
	oriDepth int
	// mq0 counts reads with mapping quality zero.
	mq0 int
	// adf / adr are per-allele forward/reverse depths (FORMAT AD support).
	adf [b2bMaxAlleles]int
	adr [b2bMaxAlleles]int
}

// bcfCallGlfgen is the Go port of bcf_call_glfgen (bam2bcf.c:250) for
// the SNP path. It collects one sample's pileup column into the packed
// base array, applies the delta_baseQ neighbour-quality cap, runs the
// errmod MAQ model, and accumulates the QS / I16 / AD annotations into
// r. ref4 is the 2-bit reference base index (0..4). It returns the
// number of bases that passed the quality filter.
func bcfCallGlfgen(pile []pileupBase, ref4 int, opts MpileupOptions, em *Errmod, r *bcfCallret) int {
	*r = bcfCallret{}
	if len(pile) == 0 {
		return 0
	}
	minBQ := int(opts.MinBQ)
	maxBQ := int(opts.MaxBQ)
	deltaBQ := opts.DeltaBQ

	bases := make([]uint16, 0, len(pile))
	for i := range pile {
		p := &pile[i]
		r.oriDepth++

		b := p.base4
		// Lowest of this base and its neighbours' qualities, capped via
		// delta_baseQ (bam2bcf.c:~427-435).
		q := int(p.rawQual)
		if p.prevQ >= 0 && q > p.prevQ+deltaBQ {
			q = p.prevQ + deltaBQ
		}
		if p.nextQ >= 0 && q > p.nextQ+deltaBQ {
			q = p.nextQ + deltaBQ
		}
		if q < minBQ {
			continue
		}
		if q > maxBQ {
			q = maxBQ
		}
		baseQ := q

		mapQ := int(p.mapq)
		if mapQ == 255 {
			mapQ = b2bDefMapQ
		}
		if mapQ == 0 {
			r.mq0++
		}
		// Cap sequence, literal port of bam2bcf.c:459-461:
		//   if (q > seqQ) q = seqQ;                       // seqQ = 99
		//   mapQ = mapQ < bca->capQ ? mapQ : bca->capQ;   // capQ = 60
		//   if (q > mapQ) q = mapQ;
		// i.e. mapQ is itself capped to 60 first, then q is capped to the
		// capped mapQ. Then q is floored into [4,63] (bam2bcf.c:462-463).
		if q > b2bSeqQ {
			q = b2bSeqQ
		}
		if mapQ > b2bCapQ {
			mapQ = b2bCapQ
		}
		if q > mapQ {
			q = mapQ
		}
		if q > 63 {
			q = 63
		}
		if q < 4 {
			q = 4
		}
		rev := 0
		if p.reverse {
			rev = 1
		}
		bases = append(bases, uint16(q<<5|rev<<4|b))

		isDiff := 0
		if !(ref4 < 4 && b == ref4) {
			isDiff = 1
		}
		if b < 4 {
			r.qs[b] += float64(q)
			if p.reverse {
				r.adr[b]++
			} else {
				r.adf[b]++
			}
		}
		// I16 annotations (compute_I16 layout).
		r.anno[0<<2|isDiff<<1|rev]++
		minDist := p.qlen - 1 - p.qpos
		if minDist > p.qpos {
			minDist = p.qpos
		}
		if minDist > b2bCapDist {
			minDist = b2bCapDist
		}
		bq := float64(baseQ)
		mq := float64(mapQ)
		md := float64(minDist)
		r.anno[1<<2|isDiff<<1|0] += bq
		r.anno[1<<2|isDiff<<1|1] += bq * bq
		r.anno[2<<2|isDiff<<1|0] += mq
		r.anno[2<<2|isDiff<<1|1] += mq * mq
		r.anno[3<<2|isDiff<<1|0] += md
		r.anno[3<<2|isDiff<<1|1] += md * md
	}

	// glfgen: errmod turns the packed base array into the 5x5 PL matrix.
	pl := r.p[:]
	em.ErrmodCal(len(bases), b2bMaxAlleles, bases, pl, nil)
	return len(bases)
}

// bcfCall is the combined multi-sample call, mirroring bcftools'
// bcf_call_t for the SNP path.
type bcfCall struct {
	// alleles holds the 2-bit base index of each allele in output order:
	// REF first, then ALTs by descending QS, then the unseen allele.
	// Entries are -1 for unused slots.
	alleles [b2bMaxAlleles]int
	// nAlleles is the number of populated allele slots.
	nAlleles int
	// unseen is the allele index that maps to the `<*>` symbolic allele,
	// or -1 if no unseen allele was appended.
	unseen int
	// qsum is the coverage-normalised QS per output allele.
	qsum [b2bMaxAlleles]float64
	// pl holds, per sample, the upper-triangle genotype-likelihood array
	// of length nAlleles*(nAlleles+1)/2.
	pl [][]int
	// qs holds, per sample, the reordered per-allele QS values.
	qs [][]int
	// adf / adr hold, per sample, the reordered per-allele fwd/rev depths.
	adf [][]int
	adr [][]int
	// anno is the combined 16-slot I16 array.
	anno [16]float64
	// depth is the post-filter total read depth.
	depth int
	// oriDepth is the pre-filter total read depth (INFO/DP).
	oriDepth int
	// mq0 is the total count of MQ0 reads.
	mq0 int
	// oriRef is the 2-bit reference base index (0..4), -1 for indels.
	oriRef int
}

// bcfCallCombine is the Go port of bcf_call_combine (bam2bcf.c:959) for
// the SNP path. It combines the per-sample calls in calls, orders the
// alleles by coverage-normalised QS sums, always appends the `<*>`
// unseen allele, builds the multi-allelic PL grid, and fills a bcfCall.
// ref4 is the 2-bit reference base index.
func bcfCallCombine(calls []bcfCallret, ref4 int) bcfCall {
	var call bcfCall
	call.oriRef = ref4

	// Coverage-normalised QS sums per base (bam2bcf.c:970-977).
	var qsum [b2bMaxAlleles]float64
	for i := range calls {
		var sum float64
		for j := 0; j < 4; j++ {
			sum += calls[i].qs[j]
		}
		if sum != 0 {
			for j := 0; j < 4; j++ {
				qsum[j] += calls[i].qs[j] / sum
			}
		}
	}

	// Sort the 5 qsum slots in ascending order, tracking original index
	// (insertion sort, matching bam2bcf.c:980-984).
	idx := [5]int{0, 1, 2, 3, 4}
	for i := 1; i < 4; i++ {
		for j := i; j > 0 && qsum[idx[j]] < qsum[idx[j-1]]; j-- {
			idx[j], idx[j-1] = idx[j-1], idx[j]
		}
	}

	for i := range call.alleles {
		call.alleles[i] = -1
	}
	call.unseen = -1
	call.alleles[0] = ref4
	// Walk the sorted qsum array from highest to lowest, assigning
	// output allele slots (bam2bcf.c:991-1001).
	j := 1
	i := 3
	for ; i >= 0; i-- {
		ipos := idx[i]
		if ipos == ref4 {
			call.qsum[0] = qsum[ipos]
		} else {
			if qsum[ipos] == 0 {
				break
			}
			call.qsum[j] = qsum[ipos]
			call.alleles[j] = ipos
			j++
		}
	}
	// Always append the `<*>` unseen allele for SNPs (bam2bcf.c:1005-1009).
	if ((ref4 < 4 && j < 4) || (ref4 == 4 && j < 5)) && i >= 0 {
		call.unseen = j
		call.alleles[j] = idx[i]
		j++
	}
	call.nAlleles = j

	// Build the upper-triangle genotype-index list g[]: for output
	// alleles a, g[z++] = a[j]*5 + a[i] for j<=i (bam2bcf.c:1025-1031).
	x := call.nAlleles * (call.nAlleles + 1) / 2
	g := make([]int, 0, x)
	for ii := 0; ii < call.nAlleles; ii++ {
		for jj := 0; jj <= ii; jj++ {
			g = append(g, call.alleles[jj]*5+call.alleles[ii])
		}
	}

	// Per-sample PL: pick out the genotypes in g[], rebase to min,
	// round, cap at 255 (bam2bcf.c:1034-1046).
	call.pl = make([][]int, len(calls))
	for s := range calls {
		r := &calls[s]
		minv := float32(math.MaxFloat32)
		for _, gi := range g {
			if r.p[gi] < minv {
				minv = r.p[gi]
			}
		}
		pl := make([]int, x)
		for k, gi := range g {
			y := int(float64(r.p[gi]-minv) + 0.499)
			if y > 255 {
				y = 255
			}
			pl[k] = y
		}
		call.pl[s] = pl
	}

	// Reorder QS and AD to match the site's allele ordering
	// (bam2bcf.c:1097-1121).
	call.qs = make([][]int, len(calls))
	call.adf = make([][]int, len(calls))
	call.adr = make([][]int, len(calls))
	for s := range calls {
		r := &calls[s]
		qs := make([]int, call.nAlleles)
		adf := make([]int, call.nAlleles)
		adr := make([]int, call.nAlleles)
		for k := 0; k < call.nAlleles; k++ {
			a := call.alleles[k]
			if a >= 0 && a < 4 {
				qs[k] = int(r.qs[a])
				adf[k] = r.adf[a]
				adr[k] = r.adr[a]
			}
		}
		call.qs[s] = qs
		call.adf[s] = adf
		call.adr[s] = adr
	}

	// Combine I16 / depth annotations (bam2bcf.c:1138-1148).
	for i := range calls {
		r := &calls[i]
		call.depth += int(r.anno[0] + r.anno[1] + r.anno[2] + r.anno[3])
		call.oriDepth += r.oriDepth
		call.mq0 += r.mq0
		for j := 0; j < 16; j++ {
			call.anno[j] += r.anno[j]
		}
	}
	return call
}
