// Per-site genotype-likelihood pipeline for `bcftools mpileup`, ported
// from bcftools' bam2bcf.c. This is slice 2 of the MAQ-model work: it
// wires the errmod error model (slice 1, pkg/htsgo/errmod) into the
// three upstream functions that turn a pileup column into a BCF/VCF
// record:
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
// SNP-only: indel branches of the upstream functions are skipped. Slice
// 4 adds the per-site bias annotations (VDB/SGB/RPBZ/MQBZ/BQBZ/MQSBZ/
// SCBZ) and MQ0F: see calcVDB, calcSegBias and calcMWUBiasZ below, and
// the bias tallies threaded through bcfCallret.
//
// BAQ realignment (slice 3) is wired separately in mpileup.go.

package bcftools

import (
	"math"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/errmod"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// b2bNpos is upstream's bca->npos: the number of read-position and
// soft-clip-length bins. get_position rescales every read to this many
// bins so 100 bp reads (the calc_vdb training length) land 1:1.
const b2bNpos = 100

// b2bNqual is upstream's bca->nqual: the number of mapping-quality and
// base-quality bins for the MWU bias tests. baseQ/mapQ are clamped to
// 59 in glfgen so all 60 bins are reachable.
const b2bNqual = 60

// b2bMaxAlleles is upstream's B2B_MAX_ALLELES: REF + up to 3 ALT + the
// `<*>` unseen allele, i.e. a 5x5 genotype matrix.
const b2bMaxAlleles = 5

// b2bNNm mirrors B2B_N_NM (bam2bcf.h:78): the number of NM-bias bins, i.e.
// the maximum number of mismatches counted before clamping. get_aux_nm
// caps inm to [0, b2bNNm-1].
const b2bNNm = 32

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
	// epos, scLen and scDist are the soft-clip-aware position annotations
	// from get_position (bam2bcf.c:146), already rescaled the way glfgen
	// expects: epos is the read-position bin in [0,b2bNpos), scLen is the
	// soft-clip-bias bin in [0,b2bNpos).
	epos  int
	scLen int
	// indel mirrors htslib's bam_pileup1_t.indel field (sam.c pileup
	// engine). It is the length of the indel that begins at the CIGAR
	// position immediately after this column: positive for an
	// insertion of that many query bases, negative for a deletion of
	// that many reference bases, 0 if the next op is a match. The
	// indel-calling helpers (bam2bcf_indel.go) consume this; the SNP
	// path ignores it.
	indel int
	// isDel mirrors htslib's bam_pileup1_t.is_del: the column lands
	// inside this read's CIGAR D op (a spanning deletion). Upstream's
	// SNP branch skips such reads (bam2bcf.c:307 `if (p->is_del &&
	// !is_indel) continue`), but the indel branch lets them through
	// with the chosen-type encoded in p.aux. The base4 slot is unused
	// when isDel=true.
	isDel bool
	// isRefskip mirrors htslib's bam_pileup1_t.is_refskip: the column
	// lands inside this read's CIGAR N op (CREF_SKIP, e.g. an RNA-seq
	// intron). Both SNP and indel branches skip these reads
	// unconditionally (bam2bcf.c:301).
	isRefskip bool
	// aux is per-read scratch space the indel core uses to thread
	// scoring state through a column. Upstream bcftools stuffs the
	// alignment-score / chosen-type bits into a 32-bit field on
	// bcf_callaux_t.bases — here we keep that state per-pileupBase.
	aux uint32
	// rec is a back-pointer to the originating BAM record. It is now
	// populated on every pileupBase: the indel-branch glfgen's
	// alt_pos/ref_pos tallies and bcfCallGapPrep's iref_*/ialt_*
	// histograms both need the read's CIGAR (for soft-clip-aware
	// get_pos) for ALL reads in the column, not just the indel-bearing
	// ones. The SNP path never reads it.
	rec *sam.Record
	// hasSoftClip mirrors htslib's PLP_HAS_SOFT_CLIP flag: true iff
	// the originating read has at least one CIGAR S op anywhere in its
	// alignment. Upstream sets the flag once in pileup_constructor
	// (mpileup.c:317-323) and consumes it in bcf_call_glfgen when
	// counting per-sample SCR (bam2bcf.c:300). We stamp the same bit
	// on every pileupBase produced from such a read so the SCR
	// accumulator can read it from the column directly.
	hasSoftClip bool
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
	// Bias-test tallies, mirroring the per-sample contribution to the
	// shared bca arrays in bcf_call_glfgen (bam2bcf.c:506-537). Upstream
	// keeps these in bcf_callaux_t and they accumulate across samples;
	// here each sample's glfgen fills its own copy and bcfCallCombine
	// sums them. refPos/altPos and refScl/altScl have b2bNpos bins;
	// refMq/altMq, refBq/altBq, fwdMqs/revMqs have b2bNqual bins.
	refPos [b2bNpos]int
	altPos [b2bNpos]int
	refScl [b2bNpos]int
	altScl [b2bNpos]int
	refMq  [b2bNqual]int
	altMq  [b2bNqual]int
	refBq  [b2bNqual]int
	altBq  [b2bNqual]int
	fwdMqs [b2bNqual]int
	revMqs [b2bNqual]int
	// scr is the per-sample count of reads at this column whose CIGAR
	// contains at least one soft-clip op. Filled by glfgen when either
	// B2BInfoSCR or B2BFmtSCR is set (matching bam2bcf.c:300). Combined
	// into bcfCall.scrTotal and bcfCall.scr by bcfCallCombine.
	scr int
	// refNm / altNm are per-sample NM-bias histograms (one bin per
	// possible mismatch count up to b2bNNm-1, matching upstream's
	// bca->ref_nm/alt_nm and r->ref_nm/alt_nm which share storage
	// via bcf_call_t). The driver's SNP and indel glfgen both fill
	// them when any of the NMBZ/NM bits are set; the SNP and indel
	// combine each sum across samples (and across the two glfgen
	// passes, mirroring upstream's shared bca-level accumulator).
	refNm [b2bNNm]int
	altNm [b2bNNm]int
}

// bcfCallGlfgen is the Go port of bcf_call_glfgen (bam2bcf.c:250) for
// the SNP path. It collects one sample's pileup column into the packed
// base array, applies the delta_baseQ neighbour-quality cap, runs the
// errmod MAQ model, and accumulates the QS / I16 / AD annotations into
// r. ref4 is the 2-bit reference base index (0..4). It returns the
// number of bases that passed the quality filter.
func bcfCallGlfgen(pile []pileupBase, ref4 int, opts MpileupOptions, em *errmod.Errmod, r *bcfCallret) int {
	return bcfCallGlfgenCore(pile, ref4, opts, em, r, false)
}

// bcfCallGlfgenIndel is the Go port of the `ref_base < 0` branch of
// bcf_call_glfgen (bam2bcf.c:~300-460). It consumes the per-read p.aux
// words produced by bcfCgpComputeIndelQ (chosen-type index in bits 16-21,
// seqQ in bits 8-15, indelQ in bits 0-7) and accumulates QS / I16 / AD /
// bias tallies into r exactly as glfgen does for SNPs, but indexed by
// indel-type slot instead of 2-bit base. Reads with no indel assignment
// (aux>>16 == 4 — the sentinel from compute_indelQ) are still ranked but
// their b<4 conditional gates the per-allele QS and AD updates.
//
// Returns the number of reads whose quality passed the min-BQ filter.
func bcfCallGlfgenIndel(pile []pileupBase, opts MpileupOptions, em *errmod.Errmod, r *bcfCallret) int {
	return bcfCallGlfgenCore(pile, -1, opts, em, r, true)
}

// bcfCallGlfgenCore is the shared SNP/indel implementation. isIndel=false
// reads b from p.base4 and q from the raw quality+delta cap; isIndel=true
// reads b from p.aux>>16 and q from p.aux&0xff (the precomputed indelQ).
func bcfCallGlfgenCore(pile []pileupBase, ref4 int, opts MpileupOptions, em *errmod.Errmod, r *bcfCallret, isIndel bool) int {
	*r = bcfCallret{}
	if len(pile) == 0 {
		return 0
	}
	minBQ := int(opts.MinBQ)
	maxBQ := int(opts.MaxBQ)
	deltaBQ := opts.DeltaBQ
	// In the indel branch ref4 is forced to 4 to mirror upstream
	// (bam2bcf.c:269 `ref4 = 4`).
	if isIndel {
		ref4 = 4
	}

	bases := make([]uint16, 0, len(pile))
	// ADR_ref_missed / ADF_ref_missed accumulate the indel-branch reads
	// whose precomputed indelQ falls below min-baseQ but which look like
	// REF (b<4, no indel in CIGAR). Upstream uses them to compensate the
	// per-allele AD when --ambig-reads is incAD / incAD0 (bam2bcf.c:540-561).
	var adrRefMissed, adfRefMissed [4]int
	// scrWanted mirrors upstream's gating on (B2B_INFO_SCR|B2B_FMT_SCR)
	// at bam2bcf.c:300. We only tally SCR in the SNP branch because the
	// indel branch reuses a separate bcfCallret and the SNP combine has
	// already emitted INFO/FMT/SCR onto the SNP row by the time the
	// indel branch runs.
	scrWanted := !isIndel && opts.FmtFlag&(B2BInfoSCR|B2BFmtSCR) != 0
	// nmWanted mirrors upstream's `fmt_flag & (B2B_FMT_NMBZ|B2B_INFO_NMBZ|
	// B2B_INFO_NM)` gate (bam2bcf.c:351, 442). When set, glfgen reads the
	// NM aux tag from each accepted read and accumulates the refNm/altNm
	// histograms used by NMBZ.
	nmWanted := opts.FmtFlag&(B2BFmtNMBZ|B2BInfoNMBZ|B2BInfoNM) != 0
	for i := range pile {
		p := &pile[i]
		// SCR accumulator runs BEFORE the is_refskip / is_unmap gate,
		// matching bam2bcf.c:300 which counts every pileup1_t entry
		// flagged with PLP_HAS_SOFT_CLIP regardless of refskip.
		if scrWanted && p.hasSoftClip {
			r.scr++
		}
		// is_refskip skips both branches unconditionally
		// (bam2bcf.c:301 `if (p->is_refskip || ...) continue;`).
		if p.isRefskip {
			continue
		}
		// is_del skips the SNP branch only; the indel branch lets the
		// read through (bam2bcf.c:307 `if (p->is_del && !is_indel) continue`).
		if p.isDel && !isIndel {
			continue
		}
		r.oriDepth++

		var b, q, baseQ, seqQ int
		if isIndel {
			// Indel branch (bam2bcf.c:312-421). Read the precomputed
			// (chosen-type, seqQ, indelQ) word produced by gap_prep.
			// Upstream's `seqQ = q = (p->aux & 0xff)` (bam2bcf.c:315):
			// both the local `seqQ` and `q` initialise from the indelQ
			// bits, NOT the saved seqQ bits — the latter are read into
			// `baseQ` later (line 420) and used for I16 only.
			b = int(p.aux>>16) & 0x3f
			q = int(p.aux) & 0xff
			seqQ = q
			N := len(pile)
			// CNS-only branch (bam2bcf.c:317-327): when --indels-cns is
			// active, the legacy REF-rescue heuristic is skipped (only the
			// !indels_v20 && !edlib path runs it). Instead, if NO read in
			// this sample has an indel in its CIGAR but the read carries
			// a non-zero indelQ from another sample's candidate, drop the
			// indelQ and use the raw base quality with seqQ = 99 (the
			// "basic sequence confidences" fallback).
			if opts.IndelsCNS {
				indelInSample := false
				for k := range pile {
					if pile[k].indel != 0 {
						indelInSample = true
						break
					}
				}
				if !indelInSample && (p.aux&0xff) != 0 {
					rawQ := 0
					if p.rec != nil && p.qpos >= 0 && p.qpos < len(p.rec.Qual) {
						rawQ = int(p.rec.Qual[p.qpos])
					}
					if rawQ > maxBQ {
						rawQ = maxBQ
					}
					q = rawQ
					seqQ = 99
				}
			} else {
				// Legacy REF-rescue heuristic for the !indels_v20 && !edlib
				// path (bam2bcf.c:338-348, originally e4e161068 fix for #1446).
				// At homopolymer / tandem-repeat sites bcf_cgp_compute_indelQ
				// often returns q == 0 for REF-type reads, so without this
				// rescue they would fail the min-baseQ gate below and be
				// undercounted in I16 / QS / AD. Upstream's heuristic: when
				// the read has no indel in its CIGAR (p.indel == 0) and
				// either the sample is shallow enough that q < _n/2 or deeper
				// than 20 reads, reclassify the read as REF (b = 0), promote
				// q to the read's raw base quality at qpos and rebuild seqQ
				// as a 3:2 blend of the old seqQ and the base quality. Cap
				// seqQ to 40 once the pile exceeds 20 reads to dampen
				// overconfident calls.
				if p.indel == 0 && (q < N/2 || N > 20) {
					rawQ := 0
					if p.rec != nil && p.qpos >= 0 && p.qpos < len(p.rec.Qual) {
						rawQ = int(p.rec.Qual[p.qpos])
					}
					b = 0
					q = rawQ
					seqQ = (3*seqQ + 2*rawQ) / 8
				}
				if N > 20 && seqQ > 40 {
					seqQ = 40
				}
			}
			// Upstream reads baseQ AFTER the heuristic from p->aux>>8&0xff
			// (bam2bcf.c:420) — i.e. the pre-heuristic seqQ bits stored by
			// bcf_cgp_compute_indelQ. The local `seqQ` variable is used
			// only to cap `q` below and to gate min-baseQ; I16 uses the
			// original aux value.
			baseQ = int(p.aux>>8) & 0xff
			if q < minBQ {
				// Low-quality indel read. Upstream's "is this a REF
				// match in disguise?" stash for --ambig-reads compensation
				// (bam2bcf.c:365-374): only reads with no indel in their
				// CIGAR (p.indel==0) and a real DNA base (b<4) qualify.
				if p.indel == 0 && b >= 0 && b < 4 {
					if p.reverse {
						adrRefMissed[b]++
					} else {
						adfRefMissed[b]++
					}
				}
				continue
			}
			// CNS-only seqQ cap and realigned-read dampener (TEST 6,
			// bam2bcf.c:386-415). Active only with --indels-cns.
			if opts.IndelsCNS {
				cap20 := N
				if cap20 > 20 {
					cap20 = 20
				}
				seqQOffset := opts.SeqQOffset
				if seqQOffset == 0 {
					seqQOffset = 120
				}
				if cs := seqQOffset - cap20*5; seqQ > cs {
					seqQ = cs
				}
				indelInSample := false
				for k := range pile {
					if pile[k].indel != 0 {
						indelInSample = true
						break
					}
				}
				if indelInSample && p.indel == 0 && b != 0 {
					alt := seqQ/2 + 5
					if alt < seqQ {
						seqQ = alt
					}
					rawQ := 0
					if p.rec != nil && p.qpos >= 0 && p.qpos < len(p.rec.Qual) {
						rawQ = int(p.rec.Qual[p.qpos])
					}
					a := rawQ/4 + 10
					bv := q/4 + 1
					min := a
					if bv < min {
						min = bv
					}
					if min < q {
						q = min
					}
				}
			}
		} else {
			b = p.base4
			// Lowest of this base and its neighbours' qualities, capped via
			// delta_baseQ (bam2bcf.c:~427-435).
			q = int(p.rawQual)
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
			baseQ = q
			// SNP branch hardcodes seqQ to 99 (bam2bcf.c:440).
			seqQ = b2bSeqQ
		}

		mapQ := int(p.mapq)
		if mapQ == 255 {
			mapQ = b2bDefMapQ
		}
		if mapQ == 0 {
			r.mq0++
		}
		// Cap sequence, literal port of bam2bcf.c:459-461:
		//   if (q > seqQ) q = seqQ;
		//   mapQ = mapQ < bca->capQ ? mapQ : bca->capQ;   // capQ = 60
		//   if (q > mapQ) q = mapQ;
		// SNP path: seqQ == 99 (b2bSeqQ). Indel path: seqQ is either the
		// stored bits 8-15 of p.aux or the REF-rescue blend above.
		if q > seqQ {
			q = seqQ
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
		if isIndel {
			// Indel branch (bam2bcf.c:350): is_diff = b ? 1 : 0. b == 0
			// is the REF indel type (length 0).
			if b != 0 {
				isDiff = 1
			}
		} else if !(ref4 < 4 && b == ref4) {
			isDiff = 1
		}
		// NM tag: upstream's get_aux_nm is called per accepted read when
		// any NMBZ/NM bit is set (bam2bcf.c:351, 442). Returns -1 if the
		// read has no NM tag; the per-read inm is used to bump the
		// per-sample refNm/altNm histograms below alongside the other
		// bias tallies. Note upstream subtracts 1 from inm for REF reads
		// (is_diff==0) and 2 for ALT reads.
		inm := -1
		if nmWanted {
			if v, ok := getAuxNm(p.rec, isDiff == 0); ok {
				inm = v
			}
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

		// Bias-test tallies (bam2bcf.c:491-537). baseQ/mapQ are clamped
		// to 59 so they index the 60-bin arrays; with nqual==60 the
		// nqual_over_60 scale factor is 1.0, hence imq==mapQ, ibq==baseQ.
		// epos / scLen were rescaled by get_position at accumulation
		// time. is_diff splits ref vs non-ref reads.
		biasBQ := baseQ
		if biasBQ > 59 {
			biasBQ = 59
		}
		biasMQ := mapQ
		if biasMQ > 59 {
			biasMQ = 59
		}
		if p.reverse {
			r.revMqs[biasMQ]++
		} else {
			r.fwdMqs[biasMQ]++
		}
		if isDiff == 0 {
			r.refPos[p.epos]++
			r.refBq[biasBQ]++
			r.refMq[biasMQ]++
			r.refScl[p.scLen]++
			if inm >= 0 {
				r.refNm[inm]++
			}
		} else {
			r.altPos[p.epos]++
			r.altBq[biasBQ]++
			r.altMq[biasMQ]++
			r.altScl[p.scLen]++
			if inm >= 0 {
				r.altNm[inm]++
			}
		}
	}

	// AD compensation for low-quality REF-looking indel reads
	// (bam2bcf.c:540-561). Only the indel branch produces ADR/ADF_ref_missed
	// entries; opts.AmbigReadsMode picks the strategy.
	if isIndel {
		switch opts.AmbigReadsMode {
		case AmbigReadsIncAD0:
			// All ambig reads are claimed as REF (allele 0).
			for j := 0; j < 4; j++ {
				r.adr[0] += adrRefMissed[j]
				r.adf[0] += adfRefMissed[j]
			}
		case AmbigReadsIncAD:
			// Distribute ambig reads across the per-allele AD slots in
			// proportion to the existing allele counts.
			dp, dpAmbig := 0, 0
			for j := 0; j < 4; j++ {
				dp += r.adr[j]
			}
			for j := 0; j < 4; j++ {
				dpAmbig += adrRefMissed[j]
			}
			if dp != 0 {
				for j := 0; j < 4; j++ {
					// upstream uses lroundf((float)dp_ambig * r->ADR[j] / dp).
					r.adr[j] += int(math.Round(float64(float32(dpAmbig) * float32(r.adr[j]) / float32(dp))))
				}
			}
			dp, dpAmbig = 0, 0
			for j := 0; j < 4; j++ {
				dp += r.adf[j]
			}
			for j := 0; j < 4; j++ {
				dpAmbig += adfRefMissed[j]
			}
			if dp != 0 {
				for j := 0; j < 4; j++ {
					r.adf[j] += int(math.Round(float64(float32(dpAmbig) * float32(r.adf[j]) / float32(dp))))
				}
			}
		}
	}

	// glfgen: errmod turns the packed base array into the 5x5 PL matrix.
	pl := r.p[:]
	em.Cal(bases, b2bMaxAlleles, pl)
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
	// hasAlt reports whether the site carries a real (non-`<*>`) ALT
	// allele; the bias annotations below are only computed when it does
	// (bam2bcf.c:1014, 1151).
	hasAlt bool
	// Bias annotations (bam2bcf.c calc_* functions). Each is paired with
	// an ok flag: upstream emits the INFO tag only when the value is not
	// HUGE_VAL, which Go represents as ok==false.
	vdb     float64
	vdbOK   bool
	segBias float64
	sgbOK   bool
	mwuPos  float64 // RPBZ
	rpbzOK  bool
	mwuMq   float64 // MQBZ
	mqbzOK  bool
	mwuBq   float64 // BQBZ
	bqbzOK  bool
	mwuMqs  float64 // MQSBZ
	mqsbzOK bool
	mwuSc   float64 // SCBZ
	scbzOK  bool
	mwuNm   float64 // NMBZ — Mann-Whitney U z-score on per-read NM tags
	nmbzOK  bool
	// scrTotal is the sum of per-sample SCR counts (INFO/SCR). scr is
	// the per-sample SCR array (FORMAT/SCR), in input sample order.
	// Both are filled by bcfCallCombine when either B2BInfoSCR or
	// B2BFmtSCR is set on the fmt_flag.
	scrTotal int
	scr      []int
	// dp4 holds the per-sample DP4 row: [fwdRef, revRef, fwdAlt, revAlt],
	// one [4]int per sample in input order. It mirrors upstream's
	// bc->DP4 (bam2bcf.c:1052-1058), sourced from calls[i].anno[0..3].
	// Always populated by bcfCallCombine — the per-sample FORMAT
	// DP/DV/DP4/SP renderers read it when the corresponding fmt_flag
	// bit is set.
	dp4 [][4]int
}

// bcfCallCombineIndel is the Go port of the `ref_base < 0` branch of
// bcf_call_combine (bam2bcf.c:959-1198). It combines the per-sample
// indel-flavored bcfCallret structs into a multi-allelic indel call:
//
//   - alleles[0] is forced to 0 (the REF indel-type index); subsequent
//     alleles are picked from indel-type indices 1..3 in descending QS,
//     stopping at the first zero-QS slot;
//   - no `<*>` unseen allele is appended;
//   - if only the REF type survives (nAlleles==1) the call is rejected
//     (returns ok=false), matching upstream's `return -1`;
//   - the multi-allelic PL grid is built from the 5x5 errmod matrix the
//     same way as for SNPs;
//   - VDB is recomputed from the per-sample altPos histogram (the indel
//     branch reused the same bca->alt_pos array, accumulating only over
//     indel-type-bearing reads);
//   - SGB is recomputed from the per-sample anno (now indel-flavored);
//   - RPBZ/MQBZ/SCBZ are recomputed from the per-call iref/ialt tallies
//     supplied via `bca`. BQBZ and MQSBZ are NOT recomputed: upstream's
//     `bcf_callaux_clean` only resets the bca histograms, not the
//     bcf_call_t.mwu_bq / mwu_mqs scalars, so those leak from the last
//     has_alt SNP combine at this or an earlier position. The caller
//     threads the leaked values in via the `leak` argument.
//   - NMBZ is recomputed from refNm/altNm summed across BOTH the SNP
//     and indel per-sample bcfCallret structs. Upstream accumulates
//     into the shared `bca->ref_nm/alt_nm` arrays in both glfgen
//     branches; `snpCalls` carries the SNP-pass contribution.
func bcfCallCombineIndel(calls []bcfCallret, snpCalls []bcfCallret, bca *bcfCallauxIndel, leak biasLeak) (bcfCall, bool) {
	var call bcfCall
	call.oriRef = -1

	// Coverage-normalised QS sums per indel-type slot (bam2bcf.c:970-977).
	// QS is indexed by indel-type index (0=REF, 1..3=ALT types in the
	// order chosen by bcfCgpComputeIndelQ).
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
	// (insertion sort, matching bam2bcf.c:980-984). Indel's "ref4" is 0.
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
	const indelRef4 = 0
	call.alleles[0] = indelRef4
	j := 1
	i := 3
	for ; i >= 0; i-- {
		ipos := idx[i]
		if ipos == indelRef4 {
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
	call.nAlleles = j
	if call.nAlleles == 1 {
		// No reliable supporting read (bam2bcf.c:1012).
		return call, false
	}
	call.hasAlt = true

	// Build the upper-triangle genotype-index list g[].
	x := call.nAlleles * (call.nAlleles + 1) / 2
	g := make([]int, 0, x)
	for ii := 0; ii < call.nAlleles; ii++ {
		for jj := 0; jj <= ii; jj++ {
			g = append(g, call.alleles[jj]*5+call.alleles[ii])
		}
	}

	// Per-sample PL.
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

	// Reorder QS and AD to match the site's allele ordering.
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
	call.dp4 = make([][4]int, len(calls))
	for i := range calls {
		r := &calls[i]
		call.depth += int(r.anno[0] + r.anno[1] + r.anno[2] + r.anno[3])
		call.oriDepth += r.oriDepth
		call.mq0 += r.mq0
		for k := 0; k < 16; k++ {
			call.anno[k] += r.anno[k]
		}
		call.dp4[i] = perSampleDP4(r)
	}

	// Sum the per-sample SNP-flavored bias tallies (refPos/altPos via
	// glfgen's indel branch — VDB is computed from alt_pos). Other
	// histograms (refMq/altMq/refBq/altBq) are accumulated but not
	// emitted on the indel row: the indel branch only consults the
	// iref_*/ialt_* histograms threaded via bca.
	var altPos [b2bNpos]int
	for i := range calls {
		r := &calls[i]
		for k := 0; k < b2bNpos; k++ {
			altPos[k] += r.altPos[k]
		}
	}

	// NM histograms for the indel-pass NMBZ. Upstream's bca->ref_nm /
	// bca->alt_nm are wiped by bcf_callaux_clean between the SNP and
	// indel passes (mpileup.c:580, 593 — the clean call before the
	// indel glfgen at line 593 clears the shared B2B_N_NM slot via the
	// `else` branch in bcf_callaux_clean, bam2bcf.c:219-223). So the
	// indel-row NMBZ at bam2bcf.c:1185 sees ONLY the indel-pass NM
	// tallies, never the SNP-pass ones — sum just `calls` here.
	var refNm, altNm [b2bNNm]int
	for i := range calls {
		r := &calls[i]
		for k := 0; k < b2bNNm; k++ {
			refNm[k] += r.refNm[k]
			altNm[k] += r.altNm[k]
		}
	}

	// SGB from the (indel-flavored) per-sample anno.
	call.segBias, call.sgbOK = calcSegBias(calls, &call)

	// RPBZ / MQBZ / SCBZ from the indel iref/ialt histograms populated by
	// bcfCallGapPrep (bam2bcf.c:1166-1171).
	call.mwuPos, call.rpbzOK = calcMWUBiasZ(bca.IrefPos, bca.IaltPos, false)
	call.mwuMq, call.mqbzOK = calcMWUBiasZ(bca.IrefMq, bca.IaltMq, true)
	call.mwuSc, call.scbzOK = calcMWUBiasZ(bca.IrefScl[:], bca.IaltScl[:], false)
	// NMBZ (bam2bcf.c:1184-1185) uses the same calc_mwu_biasZ helper on
	// the cross-pass refNm/altNm sum computed above. Note upstream does
	// NOT skip this for indel rows.
	call.mwuNm, call.nmbzOK = calcMWUBiasZ(refNm[:], altNm[:], false)

	// BQBZ / MQSBZ leak from the last has_alt SNP combine, per upstream's
	// bcf_callaux_clean which only resets the bca arrays.
	call.mwuBq, call.bqbzOK = leak.bq, leak.bqOK
	call.mwuMqs, call.mqsbzOK = leak.mqs, leak.mqsOK

	// VDB from the indel-branch alt_pos (bam2bcf.c:1195).
	call.vdb, call.vdbOK = calcVDB(altPos[:])

	// SCR: the indel pass does not tally soft-clipped reads (only the
	// SNP branch does, see bcfCallGlfgen's scrWanted gate). Upstream
	// emits SCR on the indel row using the same per-column tally as
	// the SNP row, which lives on the shared bcf_call_t.scr / .scrTotal.
	// Reuse the SNP-pass per-sample counts so both rows agree.
	call.scr = make([]int, len(snpCalls))
	for i := range snpCalls {
		call.scr[i] = snpCalls[i].scr
		call.scrTotal += snpCalls[i].scr
	}

	return call, true
}

// biasLeak carries the BQBZ / MQSBZ scalars that upstream's bcf_call_t
// retains across positions: bcf_callaux_clean wipes the bca histograms
// but not the call's mwu_bq / mwu_mqs, so an indel record at a position
// where the SNP row had no real ALT inherits whichever value was last
// computed by a has_alt SNP combine.
type biasLeak struct {
	bq, mqs     float64
	bqOK, mqsOK bool
}

// updateBiasLeak captures the BQBZ / MQSBZ scalars from a SNP combine
// for later reuse by the indel combine at the same or a subsequent
// position. Only has_alt SNP combines overwrite the slot — upstream's
// `if (!has_alt) return 0;` (bam2bcf.c:1151) preserves the previous
// value otherwise.
func (l *biasLeak) update(call *bcfCall) {
	if !call.hasAlt {
		return
	}
	l.bq, l.bqOK = call.mwuBq, call.bqbzOK
	l.mqs, l.mqsOK = call.mwuMqs, call.mqsbzOK
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
	call.dp4 = make([][4]int, len(calls))
	for i := range calls {
		r := &calls[i]
		call.depth += int(r.anno[0] + r.anno[1] + r.anno[2] + r.anno[3])
		call.oriDepth += r.oriDepth
		call.mq0 += r.mq0
		for j := 0; j < 16; j++ {
			call.anno[j] += r.anno[j]
		}
		call.dp4[i] = perSampleDP4(r)
	}

	// Fold per-sample SCR tallies into the call (bam2bcf.c:1060-1066).
	// Upstream stores SCR as an int32_t[ncall+1]: index 0 holds the
	// total (INFO/SCR), indices 1..n hold the per-sample counts
	// (FORMAT/SCR). Here we keep the total scalar separately and a
	// per-sample slice.
	call.scr = make([]int, len(calls))
	for i := range calls {
		call.scr[i] = calls[i].scr
		call.scrTotal += calls[i].scr
	}

	// has_alt: a real ALT exists unless the only non-REF allele is `<*>`
	// (bam2bcf.c:1014). The bias tests are skipped otherwise.
	call.hasAlt = !(call.nAlleles == 2 && call.unseen != -1)
	if !call.hasAlt {
		return call
	}

	// Sum the per-sample bias tallies into the combined arrays, matching
	// upstream's shared bca accumulation.
	var refPos, altPos, refScl, altScl [b2bNpos]int
	var refMq, altMq, refBq, altBq, fwdMqs, revMqs [b2bNqual]int
	var refNm, altNm [b2bNNm]int
	for i := range calls {
		r := &calls[i]
		for k := 0; k < b2bNpos; k++ {
			refPos[k] += r.refPos[k]
			altPos[k] += r.altPos[k]
			refScl[k] += r.refScl[k]
			altScl[k] += r.altScl[k]
		}
		for k := 0; k < b2bNqual; k++ {
			refMq[k] += r.refMq[k]
			altMq[k] += r.altMq[k]
			refBq[k] += r.refBq[k]
			altBq[k] += r.altBq[k]
			fwdMqs[k] += r.fwdMqs[k]
			revMqs[k] += r.revMqs[k]
		}
		for k := 0; k < b2bNNm; k++ {
			refNm[k] += r.refNm[k]
			altNm[k] += r.altNm[k]
		}
	}

	// SGB — segregation bias (bam2bcf.c:1158, calc_SegBias).
	call.segBias, call.sgbOK = calcSegBias(calls, &call)

	// Mann-Whitney U z-scores (bam2bcf.c:1172-1183). RPBZ uses the
	// read-position bins, MQBZ the mapping-quality bins (one-sided),
	// BQBZ the base-quality bins, MQSBZ the fwd/rev mapping-quality
	// bins, SCBZ the soft-clip-length bins.
	call.mwuPos, call.rpbzOK = calcMWUBiasZ(refPos[:], altPos[:], false)
	call.mwuMq, call.mqbzOK = calcMWUBiasZ(refMq[:], altMq[:], true)
	call.mwuBq, call.bqbzOK = calcMWUBiasZ(refBq[:], altBq[:], false)
	call.mwuMqs, call.mqsbzOK = calcMWUBiasZ(fwdMqs[:], revMqs[:], false)
	call.mwuSc, call.scbzOK = calcMWUBiasZ(refScl[:], altScl[:], false)
	// NMBZ — NM-bias z-score on the per-read NM-tag histograms
	// (bam2bcf.c:1184-1185). The same calc_mwu_biasZ helper applies.
	call.mwuNm, call.nmbzOK = calcMWUBiasZ(refNm[:], altNm[:], false)

	// VDB — variant distance bias (bam2bcf.c:1194, calc_vdb).
	call.vdb, call.vdbOK = calcVDB(altPos[:])

	return call
}

// getPositionResult bundles the soft-clip-aware position annotations
// computed by getPosition for one read base.
type getPositionResult struct {
	pos    int // edist: query position counted from the unclipped read start
	length int // unclipped read length (l_qseq minus soft-clips)
	scLen  int // length of the nearer soft-clip, 0 if none
	scDist int // distance from this base to the nearer soft-clip
}

// getPosition is the Go port of get_position (bam2bcf.c:146). cigarOps
// and cigarLens describe rec's CIGAR (one entry per op); qpos is the
// 0-based query position of the base; lQseq is the read's sequence
// length. The CIGAR op codes use the BAM convention 0=M..9, with 4=S
// (soft-clip) and 5=H (hard-clip).
func getPosition(cigarOps, cigarLens []int, qpos, lQseq int) getPositionResult {
	const (
		cigarSoftClip = 4
		cigarHardClip = 5
	)
	edist := qpos + 1
	scLeft, scRight := 0, 0
	scLeftDist, scRightDist := -1, -1

	// Leading soft-clip run.
	i := 0
	for ; i < len(cigarOps); i++ {
		switch cigarOps[i] {
		case cigarHardClip:
			continue
		case cigarSoftClip:
			scLeft += cigarLens[i]
		default:
			goto leftDone
		}
	}
leftDone:
	if scLeft != 0 {
		scLeftDist = qpos + 1 - scLeft
	}
	edist -= scLeft

	// Trailing soft-clip run (down to i, the first non-leading-clip op).
trailing:
	for j := len(cigarOps) - 1; j >= i; j-- {
		switch cigarOps[j] {
		case cigarHardClip:
			continue
		case cigarSoftClip:
			scRight += cigarLens[j]
		default:
			break trailing
		}
	}
	if scRight != 0 {
		scRightDist = lQseq - scRight - qpos
	}

	res := getPositionResult{pos: edist, length: lQseq - scLeft - scRight}
	switch {
	case scLeftDist >= 0:
		if scRightDist < 0 || scLeftDist < scRightDist {
			res.scLen = scLeft
			res.scDist = scLeftDist
		}
	case scRightDist >= 0:
		res.scLen = scRight
		res.scDist = scRightDist
	}
	return res
}

// biasPositionBins reduces getPosition's raw output to the rescaled
// (epos, scLen) bins glfgen's bias tallies index, exactly as the inline
// code at bam2bcf.c:497-504 does.
func biasPositionBins(r getPositionResult) (epos, scLen int) {
	epos = int(float64(r.pos) / float64(r.length+1) * float64(b2bNpos-1))
	if r.scLen != 0 {
		scLen = int(15.0 * float64(r.scLen) / float64(r.scDist+1))
		if scLen > 99 {
			scLen = 99
		}
	}
	if epos < 0 {
		epos = 0
	} else if epos >= b2bNpos {
		epos = b2bNpos - 1
	}
	return epos, scLen
}

// getAuxNm is the Go port of get_aux_nm (bam2bcf.c:96). It reads the
// BAM NM:i: aux tag, normalises it by treating each indel as a single
// event (subtracting len-1 for indel ops longer than 1) and adding the
// soft-clip lengths as mismatches, then subtracts 1 for a REF read or
// 2 for an ALT read (mirroring upstream's MNP-aware adjustment), and
// clamps the result to [0, b2bNNm-1]. Returns ok=false when the read
// has no NM tag (mirroring upstream's "-1" sentinel).
//
// rec must be the originating BAM record; isRef is true when this read
// supports the reference allele (b == ref4 for SNPs, b == 0 for indels).
// The upstream qpos argument is unused by upstream — the cache it would
// normally key on is per-read, and this port computes everything fresh.
func getAuxNm(rec *sam.Record, isRef bool) (int, bool) {
	if rec == nil {
		return 0, false
	}
	aux, ok := rec.GetAux("NM")
	if !ok {
		return 0, false
	}
	v, ok := aux.Int()
	if !ok {
		return 0, false
	}
	nm := int(v)
	// Adjust: indels collapse to single events (upstream subtracts
	// len-1 per indel op of length > 1), soft-clips count as
	// mismatches (added to nm).
	for _, op := range rec.Cigar {
		o := op.Op()
		l := int(op.Length())
		switch o {
		case sam.CigarSoftClip:
			nm += l
		case sam.CigarInsertion, sam.CigarDeletion:
			if l > 1 {
				nm -= l - 1
			}
		}
	}
	// MNP-aware subtraction (bam2bcf.c:135-137): REF reads have one
	// "expected" mismatch nearby on average, ALT reads two.
	if isRef {
		nm--
	} else {
		nm -= 2
	}
	if nm < 0 {
		nm = 0
	}
	if nm >= b2bNNm {
		nm = b2bNNm - 1
	}
	return nm, true
}

// kfErfc is the Go port of htslib's kf_erfc (kfunc.c:58). bcftools uses
// this rational approximation, not the libc erfc, so calc_vdb must call
// it for byte-for-byte parity. Note kf_erfc(x) effectively evaluates
// erfc(x*sqrt2) — the M_SQRT2 scaling is part of the function.
func kfErfc(x float64) float64 {
	const (
		p0 = 220.2068679123761
		p1 = 221.2135961699311
		p2 = 112.0792914978709
		p3 = 33.912866078383
		p4 = 6.37396220353165
		p5 = .7003830644436881
		p6 = .03526249659989109
		q0 = 440.4137358247522
		q1 = 793.8265125199484
		q2 = 637.3336333788311
		q3 = 296.5642487796737
		q4 = 86.78073220294608
		q5 = 16.06417757920695
		q6 = 1.755667163182642
		q7 = .08838834764831844
	)
	z := math.Abs(x) * math.Sqrt2
	if z > 37 {
		if x > 0 {
			return 0
		}
		return 2
	}
	expntl := math.Exp(z * z * -.5)
	var p float64
	if z < 10/math.Sqrt2 {
		p = expntl * ((((((p6*z+p5)*z+p4)*z+p3)*z+p2)*z+p1)*z + p0) /
			(((((((q7*z+q6)*z+q5)*z+q4)*z+q3)*z+q2)*z+q1)*z + q0)
	} else {
		p = expntl / 2.506628274631001 / (z + 1/(z+2/(z+3/(z+4/(z+.65)))))
	}
	if x > 0 {
		return 2 * p
	}
	return 2 * (1 - p)
}

// calcVDB is the Go port of calc_vdb (bam2bcf.c:600). It returns the
// variant distance bias (a value in [0,1], smaller meaning more biased)
// for the alt-read position histogram. ok is false when VDB cannot be
// computed because fewer than two variant reads were observed (upstream
// returns HUGE_VAL).
//
// Faithfulness note: upstream declares mean_pos and mean_diff as C
// `float`, so the running sums and the final mean are single precision.
// The port keeps that — calc_vdb is sensitive to it via the `int ipos`
// truncation and the dp==2 exact formula.
func calcVDB(pos []int) (float64, bool) {
	const readlen = 100
	// param is the upstream nparam-by-3 fitting table {dp, pscale, pshift}.
	param := [15][3]float32{
		{3, 0.079, 18}, {4, 0.09, 19.8}, {5, 0.1, 20.5}, {6, 0.11, 21.5},
		{7, 0.125, 21.6}, {8, 0.135, 22}, {9, 0.14, 22.2}, {10, 0.153, 22.3},
		{15, 0.19, 22.8}, {20, 0.22, 23.2}, {30, 0.26, 23.4}, {40, 0.29, 23.5},
		{50, 0.35, 23.65}, {100, 0.5, 23.7}, {200, 0.7, 23.7},
	}
	const nparam = 15

	dp := 0
	var meanPos float32
	for i := range pos {
		if pos[i] == 0 {
			continue
		}
		dp += pos[i]
		meanPos += float32(pos[i] * i)
	}
	if dp < 2 {
		return 0, false // one or zero reads can be placed anywhere
	}
	meanPos /= float32(dp)

	var meanDiff float32
	for i := range pos {
		if pos[i] == 0 {
			continue
		}
		// abs(i - mean_pos) in single precision, matching the C float math.
		meanDiff += float32(pos[i]) * float32(math.Abs(float64(float32(i)-meanPos)))
	}
	meanDiff /= float32(dp)

	if dp == 2 {
		ipos := int(meanDiff) // C float-to-int truncation
		// Upstream: (int expr)/(readlen-1) is integer division, only the
		// final /(readlen*0.5) is floating point.
		num := (2*readlen - 2*(ipos+1) - 1) * (ipos + 1)
		return float64(num/(readlen-1)) / (readlen * 0.5), true
	}

	var i int
	if dp >= 200 {
		i = nparam
	} else {
		for i = 0; i < nparam; i++ {
			if param[i][0] >= float32(dp) {
				break
			}
		}
	}
	var pscale, pshift float32
	switch {
	case i == nparam:
		pscale = param[nparam-1][1]
		pshift = param[nparam-1][2]
	case i > 0 && param[i][0] != float32(dp):
		pscale = (param[i-1][1] + param[i][1]) * 0.5
		pshift = (param[i-1][2] + param[i][2]) * 0.5
	default:
		pscale = param[i][1]
		pshift = param[i][2]
	}
	return 0.5 * kfErfc(-float64((meanDiff-pshift)*pscale)), true
}

// logsumexp2 is the Go port of bam2bcf.c's logsumexp2 helper.
func logsumexp2(a, b float64) float64 {
	if a > b {
		return math.Log(1+math.Exp(b-a)) + a
	}
	return math.Log(1+math.Exp(a-b)) + b
}

// calcSegBias is the Go port of calc_SegBias (bam2bcf.c:895). It returns
// the segregation-bias score and ok=false (HUGE_VAL upstream) when no
// non-reference reads were observed.
func calcSegBias(calls []bcfCallret, call *bcfCall) (float64, bool) {
	n := len(calls)
	if n == 0 {
		return 0, false
	}
	nr := int(call.anno[2] + call.anno[3]) // non-reference reads
	if nr == 0 {
		return 0, false
	}
	avgDp := int(call.anno[0]+call.anno[1]) + nr
	avgDp /= n
	if avgDp == 0 {
		// Guard against a divide-by-zero that C would not hit because
		// nr>0 implies the total is positive; keep parity by flooring.
		avgDp = 1
	}
	m := math.Floor(float64(nr)/float64(avgDp) + 0.5)
	if m > float64(n) {
		m = float64(n)
	} else if m == 0 {
		m = 1
	}
	f := m / 2.0 / float64(n)
	p := float64(nr) / float64(n)
	q := float64(nr) / m
	var sum float64
	const log2 = math.Ln2
	for i := range calls {
		oi := int(calls[i].anno[2] + calls[i].anno[3])
		var tmp float64
		if oi != 0 {
			tmp = logsumexp2(math.Log(2*(1-f)), math.Log(f)+float64(oi)*log2-q)
			tmp += math.Log(f) + float64(oi)*math.Log(q/p) - q + p
		} else {
			tmp = math.Log(2*f*(1-f)*math.Exp(-q)+f*f*math.Exp(-2*q)+(1-f)*(1-f)) + p
		}
		sum += tmp
	}
	return sum, true
}

// calcMWUBiasZ is the Go port of calc_mwu_biasZ (bam2bcf.c:817) with
// do_Z=1, i.e. the standard-deviation-normalised Mann-Whitney U z-score
// used by RPBZ/MQBZ/BQBZ/MQSBZ/SCBZ. leftOnly mirrors the upstream
// left_only argument (set for MQBZ); with do_Z it does not change the
// result but is kept for fidelity. ok=false signals HUGE_VAL (one of
// the two histograms is empty).
func calcMWUBiasZ(a, b []int, leftOnly bool) (float64, bool) {
	n := len(a)
	_ = leftOnly // unused with do_Z=1, kept to mirror the C signature

	// Optimisation: detect an all-zero b array.
	bEmpty := true
	for i := 0; i < n; i++ {
		if b[i] != 0 {
			bEmpty = false
			break
		}
	}

	var e, l, na, nb int
	var t int64
	if bEmpty {
		for i := n - 1; i >= 0; i-- {
			na += a[i]
			ai := int64(a[i])
			t += (ai*ai - 1) * ai
		}
	} else {
		for i := n - 1; i >= 0; i-- {
			e += a[i] * b[i]
			l += a[i] * nb
			na += a[i]
			nb += b[i]
			pp := int64(a[i] + b[i])
			t += (pp*pp - 1) * pp
		}
	}
	if na == 0 || nb == 0 {
		return 0, false
	}

	u := float64(l) + float64(e)*0.5
	m := float64(na) * float64(nb) / 2.0
	var2 := float64(na*nb) / 12.0 *
		(float64(na+nb+1) - float64(t)/float64((na+nb)*(na+nb-1)))
	if var2 <= 0 {
		return 0, true
	}
	return (u - m) / math.Sqrt(var2), true
}
