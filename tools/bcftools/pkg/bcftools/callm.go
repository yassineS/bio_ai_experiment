package bcftools

import (
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// This file is a faithful Go port of the multiallelic genotype-likelihood
// caller in reference_code/bcftools/mcall.c (the `bcftools call -m`
// algorithm). It is exercised when the mpileup INFO/QS annotation is
// present on the record; the synthetic PL-only fixtures used by the older
// tests fall through to the heuristic path in call.go.
//
// The port covers, for the single-sample-group case that the
// `mpileup | call -m` pipeline produces:
//
//   - set_pdg: PL -> genotype probabilities (call_init_pl2p table).
//   - mcall_find_best_alleles: the EM-style allele-frequency search that
//     scores every 1/2/3-allele combination against the per-genotype
//     probabilities weighted by QS-derived allele frequencies.
//   - the per-site QUAL from the ref-vs-nonref probability
//     (-4.343 * (ref_lk - logsumexp2(lk_sum, ref_lk))).
//   - mcall_set_ref_genotypes / mcall_call_genotypes: per-sample GT.
//   - the allele-trimming maps that drop the unseen <*> allele and any
//     unsupported ALTs, re-indexing PL (Number=G) and AD (Number=R).
//   - the INFO rewrite: AN/AC from the called genotypes, DP4/MQ from the
//     mpileup INFO/I16, and the removal of I16/QS.

// pl2pTable mirrors call_init_pl2p: pl2p[i] = 10^(-i/10) for i in [0,256).
var pl2pTable = func() [256]float64 {
	var t [256]float64
	for i := 0; i < 256; i++ {
		t[i] = math.Pow(10, -float64(i)/10.0)
	}
	return t
}()

func pl2p(pl int) float64 {
	if pl < 256 {
		return pl2pTable[pl]
	}
	return math.Pow(10, -float64(pl)/10.0)
}

// (logsumexp2 lives in bam2bcf.go and is shared.)

// mcallGroup holds one sample-group's per-record state: the sample
// indices that participate, the per-group qsum, and the per-group
// best-allele results that mcallBestAlleles fills in (mirrors
// upstream smpl_grp_t). A solo group covering every sample is the
// default when -G is not supplied.
type mcallGroup struct {
	samples []int
	qsum    []float64
	// Filled in by mcallBestAlleles.
	alsMask int
	refLk   float64
	lkSum   float64
	maxLk   float64
}

// mcallTin groups the parsed inputs the caller needs per record.
type mcallTin struct {
	nals   int // number of alleles including REF (and <*> if present)
	ngts   int // nals*(nals+1)/2
	ploidy int // representative ploidy (max across samples) — used
	//  to choose the EM formula and the per-genotype
	//  HWE weighting on the dominant branch. With
	//  per-sample ploidy (e.g. --ploidy GRCh37 on chrY)
	//  smplPloidy[i] is the authoritative value for
	//  sample i.
	nsmpl      int         // number of samples
	pdg        [][]float64 // per-sample probability vector, length ngts
	qsum       []float64   // normalized allele frequencies (pooled-group view)
	groups     []mcallGroup
	unseen     int     // index of the <*> allele, or -1
	thetaLn    float64 // log(theta * Watterson factor)
	smplPloidy []int   // per-sample ploidy (0/1/2). When nil, every
	// sample uses the global `ploidy` field.
}

// computeTheta replicates mcall.c init: theta scaled by the Watterson
// factor aM over the total number of alleles, then logged.
func computeTheta(prior float64, ploidy, nsmpl int) float64 {
	return computeThetaN(prior, ploidy*nsmpl)
}

// computeThetaN is the n-aware variant used when per-sample ploidy
// may differ (n = sum of per-sample ploidies).
func computeThetaN(prior float64, n int) float64 {
	if prior <= 0 {
		return 0
	}
	aM := 1.0
	for i := 2; i < n; i++ {
		aM += 1.0 / float64(i)
	}
	t := prior * aM
	if t >= 1 {
		t = 0.99
	}
	return math.Log(t)
}

// gtIndex returns the diploid PL index for the unordered genotype (a,b).
func gtIndex(a, b int) int {
	if a < b {
		a, b = b, a
	}
	return a*(a+1)/2 + b
}

// hasQS reports whether the record carries the mpileup INFO/QS annotation
// that the mcall path requires.
func hasQS(v *vcf.Variant) bool {
	_, ok := v.Info["QS"]
	return ok
}

// parseMcallInputs reads PL, QS, the allele list, and the unseen-allele
// index off v. It returns nil when the record lacks the data the mcall
// path needs (e.g. no PL or no QS), signalling the caller to fall back.
//
// samplePloidy, when non-nil, overrides opts.Ploidy on a per-sample
// basis. The site-level `ploidy` field is set to the max of those
// values, which selects which EM branch (haploid vs diploid HWE) runs
// — the per-sample value is consulted in genotype assignment, AC/AN,
// and the per-sample HWE weighting via in.smplPloidy.
func parseMcallInputs(v *vcf.Variant, opts CallOptions, samplePloidy []int) *mcallTin {
	ploidy := 2
	if opts.Ploidy == PloidyHaploid {
		ploidy = 1
	}
	if len(samplePloidy) > 0 {
		// The site ploidy is the maximum across samples — that picks
		// the diploid EM branch when any sample is diploid, which
		// matches mcall.c's per-site uniform formula choice on the
		// pdg/qsum loops. Per-sample ploidy is still honoured below
		// (smplPloidy + genotype calling).
		mx := 0
		for _, p := range samplePloidy {
			if p > mx {
				mx = p
			}
		}
		if mx > 0 {
			ploidy = mx
		}
	}
	nAlt := len(v.Alt)
	if nAlt == 1 && v.Alt[0] == "." {
		nAlt = 0
	}
	nals := nAlt + 1
	ngts := nals * (nals + 1) / 2

	// Locate the unseen <*> / <X> allele.
	unseen := -1
	for i, a := range v.Alt {
		if a == "<*>" || a == "<X>" || a == "X" || a == "*" {
			unseen = i + 1 // +1: REF is allele 0
			break
		}
	}

	in := &mcallTin{
		nals:       nals,
		ngts:       ngts,
		ploidy:     ploidy,
		nsmpl:      len(v.Samples),
		unseen:     unseen,
		pdg:        make([][]float64, len(v.Samples)),
		smplPloidy: samplePloidy,
	}
	// Watterson aM in mcall.c uses the per-sample ploidy at INIT time
	// (which is `ploidy_max(table)` across every sample, not the
	// per-record value). With a PloidyTable that's ploidy_max; without
	// it that's the global Ploidy.
	totAlleles := ploidy * len(v.Samples)
	if opts.PloidyTable != nil {
		totAlleles = opts.PloidyTable.MaxPloidy() * len(v.Samples)
	}
	in.thetaLn = computeThetaN(opts.Prior, totAlleles)

	for i, s := range v.Samples {
		pls, ok := decodePLInts(s.Data["PL"], ngts)
		if !ok {
			in.pdg[i] = make([]float64, ngts) // all-zero -> treated as missing
			continue
		}
		in.pdg[i] = setPdg(pls, ngts, nals, in.unseen)
	}

	// INFO/QS -> qsum (raw allele frequencies). Missing trailing alleles
	// (typical ref-only <*> site) get qsum=0. This pooled view is kept
	// as the always-available default; per-group qsums (when -G is
	// active) replace it via buildGroupQsums below.
	qs := parseFloatList(v.Info["QS"])
	in.qsum = make([]float64, nals)
	for i := 0; i < nals && i < len(qs); i++ {
		in.qsum[i] = qs[i]
	}
	// normalize so qsum sums to 1
	sum := 0.0
	for _, q := range in.qsum {
		sum += q
	}
	if sum != 0 {
		for i := range in.qsum {
			in.qsum[i] /= sum
		}
	}
	// Build per-group qsums. When no resolved groups are supplied the
	// caller still sees a single group covering every sample, with
	// the pooled qsum cloned in. Mirrors mcall.c's nsmpl_grp==1
	// fast path.
	in.groups = buildMcallGroups(v, opts, in.nals, in.qsum)
	// -F AN,AC: when prior allele frequencies are supplied via INFO
	// tags, reweight every group's qsum by (qsum + 0.5*ac) / (nsmpl
	// + 0.5*an), then re-normalise so each group's qsum sums to 1.
	// Matches mcall.c:1498-1528.
	if opts.PriorAN != "" && opts.PriorAC != "" {
		applyPriorFreqs(v, in, opts.PriorAN, opts.PriorAC)
	}
	return in
}

// applyPriorFreqs reweights every group's qsum using the INFO/<AN>
// and INFO/<AC> values supplied via `-F AN,AC`. The reweighting
// formula mirrors upstream mcall.c:1511-1518.
func applyPriorFreqs(v *vcf.Variant, in *mcallTin, anTag, acTag string) {
	if v == nil || in == nil {
		return
	}
	anStr, ok := v.Info[anTag]
	if !ok {
		return
	}
	an, err := strconv.Atoi(strings.TrimSpace(anStr))
	if err != nil || an <= 0 {
		return
	}
	acStr, ok := v.Info[acTag]
	if !ok {
		return
	}
	ac := parseIntList(acStr)
	if len(ac) != in.nals-1 {
		return
	}
	ac0 := an
	for i := 0; i < in.nals-1; i++ {
		ac0 -= ac[i]
	}
	if ac0 < 0 {
		return
	}
	for gi := range in.groups {
		grp := &in.groups[gi]
		nsmpl := float64(len(grp.samples))
		denom := nsmpl + 0.5*float64(an)
		if denom == 0 || len(grp.qsum) < in.nals {
			continue
		}
		for i := 0; i < in.nals-1; i++ {
			grp.qsum[i+1] = (grp.qsum[i+1] + 0.5*float64(ac[i])) / denom
		}
		grp.qsum[0] = (grp.qsum[0] + 0.5*float64(ac0)) / denom
		// Re-normalise (mcall.c:1521-1528).
		s := 0.0
		for _, q := range grp.qsum {
			s += q
		}
		if s != 0 {
			for j := range grp.qsum {
				grp.qsum[j] /= s
			}
		}
	}
}

// buildMcallGroups partitions samples per -G and recomputes the per-
// group qsum from the configured tag (default: AD; QS when present).
// When opts.sampleGroups is nil, a single group covering every
// sample with the pooled qsum is returned — keeping the original
// single-group code path bit-equal.
func buildMcallGroups(v *vcf.Variant, opts CallOptions, nals int, pooledQS []float64) []mcallGroup {
	if opts.sampleGroups == nil {
		g := mcallGroup{samples: make([]int, len(v.Samples))}
		for i := range v.Samples {
			g.samples[i] = i
		}
		g.qsum = append([]float64(nil), pooledQS...)
		return []mcallGroup{g}
	}
	resolved, err := opts.sampleGroups.Resolve(headerSamplesFromVariant(v))
	if err != nil {
		// Resolution errors should have been caught at Call() init;
		// fall back to single-group to avoid panic in the hot loop.
		g := mcallGroup{samples: make([]int, len(v.Samples))}
		for i := range v.Samples {
			g.samples[i] = i
		}
		g.qsum = append([]float64(nil), pooledQS...)
		return []mcallGroup{g}
	}
	// nsmpl_grp == 1 is mcall.c's fast path: even with -G supplied,
	// when the file resolves to a single group upstream skips the
	// AD-based recomputation and uses INFO/QS directly (mcall.c:1446).
	// Match the same behaviour so single-group -G byte-equals the
	// no-group default.
	if len(resolved) == 1 {
		g := mcallGroup{samples: append([]int(nil), resolved[0].Indices...)}
		g.qsum = append([]float64(nil), pooledQS...)
		return []mcallGroup{g}
	}
	// Choose the per-sample tag: explicit override wins, else QS, else AD.
	tag := opts.GroupSamplesTag
	if tag == "" {
		tag = opts.sampleGroups.Tag
	}
	if tag == "" {
		// QS is only emitted at the INFO level by mpileup; the
		// FORMAT-level fallback is AD. When -G is active mpileup
		// must have been run with `-a AD` (or `-a QS`).
		tag = "AD"
	}
	out := make([]mcallGroup, len(resolved))
	for gi, rg := range resolved {
		g := mcallGroup{samples: append([]int(nil), rg.Indices...)}
		g.qsum = make([]float64, nals)
		for _, isidx := range rg.Indices {
			if isidx >= len(v.Samples) {
				continue
			}
			raw, ok := v.Samples[isidx].Data[tag]
			if !ok {
				continue
			}
			vals := parseFloatList(raw)
			sum := 0.0
			for _, x := range vals {
				if x > 0 {
					sum += x
				}
			}
			if sum == 0 {
				continue
			}
			for j := 0; j < nals && j < len(vals); j++ {
				if vals[j] > 0 {
					g.qsum[j] += vals[j] / sum
				}
			}
		}
		// Normalize per-group so qsum sums to 1 (mcall.c:1523-1528).
		s := 0.0
		for _, q := range g.qsum {
			s += q
		}
		if s != 0 {
			for j := range g.qsum {
				g.qsum[j] /= s
			}
		}
		out[gi] = g
	}
	return out
}

// headerSamplesFromVariant returns the sample names declared on v.
func headerSamplesFromVariant(v *vcf.Variant) []string {
	out := make([]string, len(v.Samples))
	for i, s := range v.Samples {
		out[i] = s.Name
	}
	return out
}

// decodePLInts parses a "0,15,100" PL string into ints. Missing entries
// ("." / "") become a sentinel matching bcf_int32_missing handling.
func decodePLInts(s string, ngts int) ([]int, bool) {
	if s == "" || s == "." {
		return nil, false
	}
	parts := strings.Split(s, ",")
	out := make([]int, len(parts))
	for i, p := range parts {
		if p == "." || p == "" {
			out[i] = plMissing
			continue
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			return nil, false
		}
		out[i] = n
	}
	return out, true
}

const plMissing = -2147483647 // stand-in for bcf_int32_missing

// setPdg ports set_pdg for the common diploid case: it converts PLs to
// probabilities and normalizes them to sum to 1. For the unseen-allele
// filling we follow the unseen>=0 branch (mpileup always provides <*>).
func setPdg(pls []int, ngts, nals, unseen int) []float64 {
	pdg := make([]float64, ngts)
	// A PL vector shorter than ngts is the text-VCF analogue of
	// bcf_int32_vector_end on a diploid record — upstream treats it
	// as "expect diploid GLs; if not diploid treat as missing", which
	// collapses the whole sample to all-zero pdg.
	if len(pls) < ngts {
		return pdg
	}
	sum := 0.0
	j := 0
	for ; j < ngts && j < len(pls); j++ {
		if pls[j] == plMissing {
			break
		}
		pdg[j] = pl2p(pls[j])
		sum += pdg[j]
	}
	if j == 0 {
		// all missing
		return make([]float64, ngts)
	}
	if j < ngts {
		// Missing values present; fill from the unseen allele's PLs.
		if unseen < 0 {
			sum = 0
			for k := 0; k < ngts; k++ {
				val := 255
				if k < len(pls) && pls[k] != plMissing {
					val = pls[k]
				}
				pdg[k] = pl2p(val)
				sum += pdg[k]
			}
		} else {
			sum = 0
			k := 0
			plsCopy := make([]int, ngts)
			copy(plsCopy, pls)
			for len(plsCopy) < ngts {
				plsCopy = append(plsCopy, plMissing)
			}
			for ia := 0; ia < nals; ia++ {
				for ib := 0; ib <= ia; ib++ {
					if plsCopy[k] == plMissing {
						idx := gtIndex(ia, unseen)
						if idx >= len(plsCopy) || plsCopy[idx] == plMissing {
							idx = gtIndex(ib, unseen)
						}
						if idx >= len(plsCopy) || plsCopy[idx] == plMissing {
							idx = gtIndex(unseen, unseen)
						}
						if idx >= len(plsCopy) || plsCopy[idx] == plMissing {
							plsCopy[k] = 255
						} else {
							plsCopy[k] = plsCopy[idx]
						}
					}
					pdg[k] = pl2p(plsCopy[k])
					sum += pdg[k]
					k++
				}
			}
		}
	}
	if sum == float64(ngts) {
		// all missing
		return make([]float64, ngts)
	}
	for k := 0; k < ngts; k++ {
		pdg[k] /= sum
	}
	return pdg
}

// mcallBestAlleles ports mcall_find_best_alleles. It returns the allele
// bitmask of the most likely combination, plus ref_lk, lk_sum, max_lk
// for the supplied sample group. When grp is nil the pooled view
// (in.groups[0] when populated, else every sample with the global
// qsum) is used.
func mcallBestAlleles(in *mcallTin, grp *mcallGroup) (alsMask int, refLk, lkSum, maxLk float64) {
	if grp == nil {
		if len(in.groups) > 0 {
			grp = &in.groups[0]
		} else {
			// Build a transient pooled view (only triggered by very
			// old callers that bypass parseMcallInputs).
			tmp := mcallGroup{samples: make([]int, in.nsmpl), qsum: in.qsum}
			for i := 0; i < in.nsmpl; i++ {
				tmp.samples[i] = i
			}
			grp = &tmp
		}
	}
	qsum := grp.qsum
	if qsum == nil {
		qsum = in.qsum
	}
	nals := in.nals
	ninf := math.Inf(-1)
	refLk, maxLk, lkSum = ninf, ninf, ninf
	alsMask = 0

	// update ports the UPDATE_MAX_LKs macro. The max-likelihood selection
	// uses lkSet (the local lk_tot_set), while the lk_sum accumulation uses
	// the separate sumFlag the macro is invoked with (for single alleles
	// this is `ia>0 && lk_tot_set`, so the reference homozygote does not
	// feed lk_sum but can still win the max comparison).
	update := func(als int, lkTot float64, lkSet, sumFlag bool) {
		if maxLk < lkTot && lkSet {
			maxLk = lkTot
			alsMask = als
		}
		if sumFlag {
			lkSum = logsumexp2(lkTot, lkSum)
		}
	}

	// smplP returns sample is' ploidy, falling back to the site-level
	// value when smplPloidy isn't supplied. Ploidy=0 samples are
	// IGNORED only in genotype assignment (mcall.c set_pdg ignores
	// ploidy entirely, and the EM scoring loops feed every sample's
	// pdg into lk_tot regardless of ploidy — only the per-sample HWE
	// formula picks haploid vs diploid).
	smplP := func(is int) int {
		if is < len(in.smplPloidy) {
			return in.smplPloidy[is]
		}
		return in.ploidy
	}

	// Single allele.
	for ia := 0; ia < nals; ia++ {
		lkTot := 0.0
		lkSet := false
		iaa := (ia+1)*(ia+2)/2 - 1
		for _, is := range grp.samples {
			p := in.pdg[is]
			if iaa < len(p) && p[iaa] != 0 {
				lkTot += math.Log(p[iaa])
				lkSet = true
			}
		}
		if ia == 0 {
			refLk = lkTot
		} else {
			lkTot += in.thetaLn
		}
		update(1<<ia, lkTot, lkSet, ia > 0 && lkSet)
	}

	// Two alleles.
	if nals > 1 {
		for ia := 0; ia < nals; ia++ {
			if qsum[ia] == 0 {
				continue
			}
			iaa := (ia+1)*(ia+2)/2 - 1
			for ib := 0; ib < ia; ib++ {
				if qsum[ib] == 0 {
					continue
				}
				lkTot := 0.0
				lkSet := false
				fa := qsum[ia] / (qsum[ia] + qsum[ib])
				fb := qsum[ib] / (qsum[ia] + qsum[ib])
				fa2, fb2, fab := fa*fa, fb*fb, 2*fa*fb
				ibb := (ib+1)*(ib+2)/2 - 1
				iab := iaa - ia + ib
				for _, is := range grp.samples {
					p := in.pdg[is]
					pl := smplP(is)
					var val float64
					if pl == 2 || (len(in.smplPloidy) == 0 && in.ploidy == 2) {
						val = fa2*p[iaa] + fb2*p[ibb] + fab*p[iab]
					} else if pl == 1 {
						val = fa*p[iaa] + fb*p[ibb]
					}
					// pl == 0 (or any other value) → val stays 0,
					// matching upstream which leaves val=0 when neither
					// the haploid nor diploid branch matches.
					if val != 0 {
						lkTot += math.Log(val)
						lkSet = true
					}
				}
				if ia != 0 {
					lkTot += in.thetaLn
				}
				if ib != 0 {
					lkTot += in.thetaLn
				}
				update(1<<ia|1<<ib, lkTot, lkSet, lkSet)
			}
		}
	}

	// Three alleles.
	if nals > 2 {
		for ia := 0; ia < nals; ia++ {
			if qsum[ia] == 0 {
				continue
			}
			iaa := (ia+1)*(ia+2)/2 - 1
			for ib := 0; ib < ia; ib++ {
				if qsum[ib] == 0 {
					continue
				}
				ibb := (ib+1)*(ib+2)/2 - 1
				iab := iaa - ia + ib
				for ic := 0; ic < ib; ic++ {
					if qsum[ic] == 0 {
						continue
					}
					lkTot := 0.0
					lkSet := false
					tot := qsum[ia] + qsum[ib] + qsum[ic]
					fa, fb, fc := qsum[ia]/tot, qsum[ib]/tot, qsum[ic]/tot
					fa2, fb2, fc2 := fa*fa, fb*fb, fc*fc
					fab, fac, fbc := 2*fa*fb, 2*fa*fc, 2*fb*fc
					icc := (ic+1)*(ic+2)/2 - 1
					iac := iaa - ia + ic
					ibc := ibb - ib + ic
					for _, is := range grp.samples {
						p := in.pdg[is]
						pl := smplP(is)
						var val float64
						if pl == 2 || (len(in.smplPloidy) == 0 && in.ploidy == 2) {
							val = fa2*p[iaa] + fb2*p[ibb] + fc2*p[icc] + fab*p[iab] + fac*p[iac] + fbc*p[ibc]
						} else if pl == 1 {
							val = fa*p[iaa] + fb*p[ibb] + fc*p[icc]
						}
						if val != 0 {
							lkTot += math.Log(val)
							lkSet = true
						}
					}
					if ia != 0 {
						lkTot += in.thetaLn
					}
					if ib != 0 {
						lkTot += in.thetaLn
					}
					if ic != 0 {
						lkTot += in.thetaLn
					}
					update(1<<ia|1<<ib|1<<ic, lkTot, lkSet, lkSet)
				}
			}
		}
	}

	return alsMask, refLk, lkSum, maxLk
}

// mcallSite runs the full multiallelic caller on v. It returns the
// rewritten record and a keep flag. ok is false when the record lacks the
// QS/PL data the algorithm needs (caller should fall back).
func mcallSite(v *vcf.Variant, opts CallOptions, samplePloidy []int) (out *vcf.Variant, keep bool, ok bool) {
	if !hasQS(v) {
		return nil, false, false
	}
	in := parseMcallInputs(v, opts, samplePloidy)
	if in == nil {
		return nil, false, false
	}

	ninf := math.Inf(-1)
	var maxQual float64 = ninf
	var refLk, lkSum float64 = ninf, ninf
	alsMask := 0
	// Per-group best-allele scan (mcall.c:1538-1554). For each group
	// the unioned mask drives the site's final allele set; the qual
	// reported in INFO is the maximum across groups, with the
	// matching group's refLk and lkSum carried along so the QUAL
	// posterior is computed from the same group's likelihoods.
	for gi := range in.groups {
		gMask, gRefLk, gLkSum, gMaxLk := mcallBestAlleles(in, &in.groups[gi])
		in.groups[gi].alsMask = gMask
		in.groups[gi].refLk = gRefLk
		in.groups[gi].lkSum = gLkSum
		in.groups[gi].maxLk = gMaxLk
		alsMask |= gMask
		if gMaxLk == ninf {
			continue
		}
		q := -4.343 * (gRefLk - logsumexp2(gLkSum, gRefLk))
		if maxQual == ninf || q > maxQual {
			maxQual = q
			refLk = gRefLk
			lkSum = gLkSum
		}
	}

	// Make sure REF is always present.
	if alsMask&1 == 0 {
		alsMask |= 1
	}

	// is_variant reflects the natural best-allele result *before* any -A
	// forcing (mcall.c line 1558). Under -A a ref-only site keeps its ALT
	// columns but every genotype is still set to reference.
	isVariant := alsMask != 1

	// Build the new allele set, dropping the unseen allele unless
	// -A/keep-unseen requests it. With keep-unseen on a ref-only
	// site upstream forcibly adds the unseen allele (mcall.c:1564-
	// 1568) so the output retains `<*>` and the per-sample PL
	// vector against the unseen genotype.
	nalsOri := in.nals
	keepUnseen := opts.KeepUnseen
	newMask := 0
	addedCount := 0
	for i := 0; i < nalsOri; i++ {
		// CALL_KEEP_UNSEEN: if we have exactly the REF in the new
		// mask and hit the unseen index, force-include it
		// (mcall.c:1564 — `i==unseen && nals_new==1`).
		if keepUnseen && i == in.unseen && addedCount == 1 {
			newMask |= 1 << i
			addedCount++
			continue
		}
		if i > 0 && i == in.unseen {
			continue
		}
		if opts.KeepAlts {
			alsMask |= 1 << i
		}
		if alsMask&(1<<i) != 0 {
			newMask |= 1 << i
			addedCount++
		}
	}
	if newMask&1 == 0 {
		newMask |= 1
	}

	// als_map: old allele index -> new index (or -1 if dropped).
	alsMap := make([]int, nalsOri)
	nalsNew := 0
	for i := 0; i < nalsOri; i++ {
		if newMask&(1<<i) != 0 {
			alsMap[i] = nalsNew
			nalsNew++
		} else {
			alsMap[i] = -1
		}
	}

	// The genotype caller runs only when the site is genuinely variant
	// (best alleles included a non-ref) and survives trimming. A ref-only
	// site (als_new==1) or an -A-forced non-variant site uses ref GTs.
	callGenotypes := isVariant && nalsNew > 1

	// Per-sample genotypes + AC.
	ac := make([]int, nalsNew)
	gts := make([][2]int, in.nsmpl)
	nAC := 0
	if !callGenotypes {
		// All-ref: GT 0/0 (or . if no data).
		for i := 0; i < in.nsmpl; i++ {
			pl := in.ploidy
			if i < len(in.smplPloidy) {
				pl = in.smplPloidy[i]
			}
			if pdgAllZero(in.pdg[i]) || pl == 0 {
				gts[i] = [2]int{-1, -1}
				if pl == 0 {
					// Ploidy-0 sample (e.g. F on chrY): emit "." with
					// no allele, matching mcall.c's ploidy==0 branch.
				}
			} else if pl == 1 {
				gts[i] = [2]int{0, -1}
				ac[0]++
			} else {
				gts[i] = [2]int{0, 0}
				ac[0] += 2
			}
		}
	} else {
		// Iterate per-group so each sample uses its group's local
		// alsMask + qsum for genotype assignment (mcall.c:1610).
		for gi := range in.groups {
			grp := &in.groups[gi]
			for _, is := range grp.samples {
				g := mcallCallGenotype(in, alsMask, alsMap, is, grp)
				gts[is] = g
				if g[0] >= 0 {
					ac[g[0]]++
				}
				if g[1] >= 0 {
					ac[g[1]]++
				}
			}
		}
		for i := 1; i < nalsNew; i++ {
			nAC += ac[i]
		}
	}

	if isVariant && opts.VariantsOnly && nAC == 0 {
		return nil, false, true
	}
	if !isVariant && opts.VariantsOnly {
		return nil, false, true
	}

	// PL is retained whenever any ALT survives on output (nals_new>1),
	// matching mcall.c which only strips PL for the pure ref-only branch.
	keepPL := nalsNew > 1
	out = mcallEmit(v, in, opts, alsMap, nalsNew, gts, ac, nAC, maxQual, lkSum, refLk, keepPL)
	return out, true, true
}

// pdgAllZero reports whether every entry of p is zero (the "missing
// sample" condition mcall uses).
func pdgAllZero(p []float64) bool {
	for _, x := range p {
		if x != 0 {
			return false
		}
	}
	return true
}

// callAnnotateFlags captures which `-a` tags the caller requested.
// Mirrors upstream's CALL_FMT_GQ / CALL_FMT_GP / CALL_FMT_PV4 bit
// flags (vcfcall.c:868-871 and call.h:37-39).
type callAnnotateFlags struct {
	gq, gp, pv4 bool
}

func parseCallAnnotateFlags(spec string) callAnnotateFlags {
	out := callAnnotateFlags{}
	if spec == "" {
		return out
	}
	for _, tok := range strings.Split(spec, ",") {
		t := strings.ToUpper(strings.TrimSpace(tok))
		switch t {
		case "GQ", "FORMAT/GQ", "FMT/GQ":
			out.gq = true
		case "GP", "FORMAT/GP", "FMT/GP":
			out.gp = true
		case "PV4", "INFO/PV4":
			out.pv4 = true
		}
	}
	return out
}

// computeCallGQ ports the upstream `-a GQ` formula: per-sample
// genotype quality is `-4.34294 * log(1 - max/sum)`, where `gps`
// are the per-genotype likelihoods (per-sample, restricted to the
// group's allowed alleles and the per-sample ploidy). The integer
// result is capped at INT8_MAX=127 (mcall.c:881).
//
// Returns "." for missing samples (no pdg signal, or all gps==0).
func computeCallGQ(in *mcallTin, grp *mcallGroup, alsMap []int, ismpl int) string {
	if in == nil || ismpl >= len(in.pdg) {
		return "."
	}
	pdg := in.pdg[ismpl]
	if pdgAllZero(pdg) {
		return "."
	}
	ploidy := in.ploidy
	if ismpl < len(in.smplPloidy) {
		ploidy = in.smplPloidy[ismpl]
	}
	if ploidy == 0 {
		return "."
	}
	qsum := in.qsum
	alsMask := -1
	if grp != nil {
		alsMask = grp.alsMask
		if grp.qsum != nil {
			qsum = grp.qsum
		}
	}
	maxVal, sumVal := 0.0, 0.0
	visit := func(g float64) {
		if g > maxVal {
			maxVal = g
		}
		sumVal += g
	}
	for ia := 0; ia < in.nals; ia++ {
		if alsMask >= 0 && alsMask&(1<<ia) == 0 {
			continue
		}
		if alsMap != nil && ia < len(alsMap) && alsMap[ia] < 0 {
			continue
		}
		iaa := (ia+1)*(ia+2)/2 - 1
		if ploidy == 2 {
			visit(pdg[iaa] * qsum[ia] * qsum[ia])
			for ib := 0; ib < ia; ib++ {
				if alsMask >= 0 && alsMask&(1<<ib) == 0 {
					continue
				}
				if alsMap != nil && ib < len(alsMap) && alsMap[ib] < 0 {
					continue
				}
				iab := iaa - ia + ib
				visit(2 * pdg[iab] * qsum[ia] * qsum[ib])
			}
		} else if ploidy == 1 {
			visit(pdg[iaa] * qsum[ia])
		}
	}
	if maxVal <= 0 || sumVal <= 0 {
		return "0"
	}
	frac := 1 - maxVal/sumVal
	// frac can underflow to 0 when one genotype completely dominates
	// (pdg ≈ delta function). Upstream's log(1 - max/sum) becomes
	// -Inf in that case; the int cast that follows then clamps via
	// INT8_MAX (mcall.c:881). Mirror with an explicit cap.
	if frac <= 0 {
		return "127"
	}
	gq := -4.34294 * math.Log(frac)
	gqi := int(gq)
	if gqi > 127 {
		gqi = 127
	}
	if gqi < 0 {
		gqi = 0
	}
	return strconv.Itoa(gqi)
}

// appendUnique appends s to keys unless it's already present.
func appendUnique(keys []string, s string) []string {
	for _, k := range keys {
		if k == s {
			return keys
		}
	}
	return append(keys, s)
}

// parseFloatList parses a comma-separated float list ("1,0" or
// "0.165116,0.834884,0"). Non-numeric entries (e.g. ".") become 0.
func parseFloatList(s string) []float64 {
	if s == "" || s == "." {
		return nil
	}
	parts := strings.Split(s, ",")
	out := make([]float64, len(parts))
	for i, p := range parts {
		if p == "" || p == "." {
			out[i] = 0
			continue
		}
		f, err := strconv.ParseFloat(p, 64)
		if err != nil {
			out[i] = 0
			continue
		}
		out[i] = f
	}
	return out
}

// parseIntList parses a comma-separated int list ("5,0,0,0"). Non-numeric
// entries become 0.
func parseIntList(s string) []int {
	if s == "" || s == "." {
		return nil
	}
	parts := strings.Split(s, ",")
	out := make([]int, len(parts))
	for i, p := range parts {
		n, err := strconv.Atoi(strings.TrimSpace(p))
		if err != nil {
			out[i] = 0
			continue
		}
		out[i] = n
	}
	return out
}

// mcallEmit assembles the output record: trimmed alleles, called GT, the
// re-indexed PL/AD, the QUAL, and the INFO rewrite (AN/AC/DP4/MQ; drop
// I16/QS). It is the Go analogue of the tail of mcall().
func mcallEmit(v *vcf.Variant, in *mcallTin, opts CallOptions, alsMap []int, nalsNew int, gts [][2]int, ac []int, nAC int, maxQual, lkSum, refLk float64, keepPL bool) *vcf.Variant {
	out := *v
	ninf := math.Inf(-1)

	// --- alleles -------------------------------------------------------
	newAlt := make([]string, 0, nalsNew-1)
	for i := 1; i < in.nals; i++ {
		if alsMap[i] >= 0 {
			newAlt = append(newAlt, v.Alt[i-1])
		}
	}
	out.Alt = newAlt
	if len(out.Alt) == 0 {
		out.Alt = []string{"."}
	}

	// --- PL re-index (Number=G), only for variant sites; ref-only sites
	// drop PL entirely (matching mcall: PL removed when als_new==1). ----
	ngtsNew := nalsNew * (nalsNew + 1) / 2
	// pl_map[new k] -> old l
	plMap := make([]int, 0, ngtsNew)
	{
		l := 0
		for i := 0; i < in.nals; i++ {
			for j := 0; j <= i; j++ {
				if newMaskHas(alsMap, i) && newMaskHas(alsMap, j) {
					plMap = append(plMap, l)
				}
				l++
			}
		}
	}

	annot := parseCallAnnotateFlags(opts.Annotate)
	// Build a per-sample → *mcallGroup index for the GQ formula
	// (each group's local qsum + alsMask drives the per-sample
	// genotype-posterior calculation).
	sampleGrp := make([]*mcallGroup, in.nsmpl)
	for gi := range in.groups {
		grp := &in.groups[gi]
		for _, is := range grp.samples {
			if is >= 0 && is < in.nsmpl {
				sampleGrp[is] = grp
			}
		}
	}
	out.Samples = make([]vcf.Sample, in.nsmpl)
	for i, s := range v.Samples {
		ns := vcf.Sample{Name: s.Name, Data: copyStringMap(s.Data)}
		// GT
		pl := in.ploidy
		if i < len(in.smplPloidy) {
			pl = in.smplPloidy[i]
		}
		ns.Data["GT"] = renderGT(gts[i], pl)
		// AD (Number=R) re-index
		if ad, okAD := s.Data["AD"]; okAD {
			ns.Data["AD"] = reindexNumberR(ad, in.nals, alsMap)
		}
		// PL re-index or drop.
		if keepPL {
			if _, okPL := s.Data["PL"]; okPL {
				pls, ok := decodePLInts(s.Data["PL"], in.ngts)
				if ok {
					nb := make([]string, len(plMap))
					for k, l := range plMap {
						if l < len(pls) && pls[l] != plMissing {
							nb[k] = strconv.Itoa(pls[l])
						} else {
							nb[k] = "."
						}
					}
					ns.Data["PL"] = strings.Join(nb, ",")
				}
			}
		} else {
			delete(ns.Data, "PL")
		}
		// Upstream emits GQ only on variant sites where PL is
		// retained; matches mcall.c:845.
		if annot.gq && keepPL {
			ns.Data["GQ"] = computeCallGQ(in, sampleGrp[i], alsMap, i)
		}
		out.Samples[i] = ns
	}

	// FORMAT order: upstream keeps the input FORMAT order with GT
	// prepended, then PL removed for ref-only sites. When -a GQ is
	// set AND PL is retained (variant site), GQ is appended after
	// PL — mirroring mcall.c which only emits GQ alongside PL.
	out.Format = rebuildFormat(v.Format, keepPL)
	if annot.gq && keepPL {
		out.Format = appendUnique(out.Format, "GQ")
	}

	// --- QUAL ----------------------------------------------------------
	if nAC != 0 {
		out.Qual = quantizeQual(maxQual)
	} else {
		if lkSum != ninf {
			out.Qual = quantizeQual(-4.343 * (lkSum - logsumexp2(lkSum, refLk)))
		} else if ac[0] != 0 {
			if in.thetaLn != 0 {
				out.Qual = quantizeQual(-4.343 * in.thetaLn)
			} else {
				out.Qual = 0
			}
		} else {
			out.Qual = -1 // missing
		}
	}

	// --- INFO rewrite --------------------------------------------------
	out.Info = copyStringMap(v.Info)
	out.InfoOrder = append([]string(nil), v.InfoOrder...)

	// Parse I16 before removing it.
	i16 := parseFloatList(v.Info["I16"])

	delInfo(&out, "I16")
	delInfo(&out, "QS")

	// AC (Number=A), only when nals_new>1.
	if nalsNew > 1 {
		acStrs := make([]string, nalsNew-1)
		for i := 1; i < nalsNew; i++ {
			acStrs[i-1] = strconv.Itoa(ac[i])
		}
		setInfo(&out, "AC", strings.Join(acStrs, ","))
	}
	// AN.
	an := nAC + ac[0]
	setInfo(&out, "AN", strconv.Itoa(an))

	// DP4 + MQ from I16.
	if len(i16) >= 16 {
		dp4 := []string{
			strconv.Itoa(int(i16[0])),
			strconv.Itoa(int(i16[1])),
			strconv.Itoa(int(i16[2])),
			strconv.Itoa(int(i16[3])),
		}
		setInfo(&out, "DP4", strings.Join(dp4, ","))
		dp4sum := int(i16[0]) + int(i16[1]) + int(i16[2]) + int(i16[3])
		mq := 0
		if dp4sum > 0 {
			mq = (int(i16[8]) + int(i16[10])) / dp4sum
		}
		setInfo(&out, "MQ", strconv.Itoa(mq))
	}

	return &out
}

// newMaskHas reports whether old allele index i survives trimming.
func newMaskHas(alsMap []int, i int) bool {
	return i < len(alsMap) && alsMap[i] >= 0
}

// renderGT renders a called (a,b) genotype in new allele numbering.
func renderGT(g [2]int, ploidy int) string {
	if g[0] < 0 {
		if ploidy == 2 {
			return "./."
		}
		return "."
	}
	if ploidy == 1 || g[1] < 0 {
		return strconv.Itoa(g[0])
	}
	a, b := g[0], g[1]
	return strconv.Itoa(a) + "/" + strconv.Itoa(b)
}

// reindexNumberR re-orders a Number=R field (e.g. AD) by als_map and drops
// dropped alleles. nalsOri is the original allele count.
func reindexNumberR(field string, nalsOri int, alsMap []int) string {
	vals := parseIntList(field)
	nNew := 0
	for _, m := range alsMap {
		if m >= 0 {
			nNew++
		}
	}
	out := make([]string, nNew)
	for k := 0; k < nalsOri && k < len(alsMap); k++ {
		l := alsMap[k]
		if l < 0 {
			continue
		}
		if k < len(vals) {
			out[l] = strconv.Itoa(vals[k])
		} else {
			out[l] = "0"
		}
	}
	return strings.Join(out, ",")
}

// rebuildFormat returns the FORMAT key order: GT first, then the input
// keys, dropping PL when keepPL is false.
func rebuildFormat(orig []string, keepPL bool) []string {
	out := []string{"GT"}
	for _, k := range orig {
		if k == "GT" {
			continue
		}
		if k == "PL" && !keepPL {
			continue
		}
		out = append(out, k)
	}
	return out
}

// delInfo removes key from the INFO map and its order slice.
func delInfo(v *vcf.Variant, key string) {
	if v.Info != nil {
		delete(v.Info, key)
	}
	for i, k := range v.InfoOrder {
		if k == key {
			v.InfoOrder = append(v.InfoOrder[:i], v.InfoOrder[i+1:]...)
			break
		}
	}
}

// quantizeQual rounds a QUAL through the float32 + %g(6) path htslib uses
// so the VCF writer's shortest-float formatting reproduces upstream
// exactly.
func quantizeQual(q float64) float64 {
	if q < 0 {
		q = 0
	}
	s := formatFloat32G(q)
	f, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return q
	}
	return f
}

// mcallCallGenotype ports mcall_call_genotypes for one sample, returning
// the called (a,b) alleles in the *new* (trimmed) allele numbering.
// Missing samples return {-1,-1}. When grp is non-nil the group's
// own alsMask and qsum gate the per-allele scoring, mirroring
// upstream's per-group genotype call; passing nil uses the supplied
// alsMask plus the pooled qsum (backwards-compat for any callers
// that haven't been updated).
func mcallCallGenotype(in *mcallTin, alsMask int, alsMap []int, ismpl int, grp *mcallGroup) [2]int {
	pdg := in.pdg[ismpl]
	ploidy := in.ploidy
	if ismpl < len(in.smplPloidy) {
		ploidy = in.smplPloidy[ismpl]
	}
	if ploidy == 0 || pdgAllZero(pdg) {
		return [2]int{-1, -1}
	}
	qsum := in.qsum
	if grp != nil {
		alsMask = grp.alsMask
		if grp.qsum != nil {
			qsum = grp.qsum
		}
	}

	// Default fallback: 0/0.
	g0, g1 := 0, -1
	if ploidy == 2 {
		g1 = 0
	}

	bestLk := 0.0
	for ia := 0; ia < in.nals; ia++ {
		if alsMask&(1<<ia) == 0 {
			continue
		}
		iaa := (ia+1)*(ia+2)/2 - 1
		var lk float64
		if ploidy == 2 {
			lk = pdg[iaa] * qsum[ia] * qsum[ia]
		} else {
			lk = pdg[iaa] * qsum[ia]
		}
		if bestLk < lk {
			bestLk = lk
			g0 = alsMap[ia]
		}
	}
	if ploidy == 2 {
		g1 = g0
		for ia := 0; ia < in.nals; ia++ {
			if alsMask&(1<<ia) == 0 {
				continue
			}
			iaa := (ia+1)*(ia+2)/2 - 1
			for ib := 0; ib < ia; ib++ {
				if alsMask&(1<<ib) == 0 {
					continue
				}
				iab := iaa - ia + ib
				lk := 2 * pdg[iab] * qsum[ia] * qsum[ib]
				if bestLk < lk {
					bestLk = lk
					g0 = alsMap[ib]
					g1 = alsMap[ia]
				}
			}
		}
	}
	return [2]int{g0, g1}
}
