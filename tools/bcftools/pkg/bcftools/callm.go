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

// mcallTin groups the parsed inputs the caller needs per record.
type mcallTin struct {
	nals    int         // number of alleles including REF (and <*> if present)
	ngts    int         // nals*(nals+1)/2
	ploidy  int         // 1 or 2
	nsmpl   int         // number of samples
	pdg     [][]float64 // per-sample probability vector, length ngts
	qsum    []float64   // normalized allele frequencies, length nals
	unseen  int         // index of the <*> allele, or -1
	thetaLn float64     // log(theta * Watterson factor)
}

// computeTheta replicates mcall.c init: theta scaled by the Watterson
// factor aM over the total number of alleles, then logged.
func computeTheta(prior float64, ploidy, nsmpl int) float64 {
	if prior <= 0 {
		return 0
	}
	n := ploidy * nsmpl
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
func parseMcallInputs(v *vcf.Variant, opts CallOptions) *mcallTin {
	ploidy := 2
	if opts.Ploidy == PloidyHaploid {
		ploidy = 1
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
		nals:   nals,
		ngts:   ngts,
		ploidy: ploidy,
		nsmpl:  len(v.Samples),
		unseen: unseen,
		pdg:    make([][]float64, len(v.Samples)),
	}
	in.thetaLn = computeTheta(opts.Prior, ploidy, len(v.Samples))

	for i, s := range v.Samples {
		pls, ok := decodePLInts(s.Data["PL"], ngts)
		if !ok {
			in.pdg[i] = make([]float64, ngts) // all-zero -> treated as missing
			continue
		}
		in.pdg[i] = setPdg(pls, ngts, nals, in.unseen)
	}

	// INFO/QS -> qsum (raw allele frequencies). Missing trailing alleles
	// (typical ref-only <*> site) get qsum=0.
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
	return in
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
// bitmask of the most likely combination, plus ref_lk, lk_sum, max_lk.
func mcallBestAlleles(in *mcallTin) (alsMask int, refLk, lkSum, maxLk float64) {
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

	// Single allele.
	for ia := 0; ia < nals; ia++ {
		lkTot := 0.0
		lkSet := false
		iaa := (ia+1)*(ia+2)/2 - 1
		for is := 0; is < in.nsmpl; is++ {
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
			if in.qsum[ia] == 0 {
				continue
			}
			iaa := (ia+1)*(ia+2)/2 - 1
			for ib := 0; ib < ia; ib++ {
				if in.qsum[ib] == 0 {
					continue
				}
				lkTot := 0.0
				lkSet := false
				fa := in.qsum[ia] / (in.qsum[ia] + in.qsum[ib])
				fb := in.qsum[ib] / (in.qsum[ia] + in.qsum[ib])
				fa2, fb2, fab := fa*fa, fb*fb, 2*fa*fb
				ibb := (ib+1)*(ib+2)/2 - 1
				iab := iaa - ia + ib
				for is := 0; is < in.nsmpl; is++ {
					p := in.pdg[is]
					var val float64
					if in.ploidy == 2 {
						val = fa2*p[iaa] + fb2*p[ibb] + fab*p[iab]
					} else {
						val = fa*p[iaa] + fb*p[ibb]
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
				update(1<<ia|1<<ib, lkTot, lkSet, lkSet)
			}
		}
	}

	// Three alleles.
	if nals > 2 {
		for ia := 0; ia < nals; ia++ {
			if in.qsum[ia] == 0 {
				continue
			}
			iaa := (ia+1)*(ia+2)/2 - 1
			for ib := 0; ib < ia; ib++ {
				if in.qsum[ib] == 0 {
					continue
				}
				ibb := (ib+1)*(ib+2)/2 - 1
				iab := iaa - ia + ib
				for ic := 0; ic < ib; ic++ {
					if in.qsum[ic] == 0 {
						continue
					}
					lkTot := 0.0
					lkSet := false
					tot := in.qsum[ia] + in.qsum[ib] + in.qsum[ic]
					fa, fb, fc := in.qsum[ia]/tot, in.qsum[ib]/tot, in.qsum[ic]/tot
					fa2, fb2, fc2 := fa*fa, fb*fb, fc*fc
					fab, fac, fbc := 2*fa*fb, 2*fa*fc, 2*fb*fc
					icc := (ic+1)*(ic+2)/2 - 1
					iac := iaa - ia + ic
					ibc := ibb - ib + ic
					for is := 0; is < in.nsmpl; is++ {
						p := in.pdg[is]
						var val float64
						if in.ploidy == 2 {
							val = fa2*p[iaa] + fb2*p[ibb] + fc2*p[icc] + fab*p[iab] + fac*p[iac] + fbc*p[ibc]
						} else {
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
func mcallSite(v *vcf.Variant, opts CallOptions) (out *vcf.Variant, keep bool, ok bool) {
	if !hasQS(v) {
		return nil, false, false
	}
	in := parseMcallInputs(v, opts)
	if in == nil {
		return nil, false, false
	}

	alsMask, refLk, lkSum, maxLk := mcallBestAlleles(in)
	ninf := math.Inf(-1)
	var maxQual float64 = ninf
	if maxLk != ninf {
		maxQual = -4.343 * (refLk - logsumexp2(lkSum, refLk))
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
	// -A/keep-unseen requests it. Without keep-unseen, the <*> allele is
	// always dropped.
	nalsOri := in.nals
	keepUnseen := false // we do not implement --keep-unseen-allele here
	newMask := 0
	for i := 0; i < nalsOri; i++ {
		if i > 0 && i == in.unseen && !keepUnseen {
			continue
		}
		if opts.KeepAlts {
			alsMask |= 1 << i
		}
		if alsMask&(1<<i) != 0 {
			newMask |= 1 << i
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
			if pdgAllZero(in.pdg[i]) || in.ploidy == 0 {
				gts[i] = [2]int{-1, -1}
			} else {
				gts[i] = [2]int{0, 0}
				ac[0] += in.ploidy
			}
		}
	} else {
		for i := 0; i < in.nsmpl; i++ {
			g := mcallCallGenotype(in, alsMask, alsMap, i)
			gts[i] = g
			if g[0] >= 0 {
				ac[g[0]]++
			}
			if g[1] >= 0 {
				ac[g[1]]++
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

	out.Samples = make([]vcf.Sample, in.nsmpl)
	for i, s := range v.Samples {
		ns := vcf.Sample{Name: s.Name, Data: copyStringMap(s.Data)}
		// GT
		ns.Data["GT"] = renderGT(gts[i], in.ploidy)
		// AD (Number=R) re-index
		if ad, okAD := s.Data["AD"]; okAD {
			ns.Data["AD"] = reindexNumberR(ad, in.nals, alsMap)
		}
		// PL re-index or drop
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
		out.Samples[i] = ns
	}

	// FORMAT order: upstream keeps the input FORMAT order with GT
	// prepended, then PL removed for ref-only sites.
	out.Format = rebuildFormat(v.Format, keepPL)

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
// Missing samples return {-1,-1}.
func mcallCallGenotype(in *mcallTin, alsMask int, alsMap []int, ismpl int) [2]int {
	pdg := in.pdg[ismpl]
	ploidy := in.ploidy
	if ploidy == 0 || pdgAllZero(pdg) {
		return [2]int{-1, -1}
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
			lk = pdg[iaa] * in.qsum[ia] * in.qsum[ia]
		} else {
			lk = pdg[iaa] * in.qsum[ia]
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
				lk := 2 * pdg[iab] * in.qsum[ia] * in.qsum[ib]
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
