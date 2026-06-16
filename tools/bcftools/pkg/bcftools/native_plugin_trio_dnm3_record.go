// Per-record orchestration for the trio-dnm3 float models: reads FORMAT/PL, AD,
// QS, QM, SP, sets up the per-member [father,mother,child] arrays, applies the
// >4-allele trimming, runs the selected model, and writes FORMAT/DNM, FORMAT/VA
// and FORMAT/VAF. This ports process_record() (the non-NAIVE path), set_trio_PL
// / set_trio_QS / set_trio_QM / set_trio_QS_noisy, many_alts_trim and the score
// -> tag conversion from plugins/trio-dnm3.c. See the libm-tolerance note in
// native_plugin_trio_dnm3_models.go.
package bcftools

import (
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// dnmScoreType mirrors the DNM_* output-tag types in trio-dnm3.c.
type dnmScoreType int

const (
	dnmScoreLog   dnmScoreType = iota // float, log-scaled (default)
	dnmScorePhred                     // int, phred quality (0-255)
	dnmScoreProb                      // float, probability
)

// runFloatConfig bundles the resolved options for the float-model path, passed
// from RunFull into runFloatModels.
type runFloatConfig struct {
	useModel      int
	scoreType     dnmScoreType
	scoreTag      string
	alleleTag     string
	vafTag        string
	strictlyNovel bool
	useDNGPriors  bool
	withPPL       bool
	withPAD       bool
	withCAD       bool
	forceAD       bool
	hasAD         bool
	hasSP         bool
	needPL        bool
	needQS        bool
	mrate         float64
	phi           float64
	noisePrior    float64
	minVAF        float64
	allelicDrop   float64
	sbCoeff       float64
	minScore      float64
	minQM         float64
	pnSNV         pnoise
	pnIndel       pnoise
	outputType    OutputFormat
}

// runFloatModels executes the DMM/ALM/DNG float models: builds the full pprob
// priors, runs the per-record float pipeline, and writes the output VCF/BCF.
func (p *trioDNM3Plugin) runFloatModels(opts PluginOptions, out io.Writer, hdr *vcf.Header, variants []*vcf.Variant, trios []dnmTrio, filter *pluginFilter, chrX *chrXMatcher, cfg runFloatConfig) error {
	pa := newDNMPriorsFull(cfg.strictlyNovel, cfg.useDNGPriors, cfg.mrate, autosomalPriors)
	pX := newDNMPriorsFull(cfg.strictlyNovel, cfg.useDNGPriors, cfg.mrate, chrXPriors)
	pXX := newDNMPriorsFull(cfg.strictlyNovel, cfg.useDNGPriors, cfg.mrate, chrXXPriors)

	st := &dnmFloatState{
		hdr: hdr, trios: trios, chrX: chrX, priors: pa, priorsX: pX, priorsXX: pXX,
		useModel: cfg.useModel, scoreType: cfg.scoreType, scoreTag: cfg.scoreTag,
		alleleTag: cfg.alleleTag, vafTag: cfg.vafTag, minScore: cfg.minScore, sbCoeff: cfg.sbCoeff,
		hasAD: cfg.hasAD, hasSP: cfg.hasSP, withPAD: cfg.withPAD, withCAD: cfg.withCAD,
		withPPL: cfg.withPPL, forceAD: cfg.forceAD, needPL: cfg.needPL, needQS: cfg.needQS,
		minQM: cfg.minQM, pnSNV: cfg.pnSNV, pnIndel: cfg.pnIndel, phi: cfg.phi,
		noisePrior: cfg.noisePrior, allelicDrop: cfg.allelicDrop, minVAF: cfg.minVAF,
		mrate: cfg.mrate, strictNovel: cfg.strictlyNovel, filter: filter,
		trioPass: make([]bool, len(trios)),
	}

	results := make([]*vcf.Variant, 0, len(variants))
	for _, v := range variants {
		nv, err := st.processRecordFloat(v)
		if err != nil {
			return err
		}
		results = append(results, nv)
	}

	outHdr := buildDNMFloatHeader(hdr, cfg)

	w, cleanup, err := openOutput(out, ViewOptions{OutputFormat: cfg.outputType, CompressLevel: opts.CompressLevel, Threads: opts.Threads}, outHdr)
	if err != nil {
		return err
	}
	if err := w.WriteHeader(); err != nil {
		cleanup()
		return err
	}
	for _, v := range results {
		if err := w.Write(v); err != nil {
			cleanup()
			return err
		}
	}
	if err := w.Flush(); err != nil {
		cleanup()
		return err
	}
	cleanup()
	return nil
}

// buildDNMFloatHeader appends the DNM/VA(/VAF) FORMAT lines to a copy of the
// input header, mirroring the bcf_hdr_printf calls in init_data() for the float
// score types.
func buildDNMFloatHeader(hdr *vcf.Header, cfg runFloatConfig) *vcf.Header {
	outHdr := &vcf.Header{Samples: hdr.Samples}
	outHdr.MetaInfo = append(outHdr.MetaInfo, hdr.MetaInfo...)
	var typeName, vcfType string
	switch cfg.scoreType {
	case dnmScoreLog:
		typeName = "log scaled value (bigger value = bigger confidence)"
		vcfType = "Float"
	case dnmScorePhred:
		typeName = "phred value (bigger value = bigger confidence)"
		vcfType = "Integer"
	case dnmScoreProb:
		typeName = "probability"
		vcfType = "Float"
	}
	outHdr.MetaInfo = append(outHdr.MetaInfo,
		fmt.Sprintf(`##FORMAT=<ID=%s,Number=1,Type=%s,Description="De-novo mutation score given as %s">`, cfg.scoreTag, vcfType, typeName))
	outHdr.MetaInfo = append(outHdr.MetaInfo,
		fmt.Sprintf(`##FORMAT=<ID=%s,Number=1,Type=Integer,Description="The de-novo allele">`, cfg.alleleTag))
	if cfg.hasAD {
		outHdr.MetaInfo = append(outHdr.MetaInfo,
			fmt.Sprintf(`##FORMAT=<ID=%s,Number=1,Type=Integer,Description="The percentage of ALT reads">`, cfg.vafTag))
	}
	return outHdr
}

// parseDNMFloat parses a numeric option value, returning a clear error tied to
// the flag name on failure.
func parseDNMFloat(flag, v string) (float64, error) {
	f, err := strconv.ParseFloat(v, 64)
	if err != nil {
		return 0, fmt.Errorf("trio-dnm3: could not parse: %s %s", flag, v)
	}
	return f, nil
}

// dnmFloatState carries the per-record float-model transform state.
type dnmFloatState struct {
	hdr               *vcf.Header
	trios             []dnmTrio
	chrX              *chrXMatcher
	priors            *dnmPriors
	priorsX, priorsXX *dnmPriors
	useModel          int
	scoreType         dnmScoreType
	scoreTag          string
	alleleTag         string
	vafTag            string
	minScore          float64
	sbCoeff           float64
	hasAD             bool
	hasSP             bool
	withPAD           bool
	withCAD           bool
	withPPL           bool
	forceAD           bool
	needPL            bool
	needQS            bool
	minQM             float64 // signed
	pnSNV, pnIndel    pnoise
	phi               float64
	noisePrior        float64
	allelicDrop       float64
	minVAF            float64
	mrate             float64
	strictNovel       bool
	filter            *pluginFilter
	trioPass          []bool
}

// testFilters applies the -i/-e per-trio filter for the float models; see
// dnmRunFilters.
func (st *dnmFloatState) testFilters(v *vcf.Variant) bool {
	return dnmRunFilters(st.filter, st.trios, st.trioPass, v)
}

// processRecordFloat annotates v with FORMAT/DNM, FORMAT/VA, FORMAT/VAF using
// the selected float model, mirroring process_record() (non-NAIVE branch). It
// returns the record (mutated in place) and an error for fatal data problems
// that upstream reports as error().
func (st *dnmFloatState) processRecordFloat(v *vcf.Variant) (*vcf.Variant, error) {
	nAllele := len(v.Alt) + 1
	if nAllele == 1 || (len(v.Alt) == 1 && (v.Alt[0] == "." || v.Alt[0] == "")) {
		return v, nil
	}
	if !st.testFilters(v) {
		return v, nil
	}

	nSmpl := len(v.Samples)

	// FORMAT/SP (strand bias), one value per sample.
	var sp []int
	nSP := st.hasSP
	if nSP && st.sbCoeff != 0 {
		sp = st.readIntField(v, "SP", nSmpl)
		if sp == nil {
			nSP = false
		}
	} else {
		nSP = false
	}

	// FORMAT/AD: n_allele per sample.
	nAD := 0
	var ad []int
	if st.hasAD {
		ad = st.readIntField(v, "AD", nSmpl*nAllele)
		if ad != nil {
			nAD = nAllele
		} else if st.forceAD {
			// best-effort: try a per-sample width inferred from the first sample
			ad = st.readIntField(v, "AD", -1)
			if ad != nil && len(ad)%nSmpl == 0 {
				nAD = len(ad) / nSmpl
			}
		}
	}
	if st.useModel == dnmUseDMM && nAD == 0 {
		return v, errDNM("the FMT/AD tag is not available at %s:%d", v.Chrom, v.Pos)
	}

	// FORMAT/PL: n_allele*(n_allele+1)/2 per sample.
	npl1 := 0
	var pl []int
	pl = st.readIntField(v, "PL", nSmpl*nAllele*(nAllele+1)/2)
	if pl != nil {
		npl1 = nAllele * (nAllele + 1) / 2
	} else if st.needPL {
		return v, errDNM("the FORMAT/PL tag not present at %s:%d", v.Chrom, v.Pos)
	}

	// FORMAT/QS (ALM, DMM, and any model when trimming >4 alleles).
	nqs1 := 0
	var qs []int
	if st.useModel == dnmUseALM || st.useModel == dnmUseDMM || nAllele > 4 {
		qs = st.readIntField(v, "QS", nSmpl*nAllele)
		if qs == nil {
			if st.needQS {
				return v, errDNM("the FMT/QS tag is not available at %s:%d", v.Chrom, v.Pos)
			}
			if nAD == 0 {
				return v, nil // cannot trim; upstream warns once and skips
			}
			// fake QS from AD assuming BQ=30 (used by --with-pAD).
			qs = make([]int, nSmpl*nAD)
			for i := range qs {
				qs[i] = 30 * ad[i]
			}
			nqs1 = nAD
		} else {
			nqs1 = nAllele
		}
	}

	// FORMAT/QM (DMM only).
	nqm1 := 0
	var qm []int
	if st.useModel == dnmUseDMM {
		if st.minQM > 0 {
			qm = st.readIntField(v, "QM", nSmpl*nAllele)
			if qm == nil {
				return v, errDNM("the FMT/QM tag is not available at %s:%d, run with negative value of --max-QM", v.Chrom, v.Pos)
			}
			nqm1 = nAllele
		} else {
			nqm1 = nAD
		}
	}

	isIndel := variantTypeMask(v)&vtINDEL != 0
	pnCur := st.pnSNV
	if isIndel {
		pnCur = st.pnIndel
	}

	isChrX := st.chrX.overlaps(v.Chrom, v.Pos)

	var outs []dnmOutSlice
	writeDNM := false

	for ti, tr := range st.trios {
		if st.filter != nil && !st.trioPass[ti] {
			continue
		}
		priors := st.priors
		if isChrX {
			if tr.isMale {
				priors = st.priorsX
			} else {
				priors = st.priorsXX
			}
		}

		// Build per-member arrays for [father, mother, child].
		var plM [3][]float64
		if npl1 > 0 {
			plM = setTrioPL(pl, tr, npl1)
		}
		var qsM [3][]float64
		var qmM [3][]float64
		if st.useModel == dnmUseDMM {
			qmM = setTrioQM(qm, ad, tr, nqm1, st.minQM)
			qsM = setTrioQS(qs, tr, nqs1)
		} else if st.useModel == dnmUseALM {
			qsM = setTrioQSNoisy(qs, ad, tr, nqs1, nAD, pnCur)
		} else if nAllele > 4 {
			qsM = setTrioQS(qs, tr, nqs1)
		}
		var adM [3][]int
		if nAD > 0 {
			adM = setTrioAD(ad, tr, nAD, nAllele)
		}

		nals := nAllele
		curNpl := npl1
		var altIdx []int
		if nAllele > 4 {
			altIdx = manyAltsTrim(&nals, &curNpl, plM[:], qsM[:], qmM[:], adM[:], nAllele)
		}

		mp := &dnmModelParams{
			phi: st.phi, minQM: st.minQM, minVAF: st.minVAF, noisePrior: st.noisePrior,
			allelicDrop: st.allelicDrop, withCAD: st.withCAD, withPPL: st.withPPL,
			strictNovel: st.strictNovel, pnCur: pnCur,
		}

		var score float64
		var al0, al1 int
		switch st.useModel {
		case dnmUseDMM:
			score, al0, al1 = processTrioDMM(mp, priors, nals, plM, adM, qmM)
		case dnmUseALM:
			score, al0, al1 = processTrioALM(mp, priors, nals, plM, adM, qsM, priors == st.priorsX)
		case dnmUseDNG:
			score, al0, al1 = processTrioDNG(priors, nals, plM)
		}
		_ = al0

		if nSP {
			s := sp[tr.child]
			if s >= 0 {
				score += st.sbCoeff * phred2log(float64(s))
			}
		}

		if nAllele > 4 {
			al0 = altIdx[al0]
			al1 = altIdx[al1]
		}

		var out dnmOutSlice
		out.sample = tr.child
		out.allele = al1
		if score >= st.minScore {
			writeDNM = true
			switch st.scoreType {
			case dnmScoreLog:
				out.isFloat = true
				out.floatVal = score
			case dnmScoreProb:
				out.isFloat = true
				out.floatVal = math.Exp(score)
			case dnmScorePhred:
				ph := log2phred(subtractLog(0, score))
				if ph > 255 {
					ph = 255
				}
				out.intVal = int(math.Round(ph))
			}
			out.writeDNMTag = true
		}

		if nAD > 0 {
			if al1 < nAD {
				out.writeVAF = true
				idxs := [3]int{tr.father, tr.mother, tr.child}
				for j := 0; j < 3; j++ {
					srcBase := idxs[j] * nAD
					adSum := 0
					for k := 0; k < nAllele; k++ {
						adSum += ad[srcBase+k]
					}
					if adSum != 0 {
						out.vaf[j] = int(math.Round(float64(ad[srcBase+al1]) * 100. / float64(adSum)))
					} else {
						out.vaf[j] = 0
					}
					out.vafSmpl[j] = idxs[j]
				}
			} else {
				out.writeMiss = true
				out.vafSmpl = [3]int{tr.father, tr.mother, tr.child}
			}
		}
		outs = append(outs, out)
	}

	if !writeDNM {
		return v, nil
	}

	st.applyOutputs(v, outs)
	return v, nil
}

// dnmOutSlice is one trio child's computed annotations for a record.
type dnmOutSlice struct {
	sample      int
	writeDNMTag bool    // the score passed --min-score, so DNM/VA are written
	isFloat     bool    // DNM written as a float (log/prob) vs phred int
	floatVal    float64 // DNM value when isFloat
	intVal      int     // DNM value when !isFloat (phred)
	allele      int     // VA, the de-novo allele (al1)
	vaf         [3]int  // [father,mother,child] VAF percentages
	vafSmpl     [3]int  // sample indices for the vaf entries
	writeVAF    bool    // VAF computed (al1 within AD range)
	writeMiss   bool    // VAF explicitly missing (al1 out of AD range)
}

// applyOutputs writes the DNM/VA/VAF FORMAT fields onto v from the per-trio
// results, mirroring the bcf_update_format_* calls. Samples not covered by a
// trio remain "." for the per-trio tags; VAF is filled for trio members only.
func (st *dnmFloatState) applyOutputs(v *vcf.Variant, outs []dnmOutSlice) {
	scoreVals := map[int]string{}
	alleleVals := map[int]string{}
	vafVals := map[int]string{}
	anyVAF := false
	for _, o := range outs {
		if o.writeDNMTag {
			if o.isFloat {
				scoreVals[o.sample] = formatVCFFloat(o.floatVal)
			} else {
				scoreVals[o.sample] = strconv.Itoa(o.intVal)
			}
			alleleVals[o.sample] = strconv.Itoa(o.allele)
		}
		if o.writeVAF {
			anyVAF = true
			for j := 0; j < 3; j++ {
				vafVals[o.vafSmpl[j]] = strconv.Itoa(o.vaf[j])
			}
		} else if o.writeMiss {
			for j := 0; j < 3; j++ {
				if _, ok := vafVals[o.vafSmpl[j]]; !ok {
					vafVals[o.vafSmpl[j]] = "."
				}
			}
		}
	}

	v.Format = append(v.Format, st.scoreTag, st.alleleTag)
	if anyVAF {
		v.Format = append(v.Format, st.vafTag)
	}
	for s := range v.Samples {
		if sv, ok := scoreVals[s]; ok {
			v.Samples[s].Data[st.scoreTag] = sv
			v.Samples[s].Data[st.alleleTag] = alleleVals[s]
		}
		if anyVAF {
			if vv, ok := vafVals[s]; ok {
				v.Samples[s].Data[st.vafTag] = vv
			}
		}
	}
}

// readIntField parses FORMAT/tag as a flat per-sample integer slice. When
// wantTotal >= 0 it requires exactly that many values (nSmpl * per-sample
// width); a mismatch returns nil so the caller can fall back. When wantTotal < 0
// it returns whatever it can parse (used by --force-AD). Missing/"." entries
// become 0 (matching htslib's bcf_get_format_int32 0-fill for present samples;
// a wholly-absent tag returns nil).
func (st *dnmFloatState) readIntField(v *vcf.Variant, tag string, wantTotal int) []int {
	nSmpl := len(v.Samples)
	// Determine per-sample width from the first present sample.
	width := -1
	anyPresent := false
	for i := 0; i < nSmpl; i++ {
		s, ok := v.Samples[i].Data[tag]
		if !ok || s == "" || s == "." {
			continue
		}
		anyPresent = true
		w := strings.Count(s, ",") + 1
		if w > width {
			width = w
		}
	}
	if !anyPresent || width <= 0 {
		return nil
	}
	out := make([]int, nSmpl*width)
	for i := 0; i < nSmpl; i++ {
		base := i * width
		s, ok := v.Samples[i].Data[tag]
		if !ok || s == "" || s == "." {
			continue // 0-filled
		}
		parts := strings.Split(s, ",")
		for k := 0; k < width && k < len(parts); k++ {
			n, err := strconv.Atoi(strings.TrimSpace(parts[k]))
			if err == nil {
				out[base+k] = n
			}
		}
	}
	if wantTotal >= 0 && len(out) != wantTotal {
		return nil
	}
	return out
}

// errDNM formats a trio-dnm3 fatal-data error.
func errDNM(format string, a ...interface{}) error {
	return fmt.Errorf("trio-dnm3: "+format, a...)
}

// setTrioPL builds the per-member normalised log-prob PL arrays, mirroring
// set_trio_PL(): each member's PLs are converted to linear probs, summed, and
// re-logged as log(prob/sum).
func setTrioPL(pl []int, tr dnmTrio, npl1 int) [3][]float64 {
	var out [3][]float64
	idxs := [3]int{tr.father, tr.mother, tr.child}
	for j := 0; j < 3; j++ {
		src := idxs[j] * npl1
		dst := make([]float64, npl1)
		sum := 0.0
		for k := 0; k < npl1; k++ {
			dst[k] = phred2num(float64(pl[src+k]))
			sum += dst[k]
		}
		for k := 0; k < npl1; k++ {
			dst[k] = math.Log(dst[k] / sum)
		}
		out[j] = dst
	}
	return out
}

// setTrioAD builds the per-member AD arrays, mirroring set_trio_AD().
func setTrioAD(ad []int, tr dnmTrio, nad1, nals int) [3][]int {
	var out [3][]int
	idxs := [3]int{tr.father, tr.mother, tr.child}
	for j := 0; j < 3; j++ {
		src := idxs[j] * nad1
		dst := make([]int, nals)
		for k := 0; k < nals && k < nad1; k++ {
			dst[k] = ad[src+k]
		}
		out[j] = dst
	}
	return out
}

// setTrioQM builds the per-member QM (per-read error) arrays, mirroring
// set_trio_QM(): phred2num(QM) floored at fabs(min_qm), or fabs(min_qm) when QM
// is absent.
func setTrioQM(qm, ad []int, tr dnmTrio, nqm1 int, minQM float64) [3][]float64 {
	var out [3][]float64
	minq := math.Abs(minQM)
	idxs := [3]int{tr.father, tr.mother, tr.child}
	for j := 0; j < 3; j++ {
		dst := make([]float64, nqm1)
		if qm == nil {
			for k := 0; k < nqm1; k++ {
				dst[k] = minq
			}
		} else {
			src := idxs[j] * nqm1
			for k := 0; k < nqm1; k++ {
				q := qm[src+k]
				if q != 0 {
					dst[k] = phred2num(float64(q))
				} else {
					dst[k] = minq
				}
				if dst[k] < minq {
					dst[k] = minq
				}
			}
		}
		out[j] = dst
	}
	return out
}

// setTrioQS builds the per-member QS log-prob arrays, mirroring set_trio_QS():
// phred2log(QS).
func setTrioQS(qs []int, tr dnmTrio, nqs1 int) [3][]float64 {
	var out [3][]float64
	idxs := [3]int{tr.father, tr.mother, tr.child}
	for j := 0; j < 3; j++ {
		src := idxs[j] * nqs1
		dst := make([]float64, nqs1)
		for k := 0; k < nqs1; k++ {
			dst[k] = phred2log(float64(qs[src+k]))
		}
		out[j] = dst
	}
	return out
}

// setTrioQSNoisy builds the noise-adjusted per-member QS arrays for ALM,
// mirroring set_trio_QS_noisy(): parental QS are reduced by a per-read noise
// allowance before phred2log.
func setTrioQSNoisy(qs, ad []int, tr dnmTrio, nqs1, nAD int, pn pnoise) [3][]float64 {
	var out [3][]float64
	useAD := nAD
	if useAD != 0 && pn.abs == 0 && pn.abs1 == 0 && pn.frac1 == 0 {
		useAD = 0
	}
	idxs := [3]int{tr.father, tr.mother, tr.child}
	var adF, adM []int
	if useAD != 0 {
		adF = ad[tr.father*nAD : tr.father*nAD+nAD]
		adM = ad[tr.mother*nAD : tr.mother*nAD+nAD]
	}
	for j := 0; j < 3; j++ {
		qsrc := idxs[j] * nqs1
		dst := make([]float64, nqs1)
		var pnv, pns float64
		if (pn.frac != 0 || pn.frac1 != 0) && j != iCHILD {
			sumQS := 0.0
			for k := 0; k < nqs1; k++ {
				sumQS += float64(qs[qsrc+k])
			}
			pnv = sumQS * pn.frac
			pns = sumQS * pn.frac1
			if useAD != 0 {
				asrc := idxs[j] * nAD
				sumAD := 0.0
				for k := 0; k < nAD; k++ {
					sumAD += float64(ad[asrc+k])
				}
				if sumAD != 0 {
					if pnv < pn.abs*sumQS/sumAD {
						pnv = pn.abs * sumQS / sumAD
					}
					if pns < pn.abs1*sumQS/sumAD {
						pns = pn.abs1 * sumQS / sumAD
					}
				}
			}
		}
		for k := 0; k < nqs1; k++ {
			val := float64(qs[qsrc+k])
			if useAD != 0 && (adF[k] == 0 || adM[k] == 0) {
				val -= pns
			} else {
				val -= pnv
			}
			if val < 0 {
				val = 0
			}
			dst[k] = phred2log(val)
		}
		out[j] = dst
	}
	return out
}

// manyAltsTrim reduces a >4-allele site to the four most-likely alleles (by
// summed QS across the trio), permuting the per-member PL/QS/QM/AD arrays in
// place, mirroring many_alts_trim(). It returns the alt_idx map (compact index
// -> original allele) used to translate the de-novo allele back afterwards.
func manyAltsTrim(nals, npl *int, plM [][]float64, qsM, qmM [][]float64, adM [][]int, origNals int) []int {
	// Sum QS across the three members per allele.
	sumQS := make([]float64, origNals)
	for i := 0; i < 3; i++ {
		for j := 0; j < origNals; j++ {
			sumQS[j] += qsM[i][j]
		}
	}
	idx := make([]int, origNals)
	for i := range idx {
		idx[i] = i
	}
	// Insertion sort keeping REF (index 0) first, ascending by sumQS for 1..n.
	for i := 2; i < origNals; i++ {
		for j := i; j > 1 && sumQS[idx[j]] < sumQS[idx[j-1]]; j-- {
			idx[j], idx[j-1] = idx[j-1], idx[j]
		}
	}

	// QS, QM, AD: upstream memcpy's only the 4 reordered values into the head of
	// each per-member array, leaving the tail (index >=4) untouched at its
	// original value. The DMM AD/QM-consuming helpers iterate over the ORIGINAL
	// allele count (nad), so they read that stale tail; we must preserve it
	// rather than slice to length 4. (ALM/QS only ever indexes i<nals=4, so the
	// tail is never read there.)
	for i := 0; i < 3; i++ {
		head := make([]float64, 4)
		for j := 0; j < 4; j++ {
			head[j] = qsM[i][idx[j]]
		}
		copy(qsM[i], head)
	}
	if qmM[0] != nil {
		for i := 0; i < 3; i++ {
			head := make([]float64, 4)
			for j := 0; j < 4; j++ {
				head[j] = qmM[i][idx[j]]
			}
			copy(qmM[i], head)
		}
	}
	if adM[0] != nil {
		for i := 0; i < 3; i++ {
			head := make([]int, 4)
			for j := 0; j < 4; j++ {
				head[j] = adM[i][idx[j]]
			}
			copy(adM[i], head)
		}
	}
	*npl = 10
	*nals = 4
	if plM[0] != nil {
		for i := 0; i < 3; i++ {
			trimmed := make([]float64, 10)
			for j := 0; j < 4; j++ {
				for k := 0; k <= j; k++ {
					idst := bcfAlleles2GT(j, k)
					isrc := bcfAlleles2GT(idx[j], idx[k])
					trimmed[idst] = plM[i][isrc]
				}
			}
			copy(plM[i], trimmed)
			plM[i] = plM[i][:10]
		}
	} else {
		*npl = 0
	}
	return idx
}

// bcfAlleles2GT returns the genotype index for alleles a,b (a>=b not required),
// mirroring htslib's bcf_alleles2gt: gt = a*(a+1)/2 + b with a>=b.
func bcfAlleles2GT(a, b int) int {
	if a < b {
		a, b = b, a
	}
	return a*(a+1)/2 + b
}
