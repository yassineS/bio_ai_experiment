// Native port of the upstream `fill-tags` plugin (plugins/fill-tags.c) and its
// deprecated subset `fill-AN-AC` (plugins/fill-AN-AC.c).
//
// fill-tags recomputes INFO/FORMAT annotations from FORMAT/GT: AN, AC, AC_Hom,
// AC_Het, AC_Hemi, AF, MAF, NS, HWE, ExcHet, END, TYPE, and FORMAT/VAF, VAF1.
// The per-allele het/hom/hemi/half counting follows process_fmt's BRANCH_INT
// loop exactly, including the --drop-missing (-d) treatment of half-missing
// genotypes. The single default population "ALL" is supported; the -S
// sample-population grouping and the experimental TAG=func(EXPR) expressions
// require the filter engine and are reported as unsupported in batch 1.
package bcftools

import (
	"fmt"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("fill-tags", func() NativePlugin { return &fillTagsPlugin{} })
	registerNativePlugin("fill-AN-AC", func() NativePlugin { return &fillTagsPlugin{anacOnly: true} })
}

// fill-tags tag selection bits, matching the SET_* flags in fill-tags.c.
const (
	setAN = 1 << iota
	setAC
	setACHom
	setACHet
	setACHemi
	setAF
	setNS
	setMAF
	setHWE
	setExcHet
	setEND
	setType
	setVAF
	setVAF1
	setFMissing
)

// fillTagsPlugin implements fill-tags / fill-AN-AC. It is per-record and
// parallel: every record is recomputed independently from its own FORMAT/GT.
type fillTagsPlugin struct {
	anacOnly    bool // true for the fill-AN-AC entry point
	tags        int
	dropMissing bool
}

// alleleCounts accumulates, per allele, the het/hom/hemi/half-missing counts
// used to derive every fill-tags annotation. It mirrors fill-tags.c counts_t.
type alleleCounts struct {
	nhom, nhet, nhemi, nac int
}

// Name returns the plugin name.
func (p *fillTagsPlugin) Name() string {
	if p.anacOnly {
		return "fill-AN-AC"
	}
	return "fill-tags"
}

// About returns the one-line description, matching the upstream about().
func (p *fillTagsPlugin) About() string {
	if p.anacOnly {
		return "Fill INFO fields AN and AC. This plugin is DEPRECATED, use fill-tags instead."
	}
	return "Set INFO tags AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, ExcHet, HWE, MAF, NS; FORMAT/VAF and more."
}

// Parallel reports true: each record is recomputed independently.
func (p *fillTagsPlugin) Parallel() bool { return true }

// Init parses the plugin arguments and appends the relevant ##INFO/##FORMAT
// header lines in upstream's fixed order.
func (p *fillTagsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	tagsStr := "all"
	if p.anacOnly {
		p.tags = setAN | setAC
	} else {
		for i := 0; i < len(args); i++ {
			a := args[i]
			switch a {
			case "-d", "--drop-missing":
				p.dropMissing = true
			case "-t", "--tags":
				if i+1 >= len(args) {
					return nil, fmt.Errorf("fill-tags: -t requires an argument")
				}
				i++
				tagsStr = args[i]
			case "-l", "--list-tags":
				return nil, fmt.Errorf("fill-tags: -l (list-tags) is not supported in the native plugin")
			case "-S", "--samples-file":
				return nil, fmt.Errorf("fill-tags: -S/--samples-file (population grouping) is not supported in the native plugin")
			default:
				if strings.HasPrefix(a, "-t") && len(a) > 2 {
					tagsStr = a[2:]
					continue
				}
				return nil, fmt.Errorf("fill-tags: unsupported option %q", a)
			}
		}
		flag, err := p.parseTags(tagsStr)
		if err != nil {
			return nil, err
		}
		p.tags = flag
	}

	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	add := func(line string) { out.MetaInfo = appendInfoHeader(out.MetaInfo, line) }
	// F_MISSING is added first (upstream adds the expression's header during
	// parse_tags, before the fixed hdr_append block below).
	if p.tags&setFMissing != 0 {
		add(`##INFO=<ID=F_MISSING,Number=1,Type=Float,Description="Added by +fill-tags expression F_MISSING:1=F_MISSING">`)
	}
	if p.tags&setAN != 0 {
		add(`##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">`)
	}
	if p.tags&setAC != 0 {
		add(`##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">`)
	}
	if p.tags&setNS != 0 {
		add(`##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">`)
	}
	if p.tags&setACHom != 0 {
		add(`##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="Allele counts in homozygous genotypes">`)
	}
	if p.tags&setACHet != 0 {
		add(`##INFO=<ID=AC_Het,Number=A,Type=Integer,Description="Allele counts in heterozygous genotypes">`)
	}
	if p.tags&setACHemi != 0 {
		add(`##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description="Allele counts in hemizygous genotypes">`)
	}
	if p.tags&setAF != 0 {
		add(`##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">`)
	}
	if p.tags&setMAF != 0 {
		add(`##INFO=<ID=MAF,Number=1,Type=Float,Description="Frequency of the second most common allele">`)
	}
	if p.tags&setHWE != 0 {
		add(`##INFO=<ID=HWE,Number=A,Type=Float,Description="HWE test (PMID:15789306); 1=good, 0=bad">`)
	}
	if p.tags&setEND != 0 {
		add(`##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">`)
	}
	if p.tags&setType != 0 {
		add(`##INFO=<ID=TYPE,Number=.,Type=String,Description="Variant type">`)
	}
	if p.tags&setExcHet != 0 {
		add(`##INFO=<ID=ExcHet,Number=A,Type=Float,Description="Test excess heterozygosity; 1=good, 0=bad">`)
	}
	if p.tags&setVAF != 0 {
		add(`##FORMAT=<ID=VAF,Number=A,Type=Float,Description="The fraction of reads with alternate allele (nALT/nSumAll)">`)
	}
	if p.tags&setVAF1 != 0 {
		add(`##FORMAT=<ID=VAF1,Number=1,Type=Float,Description="The fraction of reads with alternate alleles (nSumALT/nSumAll)">`)
	}
	return out, nil
}

// parseTags converts the comma-separated -t list into the SET_* bitmask. The
// "all" keyword selects every supported tag except END/TYPE (matching upstream,
// where 'all' excludes SET_END|SET_TYPE), but unlike upstream it does not pull
// in the F_MISSING expression (which needs the filter engine).
func (p *fillTagsPlugin) parseTags(str string) (int, error) {
	flag := 0
	for _, raw := range strings.Split(str, ",") {
		t := strings.TrimSpace(raw)
		t = strings.TrimPrefix(t, "INFO/")
		t = strings.TrimPrefix(t, "FORMAT/")
		switch strings.ToLower(t) {
		case "all":
			flag |= setAN | setAC | setACHom | setACHet | setACHemi | setAF | setNS | setMAF | setHWE | setExcHet | setVAF | setVAF1 | setFMissing
		case "an":
			flag |= setAN
		case "ac":
			flag |= setAC
		case "ns":
			flag |= setNS
		case "ac_hom":
			flag |= setACHom
		case "ac_het":
			flag |= setACHet
		case "ac_hemi":
			flag |= setACHemi
		case "af":
			flag |= setAF
		case "maf":
			flag |= setMAF
		case "hwe":
			flag |= setHWE
		case "exchet":
			flag |= setExcHet
		case "end":
			flag |= setEND
		case "type":
			flag |= setType
		case "vaf":
			flag |= setVAF
		case "vaf1":
			flag |= setVAF1
		case "f_missing":
			flag |= setFMissing
		default:
			if strings.Contains(t, "=") {
				return 0, fmt.Errorf("fill-tags: custom expression %q is not supported in the native plugin", t)
			}
			return 0, fmt.Errorf("fill-tags: the tag %q is not supported", t)
		}
	}
	return flag, nil
}

// Process recomputes the requested annotations for a single record.
func (p *fillTagsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	nals := 1 + len(v.Alt)

	// Per-allele het/hom/hemi/half counts and the number of samples with data.
	counts := make([]alleleCounts, nals)
	ns := 0
	hasGT := formatHasTag(v, "GT")
	if hasGT && (p.tags&(setAN|setAC|setACHom|setACHet|setACHemi|setAF|setMAF|setNS|setHWE|setExcHet) != 0) {
		ns = p.countGenotypes(v, counts, nals)
	}

	// AN: total called alleles across all alleles.
	an := 0
	for j := 0; j < nals; j++ {
		an += counts[j].nhet + counts[j].nhom + counts[j].nhemi + counts[j].nac
	}

	// F_MISSING is filled before the FORMAT-derived tags, matching upstream
	// where the expression functions run before process_fmt.
	if p.tags&setFMissing != 0 {
		setInfo(v, "F_MISSING", formatVCFFloat(fractionMissing(v)))
	}
	if p.tags&setNS != 0 && hasGT {
		setInfo(v, "NS", strconv.Itoa(ns))
	}
	if p.tags&setAN != 0 && hasGT {
		setInfo(v, "AN", strconv.Itoa(an))
	}
	if p.tags&(setAF|setMAF) != 0 && hasGT {
		p.fillAFMAF(v, counts, nals, an)
	}
	if p.tags&setAC != 0 && hasGT {
		p.fillAC(v, counts, nals)
	}
	if p.tags&setACHet != 0 && hasGT {
		p.fillPerAlt(v, "AC_Het", counts, nals, func(c alleleCounts) int { return c.nhet })
	}
	if p.tags&setACHom != 0 && hasGT {
		p.fillPerAlt(v, "AC_Hom", counts, nals, func(c alleleCounts) int { return c.nhom })
	}
	if p.tags&setACHemi != 0 && nals > 1 && hasGT {
		p.fillPerAlt(v, "AC_Hemi", counts, nals, func(c alleleCounts) int { return c.nhemi })
	}
	if p.tags&(setHWE|setExcHet) != 0 && hasGT {
		p.fillHWE(v, counts, nals)
	}
	if p.tags&(setVAF|setVAF1) != 0 {
		p.fillVAF(v, nals)
	}
	if p.tags&setEND != 0 {
		setInfo(v, "END", strconv.Itoa(v.Pos+refLen(v)-1))
	}
	if p.tags&setType != 0 {
		setInfo(v, "TYPE", typeMacroValue(v))
	}
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *fillTagsPlugin) Destroy() error { return nil }

// countGenotypes fills the per-allele counts from FORMAT/GT, mirroring the
// BRANCH_INT loop of process_fmt: it classifies each sample's genotype as
// het/hom/hemi/half and increments the relevant per-allele counter. It returns
// the number of samples with at least one called allele (NS).
func (p *fillTagsPlugin) countGenotypes(v *vcf.Variant, counts []alleleCounts, nals int) int {
	ns := 0
	for i := range v.Samples {
		gt, ok := sampleGT(v, i)
		if !ok {
			continue
		}
		// Collect distinct present alleles and total called allele count.
		present := make([]bool, nals)
		nbits := 0
		ncalled := 0
		for _, a := range gt.alleles {
			if a == missingAllele {
				continue
			}
			if a < 0 || a >= nals {
				continue
			}
			ncalled++
			if !present[a] {
				nbits++
				present[a] = true
			}
		}
		if ncalled == 0 {
			continue // fully missing genotype
		}
		isHom := nbits == 1
		var isHemi, isHalf bool
		switch {
		case ncalled != gt.ploidy():
			// some alleles missing within the genotype
			if p.dropMissing {
				isHemi, isHalf = false, true
			} else {
				isHemi, isHalf = true, false
			}
		case ncalled == 1:
			isHemi, isHalf = true, false
		default:
			isHemi, isHalf = false, false
		}
		for a := 0; a < nals; a++ {
			if !present[a] {
				continue
			}
			switch {
			case isHalf:
				counts[a].nac++
			case !isHom:
				counts[a].nhet++
			case !isHemi:
				counts[a].nhom += 2
			default:
				counts[a].nhemi++
			}
		}
		ns++
	}
	return ns
}

// alleleTotal is the total called count of allele j (het+hom+hemi+half).
func alleleTotal(c alleleCounts) int { return c.nhet + c.nhom + c.nhemi + c.nac }

func (p *fillTagsPlugin) fillAC(v *vcf.Variant, counts []alleleCounts, nals int) {
	parts := make([]string, 0, nals-1)
	for j := 1; j < nals; j++ {
		parts = append(parts, strconv.Itoa(alleleTotal(counts[j])))
	}
	setInfo(v, "AC", strings.Join(parts, ","))
}

func (p *fillTagsPlugin) fillPerAlt(v *vcf.Variant, key string, counts []alleleCounts, nals int, sel func(alleleCounts) int) {
	parts := make([]string, 0, nals-1)
	for j := 1; j < nals; j++ {
		parts = append(parts, strconv.Itoa(sel(counts[j])))
	}
	setInfo(v, key, strings.Join(parts, ","))
}

// fillAFMAF computes per-allele frequencies. AF is the per-ALT frequency;
// MAF is the frequency of the second most common allele over all alleles.
func (p *fillTagsPlugin) fillAFMAF(v *vcf.Variant, counts []alleleCounts, nals int, an int) {
	freq := make([]float64, nals)
	if nals > 1 {
		for j := 0; j < nals; j++ {
			freq[j] = float64(alleleTotal(counts[j]))
		}
		if an != 0 {
			for j := 0; j < nals; j++ {
				freq[j] /= float64(an)
			}
		}
	}
	if p.tags&setAF != 0 {
		parts := make([]string, 0, nals-1)
		for j := 1; j < nals; j++ {
			if nals > 1 && an == 0 {
				parts = append(parts, ".")
			} else {
				parts = append(parts, formatVCFFloat(freq[j]))
			}
		}
		setInfo(v, "AF", strings.Join(parts, ","))
	}
	if nals > 1 && p.tags&setMAF != 0 {
		sorted := append([]float64(nil), freq...)
		if an != 0 {
			sort.Sort(sort.Reverse(sort.Float64Slice(sorted)))
		}
		// MAF = second most common allele frequency = sorted[1].
		maf := sorted[1]
		if an == 0 {
			setInfo(v, "MAF", ".")
		} else {
			setInfo(v, "MAF", formatVCFFloat(maf))
		}
	}
}

// fillHWE computes the HWE and ExcHet p-values per ALT allele, porting calc_hwe
// (Wigginton 2005, PMID 15789306).
func (p *fillTagsPlugin) fillHWE(v *vcf.Variant, counts []alleleCounts, nals int) {
	hwe := make([]string, 0, nals-1)
	exc := make([]string, 0, nals-1)
	if nals > 1 {
		nrefTot := counts[0].nhom
		for j := 0; j < nals; j++ {
			nrefTot += counts[j].nhet
		}
		for j := 1; j < nals; j++ {
			nref := nrefTot - counts[j].nhet
			nalt := counts[j].nhet + counts[j].nhom
			nhet := counts[j].nhet
			var ph, pe float64 = 1, 1
			if nref > 0 && nalt > 0 {
				ph, pe = calcHWE(nref, nalt, nhet)
			}
			hwe = append(hwe, formatVCFFloat(ph))
			exc = append(exc, formatVCFFloat(pe))
		}
	}
	if p.tags&setHWE != 0 {
		setInfo(v, "HWE", strings.Join(hwe, ","))
	}
	if p.tags&setExcHet != 0 {
		setInfo(v, "ExcHet", strings.Join(exc, ","))
	}
}

// fillVAF computes FORMAT/VAF and VAF1 from FORMAT/AD, porting process_vaf.
func (p *fillTagsPlugin) fillVAF(v *vcf.Variant, nals int) {
	if nals <= 1 {
		return
	}
	if !formatHasTag(v, "AD") {
		return
	}
	doVAF := p.tags&setVAF != 0
	doVAF1 := p.tags&setVAF1 != 0
	var vafVals, vaf1Vals []string
	for i := range v.Samples {
		ad, ok := v.Samples[i].Data["AD"]
		adVals, parseOK := parseADList(ad)
		valid := ok && parseOK && len(adVals) == nals
		sum := 0
		if valid {
			for _, x := range adVals {
				sum += x
			}
		}
		if doVAF {
			if !valid {
				// Upstream sets dst[0]=missing and the rest to vector-end, which
				// renders as a single "." in VCF text.
				vafVals = append(vafVals, ".")
			} else {
				per := make([]string, nals-1)
				for j := 1; j < nals; j++ {
					if sum != 0 {
						per[j-1] = formatVCFFloat(float64(adVals[j]) / float64(sum))
					} else {
						per[j-1] = "0"
					}
				}
				vafVals = append(vafVals, strings.Join(per, ","))
			}
		}
		if doVAF1 {
			if !valid {
				vaf1Vals = append(vaf1Vals, ".")
			} else if sum != 0 {
				vaf1Vals = append(vaf1Vals, formatVCFFloat(float64(sum-adVals[0])/float64(sum)))
			} else {
				vaf1Vals = append(vaf1Vals, "0")
			}
		}
	}
	if doVAF {
		ensureFormatTag(v, "VAF")
		for i := range v.Samples {
			v.Samples[i].Data["VAF"] = vafVals[i]
		}
	}
	if doVAF1 {
		ensureFormatTag(v, "VAF1")
		for i := range v.Samples {
			v.Samples[i].Data["VAF1"] = vaf1Vals[i]
		}
	}
}

// calcHWE ports fill-tags.c calc_hwe: it returns (p_hwe, p_exc_het) for a
// biallelic configuration with nref reference alleles, nalt alt alleles, and
// nhet heterozygous genotypes (assuming diploid: total genotypes (nref+nalt)/2).
func calcHWE(nref, nalt, nhet int) (float64, float64) {
	ngt := (nref + nalt) / 2
	nrare := nref
	if nalt < nref {
		nrare = nalt
	}
	probs := make([]float64, nrare+1)

	mid := int(float64(nrare) * float64(nref+nalt-nrare) / float64(nref+nalt))
	if (nrare&1)^(mid&1) != 0 {
		mid++
	}
	het := mid
	homR := (nrare - mid) / 2
	homC := ngt - het - homR
	sum := 1.0
	probs[mid] = 1.0

	for het = mid; het > 1; het -= 2 {
		probs[het-2] = probs[het] * float64(het) * float64(het-1) / (4.0 * float64(homR+1) * float64(homC+1))
		sum += probs[het-2]
		homR++
		homC++
	}

	het = mid
	homR = (nrare - mid) / 2
	homC = ngt - het - homR
	for het = mid; het <= nrare-2; het += 2 {
		probs[het+2] = probs[het] * 4.0 * float64(homR) * float64(homC) / (float64(het+2) * float64(het+1))
		sum += probs[het+2]
		homR--
		homC--
	}

	for h := 0; h < nrare+1; h++ {
		probs[h] /= sum
	}

	pExcHet := probs[nhet]
	for h := nhet + 1; h <= nrare; h++ {
		pExcHet += probs[h]
	}

	pHWE := 0.0
	for h := 0; h <= nrare; h++ {
		if probs[h] > probs[nhet] {
			continue
		}
		pHWE += probs[h]
	}
	if pHWE > 1 {
		pHWE = 1
	}
	return pHWE, pExcHet
}

// fractionMissing returns the fraction of samples whose genotype has at least
// one missing allele, matching the upstream F_MISSING expression
// F_PASS(GT="mis"): in bcftools the "mis" genotype class matches a sample with
// any missing allele (e.g. both "./." and the partial "./1"). A sample with no
// GT field counts as missing; with zero samples it is 0.
func fractionMissing(v *vcf.Variant) float64 {
	if len(v.Samples) == 0 {
		return 0
	}
	nmiss := 0
	for i := range v.Samples {
		gt, ok := sampleGT(v, i)
		if !ok || gt.nMissing() > 0 {
			nmiss++
		}
	}
	return float64(nmiss) / float64(len(v.Samples))
}

// formatHasTag reports whether tag appears in v.Format.
func formatHasTag(v *vcf.Variant, tag string) bool {
	for _, f := range v.Format {
		if f == tag {
			return true
		}
	}
	return false
}

// parseADList parses a comma-separated FORMAT/AD integer list. A "." entry (or
// an empty/missing field) marks the value as not fully present (false), to
// match upstream's "all values present or skip the sample" handling for AD.
func parseADList(s string) ([]int, bool) {
	if s == "" || s == "." {
		return nil, false
	}
	parts := strings.Split(s, ",")
	out := make([]int, len(parts))
	for i, p := range parts {
		if p == "." {
			return nil, false
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			return nil, false
		}
		out[i] = n
	}
	return out, true
}

// refLen returns the reference-span length of the record (len(REF) for simple
// records). This is used to compute INFO/END = POS + rlen - 1.
func refLen(v *vcf.Variant) int {
	if end, ok := v.Info["END"]; ok {
		if n, err := strconv.Atoi(end); err == nil {
			return n - v.Pos + 1
		}
	}
	return len(v.Ref)
}

// appendInfoHeader inserts an ##INFO/##FORMAT line into the meta lines if a
// definition for the same ID is not already present, preserving upstream's
// behaviour of not duplicating an existing header definition. New lines are
// appended after the last existing ##INFO/##FORMAT line group, mirroring how
// bcf_hdr_append grows the header.
func appendInfoHeader(meta []string, line string) []string {
	id := headerID(line)
	if id != "" {
		for _, m := range meta {
			if headerID(m) == id && headerKind(m) == headerKind(line) {
				return meta // already defined; do not duplicate
			}
		}
	}
	return append(meta, line)
}

// headerID extracts the ID=... value from a structured header line.
func headerID(line string) string {
	i := strings.Index(line, "ID=")
	if i < 0 {
		return ""
	}
	rest := line[i+3:]
	end := strings.IndexAny(rest, ",>")
	if end < 0 {
		return rest
	}
	return rest[:end]
}

// headerKind returns the header category prefix (e.g. "##INFO", "##FORMAT").
func headerKind(line string) string {
	if i := strings.Index(line, "=<"); i >= 0 {
		return line[:i]
	}
	return line
}
