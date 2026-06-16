// Text-model port of htslib's bcf_remove_allele_set (vcfutils.c), used by the
// ad-bias --clean-vcf mode to drop ALT alleles that fail the Fisher threshold.
// It rewrites the ALT list and remaps the Number=A / Number=R / Number=G INFO
// and FORMAT fields (and reindexes FORMAT/GT) for the retained alleles, exactly
// as htslib does for the diploid/haploid common case. REF (allele 0) is always
// kept.
package bcftools

import (
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// removeAlleleSet removes the ALT alleles flagged in rm (indexed by allele
// number, 0==REF) from v in place. rm[0] must be false (the reference allele
// cannot be removed). When no allele is removed it is a no-op. It mirrors
// bcf_remove_allele_set for the Number=A/R/G INFO and FORMAT fields plus the
// FORMAT/GT reindexing; symbolic/structural-variant special cases are not
// needed by ad-bias and are left untouched.
func removeAlleleSet(hdr *vcf.Header, v *vcf.Variant, rm []bool) error {
	nROri := 1 + len(v.Alt)
	if len(rm) < nROri {
		grown := make([]bool, nROri)
		copy(grown, rm)
		rm = grown
	}

	// Build the old->new allele-index map and the new ALT list.
	mapIdx := make([]int, nROri)
	mapIdx[0] = 0
	newAlt := make([]string, 0, len(v.Alt))
	nrm := 0
	j := 1
	for i := 1; i < nROri; i++ {
		if rm[i] {
			mapIdx[i] = -1
			nrm++
			continue
		}
		mapIdx[i] = j
		j++
		newAlt = append(newAlt, v.Alt[i-1])
	}
	if nrm == 0 {
		return nil
	}

	// Remap INFO Number=A/R/G fields, then the ALT list.
	numbers := headerNumberMapsFrom(hdr.MetaInfo)
	for _, key := range v.InfoOrder {
		val, ok := v.Info[key]
		if !ok {
			continue
		}
		num := numbers.info[key]
		newVal, changed := subsetNumberedList(val, num, rm, mapIdx, nROri)
		if changed {
			v.Info[key] = newVal
		}
	}

	v.Alt = newAlt

	// Reindex FORMAT/GT (allele indices change; removed alleles become missing,
	// preserving phasing).
	if formatHasTag(v, "GT") {
		needRemap := false
		for i := 1; i < nROri; i++ {
			if mapIdx[i] != i {
				needRemap = true
				break
			}
		}
		if needRemap {
			for s := range v.Samples {
				gt, ok := parseGT(v.Samples[s].Data["GT"])
				if !ok {
					continue
				}
				for k, al := range gt.alleles {
					if al == missingAllele || al < 0 || al >= nROri {
						continue
					}
					if mapIdx[al] < 0 {
						gt.alleles[k] = missingAllele
					} else {
						gt.alleles[k] = mapIdx[al]
					}
				}
				v.Samples[s].Data["GT"] = gt.String()
			}
		}
	}

	// Remap FORMAT Number=A/R/G fields (skip GT, handled above).
	for _, tag := range v.Format {
		if tag == "GT" {
			continue
		}
		num := numbers.format[tag]
		if num != "A" && num != "R" && num != "G" {
			continue
		}
		for s := range v.Samples {
			val, ok := v.Samples[s].Data[tag]
			if !ok {
				continue
			}
			newVal, changed := subsetNumberedList(val, num, rm, mapIdx, nROri)
			if changed {
				v.Samples[s].Data[tag] = newVal
			}
		}
	}

	return nil
}

// subsetNumberedList rewrites a comma-separated A/R/G value list for the
// retained alleles, mirroring the Number=A/R/G handling in
// bcf_remove_allele_set. A non-A/R/G number, a missing ("." / "") value, or a
// value whose element count does not match the expected cardinality is returned
// unchanged (changed=false), matching upstream's "leave it alone on a count
// mismatch / missing" behaviour.
func subsetNumberedList(val, number string, rm []bool, mapIdx []int, nROri int) (string, bool) {
	if number != "A" && number != "R" && number != "G" {
		return val, false
	}
	if val == "" || val == "." {
		return val, false
	}
	parts := strings.Split(val, ",")
	switch number {
	case "A":
		// One value per ALT allele.
		if len(parts) != nROri-1 {
			return val, false
		}
		out := make([]string, 0, len(parts))
		for i := 1; i < nROri; i++ {
			if rm[i] {
				continue
			}
			out = append(out, parts[i-1])
		}
		return strings.Join(out, ","), true
	case "R":
		// One value per allele (REF + ALTs).
		if len(parts) != nROri {
			return val, false
		}
		out := make([]string, 0, len(parts))
		for i := 0; i < nROri; i++ {
			if i > 0 && rm[i] {
				continue
			}
			out = append(out, parts[i])
		}
		return strings.Join(out, ","), true
	default: // "G"
		return subsetGenotypeList(parts, rm, nROri)
	}
}

// subsetGenotypeList rewrites a Number=G value list (diploid or haploid) for the
// retained alleles, porting the Number=G branches of bcf_remove_allele_set. For
// the diploid layout the value at genotype index b*(b+1)/2 + a is kept only when
// both alleles a and b survive; for the haploid layout (one value per allele)
// the per-allele value is kept when its allele survives. A count that matches
// neither layout leaves the value unchanged.
func subsetGenotypeList(parts []string, rm []bool, nROri int) (string, bool) {
	nGOri := nROri * (nROri + 1) / 2
	if len(parts) == nGOri {
		out := make([]string, 0, len(parts))
		idx := 0
		for b := 0; b < nROri; b++ {
			for a := 0; a <= b; a++ {
				if idx >= len(parts) {
					break
				}
				if !(rmAt(rm, a) || rmAt(rm, b)) {
					out = append(out, parts[idx])
				}
				idx++
			}
		}
		return strings.Join(out, ","), true
	}
	if len(parts) == nROri { // haploid
		out := make([]string, 0, len(parts))
		for a := 0; a < nROri; a++ {
			if !rmAt(rm, a) {
				out = append(out, parts[a])
			}
		}
		return strings.Join(out, ","), true
	}
	return strings.Join(parts, ","), false
}

// rmAt reports whether allele a is flagged for removal, treating REF (a==0) and
// out-of-range indices as retained.
func rmAt(rm []bool, a int) bool {
	if a <= 0 || a >= len(rm) {
		return false
	}
	return rm[a]
}

// headerNumberMaps indexes the Number= attribute of every INFO and FORMAT
// header line, so removeAlleleSet can decide which fields are A/R/G-cardinal.
type headerNumberMaps struct {
	info   map[string]string
	format map[string]string
}

// headerNumberMapsFrom builds the INFO/FORMAT Number= lookup tables from a
// header's meta lines.
func headerNumberMapsFrom(meta []string) headerNumberMaps {
	m := headerNumberMaps{info: map[string]string{}, format: map[string]string{}}
	for _, line := range meta {
		switch headerKind(line) {
		case "##INFO":
			if id := headerID(line); id != "" {
				m.info[id] = headerField(line, "Number")
			}
		case "##FORMAT":
			if id := headerID(line); id != "" {
				m.format[id] = headerField(line, "Number")
			}
		}
	}
	return m
}
