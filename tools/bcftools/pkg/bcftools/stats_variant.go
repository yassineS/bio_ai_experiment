package bcftools

import (
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// VCF variant-type bit flags, mirroring htslib's bcf_get_variant_type. Only the
// subset that `bcftools stats` consumes is modelled here.
const (
	vcfRef    = 0
	vcfSNP    = 1 << 0
	vcfMNP    = 1 << 1
	vcfIndel  = 1 << 2
	vcfOther  = 1 << 3
	vcfBND    = 1 << 4
	vcfOvrlap = 1 << 5
)

// statsVariantType classifies a single REF/ALT allele pair, returning the
// htslib variant-type bitmask and the signed length difference (positive for
// insertions, negative for deletions, as bcf_variant_length reports). It mirrors
// htslib's bcf_set_variant_type so that stats counters bucket exactly as
// upstream does.
func statsVariantType(ref, alt string) (typ int, n int) {
	if alt == "*" {
		return vcfOvrlap, 0
	}
	if alt == "" || alt == "." {
		return vcfRef, 0
	}
	if len(ref) == 1 && len(alt) == 1 {
		if alt[0] == '.' || equalFold1(ref[0], alt[0]) {
			return vcfRef, 0
		}
		if alt[0] == 'X' || alt[0] == 'x' {
			return vcfRef, 0
		}
		return vcfSNP, 1
	}
	if alt[0] == '<' {
		return vcfOther, 0
	}
	if alt[0] == ']' || alt[0] == '[' {
		return vcfBND, 0
	}

	// Trim the common prefix.
	r, a := 0, 0
	for r < len(ref) && a < len(alt) && toUpperByte(ref[r]) == toUpperByte(alt[a]) {
		r++
		a++
	}
	switch {
	case a < len(alt) && r >= len(ref):
		// ALT longer: insertion.
		if alt[len(alt)-1] == ']' || alt[len(alt)-1] == '[' {
			return vcfBND, 0
		}
		return vcfIndel, len(alt) - len(ref)
	case r < len(ref) && a >= len(alt):
		// REF longer: deletion.
		return vcfIndel, len(alt) - len(ref)
	case r >= len(ref) && a >= len(alt):
		return vcfRef, 0
	}

	// Trim the common suffix.
	re, ae := len(ref)-1, len(alt)-1
	if alt[ae] == ']' || alt[ae] == '[' {
		return vcfBND, 0
	}
	for re > r && ae > a && toUpperByte(ref[re]) == toUpperByte(alt[ae]) {
		re--
		ae--
	}
	switch {
	case ae == a:
		if re == r {
			return vcfSNP, 1
		}
		if toUpperByte(ref[re]) == toUpperByte(alt[ae]) {
			return vcfIndel, -(re - r)
		}
		return vcfOther, -(re - r)
	case re == r:
		if toUpperByte(ref[re]) == toUpperByte(alt[ae]) {
			return vcfIndel, ae - a
		}
		return vcfOther, ae - a
	}
	if re-r == ae-a {
		return vcfMNP, 0
	}
	return vcfOther, 0
}

// lineVariantTypes returns the union of per-allele variant types for a record,
// matching htslib's bcf_get_variant_types: REF is its own zero-valued type and
// is only reported when no alternate allele is present.
func lineVariantTypes(v *vcf.Variant) int {
	typ := 0
	any := false
	for _, alt := range v.Alt {
		t, _ := statsVariantType(v.Ref, alt)
		if t == vcfRef {
			continue
		}
		typ |= t
		any = true
	}
	if !any {
		return vcfRef
	}
	return typ
}

// genotype types, mirroring htslib's vcfutils.h GT_* constants.
const (
	gtHomRR = iota
	gtHomAA
	gtHetRA
	gtHetAA
	gtHaplR
	gtHaplA
	gtUnkn
)

// classifyGT parses a textual GT string the way htslib's bcf_gt_type does,
// returning the genotype type and the (0-based) indices of the first and second
// non-reference alleles. ial/jal index into REF(0)/ALT(1..) just like
// bcf_gt_type's outputs.
func classifyGT(gt string) (typ, ial, jal int) {
	if gt == "" {
		return gtUnkn, 0, 0
	}
	fields := strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' })
	if len(fields) == 0 {
		return gtUnkn, 0, 0
	}
	hasRef, hasAlt := false, 0
	i1, j1 := 0, 0 // 1-based allele values (>1 means an ALT)
	nals := 0
	for _, f := range fields {
		if f == "." {
			return gtUnkn, 0, 0
		}
		val, err := strconv.Atoi(f)
		if err != nil {
			return gtUnkn, 0, 0
		}
		tmp := val + 1 // htslib encodes allele a as (a+1), so REF=>1, ALT1=>2 ...
		if tmp > 1 {
			if i1 == 0 {
				i1 = tmp
				hasAlt = 1
			} else if tmp != i1 {
				if tmp < i1 {
					j1 = i1
					i1 = tmp
				} else {
					j1 = tmp
				}
				hasAlt = 2
			}
		} else {
			hasRef = true
		}
		nals++
	}
	ial = 0
	jal = 0
	if i1 > 0 {
		ial = i1 - 1
	}
	if j1 > 0 {
		jal = j1 - 1
	}
	switch {
	case nals == 0:
		return gtUnkn, ial, jal
	case nals == 1:
		if hasRef {
			return gtHaplR, ial, jal
		}
		return gtHaplA, ial, jal
	case !hasRef:
		if hasAlt == 1 {
			return gtHomAA, ial, jal
		}
		return gtHetAA, ial, jal
	case hasAlt == 0:
		return gtHomRR, ial, jal
	default:
		return gtHetRA, ial, jal
	}
}

// statsAC returns the per-allele counts the way bcftools stats derives them for
// AF binning: it prefers INFO/AC together with INFO/AN (ac[0] = AN - sum(AC)),
// and only falls back to recomputing from genotypes when samples were requested
// and the INFO tags are absent. The returned slice is indexed by allele
// (0=REF). ok is false when no count basis is available.
func (r *statsResult) statsAC(v *vcf.Variant) (ac []int, ok bool) {
	useSamples := r.useSamples
	nAllele := len(v.Alt) + 1
	r.acBuf = intBuf(r.acBuf, nAllele)
	ac = r.acBuf
	acRaw, hasAC := v.Info["AC"]
	anRaw, hasAN := v.Info["AN"]
	if hasAC && hasAN {
		an, errAN := strconv.Atoi(strings.TrimSpace(anRaw))
		parts := strings.Split(acRaw, ",")
		if errAN == nil && len(parts) == nAllele-1 {
			nac := 0
			good := true
			for i, p := range parts {
				val, err := strconv.Atoi(strings.TrimSpace(p))
				if err != nil {
					good = false
					break
				}
				ac[i+1] = val
				nac += val
			}
			if good && an >= nac {
				ac[0] = an - nac
				return ac, true
			}
			// Partial fill from a failed/short AC parse must not leak into the
			// genotype-fallback path below; reset before reusing the buffer.
			for i := range ac {
				ac[i] = 0
			}
		}
	}
	if useSamples {
		alts := computeACFromGT(v, nAllele-1)
		an := totalAN(v)
		nac := 0
		for i, c := range alts {
			ac[i+1] = c
			nac += c
		}
		if an >= nac {
			ac[0] = an - nac
		}
		return ac, true
	}
	return ac, false
}

// toUpperByte upper-cases an ASCII byte.
func toUpperByte(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 32
	}
	return b
}
