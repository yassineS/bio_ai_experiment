// Shared genotype and value helpers for native bcftools plugins. These mirror
// the relevant pieces of htslib's GT encoding (bcf_gt_allele, phasing,
// bcf_gt_missing) but operate on the textual VCF representation used by the
// pkg/htsgo/vcf record model.
package bcftools

import (
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// missingAllele is the sentinel for a missing GT allele ("." in a genotype).
const missingAllele = -1

// genotype is a parsed FORMAT/GT value: one allele index per ploidy (with
// missingAllele for "."), plus a parallel phased flag for each allele after
// the first. The phase boundary preceding allele i (i>=1) is phased[i].
type genotype struct {
	alleles []int  // allele indices; missingAllele (-1) for "."
	phased  []bool // phased[i] true if the separator before alleles[i] was '|'
}

// parseGT parses a textual GT such as "0/1", "0|1", "./.", "1", or "0/0/1".
// An empty or "." whole field yields a single missing allele. The boolean
// reports whether the field was a parseable genotype (false for fields that
// are not GT-shaped, e.g. ".").
func parseGT(s string) (genotype, bool) {
	if s == "" {
		return genotype{}, false
	}
	gt := genotype{}
	// Split on / or | while remembering the separators for phasing.
	start := 0
	phasedSep := false // separator preceding the current allele
	for i := 0; i <= len(s); i++ {
		if i == len(s) || s[i] == '/' || s[i] == '|' {
			tok := s[start:i]
			var allele int
			if tok == "." || tok == "" {
				allele = missingAllele
			} else {
				n, err := strconv.Atoi(tok)
				if err != nil {
					return genotype{}, false
				}
				allele = n
			}
			gt.alleles = append(gt.alleles, allele)
			gt.phased = append(gt.phased, phasedSep)
			if i < len(s) {
				phasedSep = s[i] == '|'
			}
			start = i + 1
		}
	}
	return gt, true
}

// String renders the genotype back to VCF text, using '|' before alleles whose
// phased flag is set and '/' otherwise. The first allele never carries a
// separator. Missing alleles render as ".".
func (g genotype) String() string {
	if len(g.alleles) == 0 {
		return "."
	}
	var b strings.Builder
	for i, a := range g.alleles {
		if i > 0 {
			if g.phased[i] {
				b.WriteByte('|')
			} else {
				b.WriteByte('/')
			}
		}
		if a == missingAllele {
			b.WriteByte('.')
		} else {
			b.WriteString(strconv.Itoa(a))
		}
	}
	return b.String()
}

// ploidy returns the number of alleles in the genotype.
func (g genotype) ploidy() int { return len(g.alleles) }

// isFullyMissing reports whether every allele is missing.
func (g genotype) isFullyMissing() bool {
	for _, a := range g.alleles {
		if a != missingAllele {
			return false
		}
	}
	return len(g.alleles) > 0
}

// nMissing returns the number of missing alleles.
func (g genotype) nMissing() int {
	n := 0
	for _, a := range g.alleles {
		if a == missingAllele {
			n++
		}
	}
	return n
}

// sampleGT returns the parsed GT for sample i, or (zero, false) if absent.
func sampleGT(v *vcf.Variant, i int) (genotype, bool) {
	if i >= len(v.Samples) {
		return genotype{}, false
	}
	s, ok := v.Samples[i].Data["GT"]
	if !ok {
		return genotype{}, false
	}
	return parseGT(s)
}

// formatVCFFloat renders a float for an INFO/FORMAT field the way htslib's
// VCF writer does: the value is stored as a 32-bit C float, then printed with
// the shortest round-tripping decimal. This matches bcftools byte-for-byte for
// the AF/MAF/HWE/ExcHet/VAF values produced by fill-tags. Whole numbers print
// without a trailing ".0" (matching kputd's integer fast path).
func formatVCFFloat(v float64) string {
	f32 := float32(v)
	if math.IsNaN(float64(f32)) {
		return "nan"
	}
	f := float64(f32)
	if f == math.Trunc(f) && !math.IsInf(f, 0) && math.Abs(f) < 1e15 {
		// Preserve negative zero, which htslib's %g prints as "-0" (e.g. GL
		// derived from PL=0 via -0.1*0 == -0.0).
		if f == 0 && math.Signbit(f) {
			return "-0"
		}
		return strconv.FormatInt(int64(f), 10)
	}
	// htslib's VCF writer prints floats with C's `%g` conversion, which uses
	// six significant digits with trailing zeros stripped (e.g. 1/6 ->
	// "0.166667", 3/4 -> "0.75"). The value is first narrowed to a 32-bit C
	// float, so the rounding matches upstream byte-for-byte.
	return strconv.FormatFloat(f, 'g', 6, 64)
}

// ensureFormatTag makes sure tag is present in v.Format (appending it if not).
// (setInfo and variantTypeMask are shared with the rest of the package and
// defined in call.go / view.go respectively.)
func ensureFormatTag(v *vcf.Variant, tag string) {
	for _, f := range v.Format {
		if f == tag {
			return
		}
	}
	v.Format = append(v.Format, tag)
}
