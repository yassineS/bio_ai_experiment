// Shared helpers for the sample/indel statistics plugins (smpl-stats,
// indel-stats). These mirror the small htslib utilities the upstream plugins
// rely on: per-sample genotype parsing (bcf_get_genotypes + the plugins'
// parse_genotype), per-allele variant typing (bcf_get_variant_type), the
// star-allele scan, and allele-count computation (bcf_calc_ac).
package bcftools

import (
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// statsReportWriter resolves the destination for a stats plugin's textual report
// given the -o/--output FILE selection and the host stdout writer. When
// outputFile is empty or "-", the report goes to the host stdout (out) and the
// returned closer is a no-op. Otherwise the file is created (truncated) and the
// returned closer closes it. This mirrors upstream's report_stats() open logic
// (`!output_fname || !strcmp("-",output_fname) ? stdout : fopen(...)`), so the
// report bytes are identical whether written to stdout or a file — the CMD line
// echoes the verbatim argv (including -o) in both cases, exactly as upstream.
func statsReportWriter(outputFile string, out io.Writer) (io.Writer, func() error, error) {
	if outputFile == "" || outputFile == "-" {
		return out, func() error { return nil }, nil
	}
	f, err := os.Create(outputFile)
	if err != nil {
		return nil, nil, fmt.Errorf("could not open the file for writing: %s", outputFile)
	}
	return f, f.Close, nil
}

// genotypeKind classifies a parsed sample genotype, mirroring the return codes
// of the plugins' parse_genotype helper.
type genotypeKind int

const (
	gtDiploid genotypeKind = iota // two present alleles
	gtHemi                        // a single present allele (haploid / vector-end)
	gtMissing                     // a missing allele was encountered
)

// parseGenotypeAlleles parses sample i's GT into a fixed two-element allele
// array and a kind, mirroring parse_genotype() in smpl-stats.c / indel-stats.c.
// For a haploid genotype, als[1] is set equal to als[0] (as upstream does) and
// the kind is gtHemi. A missing first or second allele yields gtMissing.
func parseGenotypeAlleles(v *vcf.Variant, i int) (als [2]int, kind genotypeKind) {
	gt, ok := sampleGT(v, i)
	if !ok || len(gt.alleles) == 0 {
		return als, gtMissing
	}
	if gt.alleles[0] == missingAllele {
		return als, gtMissing
	}
	als[0] = gt.alleles[0]
	if len(gt.alleles) == 1 {
		als[1] = als[0]
		return als, gtHemi
	}
	if gt.alleles[1] == missingAllele {
		return als, gtMissing
	}
	als[1] = gt.alleles[1]
	return als, gtDiploid
}

// starAlleleIndex returns the index of the '*' ALT allele (1-based into the
// allele list), or -1 if none, matching the star_allele scan in the plugins.
func starAlleleIndex(v *vcf.Variant) int {
	for i, alt := range v.Alt {
		if alt == "*" {
			return i + 1
		}
	}
	return -1
}

// altVariantType returns the variant-type bitmask of allele a (a>=1 indexes
// v.Alt[a-1]; a==0 is the reference and types as 0), mirroring htslib's
// bcf_get_variant_type for a single allele.
func altVariantType(v *vcf.Variant, a int) int {
	if a <= 0 || a-1 >= len(v.Alt) {
		return 0
	}
	return variantTypeBit(v.Ref, v.Alt[a-1])
}

// indelAlleleLen returns the net length change of allele a relative to the
// reference (len(alt)-len(ref)), mirroring htslib's rec->d.var[a].n for indels.
// It is only meaningful when allele a types as an indel.
func indelAlleleLen(v *vcf.Variant, a int) int {
	if a <= 0 || a-1 >= len(v.Alt) {
		return 0
	}
	return len(v.Alt[a-1]) - len(v.Ref)
}

// computeACWithRef returns the allele-count array (length nAllele, ac[0] is the
// reference count), preferring INFO/AC+AN and falling back to the genotypes,
// mirroring bcf_calc_ac(BCF_UN_INFO|BCF_UN_FMT). When INFO/AC is present it
// provides the ALT counts (ac[1..]) and ac[0] is AN minus their sum; otherwise
// every count is derived from the called genotypes.
func computeACWithRef(v *vcf.Variant, nAllele int) []int {
	ac := make([]int, nAllele)
	if acStr, ok := v.Info["AC"]; ok {
		parts := strings.Split(acStr, ",")
		good := len(parts) == nAllele-1
		sum := 0
		if good {
			for k, s := range parts {
				n, err := strconv.Atoi(s)
				if err != nil {
					good = false
					break
				}
				ac[k+1] = n
				sum += n
			}
		}
		if good {
			if anStr, anOK := v.Info["AN"]; anOK {
				if an, err := strconv.Atoi(anStr); err == nil {
					ac[0] = an - sum
					return ac
				}
			}
		}
		// Fall through to GT-based counting if AC/AN were not usable.
		for i := range ac {
			ac[i] = 0
		}
	}
	for i := range v.Samples {
		gt, ok := sampleGT(v, i)
		if !ok {
			continue
		}
		for _, a := range gt.alleles {
			if a == missingAllele || a < 0 || a >= nAllele {
				continue
			}
			ac[a]++
		}
	}
	return ac
}
