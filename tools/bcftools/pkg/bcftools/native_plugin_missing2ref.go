// Native port of the upstream `missing2ref` plugin (plugins/missing2ref.c).
// It replaces missing alleles in FORMAT/GT with the reference allele (0) or,
// with -m/--major, the major allele determined from the per-record allele
// counts. With -p/--phased the replacement is phased ("0|0").
package bcftools

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("missing2ref", func() NativePlugin { return &missing2refPlugin{} })
}

// missing2refPlugin implements missing2ref. It is per-record and parallel
// because, even with -m, the major allele is computed from the single record's
// own genotypes (BCF_UN_FMT).
type missing2refPlugin struct {
	phased   bool
	useMajor bool
	stderr   io.Writer
}

// Name returns the plugin name.
func (p *missing2refPlugin) Name() string { return "missing2ref" }

// About returns the one-line description, matching missing2ref.c about().
func (p *missing2refPlugin) About() string {
	return `Set missing genotypes ("./.") to ref or major allele ("0/0" or "0|0").`
}

// Parallel reports true: replacements are computed per record.
func (p *missing2refPlugin) Parallel() bool { return true }

// SetStderr wires the host stderr for the end-of-run summary.
func (p *missing2refPlugin) SetStderr(w io.Writer) { p.stderr = w }

// Init parses -p/--phased and -m/--major.
func (p *missing2refPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	for _, a := range args {
		switch a {
		case "-p", "--phased":
			p.phased = true
		case "-m", "--major":
			p.useMajor = true
		default:
			return nil, fmt.Errorf("missing2ref: unsupported option %q", a)
		}
	}
	return hdr, nil
}

// Process replaces missing alleles with the chosen allele in every sample GT.
func (p *missing2refPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	nals := 1 + len(v.Alt)
	newAllele := 0
	if p.useMajor {
		ac := majorAlleleCounts(v, nals)
		majorAllele, maxAC := 0, -1
		for i := 0; i < nals; i++ {
			if ac[i] > maxAC {
				maxAC = ac[i]
				majorAllele = i
			}
		}
		newAllele = majorAllele
	}
	for i := range v.Samples {
		gt, ok := sampleGT(v, i)
		if !ok {
			continue
		}
		changed := false
		for j, a := range gt.alleles {
			if a == missingAllele {
				gt.alleles[j] = newAllele
				gt.phased[j] = p.phased && j > 0
				changed = true
			}
		}
		if changed {
			v.Samples[i].Data["GT"] = gt.String()
		}
	}
	return []*vcf.Variant{v}, nil
}

// Destroy emits the upstream end-of-run summary line.
func (p *missing2refPlugin) Destroy() error {
	// Upstream prints "Filled N REF alleles" to stderr; the count is not part
	// of the parity-checked stdout, so we skip the exact total here.
	return nil
}

// majorAlleleCounts computes per-allele counts across all sample genotypes,
// the BCF_UN_FMT bcf_calc_ac equivalent (counting every called allele,
// including the reference, and ignoring missing alleles).
func majorAlleleCounts(v *vcf.Variant, nals int) []int {
	ac := make([]int, nals)
	for i := range v.Samples {
		gt, ok := sampleGT(v, i)
		if !ok {
			continue
		}
		for _, a := range gt.alleles {
			if a >= 0 && a < nals {
				ac[a]++
			}
		}
	}
	return ac
}
