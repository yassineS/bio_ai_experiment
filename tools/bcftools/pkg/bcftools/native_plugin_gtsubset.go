// Native port of the upstream `GTsubset` plugin (plugins/GTsubset.c). It outputs
// only those sites where every requested sample exclusively shares one genotype:
// all selected samples must carry the same GT (both alleles), and none of the
// unselected samples may carry that genotype. Missing genotypes always pass.
// The plugin is a per-record generic init/process plugin that emits the matching
// records unchanged.
package bcftools

import (
	"fmt"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("GTsubset", func() NativePlugin { return &gtSubsetPlugin{} }) }

// gtSubsetPlugin implements the `GTsubset` plugin. It keeps the per-sample
// selection mask across records but holds no cross-record accumulator, so it is
// treated as serial (its Process compares the whole sample row at once).
type gtSubsetPlugin struct {
	selected []bool // selected[i] true for samples named in -s
	nSel     int
}

// Name returns the plugin name.
func (p *gtSubsetPlugin) Name() string { return "GTsubset" }

// About returns the one-line description, matching GTsubset.c about().
func (p *gtSubsetPlugin) About() string {
	return "Output only sites where the requested samples all exclusively share a genotype (GT).\n"
}

// Parallel reports false: Process inspects every sample of a record together and
// the plugin is kept off the parallel pool for simplicity and determinism.
func (p *gtSubsetPlugin) Parallel() bool { return false }

// Init parses -s/--sample-list, validates GT presence, and builds the selection
// mask. The sample list is comma-separated; unknown samples are an error,
// mirroring upstream.
func (p *gtSubsetPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	var sampleList string
	for i := 0; i < len(args); i++ {
		a := args[i]
		switch a {
		case "-s", "--sample-list":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("GTsubset: %s requires an argument", a)
			}
			i++
			sampleList = args[i]
		case "-h", "--help":
			return nil, fmt.Errorf("GTsubset: help requested")
		default:
			return nil, fmt.Errorf("GTsubset: unsupported option %q", a)
		}
	}
	if len(hdr.Samples) == 0 {
		return nil, fmt.Errorf("GTsubset: no samples in input file")
	}
	if !hasFormatHeader(hdr.MetaInfo, "GT") {
		return nil, fmt.Errorf("GTsubset: GT not present in the header")
	}
	names := splitSampleList(sampleList)
	if len(names) == 0 {
		return nil, fmt.Errorf("GTsubset: sample specification not valid")
	}
	idx := make(map[string]int, len(hdr.Samples))
	for i, s := range hdr.Samples {
		idx[s] = i
	}
	p.selected = make([]bool, len(hdr.Samples))
	for _, n := range names {
		i, ok := idx[n]
		if !ok {
			return nil, fmt.Errorf("GTsubset: sample '%s' not in input vcf file", n)
		}
		p.selected[i] = true
		p.nSel++
	}
	return hdr, nil
}

// Process emits the record only when the selected samples exclusively share a
// genotype. The matching mirrors GTsubset.c process() operating directly on the
// raw htslib GT integer encoding (allele index and phase bit), because upstream
// compares those raw ints — including the phase bit and the literal "==0"
// missing test, so a phased and an unphased copy of the same genotype are
// treated as different. The reference genotype (a1,a2) is the first selected
// sample whose both raw values are non-zero; every selected sample must equal
// it and every unselected sample must differ. A sample with a raw "0" in either
// slot (an unphased missing allele) always passes.
func (p *gtSubsetPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	nsmp := len(v.Samples)

	a1, a2 := 0, 0
	for gt := 0; a1 == 0 || a2 == 0; gt++ {
		if gt == nsmp {
			break
		}
		if !p.selected[gt] {
			continue
		}
		a1, a2 = gtRawPair(v, gt)
	}

	allPass := true
	for i := 0; i < nsmp; i++ {
		b1, b2 := gtRawPair(v, i)
		// A raw 0 in either slot (unphased missing) always passes.
		if b1 == 0 || b2 == 0 {
			continue
		}
		if p.selected[i] {
			if b1 == a1 && b2 == a2 {
				continue
			}
			allPass = false
			break
		}
		// Unselected sample must differ from the shared genotype.
		if b1 != a1 || b2 != a2 {
			continue
		}
		allPass = false
		break
	}
	if allPass {
		return []*vcf.Variant{v}, nil
	}
	return nil, nil
}

// Destroy is a no-op for GTsubset.
func (p *gtSubsetPlugin) Destroy() error { return nil }

// gtRawPair returns the raw htslib GT integers for sample i's first two allele
// slots. The encoding matches bcf_gt_*: an allele index a maps to (a+1)<<1, the
// phase bit (1) is OR-ed in when the genotype is phased (any "|" separator), a
// missing allele decodes to index -1 (so unphased missing = 0, phased missing =
// 1), and a haploid genotype's second slot is the vector-end sentinel (a large
// distinct value that is never 0 and never equals a real allele).
func gtRawPair(v *vcf.Variant, i int) (int, int) {
	gt, ok := sampleGT(v, i)
	if !ok || len(gt.alleles) == 0 {
		return 0, 0
	}
	phased := 0
	for _, p := range gt.phased {
		if p {
			phased = 1
			break
		}
	}
	raw := func(a int) int { return ((a + 1) << 1) | phased }
	b1 := raw(gt.alleles[0])
	if len(gt.alleles) == 1 {
		return b1, gtVectorEndRaw
	}
	return b1, raw(gt.alleles[1])
}

// gtVectorEndRaw is the sentinel for a haploid genotype's absent second slot,
// mirroring bcf_int32_vector_end. It is distinct from 0 and from every real
// allele encoding so a haploid genotype only matches another identical haploid.
const gtVectorEndRaw = 1 << 30

// splitSampleList splits a comma-separated sample list, trimming spaces and
// dropping empty entries.
func splitSampleList(s string) []string {
	if s == "" {
		return nil
	}
	var out []string
	start := 0
	for i := 0; i <= len(s); i++ {
		if i == len(s) || s[i] == ',' {
			tok := s[start:i]
			if tok != "" {
				out = append(out, tok)
			}
			start = i + 1
		}
	}
	return out
}
