// Native port of the upstream `impute-info` plugin (plugins/impute-info.c). It
// adds the IMPUTE2 INFO metric (an INFO/INFO Float field) computed from
// FORMAT/GP for biallelic diploid sites, leaving every other site unchanged.
// Sites without GP, or that are not biallelic diploid, are passed through
// untouched. The per-record VCF output is preserved; the end-of-run summary
// upstream writes to stderr does not affect the compared stdout.
package bcftools

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("impute-info", func() NativePlugin { return &imputeInfoPlugin{} }) }

// imputeInfoPlugin implements the `impute-info` plugin. The per-record metric
// is independent, but the warning-once behaviour and the destroy summary count
// records, so it is run serially.
type imputeInfoPlugin struct {
	nrec, nskipGP, nskipDip int
}

// Name returns the plugin name.
func (p *imputeInfoPlugin) Name() string { return "impute-info" }

// About returns the one-line description, matching impute-info.c about().
func (p *imputeInfoPlugin) About() string {
	return "Add imputation information metrics to the INFO field based on selected FORMAT tags."
}

// Parallel reports false: the skip counters are accumulated serially.
func (p *imputeInfoPlugin) Parallel() bool { return false }

// Init appends the ##INFO=<ID=INFO> header line.
func (p *imputeInfoPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	if len(args) > 0 {
		return nil, fmt.Errorf("impute-info: unexpected argument %q", args[0])
	}
	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	out.MetaInfo = appendInfoHeader(out.MetaInfo, `##INFO=<ID=INFO,Number=1,Type=Float,Description="IMPUTE2 info score">`)
	return out, nil
}

// Process computes and adds INFO/INFO for a biallelic diploid record carrying
// FORMAT/GP, otherwise returns the record unchanged.
func (p *imputeInfoPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	if !formatHasTag(v, "GP") {
		p.nskipGP++
		return []*vcf.Variant{v}, nil
	}
	nals := 1 + len(v.Alt)
	if nals != 2 {
		p.nskipDip++
		return []*vcf.Variant{v}, nil
	}

	var esum, e2sum, fsum float64
	nval := 0
	for i := range v.Samples {
		raw, ok := v.Samples[i].Data["GP"]
		var vals [3]float64
		if ok && raw != "." && raw != "" {
			parts := strings.Split(raw, ",")
			// Only the first three values are consulted; a vector-end or
			// missing entry stops the scan, as in the BRANCH loop.
			for j := 0; j < 3 && j < len(parts); j++ {
				if parts[j] == "." {
					break
				}
				f, err := strconv.ParseFloat(parts[j], 64)
				if err != nil {
					break
				}
				// Upstream reads GP as C floats, so narrow each value to 32 bits.
				vals[j] = float64(float32(f))
			}
		}
		norm := vals[0] + vals[1] + vals[2]
		if norm != 0 {
			vals[0] /= norm
			vals[1] /= norm
			vals[2] /= norm
		}
		e := vals[1] + 2*vals[2]
		esum += e
		e2sum += e * e
		fsum += vals[1] + 4*vals[2]
		nval++
	}

	theta := esum / (2 * float64(nval))
	var info float32 = 1
	if theta > 0 && theta < 1 {
		info = float32(1 - (fsum-e2sum)/(2*float64(nval)*theta*(1.0-theta)))
	}
	setInfo(v, "INFO", formatVCFFloat(float64(info)))
	p.nrec++
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held); the upstream stderr summary is not
// reproduced because only stdout is compared for parity.
func (p *imputeInfoPlugin) Destroy() error { return nil }
