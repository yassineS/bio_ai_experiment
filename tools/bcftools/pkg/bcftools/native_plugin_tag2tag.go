// Native port of the upstream `tag2tag` plugin (plugins/tag2tag.c) for the
// common GL/PL/GP/GT conversions selected with --SRC-to-DST. These convert
// between the similar genotype-likelihood tags:
//
//	GL  : log10 likelihoods (Float, Number=G)
//	PL  : phred-scaled likelihoods (Integer, Number=G)
//	GP  : genotype probabilities (Float, Number=G)
//	GT  : hard-called genotype (String, Number=1)
//
// Options -r/--replace (drop the source tag) and -t/--threshold (GP->GT
// hard-call threshold) are supported. The localized-allele family
// (--XX-to-LXX, --LXX-to-XX) and the QR,QA->QS conversion require the htslib
// localized-allele machinery and are reported as unsupported in batch 1.
package bcftools

import (
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("tag2tag", func() NativePlugin { return &tag2tagPlugin{} }) }

// tag2tagPlugin implements the GL/PL/GP/GT conversions. It is per-record and
// parallel.
type tag2tagPlugin struct {
	src, dst string
	dropSrc  bool
	gpThresh float64
}

// Name returns the plugin name.
func (p *tag2tagPlugin) Name() string { return "tag2tag" }

// About returns the one-line description, matching tag2tag.c about().
func (p *tag2tagPlugin) About() string {
	return "Convert between similar tags, such as GL,PL,GP or QR,QA,QS or localized alleles, eg PL,LPL."
}

// Parallel reports true.
func (p *tag2tagPlugin) Parallel() bool { return true }

// Init parses the --SRC-to-DST selector and -r/-t options.
func (p *tag2tagPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	// Upstream tag2tag.c leaves gp_th calloc'd to 0 by default (unlike the
	// usage text's "[0.1]"), so the GP->GT hard call requires GP==1 unless -t
	// is given explicitly.
	p.gpThresh = 0
	for i := 0; i < len(args); i++ {
		a := args[i]
		switch {
		case a == "-r" || a == "--replace":
			p.dropSrc = true
		case a == "-t" || a == "--threshold":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("tag2tag: -t requires an argument")
			}
			i++
			v, err := strconv.ParseFloat(args[i], 64)
			if err != nil || v < 0 || v > 1 {
				return nil, fmt.Errorf("tag2tag: expected value between 0 and 1 for -t, got %q", args[i])
			}
			p.gpThresh = v
		case a == "-s" || a == "--skip-nalt" || a == "-d" || a == "--defaults":
			return nil, fmt.Errorf("tag2tag: option %q applies only to the localized-allele modes, which are not supported in the native plugin", a)
		case strings.HasPrefix(a, "--") && strings.Contains(strings.ToUpper(a), "-TO-"):
			if err := p.parseSelector(a); err != nil {
				return nil, err
			}
		default:
			return nil, fmt.Errorf("tag2tag: unsupported option %q", a)
		}
	}
	if p.src == "" {
		return nil, fmt.Errorf("tag2tag: a conversion such as --PL-to-GT must be given")
	}
	if !hasFormatHeader(hdr.MetaInfo, p.src) {
		return nil, fmt.Errorf("tag2tag: the source tag does not exist: %s", p.src)
	}

	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	if p.dropSrc {
		out.MetaInfo = removeFormatHeader(out.MetaInfo, p.src)
	}
	out.MetaInfo = appendInfoHeader(out.MetaInfo, dstHeaderLine(p.dst))
	return out, nil
}

// parseSelector parses a --SRC-to-DST option into src/dst, rejecting the
// localized-allele and QR/QA selectors.
func (p *tag2tagPlugin) parseSelector(opt string) error {
	up := strings.ToUpper(opt[2:]) // strip leading --
	switch up {
	case "XX-TO-LXX", "LXX-TO-XX", "LPL-TO-PL", "LAD-TO-AD",
		"PL-TO-LPL", "AD-TO-LAD", "QR-QA-TO-QS":
		return fmt.Errorf("tag2tag: conversion %s (localized alleles / QR,QA) is not supported in the native plugin", opt)
	}
	parts := strings.Split(up, "-TO-")
	if len(parts) != 2 {
		return fmt.Errorf("tag2tag: could not parse conversion %q", opt)
	}
	src, dst := parts[0], parts[1]
	valid := map[string]bool{"GL": true, "PL": true, "GP": true, "GT": true}
	if !valid[src] || !valid[dst] || src == dst {
		return fmt.Errorf("tag2tag: the conversion is not supported: %s", opt)
	}
	if src == "GT" {
		return fmt.Errorf("tag2tag: GT cannot be a source tag")
	}
	p.src, p.dst = src, dst
	return nil
}

// Process converts the source tag to the destination tag for one record.
func (p *tag2tagPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	if !formatHasTag(v, p.src) {
		return []*vcf.Variant{v}, nil
	}
	nals := 1 + len(v.Alt)

	for i := range v.Samples {
		raw, ok := v.Samples[i].Data[p.src]
		if !ok || raw == "." || raw == "" {
			continue
		}
		// Convert the per-sample source values to GL (log10 likelihoods).
		gl, ok := p.toGL(raw)
		if !ok {
			continue
		}
		switch p.dst {
		case "GL":
			v.Samples[i].Data["GL"] = formatGLList(gl)
		case "PL":
			v.Samples[i].Data["PL"] = glToPL(gl)
		case "GP":
			v.Samples[i].Data["GP"] = glToGP(gl)
		case "GT":
			v.Samples[i].Data["GT"] = glToGT(gl, nals, p.gpThresh)
		}
	}

	for _, tag := range []string{"GL", "PL", "GP", "GT"} {
		if tag == p.dst {
			ensureFormatTag(v, tag)
		}
	}
	if p.dropSrc {
		removeFormatTag(v, p.src)
	}
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *tag2tagPlugin) Destroy() error { return nil }

// toGL parses the source value and returns it as GL (log10 likelihoods),
// mirroring the "convert to GL" step in tag2tag.c process(). GP values are
// log10'd (0 -> -99); PL values are scaled by -0.1; GL is used as-is.
func (p *tag2tagPlugin) toGL(raw string) ([]float64, bool) {
	parts := strings.Split(raw, ",")
	out := make([]float64, len(parts))
	for i, s := range parts {
		if s == "." {
			out[i] = math.NaN()
			continue
		}
		f, err := strconv.ParseFloat(s, 64)
		if err != nil {
			return nil, false
		}
		switch p.src {
		case "GP":
			if f != 0 {
				out[i] = math.Log10(f)
			} else {
				out[i] = -99
			}
		case "PL":
			out[i] = -0.1 * f
		case "GL":
			out[i] = f
		}
	}
	return out, true
}

// formatGLList renders GL values, matching htslib's float %g (6 sig digits).
func formatGLList(gl []float64) string {
	parts := make([]string, len(gl))
	for i, x := range gl {
		if math.IsNaN(x) {
			parts[i] = "."
		} else {
			parts[i] = formatVCFFloat(x)
		}
	}
	return strings.Join(parts, ",")
}

// glToPL converts GL back to PL: PL = round(-10*GL), mirroring lroundf(-10*gl).
func glToPL(gl []float64) string {
	parts := make([]string, len(gl))
	for i, x := range gl {
		if math.IsNaN(x) {
			parts[i] = "."
		} else {
			parts[i] = strconv.Itoa(int(math.Round(-10 * x)))
		}
	}
	return strings.Join(parts, ",")
}

// glToGP converts GL to normalized GP: 10^GL then divide by the sum.
func glToGP(gl []float64) string {
	gp := make([]float64, len(gl))
	sum := 0.0
	for i, x := range gl {
		if math.IsNaN(x) {
			gp[i] = math.NaN()
			continue
		}
		gp[i] = math.Pow(10, x)
		sum += gp[i]
	}
	parts := make([]string, len(gl))
	for i := range gp {
		if math.IsNaN(gp[i]) {
			parts[i] = "."
		} else if sum > 0 {
			parts[i] = formatVCFFloat(gp[i] / sum)
		} else {
			parts[i] = formatVCFFloat(gp[i])
		}
	}
	return strings.Join(parts, ",")
}

// glToGT hard-calls a genotype from GL, mirroring the GP->GT branch of
// tag2tag.c: convert GL to normalized GP, pick the max, and apply the
// threshold (1 - gpThresh) below which the call becomes missing.
func glToGT(gl []float64, nals int, gpThresh float64) string {
	// Convert to GP (normalized 10^GL). htslib stores the likelihoods as C
	// floats, so the pow/sum/divide all happen in 32-bit precision; doing the
	// same here is what makes a near-certain call (e.g. PL 120,0,130) round to
	// GP==1 and get hard-called rather than falling just under the threshold.
	gp := make([]float64, len(gl))
	var sum float32
	for i, x := range gl {
		if math.IsNaN(x) {
			gp[i] = math.NaN()
			continue
		}
		f := float32(math.Pow(10, x))
		gp[i] = float64(f)
		sum += f
	}
	if sum > 0 {
		for i := range gp {
			if !math.IsNaN(gp[i]) {
				gp[i] = float64(float32(gp[i]) / sum)
			}
		}
	}
	if len(gp) == 0 || math.IsNaN(gp[0]) {
		return "./."
	}
	jmax := 0
	n1 := len(gp)
	for j := 1; j < n1; j++ {
		if math.IsNaN(gp[j]) {
			n1 = j
			break
		}
		if gp[j] > gp[jmax] {
			jmax = j
		}
	}

	// Haploid genotype: number of GP values == number of alleles.
	if n1 == nals {
		if gp[jmax] < 1-gpThresh {
			return "."
		}
		return strconv.Itoa(jmax)
	}
	if gp[jmax] < 1-gpThresh {
		return "./."
	}
	if jmax == 0 {
		return "0/0"
	}
	a, b := gt2alleles(jmax)
	return strconv.Itoa(a) + "/" + strconv.Itoa(b)
}

// gt2alleles maps the diploid GP/PL index back to its two alleles, mirroring
// htslib bcf_gt2alleles: index = b*(b+1)/2 + a with a <= b.
func gt2alleles(idx int) (int, int) {
	b := 0
	for (b+1)*(b+2)/2 <= idx {
		b++
	}
	a := idx - b*(b+1)/2
	return a, b
}

// dstHeaderLine returns the ##FORMAT header line for a destination tag.
func dstHeaderLine(tag string) string {
	switch tag {
	case "GP":
		return `##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype probabilities">`
	case "GL":
		return `##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods">`
	case "PL":
		return `##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">`
	case "GT":
		return `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`
	}
	return ""
}

// removeFormatHeader drops the ##FORMAT line for id from meta.
func removeFormatHeader(meta []string, id string) []string {
	out := meta[:0:0]
	for _, m := range meta {
		if headerKind(m) == "##FORMAT" && headerID(m) == id {
			continue
		}
		out = append(out, m)
	}
	return out
}

// removeFormatTag drops a FORMAT tag (and its per-sample values) from v.
func removeFormatTag(v *vcf.Variant, tag string) {
	idx := -1
	for i, f := range v.Format {
		if f == tag {
			idx = i
			break
		}
	}
	if idx < 0 {
		return
	}
	v.Format = append(v.Format[:idx], v.Format[idx+1:]...)
	for i := range v.Samples {
		delete(v.Samples[i].Data, tag)
	}
}
