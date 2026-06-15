// Native registration stubs for the batch-6 trio plugins whose output is not
// bit-reproducible against upstream and is therefore deliberately reported as
// unsupported (rather than emitting silently wrong results). Each is a
// recognised native plugin that fails cleanly from Init, mirroring how gvcfz,
// fixref's id mode and the vrfs plugin report the parts that genuinely need
// htslib internals or non-reproducible floating-point machinery.
//
//   - parental-origin: the predicted origin (paternal / maternal / uncertain)
//     is sign-robust, but its primary numeric output — the `quality` column
//     (4.3429 * |log(ppat) - log(pmat)|) and the per-site DBG probabilities
//     printed with `%e` — is a sum of pow(10,...)/log() terms and the
//     incomplete-beta tail (kf_betai); these are not bit-identical between
//     libm's pow/log/betai and Go's math, so byte parity cannot be guaranteed.
//
//   - color-chrs: writes its segmentation to a separate `<prefix>.dat` file
//     (not stdout) and decodes shared haplotype blocks with an HMM Viterbi pass
//     over a chain of up to 10,000 pre-multiplied transition matrices. The
//     repeated matrix-power floating-point accumulation is not guaranteed to be
//     bit-identical to the C HMM, and the file-output contract sits outside the
//     stdout-oriented native plugin pipeline.
//
//   - trio-dnm3: a full de-novo mutation caller (~2400 lines) over FORMAT/AD,
//     QM and SP, scoring every trio with a beta-binomial / DNG likelihood model
//     whose primary output is a phred-scaled FORMAT/DNM float annotation. That
//     output is fundamentally libm-precision-dependent (pow/log/lgamma/beta),
//     so it cannot be reproduced byte-for-byte; porting it is also out of scope
//     for the tractable-scoring criterion of this batch.
package bcftools

import (
	"fmt"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("parental-origin", func() NativePlugin { return &parentalOriginPlugin{} })
	registerNativePlugin("color-chrs", func() NativePlugin { return &colorChrsPlugin{} })
	registerNativePlugin("trio-dnm3", func() NativePlugin { return &trioDNM3Plugin{} })
}

// parentalOriginPlugin is a registration stub: its numeric output depends on
// libm pow/log and the incomplete-beta tail, which are not bit-reproducible.
type parentalOriginPlugin struct{}

// Name returns the plugin name.
func (p *parentalOriginPlugin) Name() string { return "parental-origin" }

// About returns the one-line description, matching parental-origin.c about().
func (p *parentalOriginPlugin) About() string {
	return "Determine parental origin of a CNV region in a trio.\n"
}

// RunStyle reports that parental-origin is a run()-style plugin.
func (p *parentalOriginPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of parental-origin's value-taking flags
// consumes the following CLI token, so the host can split the input file out
// even though Init then reports the mode unsupported.
func (p *parentalOriginPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-b", "--min-binom-prob", "-e", "--exclude", "-i", "--include",
		"-p", "--pfm", "-r", "--region", "-t", "--type", "-v", "--verbosity":
		return true
	}
	return false
}

// Init reports that parental-origin's numeric output is not bit-reproducible.
func (p *parentalOriginPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	return nil, fmt.Errorf("parental-origin: the quality score and DBG probabilities are sums of libm pow/log and incomplete-beta terms that are not bit-identical to Go's math, so byte parity cannot be guaranteed; run upstream bcftools for +parental-origin")
}

// Process is never reached (Init fails) but satisfies NativePlugin.
func (p *parentalOriginPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *parentalOriginPlugin) Destroy() error { return nil }

// colorChrsPlugin is a registration stub: it writes to a separate `<prefix>.dat`
// file and decodes via an HMM matrix-power Viterbi pass, neither of which fits
// the stdout-oriented, bit-reproducible native plugin pipeline.
type colorChrsPlugin struct{}

// Name returns the plugin name.
func (p *colorChrsPlugin) Name() string { return "color-chrs" }

// About returns the one-line description, matching color-chrs.c about().
func (p *colorChrsPlugin) About() string {
	return "Color shared chromosomal segments, requires phased GTs.\n"
}

// Init reports that color-chrs is unsupported by the native plugin. It is a
// generic init/process plugin (options follow `--`), so no run-style splitting
// is needed.
func (p *colorChrsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	return nil, fmt.Errorf("color-chrs: writes a separate <prefix>.dat file and decodes shared haplotype blocks with an HMM matrix-power Viterbi pass whose floating-point accumulation is not guaranteed bit-identical to the C HMM; run upstream bcftools for +color-chrs")
}

// Process is never reached (Init fails) but satisfies NativePlugin.
func (p *colorChrsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *colorChrsPlugin) Destroy() error { return nil }

// trioDNM3Plugin is a registration stub: its de-novo likelihood scoring produces
// a phred-scaled FORMAT/DNM float that is libm-precision-dependent.
type trioDNM3Plugin struct{}

// Name returns the plugin name.
func (p *trioDNM3Plugin) Name() string { return "trio-dnm3" }

// About returns the one-line description, matching trio-dnm3.c about().
func (p *trioDNM3Plugin) About() string {
	return "Screen variants for possible de-novo mutations in trios.\n"
}

// RunStyle reports that trio-dnm3 is a run()-style plugin.
func (p *trioDNM3Plugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of trio-dnm3's common value-taking flags
// consumes the following CLI token, so the host can split the input file out
// even though Init then reports the mode unsupported.
func (p *trioDNM3Plugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-p", "--pfm", "-P", "--ped", "-o", "--output", "-O", "--output-type",
		"-e", "--exclude", "-i", "--include", "-r", "--regions", "-R", "--regions-file",
		"-t", "--targets", "-T", "--targets-file", "--dnm-tag", "--use",
		"--mrate", "--pn", "--pns", "--ppl", "--min-score", "-v", "--verbosity":
		return true
	}
	return false
}

// Init reports that trio-dnm3's de-novo scoring is not bit-reproducible.
func (p *trioDNM3Plugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	return nil, fmt.Errorf("trio-dnm3: the de-novo FORMAT/DNM score is computed from a beta-binomial / DNG likelihood model (libm pow/log/lgamma/beta) that is not bit-identical to Go's math, so byte parity cannot be guaranteed; run upstream bcftools for +trio-dnm3")
}

// Process is never reached (Init fails) but satisfies NativePlugin.
func (p *trioDNM3Plugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *trioDNM3Plugin) Destroy() error { return nil }
