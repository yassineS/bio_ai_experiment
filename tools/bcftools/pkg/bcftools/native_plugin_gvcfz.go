// Native registration of the upstream `gvcfz` plugin (plugins/gvcfz.c).
//
// gvcfz resizes gVCF blocks according to -g/--group-by expressions such as
// `PASS:GQ>60 & DP<20; Flt1:QG>20`. Each group is a full bcftools FORMAT/GT
// filter expression evaluated per record (filter_init/filter_test in htslib),
// covering FORMAT fields (GQ, DP, MIN_DP), genotype predicates (GT="alt"), and
// boolean combinators. The native plugin framework does not expose that filter
// engine (the in-tree filter only evaluates INFO/site-level fields, not
// FORMAT/GT predicates), so producing the block grouping would require silently
// approximating the expression semantics. Rather than emit subtly wrong blocks,
// gvcfz is registered as a recognised native plugin that fails cleanly from
// Init, mirroring how fixref's id mode and the vrfs plugin report the parts that
// genuinely need htslib internals.
package bcftools

import (
	"fmt"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("gvcfz", func() NativePlugin { return &gvcfzPlugin{} }) }

// gvcfzPlugin is a registration stub for the gvcfz plugin. Its grouping
// requires the htslib FORMAT/GT filter expression engine, which the native
// pipeline does not provide.
type gvcfzPlugin struct{}

// Name returns the plugin name.
func (p *gvcfzPlugin) Name() string { return "gvcfz" }

// About returns the one-line description, matching gvcfz.c about().
func (p *gvcfzPlugin) About() string {
	return "Compress gVCF file by resizing gVCF blocks according to specified criteria."
}

// RunStyle reports that gvcfz is a run()-style plugin (it exports a `run`
// symbol), so the host parses its options in the run-style form.
func (p *gvcfzPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of gvcfz's own value-taking flags consumes
// the following CLI token, so the host can split the input file from the
// options even though Init then reports the mode unsupported.
func (p *gvcfzPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-g", "--group-by", "-i", "--include", "-e", "--exclude",
		"-o", "--output", "-O", "--output-type", "-v", "--verbosity":
		return true
	}
	return false
}

// Init reports that gvcfz's block grouping is unsupported by the native plugin
// because it needs the htslib FORMAT/GT filter expression engine.
func (p *gvcfzPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	return nil, fmt.Errorf("gvcfz: block grouping requires the htslib FORMAT/GT filter expression engine, which the native plugin does not provide; run upstream bcftools for +gvcfz")
}

// Process is never reached (Init fails) but satisfies NativePlugin.
func (p *gvcfzPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *gvcfzPlugin) Destroy() error { return nil }
