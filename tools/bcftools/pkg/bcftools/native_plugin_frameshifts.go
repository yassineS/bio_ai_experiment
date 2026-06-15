// Native registration for the upstream `frameshifts` plugin
// (plugins/frameshifts.c). The plugin annotates frameshift indels (INFO/OOF)
// using an exons BED via htslib's bcf_sr_regions_t cursor: after each
// bcf_sr_regions_overlap call it reads the matched exon's reg->start/reg->end
// and trims the inserted/deleted length against that exon to decide in-frame
// vs out-of-frame. Reproducing this byte-for-byte requires porting the
// synced-reader BED region cursor (seek, prev_seq/prev_start advancement, the
// is_bed from++ adjustment and the start=from-1/end=to-1 storage) whose state
// is consumed indirectly through reg->start/reg->end — machinery the native
// plugin framework does not yet expose. Rather than silently diverge, the
// plugin is registered but reports a clean "unsupported" error from Init.
package bcftools

import (
	"fmt"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("frameshifts", func() NativePlugin { return &frameshiftsPlugin{} }) }

// frameshiftsPlugin is the unsupported native entry for `frameshifts`.
type frameshiftsPlugin struct{}

// Name returns the plugin name.
func (p *frameshiftsPlugin) Name() string { return "frameshifts" }

// About returns the one-line description, matching frameshifts.c about().
func (p *frameshiftsPlugin) About() string { return "Annotate frameshift indels." }

// Init reports that the native port is unsupported and why.
func (p *frameshiftsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	return nil, fmt.Errorf("frameshifts: the native plugin is not supported; it requires htslib's bcf_sr_regions_t BED cursor (exons->start/exons->end after bcf_sr_regions_overlap) to trim indel lengths against exons, which is not available in the native plugin framework")
}

// Process is never reached (Init fails first).
func (p *frameshiftsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *frameshiftsPlugin) Destroy() error { return nil }
