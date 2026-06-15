// Native port of the upstream `add-variantkey` plugin
// (plugins/add-variantkey.c). It adds two INFO/String annotations to every
// record:
//
//   - VKX: the 16-character hexadecimal 64-bit VariantKey computed from CHROM,
//     0-based POS, REF and the FIRST ALT allele (upstream uses rec->d.allele[1]).
//   - RSX: the 8-character hexadecimal rendering of the variant ID with its
//     first two characters dropped (the "rs" prefix) and parsed as a uint32.
//
// It is a generic per-record, parallel plugin: each record is annotated
// independently from its own CHROM/POS/REF/ALT and ID.
package bcftools

import (
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("add-variantkey", func() NativePlugin { return &addVariantKeyPlugin{} })
}

// addVariantKeyPlugin implements add-variantkey.
type addVariantKeyPlugin struct{}

// Name returns the plugin name.
func (p *addVariantKeyPlugin) Name() string { return "add-variantkey" }

// About returns the one-line description, matching add-variantkey.c about().
func (p *addVariantKeyPlugin) About() string {
	return "Add VariantKey INFO fields VKX and RSX.\n"
}

// Parallel reports true: each record is annotated independently.
func (p *addVariantKeyPlugin) Parallel() bool { return true }

// Init appends the VKX and RSX INFO header lines, matching upstream's two
// bcf_hdr_append calls in order.
func (p *addVariantKeyPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	out.MetaInfo = appendInfoHeader(out.MetaInfo, `##INFO=<ID=VKX,Number=1,Type=String,Description="Hexadecimal representation of 64 bit VariantKey">`)
	out.MetaInfo = appendInfoHeader(out.MetaInfo, `##INFO=<ID=RSX,Number=1,Type=String,Description="Hexadecimal representation of ID minus the 'rs' prefix (32bit)">`)
	return out, nil
}

// Process annotates a single record with VKX and RSX.
func (p *addVariantKeyPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	alt := ""
	if len(v.Alt) > 0 {
		alt = v.Alt[0]
	}
	// Upstream passes rec->pos (0-based). v.Pos is 1-based VCF coordinate.
	vk := variantKey(v.Chrom, uint32(v.Pos-1), v.Ref, alt)
	setInfo(v, "VKX", variantKeyHex(vk))
	setInfo(v, "RSX", formatRSX(v.ID))
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *addVariantKeyPlugin) Destroy() error { return nil }

// formatRSX reproduces upstream's RSX computation: take the ID, advance the
// pointer past the first two characters (the "rs" prefix), parse the leading
// decimal digits as a uint32 (C strtoul base 10, which stops at the first
// non-digit and yields 0 when none are present), and render it as %08x. When
// the ID is shorter than two characters (e.g. "."), there is nothing left to
// parse and the result is "00000000", matching the C behaviour for those IDs.
func formatRSX(id string) string {
	var rs uint32
	if len(id) > 2 {
		s := id[2:]
		i := 0
		for i < len(s) && s[i] >= '0' && s[i] <= '9' {
			rs = rs*10 + uint32(s[i]-'0')
			i++
		}
	}
	return rsHex(rs)
}

// rsHex renders a uint32 as an 8-character lowercase hexadecimal string,
// matching the C "%08" PRIx32 formatting.
func rsHex(rs uint32) string {
	var buf [8]byte
	for i := 7; i >= 0; i-- {
		buf[i] = vkHexDigits[rs&0xF]
		rs >>= 4
	}
	return string(buf[:])
}
