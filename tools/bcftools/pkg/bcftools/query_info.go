package bcftools

import (
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// formatWholeInfo re-serialises the entire INFO column from a Variant,
// preserving the original source order via InfoOrder. Flag tags (empty
// value) are emitted as the bare key. An empty INFO column renders as ".".
// Mirrors upstream convert.c::process_info with fmt->key == NULL.
//
// Implements the bare `%INFO` token of `bcftools query -f`. The lookup
// path is: source-order via v.InfoOrder when present (the parity-correct
// branch for VCF inputs), else stable alphabetical for hand-built Variants
// so Go's randomised map iteration cannot produce nondeterministic output.
func formatWholeInfo(v *vcf.Variant) string {
	if len(v.Info) == 0 {
		return "."
	}
	order := v.InfoOrder
	if len(order) == 0 {
		keys := make([]string, 0, len(v.Info))
		for k := range v.Info {
			keys = append(keys, k)
		}
		sort.Strings(keys)
		order = keys
	}
	var sb strings.Builder
	first := true
	for _, k := range order {
		val, ok := v.Info[k]
		if !ok {
			continue
		}
		if !first {
			sb.WriteByte(';')
		}
		first = false
		sb.WriteString(k)
		if val != "" {
			sb.WriteByte('=')
			sb.WriteString(val)
		}
	}
	if first {
		return "."
	}
	return sb.String()
}
