package bcftools

import (
	"fmt"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// applySetID rewrites the ID column of each record using the query-style
// format string. A leading '+' means "only set the ID when it is missing
// (./empty)"; without it the value is replaced unconditionally. Mirrors
// upstream vcfannotate.c:3250-3253 and the per-record block at
// vcfannotate.c:3840-3850.
func applySetID(recs []*vcf.Variant, spec string) error {
	if spec == "" {
		return nil
	}
	replace := true
	if spec[0] == '+' {
		replace = false
		spec = spec[1:]
	}
	tokens, err := ParseFormatString(spec)
	if err != nil {
		return fmt.Errorf("parse format %q: %w", spec, err)
	}
	for _, v := range recs {
		if !replace && v.ID != "" && v.ID != "." {
			continue
		}
		newID := formatSetIDTokens(tokens, v)
		// Strip a trailing '\n' if the user-supplied format happened to
		// include one — the ID column must stay on one line. Upstream's
		// convert_line also leaves the trailing newline off when used in
		// this context (vcfannotate.c:3845).
		newID = strings.TrimRight(newID, "\n")
		v.ID = newID
	}
	return nil
}

// formatSetIDTokens runs the parsed format tokens over v in the outer
// scope (no sample loop): the bcftools -I/--set-id path never expands
// `[ ... ]` groups because there is no notion of an "ID per sample".
// Sample tokens, if present, are silently skipped.
func formatSetIDTokens(tokens []FormatToken, v *vcf.Variant) string {
	var sb strings.Builder
	for _, t := range tokens {
		switch t.Kind {
		case TokenLiteral:
			sb.WriteString(t.Text)
		case TokenPlaceholder:
			sb.WriteString(formatPlaceholder(t.Text, v, -1))
		case TokenSample:
			// Upstream treats a sample-loop inside --set-id as a no-op
			// because the ID column is one-per-record; we follow suit.
		}
	}
	return sb.String()
}
