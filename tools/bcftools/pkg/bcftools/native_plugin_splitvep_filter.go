// split-vep -i/-e filter wiring: the -i/-e expression evaluates against the
// expanded per-transcript CSQ subfields, which upstream registers as synthetic
// INFO tags on the output header before compiling the filter (filter_init on
// args->hdr_out, after parse_filter_str/parse_column_str). This file reproduces
// that: it auto-registers any CSQ subfield referenced by the expression as an
// extra column (parse_filter_str's "add the undefined tags to the -c string"
// step) and compiles the filter against the augmented header, reusing the shared
// filter engine. See native_plugin_splitvep_proc.go for where the compiled
// filter is applied per record (filter_and_output).
package bcftools

import (
	"fmt"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// buildFilter compiles the -i/-e expression (if any) against the split-vep
// output header so bare identifiers resolve to the registered CSQ-subfield INFO
// tags, mirroring upstream filter_init(args->hdr_out). A nil result means "no
// filter" (every record passes).
func (p *splitVepPlugin) buildFilter(outHdr *vcf.Header) (*pluginFilter, error) {
	if !p.filterSet || p.filterStr == "" {
		return nil, nil
	}
	f, err := newPluginFilterWithHeader(p.filterStr, p.filterExclude, outHdr)
	if err != nil {
		return nil, fmt.Errorf("split-vep: -i/-e expression: %w", err)
	}
	return f, nil
}

// addFilterColumns scans the -i/-e expression for CSQ subfield names and adds
// any that are not already requested columns to p.annots, so they are registered
// on the output header and populated as INFO tags before the filter runs. This
// mirrors parse_filter_str, which appends the expression's undefined-but-valid
// subfields to the --columns string. A referenced identifier that is neither a
// builtin column, an existing INFO tag in the input header, nor a CSQ subfield
// is reported as an error, matching upstream's
// "the tag ... is not defined in the VCF header or in INFO/CSQ".
func (p *splitVepPlugin) addFilterColumns() error {
	if !p.filterSet || p.filterStr == "" {
		return nil
	}
	have := map[string]bool{}
	for _, a := range p.annots {
		have[a.tag] = true
	}
	for _, tok := range filterIdentifiers(p.filterStr) {
		// The raw CSQ/BCSQ tag itself is handled separately (raw_vep_request);
		// it is always defined in the input header, so it is not auto-added.
		if tok == p.vepTag {
			continue
		}
		idx, field, ok := p.lookupField(tok)
		if !ok {
			// Not a CSQ subfield. If it is a builtin column or already declared in
			// the input header it is fine (the filter resolves it directly);
			// otherwise upstream errors out.
			if isBuiltinColumn(tok) || infoHeaderExists(p.inHdr, tok) {
				continue
			}
			return fmt.Errorf("Error: the tag \"%s\" is not defined in the VCF header or in INFO/%s", tok, p.vepTag)
		}
		if have[field] {
			continue
		}
		have[field] = true
		p.annots = append(p.annots, svAnnot{
			field: field,
			idx:   idx,
			tag:   field,
			typ:   p.defaultColumnType(field),
		})
	}
	return nil
}

// lookupField resolves an identifier to a CSQ subfield index and its (prefixed,
// sanitized) field name, trying the verbatim name first and then the sanitized
// form, exactly as parse_column_str does.
func (p *splitVepPlugin) lookupField(name string) (int, string, bool) {
	if idx, ok := p.field2idx[name]; ok {
		return idx, p.fields[idx], true
	}
	if sname := p.sanitizeField(strings.TrimPrefix(name, p.annotPrefix)); sname != name {
		if idx, ok := p.field2idx[sname]; ok {
			return idx, p.fields[idx], true
		}
	}
	return 0, "", false
}

// filterIdentifiers extracts the bare identifier tokens from a filter
// expression, skipping quoted string literals (single or double quoted) and the
// INFO/, FMT/ and FORMAT/ prefixes so only the resolvable tag names are
// returned. Numbers and operators are ignored. Identifiers are returned in order
// of first appearance, without duplicates.
func filterIdentifiers(expr string) []string {
	var out []string
	seen := map[string]bool{}
	i := 0
	for i < len(expr) {
		c := expr[i]
		switch {
		case c == '"' || c == '\'':
			// Skip the quoted literal.
			q := c
			i++
			for i < len(expr) && expr[i] != q {
				i++
			}
			if i < len(expr) {
				i++ // closing quote
			}
		case isIdentStart(c):
			start := i
			for i < len(expr) && (isIdentPart(expr[i]) || expr[i] == '/') {
				i++
			}
			tok := expr[start:i]
			// Strip a leading namespace prefix; a bare INFO/X means subfield X.
			for _, pfx := range []string{"INFO/", "FMT/", "FORMAT/"} {
				if strings.HasPrefix(tok, pfx) {
					tok = tok[len(pfx):]
					break
				}
			}
			if tok == "" {
				continue
			}
			// Skip a leading digit token (e.g. part of a number caught by the
			// identifier scan is impossible here, but keywords like true/false are
			// harmless to skip).
			switch tok {
			case "true", "false":
				continue
			}
			if !seen[tok] {
				seen[tok] = true
				out = append(out, tok)
			}
		default:
			i++
		}
	}
	return out
}
