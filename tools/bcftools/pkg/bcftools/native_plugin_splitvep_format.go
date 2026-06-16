// split-vep -f/--format query engine: a focused re-implementation of the
// `bcftools query`-style format string covering the directives split-vep
// commonly uses. See native_plugin_splitvep.go for the lifecycle.
package bcftools

import (
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// svFmtKind enumerates the kinds of format tokens.
type svFmtKind int

const (
	svFmtLiteral svFmtKind = iota
	svFmtChrom
	svFmtPos
	svFmtID
	svFmtRef
	svFmtAlt
	svFmtQual
	svFmtFilter
	svFmtInfo // a plain INFO tag (%INFO/X or %X that is an INFO key)
	svFmtCsqField
)

// svFmtItem is one parsed token of the format string.
type svFmtItem struct {
	kind    svFmtKind
	literal string // for literal tokens
	name    string // INFO tag name
	csqIdx  int    // for svFmtCsqField: subfield index
	csqType int    // for svFmtCsqField: value type for re-rendering
}

// resolveFormat parses the (possibly %CSQ-expanded) format string into tokens.
func (p *splitVepPlugin) resolveFormat() error {
	format := p.formatStr
	if p.allFields != "" {
		expanded, err := p.expandCsqExpression(format)
		if err != nil {
			return err
		}
		format = expanded
	}
	items, err := p.parseFormatTokens(format)
	if err != nil {
		return err
	}
	p.formatItems = items
	// Register the CSQ subfields referenced by the format string as columns, in
	// subfield (header) order, mirroring parse_format_str's append-to-column_str.
	// In text mode the INFO tags are not printed, but they make filterAndOutput's
	// "nannot" branch (and any -i/-e filter) see the same set upstream does.
	referenced := make([]bool, len(p.fields))
	for _, it := range items {
		if it.kind == svFmtCsqField && it.csqIdx >= 0 && it.csqIdx < len(referenced) {
			referenced[it.csqIdx] = true
		}
	}
	seen := map[string]bool{}
	for idx, ref := range referenced {
		if !ref {
			continue
		}
		field := p.fields[idx]
		if seen[field] {
			continue
		}
		seen[field] = true
		p.annots = append(p.annots, svAnnot{
			field: field,
			idx:   idx,
			tag:   field,
			typ:   p.defaultColumnType(field),
		})
	}
	return nil
}

// expandCsqExpression replaces a standalone %CSQ (i.e. %<vepTag>) token with the
// delimiter-joined list of every subfield, matching expand_csq_expression. The
// delimiter is "tab"->"\t", "space"->" ", else the literal string.
func (p *splitVepPlugin) expandCsqExpression(format string) (string, error) {
	delim := p.allFields
	switch delim {
	case "tab":
		delim = "\t"
	case "space":
		delim = " "
	}
	repl := make([]string, len(p.fields))
	for i, f := range p.fields {
		repl[i] = "%" + f
	}
	expansion := strings.Join(repl, delim)
	token := "%" + p.vepTag
	var b strings.Builder
	i := 0
	for i < len(format) {
		if strings.HasPrefix(format[i:], token) {
			after := i + len(token)
			if after >= len(format) || !isTagChar(format[after]) {
				b.WriteString(expansion)
				i = after
				continue
			}
		}
		b.WriteByte(format[i])
		i++
	}
	return b.String(), nil
}

// parseFormatTokens tokenizes a format string into literals and %DIRECTIVE
// tokens, resolving \t and \n escapes and rejecting unsupported directives.
func (p *splitVepPlugin) parseFormatTokens(format string) ([]svFmtItem, error) {
	var items []svFmtItem
	var lit strings.Builder
	flush := func() {
		if lit.Len() > 0 {
			items = append(items, svFmtItem{kind: svFmtLiteral, literal: lit.String()})
			lit.Reset()
		}
	}
	i := 0
	for i < len(format) {
		c := format[i]
		switch c {
		case '\\':
			if i+1 < len(format) {
				switch format[i+1] {
				case 'n':
					lit.WriteByte('\n')
				case 't':
					lit.WriteByte('\t')
				case '\\':
					lit.WriteByte('\\')
				default:
					lit.WriteByte(format[i+1])
				}
				i += 2
				continue
			}
			lit.WriteByte('\\')
			i++
		case '%':
			// Read the directive name (letters, digits, _, ., and an optional
			// INFO/ prefix).
			j := i + 1
			infoPrefix := false
			if strings.HasPrefix(format[j:], "INFO/") {
				infoPrefix = true
				j += len("INFO/")
			}
			start := j
			for j < len(format) && isTagChar(format[j]) {
				j++
			}
			name := format[start:j]
			if name == "" {
				lit.WriteByte('%')
				i++
				continue
			}
			flush()
			item, err := p.resolveFormatDirective(name, infoPrefix)
			if err != nil {
				return nil, err
			}
			items = append(items, item)
			i = j
		default:
			lit.WriteByte(c)
			i++
		}
	}
	flush()
	return items, nil
}

// isTagChar reports whether c may appear in a tag/field name.
func isTagChar(c byte) bool {
	return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' || c == '.'
}

// resolveFormatDirective maps a %NAME directive to a format item. A CSQ subfield
// name takes precedence over a same-named INFO tag (matching upstream, which
// warns and prefers the subfield), unless the explicit INFO/ prefix is used.
func (p *splitVepPlugin) resolveFormatDirective(name string, infoPrefix bool) (svFmtItem, error) {
	switch name {
	case "CHROM":
		return svFmtItem{kind: svFmtChrom}, nil
	case "POS":
		return svFmtItem{kind: svFmtPos}, nil
	case "ID":
		return svFmtItem{kind: svFmtID}, nil
	case "REF":
		return svFmtItem{kind: svFmtRef}, nil
	case "ALT":
		return svFmtItem{kind: svFmtAlt}, nil
	case "QUAL":
		return svFmtItem{kind: svFmtQual}, nil
	case "FILTER":
		return svFmtItem{kind: svFmtFilter}, nil
	}
	if !infoPrefix {
		if idx, ok := p.field2idx[name]; ok {
			return svFmtItem{kind: svFmtCsqField, csqIdx: idx, csqType: p.defaultColumnType(p.fields[idx])}, nil
		}
		if sname := p.sanitizeField(strings.TrimPrefix(name, p.annotPrefix)); sname != name {
			if idx, ok := p.field2idx[sname]; ok {
				return svFmtItem{kind: svFmtCsqField, csqIdx: idx, csqType: p.defaultColumnType(p.fields[idx])}, nil
			}
		}
	}
	if name == p.vepTag {
		return svFmtItem{}, fmt.Errorf("split-vep: the bare %%%s directive is only supported with -A/--all-fields", name)
	}
	// Treat as a plain INFO tag.
	return svFmtItem{kind: svFmtInfo, name: name}, nil
}

// runFormat implements the text (-f) output mode. For each record it splits the
// CSQ, selects and severity-filters transcripts, and emits one line per
// transcript (-d) or one collapsed line per record. Records with no passing
// consequence are dropped when drop_sites is set. When -i/-e is given the
// expression's CSQ subfields are populated as INFO tags and the filter gates
// each emitted line, exactly where upstream's filter_and_output evaluates it for
// the text writer.
func (p *splitVepPlugin) runFormat(variants []*vcf.Variant, out io.Writer) error {
	filter, err := p.buildFilter(p.annotateHeader(p.inHdr))
	if err != nil {
		return err
	}
	var b strings.Builder
	emit := func(rec *vcf.Variant, current [][]string) error {
		// drop_sites==0 with no passing consequence renders a bare line; upstream
		// reaches writeFormatLine with an empty transcript set in that case.
		p.writeFormatLine(&b, rec, current)
		return nil
	}
	for _, v := range variants {
		if err := p.processVariant(v, filter, emit); err != nil {
			return err
		}
	}
	_, err = io.WriteString(out, b.String())
	return err
}

// writeFormatLine renders the format string once for the given set of
// transcripts. CSQ-subfield tokens are comma-joined across the transcripts;
// site-level tokens are constant. The convert engine always forces a trailing
// newline, which the format string itself usually supplies.
func (p *splitVepPlugin) writeFormatLine(b *strings.Builder, v *vcf.Variant, transcripts [][]string) {
	for _, item := range p.formatItems {
		switch item.kind {
		case svFmtLiteral:
			b.WriteString(item.literal)
		case svFmtChrom:
			b.WriteString(v.Chrom)
		case svFmtPos:
			fmt.Fprintf(b, "%d", v.Pos)
		case svFmtID:
			b.WriteString(v.ID)
		case svFmtRef:
			b.WriteString(v.Ref)
		case svFmtAlt:
			b.WriteString(strings.Join(v.Alt, ","))
		case svFmtQual:
			b.WriteString(formatQual(v.Qual))
		case svFmtFilter:
			if len(v.Filter) == 0 {
				b.WriteByte('.')
			} else {
				b.WriteString(strings.Join(v.Filter, ";"))
			}
		case svFmtInfo:
			if val, ok := v.Info[item.name]; ok {
				b.WriteString(val)
			} else {
				b.WriteByte('.')
			}
		case svFmtCsqField:
			vals := make([]string, 0, len(transcripts))
			for _, tr := range transcripts {
				raw := trField(tr, item.csqIdx)
				if raw == "" {
					vals = append(vals, ".")
				} else {
					vals = append(vals, p.renderTyped(item.csqType, raw))
				}
			}
			if len(vals) == 0 {
				b.WriteByte('.')
			} else {
				b.WriteString(strings.Join(vals, ","))
			}
		}
	}
}

// formatQual renders the QUAL value as bcftools does: "." for missing, else the
// %g-formatted number.
func formatQual(q float64) string {
	if q < 0 {
		return "."
	}
	return formatVCFFloat(q)
}
