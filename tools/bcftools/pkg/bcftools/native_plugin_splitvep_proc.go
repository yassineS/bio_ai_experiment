// split-vep record processing: severity scale, the query format engine, the
// annotate engine, and the shared per-record CSQ split. See
// native_plugin_splitvep.go for the plugin lifecycle and option parsing.
package bcftools

import (
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// svDefaultScale is the upstream default consequence severity scale, ascending
// (tier 0 = least severe). Tokens on the same line share a tier; they are
// matched as case-insensitive substrings of a consequence term. This mirrors
// default_severity() in split-vep.c byte-for-byte.
var svDefaultScale = [][]string{
	{"intergenic"},
	{"feature_truncation", "feature_elongation"},
	{"regulatory"},
	{"tf_binding_site", "tfbs"},
	{"downstream", "upstream"},
	{"non_coding_transcript", "non_coding"},
	{"intron", "nmd_transcript"},
	{"non_coding_transcript_exon"},
	{"5_prime_utr", "3_prime_utr"},
	{"coding_sequence", "mature_mirna"},
	{"stop_retained", "start_retained", "synonymous"},
	{"incomplete_terminal_codon"},
	{"splice_region"},
	{"missense", "inframe", "protein_altering"},
	{"transcript_amplification"},
	{"exon_loss"},
	{"disruptive"},
	{"start_lost", "stop_lost", "stop_gained", "frameshift"},
	{"splice_acceptor", "splice_donor"},
	{"transcript_ablation"},
}

// initSeverityScale loads the default scale into the substring list and the
// token->tier map, mirroring init_severity().
func (p *splitVepPlugin) initSeverityScale() {
	p.csq2sev = map[string]int{}
	p.scale = nil
	for tier, tokens := range svDefaultScale {
		for _, tok := range tokens {
			p.scale = append(p.scale, tok)
			if _, ok := p.csq2sev[tok]; !ok {
				p.csq2sev[tok] = tier
			}
		}
	}
}

// csqToSeverity returns the (min,max) severity tiers spanned by a consequence
// string (the '&'-joined consequence terms), porting csq_to_severity. New
// tokens are resolved by substring against the scale (else assigned a high tier)
// and cached. When exact >= 0 it returns early on the first token whose tier
// equals exact.
func (p *splitVepPlugin) csqToSeverity(csq string, exact int) (int, int) {
	minSev := int(^uint(0) >> 1)
	maxSev := -1
	for _, tokRaw := range strings.Split(csq, "&") {
		tok := strings.ToLower(tokRaw)
		sev, ok := p.csq2sev[tok]
		if !ok {
			idx := -1
			for i, s := range p.scale {
				if strings.Contains(tok, s) {
					idx = i
					break
				}
			}
			if idx >= 0 {
				sev = p.csq2sev[p.scale[idx]]
			} else {
				sev = len(p.scale) + 1
			}
			p.scale = append(p.scale, tok)
			p.csq2sev[tok] = sev
		}
		if exact < 0 {
			if minSev > sev {
				minSev = sev
			}
			if maxSev < sev {
				maxSev = sev
			}
		} else if sev == exact {
			return sev, sev
		}
	}
	return minSev, maxSev
}

// csqSeverityPass reports whether a consequence string passes the configured
// severity range, porting csq_severity_pass.
func (p *splitVepPlugin) csqSeverityPass(csq string) bool {
	if p.minSev == p.maxSev && p.minSev == svSelectAny {
		return true
	}
	exact := -1
	if p.minSev == p.maxSev {
		exact = p.minSev
	}
	minSev, maxSev := p.csqToSeverity(csq, exact)
	if maxSev < p.minSev {
		return false
	}
	if minSev > p.maxSev {
		return false
	}
	return true
}

// worstTranscript returns the index of the transcript with the greatest maximum
// severity, porting get_worst_transcript (ties keep the first).
func (p *splitVepPlugin) worstTranscript(transcripts [][]string) int {
	maxSeverity := -1
	imax := 0
	for i, tr := range transcripts {
		if p.csqIdx >= len(tr) {
			continue
		}
		_, max := p.csqToSeverity(tr[p.csqIdx], -1)
		if maxSeverity < max {
			imax = i
			maxSeverity = max
		}
	}
	return imax
}

// splitCsq splits a record's CSQ INFO value into per-transcript subfield slices
// (comma-separated transcripts, each pipe-separated). It returns nil when the
// tag is absent.
func (p *splitVepPlugin) splitCsq(v *vcf.Variant) [][]string {
	val, ok := v.Info[p.vepTag]
	if !ok || val == "" {
		return nil
	}
	var out [][]string
	for _, tr := range strings.Split(val, ",") {
		out = append(out, strings.Split(tr, "|"))
	}
	return out
}

// selectTranscripts returns the indices of the transcripts selected by the
// transcript-selection mode.
func (p *splitVepPlugin) selectTranscripts(transcripts [][]string) []int {
	if len(transcripts) == 0 {
		return nil
	}
	switch p.selTr {
	case svTrWorst:
		return []int{p.worstTranscript(transcripts)}
	default: // svTrAll
		idx := make([]int, len(transcripts))
		for i := range transcripts {
			idx[i] = i
		}
		return idx
	}
}

// trField returns the value of subfield idx in transcript tr, or "" when the
// transcript has too few columns.
func trField(tr []string, idx int) string {
	if idx < 0 || idx >= len(tr) {
		return ""
	}
	return tr[idx]
}

// runAnnotate implements the default (VCF/BCF) output mode: each record's CSQ is
// split, transcripts are selected and severity-filtered, and the requested
// columns are written as INFO tags. Records are then re-emitted in the requested
// container.
func (p *splitVepPlugin) runAnnotate(opts PluginOptions, hdr *vcf.Header, variants []*vcf.Variant, out io.Writer) error {
	outHdr := p.annotateHeader(hdr)

	// Compile the -i/-e expression against the augmented output header so it can
	// reference the per-transcript CSQ subfield columns that annotateHeader just
	// registered as synthetic INFO tags, mirroring upstream split-vep.c, which
	// calls filter_init(args->hdr_out) after registering those columns.
	var filter *pluginFilter
	if p.filterStr != "" {
		f, ferr := newPluginFilterWithHeader(p.filterStr, p.filterExclude, outHdr)
		if ferr != nil {
			return fmt.Errorf("split-vep: -i/-e expression: %w", ferr)
		}
		filter = f
	}

	w, cleanup, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
		Threads:       opts.Threads,
	}, outHdr)
	if err != nil {
		return err
	}
	if err := w.WriteHeader(); err != nil {
		cleanup()
		return err
	}
	for _, v := range variants {
		emit := p.annotateRecord(v)
		if !emit {
			continue
		}
		// Apply the -i/-e filter after the derived columns have been written as
		// INFO tags, exactly where upstream's filter_and_output evaluates it.
		if !filter.testSite(v) {
			continue
		}
		if err := w.Write(v); err != nil {
			cleanup()
			return err
		}
	}
	if err := w.Flush(); err != nil {
		cleanup()
		return err
	}
	cleanup()
	return nil
}

// annotateHeader appends the INFO header lines for the requested columns, not
// overwriting existing definitions, matching upstream.
func (p *splitVepPlugin) annotateHeader(hdr *vcf.Header) *vcf.Header {
	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	for _, a := range p.annots {
		typ := "String"
		switch a.typ {
		case svTypeInt:
			typ = "Integer"
		case svTypeReal:
			typ = "Float"
		}
		line := fmt.Sprintf(`##INFO=<ID=%s,Number=.,Type=%s,Description="The %s field from INFO/%s">`,
			a.tag, typ, a.field, p.vepTag)
		out.MetaInfo = appendInfoHeader(out.MetaInfo, line)
	}
	return out
}

// annotateRecord splits and annotates a single record in place. It returns
// whether the record should be emitted (drop_sites controls dropping of records
// with no passing consequence).
func (p *splitVepPlugin) annotateRecord(v *vcf.Variant) bool {
	transcripts := p.splitCsq(v)
	if len(transcripts) == 0 {
		// No CSQ: pass through unannotated unless drop_sites.
		return p.dropSites == 0
	}
	selected := p.selectTranscripts(transcripts)
	values := make([][]string, len(p.annots))
	severityPass := false
	for _, ti := range selected {
		tr := transcripts[ti]
		if p.csqIdx >= 0 && p.csqIdx < len(tr) {
			if !p.csqSeverityPass(tr[p.csqIdx]) {
				continue
			}
		} else if p.minSev != svSelectAny {
			continue
		}
		severityPass = true
		for ai, a := range p.annots {
			var val string
			if a.idx == -1 {
				val = strings.Join(tr, "|")
			} else {
				val = trField(tr, a.idx)
			}
			values[ai] = append(values[ai], val)
		}
	}
	updated := false
	for ai, a := range p.annots {
		if len(values[ai]) == 0 {
			continue
		}
		joined := strings.Join(values[ai], ",")
		out := p.renderTyped(a.typ, joined)
		setInfo(v, a.tag, out)
		updated = true
	}
	if p.dropSites == 1 && (!updated || !severityPass) {
		return false
	}
	return true
}

// renderTyped reformats a comma-separated value list according to the column
// type, matching upstream's parse_array_real / parse_array_int32 (unparseable
// tokens become ".").
func (p *splitVepPlugin) renderTyped(typ int, s string) string {
	if typ == svTypeStr {
		return s
	}
	parts := strings.Split(s, ",")
	out := make([]string, len(parts))
	for i, tok := range parts {
		switch typ {
		case svTypeReal:
			f, err := strconv.ParseFloat(tok, 64)
			if err != nil {
				out[i] = "."
			} else {
				out[i] = formatVCFFloat(f)
			}
		case svTypeInt:
			n, err := strconv.ParseInt(leadingInt(tok), 10, 64)
			if err != nil {
				out[i] = "."
			} else {
				out[i] = strconv.FormatInt(n, 10)
			}
		default:
			out[i] = tok
		}
	}
	return strings.Join(out, ",")
}

// leadingInt returns the leading integer portion of s (matching C strtol, which
// stops at the first non-digit). An empty result signals a parse failure.
func leadingInt(s string) string {
	i := 0
	if i < len(s) && (s[i] == '+' || s[i] == '-') {
		i++
	}
	start := i
	for i < len(s) && s[i] >= '0' && s[i] <= '9' {
		i++
	}
	if i == start {
		return ""
	}
	return s[:i]
}

// infoHeaderExists reports whether an ##INFO line for id is present.
func infoHeaderExists(hdr *vcf.Header, id string) bool {
	for _, m := range hdr.MetaInfo {
		if strings.HasPrefix(m, "##INFO=") && headerID(m) == id {
			return true
		}
	}
	return false
}

// infoHeaderDescription returns the Description="..." value of the ##INFO line
// for id (with surrounding quotes retained on the value as upstream sees them),
// or "" if not found.
func infoHeaderDescription(hdr *vcf.Header, id string) string {
	for _, m := range hdr.MetaInfo {
		if !strings.HasPrefix(m, "##INFO=") || headerID(m) != id {
			continue
		}
		i := strings.Index(m, "Description=")
		if i < 0 {
			return ""
		}
		rest := m[i+len("Description="):]
		rest = strings.TrimPrefix(rest, `"`)
		// Trim the trailing '">' of the structured line.
		rest = strings.TrimSuffix(rest, ">")
		rest = strings.TrimSuffix(rest, `"`)
		return rest
	}
	return ""
}
