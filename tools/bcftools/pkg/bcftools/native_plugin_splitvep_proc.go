// split-vep record processing: severity scale, the query format engine, the
// annotate engine, and the shared per-record CSQ split. See
// native_plugin_splitvep.go for the plugin lifecycle and option parsing.
package bcftools

import (
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// initSeverityScale loads the consequence severity scale into the substring list
// (p.scale) and the token->tier map (p.csq2sev), mirroring the severity-scale
// block of init_data. The source is the -S FILE when given, otherwise the
// built-in default_severity() text (svDefaultSeverityText). Each line is a tier
// (tokens on the same line, in arbitrary order, share that tier); '#' comment
// lines are skipped; tokens are lowercased. The walk is a faithful port of
// upstream's token loop, where the tier advances on the '\n' that terminates a
// token, so the tier counts newlines rather than lines.
func (p *splitVepPlugin) initSeverityScale() error {
	text := svDefaultSeverityText
	if p.severity != "" && p.severity != "-" && p.severity != "?" {
		data, err := os.ReadFile(p.severity)
		if err != nil {
			return fmt.Errorf("Cannot read %s", p.severity)
		}
		// Upstream re-joins each read line with a trailing '\n', so a file lacking a
		// final newline still gets one. Normalise to that form.
		text = string(data)
		if text != "" && !strings.HasSuffix(text, "\n") {
			text += "\n"
		}
	}
	p.csq2sev = map[string]int{}
	p.scale = nil
	severity := 0
	i := 0
	for i < len(text) {
		c := text[i]
		if c == '#' {
			for i < len(text) && text[i] != '\n' {
				i++
			}
			if i < len(text) {
				i++ // consume the '\n'
			}
			continue
		}
		// Read one token (a maximal run of non-space characters), lowercased.
		start := i
		for i < len(text) && !isSpaceByte(text[i]) {
			i++
		}
		tok := strings.ToLower(text[start:i])
		p.scale = append(p.scale, tok)
		if _, ok := p.csq2sev[tok]; !ok {
			p.csq2sev[tok] = severity
		}
		if i >= len(text) {
			break
		}
		if text[i] == '\n' {
			severity++
		}
		i++
		// Skip any whitespace before the next token (spaces between same-tier
		// tokens, and blank lines, which do not further bump the tier).
		for i < len(text) && isSpaceByte(text[i]) {
			i++
		}
	}
	return nil
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

// csqRewriteWorst returns the single most-severe term of a '&'-joined
// consequence string, porting csq_rewrite_worst. The per-term severity is an
// EXACT lookup into csq2severity (no substring fallback and no lazy insertion,
// unlike csqToSeverity); a term absent from the scale therefore scores -1, and on
// a tie the first term wins. When the string has only one term it is returned
// unchanged.
func (p *splitVepPlugin) csqRewriteWorst(s string) string {
	terms := strings.Split(s, "&")
	if len(terms) <= 1 {
		return s
	}
	imax, smax := 0, -1
	for i, t := range terms {
		sev := -1
		if v, ok := p.csq2sev[t]; ok {
			sev = v
		}
		if smax < sev {
			smax = sev
			imax = i
		}
	}
	return terms[imax]
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

// selectTranscripts returns the indices (into transcripts) selected by the
// transcript-selection mode, porting the select-transcripts block of
// process_record. n is the retained transcript count after -g restriction (so
// only transcripts[:n] are considered). The EXPRESSION mode returns the matching
// transcripts; worst returns the single most-severe; all returns every retained
// transcript.
func (p *splitVepPlugin) selectTranscripts(transcripts [][]string, n int) ([]int, error) {
	if n <= 0 {
		return nil, nil
	}
	retained := transcripts[:n]
	switch p.selTr {
	case svTrExpr:
		return p.matchingTranscripts(retained)
	case svTrWorst:
		return []int{p.worstTranscript(retained)}, nil
	default: // svTrAll
		idx := make([]int, n)
		for i := range retained {
			idx[i] = i
		}
		return idx, nil
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
// container. With -d/--duplicate one record is emitted per passing transcript;
// otherwise the transcripts' values are collapsed into one record. The -i/-e
// filter, when present, gates each emitted record after its INFO tags are set,
// exactly where upstream's filter_and_output evaluates it.
func (p *splitVepPlugin) runAnnotate(opts PluginOptions, hdr *vcf.Header, variants []*vcf.Variant, out io.Writer) error {
	outHdr := p.annotateHeader(hdr)

	filter, err := p.buildFilter(outHdr)
	if err != nil {
		return err
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
		writeErr := p.processVariant(v, filter, func(rec *vcf.Variant, _ [][]string) error {
			return w.Write(rec)
		})
		if writeErr != nil {
			cleanup()
			return writeErr
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

// processVariant runs the per-record logic shared by the annotate (-c) and text
// (-f) modes, mirroring upstream's process_record + filter_and_output: it splits
// the CSQ, selects and severity-filters transcripts, sets the requested columns
// as INFO tags (collapsed across transcripts, or one record per transcript with
// -d), applies the -i/-e filter, honours drop_sites, and invokes emit for each
// record that survives. emit receives the annotated *vcf.Variant and the slice
// of transcripts that produced it (a single transcript with -d, all passing
// transcripts otherwise) so the text renderer can format the same view the
// filter saw.
func (p *splitVepPlugin) processVariant(v *vcf.Variant, filter *pluginFilter, emit func(*vcf.Variant, [][]string) error) error {
	transcripts := p.splitCsq(v)
	if len(transcripts) == 0 {
		// No CSQ: pass through unannotated unless drop_sites.
		if p.dropSites != 0 {
			return nil
		}
		p.resetAnnots()
		return p.filterAndOutput(v, filter, true, true, nil, emit)
	}
	// -g gene restriction reorders/truncates the transcript list before selection.
	ntr := len(transcripts)
	if p.genes != nil {
		transcripts, ntr = p.restrictCsqsToGenes(transcripts)
	}
	selected, err := p.selectTranscripts(transcripts, ntr)
	if err != nil {
		return err
	}
	p.resetAnnots()
	severityPass := false
	allMissing := true
	var current [][]string
	for _, ti := range selected {
		tr := transcripts[ti]
		if p.csqIdx >= 0 && p.csqIdx < len(tr) {
			if !p.csqSeverityPass(trField(tr, p.csqIdx)) {
				continue
			}
		} else if p.minSev != svSelectAny {
			continue
		}
		severityPass = true
		current = append(current, tr)
		for ai, a := range p.annots {
			val := "."
			if a.idx == -1 {
				val = strings.Join(tr, "|")
				allMissing = false
			} else if s := trField(tr, a.idx); s != "" {
				// PRN :worst rewrites the Consequence subfield to its single worst
				// term (split on '&', exact severity lookup), matching csq_rewrite_worst.
				if a.idx == p.csqIdx && p.prnCsq == svPrnWorst {
					s = p.csqRewriteWorst(s)
				}
				val = s
				allMissing = false
			}
			p.annotVals[ai] = append(p.annotVals[ai], val)
		}
		if p.duplicate {
			if err := p.filterAndOutput(v, filter, severityPass, allMissing, current, emit); err != nil {
				return err
			}
			p.resetAnnots()
			current = current[:0]
			allMissing = true
			severityPass = false
		}
	}
	if !severityPass && p.dropSites != 0 {
		return nil
	}
	if !p.duplicate {
		return p.filterAndOutput(v, filter, severityPass, allMissing, current, emit)
	}
	return nil
}

// resetAnnots clears the per-transcript value accumulator, mirroring
// annot_reset.
func (p *splitVepPlugin) resetAnnots() {
	if p.annotVals == nil {
		p.annotVals = make([][]string, len(p.annots))
	}
	for i := range p.annotVals {
		p.annotVals[i] = p.annotVals[i][:0]
	}
}

// filterAndOutput sets the accumulated annot values as INFO tags on v, applies
// the -i/-e filter, and emits the record when it passes. It mirrors upstream's
// filter_and_output for the VCF (annotate) output mode: tags with at least one
// accumulated value are written; the filter then gates the record. drop_sites
// suppresses sites that did not update any tag (all_missing). v is annotated in
// place; with -d the caller resets between transcripts so each emitted record
// carries only that transcript's values.
func (p *splitVepPlugin) filterAndOutput(v *vcf.Variant, filter *pluginFilter, severityPass, allMissing bool, current [][]string, emit func(*vcf.Variant, [][]string) error) error {
	updated := false
	for ai, a := range p.annots {
		if len(p.annotVals[ai]) == 0 {
			continue
		}
		joined := strings.Join(p.annotVals[ai], ",")
		setInfo(v, a.tag, p.renderTyped(a.typ, joined))
		updated = true
	}
	if !filter.testSite(v) {
		return nil
	}
	if len(p.annots) > 0 {
		if p.dropSites != 0 && (!updated || allMissing) {
			return nil
		}
	} else if !severityPass {
		return nil
	}
	return emit(v, current)
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
