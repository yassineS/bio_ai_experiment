// Shared -i/--include / -e/--exclude pre-filter plumbing for the native plugins
// that use upstream's filter_init/filter_test as a record (and, where the
// expression names FORMAT fields, per-sample) gate before the plugin processes
// a record. It re-implements the include/exclude bookkeeping that every such
// plugin's .c repeats verbatim (see plugins/smpl-stats.c, indel-stats.c,
// trio-stats.c, guess-ploidy.c and contrast.c): compile one expression, evaluate
// it per record, and — for an FLT_EXCLUDE expression — invert the verdict, with
// the subtle per-sample inversion semantics upstream applies when the expression
// names a FORMAT field.
package bcftools

import "github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"

// pluginFilter is a compiled -i/-e expression plus its include/exclude logic,
// shared by the stats and contrast native plugins. The zero value (and a nil
// *pluginFilter) is the "no filter" case: every record and every sample passes.
type pluginFilter struct {
	filter  *Filter
	exclude bool // true for -e/--exclude, false for -i/--include
}

// newPluginFilter compiles expr into a pluginFilter without header knowledge.
// exclude selects the FLT_EXCLUDE logic (-e/--exclude); otherwise FLT_INCLUDE
// (-i/--include) is used. An empty expression yields a nil filter (no-op). It
// returns an error if the expression does not compile. Prefer
// newPluginFilterWithHeader when the plugin's header is available so bare
// FORMAT-only tags resolve to FORMAT exactly as upstream does.
func newPluginFilter(expr string, exclude bool) (*pluginFilter, error) {
	return newPluginFilterWithHeader(expr, exclude, nil)
}

// newPluginFilterWithHeader is newPluginFilter with header-aware bare-name
// resolution: a bare tag declared only as FORMAT becomes a per-sample FORMAT
// term and a tag declared as both INFO and FORMAT is rejected as ambiguous,
// matching upstream filter.c. A nil header behaves like newPluginFilter.
func newPluginFilterWithHeader(expr string, exclude bool, hdr *vcf.Header) (*pluginFilter, error) {
	if expr == "" {
		return nil, nil
	}
	f, err := CompileFilterWithHeader(expr, hdr)
	if err != nil {
		return nil, err
	}
	return &pluginFilter{filter: f, exclude: exclude}, nil
}

// testSite returns the record-level verdict for a site-only filter, mirroring
//
//	int pass = filter_test(filter, rec, NULL);
//	if ( filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
//
// used by contrast and split, where smpl_pass is never requested. A nil
// pluginFilter always passes.
func (pf *pluginFilter) testSite(v *vcf.Variant) bool {
	if pf == nil || pf.filter == nil {
		return true
	}
	pass := pf.filter.Eval(v)
	if pf.exclude {
		return !pass
	}
	return pass
}

// testSamples returns the record-level verdict and a per-sample include mask,
// faithfully reproducing the FLT_INCLUDE / FLT_EXCLUDE bookkeeping the
// per-sample stats plugins repeat verbatim around
// filter_test(filter, rec, &smpl_pass):
//
//   - For a site-only expression (no FORMAT term, mask == nil) the verdict is
//     the site pass (inverted for exclude) and mask stays nil, so the caller
//     processes every sample when the site passes.
//   - For a per-sample expression the returned mask has one bool per sample
//     (true == keep). For include the mask is the raw filter mask; for exclude
//     each sample bit is inverted (a sample that matched the expression is
//     dropped). When the raw site verdict is false (nothing matched), exclude
//     keeps the whole record with an all-true mask, matching upstream's
//     "for (i...) smpl_pass[i] = 1" branch.
//
// When passSite is false the caller must skip the record entirely; mask is then
// nil. A nil pluginFilter passes every record with a nil mask.
func (pf *pluginFilter) testSamples(v *vcf.Variant) (passSite bool, mask []bool) {
	if pf == nil || pf.filter == nil {
		return true, nil
	}
	siteMatch, raw := pf.filter.EvalSamples(v)
	if !pf.exclude {
		// FLT_INCLUDE: drop the record only when nothing matched. The mask (if
		// any) selects the matching samples directly.
		if raw == nil {
			return siteMatch, nil
		}
		if !siteMatch {
			return false, nil
		}
		return true, raw
	}
	// FLT_EXCLUDE.
	if raw == nil {
		// Site-level exclude: invert the site verdict, no per-sample mask.
		return !siteMatch, nil
	}
	if !siteMatch {
		// Nothing matched the expression: keep the record with all samples.
		return true, allTrue(len(raw))
	}
	// Something matched: invert each sample bit and keep the record only if at
	// least one sample is now retained.
	out := make([]bool, len(raw))
	any := false
	for i, m := range raw {
		out[i] = !m
		if out[i] {
			any = true
		}
	}
	if !any {
		return false, nil
	}
	return true, out
}

// rawSamples returns the unmodified filter_test result: the raw site-match flag
// and the raw per-sample match mask (nil for a site-only expression), BEFORE any
// FLT_EXCLUDE inversion. It is used by the per-trio plugins (trio-stats and
// indel-stats with -p), which fold the per-sample mask into a per-trio verdict
// and so need the raw bits plus this filter's exclude flag, rather than the
// per-sample inversion testSamples performs. A nil pluginFilter reports a
// matching site with no mask.
func (pf *pluginFilter) rawSamples(v *vcf.Variant) (siteMatch bool, mask []bool, exclude bool) {
	if pf == nil || pf.filter == nil {
		return true, nil, false
	}
	m, raw := pf.filter.EvalSamples(v)
	return m, raw, pf.exclude
}

// allTrue returns a length-n slice of true values, the all-samples-pass mask.
func allTrue(n int) []bool {
	m := make([]bool, n)
	for i := range m {
		m[i] = true
	}
	return m
}
