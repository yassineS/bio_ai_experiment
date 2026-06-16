package bcftools

import (
	"regexp"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// This file extends the filter engine with FORMAT/sample-level support,
// mirroring the relevant semantics of upstream bcftools filter.c. The
// scalar node.eval(v) tree in filter.go remains the home for site/INFO-level
// expressions; here we add a parallel token-based evaluator that produces a
// per-sample mask so that an expression containing FORMAT or GT terms can be
// evaluated for every sample, with the site passing when ANY sample matches
// (upstream's TOK_AND_VEC / TOK_OR_VEC "can be true in different samples"
// rule, filter.c vector_logic_and / vector_logic_or).

// tokResult is the result of evaluating one (sub)expression against a record,
// modelled on upstream filter.c's token_t pass_site / pass_samples pair.
//
//   - hasSamples is false for purely site-level terms (INFO tags, QUAL, FILTER,
//     literals). Only passSite is meaningful for those.
//   - hasSamples is true for terms that touched a FORMAT/GT tag; perSample then
//     holds one boolean per sample and passSite is the OR across samples.
type tokResult struct {
	passSite   bool
	hasSamples bool
	perSample  []bool
}

// fmtNode is implemented by AST nodes that participate in per-sample
// evaluation. Every node in the tree implements it; the site-level nodes
// defined in filter.go gain a default implementation via the scalarNode
// adapter (see evalTok below).
type fmtNode interface {
	node
	// evalTok evaluates the node against v with nsmpl samples, returning a
	// tokResult.
	evalTok(v *vcf.Variant, nsmpl int) tokResult
}

// EvalSamples evaluates the filter against v and returns the per-sample match
// mask alongside the site-level pass flag. For a site-level (INFO/QUAL/FILTER)
// expression hasSamples is false and the mask is nil; callers that only need
// the site verdict can keep using Eval. The mask has one entry per sample in
// v.Samples and is true for samples that satisfy the expression, matching the
// semantics of upstream filter_test(filter, line, &smpl_pass).
func (f *Filter) EvalSamples(v *vcf.Variant) (passSite bool, mask []bool) {
	if f == nil || f.root == nil {
		return true, nil
	}
	res := evalTok(f.root, v, len(v.Samples))
	if res.hasSamples {
		return res.passSite, res.perSample
	}
	return res.passSite, nil
}

// Eval returns true if v matches the compiled expression at the site level.
// An expression containing FORMAT/GT terms passes the site when ANY sample
// satisfies it, matching upstream bcftools -i semantics.
//
// This overrides the simpler scalar Eval in filter.go: it routes through the
// token evaluator so FORMAT/GT terms are honoured while site/INFO-only
// expressions evaluate identically to before.
func (f *Filter) Eval(v *vcf.Variant) bool {
	if f == nil || f.root == nil {
		return true
	}
	return evalTok(f.root, v, len(v.Samples)).passSite
}

// evalTok dispatches to a node's per-sample evaluator when it implements
// fmtNode, otherwise it treats the node as a site-level scalar (the literal /
// INFO / FILTER / field nodes from filter.go) and lifts its truthiness into a
// site-only tokResult.
func evalTok(n node, v *vcf.Variant, nsmpl int) tokResult {
	if fn, ok := n.(fmtNode); ok {
		return fn.evalTok(v, nsmpl)
	}
	return tokResult{passSite: truthy(n.eval(v))}
}

// fullMask returns a per-sample slice with every entry set to b.
func fullMask(nsmpl int, b bool) []bool {
	m := make([]bool, nsmpl)
	if b {
		for i := range m {
			m[i] = true
		}
	}
	return m
}

// orMask returns true if any sample bit is set.
func orMask(m []bool) bool {
	for _, b := range m {
		if b {
			return true
		}
	}
	return false
}

// --- logic combiners (binOpNode for && / ||) -------------------------------

// evalTok for binOpNode implements both the boolean logic operators (&& / ||)
// and the comparison operators, with per-sample propagation following
// filter.c's vector_logic_and / vector_logic_or.
func (n *binOpNode) evalTok(v *vcf.Variant, nsmpl int) tokResult {
	switch n.op {
	case "&&":
		return logicAnd(evalTok(n.lhs, v, nsmpl), evalTok(n.rhs, v, nsmpl), nsmpl)
	case "||":
		return logicOr(evalTok(n.lhs, v, nsmpl), evalTok(n.rhs, v, nsmpl), nsmpl)
	case "==", "!=", "<", "<=", ">", ">=", "~", "!~":
		return compareTok(n.op, n.lhs, n.rhs, v, nsmpl)
	}
	return tokResult{passSite: false}
}

// logicAnd mirrors vector_logic_and with the && (TOK_AND_VEC) operator: the
// site fails unless both sides pass, and when both sides are per-sample the
// result mask is the OR of the two masks ("can be true in different samples").
// A mix of one site-level and one per-sample side carries the per-sample mask
// through unchanged (gated on the site side passing).
func logicAnd(a, b tokResult, nsmpl int) tokResult {
	if !a.passSite || !b.passSite {
		return tokResult{passSite: false, hasSamples: a.hasSamples || b.hasSamples, perSample: fullMask(nsmpl, false)}
	}
	if !a.hasSamples && !b.hasSamples {
		return tokResult{passSite: true}
	}
	if a.hasSamples != b.hasSamples {
		s := a
		if b.hasSamples {
			s = b
		}
		return tokResult{passSite: orMask(s.perSample), hasSamples: true, perSample: append([]bool(nil), s.perSample...)}
	}
	m := make([]bool, nsmpl)
	for i := range m {
		m[i] = a.perSample[i] || b.perSample[i]
	}
	return tokResult{passSite: true, hasSamples: true, perSample: m}
}

// logicOr mirrors vector_logic_or with the || (TOK_OR_VEC) operator. When both
// sides are per-sample the mask is the OR of the two; when one side is a
// passing site-level term, upstream sets every used sample (selects all
// samples if one is true), except for the QUAL>30 || FMT/GQ>30 guard where a
// failing site-level side leaves only the per-sample matches.
func logicOr(a, b tokResult, nsmpl int) tokResult {
	if !a.passSite && !b.passSite {
		return tokResult{passSite: false, hasSamples: a.hasSamples || b.hasSamples, perSample: fullMask(nsmpl, false)}
	}
	if !a.hasSamples && !b.hasSamples {
		return tokResult{passSite: true}
	}
	if a.hasSamples != b.hasSamples {
		smp, site := a, b
		if b.hasSamples {
			smp, site = b, a
		}
		// Guard (filter.c): if the site-level side is false, carry only the
		// per-sample matches; otherwise the passing site selects all samples.
		if !site.passSite {
			return tokResult{passSite: orMask(smp.perSample), hasSamples: true, perSample: append([]bool(nil), smp.perSample...)}
		}
		return tokResult{passSite: true, hasSamples: true, perSample: fullMask(nsmpl, true)}
	}
	m := make([]bool, nsmpl)
	for i := range m {
		m[i] = a.perSample[i] || b.perSample[i]
	}
	return tokResult{passSite: orMask(m), hasSamples: true, perSample: m}
}

// evalTok for notNode negates the site verdict. Negation collapses to the site
// level (upstream does not propagate a per-sample mask through `!`), so a
// negated per-sample expression contributes only its site truthiness.
func (n *notNode) evalTok(v *vcf.Variant, nsmpl int) tokResult {
	inner := evalTok(n.inner, v, nsmpl)
	return tokResult{passSite: !inner.passSite}
}

// --- comparison -------------------------------------------------------------

// compareTok evaluates a comparison whose operands may be site-level or
// per-sample (FORMAT/GT) operands. If neither operand is per-sample it reduces
// to the scalar compare from filter.go. Otherwise the comparison is applied
// per sample and the site passes if any sample matches.
func compareTok(op string, lhs, rhs node, v *vcf.Variant, nsmpl int) tokResult {
	lp, lok := lhs.(perSampleNode)
	rp, rok := rhs.(perSampleNode)
	if !lok && !rok {
		return tokResult{passSite: cmpTokScalar(op, lhs.eval(v), rhs.eval(v))}
	}

	// Determine the per-sample operand and the (possibly scalar) other side.
	mask := make([]bool, nsmpl)
	for i := 0; i < nsmpl; i++ {
		var lv, rv any
		var lpres, rpres bool
		if lok {
			lv, lpres = lp.sampleValue(v, i)
		} else {
			lv, lpres = lhs.eval(v), true
		}
		if rok {
			rv, rpres = rp.sampleValue(v, i)
		} else {
			rv, rpres = rhs.eval(v), true
		}
		if !lpres || !rpres {
			// Missing value: comparison is false for that sample (mirrors the
			// default missing_logic where == / != against present values fail
			// unless explicitly testing equality to ".").
			mask[i] = cmpMissing(op, lv, lpres, rv, rpres)
			continue
		}
		mask[i] = cmpTokScalar(op, lv, rv)
	}
	return tokResult{passSite: orMask(mask), hasSamples: true, perSample: mask}
}

// cmpMissing handles a per-sample comparison where one operand is a missing
// FORMAT value. Equality to a literal "." matches; inequality to "." matches a
// present value; all other comparisons against a missing value are false.
func cmpMissing(op string, lv any, lpres bool, rv any, rpres bool) bool {
	// Identify the literal side, if any.
	var lit string
	var litSet bool
	if lpres {
		lit, litSet = asString(lv), true
	} else if rpres {
		lit, litSet = asString(rv), true
	}
	switch op {
	case "==":
		return litSet && lit == "."
	case "!=":
		return litSet && lit != "."
	}
	return false
}

// cmpTokScalar performs a single comparison, supporting the regex operators
// (~ / !~) in addition to the arithmetic/relational operators handled by
// compare(). Multi-valued fields (comma-separated, e.g. an INFO Number=A tag
// such as AC=2,1) are matched element-wise with ANY semantics: the comparison
// passes if any element of one side compares true against any element of the
// other, mirroring upstream filter.c's vector comparison (a Number=A/G/. tag
// passes the site when any of its values satisfies the operator). This applies
// uniformly to every operator, including != and the regex operators.
func cmpTokScalar(op string, a, b any) bool {
	as, aMulti := splitMultiValue(a)
	bs, bMulti := splitMultiValue(b)
	if !aMulti && !bMulti {
		return cmpScalarOne(op, a, b)
	}
	for _, ae := range as {
		for _, be := range bs {
			if cmpScalarOne(op, ae, be) {
				return true
			}
		}
	}
	return false
}

// splitMultiValue splits v into its comma-separated elements when v is a
// multi-valued string (contains a comma). It returns the element slice and
// whether v was multi-valued. A scalar value yields a single-element slice and
// false. Only string values are split; numeric/boolean values are scalar.
func splitMultiValue(v any) ([]any, bool) {
	s, ok := v.(string)
	if !ok || !strings.Contains(s, ",") {
		return []any{v}, false
	}
	parts := strings.Split(s, ",")
	out := make([]any, len(parts))
	for i, p := range parts {
		out[i] = p
	}
	return out, true
}

// cmpScalarOne performs a single scalar comparison on already-singular
// operands, dispatching the regex operators (~ / !~) and otherwise delegating
// to compare().
func cmpScalarOne(op string, a, b any) bool {
	switch op {
	case "~", "!~":
		matched := regexMatch(asString(a), asString(b))
		if op == "!~" {
			return !matched
		}
		return matched
	}
	return compare(op, a, b)
}

// regexMatch reports whether s matches the POSIX extended regular expression
// pattern, mirroring upstream's regexec-based ~ operator. A pattern that fails
// to compile never matches.
func regexMatch(s, pattern string) bool {
	re, err := regexp.CompilePOSIX(pattern)
	if err != nil {
		re, err = regexp.Compile(pattern)
		if err != nil {
			return false
		}
	}
	return re.MatchString(s)
}

// perSampleNode is implemented by value nodes that resolve to a per-sample
// value (FORMAT tags and GT). sampleValue returns the value for sample i and
// whether it is present (false for a missing FORMAT field).
type perSampleNode interface {
	node
	sampleValue(v *vcf.Variant, i int) (any, bool)
}

// --- FORMAT tag node --------------------------------------------------------

// formatNode resolves a FORMAT/<tag> (or bare FORMAT-only tag) per sample.
type formatNode struct {
	key string
}

// eval on a formatNode returns the first sample's value so that the node can
// still participate as a scalar where a per-sample mask is not requested.
func (n *formatNode) eval(v *vcf.Variant) any {
	if len(v.Samples) == 0 {
		return nil
	}
	val, ok := n.sampleValue(v, 0)
	if !ok {
		return nil
	}
	return val
}

// sampleValue returns the FORMAT value for sample i. Multi-value fields keep
// their comma-separated text; numeric coercion happens lazily in compare().
func (n *formatNode) sampleValue(v *vcf.Variant, i int) (any, bool) {
	if i >= len(v.Samples) {
		return nil, false
	}
	val, ok := v.Samples[i].Data[n.key]
	if !ok || val == "" || val == "." {
		return nil, false
	}
	return val, true
}

// evalTok evaluates a bare FORMAT tag used on its own (truthiness per sample).
func (n *formatNode) evalTok(v *vcf.Variant, nsmpl int) tokResult {
	mask := make([]bool, nsmpl)
	for i := 0; i < nsmpl; i++ {
		val, ok := n.sampleValue(v, i)
		mask[i] = ok && truthy(val)
	}
	return tokResult{passSite: orMask(mask), hasSamples: true, perSample: mask}
}

// --- GT node ----------------------------------------------------------------

// gtNode resolves FORMAT/GT for genotype-class and exact-string tests. It is
// only ever the left-hand side of a comparison whose right side is a literal
// string (a class keyword such as "het" or an exact genotype such as "0/1").
type gtNode struct{}

// eval returns the first sample's raw GT string for scalar use.
func (n *gtNode) eval(v *vcf.Variant) any {
	if len(v.Samples) == 0 {
		return nil
	}
	gt, ok := v.Samples[0].Data["GT"]
	if !ok {
		return nil
	}
	return gt
}

// sampleValue returns the raw GT string for sample i.
func (n *gtNode) sampleValue(v *vcf.Variant, i int) (any, bool) {
	if i >= len(v.Samples) {
		return nil, false
	}
	gt, ok := v.Samples[i].Data["GT"]
	if !ok {
		return nil, false
	}
	return gt, true
}

// evalTok of a bare GT term is truthiness per sample (a non-missing genotype).
func (n *gtNode) evalTok(v *vcf.Variant, nsmpl int) tokResult {
	mask := make([]bool, nsmpl)
	for i := 0; i < nsmpl; i++ {
		gt, ok := n.sampleValue(v, i)
		if !ok {
			continue
		}
		g, ok := parseGT(asString(gt))
		mask[i] = ok && !g.isFullyMissing()
	}
	return tokResult{passSite: orMask(mask), hasSamples: true, perSample: mask}
}

// gtCompareNode handles a comparison whose left side is GT and right side a
// string literal. It dispatches between genotype-class keywords and exact
// genotype-string matching, honouring =/==, != and ~/!~.
type gtCompareNode struct {
	op  string
	rhs string // literal comparand
}

func (n *gtCompareNode) eval(v *vcf.Variant) any {
	return n.evalTok(v, len(v.Samples)).passSite
}

// evalTok evaluates the GT comparison per sample.
func (n *gtCompareNode) evalTok(v *vcf.Variant, nsmpl int) tokResult {
	mask := make([]bool, nsmpl)
	cls, isClass := gtClassKeyword(n.rhs)
	for i := 0; i < nsmpl; i++ {
		raw, ok := gtRaw(v, i)
		if !ok {
			continue
		}
		var matched bool
		if isClass {
			matched = gtMatchesClass(raw, cls)
		} else {
			matched = gtMatchesString(raw, n.rhs, n.op)
		}
		switch n.op {
		case "!=", "!~":
			mask[i] = !matched
		default:
			mask[i] = matched
		}
	}
	return tokResult{passSite: orMask(mask), hasSamples: true, perSample: mask}
}

// gtRaw returns the raw GT string of sample i.
func gtRaw(v *vcf.Variant, i int) (string, bool) {
	if i >= len(v.Samples) {
		return "", false
	}
	gt, ok := v.Samples[i].Data["GT"]
	return gt, ok
}

// gtClassKeyword reports whether s is a recognised GT class keyword, returning
// the canonical lower-case form. The keywords are case-insensitive, matching
// upstream's strcasecmp dispatch in filter.c (~line 4150): the type-3 set
// (hom/het/hap), the type-4 set (mis/alt/ref), and the type-2 set
// (rr/ra/ar/aa/aA/r/a). The aA / Aa token is the only case-sensitive one
// (filter.c uses strcmp for it), so it is matched before the case-insensitive
// fold. Tokens that are not class keywords (e.g. "phased", or an exact
// genotype such as "0/1") fall through to exact-string matching, exactly as
// upstream does — upstream has no phased/unphased pseudo-class.
func gtClassKeyword(s string) (string, bool) {
	if s == "aA" || s == "Aa" {
		return "aA", true
	}
	switch strings.ToLower(s) {
	case "het", "hom", "hap", "mis", "alt", "ref",
		"rr", "ra", "ar", "aa", "r", "a":
		return strings.ToLower(s), true
	}
	return "", false
}

// gtMatchesClass reports whether the raw genotype string matches the class
// keyword, reproducing filter.c's _filters_set_genotype classification
// (filter.c:1222-1274):
//
//   - type 4 (mis/alt/ref): missing genotype => mis; no ref allele => alt;
//     ref-only (homozygous ref or haploid ref) => ref; heterozygous with a ref
//     allele => alt.
//   - type 3 (hap/hom/het): haploid (ploidy 1, present) => hap; all alleles
//     equal => hom; otherwise het.
//   - type 2 (rr/ra/aa/aA/r/a): finer split by ref presence and ploidy.
func gtMatchesClass(raw, cls string) bool {
	g, ok := parseGT(raw)
	if !ok {
		return false
	}
	ploidy := g.ploidy()
	missing := g.nMissing() > 0 || ploidy == 0
	hasRef := false
	isHet := false
	for j, a := range g.alleles {
		if a == 0 {
			hasRef = true
		}
		if j > 0 && a != g.alleles[j-1] {
			isHet = true
		}
	}

	switch cls {
	// type 4
	case "mis":
		return missing
	case "alt":
		return !missing && (!hasRef || isHet)
	case "ref":
		return !missing && hasRef && !isHet
	// type 3
	case "hap":
		return !missing && ploidy == 1
	case "hom":
		return !missing && ploidy != 1 && !isHet
	case "het":
		return !missing && isHet
	// type 2 (haploid)
	case "r":
		return !missing && ploidy == 1 && hasRef
	case "a":
		return !missing && ploidy == 1 && !hasRef
	// type 2 (diploid+)
	case "rr":
		return !missing && ploidy != 1 && !isHet && hasRef
	case "aa":
		return !missing && ploidy != 1 && !isHet && !hasRef
	case "ra", "ar":
		return !missing && ploidy != 1 && isHet && hasRef
	case "aA":
		return !missing && ploidy != 1 && isHet && !hasRef
	}
	return false
}

// gtMatchesString matches an exact genotype string such as "0/1" or "0|1".
// For == / != upstream compares the formatted genotype text; phase and allele
// order are significant exactly as written. The ~ / !~ operators apply the
// literal as a regular expression over the genotype text.
func gtMatchesString(raw, lit, op string) bool {
	switch op {
	case "~", "!~":
		return regexMatch(raw, lit)
	}
	return raw == lit
}
