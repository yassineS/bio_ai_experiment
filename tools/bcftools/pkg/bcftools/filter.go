package bcftools

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Filter is a compiled boolean expression that can be evaluated against a
// vcf.Variant. The supported grammar (recursive descent):
//
//	expr     := or_expr
//	or_expr  := and_expr ("||" and_expr)*
//	and_expr := unary ("&&" unary)*
//	unary    := "!" unary | primary
//	primary  := "(" expr ")" | comparison | value
//	comparison := value op value
//	value    := tag | function | NUMBER | STRING
//	tag      := ["INFO/"|"FMT/"|"FORMAT/"] IDENT ["[" INT "]"]
//	function := "TYPE" | "N_ALT" | "QUAL" | "FILTER" | "N_PASS" | ...
//	op       := "==" | "=" | "!=" | "<" | "<=" | ">" | ">="
//
// Tag resolution follows htslib filter.c: a bare INFO-defined name resolves to
// INFO; the bare name `DP` (and other FORMAT-defined tags) resolves to FORMAT;
// `INFO/X` and `FMT/X` force the source explicitly. Multi-value INFO fields use
// "any element satisfies" comparison semantics, and `TAG[i]` selects a single
// element. `-e EXPR` is the exact complement of `-i EXPR`.
type Filter struct {
	root node
}

// CompileFilter parses an expression string into a reusable Filter without
// header type information. Bare tags are resolved with the observable
// fallback rule (INFO when present on the record, else treated as FORMAT).
// Prefer CompileFilterWithHeader so bare-tag resolution matches htslib
// exactly.
func CompileFilter(expr string) (*Filter, error) {
	return CompileFilterWithHeader(expr, nil)
}

// CompileFilterWithHeader parses an expression string into a reusable Filter,
// resolving bare (unprefixed) tags against the VCF header exactly as htslib's
// filter.c does: a bare name defined only in INFO resolves to INFO, a name
// defined only in FORMAT resolves to FORMAT, and a name defined in *both* is
// ambiguous and rejected with an error (mirroring filter.c:3434-3437). When
// hdr is nil the observable fallback rule is used instead.
func CompileFilterWithHeader(expr string, hdr *vcf.Header) (*Filter, error) {
	p := &parser{src: expr, hdr: newHeaderTags(hdr)}
	root, err := p.parseExpr()
	if err != nil {
		return nil, err
	}
	p.skipSpace()
	if p.pos != len(p.src) {
		return nil, fmt.Errorf("bcftools: trailing tokens in expression at %d", p.pos)
	}
	if p.tagErr != nil {
		return nil, p.tagErr
	}
	return &Filter{root: root}, nil
}

// headerTags records which tag names are defined in the INFO and FORMAT
// sections of a VCF header, so bare tags can be resolved the way htslib does.
type headerTags struct {
	info map[string]bool
	fmt  map[string]bool
}

// newHeaderTags scans hdr's meta lines for ##INFO/##FORMAT definitions. A nil
// header (or one with no structured lines) yields a nil-map headerTags, which
// callers treat as "no header information available".
func newHeaderTags(hdr *vcf.Header) *headerTags {
	if hdr == nil {
		return nil
	}
	ht := &headerTags{info: map[string]bool{}, fmt: map[string]bool{}}
	for _, m := range hdr.MetaInfo {
		kind, id := structuredID(m)
		switch kind {
		case "INFO":
			ht.info[id] = true
		case "FORMAT":
			ht.fmt[id] = true
		}
	}
	return ht
}

// resolve returns the tagSource for a bare tag name using header definitions.
// It mirrors filter.c: INFO-only -> tagInfo, FORMAT-only -> tagFormat,
// defined in both -> ambiguity error, defined in neither -> tagInfo (so an
// unknown bare tag is still looked up against INFO, matching htslib's default
// of is_fmt=0). A nil headerTags falls back to tagAuto (observable rule).
func (ht *headerTags) resolve(name string) (tagSource, error) {
	if ht == nil {
		return tagAuto, nil
	}
	inInfo := ht.info[name]
	inFmt := ht.fmt[name]
	switch {
	case inInfo && inFmt:
		return tagInfo, fmt.Errorf("Error: ambiguous filtering expression, both INFO/%s and FORMAT/%s are defined in the VCF header.", name, name)
	case inFmt:
		return tagFormat, nil
	default:
		// INFO-only or unknown: htslib defaults is_fmt to 0 (INFO).
		return tagInfo, nil
	}
}

// Eval returns true if v matches the compiled expression.
func (f *Filter) Eval(v *vcf.Variant) bool {
	if f == nil || f.root == nil {
		return true
	}
	res := f.root.eval(v)
	return truthy(res)
}

// node is one AST node of a compiled expression.
type node interface {
	eval(v *vcf.Variant) any
}

type binOpNode struct {
	op       string
	lhs, rhs node
}

func (n *binOpNode) eval(v *vcf.Variant) any {
	switch n.op {
	case "&&":
		return truthy(n.lhs.eval(v)) && truthy(n.rhs.eval(v))
	case "||":
		return truthy(n.lhs.eval(v)) || truthy(n.rhs.eval(v))
	case "==", "!=", "<", "<=", ">", ">=":
		// Per-sample comparison (FORMAT tags, GT) uses any-sample-by-default
		// semantics: the site passes if any used sample satisfies the operator.
		if mask, ok := n.sampleMask(v); ok {
			for _, p := range mask {
				if p {
					return true
				}
			}
			return false
		}
		return compare(n.op, n.lhs.eval(v), n.rhs.eval(v))
	}
	return false
}

// sampleMask evaluates a comparison per sample when either operand is a
// per-sample producer (FORMAT tag or GT). It returns the per-sample pass mask
// and ok=true; ok=false means neither operand is per-sample and the caller
// should fall back to scalar/site-level comparison. The mask drives both the
// any-sample default in eval and the aggregation functions (N_PASS/COUNT/...).
func (n *binOpNode) sampleMask(v *vcf.Variant) ([]bool, bool) {
	// GT genotype-class comparison: GT="het" / GT="RR" / GT="0/1" etc.
	if gt, ok := n.lhs.(*gtNode); ok {
		if mask, ok := gtCompareMask(n.op, gt, n.rhs, v); ok {
			return mask, true
		}
	}
	if gt, ok := n.rhs.(*gtNode); ok {
		if mask, ok := gtCompareMask(n.op, gt, n.lhs, v); ok {
			return mask, true
		}
	}

	lhsS, lhsOK := n.lhs.(sampleProducer)
	rhsS, rhsOK := n.rhs.(sampleProducer)
	if !lhsOK && !rhsOK {
		return nil, false
	}
	mask := make([]bool, len(v.Samples))
	switch {
	case lhsOK && rhsOK:
		lv := lhsS.evalSamples(v)
		rv := rhsS.evalSamples(v)
		for i := range mask {
			mask[i] = compare(n.op, lv[i], rv[i])
		}
	case lhsOK:
		lv := lhsS.evalSamples(v)
		rval := n.rhs.eval(v)
		for i := range mask {
			mask[i] = compare(n.op, lv[i], rval)
		}
	default:
		lval := n.lhs.eval(v)
		rv := rhsS.evalSamples(v)
		for i := range mask {
			mask[i] = compare(n.op, lval, rv[i])
		}
	}
	return mask, true
}

// gtCompareMask evaluates a GT comparison against a string literal per sample.
// The literal selects the classification mode (filter.c: hom/het/hap ->
// genotype3, mis/alt/ref -> genotype4, RR/RA/AA/aA/... -> genotype2) or, when
// it is an explicit genotype like "0/1", a direct string match against the
// formatted GT. Returns ok=false when other isn't a usable string literal.
func gtCompareMask(op string, gt *gtNode, other node, v *vcf.Variant) ([]bool, bool) {
	if op != "==" && op != "!=" {
		return nil, false
	}
	lit, ok := stringLiteral(other)
	if !ok {
		return nil, false
	}
	mode, target, isClass := gtClassTarget(lit)
	samples := gt.evalSamples(v)
	mask := make([]bool, len(samples))
	for i, sv := range samples {
		raw, _ := sv.(string)
		var eq bool
		if isClass {
			eq = gtClassify(raw, mode) == target
		} else {
			eq = gtMatchExplicit(raw, lit)
		}
		if op == "!=" {
			// filter.c only sets pass for samples with a genotype present.
			if raw == "" {
				mask[i] = false
				continue
			}
			eq = !eq
		}
		mask[i] = eq
	}
	return mask, true
}

// stringLiteral returns the string value of a node that is a string literal or
// a bare identifier used as a literal (e.g. unquoted het).
func stringLiteral(n node) (string, bool) {
	switch x := n.(type) {
	case *literalNode:
		if s, ok := x.value.(string); ok {
			return s, true
		}
	case *identNode:
		return x.name, true
	}
	return "", false
}

// gtClassTarget maps a GT class keyword to (mode, lowercased-target, true).
// Modes mirror filter.c: 3 = {hom,het,hap}, 4 = {mis,alt,ref}, 2 = {rr,ra,aa,
// aA,r,a}. Unknown keywords are treated as explicit genotype strings
// (isClass=false).
func gtClassTarget(lit string) (mode int, target string, isClass bool) {
	// "aA"/"Aa" are case-sensitive het-alt; check before the case-insensitive
	// switch so they are not swallowed by the lowercase "aa" (hom-alt) case.
	if lit == "aA" || lit == "Aa" {
		return 2, "aA", true
	}
	switch strings.ToLower(lit) {
	case "hom":
		return 3, "hom", true
	case "het":
		return 3, "het", true
	case "hap":
		return 3, "hap", true
	case "mis":
		return 4, "mis", true
	case "alt":
		return 4, "alt", true
	case "ref":
		return 4, "ref", true
	case "rr":
		return 2, "rr", true
	case "ra", "ar":
		return 2, "ra", true
	case "aa":
		return 2, "aa", true
	case "r":
		return 2, "r", true
	case "a":
		return 2, "a", true
	}
	return 0, "", false
}

// gtClassify classifies a raw GT string (e.g. "0/1", "1|1", "./.") into the
// class for the requested mode, matching _filters_set_genotype's logic. An
// empty/fully-missing genotype classifies as "mis" (mode 4) or "." otherwise.
func gtClassify(raw string, mode int) string {
	alleles, missing, ploidy := parseGT(raw)
	if ploidy == 0 || missing {
		if mode == 4 {
			return "mis"
		}
		return "."
	}
	hasRef := false
	isHet := false
	for i, a := range alleles {
		if a == 0 {
			hasRef = true
		}
		if i > 0 && a != alleles[i-1] {
			isHet = true
		}
	}
	switch mode {
	case 4:
		switch {
		case !hasRef:
			return "alt"
		case !isHet:
			return "ref"
		default:
			return "alt"
		}
	case 3:
		switch {
		case ploidy == 1:
			return "hap"
		case !isHet:
			return "hom"
		default:
			return "het"
		}
	default: // mode 2
		if ploidy == 1 {
			if hasRef {
				return "r"
			}
			return "a"
		}
		if !isHet {
			if hasRef {
				return "rr"
			}
			return "aa"
		}
		if hasRef {
			return "ra"
		}
		return "aA"
	}
}

// parseGT parses a VCF GT string into per-allele integer indices. It reports
// missing=true if any allele is missing ('.') and the ploidy (number of
// alleles). Phasing separators '|' and '/' are both accepted.
func parseGT(raw string) (alleles []int, missing bool, ploidy int) {
	if raw == "" || raw == "." {
		return nil, true, 0
	}
	parts := strings.FieldsFunc(raw, func(r rune) bool { return r == '/' || r == '|' })
	if len(parts) == 0 {
		return nil, true, 0
	}
	alleles = make([]int, 0, len(parts))
	for _, p := range parts {
		if p == "." {
			return nil, true, len(parts)
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			return nil, true, len(parts)
		}
		alleles = append(alleles, n)
	}
	return alleles, false, len(alleles)
}

// gtMatchExplicit reports whether a sample's GT matches an explicit genotype
// literal like "0/1" or "1|1". Matching is unphased-aware: the separator in the
// literal is honored, but allele order is compared as written (matching
// filter.c's string comparison of bcf_format_gt output). To be tolerant of the
// phasing separator, "0/1" also matches "0|1" only when the literal uses '/'
// — filter.c compares the formatted string verbatim, so we mirror that and only
// normalize when the literal omits phasing by using '/'.
func gtMatchExplicit(raw, lit string) bool {
	if raw == "" {
		return false
	}
	if raw == lit {
		return true
	}
	// filter.c compares bcf_format_gt output verbatim; an unphased literal
	// ("0/1") will not match a phased genotype ("0|1"). We compare verbatim.
	return false
}

type notNode struct{ inner node }

func (n *notNode) eval(v *vcf.Variant) any { return !truthy(n.inner.eval(v)) }

type literalNode struct{ value any }

func (n *literalNode) eval(*vcf.Variant) any { return n.value }

// tagSource selects where a tag node reads its value from.
type tagSource int

const (
	tagAuto   tagSource = iota // bare name: resolve INFO vs FORMAT by definition
	tagInfo                    // INFO/X
	tagFormat                  // FMT/X or FORMAT/X
)

// tagNode reads a named VCF field. index >= 0 selects a single element of a
// multi-value field; index < 0 means "the whole vector" (any-element compare).
type tagNode struct {
	name   string
	source tagSource
	index  int
}

func (n *tagNode) eval(v *vcf.Variant) any {
	switch n.source {
	case tagInfo:
		return n.indexed(infoValue(v, n.name))
	case tagFormat:
		// Per-sample FORMAT tags are not evaluated at the site level here; a
		// site-level test against a FORMAT tag matches nothing (mirrors the
		// observable behaviour of bare `DP` resolving to FORMAT/DP).
		return nil
	default: // tagAuto
		// htslib resolves a bare tag by header definition: prefer INFO unless
		// the tag is only defined in FORMAT (e.g. DP -> FORMAT/DP). We don't
		// have header type info threaded here, so we apply the observable
		// rule: a bare tag resolves to INFO when an INFO value is present,
		// otherwise it is treated as a FORMAT tag (-> no site-level match).
		if _, ok := v.Info[n.name]; ok {
			return n.indexed(infoValue(v, n.name))
		}
		return nil
	}
}

// indexed applies the optional [i] selector to a resolved value.
func (n *tagNode) indexed(val any) any {
	if n.index < 0 {
		return val
	}
	vec, ok := val.(vecValue)
	if !ok {
		if n.index == 0 {
			return val
		}
		return nil
	}
	if n.index >= len(vec) {
		return nil
	}
	return vec[n.index]
}

// infoValue returns the parsed value of an INFO tag: nil when absent, true for
// a present flag (empty value), a scalar for a single value, or a vecValue for
// a comma-separated multi-value field.
func infoValue(v *vcf.Variant, key string) any {
	raw, ok := v.Info[key]
	if !ok {
		return nil
	}
	if raw == "" {
		// Flag-style INFO: present means true.
		return true
	}
	if !strings.Contains(raw, ",") {
		return coerce(raw)
	}
	parts := strings.Split(raw, ",")
	vec := make(vecValue, len(parts))
	for i, p := range parts {
		vec[i] = coerce(p)
	}
	return vec
}

// coerce parses a raw INFO token into a float64 when numeric, else the string.
func coerce(s string) any {
	if s == "" || s == "." {
		return s
	}
	if f, err := strconv.ParseFloat(s, 64); err == nil {
		return f
	}
	return s
}

// vecValue is a multi-element field value; comparisons use any-element
// semantics (true when any element satisfies the operator).
type vecValue []any

type qualNode struct{}

func (n *qualNode) eval(v *vcf.Variant) any {
	if v.Qual < 0 {
		// -1 is our missing-QUAL sentinel.
		return nil
	}
	return v.Qual
}

type nAltNode struct{}

func (n *nAltNode) eval(v *vcf.Variant) any { return float64(len(v.Alt)) }

// typeNode evaluates to the set of variant type names of the record, using
// htslib's singular spellings (snp, indel, mnp, bnd, ref, other). The set is a
// vecValue so that `TYPE="snp"` is true when any allele is a SNP.
type typeNode struct{}

func (n *typeNode) eval(v *vcf.Variant) any {
	types := variantTypeSet(v)
	vec := make(vecValue, len(types))
	for i, t := range types {
		vec[i] = t
	}
	return vec
}

// variantTypeSet returns the per-allele variant types using the singular
// htslib filter spellings. Mirrors bcf_get_variant_types: each ALT is
// classified independently; REF-only records are "ref".
func variantTypeSet(v *vcf.Variant) []string {
	if len(v.Alt) == 0 {
		return []string{"ref"}
	}
	refLen := len(v.Ref)
	out := make([]string, 0, len(v.Alt))
	for _, a := range v.Alt {
		switch {
		case a == "" || a == "." || a == "*":
			out = append(out, "other")
		case strings.ContainsAny(a, "[]"):
			out = append(out, "bnd")
		case strings.HasPrefix(a, "<") && strings.HasSuffix(a, ">"):
			out = append(out, "other")
		case len(a) == 1 && refLen == 1:
			out = append(out, "snp")
		case len(a) == refLen:
			out = append(out, "mnp")
		default:
			out = append(out, "indel")
		}
	}
	return out
}

type filterNode struct{}

func (n *filterNode) eval(v *vcf.Variant) any {
	if len(v.Filter) == 0 {
		return "."
	}
	return strings.Join(v.Filter, ";")
}

// chromNode evaluates to the record's CHROM column (a string).
type chromNode struct{}

func (n *chromNode) eval(v *vcf.Variant) any { return v.Chrom }

// posNode evaluates to the record's POS column (a 1-based numeric position).
type posNode struct{}

func (n *posNode) eval(v *vcf.Variant) any { return float64(v.Pos) }

// idNode evaluates to the record's ID column. Like filter.c, a missing ID (".")
// compares as the literal ".".
type idNode struct{}

func (n *idNode) eval(v *vcf.Variant) any {
	if v.ID == "" {
		return "."
	}
	return v.ID
}

// refNode evaluates to the record's REF allele string.
type refNode struct{}

func (n *refNode) eval(v *vcf.Variant) any { return v.Ref }

// altNode evaluates to the record's ALT allele(s). With index < 0 it is the set
// of all ALT alleles (any-element compare); index >= 0 selects a single ALT
// (0-based), mirroring filter.c's ALT[i]. A record with no ALT yields ".".
type altNode struct{ index int }

func (n *altNode) eval(v *vcf.Variant) any {
	if n.index >= 0 {
		if n.index < len(v.Alt) {
			return v.Alt[n.index]
		}
		return "."
	}
	if len(v.Alt) == 0 {
		return "."
	}
	if len(v.Alt) == 1 {
		return v.Alt[0]
	}
	vec := make(vecValue, len(v.Alt))
	for i, a := range v.Alt {
		vec[i] = a
	}
	return vec
}

// sampleProducer is implemented by nodes that yield one value per sample
// (FORMAT tags and GT). evalSamples returns a slice with len(v.Samples)
// entries; an entry is nil when that sample has no value (e.g. missing GT).
// This lets comparison and aggregation nodes evaluate per-sample with the
// any-sample-by-default semantics of filter.c.
type sampleProducer interface {
	evalSamples(v *vcf.Variant) []any
}

// gtNode reads the per-sample GT genotype. As a bare site-level value it has no
// scalar meaning; it is only useful inside a comparison (e.g. GT="het") or an
// aggregation, where its per-sample classification/strings are used.
type gtNode struct{}

// eval makes gtNode usable as a truthy operand (rare); it reports whether any
// sample has a non-missing genotype.
func (n *gtNode) eval(v *vcf.Variant) any {
	for _, s := range v.Samples {
		if gt, ok := s.Data["GT"]; ok && gt != "" && gt != "." && gt != "./." && gt != ".|." {
			return true
		}
	}
	return false
}

// evalSamples returns the raw per-sample GT strings (e.g. "0/1"), or nil for a
// sample lacking GT. Genotype-class matching ("hom"/"het"/"RR"/...) is applied
// in the comparison logic, which recognizes the class keywords on the literal.
func (n *gtNode) evalSamples(v *vcf.Variant) []any {
	out := make([]any, len(v.Samples))
	for i, s := range v.Samples {
		gt, ok := s.Data["GT"]
		if !ok {
			out[i] = nil
			continue
		}
		out[i] = gt
	}
	return out
}

// formatNode reads a per-sample FORMAT tag (FMT/X, FORMAT/X). index >= 0
// selects a single element of a multi-value field per sample.
type formatNode struct {
	name  string
	index int
}

// eval makes a bare FORMAT tag usable as a truthy site-level operand: it is
// true when any sample has a truthy value. Comparisons use evalSamples.
func (n *formatNode) eval(v *vcf.Variant) any {
	for _, val := range n.evalSamples(v) {
		if truthy(val) {
			return true
		}
	}
	return false
}

// evalSamples returns the per-sample value of the FORMAT tag.
func (n *formatNode) evalSamples(v *vcf.Variant) []any {
	out := make([]any, len(v.Samples))
	for i, s := range v.Samples {
		raw, ok := s.Data[n.name]
		if !ok {
			out[i] = nil
			continue
		}
		out[i] = n.indexed(formatValue(raw))
	}
	return out
}

// indexed applies the optional [i] element selector to a per-sample value.
func (n *formatNode) indexed(val any) any {
	if n.index < 0 {
		return val
	}
	vec, ok := val.(vecValue)
	if !ok {
		if n.index == 0 {
			return val
		}
		return nil
	}
	if n.index >= len(vec) {
		return nil
	}
	return vec[n.index]
}

// formatValue parses a per-sample FORMAT token: nil for a missing field ("."),
// a scalar for a single value, or a vecValue for a comma-separated list.
func formatValue(raw string) any {
	if raw == "" || raw == "." {
		return nil
	}
	if !strings.Contains(raw, ",") {
		return coerce(raw)
	}
	parts := strings.Split(raw, ",")
	vec := make(vecValue, len(parts))
	for i, p := range parts {
		vec[i] = coerce(p)
	}
	return vec
}

// absNode implements ABS(x): the absolute value of a numeric field. Vectors are
// mapped element-wise.
type absNode struct{ inner node }

func (n *absNode) eval(v *vcf.Variant) any {
	val := n.inner.eval(v)
	switch x := val.(type) {
	case nil:
		return nil
	case vecValue:
		out := make(vecValue, len(x))
		for i, e := range x {
			out[i] = absOf(e)
		}
		return out
	default:
		return absOf(x)
	}
}

func absOf(v any) any {
	if f, ok := asFloat(v); ok {
		if f < 0 {
			return -f
		}
		return f
	}
	return v
}

// aggNode implements the sample-aggregating functions over an inner expression:
// N_PASS/F_PASS (count/fraction of samples whose per-sample comparison passed),
// COUNT (number of passing samples, or number of present values for a bare
// FORMAT tag), and SUM/MAX/MIN/AVG over per-sample numeric values. It mirrors
// filter.c's func_npass / func_count / func_sum family.
type aggNode struct {
	kind  string
	inner node
}

func (n *aggNode) eval(v *vcf.Variant) any {
	switch n.kind {
	case "N_PASS", "F_PASS", "COUNT":
		npass, nsmpl := n.passCounts(v)
		switch n.kind {
		case "N_PASS":
			return float64(npass)
		case "F_PASS":
			if nsmpl == 0 {
				return float64(0)
			}
			return float64(npass) / float64(nsmpl)
		default: // COUNT
			return float64(npass)
		}
	case "SUM", "MAX", "MIN", "AVG":
		return n.reduce(v)
	}
	return nil
}

// passCounts returns the number of samples whose inner per-sample comparison
// passed and the number of used samples. For a bare FORMAT tag (COUNT(FMT/X))
// it counts samples with a present value.
func (n *aggNode) passCounts(v *vcf.Variant) (npass, nsmpl int) {
	if cmp, ok := n.inner.(*binOpNode); ok {
		if mask, ok := cmp.sampleMask(v); ok {
			for _, p := range mask {
				nsmpl++
				if p {
					npass++
				}
			}
			return npass, nsmpl
		}
	}
	if sp, ok := n.inner.(sampleProducer); ok {
		// COUNT(FMT/TAG): number of present (non-missing) per-sample values.
		for _, val := range sp.evalSamples(v) {
			nsmpl++
			if val != nil {
				countNonMissing(val, &npass)
			}
		}
		return npass, nsmpl
	}
	// Fall back: a truthy site-level inner counts as 1.
	if truthy(n.inner.eval(v)) {
		return 1, 1
	}
	return 0, 1
}

// countNonMissing adds the number of non-missing scalar values in val to *cnt.
func countNonMissing(val any, cnt *int) {
	if vec, ok := val.(vecValue); ok {
		for _, e := range vec {
			if e != nil {
				*cnt++
			}
		}
		return
	}
	*cnt++
}

// reduce computes SUM/MAX/MIN/AVG over the inner expression's per-sample (or
// site-level vector) numeric values.
func (n *aggNode) reduce(v *vcf.Variant) any {
	var vals []float64
	if sp, ok := n.inner.(sampleProducer); ok {
		for _, sv := range sp.evalSamples(v) {
			collectFloats(sv, &vals)
		}
	} else {
		collectFloats(n.inner.eval(v), &vals)
	}
	if len(vals) == 0 {
		return nil
	}
	switch n.kind {
	case "SUM":
		s := 0.0
		for _, x := range vals {
			s += x
		}
		return s
	case "MAX":
		m := vals[0]
		for _, x := range vals[1:] {
			if x > m {
				m = x
			}
		}
		return m
	case "MIN":
		m := vals[0]
		for _, x := range vals[1:] {
			if x < m {
				m = x
			}
		}
		return m
	default: // AVG
		s := 0.0
		for _, x := range vals {
			s += x
		}
		return s / float64(len(vals))
	}
}

// collectFloats appends the numeric values found in val (scalar or vector) to
// *out, skipping nil and non-numeric entries.
func collectFloats(val any, out *[]float64) {
	switch x := val.(type) {
	case nil:
		return
	case vecValue:
		for _, e := range x {
			if f, ok := asFloat(e); ok {
				*out = append(*out, f)
			}
		}
	default:
		if f, ok := asFloat(x); ok {
			*out = append(*out, f)
		}
	}
}

// strlenNode implements STRLEN(field): the byte length of a string field. A
// missing value (".") has length 0, matching filter.c's func_strlen. When the
// inner field is multi-valued (e.g. STRLEN(ALT)) it yields a vector of lengths
// with any-element compare semantics.
type strlenNode struct{ inner node }

func (n *strlenNode) eval(v *vcf.Variant) any {
	val := n.inner.eval(v)
	switch x := val.(type) {
	case nil:
		return nil
	case vecValue:
		out := make(vecValue, 0, len(x))
		for _, e := range x {
			out = append(out, float64(strlenOf(e)))
		}
		return out
	default:
		return float64(strlenOf(x))
	}
}

// strlenOf returns the STRLEN of a single value: the byte length of its string
// form, with "." treated as length 0.
func strlenOf(v any) int {
	s := asString(v)
	if s == "." {
		return 0
	}
	return len(s)
}

// truthy turns a Go value into a boolean using the rules bcftools follows:
// strings are true when non-empty and not ".", numbers are true when
// nonzero, booleans pass through, vectors are true when any element is, and
// nil is false.
func truthy(v any) bool {
	switch x := v.(type) {
	case nil:
		return false
	case bool:
		return x
	case string:
		return x != "" && x != "."
	case float64:
		return x != 0
	case int:
		return x != 0
	case vecValue:
		for _, e := range x {
			if truthy(e) {
				return true
			}
		}
		return false
	}
	return false
}

// compare evaluates a comparison operator. Vectors compare with any-element
// semantics; a nil operand (missing field) never matches.
func compare(op string, a, b any) bool {
	if a == nil || b == nil {
		return false
	}
	av, aIsVec := a.(vecValue)
	bv, bIsVec := b.(vecValue)
	switch {
	case aIsVec && bIsVec:
		for _, ae := range av {
			for _, be := range bv {
				if compareScalar(op, ae, be) {
					return true
				}
			}
		}
		return false
	case aIsVec:
		for _, ae := range av {
			if compareScalar(op, ae, b) {
				return true
			}
		}
		return false
	case bIsVec:
		for _, be := range bv {
			if compareScalar(op, a, be) {
				return true
			}
		}
		return false
	default:
		return compareScalar(op, a, b)
	}
}

// compareScalar compares two non-vector values, preferring numeric comparison
// and falling back to lexical comparison when either side is non-numeric.
func compareScalar(op string, a, b any) bool {
	if a == nil || b == nil {
		return false
	}
	if af, aok := asFloat(a); aok {
		if bf, bok := asFloat(b); bok {
			return cmpFloat(op, af, bf)
		}
	}
	return cmpString(op, asString(a), asString(b))
}

func cmpFloat(op string, a, b float64) bool {
	switch op {
	case "==":
		return a == b
	case "!=":
		return a != b
	case "<":
		return a < b
	case "<=":
		return a <= b
	case ">":
		return a > b
	case ">=":
		return a >= b
	}
	return false
}

func cmpString(op string, a, b string) bool {
	switch op {
	case "==":
		return a == b
	case "!=":
		return a != b
	case "<":
		return a < b
	case "<=":
		return a <= b
	case ">":
		return a > b
	case ">=":
		return a >= b
	}
	return false
}

func asFloat(v any) (float64, bool) {
	switch x := v.(type) {
	case float64:
		return x, true
	case int:
		return float64(x), true
	case string:
		f, err := strconv.ParseFloat(x, 64)
		if err != nil {
			return 0, false
		}
		return f, true
	}
	return 0, false
}

func asString(v any) string {
	switch x := v.(type) {
	case string:
		return x
	case bool:
		if x {
			return "1"
		}
		return "0"
	case float64:
		return strconv.FormatFloat(x, 'g', -1, 64)
	case int:
		return strconv.Itoa(x)
	}
	return ""
}

// parser is the recursive-descent parser for filter expressions.
type parser struct {
	src    string
	pos    int
	hdr    *headerTags // header tag definitions for bare-tag resolution (may be nil)
	tagErr error       // first bare-tag resolution error (e.g. ambiguity); surfaced by Compile
}

func (p *parser) skipSpace() {
	for p.pos < len(p.src) && (p.src[p.pos] == ' ' || p.src[p.pos] == '\t') {
		p.pos++
	}
}

func (p *parser) peek() byte {
	if p.pos >= len(p.src) {
		return 0
	}
	return p.src[p.pos]
}

func (p *parser) match(s string) bool {
	if strings.HasPrefix(p.src[p.pos:], s) {
		p.pos += len(s)
		return true
	}
	return false
}

func (p *parser) parseExpr() (node, error) {
	return p.parseOr()
}

func (p *parser) parseOr() (node, error) {
	left, err := p.parseAnd()
	if err != nil {
		return nil, err
	}
	for {
		p.skipSpace()
		if !p.match("||") {
			return left, nil
		}
		right, err := p.parseAnd()
		if err != nil {
			return nil, err
		}
		left = &binOpNode{op: "||", lhs: left, rhs: right}
	}
}

func (p *parser) parseAnd() (node, error) {
	left, err := p.parseUnary()
	if err != nil {
		return nil, err
	}
	for {
		p.skipSpace()
		if !p.match("&&") {
			return left, nil
		}
		right, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		left = &binOpNode{op: "&&", lhs: left, rhs: right}
	}
}

func (p *parser) parseUnary() (node, error) {
	p.skipSpace()
	if p.match("!") {
		// "!=" is a comparison operator, not unary negation. We just consumed
		// the '!' though, so peek the next char.
		if p.peek() == '=' {
			// rewind one char so the comparison parser sees the full token.
			p.pos--
		} else {
			inner, err := p.parseUnary()
			if err != nil {
				return nil, err
			}
			return &notNode{inner: inner}, nil
		}
	}
	return p.parsePrimary()
}

func (p *parser) parsePrimary() (node, error) {
	p.skipSpace()
	if p.pos >= len(p.src) {
		return nil, fmt.Errorf("bcftools: unexpected end of expression")
	}
	if p.match("(") {
		inner, err := p.parseExpr()
		if err != nil {
			return nil, err
		}
		p.skipSpace()
		if !p.match(")") {
			return nil, fmt.Errorf("bcftools: missing ')' at %d", p.pos)
		}
		return inner, nil
	}
	left, err := p.parseValue()
	if err != nil {
		return nil, err
	}
	// If a comparison operator follows, build a comparison node.
	op := p.matchCompareOp()
	if op == "" {
		return left, nil
	}
	right, err := p.parseValue()
	if err != nil {
		return nil, err
	}
	return &binOpNode{op: op, lhs: left, rhs: right}, nil
}

func (p *parser) matchCompareOp() string {
	p.skipSpace()
	// Order matters: match the longer operators first so that "==" is not
	// accidentally consumed as two "=" tokens, and "<=" / ">=" are seen
	// before bare "<" / ">".
	switch {
	case p.match("=="):
		return "=="
	case p.match("!="):
		return "!="
	case p.match("<="):
		return "<="
	case p.match(">="):
		return ">="
	case p.match("="):
		// bcftools accepts the single-equal spelling as a synonym for ==.
		return "=="
	case p.match("<"):
		return "<"
	case p.match(">"):
		return ">"
	}
	return ""
}

func (p *parser) parseValue() (node, error) {
	p.skipSpace()
	if p.pos >= len(p.src) {
		return nil, fmt.Errorf("bcftools: unexpected end of expression")
	}
	c := p.src[p.pos]
	switch {
	case c == '"' || c == '\'':
		return p.parseString(c)
	case c == '-' || c == '+' || (c >= '0' && c <= '9') || c == '.':
		return p.parseNumber()
	case isIdentStart(c):
		return p.parseIdent()
	}
	return nil, fmt.Errorf("bcftools: unexpected character %q at %d", c, p.pos)
}

func (p *parser) parseString(quote byte) (node, error) {
	p.pos++ // skip opening quote
	start := p.pos
	for p.pos < len(p.src) && p.src[p.pos] != quote {
		p.pos++
	}
	if p.pos >= len(p.src) {
		return nil, fmt.Errorf("bcftools: unterminated string literal")
	}
	s := p.src[start:p.pos]
	p.pos++ // skip closing quote
	return &literalNode{value: s}, nil
}

func (p *parser) parseNumber() (node, error) {
	start := p.pos
	if p.src[p.pos] == '-' || p.src[p.pos] == '+' {
		p.pos++
	}
	for p.pos < len(p.src) && (isDigit(p.src[p.pos]) || p.src[p.pos] == '.' || p.src[p.pos] == 'e' || p.src[p.pos] == 'E' || p.src[p.pos] == '-' || p.src[p.pos] == '+') {
		// Allow only one minus sign after 'e'/'E'. Simpler: stop at the
		// first '-' that doesn't follow an exponent letter.
		if (p.src[p.pos] == '-' || p.src[p.pos] == '+') && p.pos > start {
			prev := p.src[p.pos-1]
			if prev != 'e' && prev != 'E' {
				break
			}
		}
		p.pos++
	}
	tok := p.src[start:p.pos]
	f, err := strconv.ParseFloat(tok, 64)
	if err != nil {
		return nil, fmt.Errorf("bcftools: bad number %q", tok)
	}
	return &literalNode{value: f}, nil
}

func (p *parser) parseIdent() (node, error) {
	start := p.pos
	for p.pos < len(p.src) && (isIdentPart(p.src[p.pos]) || p.src[p.pos] == '/') {
		p.pos++
	}
	tok := p.src[start:p.pos]

	// Function call: NAME( ... ). Recognized before the [i] selector so that
	// e.g. N_PASS(...) and STRLEN(...) are dispatched to the function parser.
	if p.peek() == '(' {
		if fn, ok, err := p.parseFunctionCall(tok); ok || err != nil {
			return fn, err
		}
	}

	// Optional [i] index selector immediately following the identifier.
	index := -1
	if p.peek() == '[' {
		save := p.pos
		p.pos++ // consume '['
		numStart := p.pos
		for p.pos < len(p.src) && isDigit(p.src[p.pos]) {
			p.pos++
		}
		if p.pos > numStart && p.peek() == ']' {
			idx, err := strconv.Atoi(p.src[numStart:p.pos])
			if err == nil {
				index = idx
				p.pos++ // consume ']'
			} else {
				p.pos = save
			}
		} else {
			p.pos = save
		}
	}

	// Explicit source prefixes.
	switch {
	case strings.HasPrefix(tok, "INFO/"):
		return &tagNode{name: tok[len("INFO/"):], source: tagInfo, index: index}, nil
	case strings.HasPrefix(tok, "FORMAT/"):
		name := tok[len("FORMAT/"):]
		if name == "GT" {
			return &gtNode{}, nil
		}
		return &formatNode{name: name, index: index}, nil
	case strings.HasPrefix(tok, "FMT/"):
		name := tok[len("FMT/"):]
		if name == "GT" {
			return &gtNode{}, nil
		}
		return &formatNode{name: name, index: index}, nil
	}

	// Special keywords / functions.
	switch tok {
	case "FILTER":
		return &filterNode{}, nil
	case "QUAL":
		return &qualNode{}, nil
	case "TYPE":
		return &typeNode{}, nil
	case "N_ALT":
		return &nAltNode{}, nil
	case "CHROM":
		return &chromNode{}, nil
	case "POS":
		return &posNode{}, nil
	case "ID":
		return &idNode{}, nil
	case "REF":
		return &refNode{}, nil
	case "ALT":
		return &altNode{index: index}, nil
	case "GT":
		// Bare GT resolves to FORMAT/GT (the genotype). See filter.c:3454.
		return &gtNode{}, nil
	case "true":
		return &literalNode{value: true}, nil
	case "false":
		return &literalNode{value: false}, nil
	}

	// A bare identifier is a tag reference. When header type information is
	// available we resolve INFO vs FORMAT now (and record the ambiguity error
	// for the compiler to surface), exactly as htslib does at parse time.
	// Without a header we defer to identNode's observable fallback rule so the
	// node can still behave as a string literal on the RHS of e.g.
	// `FILTER="PASS"` written without quotes.
	src, err := p.hdr.resolve(tok)
	if err != nil && p.tagErr == nil {
		p.tagErr = err
	}
	if p.hdr != nil {
		return &tagNode{name: tok, source: src, index: index}, nil
	}
	if index >= 0 {
		return &tagNode{name: tok, source: tagAuto, index: index}, nil
	}
	return &identNode{name: tok}, nil
}

// parseFunctionCall parses a recognized function-call form NAME(expr). It is
// called with the consumed identifier name and the parser positioned at '('. It
// returns ok=false (without consuming the '(') for names that are not known
// functions, so the caller can fall back to treating NAME as a tag. Function
// names are matched case-insensitively to mirror filter.c's strncasecmp checks.
func (p *parser) parseFunctionCall(name string) (node, bool, error) {
	upper := strings.ToUpper(name)
	var kind string
	switch upper {
	case "N_PASS", "F_PASS", "COUNT", "SUM", "MAX", "MIN", "AVG", "MEAN", "STRLEN", "ABS":
		kind = upper
	default:
		return nil, false, nil
	}
	if p.peek() != '(' {
		return nil, false, nil
	}
	p.pos++ // consume '('
	inner, err := p.parseExpr()
	if err != nil {
		return nil, true, err
	}
	p.skipSpace()
	if !p.match(")") {
		return nil, true, fmt.Errorf("bcftools: missing ')' in %s() at %d", name, p.pos)
	}
	switch kind {
	case "STRLEN":
		return &strlenNode{inner: inner}, true, nil
	case "ABS":
		return &absNode{inner: inner}, true, nil
	case "N_PASS":
		return &aggNode{kind: "N_PASS", inner: inner}, true, nil
	case "F_PASS":
		return &aggNode{kind: "F_PASS", inner: inner}, true, nil
	case "COUNT":
		return &aggNode{kind: "COUNT", inner: inner}, true, nil
	case "SUM":
		return &aggNode{kind: "SUM", inner: inner}, true, nil
	case "MAX":
		return &aggNode{kind: "MAX", inner: inner}, true, nil
	case "MIN":
		return &aggNode{kind: "MIN", inner: inner}, true, nil
	case "AVG", "MEAN":
		return &aggNode{kind: "AVG", inner: inner}, true, nil
	}
	return nil, false, nil
}

// identNode is a bare identifier whose meaning depends on context: when used as
// the left side of a comparison it resolves as a tag; otherwise it is a string
// literal (so `FILTER="PASS"` / `TYPE="snp"` compare against the literal). The
// disambiguation happens in eval: identNode resolves to a tag value when one
// exists for the record, otherwise to the literal name. To keep `-e` the exact
// complement of `-i`, the same resolution is used in every position.
type identNode struct{ name string }

func (n *identNode) eval(v *vcf.Variant) any {
	if val := (&tagNode{name: n.name, source: tagAuto, index: -1}).eval(v); val != nil {
		return val
	}
	// Fall back to the literal string (handles RHS of `=="PASS"`-style
	// comparisons written without quotes, and unknown tags).
	return n.name
}

func isIdentStart(c byte) bool {
	return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_'
}

func isIdentPart(c byte) bool {
	return isIdentStart(c) || isDigit(c)
}

func isDigit(c byte) bool { return c >= '0' && c <= '9' }
