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
		return compare(n.op, n.lhs.eval(v), n.rhs.eval(v))
	}
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
		return &tagNode{name: tok[len("FORMAT/"):], source: tagFormat, index: index}, nil
	case strings.HasPrefix(tok, "FMT/"):
		return &tagNode{name: tok[len("FMT/"):], source: tagFormat, index: index}, nil
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
