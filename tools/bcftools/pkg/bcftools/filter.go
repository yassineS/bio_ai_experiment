package bcftools

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// Filter is a compiled boolean expression that can be evaluated against a
// vcf.Variant. The supported grammar (recursive descent):
//
//	expr     := or_expr
//	or_expr  := and_expr ("||" and_expr)*
//	and_expr := unary ("&&" unary)*
//	unary    := "!" unary | primary
//	primary  := "(" expr ")" | comparison | identifier
//	comparison := value op value
//	value    := "INFO/" IDENT | "FILTER" | NUMBER | STRING | IDENT
//	op       := "==" | "!=" | "<" | "<=" | ">" | ">="
//
// Identifiers may also be the bare keyword `FILTER` (compares against the
// joined FILTER column). Quoted strings use double quotes.
type Filter struct {
	root node
}

// CompileFilter parses an expression string into a reusable Filter.
func CompileFilter(expr string) (*Filter, error) {
	p := &parser{src: expr}
	root, err := p.parseExpr()
	if err != nil {
		return nil, err
	}
	p.skipSpace()
	if p.pos != len(p.src) {
		return nil, fmt.Errorf("bcftools: trailing tokens in expression at %d", p.pos)
	}
	return &Filter{root: root}, nil
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

type infoNode struct{ key string }

func (n *infoNode) eval(v *vcf.Variant) any {
	val, ok := v.Info[n.key]
	if !ok {
		return nil
	}
	if val == "" {
		// Flag-style INFO: present means true.
		return true
	}
	// Try numeric coercion lazily; comparisons handle string/number mix.
	return val
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
// nonzero, booleans pass through, and nil is false.
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
	}
	return false
}

// compare evaluates a comparison operator on two values. It tries numeric
// comparison first (coercing strings to float64); falling back to lexical
// comparison if either side is non-numeric.
func compare(op string, a, b any) bool {
	if af, aok := asFloat(a); aok {
		if bf, bok := asFloat(b); bok {
			return cmpFloat(op, af, bf)
		}
	}
	as, bs := asString(a), asString(b)
	return cmpString(op, as, bs)
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
		// Multi-value INFO fields use commas; take the first numeric.
		s := x
		if i := strings.IndexByte(s, ','); i >= 0 {
			s = s[:i]
		}
		f, err := strconv.ParseFloat(s, 64)
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
	src string
	pos int
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
	switch {
	case strings.HasPrefix(tok, "INFO/"):
		return &infoNode{key: tok[len("INFO/"):]}, nil
	case tok == "FILTER":
		return &filterNode{}, nil
	case tok == "true":
		return &literalNode{value: true}, nil
	case tok == "false":
		return &literalNode{value: false}, nil
	}
	// Bare identifier: treat as a string literal so `FILTER=PASS` works.
	return &literalNode{value: tok}, nil
}

func isIdentStart(c byte) bool {
	return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_'
}

func isIdentPart(c byte) bool {
	return isIdentStart(c) || isDigit(c)
}

func isDigit(c byte) bool { return c >= '0' && c <= '9' }
