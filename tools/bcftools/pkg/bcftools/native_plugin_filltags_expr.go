// Expression evaluator for the fill-tags custom-expression surface
// (TAG[:Number]=[int|integer|float](EXPR)). This is a self-contained port of
// the slice of the upstream bcftools filter engine (filter.c) that the
// fill-tags plugin exercises: numeric INFO/FORMAT tag references, the
// aggregation functions (SUM/AVG|MEAN/MAX/MIN/MEDIAN/STDEV and their
// per-sample SMPL_* / sXXX variants), arithmetic (+ - * /), unary minus, ABS,
// PHRED, and the genotype-counting reductions F_MISSING, F_PASS, N_PASS,
// N_MISSING.
//
// Evaluation produces, for a record, either a single site-level slice of
// doubles (Number=.) or one slice per sample (per-sample / FORMAT result),
// exactly as filter_get_doubles does upstream. Missing/vector-end values are
// represented by NaN and excluded from aggregations, matching
// bcf_double_is_missing_or_vector_end.
package bcftools

import (
	"fmt"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// exprMissing is the sentinel double used for a missing or vector-end value,
// mirroring htslib's bcf_double_is_missing_or_vector_end.
var exprMissing = math.NaN()

func isExprMissing(x float64) bool { return math.IsNaN(x) }

// exprResult is the outcome of evaluating a (sub)expression against a record.
// When perSample is true, values holds nsmpl slices (one per sample, indexed
// by sample); the i-th sample's vector is values[i]. When perSample is false
// it is a single site-level vector in site.
type exprResult struct {
	perSample bool
	site      []float64   // site-level vector (perSample==false)
	values    [][]float64 // per-sample vectors (perSample==true)
}

// fillExpr is a compiled fill-tags expression.
type fillExpr struct {
	root    exprNode
	nsmpl   int
	exprStr string
}

// exprNode is a node in the parsed expression tree.
type exprNode interface {
	eval(ctx *exprCtx) exprResult
}

// exprCtx carries the per-record evaluation state.
type exprCtx struct {
	v     *vcf.Variant
	nsmpl int
	usmpl []bool // active sample mask (population restriction); nil means all
}

// compileFillExpr parses the EXPR portion of a TAG=EXPR specification against
// the header so INFO vs FORMAT references resolve correctly.
func compileFillExpr(exprStr string, hdr *vcf.Header) (*fillExpr, error) {
	p := &exprParser{src: exprStr, hdr: hdr}
	node, err := p.parse()
	if err != nil {
		return nil, err
	}
	return &fillExpr{root: node, nsmpl: len(hdr.Samples), exprStr: exprStr}, nil
}

// evaluate runs the expression against v with the given active-sample mask
// (nil = all samples active).
func (e *fillExpr) evaluate(v *vcf.Variant, usmpl []bool) exprResult {
	ctx := &exprCtx{v: v, nsmpl: len(v.Samples), usmpl: usmpl}
	return e.root.eval(ctx)
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

type exprParser struct {
	src string
	pos int
	hdr *vcf.Header
}

func (p *exprParser) parse() (exprNode, error) {
	n, err := p.parseAddSub()
	if err != nil {
		return nil, err
	}
	p.skipSpace()
	if p.pos != len(p.src) {
		return nil, fmt.Errorf("fill-tags: trailing characters in expression at %d: %q", p.pos, p.src[p.pos:])
	}
	return n, nil
}

func (p *exprParser) skipSpace() {
	for p.pos < len(p.src) && (p.src[p.pos] == ' ' || p.src[p.pos] == '\t') {
		p.pos++
	}
}

func (p *exprParser) parseAddSub() (exprNode, error) {
	left, err := p.parseMulDiv()
	if err != nil {
		return nil, err
	}
	for {
		p.skipSpace()
		if p.pos >= len(p.src) {
			break
		}
		op := p.src[p.pos]
		if op != '+' && op != '-' {
			break
		}
		p.pos++
		right, err := p.parseMulDiv()
		if err != nil {
			return nil, err
		}
		left = &arithNode{op: rune(op), lhs: left, rhs: right}
	}
	return left, nil
}

func (p *exprParser) parseMulDiv() (exprNode, error) {
	left, err := p.parseUnary()
	if err != nil {
		return nil, err
	}
	for {
		p.skipSpace()
		if p.pos >= len(p.src) {
			break
		}
		op := p.src[p.pos]
		if op != '*' && op != '/' {
			break
		}
		p.pos++
		right, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		left = &arithNode{op: rune(op), lhs: left, rhs: right}
	}
	return left, nil
}

func (p *exprParser) parseUnary() (exprNode, error) {
	p.skipSpace()
	if p.pos < len(p.src) && p.src[p.pos] == '-' {
		p.pos++
		inner, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return &negNode{inner: inner}, nil
	}
	return p.parsePrimary()
}

func (p *exprParser) parsePrimary() (exprNode, error) {
	p.skipSpace()
	if p.pos >= len(p.src) {
		return nil, fmt.Errorf("fill-tags: unexpected end of expression")
	}
	c := p.src[p.pos]
	if c == '(' {
		p.pos++
		inner, err := p.parseAddSub()
		if err != nil {
			return nil, err
		}
		p.skipSpace()
		if p.pos >= len(p.src) || p.src[p.pos] != ')' {
			return nil, fmt.Errorf("fill-tags: missing ')' in expression")
		}
		p.pos++
		return inner, nil
	}
	if c >= '0' && c <= '9' || c == '.' {
		return p.parseNumber()
	}
	return p.parseIdent()
}

func (p *exprParser) parseNumber() (exprNode, error) {
	start := p.pos
	for p.pos < len(p.src) {
		ch := p.src[p.pos]
		if (ch >= '0' && ch <= '9') || ch == '.' || ch == 'e' || ch == 'E' {
			p.pos++
			continue
		}
		if (ch == '+' || ch == '-') && p.pos > start {
			prev := p.src[p.pos-1]
			if prev == 'e' || prev == 'E' {
				p.pos++
				continue
			}
		}
		break
	}
	f, err := strconv.ParseFloat(p.src[start:p.pos], 64)
	if err != nil {
		return nil, fmt.Errorf("fill-tags: bad number %q", p.src[start:p.pos])
	}
	return &constNode{val: f}, nil
}

func (p *exprParser) parseIdent() (exprNode, error) {
	start := p.pos
	for p.pos < len(p.src) {
		ch := p.src[p.pos]
		if ch == '/' {
			// A '/' continues the identifier only when it forms one of the
			// INFO/, FORMAT/ or FMT/ tag-source prefixes; otherwise it is the
			// division operator and terminates the token.
			prefix := strings.ToUpper(p.src[start:p.pos])
			if prefix == "INFO" || prefix == "FORMAT" || prefix == "FMT" {
				p.pos++
				continue
			}
			break
		}
		if ch == '(' || ch == ')' || ch == '+' || ch == '-' || ch == '*' || ch == ' ' || ch == '\t' || ch == ',' {
			break
		}
		p.pos++
	}
	name := p.src[start:p.pos]
	p.skipSpace()
	if p.pos < len(p.src) && p.src[p.pos] == '(' {
		return p.parseCall(name)
	}
	// Bare F_MISSING / N_MISSING are built-in genotype-missingness reductions
	// (filter.c TOK_FUNC), equivalent to F_PASS(GT="mis") / N_PASS(GT="mis").
	switch strings.ToUpper(name) {
	case "F_MISSING", "N_MISSING":
		return newGTReduceNode(strings.ToUpper(name), "", p.hdr)
	}
	return p.makeTagNode(name)
}

func (p *exprParser) parseCall(name string) (exprNode, error) {
	p.pos++ // consume '('
	upper := strings.ToUpper(name)

	switch upper {
	case "ABS":
		inner, err := p.parseSingleArg(name)
		if err != nil {
			return nil, err
		}
		return &absNode{inner: inner}, nil
	case "PHRED":
		inner, err := p.parseSingleArg(name)
		if err != nil {
			return nil, err
		}
		return &phredNode{inner: inner}, nil
	case "F_MISSING", "F_PASS", "N_PASS", "N_MISSING":
		inner, err := p.parseRawParen()
		if err != nil {
			return nil, err
		}
		return newGTReduceNode(upper, inner, p.hdr)
	}

	agg, perSample, ok := aggKind(upper)
	if !ok {
		return nil, fmt.Errorf("fill-tags: unsupported function %q in expression", name)
	}
	inner, err := p.parseSingleArg(name)
	if err != nil {
		return nil, err
	}
	return &aggNode{kind: agg, perSample: perSample, inner: inner}, nil
}

// parseSingleArg parses a single inner expression up to the matching ')'.
func (p *exprParser) parseSingleArg(name string) (exprNode, error) {
	inner, err := p.parseAddSub()
	if err != nil {
		return nil, err
	}
	p.skipSpace()
	if p.pos >= len(p.src) || p.src[p.pos] != ')' {
		return nil, fmt.Errorf("fill-tags: missing ')' after %s(", name)
	}
	p.pos++
	return inner, nil
}

// parseRawParen returns the raw text inside the current parentheses (the '('
// was already consumed), balancing nested parentheses.
func (p *exprParser) parseRawParen() (string, error) {
	depth := 1
	start := p.pos
	for p.pos < len(p.src) {
		switch p.src[p.pos] {
		case '(':
			depth++
		case ')':
			depth--
			if depth == 0 {
				inner := p.src[start:p.pos]
				p.pos++ // consume ')'
				return inner, nil
			}
		}
		p.pos++
	}
	return "", fmt.Errorf("fill-tags: missing ')' in expression")
}

// makeTagNode builds a reference to an INFO or FORMAT tag (or a bare name,
// resolved against the header: FORMAT preferred when it exists, matching the
// upstream filter engine).
func (p *exprParser) makeTagNode(name string) (exprNode, error) {
	if name == "" {
		return nil, fmt.Errorf("fill-tags: empty operand in expression")
	}
	isFmt := false
	tag := name
	switch {
	case strings.HasPrefix(strings.ToUpper(name), "INFO/"):
		tag = name[5:]
	case strings.HasPrefix(strings.ToUpper(name), "FORMAT/"):
		tag = name[7:]
		isFmt = true
	case strings.HasPrefix(strings.ToUpper(name), "FMT/"):
		tag = name[4:]
		isFmt = true
	default:
		if headerHasFormat(p.hdr, tag) {
			isFmt = true
		}
	}
	return &tagNode{tag: tag, isFormat: isFmt}, nil
}

// ---------------------------------------------------------------------------
// Nodes
// ---------------------------------------------------------------------------

type constNode struct{ val float64 }

func (n *constNode) eval(*exprCtx) exprResult {
	return exprResult{perSample: false, site: []float64{n.val}}
}

type tagNode struct {
	tag      string
	isFormat bool
}

func (n *tagNode) eval(ctx *exprCtx) exprResult {
	if n.isFormat {
		vals := make([][]float64, ctx.nsmpl)
		for i := range ctx.v.Samples {
			vals[i] = parseExprFloatList(ctx.v.Samples[i].Data[n.tag])
		}
		return exprResult{perSample: true, values: vals}
	}
	return exprResult{perSample: false, site: parseExprFloatList(ctx.v.Info[n.tag])}
}

type negNode struct{ inner exprNode }

func (n *negNode) eval(ctx *exprCtx) exprResult {
	return mapResult(n.inner.eval(ctx), func(x float64) float64 {
		if isExprMissing(x) {
			return x
		}
		return -x
	})
}

type absNode struct{ inner exprNode }

func (n *absNode) eval(ctx *exprCtx) exprResult {
	return mapResult(n.inner.eval(ctx), func(x float64) float64 {
		if isExprMissing(x) {
			return x
		}
		return math.Abs(x)
	})
}

type phredNode struct{ inner exprNode }

func (n *phredNode) eval(ctx *exprCtx) exprResult {
	return mapResult(n.inner.eval(ctx), func(x float64) float64 {
		if isExprMissing(x) {
			return x
		}
		// upstream PHRED uses -4.34294481903251 * log(x) (i.e. -10*log10(x)).
		return -4.34294481903251 * math.Log(x)
	})
}

type arithNode struct {
	op       rune
	lhs, rhs exprNode
}

func (n *arithNode) eval(ctx *exprCtx) exprResult {
	a := n.lhs.eval(ctx)
	b := n.rhs.eval(ctx)
	return combineResults(a, b, ctx.nsmpl, func(x, y float64) float64 {
		if isExprMissing(x) || isExprMissing(y) {
			return exprMissing
		}
		switch n.op {
		case '+':
			return x + y
		case '-':
			return x - y
		case '*':
			return x * y
		case '/':
			if y == 0 {
				return exprMissing
			}
			return x / y
		}
		return exprMissing
	})
}

// aggregation kinds.
type aggOp int

const (
	aggSum aggOp = iota
	aggAvg
	aggMax
	aggMin
	aggMedian
	aggStdev
)

func aggKind(upper string) (aggOp, bool, bool) {
	switch upper {
	case "SUM":
		return aggSum, false, true
	case "AVG", "MEAN":
		return aggAvg, false, true
	case "MAX":
		return aggMax, false, true
	case "MIN":
		return aggMin, false, true
	case "MEDIAN":
		return aggMedian, false, true
	case "STDEV":
		return aggStdev, false, true
	case "SMPL_SUM", "SSUM":
		return aggSum, true, true
	case "SMPL_AVG", "SMPL_MEAN", "SAVG", "SMEAN":
		return aggAvg, true, true
	case "SMPL_MAX", "SMAX":
		return aggMax, true, true
	case "SMPL_MIN", "SMIN":
		return aggMin, true, true
	case "SMPL_MEDIAN", "SMEDIAN":
		return aggMedian, true, true
	case "SMPL_STDEV", "SSTDEV":
		return aggStdev, true, true
	}
	return 0, false, false
}

type aggNode struct {
	kind      aggOp
	perSample bool
	inner     exprNode
}

func (n *aggNode) eval(ctx *exprCtx) exprResult {
	r := n.inner.eval(ctx)
	if n.perSample && r.perSample {
		out := make([][]float64, ctx.nsmpl)
		for i := 0; i < ctx.nsmpl; i++ {
			if ctx.usmpl != nil && !ctx.usmpl[i] {
				continue
			}
			val, ok := reduce(n.kind, r.values[i])
			if ok {
				out[i] = []float64{val}
			} else {
				out[i] = []float64{exprMissing}
			}
		}
		return exprResult{perSample: true, values: out}
	}
	var flat []float64
	if r.perSample {
		for i := 0; i < ctx.nsmpl; i++ {
			if ctx.usmpl != nil && !ctx.usmpl[i] {
				continue
			}
			flat = append(flat, r.values[i]...)
		}
	} else {
		flat = r.site
	}
	val, ok := reduce(n.kind, flat)
	if !ok {
		return exprResult{perSample: false, site: nil}
	}
	return exprResult{perSample: false, site: []float64{val}}
}

// reduce applies the aggregation, excluding missing values. ok is false when
// there are no non-missing values.
func reduce(op aggOp, vals []float64) (float64, bool) {
	clean := make([]float64, 0, len(vals))
	for _, x := range vals {
		if !isExprMissing(x) {
			clean = append(clean, x)
		}
	}
	if len(clean) == 0 {
		return 0, false
	}
	switch op {
	case aggSum:
		s := 0.0
		for _, x := range clean {
			s += x
		}
		return s, true
	case aggAvg:
		s := 0.0
		for _, x := range clean {
			s += x
		}
		return s / float64(len(clean)), true
	case aggMax:
		m := clean[0]
		for _, x := range clean[1:] {
			if x > m {
				m = x
			}
		}
		return m, true
	case aggMin:
		m := clean[0]
		for _, x := range clean[1:] {
			if x < m {
				m = x
			}
		}
		return m, true
	case aggMedian:
		if len(clean) == 1 {
			return clean[0], true
		}
		sort.Float64s(clean)
		nn := len(clean)
		if nn%2 == 1 {
			return clean[nn/2], true
		}
		return (clean[nn/2-1] + clean[nn/2]) * 0.5, true
	case aggStdev:
		if len(clean) == 1 {
			return 0, true
		}
		avg := 0.0
		for _, x := range clean {
			avg += x
		}
		avg /= float64(len(clean))
		sd := 0.0
		for _, x := range clean {
			sd += (x - avg) * (x - avg)
		}
		return math.Sqrt(sd / float64(len(clean))), true
	}
	return 0, false
}

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

// mapResult applies f elementwise, preserving the result shape.
func mapResult(r exprResult, f func(float64) float64) exprResult {
	if r.perSample {
		out := make([][]float64, len(r.values))
		for i, v := range r.values {
			if v == nil {
				continue
			}
			nv := make([]float64, len(v))
			for j, x := range v {
				nv[j] = f(x)
			}
			out[i] = nv
		}
		return exprResult{perSample: true, values: out}
	}
	nv := make([]float64, len(r.site))
	for j, x := range r.site {
		nv[j] = f(x)
	}
	return exprResult{perSample: false, site: nv}
}

// combineResults applies a binary op between two results, broadcasting a
// site-level scalar against a per-sample / vector operand as upstream does.
func combineResults(a, b exprResult, nsmpl int, f func(x, y float64) float64) exprResult {
	if !a.perSample && !b.perSample {
		n := len(a.site)
		if len(b.site) > n {
			n = len(b.site)
		}
		out := make([]float64, n)
		for i := 0; i < n; i++ {
			out[i] = f(at(a.site, i), at(b.site, i))
		}
		return exprResult{perSample: false, site: out}
	}
	out := make([][]float64, nsmpl)
	for i := 0; i < nsmpl; i++ {
		av := operandVec(a, i)
		bv := operandVec(b, i)
		n := len(av)
		if len(bv) > n {
			n = len(bv)
		}
		nv := make([]float64, n)
		for j := 0; j < n; j++ {
			nv[j] = f(at(av, j), at(bv, j))
		}
		out[i] = nv
	}
	return exprResult{perSample: true, values: out}
}

func operandVec(r exprResult, sample int) []float64 {
	if r.perSample {
		if sample < len(r.values) {
			return r.values[sample]
		}
		return nil
	}
	return r.site
}

// at returns the i-th value of v, broadcasting a single-element vector and
// returning missing past the end.
func at(v []float64, i int) float64 {
	if len(v) == 1 {
		return v[0]
	}
	if i < len(v) {
		return v[i]
	}
	return exprMissing
}

// parseExprFloatList parses a comma-separated numeric field into doubles,
// mapping "." entries to the expression missing sentinel (NaN).
func parseExprFloatList(s string) []float64 {
	if s == "" {
		return nil
	}
	parts := strings.Split(s, ",")
	out := make([]float64, len(parts))
	for i, p := range parts {
		if p == "." || p == "" {
			out[i] = exprMissing
			continue
		}
		f, err := strconv.ParseFloat(p, 64)
		if err != nil {
			out[i] = exprMissing
			continue
		}
		out[i] = f
	}
	return out
}
