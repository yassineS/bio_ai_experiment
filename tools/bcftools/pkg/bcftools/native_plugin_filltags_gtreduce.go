// Genotype-count reduction nodes for the fill-tags expression evaluator:
// F_MISSING, F_PASS, N_PASS, N_MISSING. These mirror filter.c's func_npass:
// the inner expression is a per-sample condition (e.g. GT="mis", GT="het",
// FMT/DP>10); N_PASS returns the number of active samples that pass, F_PASS
// returns that count divided by the number of active samples. N_MISSING /
// F_MISSING are the genotype-missingness counterparts (npass over samples
// whose genotype has a missing allele). The active-sample set honours the
// -S population restriction passed in via the evaluation context.
package bcftools

import (
	"fmt"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// gtReduceNode evaluates a per-sample condition and reduces it to a single
// site-level count (N_*) or fraction (F_*).
type gtReduceNode struct {
	name     string // F_MISSING, F_PASS, N_PASS, N_MISSING
	fraction bool   // true for the F_* forms
	missing  bool   // true for the *_MISSING forms (count missing genotypes)
	cond     *Filter
}

// newGTReduceNode compiles a genotype-count reduction. For the *_MISSING forms
// the inner text (if any) is ignored: upstream defines F_MISSING as the
// synonym F_PASS(GT="mis"), so the reduction simply counts samples with a
// missing genotype. For F_PASS / N_PASS the inner text is the per-sample
// condition, compiled with the shared native filter engine.
func newGTReduceNode(name, inner string, hdr *vcf.Header) (exprNode, error) {
	n := &gtReduceNode{name: name}
	switch name {
	case "F_MISSING":
		n.fraction, n.missing = true, true
	case "N_MISSING":
		n.fraction, n.missing = false, true
	case "F_PASS":
		n.fraction = true
	case "N_PASS":
		n.fraction = false
	default:
		return nil, fmt.Errorf("fill-tags: unsupported reduction %q", name)
	}
	if !n.missing {
		trimmed := strings.TrimSpace(inner)
		if trimmed == "" {
			return nil, fmt.Errorf("fill-tags: %s() requires a condition", name)
		}
		f, err := CompileFilterWithHeader(trimmed, hdr)
		if err != nil {
			return nil, fmt.Errorf("fill-tags: %s(%s): %w", name, inner, err)
		}
		n.cond = f
	} else {
		// F_MISSING accepts the explicit GT="mis" form too; if a condition was
		// given that is not the missing shorthand, compile it as a filter.
		trimmed := strings.TrimSpace(inner)
		if trimmed != "" && !isMissingShorthand(trimmed) {
			f, err := CompileFilterWithHeader(trimmed, hdr)
			if err != nil {
				return nil, fmt.Errorf("fill-tags: %s(%s): %w", name, inner, err)
			}
			n.cond = f
			n.missing = false
		}
	}
	return n, nil
}

// isMissingShorthand reports whether the condition is the GT="mis" shorthand
// that F_MISSING/N_MISSING expand from.
func isMissingShorthand(s string) bool {
	s = strings.ToLower(strings.ReplaceAll(s, " ", ""))
	return s == `gt="mis"` || s == `gt='mis'` || s == "fmt/gt=\"mis\"" || s == "format/gt=\"mis\""
}

func (n *gtReduceNode) eval(ctx *exprCtx) exprResult {
	nsmpl := 0
	npass := 0
	var mask []bool
	if n.cond != nil {
		_, mask = n.cond.EvalSamples(ctx.v)
	}
	for i := 0; i < ctx.nsmpl; i++ {
		if ctx.usmpl != nil && !ctx.usmpl[i] {
			continue
		}
		nsmpl++
		if n.missing {
			gt, ok := sampleGT(ctx.v, i)
			if !ok || gt.nMissing() > 0 {
				npass++
			}
			continue
		}
		if mask != nil && i < len(mask) && mask[i] {
			npass++
		}
	}
	var val float64
	if n.fraction {
		if nsmpl != 0 {
			val = float64(npass) / float64(nsmpl)
		}
	} else {
		val = float64(npass)
	}
	return exprResult{perSample: false, site: []float64{val}}
}
