// Curly-brace multi-threshold expansion for the native stats plugins'
// -i/--include / -e/--exclude expressions, matching the parse_filters() routine
// that smpl-stats.c, indel-stats.c and trio-stats.c repeat verbatim. Upstream
// lets a single expression scan a range of thresholds at once, e.g.
// `-i 'GQ>{10,20,30}'`, which is expanded into one filter expression per
// comma-separated element (the braces replaced by that element), each tallied
// into its own FLT* report section. Multiple `{...}` groups combine as a
// cartesian product, with the exact ordering the C code produces.
package bcftools

import (
	"fmt"
	"strings"
)

// expandPluginFilterExpr expands the curly-brace lists in a single -i/-e filter
// expression into the list of concrete filter expressions, replicating the
// upstream parse_filters() algorithm (plugins/{smpl,indel,trio}-stats.c)
// byte-for-byte, including its expansion order for one or more `{...}` groups.
//
// An empty expr (no -i/-e given) yields a nil slice. An expression with no `{`
// yields a single-element slice ([]string{expr}). An expression whose braces
// collapse to nothing (e.g. "GQ>{}") can yield an empty, non-nil slice, exactly
// as upstream where nflt_str becomes 0 and the plugin falls back to the single
// "all" filter. An unmatched `{` (no `}`) is a hard error, matching upstream's
// "Could not parse the expression" error.
//
// The implementation mirrors the C loop verbatim: it walks the working set
// backwards (for i=n-1; i>=0; i--), and when it finds a brace group in entry i
// it appends one new entry per comma-separated element to the END of the set
// and then deletes entry i (shifting later entries left by one). That tail-
// append + in-place-delete is what produces upstream's exact FLT ordering for
// nested/multiple groups, so it is reproduced rather than a simpler in-place
// substitution.
func expandPluginFilterExpr(expr string) ([]string, error) {
	if expr == "" {
		return nil, nil
	}
	fltStr := []string{expr}
	for {
		expanded := false
		for i := len(fltStr) - 1; i >= 0; i-- {
			s := fltStr[i]
			begPos := strings.IndexByte(s, '{')
			if begPos < 0 {
				continue
			}
			endPos := strings.IndexByte(s[begPos+1:], '}')
			if endPos < 0 {
				return nil, fmt.Errorf("Could not parse the expression: %s", expr)
			}
			endPos += begPos + 1 // index of '}' in s

			prefix := s[:begPos]
			suffix := s[endPos+1:]
			inner := s[begPos+1 : endPos] // text between the braces

			// Split the inner text on commas, appending one new expression per
			// element to the tail. Upstream advances beg past each comma and
			// stops when beg reaches the closing brace, so a trailing comma
			// (e.g. "{10,}") does NOT yield an empty final element, and "{}"
			// yields no elements at all.
			beg := 0
			for beg < len(inner) {
				mid := beg
				for mid < len(inner) && inner[mid] != ',' {
					mid++
				}
				fltStr = append(fltStr, prefix+inner[beg:mid]+suffix)
				beg = mid + 1
			}

			// Delete entry i (shift later entries left by one), mirroring the C
			// memmove + nflt_str-- after the tail append.
			fltStr = append(fltStr[:i], fltStr[i+1:]...)
			expanded = true
		}
		if !expanded {
			break
		}
	}
	return fltStr, nil
}

// pluginExprLabel renders the DEF-line label for an expanded filter expression,
// replacing tab characters with spaces so the tab-separated report stays
// parsable — matching the tab-to-space rewrite the stats plugins apply to
// flt->expr in init_data().
func pluginExprLabel(expr string) string {
	return strings.ReplaceAll(expr, "\t", " ")
}
