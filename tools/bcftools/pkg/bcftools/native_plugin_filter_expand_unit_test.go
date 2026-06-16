// Binary-free unit tests for expandPluginFilterExpr and pluginExprLabel
// (native_plugin_filter_expand.go). These pin the curly-brace expansion order
// directly against the algorithm documented in the source, so the logic is
// verifiable on a clean/offline checkout with no upstream binary present.
package bcftools

import (
	"reflect"
	"testing"
)

// TestUnitExpandPluginFilterExpr covers the documented expansion semantics:
// single group, trailing comma, empty braces, no braces, empty input, and the
// unbalanced-brace error.
func TestUnitExpandPluginFilterExpr(t *testing.T) {
	tests := []struct {
		name    string
		expr    string
		want    []string
		wantErr bool
		// wantNil distinguishes a nil result (empty input) from an empty slice.
		wantNil bool
	}{
		{
			name: "empty input yields nil",
			expr: "",
			// expr == "" returns (nil, nil).
			wantNil: true,
		},
		{
			name: "no brace yields single element",
			expr: "GQ>10",
			want: []string{"GQ>10"},
		},
		{
			// One group of three: the source walks the working set backwards,
			// appends one entry per comma-element to the tail, then deletes the
			// source entry. With a single starting entry the appended order is
			// preserved: 1, 2, 3.
			name: "single group preserves element order",
			expr: "X>{1,2,3}",
			want: []string{"X>1", "X>2", "X>3"},
		},
		{
			// A trailing comma does NOT produce an empty final element: beg stops
			// once it passes the last real element (beg advances past the comma to
			// len(inner)), so "{10,}" yields only "10".
			name: "trailing comma drops empty element",
			expr: "X>{10,}",
			want: []string{"X>10"},
		},
		{
			// "{}" has empty inner text, so the inner loop never runs: no entries
			// are appended and the source entry is deleted, collapsing to an empty
			// (non-nil) slice.
			name: "empty braces collapse to empty slice",
			expr: "GQ>{}",
			want: []string{},
		},
		{
			name:    "unbalanced open brace errors",
			expr:    "X>{1,2",
			wantErr: true,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := expandPluginFilterExpr(tt.expr)
			if tt.wantErr {
				if err == nil {
					t.Fatalf("expandPluginFilterExpr(%q) = %v, want error", tt.expr, got)
				}
				return
			}
			if err != nil {
				t.Fatalf("expandPluginFilterExpr(%q) unexpected error: %v", tt.expr, err)
			}
			if tt.wantNil {
				if got != nil {
					t.Fatalf("expandPluginFilterExpr(%q) = %v, want nil", tt.expr, got)
				}
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Fatalf("expandPluginFilterExpr(%q) = %v, want %v", tt.expr, got, tt.want)
			}
		})
	}
}

// TestUnitExpandPluginFilterExprMultiGroup pins the cartesian-product ordering
// for multiple groups. The algorithm finds the FIRST '{' in each entry,
// tail-appends one entry per element, and deletes the source; the net effect is
// that the LAST group in the string varies fastest and the FIRST group's
// elements appear in reverse order. The expected vectors below were derived by
// tracing the source loop directly (and cross-checked against the algorithm),
// not by guessing at upstream behavior.
func TestUnitExpandPluginFilterExprMultiGroup(t *testing.T) {
	tests := []struct {
		name string
		expr string
		want []string
	}{
		{
			name: "two groups",
			expr: "a{1,2}b{3,4}c",
			want: []string{"a2b3c", "a2b4c", "a1b3c", "a1b4c"},
		},
		{
			name: "three-by-two groups",
			expr: "p{1,2,3}q{9,8}r",
			want: []string{"p3q9r", "p3q8r", "p2q9r", "p2q8r", "p1q9r", "p1q8r"},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := expandPluginFilterExpr(tt.expr)
			if err != nil {
				t.Fatalf("expandPluginFilterExpr(%q) unexpected error: %v", tt.expr, err)
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Fatalf("expandPluginFilterExpr(%q) = %v, want %v", tt.expr, got, tt.want)
			}
		})
	}
}

// TestUnitExpandPluginFilterPluginExprLabel covers the tab->space rewrite that
// keeps the tab-separated report parsable.
func TestUnitExpandPluginFilterPluginExprLabel(t *testing.T) {
	tests := []struct {
		in, want string
	}{
		{"GQ>10", "GQ>10"},
		{"A\tB", "A B"},
		{"\t\t", "  "},
		{"no tabs here", "no tabs here"},
	}
	for _, tt := range tests {
		if got := pluginExprLabel(tt.in); got != tt.want {
			t.Errorf("pluginExprLabel(%q) = %q, want %q", tt.in, got, tt.want)
		}
	}
}
