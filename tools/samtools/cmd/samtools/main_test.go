package main

import (
	"reflect"
	"testing"
)

// TestExpandShortAA pins the `-aa` → `--all-positions-all-chroms`
// rewrite that bridges upstream samtools' fused-short syntax to Go's
// flag package. Previously a post-parse scan tried to handle it but
// `flag.Parse` rejected `-aa` first, leaving the headline `-aa`
// behaviour broken at the CLI even though library-level tests passed.
func TestExpandShortAA(t *testing.T) {
	cases := []struct {
		name string
		in   []string
		want []string
	}{
		{
			name: "no_aa",
			in:   []string{"-a", "in.bam"},
			want: []string{"-a", "in.bam"},
		},
		{
			name: "rewrite_lone",
			in:   []string{"-aa", "in.bam"},
			want: []string{"--all-positions-all-chroms", "in.bam"},
		},
		{
			name: "rewrite_mid",
			in:   []string{"-r", "chr1", "-aa", "-f", "ref.fa", "in.bam"},
			want: []string{"-r", "chr1", "--all-positions-all-chroms", "-f", "ref.fa", "in.bam"},
		},
		{
			name: "dash_dash_stops_rewriting",
			in:   []string{"--", "-aa"}, // positional arg literally named "-aa"
			want: []string{"--", "-aa"},
		},
		{
			name: "preserves_double_a_long",
			in:   []string{"-a", "-a", "in.bam"},
			want: []string{"-a", "-a", "in.bam"},
		},
		{
			name: "empty",
			in:   nil,
			want: []string{},
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := expandShortAA(tc.in)
			if !reflect.DeepEqual(got, tc.want) {
				t.Fatalf("expandShortAA(%v) = %v, want %v", tc.in, got, tc.want)
			}
		})
	}
}

// TestParseOutputFmtOptions checks the `--output-fmt-option` KEY=VALUE
// parser: a valid qbin value is extracted, last-wins on repeats, and a
// malformed entry or unknown key is an error.
func TestParseOutputFmtOptions(t *testing.T) {
	if got, err := parseOutputFmtOptions(nil); err != nil || got != "" {
		t.Errorf("parseOutputFmtOptions(nil) = (%q, %v), want (\"\", nil)", got, err)
	}
	if got, err := parseOutputFmtOptions(multiString{"qbin=8"}); err != nil || got != "8" {
		t.Errorf("parseOutputFmtOptions(qbin=8) = (%q, %v), want (\"8\", nil)", got, err)
	}
	// Last value wins for a repeated key.
	if got, err := parseOutputFmtOptions(multiString{"qbin=8", "qbin=4"}); err != nil || got != "4" {
		t.Errorf("repeated qbin = (%q, %v), want (\"4\", nil)", got, err)
	}
	if _, err := parseOutputFmtOptions(multiString{"qbin"}); err == nil {
		t.Error("parseOutputFmtOptions(\"qbin\") should fail (no '=')")
	}
	if _, err := parseOutputFmtOptions(multiString{"bogus=1"}); err == nil {
		t.Error("parseOutputFmtOptions with unknown key should fail")
	}
}
