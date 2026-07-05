package samtools

import (
	"bytes"
	"strings"
	"testing"
)

// TestDict_FieldOrder pins the @SQ field order to upstream `samtools dict`:
// SN, LN, M5, [AN], [UR], [AS], [SP]. A regression here (e.g. AS emitted
// before UR) breaks byte-exact parity on any multi-field invocation.
func TestDict_FieldOrder(t *testing.T) {
	in := strings.NewReader(">chr1\nACGTACGT\n")
	var out bytes.Buffer
	if err := Dict(in, &out, DictOptions{
		Assembly: "GRCh38",
		Species:  "Homo sapiens",
		URI:      "file:///ref.fa",
		NoHeader: true,
	}); err != nil {
		t.Fatalf("Dict: %v", err)
	}
	got := strings.TrimRight(out.String(), "\n")
	// Upstream emits M5 before UR before AS before SP.
	order := []string{"M5:", "UR:file:///ref.fa", "AS:GRCh38", "SP:Homo sapiens"}
	prev := -1
	for _, tok := range order {
		idx := strings.Index(got, tok)
		if idx < 0 {
			t.Fatalf("missing %q in %q", tok, got)
		}
		if idx < prev {
			t.Fatalf("field %q out of order in %q", tok, got)
		}
		prev = idx
	}
}

// TestDict_AliasChrPrefix verifies -A/--alias reproduces upstream's
// add/remove-"chr" semantics as a single comma-joined AN tag, including
// the chrM/chrMT mitochondrial synonyms, emitted before UR/AS/SP.
func TestDict_AliasChrPrefix(t *testing.T) {
	cases := []struct {
		name string
		want string
	}{
		{"chr1", "AN:1"},
		{"1", "AN:chr1"},
		{"chrM", "AN:M,chrMT,MT"},
		{"M", "AN:chrM,chrMT,MT"},
		{"chrMT", "AN:MT,chrM,M"},
		{"MT", "AN:chrMT,chrM,M"},
		{"chrX", "AN:X"},
	}
	for _, tc := range cases {
		in := strings.NewReader(">" + tc.name + "\nACGT\n")
		var out bytes.Buffer
		if err := Dict(in, &out, DictOptions{AliasFromHeader: true, NoHeader: true}); err != nil {
			t.Fatalf("Dict(%s): %v", tc.name, err)
		}
		got := strings.TrimRight(out.String(), "\n")
		if !strings.Contains(got, "\t"+tc.want+"\t") && !strings.HasSuffix(got, "\t"+tc.want) {
			t.Errorf("alias for %q: got %q, want field %q", tc.name, got, tc.want)
		}
		// AN must precede any UR (none here) — assert it sits right after M5.
		anIdx := strings.Index(got, "\tAN:")
		m5Idx := strings.Index(got, "\tM5:")
		if anIdx < 0 || m5Idx < 0 || anIdx < m5Idx {
			t.Errorf("alias for %q: AN not positioned after M5 in %q", tc.name, got)
		}
	}
}
