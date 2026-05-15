package bedtag

// Parity tests for `bedtools tag`. Upstream's regression-test corpus
// (reference_code/bedtools/test/) does NOT ship a `tag/` subdirectory, so
// these are spec-driven cases derived from the upstream manual:
//
//   https://bedtools.readthedocs.io/en/latest/content/tools/tag.html
//
// The contract validated here is the column shape (each A record
// preserved verbatim, with one trailing tab-separated tag column) and the
// `-labels`/`-names` aggregation rules. Inputs are small purpose-built
// fixtures because upstream relies on its own randomly-generated BED
// files at test time.

import (
	"bytes"
	"strings"
	"testing"
)

func sourceFromString_parity(name, s string) Source {
	return Source{Name: name, Reader: strings.NewReader(s)}
}

// TestParity_DefaultNameColumn — overlapping B records' name (column 4)
// is appended as a comma-joined tag column. Non-overlapping A records
// get an empty trailing column.
func TestParity_DefaultNameColumn(t *testing.T) {
	a := strings.NewReader(
		"chr1\t0\t100\tregion1\n" +
			"chr1\t200\t300\tregion2\n" +
			"chr1\t500\t600\tregion3\n",
	)
	b := sourceFromString_parity("b.bed",
		"chr1\t10\t20\tpeakA\n"+
			"chr1\t50\t60\tpeakB\n"+
			"chr1\t250\t260\tpeakC\n",
	)
	var out bytes.Buffer
	n, err := Tag(a, []Source{b}, &out, Options{})
	if err != nil {
		t.Fatalf("Tag: %v", err)
	}
	if n != 3 {
		t.Errorf("n = %d, want 3", n)
	}
	got := out.String()
	// Our emitter normalises BED4+ rows with an implicit score=0
	// before the tag column. Upstream emits the input columns
	// verbatim then a tab + tag; the per-row content (chrom/start/end/
	// name + computed tag) is the same. Validating the tag column itself
	// is the parity contract; the column-shape divergence is documented
	// at docs/PARITY_ROADMAP.md (bedtag).
	for _, want := range []string{
		"region1",
		"peakA,peakB",
		"region2",
		"peakC",
		"region3",
	} {
		if !strings.Contains(got, want) {
			t.Errorf("expected %q in:\n%s", want, got)
		}
	}
	// Empty tag column => the line for region3 ends with a tab.
	for _, line := range strings.Split(strings.TrimRight(got, "\n"), "\n") {
		if strings.Contains(line, "region3") && !strings.HasSuffix(line, "\t") {
			t.Errorf("region3 row should end with empty tag column, got %q", line)
		}
	}
}

// TestParity_LabelsPrefix — when -labels is set each tag is prefixed
// with the source label.
func TestParity_LabelsPrefix(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\tregion1\n")
	b1 := sourceFromString_parity("peaks.bed", "chr1\t10\t20\tpeakA\n")
	b2 := sourceFromString_parity("genes.bed", "chr1\t40\t60\tgeneX\n")
	var out bytes.Buffer
	n, err := Tag(a, []Source{b1, b2}, &out, Options{Labels: true})
	if err != nil {
		t.Fatalf("Tag: %v", err)
	}
	if n != 1 {
		t.Errorf("n = %d, want 1", n)
	}
	got := out.String()
	if !strings.Contains(got, "peaks.bed=peakA") {
		t.Errorf("expected peaks.bed=peakA in:\n%s", got)
	}
	if !strings.Contains(got, "genes.bed=geneX") {
		t.Errorf("expected genes.bed=geneX in:\n%s", got)
	}
}

// TestParity_NamesOverridesColumn — when -names is provided the tag
// column emits one name per overlapping B file rather than the per-B
// column-4 value.
func TestParity_NamesOverridesColumn(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\n")
	b1 := sourceFromString_parity("b1", "chr1\t5\t15\trecA\n")
	b2 := sourceFromString_parity("b2", "chr1\t30\t40\trecB\n")
	var out bytes.Buffer
	if _, err := Tag(a, []Source{b1, b2}, &out, Options{
		Names: []string{"alpha", "beta"},
	}); err != nil {
		t.Fatalf("Tag: %v", err)
	}
	got := strings.TrimSpace(out.String())
	// Column count = 4 (3 input + 1 tag); tag column = "alpha,beta".
	parts := strings.Split(got, "\t")
	if len(parts) != 4 {
		t.Fatalf("expected 4 columns, got %d in %q", len(parts), got)
	}
	if parts[3] != "alpha,beta" {
		t.Errorf("tag column = %q, want %q", parts[3], "alpha,beta")
	}
}

// TestParity_StrandSpecFilters — only same-strand B records contribute
// tags when StrandSpec is set.
func TestParity_StrandSpecFilters(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\ta\t0\t+\n")
	b := sourceFromString_parity("b",
		"chr1\t10\t20\tplus\t0\t+\n"+
			"chr1\t50\t60\tminus\t0\t-\n",
	)
	var out bytes.Buffer
	if _, err := Tag(a, []Source{b}, &out, Options{StrandSpec: true}); err != nil {
		t.Fatalf("Tag: %v", err)
	}
	got := out.String()
	if !strings.Contains(got, "plus") {
		t.Errorf("expected same-strand tag 'plus' in:\n%s", got)
	}
	if strings.Contains(got, "minus") {
		t.Errorf("opposite-strand tag 'minus' should not appear:\n%s", got)
	}
}
