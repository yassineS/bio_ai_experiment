package bedwindow

// Parity tests for `bedtools window`. The upstream bedtools test corpus
// (reference_code/bedtools/test/) does NOT ship a `window/` subdirectory,
// so these are spec-driven cases derived from the upstream manual:
//
//   https://bedtools.readthedocs.io/en/latest/content/tools/window.html
//
// What is validated:
//   - default writer mode = `A<TAB>B` for each overlap;
//   - -w (symmetric expansion of B) creates new overlaps;
//   - -c (count-only) returns one row per A;
//   - -v (invert) emits A records with no B overlap;
//   - -l/-r asymmetric expansion;
//   - clipping at 0 on the low end.

import (
	"bytes"
	"strings"
	"testing"
)

// TestParity_NoExpansionDefault — default writer is A<TAB>B for each
// overlap (-wa/-wb both implicit).
func TestParity_NoExpansionDefault(t *testing.T) {
	a := strings.NewReader("chr1\t0\t10\n")
	b := strings.NewReader("chr1\t5\t15\n")
	var out bytes.Buffer
	n, err := Window(a, b, &out, Options{})
	if err != nil {
		t.Fatalf("Window: %v", err)
	}
	if n != 1 {
		t.Errorf("n = %d, want 1", n)
	}
	if !strings.Contains(out.String(), "chr1\t0\t10\tchr1\t5\t15") {
		t.Errorf("expected A<TAB>B in output, got %q", out.String())
	}
}

// TestParity_WindowExpandsB — symmetric -w 100 expansion of B should
// pull a non-touching B into the overlap set.
func TestParity_WindowExpandsB(t *testing.T) {
	a := strings.NewReader("chr1\t300\t310\n")
	b := strings.NewReader("chr1\t200\t210\n")
	var out bytes.Buffer
	n, err := Window(a, b, &out, Options{Left: 100, Right: 100})
	if err != nil {
		t.Fatalf("Window: %v", err)
	}
	if n != 1 {
		t.Errorf("n = %d, want 1 (B expanded to [100,310))", n)
	}
}

// TestParity_CountMode — -c emits A<TAB>count, one row per A even when
// there is no overlap (count = 0).
func TestParity_CountMode(t *testing.T) {
	a := strings.NewReader("chr1\t0\t10\nchr2\t0\t10\n")
	b := strings.NewReader("chr1\t5\t6\nchr1\t8\t9\n")
	var out bytes.Buffer
	if _, err := Window(a, b, &out, Options{Count: true}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	got := out.String()
	if !strings.Contains(got, "chr1\t0\t10\t2\n") {
		t.Errorf("expected chr1 count=2 line, got:\n%s", got)
	}
	if !strings.Contains(got, "chr2\t0\t10\t0\n") {
		t.Errorf("expected chr2 count=0 line, got:\n%s", got)
	}
}

// TestParity_InvertMode — -v emits only A records with NO B overlap.
func TestParity_InvertMode(t *testing.T) {
	a := strings.NewReader("chr1\t0\t10\nchr2\t0\t10\n")
	b := strings.NewReader("chr1\t5\t6\n")
	var out bytes.Buffer
	if _, err := Window(a, b, &out, Options{Invert: true}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	got := strings.TrimSpace(out.String())
	if !strings.Contains(got, "chr2\t0\t10") {
		t.Errorf("expected chr2 in invert output, got %q", got)
	}
	if strings.Contains(got, "chr1\t0\t10") {
		t.Errorf("chr1 has overlap; should NOT appear in invert output: %q", got)
	}
}

// TestParity_AsymmetricExpansion — -l 0 -r 100: B is extended only to
// the right, so a B that is upstream of A still doesn't overlap.
func TestParity_AsymmetricExpansion(t *testing.T) {
	a := strings.NewReader("chr1\t300\t310\n")
	b := strings.NewReader("chr1\t100\t150\n")
	var out bytes.Buffer
	n, err := Window(a, b, &out, Options{Left: 0, Right: 100})
	if err != nil {
		t.Fatalf("Window: %v", err)
	}
	if n != 0 {
		t.Errorf("n = %d, want 0 (right-only expansion shouldn't pull in upstream B)", n)
	}
}

// TestParity_LowClipping — -l 1000 on a B at start=50 should clip to 0
// rather than going negative.
func TestParity_LowClipping(t *testing.T) {
	a := strings.NewReader("chr1\t10\t20\n")
	b := strings.NewReader("chr1\t50\t60\n")
	var out bytes.Buffer
	if _, err := Window(a, b, &out, Options{Left: 1000, Right: 0}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	if !strings.Contains(out.String(), "chr1\t10\t20") {
		t.Errorf("expected overlap after low-clipping, got %q", out.String())
	}
}
