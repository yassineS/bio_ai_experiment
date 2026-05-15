package bedwindow

import (
	"bytes"
	"strings"
	"testing"
)

func TestWindow_BasicNoExpansion(t *testing.T) {
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
	// Default mode = A<TAB>B
	if !strings.Contains(out.String(), "chr1\t0\t10\tchr1\t5\t15") {
		t.Errorf("expected A<TAB>B, got: %q", out.String())
	}
}

func TestWindow_ExpandsB(t *testing.T) {
	// A=[100,110), B=[200,210). Without expansion, no overlap. With -l 100,
	// B becomes [100,210), which overlaps A.
	a := strings.NewReader("chr1\t100\t110\n")
	b := strings.NewReader("chr1\t200\t210\n")
	var out bytes.Buffer
	n, err := Window(a, b, &out, Options{Left: 100, Right: 0})
	if err != nil {
		t.Fatalf("Window: %v", err)
	}
	if n != 1 {
		t.Errorf("n = %d, want 1 with left=100", n)
	}
}

func TestWindow_BothSidesExtension(t *testing.T) {
	// A=[300,310). B=[200,210). With -w 100, B becomes [100,310) — overlaps A.
	a := strings.NewReader("chr1\t300\t310\n")
	b := strings.NewReader("chr1\t200\t210\n")
	var out bytes.Buffer
	if _, err := Window(a, b, &out, Options{Left: 100, Right: 100}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	if !strings.Contains(out.String(), "chr1\t300\t310") {
		t.Errorf("expected overlap output, got: %q", out.String())
	}
}

func TestWindow_NoExpansionNoHit(t *testing.T) {
	a := strings.NewReader("chr1\t100\t110\n")
	b := strings.NewReader("chr1\t200\t210\n")
	var out bytes.Buffer
	n, err := Window(a, b, &out, Options{})
	if err != nil {
		t.Fatalf("Window: %v", err)
	}
	if n != 0 {
		t.Errorf("n = %d, want 0", n)
	}
	if out.Len() != 0 {
		t.Errorf("expected empty output, got: %q", out.String())
	}
}

func TestWindow_Count(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\n")
	b := strings.NewReader("chr1\t10\t20\nchr1\t50\t60\n")
	var out bytes.Buffer
	if _, err := Window(a, b, &out, Options{Count: true}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	if !strings.Contains(out.String(), "chr1\t0\t100\t2") {
		t.Errorf("expected count=2, got: %q", out.String())
	}
}

func TestWindow_Invert(t *testing.T) {
	a := strings.NewReader("chr1\t0\t10\nchr1\t100\t110\n")
	b := strings.NewReader("chr1\t5\t15\n")
	var out bytes.Buffer
	if _, err := Window(a, b, &out, Options{Invert: true}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	if !strings.Contains(out.String(), "chr1\t100\t110") {
		t.Errorf("expected the un-overlapping record, got: %q", out.String())
	}
	if strings.Contains(out.String(), "chr1\t0\t10") {
		t.Errorf("did not want the overlapping record, got: %q", out.String())
	}
}

func TestWindow_WriteA(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\tA1\n")
	b := strings.NewReader("chr1\t10\t20\tB1\n")
	var out bytes.Buffer
	if _, err := Window(a, b, &out, Options{WriteA: true}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	got := strings.TrimRight(out.String(), "\n")
	if got != "chr1\t0\t100\tA1" {
		t.Errorf("expected A only, got: %q", got)
	}
}

func TestWindow_WriteB(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\tA1\n")
	b := strings.NewReader("chr1\t10\t20\tB1\n")
	var out bytes.Buffer
	if _, err := Window(a, b, &out, Options{WriteB: true}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	got := strings.TrimRight(out.String(), "\n")
	if got != "chr1\t10\t20\tB1" {
		t.Errorf("expected B only, got: %q", got)
	}
}

func TestWindow_StrandSpec(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\tA\t0\t+\n")
	b := strings.NewReader(
		"chr1\t10\t20\tplus\t0\t+\n" +
			"chr1\t50\t60\tminus\t0\t-\n",
	)
	var out bytes.Buffer
	if _, err := Window(a, b, &out, Options{StrandSpec: true, WriteB: true}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	if !strings.Contains(out.String(), "plus") || strings.Contains(out.String(), "minus") {
		t.Errorf("expected only plus, got: %q", out.String())
	}
}

func TestWindow_ClipsLeftAtZero(t *testing.T) {
	// B at [0,10), Left=100. Expanded becomes [-100,10), clipped to [0,10).
	// A at [5,15) overlaps the original B; result should still be one hit.
	a := strings.NewReader("chr1\t5\t15\n")
	b := strings.NewReader("chr1\t0\t10\n")
	var out bytes.Buffer
	n, err := Window(a, b, &out, Options{Left: 100})
	if err != nil {
		t.Fatalf("Window: %v", err)
	}
	if n != 1 {
		t.Errorf("n = %d, want 1", n)
	}
}

func TestWindow_StrandAndInverse_Error(t *testing.T) {
	if _, err := Window(strings.NewReader(""), strings.NewReader(""),
		&bytes.Buffer{}, Options{StrandSpec: true, InverseStrand: true}); err == nil {
		t.Error("expected error: -sm and -sw are mutually exclusive")
	}
}
