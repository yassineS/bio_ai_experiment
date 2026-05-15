package bedtag

import (
	"bytes"
	"strings"
	"testing"
)

func sourceFromString(name, s string) Source {
	return Source{Name: name, Reader: strings.NewReader(s)}
}

func TestTag_Basic_NameColumn(t *testing.T) {
	a := strings.NewReader(
		"chr1\t0\t100\tregion1\n" +
			"chr1\t200\t300\tregion2\n",
	)
	b := sourceFromString("b.bed",
		"chr1\t10\t20\tpeak1\n"+
			"chr1\t50\t60\tpeak2\n"+
			"chr1\t250\t260\tpeak3\n",
	)
	var out bytes.Buffer
	n, err := Tag(a, []Source{b}, &out, Options{})
	if err != nil {
		t.Fatalf("Tag: %v", err)
	}
	if n != 2 {
		t.Errorf("n = %d, want 2", n)
	}
	got := out.String()
	if !strings.Contains(got, "region1") || !strings.Contains(got, "peak1,peak2") {
		t.Errorf("expected region1 tagged with peak1,peak2; got:\n%s", got)
	}
	if !strings.Contains(got, "region2") || !strings.Contains(got, "peak3") {
		t.Errorf("expected region2 tagged with peak3; got:\n%s", got)
	}
}

func TestTag_NoOverlapEmptyTagColumn(t *testing.T) {
	a := strings.NewReader("chr1\t0\t10\tlonely\n")
	b := sourceFromString("b.bed", "chr2\t0\t10\tother\n")
	var out bytes.Buffer
	if _, err := Tag(a, []Source{b}, &out, Options{}); err != nil {
		t.Fatalf("Tag: %v", err)
	}
	// Trailing empty column => the line ends with a tab.
	if !strings.HasSuffix(strings.TrimRight(out.String(), "\n"), "\t") {
		t.Errorf("expected trailing empty tag column, got: %q", out.String())
	}
}

func TestTag_StrandSpec_Filters(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\ta\t0\t+\n")
	b := sourceFromString("b.bed",
		"chr1\t10\t20\tplus\t0\t+\n"+
			"chr1\t50\t60\tminus\t0\t-\n",
	)
	var out bytes.Buffer
	if _, err := Tag(a, []Source{b}, &out, Options{StrandSpec: true}); err != nil {
		t.Fatalf("Tag: %v", err)
	}
	if !strings.Contains(out.String(), "plus") {
		t.Errorf("expected plus tag, got: %q", out.String())
	}
	if strings.Contains(out.String(), "minus") {
		t.Errorf("did not want minus tag with -s, got: %q", out.String())
	}
}

func TestTag_InverseStrand_Filters(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\ta\t0\t+\n")
	b := sourceFromString("b.bed",
		"chr1\t10\t20\tplus\t0\t+\n"+
			"chr1\t50\t60\tminus\t0\t-\n",
	)
	var out bytes.Buffer
	if _, err := Tag(a, []Source{b}, &out, Options{InverseStrand: true}); err != nil {
		t.Fatalf("Tag: %v", err)
	}
	if !strings.Contains(out.String(), "minus") {
		t.Errorf("expected minus tag, got: %q", out.String())
	}
	if strings.Contains(out.String(), "plus") {
		t.Errorf("did not want plus tag with -S, got: %q", out.String())
	}
}

func TestTag_FractionA(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\tregion\n")
	b := sourceFromString("b.bed", "chr1\t0\t10\tsmall\n")
	// 10bp / 100bp = 0.1 < 0.5 -> should NOT tag.
	var out bytes.Buffer
	if _, err := Tag(a, []Source{b}, &out, Options{FractionA: 0.5}); err != nil {
		t.Fatalf("Tag: %v", err)
	}
	if strings.Contains(out.String(), "small") {
		t.Errorf("expected no tag (FractionA filter), got: %q", out.String())
	}
}

func TestTag_NamesReplacement(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\n")
	b1 := sourceFromString("b1.bed", "chr1\t10\t20\torig1\n")
	b2 := sourceFromString("b2.bed", "chr1\t30\t40\torig2\n")
	var out bytes.Buffer
	opts := Options{Names: []string{"FILE1", "FILE2"}}
	if _, err := Tag(a, []Source{b1, b2}, &out, opts); err != nil {
		t.Fatalf("Tag: %v", err)
	}
	got := out.String()
	if !strings.Contains(got, "FILE1") || !strings.Contains(got, "FILE2") {
		t.Errorf("expected FILE1,FILE2 tags, got: %q", got)
	}
	if strings.Contains(got, "orig1") || strings.Contains(got, "orig2") {
		t.Errorf("Names should override original names; got: %q", got)
	}
}

func TestTag_LabelsPrefix(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\n")
	b := sourceFromString("peaks.bed", "chr1\t10\t20\tpeak1\n")
	var out bytes.Buffer
	if _, err := Tag(a, []Source{b}, &out, Options{Labels: true}); err != nil {
		t.Fatalf("Tag: %v", err)
	}
	if !strings.Contains(out.String(), "peaks.bed=peak1") {
		t.Errorf("expected `peaks.bed=peak1` prefix, got: %q", out.String())
	}
}

func TestTag_TagColumn_NonName(t *testing.T) {
	// Use B's column 5 (score) as the tag.
	a := strings.NewReader("chr1\t0\t100\n")
	b := sourceFromString("b.bed", "chr1\t10\t20\tname\t999\t+\n")
	var out bytes.Buffer
	if _, err := Tag(a, []Source{b}, &out, Options{TagColumn: 5}); err != nil {
		t.Fatalf("Tag: %v", err)
	}
	if !strings.Contains(out.String(), "999") {
		t.Errorf("expected score 999 as tag, got: %q", out.String())
	}
}

func TestTag_MultiSource_OrderPreserved(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\n")
	b1 := sourceFromString("b1.bed", "chr1\t10\t20\tFROM_B1\n")
	b2 := sourceFromString("b2.bed", "chr1\t30\t40\tFROM_B2\n")
	var out bytes.Buffer
	if _, err := Tag(a, []Source{b1, b2}, &out, Options{}); err != nil {
		t.Fatalf("Tag: %v", err)
	}
	// Tags should appear in the order of the source files.
	idx1 := strings.Index(out.String(), "FROM_B1")
	idx2 := strings.Index(out.String(), "FROM_B2")
	if idx1 < 0 || idx2 < 0 || idx1 > idx2 {
		t.Errorf("source order not preserved, got: %q", out.String())
	}
}

func TestTag_Mismatched_NamesLength_Error(t *testing.T) {
	a := strings.NewReader("chr1\t0\t100\n")
	b := sourceFromString("b.bed", "chr1\t0\t10\tx\n")
	if _, err := Tag(a, []Source{b}, &bytes.Buffer{}, Options{Names: []string{"X", "Y"}}); err == nil {
		t.Error("expected error for mismatched Names length")
	}
}

func TestTag_StrandSpec_And_Inverse_Error(t *testing.T) {
	if _, err := Tag(strings.NewReader(""), nil, &bytes.Buffer{},
		Options{StrandSpec: true, InverseStrand: true}); err == nil {
		t.Error("expected error: StrandSpec and InverseStrand are mutually exclusive")
	}
}
