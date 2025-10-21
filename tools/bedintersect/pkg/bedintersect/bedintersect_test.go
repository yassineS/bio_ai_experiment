package bedintersect

import (
	"bytes"
	"strings"
	"testing"
)

func TestIntersectBasic(t *testing.T) {
	fileA := `chr1	100	200
chr1	300	400`

	fileB := `chr1	150	250
chr1	350	450`

	expected := `chr1	150	200
chr1	350	400
`

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{MinOverlap: 1})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 overlaps, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectNoOverlap(t *testing.T) {
	fileA := `chr1	100	200
chr1	300	400`

	fileB := `chr1	250	280
chr1	450	500`

	expected := ``

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{MinOverlap: 1})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 0 {
		t.Errorf("Expected 0 overlaps, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectReportNoOverlap(t *testing.T) {
	fileA := `chr1	100	200
chr1	300	400
chr1	500	600`

	fileB := `chr1	150	250`

	expected := `chr1	300	400
chr1	500	600
`

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{MinOverlap: 1, NoOverlap: true})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 non-overlapping intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectCount(t *testing.T) {
	fileA := `chr1	100	500`

	fileB := `chr1	150	200
chr1	250	300
chr1	400	450`

	expected := `chr1	100	500	3
`

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{MinOverlap: 1, Count: true})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 1 {
		t.Errorf("Expected 1 output, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectWriteA(t *testing.T) {
	fileA := `chr1	100	200`

	fileB := `chr1	150	250`

	expected := `chr1	100	200
`

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{MinOverlap: 1, WriteA: true})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 1 {
		t.Errorf("Expected 1 overlap, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectWriteB(t *testing.T) {
	fileA := `chr1	100	200`

	fileB := `chr1	150	250`

	expected := `chr1	150	250
`

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{MinOverlap: 1, WriteB: true})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 1 {
		t.Errorf("Expected 1 overlap, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectMinOverlap(t *testing.T) {
	fileA := `chr1	100	200`

	fileB := `chr1	199	250`

	expected := ``

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	// Only 1bp overlap, require 10bp minimum
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{MinOverlap: 10})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 0 {
		t.Errorf("Expected 0 overlaps (below minimum), got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectMultipleChromosomes(t *testing.T) {
	fileA := `chr1	100	200
chr2	100	200`

	fileB := `chr1	150	250
chr2	150	250`

	expected := `chr1	150	200
chr2	150	200
`

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{MinOverlap: 1})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 overlaps, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectWithStats(t *testing.T) {
	fileA := `chr1	100	200
chr1	300	400
chr1	500	600`

	fileB := `chr1	150	250
chr1	350	450`

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	stats, err := IntersectWithStats(readerA, readerB, &buf, IntersectOptions{MinOverlap: 1})
	if err != nil {
		t.Fatalf("IntersectWithStats failed: %v", err)
	}

	if stats.IntervalsA != 3 {
		t.Errorf("Expected 3 A intervals, got %d", stats.IntervalsA)
	}

	if stats.IntervalsB != 2 {
		t.Errorf("Expected 2 B intervals, got %d", stats.IntervalsB)
	}

	if stats.IntervalsAHit != 2 {
		t.Errorf("Expected 2 A intervals with hits, got %d", stats.IntervalsAHit)
	}

	if stats.IntervalsAMiss != 1 {
		t.Errorf("Expected 1 A interval with no hits, got %d", stats.IntervalsAMiss)
	}

	if stats.Overlaps != 2 {
		t.Errorf("Expected 2 total overlaps, got %d", stats.Overlaps)
	}
}

func TestIntersectFractionA(t *testing.T) {
	fileA := `chr1	100	200`

	fileB := `chr1	150	200`

	expected := ``

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	// 50bp overlap / 100bp length = 0.5, require 0.75
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{MinOverlap: 1, FractionA: 0.75})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 0 {
		t.Errorf("Expected 0 overlaps (below fraction threshold), got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectFractionB(t *testing.T) {
	fileA := `chr1	100	200`

	fileB := `chr1	150	300`

	expected := `chr1	150	200
`

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	// 50bp overlap / 150bp B length = 0.33, require 0.25
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{MinOverlap: 1, FractionB: 0.25})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 1 {
		t.Errorf("Expected 1 overlap, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectReciprocal(t *testing.T) {
	fileA := `chr1	100	200`

	fileB := `chr1	150	300`

	expected := ``

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	// Overlap: 50bp
	// Fraction of A: 50/100 = 0.5 (meets 0.5 threshold)
	// Fraction of B: 50/150 = 0.33 (does NOT meet 0.5 threshold)
	// With reciprocal mode, both must be satisfied
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{
		MinOverlap: 1,
		FractionA:  0.5,
		FractionB:  0.5,
		Reciprocal: true,
	})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 0 {
		t.Errorf("Expected 0 overlaps (reciprocal not met), got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectReciprocalMet(t *testing.T) {
	fileA := `chr1	100	200`

	fileB := `chr1	125	175`

	expected := `chr1	125	175
`

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	// Overlap: 50bp
	// Fraction of A: 50/100 = 0.5 (meets 0.5 threshold)
	// Fraction of B: 50/50 = 1.0 (meets 0.5 threshold)
	// Both satisfied
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{
		MinOverlap: 1,
		FractionA:  0.5,
		FractionB:  0.5,
		Reciprocal: true,
	})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 1 {
		t.Errorf("Expected 1 overlap (reciprocal met), got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectDistance(t *testing.T) {
	fileA := `chr1	100	200
chr1	300	400
chr1	500	600`

	fileB := `chr1	250	280`

	expected := `chr1	100	200	50
chr1	300	400	20
chr1	500	600	220
`

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{MinOverlap: 1, Distance: true})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 3 {
		t.Errorf("Expected 3 results, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectClosest(t *testing.T) {
	fileA := `chr1	100	200
chr1	500	600`

	fileB := `chr1	250	280
chr1	350	380`

	expected := `chr1	250	280
chr1	350	380
`

	readerA := strings.NewReader(fileA)
	readerB := strings.NewReader(fileB)
	var buf bytes.Buffer
	
	count, err := Intersect(readerA, readerB, &buf, IntersectOptions{MinOverlap: 1, Closest: true})
	if err != nil {
		t.Fatalf("Intersect failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 results, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestIntersectWithIntervalTree(t *testing.T) {
	fileA := `chr1	100	200
chr1	300	400
chr2	100	200`

	fileB := `chr1	150	250
chr1	350	450
chr2	150	250`

	// Test with linear search
	readerA1 := strings.NewReader(fileA)
	readerB1 := strings.NewReader(fileB)
	var buf1 bytes.Buffer
	
	count1, err := Intersect(readerA1, readerB1, &buf1, IntersectOptions{MinOverlap: 1, UseTree: false})
	if err != nil {
		t.Fatalf("Intersect (linear) failed: %v", err)
	}

	// Test with interval tree
	readerA2 := strings.NewReader(fileA)
	readerB2 := strings.NewReader(fileB)
	var buf2 bytes.Buffer
	
	count2, err := Intersect(readerA2, readerB2, &buf2, IntersectOptions{MinOverlap: 1, UseTree: true})
	if err != nil {
		t.Fatalf("Intersect (tree) failed: %v", err)
	}

	// Results should be identical
	if count1 != count2 {
		t.Errorf("Count mismatch: linear=%d, tree=%d", count1, count2)
	}

	if buf1.String() != buf2.String() {
		t.Errorf("Output mismatch.\nLinear:\n%s\nTree:\n%s", buf1.String(), buf2.String())
	}
}
