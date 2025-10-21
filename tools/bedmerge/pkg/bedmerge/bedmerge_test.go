package bedmerge

import (
	"bytes"
	"strings"
	"testing"
)

func TestMergeBasic(t *testing.T) {
	input := `chr1	100	200
chr1	150	250
chr1	300	400`

	expected := `chr1	100	250
chr1	300	400
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	count, err := Merge(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeAdjacent(t *testing.T) {
	input := `chr1	100	200
chr1	200	300
chr1	300	400`

	expected := `chr1	100	400
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	count, err := Merge(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 1 {
		t.Errorf("Expected 1 merged interval, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeNonOverlapping(t *testing.T) {
	input := `chr1	100	200
chr1	300	400
chr1	500	600`

	expected := `chr1	100	200
chr1	300	400
chr1	500	600
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	count, err := Merge(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 3 {
		t.Errorf("Expected 3 intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeMaxDistance(t *testing.T) {
	input := `chr1	100	200
chr1	250	300
chr1	500	600`

	expected := `chr1	100	300
chr1	500	600
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	// Merge intervals within 50bp
	count, err := Merge(reader, &buf, MergeOptions{MaxDistance: 50})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeMultipleChromosomes(t *testing.T) {
	input := `chr1	100	200
chr1	150	250
chr2	100	200
chr2	150	250
chr3	100	200`

	expected := `chr1	100	250
chr2	100	250
chr3	100	200
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	count, err := Merge(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 3 {
		t.Errorf("Expected 3 intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeEmpty(t *testing.T) {
	input := ``

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	count, err := Merge(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 0 {
		t.Errorf("Expected 0 intervals, got %d", count)
	}
}

func TestMergeWithStats(t *testing.T) {
	input := `chr1	100	200
chr1	150	250
chr1	300	400
chr1	350	450`

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	stats, err := MergeWithStats(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("MergeWithStats failed: %v", err)
	}

	if stats.InputIntervals != 4 {
		t.Errorf("Expected 4 input intervals, got %d", stats.InputIntervals)
	}

	if stats.OutputIntervals != 2 {
		t.Errorf("Expected 2 output intervals, got %d", stats.OutputIntervals)
	}

	if stats.MergedCount != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", stats.MergedCount)
	}
}

func TestMergeUnsorted(t *testing.T) {
	// Input is not sorted - should still work
	input := `chr1	300	400
chr1	100	200
chr1	150	250`

	expected := `chr1	100	250
chr1	300	400
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	count, err := Merge(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeWithCount(t *testing.T) {
	input := `chr1	100	200
chr1	150	250
chr1	300	400`

	expected := `chr1	100	250	2
chr1	300	400	1
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	opts := MergeOptions{
		OutputFields: OutputFields{
			Count: true,
		},
	}
	
	count, err := Merge(reader, &buf, opts)
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeBedGraph(t *testing.T) {
	input := `chr1	100	200	10
chr1	150	250	20
chr1	300	400	30`

	expected := `chr1	100	250	10
chr1	300	400	30
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	opts := MergeOptions{
		OutputFields: OutputFields{
			BedGraph: true,
		},
	}
	
	count, err := Merge(reader, &buf, opts)
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeStreaming(t *testing.T) {
	input := `chr1	100	200
chr1	150	250
chr2	100	200
chr2	150	250
chr3	100	200`

	expected := `chr1	100	250
chr2	100	250
chr3	100	200
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	opts := MergeOptions{
		Streaming: true,
	}
	
	count, err := Merge(reader, &buf, opts)
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 3 {
		t.Errorf("Expected 3 intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeStreamingWithCount(t *testing.T) {
	input := `chr1	100	200
chr1	150	250
chr2	100	200
chr2	150	250`

	expected := `chr1	100	250	2
chr2	100	250	2
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	opts := MergeOptions{
		Streaming: true,
		OutputFields: OutputFields{
			Count: true,
		},
	}
	
	count, err := Merge(reader, &buf, opts)
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeOutputFields(t *testing.T) {
	input := `chr1	100	200	name1	100	+
chr1	150	250	name2	200	+
chr1	300	400	name3	300	-`

	// Test with name and strand
	expected := `chr1	100	250	name1	100	+
chr1	300	400	name3	300	-
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer
	
	opts := MergeOptions{
		OutputFields: OutputFields{
			Name:   true,
			Score:  true,
			Strand: true,
		},
	}
	
	count, err := Merge(reader, &buf, opts)
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}
