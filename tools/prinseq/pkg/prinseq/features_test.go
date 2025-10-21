package prinseq

import (
	"bytes"
	"encoding/json"
	"strings"
	"testing"
)

func TestEnhancedStats(t *testing.T) {
	input := `>seq1
ACGTACGTACGT
>seq2
GGCCGGCCGGCCGGCC
>seq3
ATATATAT
`

	stats, err := CalculateEnhancedStats(strings.NewReader(input), false)
	if err != nil {
		t.Fatalf("CalculateEnhancedStats failed: %v", err)
	}

	if stats.NumReads != 3 {
		t.Errorf("Expected 3 reads, got %d", stats.NumReads)
	}

	// Check length distribution
	if len(stats.LengthDistribution) == 0 {
		t.Error("Expected length distribution to be populated")
	}

	if stats.LengthDistribution[12] != 1 {
		t.Errorf("Expected 1 sequence of length 12, got %d", stats.LengthDistribution[12])
	}

	// Check base composition
	if len(stats.BaseComposition) == 0 {
		t.Error("Expected base composition to be populated")
	}

	// Check dinucleotides
	if len(stats.Dinucleotides) == 0 {
		t.Error("Expected dinucleotides to be populated")
	}
}

func TestEnhancedStatsFastq(t *testing.T) {
	input := `@seq1
ACGTACGTACGT
+
IIIIIIIIIIII
@seq2
GGCCGGCCGGCCGGCC
+
HHHHHHHHHHHHHHHH
`

	stats, err := CalculateEnhancedStats(strings.NewReader(input), true)
	if err != nil {
		t.Fatalf("CalculateEnhancedStats failed: %v", err)
	}

	if stats.NumReads != 2 {
		t.Errorf("Expected 2 reads, got %d", stats.NumReads)
	}

	// Check quality distribution
	if len(stats.QualityDistribution) == 0 {
		t.Error("Expected quality distribution to be populated")
	}

	// Check positional quality
	if len(stats.PositionalQuality) == 0 {
		t.Error("Expected positional quality to be populated")
	}
}

func TestJSONMarshal(t *testing.T) {
	input := `>seq1
ACGTACGT
`

	stats, err := CalculateEnhancedStats(strings.NewReader(input), false)
	if err != nil {
		t.Fatalf("CalculateEnhancedStats failed: %v", err)
	}

	// Test JSON marshaling
	jsonData, err := json.Marshal(stats)
	if err != nil {
		t.Fatalf("JSON marshal failed: %v", err)
	}

	// Test JSON unmarshaling
	var stats2 Stats
	if err := json.Unmarshal(jsonData, &stats2); err != nil {
		t.Fatalf("JSON unmarshal failed: %v", err)
	}

	if stats2.NumReads != stats.NumReads {
		t.Errorf("NumReads mismatch after JSON round-trip: %d vs %d", stats2.NumReads, stats.NumReads)
	}
}

func TestGenerateGraph(t *testing.T) {
	input := `>seq1
ACGTACGTACGT
>seq2
GGCCGGCCGGCCGGCC
>seq3
ATATATAT
`

	stats, err := CalculateEnhancedStats(strings.NewReader(input), false)
	if err != nil {
		t.Fatalf("CalculateEnhancedStats failed: %v", err)
	}

	// Test length graph
	var buf bytes.Buffer
	if err := GenerateGraph(stats, GraphTypeLength, &buf); err != nil {
		t.Errorf("GenerateGraph failed: %v", err)
	}

	output := buf.String()
	if !strings.Contains(output, "Length Distribution") {
		t.Error("Expected 'Length Distribution' in graph output")
	}

	// Test GC graph
	buf.Reset()
	if err := GenerateGraph(stats, GraphTypeGC, &buf); err != nil {
		t.Errorf("GenerateGraph GC failed: %v", err)
	}

	// Test dinucleotide graph
	buf.Reset()
	if err := GenerateGraph(stats, GraphTypeDinuc, &buf); err != nil {
		t.Errorf("GenerateGraph dinucleotides failed: %v", err)
	}
}

func TestGenerateSVG(t *testing.T) {
	input := `@seq1
ACGTACGTACGT
+
IIIIIIIIIIII
`

	stats, err := CalculateEnhancedStats(strings.NewReader(input), true)
	if err != nil {
		t.Fatalf("CalculateEnhancedStats failed: %v", err)
	}

	var buf bytes.Buffer
	if err := GenerateSVG(stats, &buf); err != nil {
		t.Errorf("GenerateSVG failed: %v", err)
	}

	output := buf.String()
	if !strings.Contains(output, "<svg") {
		t.Error("Expected SVG tag in output")
	}
	if !strings.Contains(output, "</svg>") {
		t.Error("Expected closing SVG tag in output")
	}
}

func TestGenerateHTMLReport(t *testing.T) {
	input := `@seq1
ACGTACGTACGT
+
IIIIIIIIIIII
@seq2
GGCCGGCCGGCCGGCC
+
HHHHHHHHHHHHHHHH
`

	stats, err := CalculateEnhancedStats(strings.NewReader(input), true)
	if err != nil {
		t.Fatalf("CalculateEnhancedStats failed: %v", err)
	}

	var buf bytes.Buffer
	if err := GenerateHTMLReport(stats, &buf); err != nil {
		t.Errorf("GenerateHTMLReport failed: %v", err)
	}

	output := buf.String()
	if !strings.Contains(output, "<!DOCTYPE html>") {
		t.Error("Expected HTML doctype in output")
	}
	if !strings.Contains(output, "PRINSEQ Quality Control Report") {
		t.Error("Expected report title in output")
	}
	if !strings.Contains(output, "<svg") {
		t.Error("Expected embedded SVG in report")
	}
}

func TestBenchmarkStats(t *testing.T) {
	input := `>seq1
ACGTACGTACGT
>seq2
GGCCGGCCGGCCGGCC
`

	result, stats, err := BenchmarkStats(strings.NewReader(input), false)
	if err != nil {
		t.Fatalf("BenchmarkStats failed: %v", err)
	}

	if result.Operation != "stats" {
		t.Errorf("Expected operation 'stats', got '%s'", result.Operation)
	}

	if result.Duration == 0 {
		t.Error("Expected non-zero duration")
	}

	if stats.NumReads != 2 {
		t.Errorf("Expected 2 reads, got %d", stats.NumReads)
	}

	if result.ThroughputMBs == 0 {
		t.Error("Expected non-zero throughput")
	}
}

func TestBenchmarkFilter(t *testing.T) {
	input := `>seq1
ACGTACGTACGT
>seq2
GGCCGGCCGGCCGGCC
>seq3
AT
`

	opts := FilterOptions{
		MinLen: 10,
	}

	result, err := BenchmarkFilter(strings.NewReader(input), false, opts)
	if err != nil {
		t.Fatalf("BenchmarkFilter failed: %v", err)
	}

	if result.Operation != "filter" {
		t.Errorf("Expected operation 'filter', got '%s'", result.Operation)
	}

	if result.Duration == 0 {
		t.Error("Expected non-zero duration")
	}
}

func TestRunBenchmarkSuite(t *testing.T) {
	input := `@seq1
ACGTACGTACGT
+
IIIIIIIIIIII
@seq2
GGCCGGCCGGCCGGCC
+
HHHHHHHHHHHHHHHH
`

	results, err := RunBenchmarkSuite(strings.NewReader(input), true)
	if err != nil {
		t.Fatalf("RunBenchmarkSuite failed: %v", err)
	}

	if len(results) == 0 {
		t.Error("Expected non-empty benchmark results")
	}

	// Check that we have multiple benchmark operations
	operations := make(map[string]bool)
	for _, r := range results {
		operations[r.Operation] = true
	}

	if !operations["stats"] {
		t.Error("Expected 'stats' operation in benchmark results")
	}
}

func TestFormatBenchmarkResults(t *testing.T) {
	results := []*BenchmarkResult{
		{
			Operation:     "test_op",
			DurationMs:    10.5,
			ThroughputMBs: 100.0,
			ReadsPerSec:   1000.0,
		},
	}

	output := FormatBenchmarkResults(results)
	if !strings.Contains(output, "Benchmark Results") {
		t.Error("Expected 'Benchmark Results' header")
	}
	if !strings.Contains(output, "test_op") {
		t.Error("Expected operation name in output")
	}
}
