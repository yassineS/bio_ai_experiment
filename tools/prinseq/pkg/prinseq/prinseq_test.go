package prinseq

import (
	"bytes"
	"strings"
	"testing"
)

func TestCalculateFastaStats(t *testing.T) {
	input := `>seq1
ACGTACGTACGT
>seq2
GGCCGGCCGGCCGGCC
>seq3
ATATATAT
`

	stats, err := CalculateStats(strings.NewReader(input), false)
	if err != nil {
		t.Fatalf("CalculateStats failed: %v", err)
	}

	if stats.NumReads != 3 {
		t.Errorf("Expected 3 reads, got %d", stats.NumReads)
	}

	expectedBases := 12 + 16 + 8
	if stats.TotalBases != expectedBases {
		t.Errorf("Expected %d total bases, got %d", expectedBases, stats.TotalBases)
	}

	if stats.MinLength != 8 {
		t.Errorf("Expected min length 8, got %d", stats.MinLength)
	}

	if stats.MaxLength != 16 {
		t.Errorf("Expected max length 16, got %d", stats.MaxLength)
	}

	expectedAvg := float64(expectedBases) / 3.0
	if stats.AvgLength != expectedAvg {
		t.Errorf("Expected avg length %.2f, got %.2f", expectedAvg, stats.AvgLength)
	}

	// seq1: 6 GC, seq2: 16 GC, seq3: 0 GC = 22 GC out of 36 bases = 61.11%
	expectedGC := 22.0 / 36.0 * 100.0
	if diff := stats.GCContent - expectedGC; diff > 0.01 || diff < -0.01 {
		t.Errorf("Expected GC content %.2f%%, got %.2f%%", expectedGC, stats.GCContent)
	}
}

func TestCalculateFastqStats(t *testing.T) {
	input := `@seq1
ACGTACGTACGT
+
IIIIIIIIIIII
@seq2
GGCCGGCCGGCCGGCC
+
HHHHHHHHHHHHHHHH
@seq3
ATATATAT
+
FFFFFFFF
`

	stats, err := CalculateStats(strings.NewReader(input), true)
	if err != nil {
		t.Fatalf("CalculateStats failed: %v", err)
	}

	if stats.NumReads != 3 {
		t.Errorf("Expected 3 reads, got %d", stats.NumReads)
	}

	expectedBases := 12 + 16 + 8
	if stats.TotalBases != expectedBases {
		t.Errorf("Expected %d total bases, got %d", expectedBases, stats.TotalBases)
	}

	if stats.AvgQuality == 0 {
		t.Error("Expected non-zero average quality for FASTQ")
	}
}

func TestCalculateStatsWithNs(t *testing.T) {
	input := `>seq1
ACGTNNNNACGT
>seq2
GGCCNNGGCCNN
`

	stats, err := CalculateStats(strings.NewReader(input), false)
	if err != nil {
		t.Fatalf("CalculateStats failed: %v", err)
	}

	expectedNs := 4 + 4
	if stats.NumNs != expectedNs {
		t.Errorf("Expected %d Ns, got %d", expectedNs, stats.NumNs)
	}
}

func TestFilterFastaByLength(t *testing.T) {
	input := `>seq1
ACGT
>seq2
ACGTACGTACGT
>seq3
ACGTACGT
`

	var output bytes.Buffer
	opts := FilterOptions{
		MinLen: 6,
		MaxLen: 10,
	}

	err := Filter(strings.NewReader(input), &output, false, opts)
	if err != nil {
		t.Fatalf("Filter failed: %v", err)
	}

	result := output.String()
	// Should keep only seq3 (length 8)
	if !strings.Contains(result, ">seq3") {
		t.Error("Expected seq3 to pass filter")
	}
	if strings.Contains(result, ">seq1") {
		t.Error("Expected seq1 to be filtered out (too short)")
	}
	if strings.Contains(result, ">seq2") {
		t.Error("Expected seq2 to be filtered out (too long)")
	}
}

func TestFilterFastaByGC(t *testing.T) {
	input := `>seq1_low_gc
AAAAAAAAAA
>seq2_medium_gc
ACGTACGTAC
>seq3_high_gc
GGGGGGGGGG
`

	var output bytes.Buffer
	opts := FilterOptions{
		MinGC: 40,
		MaxGC: 60,
	}

	err := Filter(strings.NewReader(input), &output, false, opts)
	if err != nil {
		t.Fatalf("Filter failed: %v", err)
	}

	result := output.String()
	// Should keep only seq2 (50% GC)
	if !strings.Contains(result, ">seq2_medium_gc") {
		t.Error("Expected seq2_medium_gc to pass filter")
	}
	if strings.Contains(result, ">seq1_low_gc") {
		t.Error("Expected seq1_low_gc to be filtered out (GC too low)")
	}
	if strings.Contains(result, ">seq3_high_gc") {
		t.Error("Expected seq3_high_gc to be filtered out (GC too high)")
	}
}

func TestFilterFastaByNs(t *testing.T) {
	input := `>seq1
ACGTACGTACGT
>seq2
ACGTNNNACGT
>seq3
ACGTNNNNNNNNACGT
`

	var output bytes.Buffer
	opts := FilterOptions{
		MaxNsN: 3,
	}

	err := Filter(strings.NewReader(input), &output, false, opts)
	if err != nil {
		t.Fatalf("Filter failed: %v", err)
	}

	result := output.String()
	// Should keep seq1 (0 Ns) and seq2 (3 Ns), filter seq3 (8 Ns)
	if !strings.Contains(result, ">seq1") {
		t.Error("Expected seq1 to pass filter")
	}
	if !strings.Contains(result, ">seq2") {
		t.Error("Expected seq2 to pass filter")
	}
	if strings.Contains(result, ">seq3") {
		t.Error("Expected seq3 to be filtered out (too many Ns)")
	}
}

func TestFilterFastqByQuality(t *testing.T) {
	input := `@seq1
ACGTACGTACGT
+
IIIIIIIIIIII
@seq2
ACGTACGTACGT
+
555555555555
@seq3
ACGTACGTACGT
+
############
`

	var output bytes.Buffer
	opts := FilterOptions{
		MinQualMean: 30, // Phred 30 = 'I' - 33 = 73 - 33 = 40 (wrong calculation, let me fix)
	}

	err := Filter(strings.NewReader(input), &output, true, opts)
	if err != nil {
		t.Fatalf("Filter failed: %v", err)
	}

	result := output.String()
	// seq1: 'I' = ASCII 73, qual = 40 (pass)
	// seq2: '5' = ASCII 53, qual = 20 (fail)
	// seq3: '#' = ASCII 35, qual = 2 (fail)
	if !strings.Contains(result, "@seq1") {
		t.Error("Expected seq1 to pass quality filter")
	}
	if strings.Contains(result, "@seq2") {
		t.Error("Expected seq2 to be filtered out (quality too low)")
	}
	if strings.Contains(result, "@seq3") {
		t.Error("Expected seq3 to be filtered out (quality too low)")
	}
}

func TestInvalidFastqFormat(t *testing.T) {
	tests := []struct {
		name  string
		input string
	}{
		{
			name: "missing plus line",
			input: `@seq1
ACGT
IIII
`,
		},
		{
			name: "mismatched quality length",
			input: `@seq1
ACGTACGT
+
III
`,
		},
		{
			name: "incomplete record",
			input: `@seq1
ACGT
+
`,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			_, err := CalculateStats(strings.NewReader(tt.input), true)
			if err == nil {
				t.Error("Expected error for invalid FASTQ format")
			}
		})
	}
}
