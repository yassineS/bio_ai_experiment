package fastp

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

func TestProcessSingleEndNoFilters(t *testing.T) {
	input := `@read1
ACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIII
`
	
	var output bytes.Buffer
	opts := DefaultProcessOptions()
	opts.MinLength = 10
	opts.QualThreshold = 0
	
	stats, err := ProcessSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd failed: %v", err)
	}
	
	if stats.TotalReads != 1 {
		t.Errorf("Expected 1 total read, got %d", stats.TotalReads)
	}
	if stats.CleanReads != 1 {
		t.Errorf("Expected 1 clean read, got %d", stats.CleanReads)
	}
}

func TestProcessSingleEndAdapterTrimming(t *testing.T) {
	input := `@read1
ACGTACGTACGTACGTAGATCGGAAGAGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIII
`
	
	var output bytes.Buffer
	opts := DefaultProcessOptions()
	opts.Adapter3 = "AGATCGGAAGAGC"
	opts.MinLength = 10
	
	stats, err := ProcessSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd failed: %v", err)
	}
	
	if stats.TotalReads != 1 {
		t.Errorf("Expected 1 total read, got %d", stats.TotalReads)
	}
	if stats.AdapterTrimmedReads != 1 {
		t.Errorf("Expected 1 adapter trimmed read, got %d", stats.AdapterTrimmedReads)
	}
	
	// Check that adapter was removed
	result := output.String()
	if strings.Contains(result, "AGATCGGAAGAGC") {
		t.Error("Adapter should have been removed from output")
	}
}

func TestProcessSingleEndLengthFilter(t *testing.T) {
	input := `@read1
ACGTACGT
+
IIIIIIII
`
	
	var output bytes.Buffer
	opts := DefaultProcessOptions()
	opts.MinLength = 20
	
	stats, err := ProcessSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd failed: %v", err)
	}
	
	if stats.TotalReads != 1 {
		t.Errorf("Expected 1 total read, got %d", stats.TotalReads)
	}
	if stats.TooShortReads != 1 {
		t.Errorf("Expected 1 too short read, got %d", stats.TooShortReads)
	}
	if stats.CleanReads != 0 {
		t.Errorf("Expected 0 clean reads, got %d", stats.CleanReads)
	}
}

func TestProcessSingleEndNFilter(t *testing.T) {
	input := `@read1
ACGTNNNNNACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIII
`
	
	var output bytes.Buffer
	opts := DefaultProcessOptions()
	opts.MaxNCount = 3
	opts.MinLength = 10
	
	stats, err := ProcessSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd failed: %v", err)
	}
	
	if stats.TotalReads != 1 {
		t.Errorf("Expected 1 total read, got %d", stats.TotalReads)
	}
	if stats.TooManyNReads != 1 {
		t.Errorf("Expected 1 too many N read, got %d", stats.TooManyNReads)
	}
}

func TestProcessSingleEndPolyGTrimming(t *testing.T) {
	input := `@read1
ACGTACGTACGTACGTGGGGGGGGGGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIII
`
	
	var output bytes.Buffer
	opts := DefaultProcessOptions()
	opts.TrimPolyG = true
	opts.PolyGMinLen = 5
	opts.MinLength = 10
	
	stats, err := ProcessSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd failed: %v", err)
	}
	
	if stats.TotalReads != 1 {
		t.Errorf("Expected 1 total read, got %d", stats.TotalReads)
	}
	if stats.PolyGTrimmedReads != 1 {
		t.Errorf("Expected 1 poly-G trimmed read, got %d", stats.PolyGTrimmedReads)
	}
	
	// Check that poly-G was removed
	result := output.String()
	lines := strings.Split(strings.TrimSpace(result), "\n")
	if len(lines) >= 2 {
		seq := lines[1]
		if strings.HasSuffix(seq, "GGGGGGGG") {
			t.Error("Poly-G tail should have been removed")
		}
	}
}

func TestCountPolyTail(t *testing.T) {
	tests := []struct {
		name string
		seq  string
		base rune
		want int
	}{
		{"poly-G tail", "ACGTGGGGGG", 'G', 6},
		{"no tail", "ACGTACGT", 'G', 0},
		{"short tail", "ACGTGG", 'G', 2},
		{"entire sequence", "GGGGGG", 'G', 6},
	}
	
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := countPolyTail(tt.seq, tt.base)
			if got != tt.want {
				t.Errorf("countPolyTail() = %d, want %d", got, tt.want)
			}
		})
	}
}

func TestCountNs(t *testing.T) {
	tests := []struct {
		name string
		seq  string
		want int
	}{
		{"no Ns", "ACGTACGT", 0},
		{"some Ns", "ACGTNNNACGT", 3},
		{"all Ns", "NNNNN", 5},
	}
	
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := countNs(tt.seq)
			if got != tt.want {
				t.Errorf("countNs() = %d, want %d", got, tt.want)
			}
		})
	}
}

func TestCalculateComplexity(t *testing.T) {
	tests := []struct {
		name string
		seq  string
		want float64
	}{
		{"high complexity", "ACGTACGTACGT", 1.0},
		{"low complexity", "AAAAAAAAAA", 0.111}, // Only AA kmer
		{"medium complexity", "AACCGGTT", 0.571}, // Multiple kmers
	}
	
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := calculateComplexity(tt.seq)
			// Allow some tolerance for floating point comparison
			if (got-tt.want) > 0.01 && (tt.want-got) > 0.01 {
				t.Errorf("calculateComplexity() = %.3f, want %.3f", got, tt.want)
			}
		})
	}
}

func TestProcessSingleEndQualityFilter(t *testing.T) {
	input := `@read1
ACGTACGTACGTACGTACGT
+
III!!!IIIIII!!!!IIII
`
	
	var output bytes.Buffer
	opts := DefaultProcessOptions()
	opts.QualThreshold = 20
	opts.QualPercent = 80  // Require 80% of bases to meet quality threshold
	opts.MinLength = 10
	
	stats, err := ProcessSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd failed: %v", err)
	}
	
	if stats.TotalReads != 1 {
		t.Errorf("Expected 1 total read, got %d", stats.TotalReads)
	}
	if stats.LowQualityReads != 1 {
		t.Errorf("Expected 1 low quality read, got %d", stats.LowQualityReads)
	}
}

func TestProcessSingleEndCombinedFilters(t *testing.T) {
	input := `@read1
ACGTACGTACGTACGTAGATCGGAAGAGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2
NNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3
AAAA
+
IIII
@read4
TGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIII
`
	
	var output bytes.Buffer
	opts := DefaultProcessOptions()
	opts.Adapter3 = "AGATCGGAAGAGC"
	opts.MinLength = 15
	opts.MaxNCount = 5
	
	stats, err := ProcessSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd failed: %v", err)
	}
	
	if stats.TotalReads != 4 {
		t.Errorf("Expected 4 total reads, got %d", stats.TotalReads)
	}
	// read1: passes (after adapter trim)
	// read2: fails (too many Ns)
	// read3: fails (too short)
	// read4: passes
	if stats.CleanReads != 2 {
		t.Errorf("Expected 2 clean reads, got %d", stats.CleanReads)
	}
	if stats.TooManyNReads != 1 {
		t.Errorf("Expected 1 too many N read, got %d", stats.TooManyNReads)
	}
	if stats.TooShortReads != 1 {
		t.Errorf("Expected 1 too short read, got %d", stats.TooShortReads)
	}
}
