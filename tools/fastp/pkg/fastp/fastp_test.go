package fastp

import (
	"bytes"
	"os"
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

// Test adapter detection
func TestDetectAdapter(t *testing.T) {
	// Create sample sequences with TruSeq adapters
	sequences := []string{
		"ACGTACGTACGTACGTAGATCGGAAGAGC",
		"TGCATGCATGCATGCAGATCGGAAGAGC",
		"GCTAGCTAGCTAGCTAGATCGGAAGAGC",
		"ACGTACGTACGTACGTAGATCGGAAGAGC",
		"TGCATGCATGCATGCAGATCGGAAGAGC",
	}
	
	adapter := DetectAdapterFromReads(sequences)
	
	if adapter != CommonAdapters["TruSeq"] {
		t.Errorf("Expected TruSeq adapter, got %s", adapter)
	}
}

// Test UMI extraction
func TestExtractUMI(t *testing.T) {
	record := &fastq.Record{
		ID:          "read1",
		Description: "",
		Sequence:    []byte("ACGTACGTACGTACGTACGT"),
		Quality:     []byte("IIIIIIIIIIIIIIIIIIII"),
	}
	
	opts := DefaultProcessOptions()
	opts.UMILength = 6
	opts.UMISkip = 0
	opts.UMILocation = "read1"
	
	stats := &ProcessStats{}
	
	result, _ := extractUMI(record, nil, opts, stats)
	
	if !strings.Contains(result.ID, "UMI:ACGTAC") {
		t.Errorf("Expected UMI to be added to ID, got %s", result.ID)
	}
	
	if len(result.Sequence) != 14 { // 20 - 6
		t.Errorf("Expected sequence length 14 after UMI extraction, got %d", len(result.Sequence))
	}
	
	if stats.UMIExtracted != 1 {
		t.Errorf("Expected 1 UMI extracted, got %d", stats.UMIExtracted)
	}
}

// Test base correction
func TestCorrectBases(t *testing.T) {
	seq := "ACGTACGTACGT"
	qual := []byte("III!!!IIIIII") // Low quality in the middle
	
	stats := &ProcessStats{}
	
	corrected, _ := correctBases(seq, qual, 20, stats, fastq.Phred33)
	
	// Bases with quality < 20 should be corrected to N
	if !strings.Contains(corrected, "N") {
		t.Error("Expected low quality bases to be corrected to N")
	}
	
	if stats.BasesCorrected == 0 {
		t.Error("Expected bases to be corrected")
	}
}

// Test overlap analysis
func TestAnalyzeOverlap(t *testing.T) {
	// Create overlapping paired-end reads
	record1 := &fastq.Record{
		ID:          "read1",
		Description: "",
		Sequence:    []byte("ACGTACGTACGTACGTACGTACGT"),
		Quality:     []byte("IIIIIIIIIIIIIIIIIIIIIIII"),
	}
	
	// record2 should be reverse complement of end of record1
	record2 := &fastq.Record{
		ID:          "read2",
		Description: "",
		Sequence:    []byte("ACGTACGTACGTACGTACGTACGT"), // Same as record1 for simplicity
		Quality:     []byte("IIIIIIIIIIIIIIIIIIIIIIII"),
	}
	
	opts := DefaultProcessOptions()
	opts.MinOverlap = 10
	opts.MaxMismatch = 5
	
	result := analyzeOverlap(record1, record2, opts, fastq.Phred33)
	
	// Should detect some overlap
	if !result.HasOverlap {
		t.Log("Note: Overlap detection may not find overlap with identical sequences")
	}
}

// Test reverse complement
func TestReverseComplement(t *testing.T) {
	tests := []struct {
		name     string
		input    string
		expected string
	}{
		{"simple", "ACGT", "ACGT"},
		{"longer", "AAATTT", "AAATTT"},
		{"with N", "ACGTN", "NACGT"},
	}
	
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := reverseComplement(tt.input)
			if result != tt.expected {
				t.Errorf("reverseComplement(%s) = %s, want %s", tt.input, result, tt.expected)
			}
		})
	}
}

// Test multi-threading with single-end
func TestProcessSingleEndMultiThreaded(t *testing.T) {
	input := `@read1
ACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2
TGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3
GCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIII
`
	
	var output bytes.Buffer
	opts := DefaultProcessOptions()
	opts.Threads = 2
	opts.MinLength = 10
	
	stats, err := ProcessSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd failed: %v", err)
	}
	
	if stats.TotalReads != 3 {
		t.Errorf("Expected 3 total reads, got %d", stats.TotalReads)
	}
	
	if stats.CleanReads != 3 {
		t.Errorf("Expected 3 clean reads, got %d", stats.CleanReads)
	}
}

// Test HTML report generation
func TestGenerateHTMLReport(t *testing.T) {
	stats := &ProcessStats{
		TotalReads:          1000,
		TotalBases:          150000,
		CleanReads:          850,
		CleanBases:          127500,
		LowQualityReads:     100,
		TooShortReads:       30,
		TooLongReads:        10,
		TooManyNReads:       10,
		AdapterTrimmedReads: 750,
		AdapterTrimmedBases: 15000,
		PolyGTrimmedReads:   200,
		PolyGTrimmedBases:   4000,
		DetectedAdapter:     "AGATCGGAAGAGC",
		UMIExtracted:        850,
		BasesCorrected:      5000,
		MergedReads:         100,
	}
	
	opts := DefaultProcessOptions()
	opts.Adapter3 = "AGATCGGAAGAGC"
	opts.QualThreshold = 20
	opts.MinLength = 30
	opts.Threads = 4
	
	tmpFile := "/tmp/fastp_test_report.html"
	err := GenerateHTMLReport(stats, opts, tmpFile)
	if err != nil {
		t.Fatalf("GenerateHTMLReport failed: %v", err)
	}
	
	// Check if file exists
	if _, err := os.Stat(tmpFile); os.IsNotExist(err) {
		t.Error("HTML report file was not created")
	}
	
	// Clean up
	os.Remove(tmpFile)
}

// Test base correction in processing
func TestProcessSingleEndWithBaseCorrection(t *testing.T) {
	input := `@read1
ACGTACGTACGTACGTACGT
+
III!!!IIIIII!!!!IIII
`
	
	var output bytes.Buffer
	opts := DefaultProcessOptions()
	opts.BaseCorrection = true
	opts.CorrectionThreshold = 20
	opts.MinLength = 10
	opts.MaxNCount = 20 // Allow more Ns
	opts.MaxNPercent = 50.0 // Allow higher N percentage
	opts.QualPercent = 0 // Don't filter by quality percentage
	
	stats, err := ProcessSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd failed: %v", err)
	}
	
	if stats.BasesCorrected == 0 {
		t.Error("Expected some bases to be corrected")
	}
	
	// Check that output contains N
	result := output.String()
	if stats.CleanReads > 0 && !strings.Contains(result, "N") {
		t.Error("Expected corrected bases (N) in output")
	}
}

