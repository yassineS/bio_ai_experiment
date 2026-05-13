package skewer

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

func TestTrimSingleEndNoAdapter(t *testing.T) {
	input := `@read1
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
`

	var output bytes.Buffer
	opts := DefaultTrimOptions()
	opts.MinLength = 10 // Lower threshold to ensure read passes

	stats, err := TrimSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd failed: %v", err)
	}

	if stats.TotalReads != 1 {
		t.Errorf("Expected 1 total read, got %d", stats.TotalReads)
	}
	if stats.TrimmedReads != 0 {
		t.Errorf("Expected 0 trimmed reads, got %d", stats.TrimmedReads)
	}
	if stats.DiscardedReads != 0 {
		t.Errorf("Expected 0 discarded reads, got %d", stats.DiscardedReads)
	}
}

func TestTrimSingleEndWith3PrimeAdapter(t *testing.T) {
	// Read with adapter at 3' end
	input := `@read1
ACGTACGTACGTAGATCGGAAGAGC
+
IIIIIIIIIIIIIIIIIIIIIIIII
`

	var output bytes.Buffer
	opts := DefaultTrimOptions()
	opts.Adapter3 = "AGATCGGAAGAGC"

	stats, err := TrimSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd failed: %v", err)
	}

	if stats.TotalReads != 1 {
		t.Errorf("Expected 1 total read, got %d", stats.TotalReads)
	}
	if stats.AdapterFound3 != 1 {
		t.Errorf("Expected 1 adapter found, got %d", stats.AdapterFound3)
	}

	// Check that adapter was removed
	result := output.String()
	if strings.Contains(result, "AGATCGGAAGAGC") {
		t.Error("Adapter should have been removed from output")
	}
}

func TestTrimSingleEndWith5PrimeAdapter(t *testing.T) {
	// Read with adapter at 5' end
	input := `@read1
AGATCGGAAGAGCACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIII
`

	var output bytes.Buffer
	opts := DefaultTrimOptions()
	opts.Adapter5 = "AGATCGGAAGAGC"

	stats, err := TrimSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd failed: %v", err)
	}

	if stats.TotalReads != 1 {
		t.Errorf("Expected 1 total read, got %d", stats.TotalReads)
	}
	if stats.AdapterFound5 != 1 {
		t.Errorf("Expected 1 adapter found, got %d", stats.AdapterFound5)
	}

	// Check that adapter was removed
	result := output.String()
	if strings.Contains(result, "AGATCGGAAGAGC") {
		t.Error("Adapter should have been removed from output")
	}
}

func TestTrimSingleEndShortRead(t *testing.T) {
	// Read too short after trimming
	input := `@read1
ACGTACGTAGATCGGAAGAGC
+
IIIIIIIIIIIIIIIIIIIII
`

	var output bytes.Buffer
	opts := DefaultTrimOptions()
	opts.Adapter3 = "AGATCGGAAGAGC"
	opts.MinLength = 18

	stats, err := TrimSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd failed: %v", err)
	}

	if stats.TotalReads != 1 {
		t.Errorf("Expected 1 total read, got %d", stats.TotalReads)
	}
	if stats.DiscardedReads != 1 {
		t.Errorf("Expected 1 discarded read, got %d", stats.DiscardedReads)
	}
}

func TestFindAdapter(t *testing.T) {
	tests := []struct {
		name       string
		seq        string
		adapter    string
		minOverlap int
		errorRate  float64
		want       int
	}{
		{
			name:       "exact match",
			seq:        "ACGTACGTAGATCGGAAGAGC",
			adapter:    "AGATCGGAAGAGC",
			minOverlap: 3,
			errorRate:  0.1,
			want:       8, // Correct position (0-indexed)
		},
		{
			name:       "no match",
			seq:        "ACGTACGTACGTACGT",
			adapter:    "AGATCGGAAGAGC",
			minOverlap: 3,
			errorRate:  0.1,
			want:       -1,
		},
		{
			name:       "partial match at end",
			seq:        "ACGTACGTACGTAGA",
			adapter:    "AGATCGGAAGAGC",
			minOverlap: 3,
			errorRate:  0.1,
			want:       12,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := findAdapter(tt.seq, tt.adapter, tt.minOverlap, tt.errorRate)
			if got != tt.want {
				t.Errorf("findAdapter() = %d, want %d", got, tt.want)
			}
		})
	}
}

func TestTrimPairedEnd(t *testing.T) {
	input1 := `@read1
ACGTACGTACGTACGTAGATCGGAAGAGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIII
`
	input2 := `@read1
TGCATGCATGCATGCAAGATCGGAAGAGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIII
`

	var output1, output2 bytes.Buffer
	opts := DefaultTrimOptions()
	opts.Adapter3 = "AGATCGGAAGAGC"

	stats, err := TrimPairedEnd(
		strings.NewReader(input1),
		strings.NewReader(input2),
		&output1,
		&output2,
		nil,
		fastq.Phred33,
		opts,
	)

	if err != nil {
		t.Fatalf("TrimPairedEnd failed: %v", err)
	}

	if stats.TotalReads != 2 {
		t.Errorf("Expected 2 total reads, got %d", stats.TotalReads)
	}
	if stats.AdapterFound3 != 2 {
		t.Errorf("Expected 2 adapters found, got %d", stats.AdapterFound3)
	}
}

func TestTrimRecordBothEnds(t *testing.T) {
	record := &fastq.Record{
		ID:          "read1",
		Description: "read1",
		Sequence:    []byte("AGATCGGAAGAGCACGTACGTACGTAGATCGGAAGAGC"),
		Quality:     []byte("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"),
	}

	opts := DefaultTrimOptions()
	opts.Adapter5 = "AGATCGGAAGAGC"
	opts.Adapter3 = "AGATCGGAAGAGC"

	stats := &TrimStats{}
	trimmed := trimRecord(record, opts, stats)

	// Should have both adapters removed
	seq := string(trimmed.Sequence)
	if strings.Contains(seq, "AGATCGGAAGAGC") {
		t.Error("Adapters should have been removed")
	}

	if stats.AdapterFound5 != 1 {
		t.Errorf("Expected 1 5' adapter found, got %d", stats.AdapterFound5)
	}
	if stats.AdapterFound3 != 1 {
		t.Errorf("Expected 1 3' adapter found, got %d", stats.AdapterFound3)
	}
}

func TestJSONOutput(t *testing.T) {
	stats := &TrimStats{
		TotalReads:     100,
		TrimmedReads:   75,
		AdapterFound3:  75,
		AdapterFound5:  0,
		DiscardedReads: 5,
		TotalBases:     10000,
		TrimmedBases:   2500,
	}

	json := stats.ToJSON()

	if !strings.Contains(json, "total_reads") {
		t.Error("JSON should contain total_reads field")
	}
	if !strings.Contains(json, "100") {
		t.Error("JSON should contain the value 100")
	}
}

func TestHTMLOutput(t *testing.T) {
	stats := &TrimStats{
		TotalReads:     100,
		TrimmedReads:   75,
		AdapterFound3:  75,
		AdapterFound5:  0,
		DiscardedReads: 5,
		TotalBases:     10000,
		TrimmedBases:   2500,
	}

	html := stats.ToHTML()

	if !strings.Contains(html, "<!DOCTYPE html>") {
		t.Error("HTML should start with DOCTYPE")
	}
	if !strings.Contains(html, "Adapter Trimming Report") {
		t.Error("HTML should contain title")
	}
	if !strings.Contains(html, "100") {
		t.Error("HTML should contain total reads value")
	}
}

func TestUMIExtraction5Prime(t *testing.T) {
	record := &fastq.Record{
		ID:          "read1",
		Description: "read1",
		Sequence:    []byte("ACGTACGTACGTACGTACGT"),
		Quality:     []byte("IIIIIIIIIIIIIIIIIIII"),
	}

	opts := DefaultTrimOptions()
	opts.UMILength = 8
	opts.UMIPosition = "5prime"

	trimmed, umi := extractUMI(record, opts)

	if umi != "ACGTACGT" {
		t.Errorf("Expected UMI 'ACGTACGT', got '%s'", umi)
	}
	if string(trimmed.Sequence) != "ACGTACGTACGT" {
		t.Errorf("Expected sequence 'ACGTACGTACGT', got '%s'", string(trimmed.Sequence))
	}
}

func TestUMIExtraction3Prime(t *testing.T) {
	record := &fastq.Record{
		ID:          "read1",
		Description: "read1",
		Sequence:    []byte("ACGTACGTACGTACGTACGT"),
		Quality:     []byte("IIIIIIIIIIIIIIIIIIII"),
	}

	opts := DefaultTrimOptions()
	opts.UMILength = 8
	opts.UMIPosition = "3prime"

	trimmed, umi := extractUMI(record, opts)

	if umi != "ACGTACGT" {
		t.Errorf("Expected UMI 'ACGTACGT', got '%s'", umi)
	}
	if string(trimmed.Sequence) != "ACGTACGTACGT" {
		t.Errorf("Expected sequence 'ACGTACGTACGT', got '%s'", string(trimmed.Sequence))
	}
}

func TestImprovedFindAdapter(t *testing.T) {
	tests := []struct {
		name       string
		seq        string
		adapter    string
		minOverlap int
		errorRate  float64
		wantPos    int
	}{
		{
			name:       "exact match",
			seq:        "ACGTACGTAGATCGGAAGAGC",
			adapter:    "AGATCGGAAGAGC",
			minOverlap: 3,
			errorRate:  0.1,
			wantPos:    8,
		},
		{
			name:       "match with error",
			seq:        "ACGTACGTAGATCAGAAGAGC",
			adapter:    "AGATCGGAAGAGC",
			minOverlap: 3,
			errorRate:  0.15,
			wantPos:    8,
		},
		{
			name:       "no match",
			seq:        "ACGTACGTACGTACGT",
			adapter:    "AGATCGGAAGAGC",
			minOverlap: 3,
			errorRate:  0.1,
			wantPos:    -1,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := improvedFindAdapter(tt.seq, tt.adapter, tt.minOverlap, tt.errorRate)
			if got != tt.wantPos {
				t.Errorf("improvedFindAdapter() = %d, want %d", got, tt.wantPos)
			}
		})
	}
}

func TestUMIStats(t *testing.T) {
	input := `@read1
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@read2
ACGTACGTTTTTTTTTTTTT
+
IIIIIIIIIIIIIIIIIIII
@read3
GGGGGGGGACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
`

	var output bytes.Buffer
	opts := DefaultTrimOptions()
	opts.UMILength = 8
	opts.UMIPosition = "5prime"

	stats, err := TrimSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd failed: %v", err)
	}

	if stats.UMIStats == nil {
		t.Fatal("UMI stats should not be nil")
	}

	if stats.UMIStats.TotalUMIs != 3 {
		t.Errorf("Expected 3 total UMIs, got %d", stats.UMIStats.TotalUMIs)
	}

	if stats.UMIStats.UniqueUMIs != 2 {
		t.Errorf("Expected 2 unique UMIs, got %d", stats.UMIStats.UniqueUMIs)
	}

	// Check that ACGTACGT appears twice
	if count, ok := stats.UMIStats.UMIDistribution["ACGTACGT"]; !ok || count != 2 {
		t.Errorf("Expected UMI 'ACGTACGT' to appear 2 times, got %d", count)
	}
}
