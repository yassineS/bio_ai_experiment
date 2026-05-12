package sickle

import (
	"bytes"
	"fmt"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

func TestTrimSingleEnd(t *testing.T) {
	// Test data with varying quality scores
	input := `@read1
ACGTACGTACGTACGT
+
IIIIIIII########
@read2
NNNNACGTACGT
+
############
@read3
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
`

	var output bytes.Buffer
	opts := DefaultTrimOptions()
	opts.QualThreshold = 30 // '#' is quality 2, 'I' is quality 40 in Phred+33

	stats, err := TrimSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd failed: %v", err)
	}

	if stats.TotalReads != 3 {
		t.Errorf("Expected 3 total reads, got %d", stats.TotalReads)
	}

	// read1 should be trimmed (low quality at end)
	// read2 should be discarded or heavily trimmed
	// read3 should pass unchanged

	result := output.String()
	if !strings.Contains(result, "read3") {
		t.Error("Expected read3 in output")
	}
}

func TestTrimSingleEndWithTruncateN(t *testing.T) {
	input := `@read1
ACGTNNNNACGT
+
IIIIIIIIIIII
@read2
ACGTACGTACGT
+
IIIIIIIIIIII
`

	var output bytes.Buffer
	opts := DefaultTrimOptions()
	opts.TruncateN = true
	opts.LengthThreshold = 4

	stats, err := TrimSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd with TruncateN failed: %v", err)
	}

	result := output.String()

	// read1 should be truncated at first N (position 4)
	if strings.Contains(result, "ACGTNNNN") {
		t.Error("Expected read1 to be truncated before N's")
	}

	// read2 should pass through
	if !strings.Contains(result, "ACGTACGTACGT") {
		t.Error("Expected read2 to pass unchanged")
	}

	if stats.TotalReads != 2 {
		t.Errorf("Expected 2 total reads, got %d", stats.TotalReads)
	}
}

func TestTrimPairedEnd(t *testing.T) {
	input1 := `@read1/1
ACGTACGTACGTACGT
+
IIIIIIII########
@read2/1
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
`

	input2 := `@read1/2
TGCATGCATGCATGCA
+
IIIIIIIIIIIIIIII
@read2/2
TGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIII
`

	var output1, output2, outputSingle bytes.Buffer
	opts := DefaultTrimOptions()
	opts.QualThreshold = 30

	stats, err := TrimPairedEnd(
		strings.NewReader(input1),
		strings.NewReader(input2),
		&output1,
		&output2,
		&outputSingle,
		fastq.Phred33,
		opts,
	)

	if err != nil {
		t.Fatalf("TrimPairedEnd failed: %v", err)
	}

	if stats.TotalReads != 4 {
		t.Errorf("Expected 4 total reads (2 pairs), got %d", stats.TotalReads)
	}

	// read2 pair should pass through
	result1 := output1.String()
	result2 := output2.String()

	if !strings.Contains(result1, "read2/1") || !strings.Contains(result2, "read2/2") {
		t.Error("Expected read2 pair in paired output")
	}
}

func TestTrimPairedEndWithSingleOutput(t *testing.T) {
	// First read pair: both pass
	// Second read pair: only forward passes
	input1 := `@read1/1
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@read2/1
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
`

	input2 := `@read1/2
TGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIII
@read2/2
NNNN
+
####
`

	var output1, output2, outputSingle bytes.Buffer
	opts := DefaultTrimOptions()
	opts.TruncateN = true
	opts.LengthThreshold = 10

	stats, err := TrimPairedEnd(
		strings.NewReader(input1),
		strings.NewReader(input2),
		&output1,
		&output2,
		&outputSingle,
		fastq.Phred33,
		opts,
	)

	if err != nil {
		t.Fatalf("TrimPairedEnd with single output failed: %v", err)
	}

	if stats.TotalReads != 4 {
		t.Errorf("Expected 4 total reads, got %d", stats.TotalReads)
	}

	// read1 pair should be in paired output
	result1 := output1.String()
	result2 := output2.String()
	if !strings.Contains(result1, "read1/1") || !strings.Contains(result2, "read1/2") {
		t.Error("Expected read1 pair in paired output")
	}

	// read2/1 should be in single output (read2/2 fails)
	resultSingle := outputSingle.String()
	if !strings.Contains(resultSingle, "read2/1") {
		t.Error("Expected read2/1 in single output")
	}
}

func TestTrimRecordBasic(t *testing.T) {
	record := &fastq.Record{
		ID:       "test",
		Sequence: []byte("ACGTACGTACGT"),
		Quality:  []byte("IIIIII######"),
	}

	opts := DefaultTrimOptions()
	opts.QualThreshold = 30

	trimmed := trimRecord(record, opts)

	// Should trim low quality end
	if len(trimmed.Sequence) >= len(record.Sequence) {
		t.Error("Expected sequence to be trimmed")
	}

	if len(trimmed.Sequence) != len(trimmed.Quality) {
		t.Error("Sequence and quality lengths don't match after trimming")
	}
}

func TestTrimRecordWithNoFivePrime(t *testing.T) {
	record := &fastq.Record{
		ID:       "test",
		Sequence: []byte("ACGTACGTACGT"),
		Quality:  []byte("####IIIIIIII"),
	}

	opts := DefaultTrimOptions()
	opts.QualThreshold = 30
	opts.NoFivePrime = true

	trimmed := trimRecord(record, opts)

	// Should not trim 5' end even though quality is low
	if len(trimmed.Sequence) == 0 {
		t.Error("Sequence was completely trimmed when NoFivePrime was set")
	}
}

func TestTrimRecordWithTruncateN(t *testing.T) {
	record := &fastq.Record{
		ID:       "test",
		Sequence: []byte("ACGTNNNNACGT"),
		Quality:  []byte("IIIIIIIIIIII"),
	}

	opts := DefaultTrimOptions()
	opts.TruncateN = true

	trimmed := trimRecord(record, opts)

	// Should truncate at first N
	if strings.Contains(string(trimmed.Sequence), "N") {
		t.Error("Sequence still contains N after truncation")
	}

	if !strings.HasPrefix(string(trimmed.Sequence), "ACGT") {
		t.Error("Expected sequence to keep bases before N")
	}
}

func TestTrimRecordHighQuality(t *testing.T) {
	record := &fastq.Record{
		ID:       "test",
		Sequence: []byte("ACGTACGTACGTACGT"),
		Quality:  []byte("IIIIIIIIIIIIIIII"),
	}

	opts := DefaultTrimOptions()

	trimmed := trimRecord(record, opts)

	// High quality read should not be trimmed
	if string(trimmed.Sequence) != string(record.Sequence) {
		t.Error("High quality sequence was unexpectedly trimmed")
	}
}

func TestTrimRecordLowQualityEntire(t *testing.T) {
	record := &fastq.Record{
		ID:       "test",
		Sequence: []byte("NNNNNNNN"),
		Quality:  []byte("########"),
	}

	opts := DefaultTrimOptions()
	opts.TruncateN = true

	trimmed := trimRecord(record, opts)

	// Entire sequence is low quality or N, should be empty
	if len(trimmed.Sequence) != 0 {
		t.Error("Expected empty sequence for entirely low quality read")
	}
}

func TestTrim5Prime(t *testing.T) {
	// '#' is quality 2, 'I' is quality 40 in Phred+33
	quality := "####IIIIIIII"
	threshold := 30
	windowSize := 3

	start := trim5Prime(quality, threshold, windowSize, len(quality))

	// Should find start position after low quality region
	if start >= len(quality)-4 {
		t.Errorf("Expected to find good quality region, got start=%d", start)
	}
}

func TestTrim3Prime(t *testing.T) {
	// '#' is quality 2, 'I' is quality 40 in Phred+33
	quality := "IIIIIIII####"
	threshold := 30
	windowSize := 3

	end := trim3Prime(quality, threshold, windowSize, 0, len(quality))

	// Should find end position before low quality region
	if end > len(quality)-3 {
		t.Errorf("Expected to trim low quality region, got end=%d", end)
	}
}

func TestTrimOptionsDefault(t *testing.T) {
	opts := DefaultTrimOptions()

	if opts.QualThreshold != 20 {
		t.Errorf("Expected default quality threshold 20, got %d", opts.QualThreshold)
	}

	if opts.LengthThreshold != 20 {
		t.Errorf("Expected default length threshold 20, got %d", opts.LengthThreshold)
	}

	if opts.WindowSize != 10 {
		t.Errorf("Expected default window size 10, got %d", opts.WindowSize)
	}

	if opts.NoFivePrime {
		t.Error("Expected NoFivePrime to be false by default")
	}

	if opts.TruncateN {
		t.Error("Expected TruncateN to be false by default")
	}
}

func TestTrimStatsAccumulation(t *testing.T) {
	input := `@read1
ACGTACGTACGTACGT
+
IIIIIIII########
@read2
ACGT
+
####
@read3
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
`

	var output bytes.Buffer
	opts := DefaultTrimOptions()
	opts.QualThreshold = 30
	opts.LengthThreshold = 10

	stats, err := TrimSingleEnd(strings.NewReader(input), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd failed: %v", err)
	}

	if stats.TotalReads != 3 {
		t.Errorf("Expected 3 total reads, got %d", stats.TotalReads)
	}

	if stats.TotalBases <= 0 {
		t.Error("Expected positive total bases")
	}

	// At least one read should be discarded (read2 is too short)
	if stats.DiscardedReads == 0 {
		t.Error("Expected at least one discarded read")
	}
}

func TestCustomWindowSize(t *testing.T) {
	input := `@read1
ACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIII############
`

	var output1, output2 bytes.Buffer

	// Test with default window size (10)
	opts1 := DefaultTrimOptions()
	opts1.QualThreshold = 30
	opts1.WindowSize = 10

	stats1, err := TrimSingleEnd(strings.NewReader(input), &output1, fastq.Phred33, opts1)
	if err != nil {
		t.Fatalf("TrimSingleEnd with default window failed: %v", err)
	}

	// Test with smaller window size (5)
	opts2 := DefaultTrimOptions()
	opts2.QualThreshold = 30
	opts2.WindowSize = 5

	stats2, err := TrimSingleEnd(strings.NewReader(input), &output2, fastq.Phred33, opts2)
	if err != nil {
		t.Fatalf("TrimSingleEnd with small window failed: %v", err)
	}

	// Both should process the same number of reads
	if stats1.TotalReads != stats2.TotalReads {
		t.Error("Window size should not affect total reads")
	}

	// Different window sizes may result in different trimming
	// Just verify both ran successfully
	if stats1.TotalBases == 0 || stats2.TotalBases == 0 {
		t.Error("Expected positive total bases for both window sizes")
	}
}

func TestProgressReporting(t *testing.T) {
	// Generate a larger input for progress testing
	var inputBuilder strings.Builder
	for i := 0; i < 100; i++ {
		inputBuilder.WriteString(fmt.Sprintf("@read%d\n", i))
		inputBuilder.WriteString("ACGTACGTACGTACGTACGTACGT\n")
		inputBuilder.WriteString("+\n")
		inputBuilder.WriteString("IIIIIIIIIIIIIIIIIIIIIIII\n")
	}

	var output bytes.Buffer
	opts := DefaultTrimOptions()
	opts.Progress = true // Enable progress reporting

	stats, err := TrimSingleEnd(strings.NewReader(inputBuilder.String()), &output, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimSingleEnd with progress failed: %v", err)
	}

	if stats.TotalReads != 100 {
		t.Errorf("Expected 100 total reads, got %d", stats.TotalReads)
	}
}

func TestQualityRecalibration(t *testing.T) {
	input := `@read1
ACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIII
`

	var output1, output2 bytes.Buffer

	// Test without recalibration
	opts1 := DefaultTrimOptions()
	opts1.Recalibrate = false

	stats1, err := TrimSingleEnd(strings.NewReader(input), &output1, fastq.Phred33, opts1)
	if err != nil {
		t.Fatalf("TrimSingleEnd without recalibration failed: %v", err)
	}

	// Test with recalibration
	opts2 := DefaultTrimOptions()
	opts2.Recalibrate = true

	stats2, err := TrimSingleEnd(strings.NewReader(input), &output2, fastq.Phred33, opts2)
	if err != nil {
		t.Fatalf("TrimSingleEnd with recalibration failed: %v", err)
	}

	// Both should process the same number of reads
	if stats1.TotalReads != stats2.TotalReads {
		t.Error("Recalibration should not affect total reads")
	}

	// High quality reads should pass regardless of recalibration
	if stats1.DiscardedReads > 0 || stats2.DiscardedReads > 0 {
		t.Error("High quality read should not be discarded")
	}
}

func TestRecalibrateRecord(t *testing.T) {
	record := &fastq.Record{
		ID:       "test",
		Sequence: []byte("AAAAACGTACGTACGTACGT"),
		Quality:  []byte("IIIIIIIIIIIIIIIIIIII"),
	}

	recalibrated := recalibrateRecord(record, fastq.Phred33)

	if len(recalibrated.Quality) != len(record.Quality) {
		t.Error("Recalibrated quality should have same length")
	}

	if string(recalibrated.Sequence) != string(record.Sequence) {
		t.Error("Recalibration should not change sequence")
	}

	// Quality values may change but should still be valid ASCII
	for i, q := range recalibrated.Quality {
		if q < 33 || q > 126 {
			t.Errorf("Invalid quality score at position %d: %d", i, q)
		}
	}
}
