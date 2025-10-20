package seqtk

import (
	"bytes"
	"os"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

func TestCalculateFastaStats(t *testing.T) {
	fasta := `>seq1
ACGTACGT
>seq2
GCGCGCGCGC
>seq3
ATATATAT
`
	r := strings.NewReader(fasta)
	stats, err := CalculateFastaStats(r)
	if err != nil {
		t.Fatalf("CalculateFastaStats failed: %v", err)
	}

	if stats.NumSequences != 3 {
		t.Errorf("Expected 3 sequences, got %d", stats.NumSequences)
	}

	if stats.TotalBases != 26 {
		t.Errorf("Expected 26 total bases, got %d", stats.TotalBases)
	}

	if stats.MinLength != 8 {
		t.Errorf("Expected min length 8, got %d", stats.MinLength)
	}

	if stats.MaxLength != 10 {
		t.Errorf("Expected max length 10, got %d", stats.MaxLength)
	}

	expectedGC := (10.0 + 4.0) / 26.0 * 100.0 // 10 GC from seq2, 4 from seq1
	if diff := stats.GCContent - expectedGC; diff > 0.01 || diff < -0.01 {
		t.Errorf("Expected GC content %.2f%%, got %.2f%%", expectedGC, stats.GCContent)
	}
}

func TestCalculateFastqStats(t *testing.T) {
	fastqData := `@read1
ACGTACGT
+
IIIIIIII
@read2
GCGCGCGC
+
IIIIIIII
`
	r := strings.NewReader(fastqData)
	stats, err := CalculateFastqStats(r, fastq.Phred33)
	if err != nil {
		t.Fatalf("CalculateFastqStats failed: %v", err)
	}

	if stats.NumSequences != 2 {
		t.Errorf("Expected 2 sequences, got %d", stats.NumSequences)
	}

	if stats.TotalBases != 16 {
		t.Errorf("Expected 16 total bases, got %d", stats.TotalBases)
	}
}

func TestConvertFastqToFasta(t *testing.T) {
	fastqData := `@read1 description
ACGTACGT
+
IIIIIIII
@read2
GCGCGCGC
+
IIIIIIII
`
	r := strings.NewReader(fastqData)
	var buf bytes.Buffer

	err := ConvertFastqToFasta(r, &buf, fastq.Phred33)
	if err != nil {
		t.Fatalf("ConvertFastqToFasta failed: %v", err)
	}

	output := buf.String()
	if !strings.Contains(output, ">read1 description") {
		t.Error("Output doesn't contain expected header")
	}
	if !strings.Contains(output, "ACGTACGT") {
		t.Error("Output doesn't contain expected sequence")
	}
	if strings.Contains(output, "IIIIIIII") {
		t.Error("Output shouldn't contain quality scores")
	}
}

func TestReverseComplementFasta(t *testing.T) {
	fasta := `>seq1
ACGT
`
	r := strings.NewReader(fasta)
	var buf bytes.Buffer

	err := ReverseComplement(r, &buf, false, fastq.Phred33)
	if err != nil {
		t.Fatalf("ReverseComplement failed: %v", err)
	}

	output := buf.String()
	if !strings.Contains(output, "ACGT") {
		t.Error("Expected reverse complement ACGT")
	}
}

func TestSampleFasta(t *testing.T) {
	fasta := `>seq1
ACGT
>seq2
GCTA
>seq3
TGCA
>seq4
CATG
`
	r := strings.NewReader(fasta)
	var buf bytes.Buffer

	err := Sample(r, &buf, 0.5, false, fastq.Phred33)
	if err != nil {
		t.Fatalf("Sample failed: %v", err)
	}

	output := buf.String()
	// Should contain at least one sequence
	if !strings.Contains(output, ">seq") {
		t.Error("Output should contain at least one sequence")
	}
}

func TestTrimQuality(t *testing.T) {
	// Quality scores: ! = 0, I = 40 (Phred+33)
	fastqData := `@read1
ACGTACGT
+
!!!!IIII
`
	r := strings.NewReader(fastqData)
	var buf bytes.Buffer

	err := TrimQuality(r, &buf, 30, fastq.Phred33)
	if err != nil {
		t.Fatalf("TrimQuality failed: %v", err)
	}

	output := buf.String()
	// Should only keep high-quality portion
	if !strings.Contains(output, "@read1") {
		t.Error("Output should contain trimmed read")
	}
}

func TestGetFileType(t *testing.T) {
	// Create temp FASTA file
	fastaContent := []byte(">seq1\nACGT\n")
	fastaFile := "/tmp/test.fasta"
	if err := os.WriteFile(fastaFile, fastaContent, 0644); err != nil {
		t.Fatalf("Failed to create test file: %v", err)
	}
	defer os.Remove(fastaFile)

	isFastq, err := GetFileType(fastaFile)
	if err != nil {
		t.Fatalf("GetFileType failed: %v", err)
	}
	if isFastq {
		t.Error("FASTA file detected as FASTQ")
	}

	// Create temp FASTQ file
	fastqContent := []byte("@seq1\nACGT\n+\nIIII\n")
	fastqFile := "/tmp/test.fastq"
	if err := os.WriteFile(fastqFile, fastqContent, 0644); err != nil {
		t.Fatalf("Failed to create test file: %v", err)
	}
	defer os.Remove(fastqFile)

	isFastq, err = GetFileType(fastqFile)
	if err != nil {
		t.Fatalf("GetFileType failed: %v", err)
	}
	if !isFastq {
		t.Error("FASTQ file not detected as FASTQ")
	}
}
