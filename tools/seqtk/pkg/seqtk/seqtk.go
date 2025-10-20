// Package seqtk provides core functionality for sequence processing.
// This is a Go reimplementation of seqtk, a fast FASTA/Q processor.
package seqtk

import (
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// Stats represents sequence statistics.
type Stats struct {
	NumSequences int
	TotalBases   int64
	MinLength    int
	MaxLength    int
	AvgLength    float64
	AvgQuality   float64 // For FASTQ only
	GCContent    float64
}

// CalculateFastaStats calculates statistics for a FASTA file.
func CalculateFastaStats(r io.Reader) (*Stats, error) {
	reader := fasta.NewReader(r)
	records, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	if len(records) == 0 {
		return &Stats{}, nil
	}

	stats := &Stats{
		NumSequences: len(records),
		MinLength:    records[0].Length(),
		MaxLength:    records[0].Length(),
	}

	var totalBases int64
	var totalGC int64

	for _, record := range records {
		length := record.Length()
		totalBases += int64(length)

		if length < stats.MinLength {
			stats.MinLength = length
		}
		if length > stats.MaxLength {
			stats.MaxLength = length
		}

		// Count GC
		for _, b := range record.Sequence {
			if b == 'G' || b == 'C' || b == 'g' || b == 'c' {
				totalGC++
			}
		}
	}

	stats.TotalBases = totalBases
	stats.AvgLength = float64(totalBases) / float64(len(records))
	stats.GCContent = float64(totalGC) / float64(totalBases) * 100

	return stats, nil
}

// CalculateFastqStats calculates statistics for a FASTQ file.
func CalculateFastqStats(r io.Reader, encoding fastq.QualityEncoding) (*Stats, error) {
	reader := fastq.NewReader(r, encoding)
	records, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	if len(records) == 0 {
		return &Stats{}, nil
	}

	stats := &Stats{
		NumSequences: len(records),
		MinLength:    records[0].Length(),
		MaxLength:    records[0].Length(),
	}

	var totalBases int64
	var totalGC int64
	var totalQuality float64

	for _, record := range records {
		length := record.Length()
		totalBases += int64(length)

		if length < stats.MinLength {
			stats.MinLength = length
		}
		if length > stats.MaxLength {
			stats.MaxLength = length
		}

		// Count GC
		for _, b := range record.Sequence {
			if b == 'G' || b == 'C' || b == 'g' || b == 'c' {
				totalGC++
			}
		}

		// Calculate average quality
		totalQuality += record.AverageQuality(encoding)
	}

	stats.TotalBases = totalBases
	stats.AvgLength = float64(totalBases) / float64(len(records))
	stats.GCContent = float64(totalGC) / float64(totalBases) * 100
	stats.AvgQuality = totalQuality / float64(len(records))

	return stats, nil
}

// ConvertFastqToFasta converts a FASTQ file to FASTA format.
func ConvertFastqToFasta(input io.Reader, output io.Writer, encoding fastq.QualityEncoding) error {
	fqReader := fastq.NewReader(input, encoding)
	faWriter := fasta.NewWriter(output, 80)

	for {
		record, err := fqReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		faRecord := &fasta.Record{
			ID:          record.ID,
			Description: record.Description,
			Sequence:    record.Sequence,
		}

		if err := faWriter.Write(faRecord); err != nil {
			return err
		}
	}

	return faWriter.Flush()
}

// ReverseComplement generates reverse complement of sequences.
func ReverseComplement(input io.Reader, output io.Writer, isFastq bool, encoding fastq.QualityEncoding) error {
	if isFastq {
		return reverseComplementFastq(input, output, encoding)
	}
	return reverseComplementFasta(input, output)
}

func reverseComplementFasta(input io.Reader, output io.Writer) error {
	reader := fasta.NewReader(input)
	writer := fasta.NewWriter(output, 80)

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		rc := record.ReverseComplement()
		if err := writer.Write(rc); err != nil {
			return err
		}
	}

	return writer.Flush()
}

func reverseComplementFastq(input io.Reader, output io.Writer, encoding fastq.QualityEncoding) error {
	reader := fastq.NewReader(input, encoding)
	writer := fastq.NewWriter(output, encoding)

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		rc := record.ReverseComplement()
		if err := writer.Write(rc); err != nil {
			return err
		}
	}

	return writer.Flush()
}

// Sample randomly samples a fraction of sequences.
func Sample(input io.Reader, output io.Writer, fraction float64, isFastq bool, encoding fastq.QualityEncoding) error {
	if fraction <= 0 || fraction > 1 {
		return fmt.Errorf("fraction must be between 0 and 1")
	}

	if isFastq {
		return sampleFastq(input, output, fraction, encoding)
	}
	return sampleFasta(input, output, fraction)
}

func sampleFasta(input io.Reader, output io.Writer, fraction float64) error {
	reader := fasta.NewReader(input)
	writer := fasta.NewWriter(output, 80)

	count := 0
	written := 0

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		count++
		// Simple deterministic sampling: write every Nth record
		if float64(written)/float64(count) < fraction {
			if err := writer.Write(record); err != nil {
				return err
			}
			written++
		}
	}

	return writer.Flush()
}

func sampleFastq(input io.Reader, output io.Writer, fraction float64, encoding fastq.QualityEncoding) error {
	reader := fastq.NewReader(input, encoding)
	writer := fastq.NewWriter(output, encoding)

	count := 0
	written := 0

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		count++
		// Simple deterministic sampling: write every Nth record
		if float64(written)/float64(count) < fraction {
			if err := writer.Write(record); err != nil {
				return err
			}
			written++
		}
	}

	return writer.Flush()
}

// TrimQuality trims sequences based on quality threshold.
func TrimQuality(input io.Reader, output io.Writer, threshold int, encoding fastq.QualityEncoding) error {
	reader := fastq.NewReader(input, encoding)
	writer := fastq.NewWriter(output, encoding)

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		trimmed := record.Trim(threshold, encoding)
		// Only write if trimmed sequence is not empty
		if len(trimmed.Sequence) > 0 {
			if err := writer.Write(trimmed); err != nil {
				return err
			}
		}
	}

	return writer.Flush()
}

// GetFileType determines if a file is FASTA or FASTQ.
func GetFileType(filename string) (bool, error) {
	file, err := os.Open(filename)
	if err != nil {
		return false, err
	}
	defer file.Close()

	buf := make([]byte, 1)
	_, err = file.Read(buf)
	if err != nil {
		return false, err
	}

	// FASTQ starts with '@', FASTA with '>'
	return buf[0] == '@', nil
}
