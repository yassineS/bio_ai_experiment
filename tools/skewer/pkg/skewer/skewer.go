// Package skewer provides adapter trimming for FASTQ files.
// It detects and removes adapter sequences from the 3' and 5' ends of reads.
package skewer

import (
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// TrimOptions contains parameters for adapter trimming.
type TrimOptions struct {
	Adapter3         string  // 3' adapter sequence
	Adapter5         string  // 5' adapter sequence
	MinLength        int     // Minimum read length after trimming
	QualThreshold    int     // Quality threshold for trimming
	MinOverlap       int     // Minimum overlap for adapter detection
	ErrorRate        float64 // Maximum error rate for adapter matching
	TrimBothEnds     bool    // Trim adapters from both ends
}

// DefaultTrimOptions returns default trimming options.
func DefaultTrimOptions() TrimOptions {
	return TrimOptions{
		Adapter3:      "",
		Adapter5:      "",
		MinLength:     18,
		QualThreshold: 0,
		MinOverlap:    3,
		ErrorRate:     0.1,
		TrimBothEnds:  false,
	}
}

// TrimStats tracks adapter trimming statistics.
type TrimStats struct {
	TotalReads      int
	TrimmedReads    int
	AdapterFound3   int
	AdapterFound5   int
	DiscardedReads  int
	TotalBases      int64
	TrimmedBases    int64
}

// TrimSingleEnd trims adapters from single-end FASTQ reads.
func TrimSingleEnd(input io.Reader, output io.Writer, encoding fastq.QualityEncoding, opts TrimOptions) (*TrimStats, error) {
	reader := fastq.NewReader(input, encoding)
	writer := fastq.NewWriter(output, encoding)
	
	stats := &TrimStats{}
	
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return stats, fmt.Errorf("error reading FASTQ: %w", err)
		}
		
		stats.TotalReads++
		originalLength := len(record.Sequence)
		stats.TotalBases += int64(originalLength)
		
		// Trim adapters
		trimmed := trimRecord(record, opts, stats)
		
		// Check if read passes length threshold
		if len(trimmed.Sequence) >= opts.MinLength {
			if err := writer.Write(trimmed); err != nil {
				return stats, fmt.Errorf("error writing FASTQ: %w", err)
			}
			if len(trimmed.Sequence) < originalLength {
				stats.TrimmedReads++
				stats.TrimmedBases += int64(originalLength - len(trimmed.Sequence))
			}
		} else {
			stats.DiscardedReads++
		}
	}
	
	// Flush writer
	if err := writer.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output: %w", err)
	}
	
	return stats, nil
}

// TrimPairedEnd trims adapters from paired-end FASTQ reads.
func TrimPairedEnd(input1, input2 io.Reader, output1, output2, outputSingle io.Writer,
	encoding fastq.QualityEncoding, opts TrimOptions) (*TrimStats, error) {
	
	reader1 := fastq.NewReader(input1, encoding)
	reader2 := fastq.NewReader(input2, encoding)
	writer1 := fastq.NewWriter(output1, encoding)
	writer2 := fastq.NewWriter(output2, encoding)
	
	var writerSingle *fastq.Writer
	if outputSingle != nil {
		writerSingle = fastq.NewWriter(outputSingle, encoding)
	}
	
	stats := &TrimStats{}
	
	for {
		record1, err1 := reader1.Read()
		record2, err2 := reader2.Read()
		
		if err1 == io.EOF || err2 == io.EOF {
			break
		}
		if err1 != nil {
			return stats, fmt.Errorf("error reading first input: %w", err1)
		}
		if err2 != nil {
			return stats, fmt.Errorf("error reading second input: %w", err2)
		}
		
		stats.TotalReads += 2
		originalLen1 := len(record1.Sequence)
		originalLen2 := len(record2.Sequence)
		stats.TotalBases += int64(originalLen1 + originalLen2)
		
		// Trim both reads
		trimmed1 := trimRecord(record1, opts, stats)
		trimmed2 := trimRecord(record2, opts, stats)
		
		// Check if both reads pass length threshold
		pass1 := len(trimmed1.Sequence) >= opts.MinLength
		pass2 := len(trimmed2.Sequence) >= opts.MinLength
		
		if pass1 && pass2 {
			// Both pass - write to paired output
			if err := writer1.Write(trimmed1); err != nil {
				return stats, fmt.Errorf("error writing first output: %w", err)
			}
			if err := writer2.Write(trimmed2); err != nil {
				return stats, fmt.Errorf("error writing second output: %w", err)
			}
			
			if len(trimmed1.Sequence) < originalLen1 || len(trimmed2.Sequence) < originalLen2 {
				stats.TrimmedReads++
				stats.TrimmedBases += int64((originalLen1 - len(trimmed1.Sequence)) +
					(originalLen2 - len(trimmed2.Sequence)))
			}
		} else if writerSingle != nil {
			// One or both fail - write survivors to single output
			if pass1 {
				if err := writerSingle.Write(trimmed1); err != nil {
					return stats, fmt.Errorf("error writing single output: %w", err)
				}
			} else {
				stats.DiscardedReads++
			}
			
			if pass2 {
				if err := writerSingle.Write(trimmed2); err != nil {
					return stats, fmt.Errorf("error writing single output: %w", err)
				}
			} else {
				stats.DiscardedReads++
			}
		} else {
			// Both discarded
			stats.DiscardedReads += 2
		}
	}
	
	// Flush writers
	if err := writer1.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing first output: %w", err)
	}
	if err := writer2.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing second output: %w", err)
	}
	if writerSingle != nil {
		if err := writerSingle.Flush(); err != nil {
			return stats, fmt.Errorf("error flushing single output: %w", err)
		}
	}
	
	return stats, nil
}

// trimRecord trims adapters from a single record.
func trimRecord(record *fastq.Record, opts TrimOptions, stats *TrimStats) *fastq.Record {
	seq := string(record.Sequence)
	qual := record.Quality
	
	start := 0
	end := len(seq)
	
	// Trim 5' adapter if specified
	if opts.Adapter5 != "" {
		pos := findAdapter(seq, opts.Adapter5, opts.MinOverlap, opts.ErrorRate)
		if pos >= 0 {
			// Found 5' adapter - trim from start to end of adapter
			start = pos + len(opts.Adapter5)
			if stats != nil {
				stats.AdapterFound5++
			}
		}
	}
	
	// Trim 3' adapter if specified
	if opts.Adapter3 != "" {
		pos := findAdapter(seq[start:], opts.Adapter3, opts.MinOverlap, opts.ErrorRate)
		if pos >= 0 {
			// Found 3' adapter - trim from adapter position to end
			end = start + pos
			if stats != nil {
				stats.AdapterFound3++
			}
		}
	}
	
	// Apply quality-based trimming if threshold is set
	if opts.QualThreshold > 0 {
		start, end = trimByQuality(qual[start:end], opts.QualThreshold, start, end)
	}
	
	// Create trimmed record
	if start >= end || end-start < opts.MinLength {
		// Return empty record if too short
		return &fastq.Record{
			ID:          record.ID,
			Description: record.Description,
			Sequence:    []byte{},
			Quality:     []byte{},
		}
	}
	
	return &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    record.Sequence[start:end],
		Quality:     record.Quality[start:end],
	}
}

// findAdapter finds the position of an adapter in a sequence with error tolerance.
// Returns -1 if adapter not found.
func findAdapter(seq string, adapter string, minOverlap int, errorRate float64) int {
	if len(adapter) == 0 {
		return -1
	}
	
	// Try full adapter first
	if strings.Contains(seq, adapter) {
		return strings.Index(seq, adapter)
	}
	
	// Try with error tolerance using simple string matching
	// For minimal implementation, we'll use exact matching for now
	maxErrors := int(float64(len(adapter)) * errorRate)
	
	for i := 0; i <= len(seq)-minOverlap; i++ {
		// Check if adapter matches at position i
		matches := 0
		compareLen := min(len(adapter), len(seq)-i)
		
		if compareLen < minOverlap {
			continue
		}
		
		for j := 0; j < compareLen; j++ {
			if seq[i+j] == adapter[j] {
				matches++
			}
		}
		
		errors := compareLen - matches
		if errors <= maxErrors && compareLen >= minOverlap {
			return i
		}
	}
	
	return -1
}

// trimByQuality trims low-quality regions from both ends.
func trimByQuality(quality []byte, threshold int, start, end int) (int, int) {
	// Trim from 3' end
	for end > start && int(quality[end-start-1]) < threshold+33 {
		end--
	}
	
	// Trim from 5' end
	for start < end && int(quality[0]) < threshold+33 {
		start++
		quality = quality[1:]
	}
	
	return start, end
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
