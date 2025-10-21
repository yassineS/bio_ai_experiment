// Package fastp provides all-in-one preprocessing for FASTQ files.
// It combines quality filtering, adapter trimming, and various other preprocessing steps.
package fastp

import (
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// ProcessOptions contains all preprocessing parameters.
type ProcessOptions struct {
	// Adapter trimming
	Adapter3         string
	Adapter5         string
	DetectAdapter    bool
	
	// Quality filtering
	QualThreshold    int
	MinLength        int
	MaxLength        int
	QualPercent      int  // Percentage of bases that must meet quality threshold
	
	// Complexity filtering
	LowComplexity    bool
	ComplexityThreshold float64
	
	// Poly-tail trimming
	TrimPolyG        bool
	TrimPolyX        bool
	PolyGMinLen      int
	
	// N filtering
	MaxNCount        int
	MaxNPercent      float64
	
	// Length filtering
	LengthRequired   int
	LengthLimit      int
}

// DefaultProcessOptions returns default processing options.
func DefaultProcessOptions() ProcessOptions {
	return ProcessOptions{
		Adapter3:            "",
		Adapter5:            "",
		DetectAdapter:       false,
		QualThreshold:       15,
		MinLength:           15,
		MaxLength:           0, // no limit
		QualPercent:         40,
		LowComplexity:       false,
		ComplexityThreshold: 0.3,
		TrimPolyG:           false,
		TrimPolyX:           false,
		PolyGMinLen:         10,
		MaxNCount:           5,
		MaxNPercent:         20.0,
		LengthRequired:      15,
		LengthLimit:         0,
	}
}

// ProcessStats tracks preprocessing statistics.
type ProcessStats struct {
	TotalReads          int
	TotalBases          int64
	CleanReads          int
	CleanBases          int64
	LowQualityReads     int
	TooShortReads       int
	TooLongReads        int
	TooManyNReads       int
	AdapterTrimmedReads int
	AdapterTrimmedBases int64
	PolyGTrimmedReads   int
	PolyGTrimmedBases   int64
}

// ProcessPairedEnd processes paired-end FASTQ reads with all filters.
func ProcessPairedEnd(input1, input2 io.Reader, output1, output2 io.Writer, encoding fastq.QualityEncoding, opts ProcessOptions) (*ProcessStats, error) {
	reader1 := fastq.NewReader(input1, encoding)
	reader2 := fastq.NewReader(input2, encoding)
	writer1 := fastq.NewWriter(output1, encoding)
	writer2 := fastq.NewWriter(output2, encoding)
	
	stats := &ProcessStats{}
	
	for {
		// Read both pairs
		record1, err1 := reader1.Read()
		record2, err2 := reader2.Read()
		
		// Check for EOF
		if err1 == io.EOF && err2 == io.EOF {
			break
		}
		if err1 == io.EOF || err2 == io.EOF {
			return stats, fmt.Errorf("paired files have different number of reads")
		}
		if err1 != nil {
			return stats, fmt.Errorf("error reading read1: %w", err1)
		}
		if err2 != nil {
			return stats, fmt.Errorf("error reading read2: %w", err2)
		}
		
		stats.TotalReads += 2
		stats.TotalBases += int64(len(record1.Sequence) + len(record2.Sequence))
		
		// Process both records
		processed1, pass1 := processRecord(record1, opts, stats, encoding)
		processed2, pass2 := processRecord(record2, opts, stats, encoding)
		
		// Both must pass for the pair to be kept
		if pass1 && pass2 {
			if err := writer1.Write(processed1); err != nil {
				return stats, fmt.Errorf("error writing read1: %w", err)
			}
			if err := writer2.Write(processed2); err != nil {
				return stats, fmt.Errorf("error writing read2: %w", err)
			}
			stats.CleanReads += 2
			stats.CleanBases += int64(len(processed1.Sequence) + len(processed2.Sequence))
		}
	}
	
	// Flush writers
	if err := writer1.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output1: %w", err)
	}
	if err := writer2.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output2: %w", err)
	}
	
	return stats, nil
}
func ProcessSingleEnd(input io.Reader, output io.Writer, encoding fastq.QualityEncoding, opts ProcessOptions) (*ProcessStats, error) {
	reader := fastq.NewReader(input, encoding)
	writer := fastq.NewWriter(output, encoding)
	
	stats := &ProcessStats{}
	
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
		
		// Process the record
		processed, pass := processRecord(record, opts, stats, encoding)
		
		// Write if passed all filters
		if pass {
			if err := writer.Write(processed); err != nil {
				return stats, fmt.Errorf("error writing FASTQ: %w", err)
			}
			stats.CleanReads++
			stats.CleanBases += int64(len(processed.Sequence))
		}
	}
	
	// Flush writer
	if err := writer.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output: %w", err)
	}
	
	return stats, nil
}

// processRecord applies all processing steps to a single record.
func processRecord(record *fastq.Record, opts ProcessOptions, stats *ProcessStats, encoding fastq.QualityEncoding) (*fastq.Record, bool) {
	seq := string(record.Sequence)
	qual := record.Quality
	
	start := 0
	end := len(seq)
	
	// Step 1: Trim adapters if specified
	if opts.Adapter5 != "" {
		pos := findAdapter(seq, opts.Adapter5)
		if pos >= 0 {
			start = pos + len(opts.Adapter5)
			stats.AdapterTrimmedReads++
			stats.AdapterTrimmedBases += int64(pos + len(opts.Adapter5))
		}
	}
	
	if opts.Adapter3 != "" {
		pos := findAdapter(seq[start:], opts.Adapter3)
		if pos >= 0 {
			end = start + pos
			stats.AdapterTrimmedReads++
			stats.AdapterTrimmedBases += int64(len(seq) - end)
		}
	}
	
	// Step 2: Trim poly-G tails if enabled
	if opts.TrimPolyG {
		polyLen := countPolyTail(seq[start:end], 'G')
		if polyLen >= opts.PolyGMinLen {
			end -= polyLen
			stats.PolyGTrimmedReads++
			stats.PolyGTrimmedBases += int64(polyLen)
		}
	}
	
	// Step 3: Trim poly-X tails if enabled
	if opts.TrimPolyX {
		for _, base := range []byte{'A', 'T', 'C'} {
			polyLen := countPolyTail(seq[start:end], rune(base))
			if polyLen >= opts.PolyGMinLen {
				end -= polyLen
			}
		}
	}
	
	// Step 4: Quality-based trimming
	if opts.QualThreshold > 0 {
		start, end = trimByQuality(qual[start:end], opts.QualThreshold, start, end, encoding)
	}
	
	// Check if read is too short after trimming
	if end-start < opts.MinLength || end-start < opts.LengthRequired {
		stats.TooShortReads++
		return nil, false
	}
	
	// Check if read is too long
	if opts.MaxLength > 0 && end-start > opts.MaxLength {
		stats.TooLongReads++
		return nil, false
	}
	
	if opts.LengthLimit > 0 && end-start > opts.LengthLimit {
		end = start + opts.LengthLimit
	}
	
	// Step 5: Check N content
	nCount := countNs(seq[start:end])
	nPercent := 100.0 * float64(nCount) / float64(end-start)
	
	if nCount > opts.MaxNCount || nPercent > opts.MaxNPercent {
		stats.TooManyNReads++
		return nil, false
	}
	
	// Step 6: Check quality (percentage of bases meeting threshold)
	if opts.QualPercent > 0 {
		qualScores := getQualityScores(qual[start:end], encoding)
		passCount := 0
		for _, q := range qualScores {
			if q >= opts.QualThreshold {
				passCount++
			}
		}
		passPercent := 100.0 * float64(passCount) / float64(len(qualScores))
		if passPercent < float64(opts.QualPercent) {
			stats.LowQualityReads++
			return nil, false
		}
	}
	
	// Step 7: Check complexity if enabled
	if opts.LowComplexity {
		complexity := calculateComplexity(seq[start:end])
		if complexity < opts.ComplexityThreshold {
			// Low complexity read - discard
			return nil, false
		}
	}
	
	// Create processed record
	return &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    record.Sequence[start:end],
		Quality:     record.Quality[start:end],
	}, true
}

// findAdapter finds the position of an adapter in a sequence.
func findAdapter(seq string, adapter string) int {
	return strings.Index(seq, adapter)
}

// countPolyTail counts the length of a poly-X tail at the end of a sequence.
func countPolyTail(seq string, base rune) int {
	count := 0
	for i := len(seq) - 1; i >= 0 && rune(seq[i]) == base; i-- {
		count++
	}
	return count
}

// countNs counts the number of N bases in a sequence.
func countNs(seq string) int {
	count := 0
	for _, b := range seq {
		if b == 'N' || b == 'n' {
			count++
		}
	}
	return count
}

// trimByQuality trims low-quality regions from both ends.
func trimByQuality(quality []byte, threshold int, start, end int, encoding fastq.QualityEncoding) (int, int) {
	offset := 33
	if encoding == fastq.Phred64 {
		offset = 64
	}
	
	// Trim from 3' end
	for end > start && int(quality[end-start-1])-offset < threshold {
		end--
	}
	
	// Trim from 5' end
	for start < end && int(quality[0])-offset < threshold {
		start++
		quality = quality[1:]
	}
	
	return start, end
}

// getQualityScores converts ASCII-encoded quality scores to numeric values.
func getQualityScores(quality []byte, encoding fastq.QualityEncoding) []int {
	scores := make([]int, len(quality))
	offset := 33
	if encoding == fastq.Phred64 {
		offset = 64
	}
	for i, q := range quality {
		scores[i] = int(q) - offset
	}
	return scores
}

// calculateComplexity calculates sequence complexity (0-1, higher is more complex).
func calculateComplexity(seq string) float64 {
	if len(seq) == 0 {
		return 0
	}
	
	// Simple complexity measure: unique 2-mers / total 2-mers
	if len(seq) < 2 {
		return 1.0
	}
	
	kmers := make(map[string]bool)
	for i := 0; i < len(seq)-1; i++ {
		kmer := seq[i : i+2]
		kmers[kmer] = true
	}
	
	return float64(len(kmers)) / float64(len(seq)-1)
}
