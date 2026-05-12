// Package sickle provides quality-based trimming for FASTQ files using a sliding window approach.
package sickle

import (
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// TrimOptions contains parameters for quality-based trimming.
type TrimOptions struct {
	QualThreshold   int  // Minimum quality score threshold
	LengthThreshold int  // Minimum length threshold after trimming
	NoFivePrime     bool // Don't trim 5' end
	TruncateN       bool // Truncate at first N
	WindowSize      int  // Size of sliding window for quality assessment
	Progress        bool // Show progress reporting
	Recalibrate     bool // Recalibrate quality scores
}

// DefaultTrimOptions returns default trimming options.
func DefaultTrimOptions() TrimOptions {
	return TrimOptions{
		QualThreshold:   20,
		LengthThreshold: 20,
		NoFivePrime:     false,
		TruncateN:       false,
		WindowSize:      10,
	}
}

// TrimStats tracks trimming statistics.
type TrimStats struct {
	TotalReads     int
	TrimmedReads   int
	DiscardedReads int
	TotalBases     int64
	TrimmedBases   int64
}

// TrimSingleEnd trims a single-end FASTQ file based on quality scores.
func TrimSingleEnd(input io.Reader, output io.Writer, encoding fastq.QualityEncoding, opts TrimOptions) (*TrimStats, error) {
	reader := fastq.NewReader(input, encoding)
	writer := fastq.NewWriter(output, encoding)

	stats := &TrimStats{}

	// Progress reporting setup
	var progressCounter int
	const progressInterval = 10000

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return stats, fmt.Errorf("error reading FASTQ: %w", err)
		}

		stats.TotalReads++
		stats.TotalBases += int64(len(record.Sequence))

		// Apply recalibration if requested
		if opts.Recalibrate {
			record = recalibrateRecord(record, encoding)
		}

		// Apply trimming
		trimmed := trimRecord(record, opts)

		// Check if read passes length threshold
		if len(trimmed.Sequence) >= opts.LengthThreshold {
			if err := writer.Write(trimmed); err != nil {
				return stats, fmt.Errorf("error writing FASTQ: %w", err)
			}
			if len(trimmed.Sequence) < len(record.Sequence) {
				stats.TrimmedReads++
				stats.TrimmedBases += int64(len(record.Sequence) - len(trimmed.Sequence))
			}
		} else {
			stats.DiscardedReads++
		}

		// Progress reporting
		if opts.Progress && stats.TotalReads%progressInterval == 0 {
			progressCounter++
			fmt.Fprintf(os.Stderr, "\rProcessed %d reads...", stats.TotalReads)
		}
	}

	// Clear progress line
	if opts.Progress {
		fmt.Fprintf(os.Stderr, "\r")
	}

	// Flush writer
	if err := writer.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output: %w", err)
	}

	return stats, nil
}

// TrimPairedEnd trims paired-end FASTQ files, maintaining synchronization.
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

	// Progress reporting setup
	const progressInterval = 10000

	for {
		record1, err1 := reader1.Read()
		record2, err2 := reader2.Read()

		if err1 == io.EOF && err2 == io.EOF {
			break
		}
		if err1 != nil && err1 != io.EOF {
			return stats, fmt.Errorf("error reading first FASTQ: %w", err1)
		}
		if err2 != nil && err2 != io.EOF {
			return stats, fmt.Errorf("error reading second FASTQ: %w", err2)
		}
		if (err1 == io.EOF) != (err2 == io.EOF) {
			return stats, fmt.Errorf("paired-end files have different number of reads")
		}

		stats.TotalReads += 2
		stats.TotalBases += int64(len(record1.Sequence) + len(record2.Sequence))

		// Apply recalibration if requested
		if opts.Recalibrate {
			record1 = recalibrateRecord(record1, encoding)
			record2 = recalibrateRecord(record2, encoding)
		}

		// Trim both reads
		trimmed1 := trimRecord(record1, opts)
		trimmed2 := trimRecord(record2, opts)

		pass1 := len(trimmed1.Sequence) >= opts.LengthThreshold
		pass2 := len(trimmed2.Sequence) >= opts.LengthThreshold

		// Both reads pass - write to paired output
		if pass1 && pass2 {
			if err := writer1.Write(trimmed1); err != nil {
				return stats, fmt.Errorf("error writing first FASTQ: %w", err)
			}
			if err := writer2.Write(trimmed2); err != nil {
				return stats, fmt.Errorf("error writing second FASTQ: %w", err)
			}

			if len(trimmed1.Sequence) < len(record1.Sequence) {
				stats.TrimmedReads++
				stats.TrimmedBases += int64(len(record1.Sequence) - len(trimmed1.Sequence))
			}
			if len(trimmed2.Sequence) < len(record2.Sequence) {
				stats.TrimmedReads++
				stats.TrimmedBases += int64(len(record2.Sequence) - len(trimmed2.Sequence))
			}
		} else if writerSingle != nil {
			// One read passes - write to single output if available
			if pass1 {
				if err := writerSingle.Write(trimmed1); err != nil {
					return stats, fmt.Errorf("error writing single FASTQ: %w", err)
				}
				if len(trimmed1.Sequence) < len(record1.Sequence) {
					stats.TrimmedReads++
					stats.TrimmedBases += int64(len(record1.Sequence) - len(trimmed1.Sequence))
				}
			}
			if pass2 {
				if err := writerSingle.Write(trimmed2); err != nil {
					return stats, fmt.Errorf("error writing single FASTQ: %w", err)
				}
				if len(trimmed2.Sequence) < len(record2.Sequence) {
					stats.TrimmedReads++
					stats.TrimmedBases += int64(len(record2.Sequence) - len(trimmed2.Sequence))
				}
			}
			if !pass1 {
				stats.DiscardedReads++
			}
			if !pass2 {
				stats.DiscardedReads++
			}
		} else {
			// Discard both if either fails and no single output
			stats.DiscardedReads += 2
		}

		// Progress reporting
		if opts.Progress && stats.TotalReads%progressInterval == 0 {
			fmt.Fprintf(os.Stderr, "\rProcessed %d reads...", stats.TotalReads)
		}
	}

	// Clear progress line
	if opts.Progress {
		fmt.Fprintf(os.Stderr, "\r")
	}

	// Flush all writers
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

// trimRecord applies quality-based trimming to a FASTQ record.
func trimRecord(record *fastq.Record, opts TrimOptions) *fastq.Record {
	seq := record.Sequence
	qual := record.Quality

	if len(seq) == 0 {
		return record
	}

	start := 0
	end := len(seq)

	// Truncate at first N if requested
	if opts.TruncateN {
		for i, base := range seq {
			if base == 'N' || base == 'n' {
				end = i
				break
			}
		}
		if end == 0 {
			// Entire sequence is N or starts with N
			return &fastq.Record{
				ID:          record.ID,
				Description: record.Description,
				Sequence:    []byte{},
				Quality:     []byte{},
			}
		}
	}

	// Trim 5' end (left side) unless disabled
	if !opts.NoFivePrime {
		start = trim5Prime(string(qual), opts.QualThreshold, opts.WindowSize, end)
	}

	// Trim 3' end (right side)
	if start < end {
		end = trim3Prime(string(qual), opts.QualThreshold, opts.WindowSize, start, end)
	}

	// Return trimmed record
	if start >= end || end <= 0 {
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
		Sequence:    seq[start:end],
		Quality:     qual[start:end],
	}
}

// trim5Prime finds the start position by trimming from 5' end.
// Assumes Phred+33 encoding for quality scores
func trim5Prime(quality string, threshold, windowSize, maxPos int) int {
	if len(quality) == 0 || maxPos == 0 {
		return 0
	}

	// Convert quality string to scores (Phred+33)
	scores := make([]int, len(quality))
	for i := 0; i < len(quality); i++ {
		scores[i] = int(quality[i]) - 33
	}

	windowSum := 0
	window := 0

	// Calculate initial window
	for i := 0; i < windowSize && i < maxPos && i < len(scores); i++ {
		windowSum += scores[i]
		window++
	}

	// Slide window from left to right
	for i := 0; i <= maxPos-window; i++ {
		if window == 0 {
			break
		}
		avgQual := windowSum / window
		if avgQual >= threshold {
			return i
		}

		// Slide window
		if i+window < maxPos && i+window < len(scores) {
			windowSum -= scores[i]
			windowSum += scores[i+window]
		} else {
			window--
		}
	}

	return 0
}

// trim3Prime finds the end position by trimming from 3' end.
// Assumes Phred+33 encoding for quality scores
func trim3Prime(quality string, threshold, windowSize, minPos, maxPos int) int {
	if len(quality) == 0 || maxPos <= minPos {
		return maxPos
	}

	// Convert quality string to scores (Phred+33)
	scores := make([]int, len(quality))
	for i := 0; i < len(quality); i++ {
		scores[i] = int(quality[i]) - 33
	}

	windowSum := 0
	window := 0

	// Calculate initial window from the right
	startIdx := maxPos - windowSize
	if startIdx < minPos {
		startIdx = minPos
	}

	for i := startIdx; i < maxPos && i < len(scores); i++ {
		windowSum += scores[i]
		window++
	}

	// Slide window from right to left
	for i := maxPos; i >= minPos+window; i-- {
		if window == 0 {
			break
		}
		avgQual := windowSum / window
		if avgQual >= threshold {
			return i
		}

		// Slide window
		if i-1 >= minPos && i-1 < len(scores) {
			windowSum -= scores[i-1]
			if i-window-1 >= minPos && i-window-1 < len(scores) {
				windowSum += scores[i-window-1]
			} else {
				window--
			}
		} else {
			break
		}
	}

	return maxPos
}

// recalibrateRecord recalibrates quality scores using empirical base quality score recalibration.
// This is a simplified version that adjusts quality scores based on sequence context.
func recalibrateRecord(record *fastq.Record, encoding fastq.QualityEncoding) *fastq.Record {
	if len(record.Quality) == 0 {
		return record
	}

	// Get quality offset
	offset := 33
	if encoding == fastq.Phred64 {
		offset = 64
	}

	// Create a copy of the quality scores
	newQuality := make([]byte, len(record.Quality))
	copy(newQuality, record.Quality)

	// Simple recalibration: adjust quality scores based on position and context
	// This is a simplified algorithm - real recalibration would use machine learning
	for i := 0; i < len(newQuality); i++ {
		currentQual := int(newQuality[i]) - offset

		// Position-based adjustment: quality tends to degrade toward read ends
		positionFactor := 1.0
		readLength := len(newQuality)
		if i < readLength/10 || i > 9*readLength/10 {
			// First and last 10% of read: slight quality penalty
			positionFactor = 0.95
		}

		// Context-based adjustment: homopolymers are more error-prone
		if i > 0 && i < len(record.Sequence)-1 {
			prevBase := record.Sequence[i-1]
			currBase := record.Sequence[i]
			nextBase := record.Sequence[i+1]

			// Penalize homopolymer runs
			if prevBase == currBase || currBase == nextBase {
				positionFactor *= 0.95
			}
		}

		// Apply adjustments
		adjustedQual := int(float64(currentQual) * positionFactor)

		// Clamp to valid range
		if adjustedQual < 0 {
			adjustedQual = 0
		}
		maxQual := 93
		if encoding == fastq.Phred64 {
			maxQual = 62
		}
		if adjustedQual > maxQual {
			adjustedQual = maxQual
		}

		newQuality[i] = byte(adjustedQual + offset)
	}

	return &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    record.Sequence,
		Quality:     newQuality,
	}
}
