// Package skewer provides adapter trimming for FASTQ files.
// It detects and removes adapter sequences from the 3' and 5' ends of reads.
package skewer

import (
	"encoding/json"
	"fmt"
	"io"
	"sync"
	"time"

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
	AutoDetect       bool    // Auto-detect adapters from reads
	AutoDetectReads  int     // Number of reads to analyze for auto-detection
	ProgressReport   bool    // Report progress during processing
	ProgressInterval int     // Report progress every N reads
	UMILength        int     // UMI length to extract (0 = disabled)
	UMIPosition      string  // UMI position: "5prime" or "3prime"
}

// DefaultTrimOptions returns default trimming options.
func DefaultTrimOptions() TrimOptions {
	return TrimOptions{
		Adapter3:         "",
		Adapter5:         "",
		MinLength:        18,
		QualThreshold:    0,
		MinOverlap:       3,
		ErrorRate:        0.1,
		TrimBothEnds:     false,
		AutoDetect:       false,
		AutoDetectReads:  1000,
		ProgressReport:   false,
		ProgressInterval: 100000,
		UMILength:        0,
		UMIPosition:      "5prime",
	}
}

// TrimStats tracks adapter trimming statistics.
type TrimStats struct {
	TotalReads       int           `json:"total_reads"`
	TrimmedReads     int           `json:"trimmed_reads"`
	AdapterFound3    int           `json:"adapter_found_3"`
	AdapterFound5    int           `json:"adapter_found_5"`
	DiscardedReads   int           `json:"discarded_reads"`
	TotalBases       int64         `json:"total_bases"`
	TrimmedBases     int64         `json:"trimmed_bases"`
	ProcessingTime   time.Duration `json:"processing_time"`
	DetectedAdapter3 string        `json:"detected_adapter_3,omitempty"`
	DetectedAdapter5 string        `json:"detected_adapter_5,omitempty"`
	UMIStats         *UMIStats     `json:"umi_stats,omitempty"`
	mu               sync.Mutex    `json:"-"`
}

// UMIStats tracks UMI/barcode statistics.
type UMIStats struct {
	TotalUMIs       int            `json:"total_umis"`
	UniqueUMIs      int            `json:"unique_umis"`
	UMIDistribution map[string]int `json:"umi_distribution,omitempty"`
}

// ToJSON exports statistics as JSON.
func (s *TrimStats) ToJSON() (string, error) {
	s.mu.Lock()
	defer s.mu.Unlock()

	data, err := json.MarshalIndent(s, "", "  ")
	if err != nil {
		return "", err
	}
	return string(data), nil
}

// ToHTML generates an HTML report of the statistics.
func (s *TrimStats) ToHTML() string {
	s.mu.Lock()
	defer s.mu.Unlock()

	html := `<!DOCTYPE html>
<html>
<head>
    <title>Skewer Trimming Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #333; }
        table { border-collapse: collapse; width: 100%; max-width: 800px; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #4CAF50; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .metric { font-weight: bold; }
        .value { color: #2196F3; }
        .percentage { color: #FF9800; font-size: 0.9em; }
    </style>
</head>
<body>
    <h1>Adapter Trimming Report</h1>
    <table>
        <tr><th>Metric</th><th>Value</th></tr>
`

	html += fmt.Sprintf("        <tr><td class='metric'>Total Reads</td><td class='value'>%d</td></tr>\n", s.TotalReads)

	if s.TotalReads > 0 {
		trimmedPct := 100.0 * float64(s.TrimmedReads) / float64(s.TotalReads)
		html += fmt.Sprintf("        <tr><td class='metric'>Trimmed Reads</td><td class='value'>%d <span class='percentage'>(%.2f%%)</span></td></tr>\n",
			s.TrimmedReads, trimmedPct)

		adapter3Pct := 100.0 * float64(s.AdapterFound3) / float64(s.TotalReads)
		html += fmt.Sprintf("        <tr><td class='metric'>3' Adapters Found</td><td class='value'>%d <span class='percentage'>(%.2f%%)</span></td></tr>\n",
			s.AdapterFound3, adapter3Pct)

		adapter5Pct := 100.0 * float64(s.AdapterFound5) / float64(s.TotalReads)
		html += fmt.Sprintf("        <tr><td class='metric'>5' Adapters Found</td><td class='value'>%d <span class='percentage'>(%.2f%%)</span></td></tr>\n",
			s.AdapterFound5, adapter5Pct)

		discardedPct := 100.0 * float64(s.DiscardedReads) / float64(s.TotalReads)
		html += fmt.Sprintf("        <tr><td class='metric'>Discarded Reads</td><td class='value'>%d <span class='percentage'>(%.2f%%)</span></td></tr>\n",
			s.DiscardedReads, discardedPct)

		keptPct := 100.0 * float64(s.TotalReads-s.DiscardedReads) / float64(s.TotalReads)
		html += fmt.Sprintf("        <tr><td class='metric'>Kept Reads</td><td class='value'>%d <span class='percentage'>(%.2f%%)</span></td></tr>\n",
			s.TotalReads-s.DiscardedReads, keptPct)
	}

	html += fmt.Sprintf("        <tr><td class='metric'>Total Bases</td><td class='value'>%d</td></tr>\n", s.TotalBases)

	if s.TotalBases > 0 {
		trimmedBasesPct := 100.0 * float64(s.TrimmedBases) / float64(s.TotalBases)
		html += fmt.Sprintf("        <tr><td class='metric'>Trimmed Bases</td><td class='value'>%d <span class='percentage'>(%.2f%%)</span></td></tr>\n",
			s.TrimmedBases, trimmedBasesPct)
	}

	html += fmt.Sprintf("        <tr><td class='metric'>Processing Time</td><td class='value'>%v</td></tr>\n", s.ProcessingTime)

	if s.DetectedAdapter3 != "" {
		html += fmt.Sprintf("        <tr><td class='metric'>Detected 3' Adapter</td><td class='value'>%s</td></tr>\n", s.DetectedAdapter3)
	}
	if s.DetectedAdapter5 != "" {
		html += fmt.Sprintf("        <tr><td class='metric'>Detected 5' Adapter</td><td class='value'>%s</td></tr>\n", s.DetectedAdapter5)
	}

	if s.UMIStats != nil {
		html += fmt.Sprintf("        <tr><td class='metric'>Total UMIs</td><td class='value'>%d</td></tr>\n", s.UMIStats.TotalUMIs)
		html += fmt.Sprintf("        <tr><td class='metric'>Unique UMIs</td><td class='value'>%d</td></tr>\n", s.UMIStats.UniqueUMIs)
	}

	html += `    </table>
</body>
</html>
`
	return html
}

// TrimSingleEnd trims adapters from single-end FASTQ reads.
func TrimSingleEnd(input io.Reader, output io.Writer, encoding fastq.QualityEncoding, opts TrimOptions) (*TrimStats, error) {
	startTime := time.Now()

	reader := fastq.NewReader(input, encoding)
	writer := fastq.NewWriter(output, encoding)

	stats := &TrimStats{}

	// buffered holds records that were pre-read for auto-detection and still
	// need to be processed; it is drained before reading further from reader.
	var buffered []*fastq.Record

	// Auto-detect adapters if requested
	if opts.AutoDetect && opts.Adapter3 == "" && opts.Adapter5 == "" {
		maxReads := opts.AutoDetectReads
		if maxReads <= 0 {
			maxReads = 1000
		}
		for i := 0; i < maxReads; i++ {
			record, err := reader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return stats, fmt.Errorf("error reading FASTQ: %w", err)
			}
			buffered = append(buffered, record)
		}
		detected := detectAdaptersFromReads(buffered)
		if detected.Adapter3 != "" {
			opts.Adapter3 = detected.Adapter3
			stats.DetectedAdapter3 = detected.Adapter3
		}
		if detected.Adapter5 != "" {
			opts.Adapter5 = detected.Adapter5
			stats.DetectedAdapter5 = detected.Adapter5
		}
	}

	// Initialize UMI stats if UMI extraction is enabled
	if opts.UMILength > 0 {
		stats.UMIStats = &UMIStats{
			UMIDistribution: make(map[string]int),
		}
	}

	readCount := 0
	for {
		var record *fastq.Record
		if len(buffered) > 0 {
			record = buffered[0]
			buffered = buffered[1:]
		} else {
			var err error
			record, err = reader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return stats, fmt.Errorf("error reading FASTQ: %w", err)
			}
		}

		stats.mu.Lock()
		stats.TotalReads++
		readCount++
		stats.mu.Unlock()

		originalLength := len(record.Sequence)
		stats.mu.Lock()
		stats.TotalBases += int64(originalLength)
		stats.mu.Unlock()

		// Extract UMI if enabled
		var umi string
		if opts.UMILength > 0 {
			record, umi = extractUMI(record, opts)
			if umi != "" {
				stats.mu.Lock()
				stats.UMIStats.TotalUMIs++
				stats.UMIStats.UMIDistribution[umi]++
				stats.mu.Unlock()
			}
		}

		// Trim adapters
		trimmed := trimRecord(record, opts, stats)

		// Add UMI to description if extracted
		if umi != "" {
			trimmed.Description = trimmed.Description + " UMI:" + umi
		}

		// Check if read passes length threshold
		if len(trimmed.Sequence) >= opts.MinLength {
			if err := writer.Write(trimmed); err != nil {
				return stats, fmt.Errorf("error writing FASTQ: %w", err)
			}
			if len(trimmed.Sequence) < originalLength {
				stats.mu.Lock()
				stats.TrimmedReads++
				stats.TrimmedBases += int64(originalLength - len(trimmed.Sequence))
				stats.mu.Unlock()
			}
		} else {
			stats.mu.Lock()
			stats.DiscardedReads++
			stats.mu.Unlock()
		}

		// Progress reporting
		if opts.ProgressReport && readCount%opts.ProgressInterval == 0 {
			fmt.Fprintf(io.Discard, "\rProcessed %d reads...", readCount)
		}
	}

	// Flush writer
	if err := writer.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output: %w", err)
	}

	// Calculate unique UMIs
	if stats.UMIStats != nil {
		stats.UMIStats.UniqueUMIs = len(stats.UMIStats.UMIDistribution)
	}

	stats.ProcessingTime = time.Since(startTime)
	return stats, nil
}

// TrimPairedEnd trims adapters from paired-end FASTQ reads.
func TrimPairedEnd(input1, input2 io.Reader, output1, output2, outputSingle io.Writer,
	encoding fastq.QualityEncoding, opts TrimOptions) (*TrimStats, error) {

	startTime := time.Now()

	reader1 := fastq.NewReader(input1, encoding)
	reader2 := fastq.NewReader(input2, encoding)
	writer1 := fastq.NewWriter(output1, encoding)
	writer2 := fastq.NewWriter(output2, encoding)

	var writerSingle *fastq.Writer
	if outputSingle != nil {
		writerSingle = fastq.NewWriter(outputSingle, encoding)
	}

	stats := &TrimStats{}

	// Initialize UMI stats if UMI extraction is enabled
	if opts.UMILength > 0 {
		stats.UMIStats = &UMIStats{
			UMIDistribution: make(map[string]int),
		}
	}

	readCount := 0
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

		stats.mu.Lock()
		stats.TotalReads += 2
		readCount += 2
		stats.mu.Unlock()

		originalLen1 := len(record1.Sequence)
		originalLen2 := len(record2.Sequence)
		stats.mu.Lock()
		stats.TotalBases += int64(originalLen1 + originalLen2)
		stats.mu.Unlock()

		// Extract UMI if enabled (from first read only)
		var umi string
		if opts.UMILength > 0 {
			record1, umi = extractUMI(record1, opts)
			if umi != "" {
				stats.mu.Lock()
				stats.UMIStats.TotalUMIs++
				stats.UMIStats.UMIDistribution[umi]++
				stats.mu.Unlock()
			}
		}

		// Trim both reads
		trimmed1 := trimRecord(record1, opts, stats)
		trimmed2 := trimRecord(record2, opts, stats)

		// Add UMI to descriptions if extracted
		if umi != "" {
			trimmed1.Description = trimmed1.Description + " UMI:" + umi
			trimmed2.Description = trimmed2.Description + " UMI:" + umi
		}

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
				stats.mu.Lock()
				stats.TrimmedReads++
				stats.TrimmedBases += int64((originalLen1 - len(trimmed1.Sequence)) +
					(originalLen2 - len(trimmed2.Sequence)))
				stats.mu.Unlock()
			}
		} else if writerSingle != nil {
			// One or both fail - write survivors to single output
			if pass1 {
				if err := writerSingle.Write(trimmed1); err != nil {
					return stats, fmt.Errorf("error writing single output: %w", err)
				}
			} else {
				stats.mu.Lock()
				stats.DiscardedReads++
				stats.mu.Unlock()
			}

			if pass2 {
				if err := writerSingle.Write(trimmed2); err != nil {
					return stats, fmt.Errorf("error writing single output: %w", err)
				}
			} else {
				stats.mu.Lock()
				stats.DiscardedReads++
				stats.mu.Unlock()
			}
		} else {
			// Both discarded
			stats.mu.Lock()
			stats.DiscardedReads += 2
			stats.mu.Unlock()
		}

		// Progress reporting
		if opts.ProgressReport && readCount%opts.ProgressInterval == 0 {
			fmt.Fprintf(io.Discard, "\rProcessed %d reads...", readCount)
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

	// Calculate unique UMIs
	if stats.UMIStats != nil {
		stats.UMIStats.UniqueUMIs = len(stats.UMIStats.UMIDistribution)
	}

	stats.ProcessingTime = time.Since(startTime)
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
// Returns -1 if adapter not found. Uses improved alignment scoring.
func findAdapter(seq string, adapter string, minOverlap int, errorRate float64) int {
	// Use improved algorithm for better accuracy
	return improvedFindAdapter(seq, adapter, minOverlap, errorRate)
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

// Common adapter sequences for auto-detection
var commonAdapters = []string{
	"AGATCGGAAGAGC",                      // Illumina TruSeq Universal
	"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", // Illumina TruSeq Full
	"CTGTCTCTTATACACATCT",                // Illumina Nextera
	"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",  // Illumina Small RNA 3'
	"TGGAATTCTCGGGTGCCAAGG",              // Illumina Small RNA 5'
	"CGCCTTGGCCGTACAGCAG",                // SOLiD
	"ATCTCGTATGCCGTCTTCTGCTTG",           // Ion Torrent
}

// detectAdapters attempts to auto-detect adapter sequences by reading up to
// maxReads records from reader.
func detectAdapters(reader *fastq.Reader, maxReads int) (TrimOptions, error) {
	var reads []*fastq.Record
	for i := 0; i < maxReads; i++ {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return TrimOptions{}, err
		}
		reads = append(reads, record)
	}
	return detectAdaptersFromReads(reads), nil
}

// detectAdaptersFromReads inspects a slice of already-read records and returns
// TrimOptions populated with any common adapters that appear in at least 10% of
// them.
func detectAdaptersFromReads(reads []*fastq.Record) TrimOptions {
	opts := TrimOptions{}

	if len(reads) == 0 {
		return opts
	}

	// Count adapter occurrences
	adapter3Counts := make(map[string]int)
	adapter5Counts := make(map[string]int)

	for _, record := range reads {
		seq := string(record.Sequence)

		// Check 3' adapters (at end of read)
		for _, adapter := range commonAdapters {
			if pos := findAdapter(seq, adapter, 5, 0.15); pos >= 0 {
				adapter3Counts[adapter]++
			}
		}

		// Check 5' adapters (at beginning of read)
		for _, adapter := range commonAdapters {
			if pos := findAdapter(seq, adapter, 5, 0.15); pos >= 0 && pos < 10 {
				adapter5Counts[adapter]++
			}
		}
	}

	// Find most common adapters (must appear in at least 10% of reads).
	// Iterate commonAdapters in declaration order so the result is
	// deterministic when several candidates tie on count.
	threshold := len(reads) / 10

	var maxCount3 int
	for _, adapter := range commonAdapters {
		if count := adapter3Counts[adapter]; count > maxCount3 && count >= threshold {
			maxCount3 = count
			opts.Adapter3 = adapter
		}
	}

	var maxCount5 int
	for _, adapter := range commonAdapters {
		if count := adapter5Counts[adapter]; count > maxCount5 && count >= threshold {
			maxCount5 = count
			opts.Adapter5 = adapter
		}
	}

	return opts
}

// extractUMI extracts UMI/barcode from read.
func extractUMI(record *fastq.Record, opts TrimOptions) (*fastq.Record, string) {
	if opts.UMILength <= 0 || opts.UMILength >= len(record.Sequence) {
		return record, ""
	}

	var umi string
	var newSeq, newQual []byte

	if opts.UMIPosition == "5prime" {
		// Extract from 5' end
		umi = string(record.Sequence[:opts.UMILength])
		newSeq = record.Sequence[opts.UMILength:]
		newQual = record.Quality[opts.UMILength:]
	} else {
		// Extract from 3' end
		umi = string(record.Sequence[len(record.Sequence)-opts.UMILength:])
		newSeq = record.Sequence[:len(record.Sequence)-opts.UMILength]
		newQual = record.Quality[:len(record.Quality)-opts.UMILength]
	}

	return &fastq.Record{
		ID:          record.ID,
		Description: record.Description,
		Sequence:    newSeq,
		Quality:     newQual,
	}, umi
}

// improvedFindAdapter uses improved Smith-Waterman-like local alignment for adapter matching.
func improvedFindAdapter(seq string, adapter string, minOverlap int, errorRate float64) int {
	if len(adapter) == 0 || len(seq) < minOverlap {
		return -1
	}

	maxErrors := int(float64(len(adapter)) * errorRate)
	bestPos := -1
	bestScore := -1

	// Scan through sequence
	for i := 0; i <= len(seq)-minOverlap; i++ {
		compareLen := min(len(adapter), len(seq)-i)

		if compareLen < minOverlap {
			continue
		}

		// Calculate alignment score
		matches := 0
		for j := 0; j < compareLen; j++ {
			if seq[i+j] == adapter[j] {
				matches++
			}
		}

		errors := compareLen - matches

		// Score: favor longer alignments with fewer errors
		score := matches*2 - errors

		if errors <= maxErrors && score > bestScore {
			bestScore = score
			bestPos = i
		}
	}

	return bestPos
}
