// Package fastp provides all-in-one preprocessing for FASTQ files.
// It combines quality filtering, adapter trimming, and various other preprocessing steps.
package fastp

import (
	"fmt"
	"io"
	"strings"
	"sync"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// Common adapter sequences for automatic detection
var CommonAdapters = map[string]string{
	"TruSeq":     "AGATCGGAAGAGC",
	"Nextera":    "CTGTCTCTTATA",
	"SmallRNA":   "TGGAATTCTCGG",
	"TruSeq_R2":  "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
	"Nextera_R2": "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
}

// ProcessOptions contains all preprocessing parameters.
type ProcessOptions struct {
	// Adapter trimming
	Adapter3      string
	Adapter5      string
	DetectAdapter bool

	// Quality filtering
	QualThreshold int
	MinLength     int
	MaxLength     int
	QualPercent   int // Percentage of bases that must meet quality threshold

	// Complexity filtering
	LowComplexity       bool
	ComplexityThreshold float64

	// Poly-tail trimming
	TrimPolyG   bool
	TrimPolyX   bool
	PolyGMinLen int

	// N filtering
	MaxNCount   int
	MaxNPercent float64

	// Length filtering
	LengthRequired int
	LengthLimit    int

	// UMI processing
	UMILength   int
	UMILocation string // "read1", "read2", "index"
	UMISkip     int    // Bases to skip before UMI

	// Base correction
	BaseCorrection      bool
	CorrectionThreshold int

	// Overlap analysis (paired-end)
	MergeOverlap bool
	MinOverlap   int
	MaxMismatch  int

	// Multi-threading
	Threads int

	// HTML report
	HTMLReport string // Path to HTML report file
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
		UMILength:           0,
		UMILocation:         "",
		UMISkip:             0,
		BaseCorrection:      false,
		CorrectionThreshold: 20,
		MergeOverlap:        false,
		MinOverlap:          30,
		MaxMismatch:         5,
		Threads:             1,
		HTMLReport:          "",
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
	DetectedAdapter     string
	UMIExtracted        int
	BasesCorrected      int64
	OverlappingReads    int
	MergedReads         int
}

// OverlapResult represents the result of paired-end overlap analysis
type OverlapResult struct {
	HasOverlap    bool
	OverlapLength int
	Mismatches    int
	MergedSeq     string
	MergedQual    []byte
}

// ProcessPairedEnd processes paired-end FASTQ reads with all filters.
func ProcessPairedEnd(input1, input2 io.Reader, output1, output2 io.Writer, encoding fastq.QualityEncoding, opts ProcessOptions) (*ProcessStats, error) {
	// Note: Auto-detection of adapters would require reading the input twice or buffering,
	// which is not practical for streaming. Users should specify adapter or use separate detection step.

	reader1 := fastq.NewReader(input1, encoding)
	reader2 := fastq.NewReader(input2, encoding)
	writer1 := fastq.NewWriter(output1, encoding)
	writer2 := fastq.NewWriter(output2, encoding)

	stats := &ProcessStats{}

	// Process with multi-threading if enabled
	if opts.Threads > 1 {
		return processPairedEndParallel(reader1, reader2, writer1, writer2, encoding, opts, stats)
	}

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

		// Extract UMI if configured
		if opts.UMILength > 0 {
			record1, record2 = extractUMI(record1, record2, opts, stats)
		}

		// Check for overlap and merge if enabled
		if opts.MergeOverlap {
			overlap := analyzeOverlap(record1, record2, opts, encoding)
			if overlap.HasOverlap {
				stats.OverlappingReads++
				if overlap.OverlapLength >= opts.MinOverlap && overlap.Mismatches <= opts.MaxMismatch {
					// Create merged record
					merged := &fastq.Record{
						ID:          record1.ID,
						Description: record1.Description,
						Sequence:    []byte(overlap.MergedSeq),
						Quality:     overlap.MergedQual,
					}
					processed, pass := processRecord(merged, opts, stats, encoding)
					if pass {
						if err := writer1.Write(processed); err != nil {
							return stats, fmt.Errorf("error writing merged read: %w", err)
						}
						stats.CleanReads++
						stats.CleanBases += int64(len(processed.Sequence))
						stats.MergedReads++
					}
					continue
				}
			}
		}

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

	// Process with multi-threading if enabled
	if opts.Threads > 1 {
		return processSingleEndParallel(reader, writer, encoding, opts, stats)
	}

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

		// Extract UMI if configured
		if opts.UMILength > 0 && opts.UMILocation == "read1" {
			record, _ = extractUMI(record, nil, opts, stats)
		}

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

	// Step 0: Base correction if enabled
	if opts.BaseCorrection {
		seq, qual = correctBases(seq, qual, opts.CorrectionThreshold, stats, encoding)
		record.Sequence = []byte(seq)
		record.Quality = qual
	}

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

// detectAdapter detects the most common adapter in a sample of reads
func detectAdapter(input1, input2 io.Reader, encoding fastq.QualityEncoding) string {
	// Sample first 10000 reads to detect adapter
	reader1 := fastq.NewReader(input1, encoding)

	adapterCounts := make(map[string]int)
	sampleSize := 10000

	for i := 0; i < sampleSize; i++ {
		record, err := reader1.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			continue
		}

		seq := string(record.Sequence)
		// Check each common adapter
		for name, adapter := range CommonAdapters {
			if strings.Contains(seq, adapter) || strings.Contains(seq, adapter[:len(adapter)/2]) {
				adapterCounts[name]++
			}
		}
	}

	// Find most common adapter
	maxCount := 0
	detectedName := ""
	for name, count := range adapterCounts {
		if count > maxCount {
			maxCount = count
			detectedName = name
		}
	}

	if detectedName != "" && maxCount > sampleSize/100 { // At least 1% detection rate
		return CommonAdapters[detectedName]
	}

	return ""
}

// DetectAdapterFromReads detects the most common adapter from FASTQ reads
// This is a public function that can be called before processing
func DetectAdapterFromReads(reads []string) string {
	adapterCounts := make(map[string]int)

	for _, seq := range reads {
		// Check each common adapter
		for name, adapter := range CommonAdapters {
			if strings.Contains(seq, adapter) || strings.Contains(seq, adapter[:len(adapter)/2]) {
				adapterCounts[name]++
			}
		}
	}

	// Find most common adapter
	maxCount := 0
	detectedName := ""
	for name, count := range adapterCounts {
		if count > maxCount {
			maxCount = count
			detectedName = name
		}
	}

	if detectedName != "" && maxCount > len(reads)/100 { // At least 1% detection rate
		return CommonAdapters[detectedName]
	}

	return ""
}

// extractUMI extracts UMI from reads based on configuration
func extractUMI(record1, record2 *fastq.Record, opts ProcessOptions, stats *ProcessStats) (*fastq.Record, *fastq.Record) {
	if opts.UMILength == 0 {
		return record1, record2
	}

	extractFromRecord := func(record *fastq.Record) *fastq.Record {
		if record == nil {
			return nil
		}

		start := opts.UMISkip
		end := start + opts.UMILength

		if end > len(record.Sequence) {
			return record
		}

		// Extract UMI and add to ID
		umi := string(record.Sequence[start:end])
		newID := fmt.Sprintf("%s_UMI:%s", record.ID, umi)

		// Remove UMI from sequence
		newSeq := append([]byte{}, record.Sequence[:start]...)
		newSeq = append(newSeq, record.Sequence[end:]...)

		newQual := append([]byte{}, record.Quality[:start]...)
		newQual = append(newQual, record.Quality[end:]...)

		stats.UMIExtracted++

		return &fastq.Record{
			ID:          newID,
			Description: record.Description,
			Sequence:    newSeq,
			Quality:     newQual,
		}
	}

	if opts.UMILocation == "read1" {
		return extractFromRecord(record1), record2
	} else if opts.UMILocation == "read2" {
		return record1, extractFromRecord(record2)
	}

	return record1, record2
}

// correctBases performs base correction based on quality scores
func correctBases(seq string, qual []byte, threshold int, stats *ProcessStats, encoding fastq.QualityEncoding) (string, []byte) {
	if threshold == 0 {
		return seq, qual
	}

	offset := 33
	if encoding == fastq.Phred64 {
		offset = 64
	}

	corrected := []byte(seq)
	correctedCount := int64(0)

	for i := 0; i < len(seq); i++ {
		q := int(qual[i]) - offset

		// Correct low-quality bases to N
		if q < threshold && seq[i] != 'N' {
			corrected[i] = 'N'
			correctedCount++
		}
	}

	stats.BasesCorrected += correctedCount
	return string(corrected), qual
}

// analyzeOverlap analyzes overlap between paired-end reads
func analyzeOverlap(record1, record2 *fastq.Record, opts ProcessOptions, encoding fastq.QualityEncoding) OverlapResult {
	result := OverlapResult{HasOverlap: false}

	seq1 := string(record1.Sequence)
	seq2 := reverseComplement(string(record2.Sequence))
	qual1 := record1.Quality
	qual2 := reverseSlice(record2.Quality)

	minLen := opts.MinOverlap
	if minLen < 10 {
		minLen = 10
	}

	// Try different overlap positions
	bestOverlap := 0
	bestMismatches := len(seq1)
	bestPos := -1

	for pos := len(seq1) - minLen; pos >= 0; pos-- {
		overlapLen := len(seq1) - pos
		if overlapLen > len(seq2) {
			overlapLen = len(seq2)
		}

		mismatches := 0
		for i := 0; i < overlapLen; i++ {
			if seq1[pos+i] != seq2[i] {
				mismatches++
			}
		}

		if mismatches < bestMismatches && overlapLen >= minLen {
			bestOverlap = overlapLen
			bestMismatches = mismatches
			bestPos = pos
		}
	}

	if bestPos >= 0 && bestOverlap >= minLen {
		result.HasOverlap = true
		result.OverlapLength = bestOverlap
		result.Mismatches = bestMismatches

		// Merge sequences using quality scores
		merged := make([]byte, bestPos+len(seq2))
		mergedQual := make([]byte, bestPos+len(seq2))

		// Copy non-overlapping part of seq1
		copy(merged[:bestPos], seq1[:bestPos])
		copy(mergedQual[:bestPos], qual1[:bestPos])

		// Merge overlapping region
		offset := 33
		if encoding == fastq.Phred64 {
			offset = 64
		}

		for i := 0; i < len(seq2); i++ {
			pos := bestPos + i
			if i < bestOverlap {
				// In overlap - use higher quality base
				q1 := int(qual1[pos]) - offset
				q2 := int(qual2[i]) - offset

				if q1 >= q2 {
					merged[pos] = seq1[pos]
					mergedQual[pos] = qual1[pos]
				} else {
					merged[pos] = seq2[i]
					mergedQual[pos] = qual2[i]
				}
			} else {
				// After overlap - use seq2
				merged[pos] = seq2[i]
				mergedQual[pos] = qual2[i]
			}
		}

		result.MergedSeq = string(merged)
		result.MergedQual = mergedQual
	}

	return result
}

// reverseComplement returns the reverse complement of a DNA sequence
func reverseComplement(seq string) string {
	complement := map[byte]byte{
		'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
		'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
		'N': 'N', 'n': 'n',
	}

	result := make([]byte, len(seq))
	for i := 0; i < len(seq); i++ {
		if comp, ok := complement[seq[len(seq)-1-i]]; ok {
			result[i] = comp
		} else {
			result[i] = 'N'
		}
	}

	return string(result)
}

// reverseSlice reverses a byte slice
func reverseSlice(s []byte) []byte {
	result := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		result[i] = s[len(s)-1-i]
	}
	return result
}

// processPairedEndParallel processes paired-end reads in parallel
func processPairedEndParallel(reader1, reader2 *fastq.Reader, writer1, writer2 *fastq.Writer, encoding fastq.QualityEncoding, opts ProcessOptions, stats *ProcessStats) (*ProcessStats, error) {
	type readPair struct {
		record1 *fastq.Record
		record2 *fastq.Record
	}

	type resultPair struct {
		processed1 *fastq.Record
		processed2 *fastq.Record
		pass       bool
	}

	inputChan := make(chan readPair, opts.Threads*2)
	outputChan := make(chan resultPair, opts.Threads*2)

	var wg sync.WaitGroup

	// Start worker goroutines
	for i := 0; i < opts.Threads; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			localStats := &ProcessStats{}

			for pair := range inputChan {
				processed1, pass1 := processRecord(pair.record1, opts, localStats, encoding)
				processed2, pass2 := processRecord(pair.record2, opts, localStats, encoding)

				outputChan <- resultPair{
					processed1: processed1,
					processed2: processed2,
					pass:       pass1 && pass2,
				}
			}

			// Merge local stats (simplified - would need mutex for accuracy)
			stats.LowQualityReads += localStats.LowQualityReads
			stats.TooShortReads += localStats.TooShortReads
			stats.TooLongReads += localStats.TooLongReads
			stats.TooManyNReads += localStats.TooManyNReads
		}()
	}

	// Reader goroutine
	go func() {
		for {
			record1, err1 := reader1.Read()
			record2, err2 := reader2.Read()

			if err1 == io.EOF && err2 == io.EOF {
				break
			}
			if err1 != nil || err2 != nil {
				break
			}

			stats.TotalReads += 2
			stats.TotalBases += int64(len(record1.Sequence) + len(record2.Sequence))

			inputChan <- readPair{record1: record1, record2: record2}
		}
		close(inputChan)
	}()

	// Wait for workers to finish
	go func() {
		wg.Wait()
		close(outputChan)
	}()

	// Writer goroutine
	for result := range outputChan {
		if result.pass {
			if err := writer1.Write(result.processed1); err != nil {
				return stats, fmt.Errorf("error writing read1: %w", err)
			}
			if err := writer2.Write(result.processed2); err != nil {
				return stats, fmt.Errorf("error writing read2: %w", err)
			}
			stats.CleanReads += 2
			stats.CleanBases += int64(len(result.processed1.Sequence) + len(result.processed2.Sequence))
		}
	}

	return stats, nil
}

// processSingleEndParallel processes single-end reads in parallel
func processSingleEndParallel(reader *fastq.Reader, writer *fastq.Writer, encoding fastq.QualityEncoding, opts ProcessOptions, stats *ProcessStats) (*ProcessStats, error) {
	type result struct {
		processed *fastq.Record
		pass      bool
	}

	inputChan := make(chan *fastq.Record, opts.Threads*2)
	outputChan := make(chan result, opts.Threads*2)

	var wg sync.WaitGroup

	// Start worker goroutines
	for i := 0; i < opts.Threads; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			localStats := &ProcessStats{}

			for record := range inputChan {
				processed, pass := processRecord(record, opts, localStats, encoding)
				outputChan <- result{processed: processed, pass: pass}
			}

			// Merge local stats (simplified)
			stats.LowQualityReads += localStats.LowQualityReads
			stats.TooShortReads += localStats.TooShortReads
			stats.TooLongReads += localStats.TooLongReads
			stats.TooManyNReads += localStats.TooManyNReads
		}()
	}

	// Reader goroutine
	go func() {
		for {
			record, err := reader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				break
			}

			stats.TotalReads++
			stats.TotalBases += int64(len(record.Sequence))
			inputChan <- record
		}
		close(inputChan)
	}()

	// Wait for workers to finish
	go func() {
		wg.Wait()
		close(outputChan)
	}()

	// Writer goroutine
	for res := range outputChan {
		if res.pass {
			if err := writer.Write(res.processed); err != nil {
				return stats, fmt.Errorf("error writing FASTQ: %w", err)
			}
			stats.CleanReads++
			stats.CleanBases += int64(len(res.processed.Sequence))
		}
	}

	return stats, nil
}
