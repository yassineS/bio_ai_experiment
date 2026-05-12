// Package seqtk provides core functionality for sequence processing.
// This is a Go reimplementation of seqtk, a fast FASTA/Q processor.
package seqtk

import (
	"compress/bzip2"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"

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
	if totalBases > 0 {
		stats.GCContent = float64(totalGC) / float64(totalBases) * 100
	}

	return stats, nil
}

// CalculateFastaStatsParallel calculates statistics for a FASTA file using parallel processing.
func CalculateFastaStatsParallel(r io.Reader, workers int) (*Stats, error) {
	reader := fasta.NewReader(r)
	records, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	if len(records) == 0 {
		return &Stats{}, nil
	}

	// For small files, use sequential processing
	if len(records) < 100 || workers <= 1 {
		return CalculateFastaStats(r)
	}

	// Split work among workers
	chunkSize := (len(records) + workers - 1) / workers
	type result struct {
		totalBases int64
		totalGC    int64
		minLen     int
		maxLen     int
	}

	resultChan := make(chan result, workers)

	for i := 0; i < workers; i++ {
		start := i * chunkSize
		end := start + chunkSize
		if end > len(records) {
			end = len(records)
		}
		if start >= len(records) {
			break
		}

		go func(chunk []*fasta.Record) {
			var r result
			r.minLen = chunk[0].Length()
			r.maxLen = chunk[0].Length()

			for _, record := range chunk {
				length := record.Length()
				r.totalBases += int64(length)

				if length < r.minLen {
					r.minLen = length
				}
				if length > r.maxLen {
					r.maxLen = length
				}

				for _, b := range record.Sequence {
					if b == 'G' || b == 'C' || b == 'g' || b == 'c' {
						r.totalGC++
					}
				}
			}
			resultChan <- r
		}(records[start:end])
	}

	// Collect results
	stats := &Stats{
		NumSequences: len(records),
		MinLength:    records[0].Length(),
		MaxLength:    records[0].Length(),
	}

	var totalBases int64
	var totalGC int64
	activeWorkers := workers
	if len(records) < workers*chunkSize {
		activeWorkers = (len(records) + chunkSize - 1) / chunkSize
	}

	for i := 0; i < activeWorkers; i++ {
		r := <-resultChan
		totalBases += r.totalBases
		totalGC += r.totalGC
		if r.minLen < stats.MinLength {
			stats.MinLength = r.minLen
		}
		if r.maxLen > stats.MaxLength {
			stats.MaxLength = r.maxLen
		}
	}

	stats.TotalBases = totalBases
	stats.AvgLength = float64(totalBases) / float64(len(records))
	if totalBases > 0 {
		stats.GCContent = float64(totalGC) / float64(totalBases) * 100
	}

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

// FilterOptions contains options for sequence filtering.
type FilterOptions struct {
	MinLength int    // Minimum sequence length (0 = no filter)
	MaxLength int    // Maximum sequence length (0 = no filter)
	Pattern   string // Pattern to match in sequence ID (empty = no filter)
}

// Filter sequences based on filter options.
func Filter(input io.Reader, output io.Writer, opts FilterOptions, isFastq bool, encoding fastq.QualityEncoding) error {
	if isFastq {
		return filterFastq(input, output, opts, encoding)
	}
	return filterFasta(input, output, opts)
}

func filterFasta(input io.Reader, output io.Writer, opts FilterOptions) error {
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

		// Apply filters
		if !passesFilter(record.ID, record.Length(), opts) {
			continue
		}

		if err := writer.Write(record); err != nil {
			return err
		}
	}

	return writer.Flush()
}

func filterFastq(input io.Reader, output io.Writer, opts FilterOptions, encoding fastq.QualityEncoding) error {
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

		// Apply filters
		if !passesFilter(record.ID, record.Length(), opts) {
			continue
		}

		if err := writer.Write(record); err != nil {
			return err
		}
	}

	return writer.Flush()
}

func passesFilter(id string, length int, opts FilterOptions) bool {
	// Check length filters
	if opts.MinLength > 0 && length < opts.MinLength {
		return false
	}
	if opts.MaxLength > 0 && length > opts.MaxLength {
		return false
	}

	// Check pattern filter
	if opts.Pattern != "" && !strings.Contains(id, opts.Pattern) {
		return false
	}

	return true
}

// Subseq extracts a subsequence from each sequence.
func Subseq(input io.Reader, output io.Writer, start, end int, isFastq bool, encoding fastq.QualityEncoding) error {
	if isFastq {
		return subseqFastq(input, output, start, end, encoding)
	}
	return subseqFasta(input, output, start, end)
}

func subseqFasta(input io.Reader, output io.Writer, start, end int) error {
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

		// Extract subsequence (1-based indexing, inclusive)
		length := record.Length()

		// Adjust negative indices (from end)
		if end < 0 {
			end = length + end + 1
		}
		if start < 0 {
			start = length + start + 1
		}

		// Convert to 0-based indexing
		startIdx := start - 1
		endIdx := end

		// Bounds checking
		if startIdx < 0 {
			startIdx = 0
		}
		if endIdx > length {
			endIdx = length
		}
		if startIdx >= length || startIdx >= endIdx {
			continue // Skip this sequence
		}

		// Create new record with subsequence
		subRecord := &fasta.Record{
			ID:          record.ID,
			Description: record.Description,
			Sequence:    record.Sequence[startIdx:endIdx],
		}

		if err := writer.Write(subRecord); err != nil {
			return err
		}
	}

	return writer.Flush()
}

func subseqFastq(input io.Reader, output io.Writer, start, end int, encoding fastq.QualityEncoding) error {
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

		// Extract subsequence (1-based indexing, inclusive)
		length := record.Length()

		// Adjust negative indices (from end)
		if end < 0 {
			end = length + end + 1
		}
		if start < 0 {
			start = length + start + 1
		}

		// Convert to 0-based indexing
		startIdx := start - 1
		endIdx := end

		// Bounds checking
		if startIdx < 0 {
			startIdx = 0
		}
		if endIdx > length {
			endIdx = length
		}
		if startIdx >= length || startIdx >= endIdx {
			continue // Skip this sequence
		}

		// Create new record with subsequence
		subRecord := &fastq.Record{
			ID:          record.ID,
			Description: record.Description,
			Sequence:    record.Sequence[startIdx:endIdx],
			Quality:     record.Quality[startIdx:endIdx],
		}

		if err := writer.Write(subRecord); err != nil {
			return err
		}
	}

	return writer.Flush()
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
// Handles compressed files (.gz, .bz2) automatically.
func GetFileType(filename string) (bool, error) {
	file, err := os.Open(filename)
	if err != nil {
		return false, err
	}
	defer file.Close()

	reader, err := DecompressReader(file, filename)
	if err != nil {
		return false, err
	}

	buf := make([]byte, 1)
	_, err = reader.Read(buf)
	if err != nil {
		return false, err
	}

	// FASTQ starts with '@', FASTA with '>'
	return buf[0] == '@', nil
}

// DecompressReader wraps a reader with decompression based on file extension.
// Supports .gz (gzip) and .bz2 (bzip2) compression.
func DecompressReader(r io.Reader, filename string) (io.Reader, error) {
	ext := strings.ToLower(filepath.Ext(filename))

	switch ext {
	case ".gz":
		return gzip.NewReader(r)
	case ".bz2":
		return bzip2.NewReader(r), nil
	default:
		return r, nil
	}
}

// CompressWriter wraps a writer with compression based on file extension.
// Supports .gz (gzip) compression.
// Returns nil if no compression is needed (caller should use original writer).
func CompressWriter(w io.Writer, filename string) (io.WriteCloser, error) {
	ext := strings.ToLower(filepath.Ext(filename))

	switch ext {
	case ".gz":
		return gzip.NewWriter(w), nil
	case ".bz2":
		return nil, fmt.Errorf("bzip2 compression not supported for writing")
	default:
		return nil, nil
	}
}

// nopCloser wraps a Writer to provide a no-op Close method
type nopCloser struct {
	io.Writer
}

func (nopCloser) Close() error { return nil }

// OpenInput opens a file for reading with automatic decompression support.
// If filename is "-", reads from stdin.
func OpenInput(filename string) (io.ReadCloser, error) {
	if filename == "-" {
		// For stdin, we can't decompress based on filename, so try to detect
		return io.NopCloser(os.Stdin), nil
	}

	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	reader, err := DecompressReader(file, filename)
	if err != nil {
		file.Close()
		return nil, err
	}

	// If reader is not the file itself, wrap in a composite closer
	if reader != file {
		return &compositeCloser{reader: reader, file: file}, nil
	}

	return file, nil
}

// compositeCloser closes both the reader and underlying file
type compositeCloser struct {
	reader io.Reader
	file   *os.File
}

func (c *compositeCloser) Read(p []byte) (n int, err error) {
	return c.reader.Read(p)
}

func (c *compositeCloser) Close() error {
	// Close file first (reader might need it)
	return c.file.Close()
}

// OpenOutput opens a file for writing with automatic compression support.
// If filename is "-" or empty, writes to stdout.
func OpenOutput(filename string) (io.WriteCloser, error) {
	if filename == "-" || filename == "" {
		return &nopCloser{os.Stdout}, nil
	}

	file, err := os.Create(filename)
	if err != nil {
		return nil, err
	}

	writer, err := CompressWriter(file, filename)
	if err != nil {
		file.Close()
		return nil, err
	}

	// If writer is nil (no compression), just use the file
	if writer == nil {
		return file, nil
	}

	// Otherwise wrap in a composite closer
	return &compositeWriter{writer: writer, file: file}, nil
}

// compositeWriter closes both the writer and underlying file
type compositeWriter struct {
	writer io.WriteCloser
	file   *os.File
}

func (c *compositeWriter) Write(p []byte) (n int, err error) {
	return c.writer.Write(p)
}

func (c *compositeWriter) Close() error {
	// Close writer first (to flush), then file
	if err := c.writer.Close(); err != nil {
		c.file.Close()
		return err
	}
	return c.file.Close()
}
