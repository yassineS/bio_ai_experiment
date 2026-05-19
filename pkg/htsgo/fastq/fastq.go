// Package fastq provides utilities for reading and writing FASTQ format files.
// FASTQ is a text-based format for storing biological sequences and their quality scores.
//
// Format specification:
//   - Line 1: '@' followed by sequence identifier and optional description
//   - Line 2: Raw sequence
//   - Line 3: '+' optionally followed by sequence identifier and description
//   - Line 4: Quality scores (ASCII-encoded, typically Phred+33 or Phred+64)
package fastq

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strings"
)

// QualityEncoding represents the quality score encoding format.
type QualityEncoding int

const (
	// Phred33 (Sanger, Illumina 1.8+) - quality scores from 0-93 using ASCII 33-126
	Phred33 QualityEncoding = iota
	// Phred64 (Illumina 1.3-1.7) - quality scores from 0-62 using ASCII 64-126
	Phred64
)

// Record represents a single FASTQ sequence record.
type Record struct {
	ID          string // Sequence identifier (first word after '@')
	Description string // Full description line (everything after '@')
	Sequence    []byte // Nucleotide sequence
	Quality     []byte // Quality scores (ASCII-encoded)
}

// Reader provides sequential access to FASTQ records.
type Reader struct {
	scanner  *bufio.Scanner
	encoding QualityEncoding
	err      error
}

// NewReader creates a new FASTQ reader from an io.Reader.
// encoding specifies the quality score encoding (Phred33 or Phred64).
func NewReader(r io.Reader, encoding QualityEncoding) *Reader {
	scanner := bufio.NewScanner(r)
	// Set larger buffer for long sequences
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024) // 10MB max token size
	return &Reader{
		scanner:  scanner,
		encoding: encoding,
	}
}

// Read reads the next FASTQ record.
// Returns io.EOF when no more records are available.
func (r *Reader) Read() (*Record, error) {
	if r.err != nil {
		return nil, r.err
	}

	// Read line 1: header (starts with '@')
	if !r.scanner.Scan() {
		if err := r.scanner.Err(); err != nil {
			r.err = err
			return nil, err
		}
		r.err = io.EOF
		return nil, io.EOF
	}
	header := r.scanner.Text()
	if !strings.HasPrefix(header, "@") {
		return nil, fmt.Errorf("expected '@' at start of FASTQ header, got: %s", header)
	}
	description := strings.TrimPrefix(header, "@")

	// Extract ID (first word)
	fields := strings.Fields(description)
	var id string
	if len(fields) > 0 {
		id = fields[0]
	}

	// Read line 2: sequence
	if !r.scanner.Scan() {
		if err := r.scanner.Err(); err != nil {
			return nil, err
		}
		return nil, fmt.Errorf("unexpected EOF reading sequence for record %s", id)
	}
	sequence := []byte(r.scanner.Text())

	// Read line 3: separator (starts with '+')
	if !r.scanner.Scan() {
		if err := r.scanner.Err(); err != nil {
			return nil, err
		}
		return nil, fmt.Errorf("unexpected EOF reading separator for record %s", id)
	}
	separator := r.scanner.Text()
	if !strings.HasPrefix(separator, "+") {
		return nil, fmt.Errorf("expected '+' separator, got: %s", separator)
	}

	// Read line 4: quality scores
	if !r.scanner.Scan() {
		if err := r.scanner.Err(); err != nil {
			return nil, err
		}
		return nil, fmt.Errorf("unexpected EOF reading quality for record %s", id)
	}
	quality := []byte(r.scanner.Text())

	// Validate that sequence and quality have same length
	if len(sequence) != len(quality) {
		return nil, fmt.Errorf("sequence and quality length mismatch for record %s: %d vs %d",
			id, len(sequence), len(quality))
	}

	return &Record{
		ID:          id,
		Description: description,
		Sequence:    sequence,
		Quality:     quality,
	}, nil
}

// ReadAll reads all FASTQ records from the reader.
func (r *Reader) ReadAll() ([]*Record, error) {
	var records []*Record
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return records, err
		}
		records = append(records, record)
	}
	return records, nil
}

// Writer provides sequential writing of FASTQ records.
type Writer struct {
	w        *bufio.Writer
	encoding QualityEncoding
}

// NewWriter creates a new FASTQ writer.
func NewWriter(w io.Writer, encoding QualityEncoding) *Writer {
	return &Writer{
		w:        bufio.NewWriter(w),
		encoding: encoding,
	}
}

// Write writes a FASTQ record.
func (fw *Writer) Write(record *Record) error {
	// Validate sequence and quality length match
	if len(record.Sequence) != len(record.Quality) {
		return fmt.Errorf("sequence and quality length mismatch: %d vs %d",
			len(record.Sequence), len(record.Quality))
	}

	// Write header
	if _, err := fmt.Fprintf(fw.w, "@%s\n", record.Description); err != nil {
		return err
	}

	// Write sequence
	if _, err := fw.w.Write(record.Sequence); err != nil {
		return err
	}
	if _, err := fw.w.WriteString("\n"); err != nil {
		return err
	}

	// Write separator
	if _, err := fw.w.WriteString("+\n"); err != nil {
		return err
	}

	// Write quality
	if _, err := fw.w.Write(record.Quality); err != nil {
		return err
	}
	if _, err := fw.w.WriteString("\n"); err != nil {
		return err
	}

	return nil
}

// Flush writes any buffered data to the underlying writer.
func (fw *Writer) Flush() error {
	return fw.w.Flush()
}

// WriteAll writes all records and flushes.
func (fw *Writer) WriteAll(records []*Record) error {
	for _, record := range records {
		if err := fw.Write(record); err != nil {
			return err
		}
	}
	return fw.Flush()
}

// Validate checks if a FASTQ record is valid.
func (r *Record) Validate() error {
	if r.ID == "" {
		return fmt.Errorf("record has empty ID")
	}
	if len(r.Sequence) == 0 {
		return fmt.Errorf("record %s has empty sequence", r.ID)
	}
	if len(r.Quality) == 0 {
		return fmt.Errorf("record %s has empty quality", r.ID)
	}
	if len(r.Sequence) != len(r.Quality) {
		return fmt.Errorf("record %s: sequence and quality length mismatch: %d vs %d",
			r.ID, len(r.Sequence), len(r.Quality))
	}
	return nil
}

// Length returns the sequence length.
func (r *Record) Length() int {
	return len(r.Sequence)
}

// QualityScores converts ASCII-encoded quality scores to numeric values.
// encoding specifies whether to use Phred33 or Phred64.
func (r *Record) QualityScores(encoding QualityEncoding) []int {
	scores := make([]int, len(r.Quality))
	offset := 33
	if encoding == Phred64 {
		offset = 64
	}
	for i, q := range r.Quality {
		scores[i] = int(q) - offset
	}
	return scores
}

// AverageQuality calculates the average quality score.
func (r *Record) AverageQuality(encoding QualityEncoding) float64 {
	if len(r.Quality) == 0 {
		return 0
	}
	scores := r.QualityScores(encoding)
	sum := 0
	for _, s := range scores {
		sum += s
	}
	return float64(sum) / float64(len(scores))
}

// MinQuality returns the minimum quality score in the record.
func (r *Record) MinQuality(encoding QualityEncoding) int {
	if len(r.Quality) == 0 {
		return 0
	}
	scores := r.QualityScores(encoding)
	min := scores[0]
	for _, s := range scores {
		if s < min {
			min = s
		}
	}
	return min
}

// ReverseComplement returns the reverse complement of the sequence with reversed quality.
func (r *Record) ReverseComplement() *Record {
	complement := map[byte]byte{
		'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
		'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
		'U': 'A', 'u': 'a',
		'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
		'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
		'D': 'H', 'H': 'D', 'N': 'N',
		'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
		'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
		'd': 'h', 'h': 'd', 'n': 'n',
		'-': '-',
	}

	revSeq := make([]byte, len(r.Sequence))
	revQual := make([]byte, len(r.Quality))

	for i, b := range r.Sequence {
		if comp, ok := complement[b]; ok {
			revSeq[len(r.Sequence)-1-i] = comp
		} else {
			revSeq[len(r.Sequence)-1-i] = b
		}
		revQual[len(r.Quality)-1-i] = r.Quality[i]
	}

	// Note: Description is preserved verbatim. We previously appended
	// " (reverse complement)" here, but that diverged from upstream seqtk
	// (which emits the original header unchanged) and broke downstream
	// parsers that key on the FASTQ description field.
	return &Record{
		ID:          r.ID,
		Description: r.Description,
		Sequence:    revSeq,
		Quality:     revQual,
	}
}

// Trim trims the sequence and quality from both ends based on quality threshold.
// threshold is the minimum quality score to keep (in Phred scale).
// Returns a new trimmed record.
func (r *Record) Trim(threshold int, encoding QualityEncoding) *Record {
	scores := r.QualityScores(encoding)

	// Find the first position >= threshold from the start
	start := 0
	for start < len(scores) && scores[start] < threshold {
		start++
	}

	// Find the last position >= threshold from the end
	end := len(scores) - 1
	for end >= 0 && scores[end] < threshold {
		end--
	}

	// If all bases are below threshold or invalid range
	if start > end {
		return &Record{
			ID:          r.ID,
			Description: r.Description,
			Sequence:    []byte{},
			Quality:     []byte{},
		}
	}

	return &Record{
		ID:          r.ID,
		Description: r.Description,
		Sequence:    r.Sequence[start : end+1],
		Quality:     r.Quality[start : end+1],
	}
}

// GCContent calculates the GC content percentage of the sequence.
func (r *Record) GCContent() float64 {
	if len(r.Sequence) == 0 {
		return 0
	}

	gcCount := 0
	for _, b := range r.Sequence {
		if b == 'G' || b == 'C' || b == 'g' || b == 'c' {
			gcCount++
		}
	}

	return float64(gcCount) / float64(len(r.Sequence)) * 100
}

// ToFasta converts a FASTQ record to FASTA format (drops quality scores).
func (r *Record) ToFasta() *bytes.Buffer {
	var buf bytes.Buffer
	fmt.Fprintf(&buf, ">%s\n", r.Description)
	buf.Write(r.Sequence)
	buf.WriteString("\n")
	return &buf
}
