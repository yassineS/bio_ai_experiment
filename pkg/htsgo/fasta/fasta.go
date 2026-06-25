// Package fasta provides utilities for reading and writing FASTA format files.
// FASTA is a text-based format for representing nucleotide or peptide sequences.
//
// Format specification:
//   - Header line begins with '>'
//   - Sequence data follows on subsequent lines
//   - Multiple sequences can be in one file
package fasta

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strings"
)

// Record represents a single FASTA sequence record.
type Record struct {
	ID          string // Sequence identifier (first word after '>')
	Description string // Full description line (everything after '>')
	Sequence    []byte // Nucleotide or amino acid sequence
}

// Reader provides sequential access to FASTA records.
//
// The Sequence of a Record returned by Read aliases an internal buffer that is
// REUSED by the next Read call, so it must be consumed (or copied) before
// calling Read again. This is what keeps the reader's peak memory at one
// record's sequence rather than reallocating a fresh multi-hundred-MB buffer
// per contig (a whole-chromosome FASTA otherwise cost ~5x upstream's RSS).
// ReadAll, which retains every record, copies each Sequence so its results stay
// independent; callers that retain a streamed Record's Sequence across a Read
// must copy it themselves.
type Reader struct {
	scanner       *bufio.Scanner
	err           error
	nextHeader    string // Buffer for the next header when we read ahead
	hasNextHeader bool
	seqBuffer     bytes.Buffer // reused sequence accumulator (see the type doc)
}

// NewReader creates a new FASTA reader from an io.Reader.
func NewReader(r io.Reader) *Reader {
	scanner := bufio.NewScanner(r)
	// Set larger buffer for long sequences
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024) // 10MB max token size
	return &Reader{scanner: scanner}
}

// Read reads the next FASTA record.
// Returns io.EOF when no more records are available.
func (r *Reader) Read() (*Record, error) {
	if r.err != nil {
		return nil, r.err
	}

	// Get header - either from buffer or scan for it
	var description string
	if r.hasNextHeader {
		description = r.nextHeader
		r.hasNextHeader = false
	} else {
		// Skip empty lines until we find a header
		for r.scanner.Scan() {
			line := r.scanner.Text()
			if strings.HasPrefix(line, ">") {
				description = strings.TrimPrefix(line, ">")
				break
			}
			// Skip empty or whitespace-only lines
			if strings.TrimSpace(line) == "" {
				continue
			}
		}

		if err := r.scanner.Err(); err != nil {
			r.err = err
			return nil, err
		}

		if description == "" {
			r.err = io.EOF
			return nil, io.EOF
		}
	}

	// Extract ID (first word) and full description
	fields := strings.Fields(description)
	var id string
	if len(fields) > 0 {
		id = fields[0]
	}

	// Read sequence lines until next header or EOF. Use scanner.Bytes() (not
	// Text()) and trim in place: a whole-chromosome record is millions of
	// lines, and allocating a string per line — as Text()+TrimSpace did —
	// churned hundreds of MB of short-lived garbage per contig, which (not the
	// sequence itself) dominated this reader's peak heap. scanner.Bytes()
	// returns the line backed by the scanner's buffer, so the header copy below
	// is taken as a string to detach it before the next Scan overwrites it.
	// Accumulate into the reader's reused buffer (Reset keeps its capacity), so
	// the next record overwrites this one rather than allocating afresh — the
	// returned Sequence therefore aliases the buffer and is valid only until the
	// next Read (see the Reader type doc).
	r.seqBuffer.Reset()
	for r.scanner.Scan() {
		line := r.scanner.Bytes()
		if len(line) > 0 && line[0] == '>' {
			// Save the header for next read.
			r.nextHeader = string(line[1:])
			r.hasNextHeader = true
			break
		}
		// Trim whitespace and append sequence (TrimSpace returns a subslice,
		// Write copies it into the buffer — no per-line allocation).
		r.seqBuffer.Write(bytes.TrimSpace(line))
	}

	if err := r.scanner.Err(); err != nil {
		r.err = err
		return nil, err
	}

	return &Record{
		ID:          id,
		Description: description,
		Sequence:    r.seqBuffer.Bytes(),
	}, nil
}

// ReadAll reads all FASTA records from the reader. Because every record is
// retained, each record's Sequence is copied out of the reader's reused buffer
// so the returned records are fully independent (Read alone aliases that
// buffer — see the Reader type doc).
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
		record.Sequence = append([]byte(nil), record.Sequence...)
		records = append(records, record)
	}
	return records, nil
}

// Writer provides sequential writing of FASTA records.
type Writer struct {
	w         *bufio.Writer
	lineWidth int // Maximum line width for sequence (0 = no wrapping)
}

// NewWriter creates a new FASTA writer.
// lineWidth specifies the maximum sequence line width (0 for no wrapping, default 80).
func NewWriter(w io.Writer, lineWidth int) *Writer {
	if lineWidth <= 0 {
		lineWidth = 80
	}
	return &Writer{
		w:         bufio.NewWriter(w),
		lineWidth: lineWidth,
	}
}

// Write writes a FASTA record.
func (fw *Writer) Write(record *Record) error {
	// Write header
	if _, err := fmt.Fprintf(fw.w, ">%s\n", record.Description); err != nil {
		return err
	}

	// Write sequence with line wrapping
	seq := record.Sequence
	for len(seq) > 0 {
		end := fw.lineWidth
		if end > len(seq) {
			end = len(seq)
		}
		if _, err := fw.w.Write(seq[:end]); err != nil {
			return err
		}
		if _, err := fw.w.WriteString("\n"); err != nil {
			return err
		}
		seq = seq[end:]
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

// Validate checks if a FASTA record is valid.
func (r *Record) Validate(allowAmbiguous bool) error {
	if r.ID == "" {
		return fmt.Errorf("record has empty ID")
	}
	if len(r.Sequence) == 0 {
		return fmt.Errorf("record %s has empty sequence", r.ID)
	}

	// Check for valid nucleotide or amino acid characters
	validNucleotides := "ACGTURYSWKMBDHVNacgturyswkmbdhvn-"
	validAminoAcids := "ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx*-"

	isNucleotide := true
	isAminoAcid := true

	for _, b := range r.Sequence {
		if !strings.ContainsRune(validNucleotides, rune(b)) {
			isNucleotide = false
		}
		if !strings.ContainsRune(validAminoAcids, rune(b)) {
			isAminoAcid = false
		}
	}

	if !isNucleotide && !isAminoAcid {
		if !allowAmbiguous {
			return fmt.Errorf("record %s contains invalid sequence characters", r.ID)
		}
	}

	return nil
}

// Length returns the sequence length.
func (r *Record) Length() int {
	return len(r.Sequence)
}

// ReverseComplement returns the reverse complement of a DNA sequence.
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

	revComp := make([]byte, len(r.Sequence))
	for i, b := range r.Sequence {
		if comp, ok := complement[b]; ok {
			revComp[len(r.Sequence)-1-i] = comp
		} else {
			revComp[len(r.Sequence)-1-i] = b
		}
	}

	// Note: Description is preserved verbatim. We previously appended
	// " (reverse complement)" here, but that diverged from upstream seqtk
	// (which emits the original header unchanged) and broke downstream
	// parsers that key on the FASTA description field.
	return &Record{
		ID:          r.ID,
		Description: r.Description,
		Sequence:    revComp,
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
