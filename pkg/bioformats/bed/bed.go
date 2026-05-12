// Package bed provides utilities for reading and writing BED (Browser Extensible Data) format files.
// BED format is used to represent genomic regions and annotations.
//
// Format specification:
//   - Tab-delimited text file
//   - Minimum 3 required fields: chrom, chromStart, chromEnd
//   - Up to 12 standard fields defined
//   - 0-based, half-open coordinates [start, end)
package bed

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Record represents a BED format record.
// BED format supports 3-12 fields. Fields beyond Name are optional.
type Record struct {
	Chrom       string   // Chromosome name
	ChromStart  int      // Start position (0-based)
	ChromEnd    int      // End position (exclusive)
	Name        string   // Name of the feature (optional, BED4+)
	Score       int      // Score (0-1000, optional, BED5+)
	Strand      string   // Strand: '+', '-', or '.' (optional, BED6+)
	ThickStart  int      // Thick start position (optional, BED7+)
	ThickEnd    int      // Thick end position (optional, BED8+)
	ItemRGB     string   // RGB color value (optional, BED9+)
	BlockCount  int      // Number of blocks (optional, BED10+)
	BlockSizes  []int    // Block sizes (optional, BED11+)
	BlockStarts []int    // Block starts (optional, BED12+)
	ExtraFields []string // Additional custom fields beyond BED12
}

// Reader provides sequential access to BED records.
type Reader struct {
	scanner *bufio.Scanner
	err     error
}

// NewReader creates a new BED reader from an io.Reader.
func NewReader(r io.Reader) *Reader {
	return &Reader{
		scanner: bufio.NewScanner(r),
	}
}

// Read reads the next BED record.
// Returns io.EOF when no more records are available.
func (r *Reader) Read() (*Record, error) {
	if r.err != nil {
		return nil, r.err
	}

	for r.scanner.Scan() {
		line := r.scanner.Text()

		// Skip empty lines and comments
		line = strings.TrimSpace(line)
		if line == "" || strings.HasPrefix(line, "#") || strings.HasPrefix(line, "track") || strings.HasPrefix(line, "browser") {
			continue
		}

		fields := strings.Split(line, "\t")
		if len(fields) < 3 {
			return nil, fmt.Errorf("BED record must have at least 3 fields, got %d", len(fields))
		}

		record := &Record{
			Chrom: fields[0],
		}

		// Parse chromStart
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("invalid chromStart %s: %v", fields[1], err)
		}
		record.ChromStart = start

		// Parse chromEnd
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return nil, fmt.Errorf("invalid chromEnd %s: %v", fields[2], err)
		}
		record.ChromEnd = end

		// Parse optional fields
		if len(fields) > 3 {
			record.Name = fields[3]
		}

		if len(fields) > 4 {
			score, err := strconv.Atoi(fields[4])
			if err != nil {
				return nil, fmt.Errorf("invalid score %s: %v", fields[4], err)
			}
			record.Score = score
		}

		if len(fields) > 5 {
			record.Strand = fields[5]
		}

		if len(fields) > 6 {
			thickStart, err := strconv.Atoi(fields[6])
			if err != nil {
				return nil, fmt.Errorf("invalid thickStart %s: %v", fields[6], err)
			}
			record.ThickStart = thickStart
		}

		if len(fields) > 7 {
			thickEnd, err := strconv.Atoi(fields[7])
			if err != nil {
				return nil, fmt.Errorf("invalid thickEnd %s: %v", fields[7], err)
			}
			record.ThickEnd = thickEnd
		}

		if len(fields) > 8 {
			record.ItemRGB = fields[8]
		}

		if len(fields) > 9 {
			blockCount, err := strconv.Atoi(fields[9])
			if err != nil {
				return nil, fmt.Errorf("invalid blockCount %s: %v", fields[9], err)
			}
			record.BlockCount = blockCount
		}

		if len(fields) > 10 {
			sizes := strings.Split(strings.TrimSuffix(fields[10], ","), ",")
			record.BlockSizes = make([]int, len(sizes))
			for i, s := range sizes {
				size, err := strconv.Atoi(s)
				if err != nil {
					return nil, fmt.Errorf("invalid block size %s: %v", s, err)
				}
				record.BlockSizes[i] = size
			}
		}

		if len(fields) > 11 {
			starts := strings.Split(strings.TrimSuffix(fields[11], ","), ",")
			record.BlockStarts = make([]int, len(starts))
			for i, s := range starts {
				start, err := strconv.Atoi(s)
				if err != nil {
					return nil, fmt.Errorf("invalid block start %s: %v", s, err)
				}
				record.BlockStarts[i] = start
			}
		}

		// Store any extra fields
		if len(fields) > 12 {
			record.ExtraFields = fields[12:]
		}

		return record, nil
	}

	if err := r.scanner.Err(); err != nil {
		r.err = err
		return nil, err
	}

	r.err = io.EOF
	return nil, io.EOF
}

// ReadAll reads all BED records from the reader.
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

// Writer provides sequential writing of BED records.
type Writer struct {
	w *bufio.Writer
}

// NewWriter creates a new BED writer.
func NewWriter(w io.Writer) *Writer {
	return &Writer{
		w: bufio.NewWriter(w),
	}
}

// Write writes a BED record.
// The number of fields written depends on which fields are populated.
func (w *Writer) Write(record *Record) error {
	fields := []string{
		record.Chrom,
		strconv.Itoa(record.ChromStart),
		strconv.Itoa(record.ChromEnd),
	}

	// Add optional fields only if they exist
	if record.Name != "" {
		fields = append(fields, record.Name)

		if record.Score != 0 {
			fields = append(fields, strconv.Itoa(record.Score))

			if record.Strand != "" {
				fields = append(fields, record.Strand)

				if record.ThickStart != 0 || record.ThickEnd != 0 {
					fields = append(fields,
						strconv.Itoa(record.ThickStart),
						strconv.Itoa(record.ThickEnd))

					if record.ItemRGB != "" {
						fields = append(fields, record.ItemRGB)

						if record.BlockCount != 0 {
							fields = append(fields, strconv.Itoa(record.BlockCount))

							if len(record.BlockSizes) > 0 {
								sizes := make([]string, len(record.BlockSizes))
								for i, size := range record.BlockSizes {
									sizes[i] = strconv.Itoa(size)
								}
								fields = append(fields, strings.Join(sizes, ","))

								if len(record.BlockStarts) > 0 {
									starts := make([]string, len(record.BlockStarts))
									for i, start := range record.BlockStarts {
										starts[i] = strconv.Itoa(start)
									}
									fields = append(fields, strings.Join(starts, ","))
								}
							}
						}
					}
				}
			}
		}
	}

	// Add extra fields
	fields = append(fields, record.ExtraFields...)

	if _, err := fmt.Fprintln(w.w, strings.Join(fields, "\t")); err != nil {
		return err
	}

	return nil
}

// Flush writes any buffered data to the underlying writer.
func (w *Writer) Flush() error {
	return w.w.Flush()
}

// WriteAll writes all records and flushes.
func (w *Writer) WriteAll(records []*Record) error {
	for _, record := range records {
		if err := w.Write(record); err != nil {
			return err
		}
	}
	return w.Flush()
}

// Length returns the length of the genomic region.
func (r *Record) Length() int {
	return r.ChromEnd - r.ChromStart
}

// Overlaps checks if this record overlaps with another record.
// Records overlap if they share any genomic positions.
func (r *Record) Overlaps(other *Record) bool {
	if r.Chrom != other.Chrom {
		return false
	}
	return r.ChromStart < other.ChromEnd && other.ChromStart < r.ChromEnd
}

// Contains checks if this record completely contains another record.
func (r *Record) Contains(other *Record) bool {
	if r.Chrom != other.Chrom {
		return false
	}
	return r.ChromStart <= other.ChromStart && r.ChromEnd >= other.ChromEnd
}

// Validate performs basic validation on the record.
func (r *Record) Validate() error {
	if r.Chrom == "" {
		return fmt.Errorf("chromosome name is empty")
	}
	if r.ChromStart < 0 {
		return fmt.Errorf("chromStart must be non-negative")
	}
	if r.ChromEnd <= r.ChromStart {
		return fmt.Errorf("chromEnd must be greater than chromStart")
	}
	if r.Strand != "" && r.Strand != "+" && r.Strand != "-" && r.Strand != "." {
		return fmt.Errorf("invalid strand: %s (must be +, -, or .)", r.Strand)
	}
	if r.Score < 0 || r.Score > 1000 {
		return fmt.Errorf("score must be between 0 and 1000")
	}
	return nil
}
