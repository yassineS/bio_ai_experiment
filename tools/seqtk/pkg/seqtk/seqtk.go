// Package seqtk provides core functionality for sequence processing.
// This is a Go reimplementation of seqtk, a fast FASTA/Q processor.
package seqtk

import (
	"bufio"
	"compress/bzip2"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
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

// subseqRegion describes a half-open [Start, End) interval (0-based) requested
// for a particular sequence. A region with Start == 0 and End == -1 (a sentinel)
// means "the whole sequence" — this is how name-list entries are represented.
type subseqRegion struct {
	Start int
	End   int // -1 means "to the end of the sequence"
}

// subseqSpec holds the parsed contents of a subseq region/name file: an ordered
// list of sequence names and, for each name, the list of requested regions.
type subseqSpec struct {
	order   []string                  // sequence names in the order first seen
	regions map[string][]subseqRegion // name -> requested regions (nil/empty for name-list mode)
	isBED   bool                      // true if the file was detected as BED
	seen    map[string]bool           // membership set for quick lookup
}

func newSubseqSpec() *subseqSpec {
	return &subseqSpec{
		regions: make(map[string][]subseqRegion),
		seen:    make(map[string]bool),
	}
}

func (s *subseqSpec) add(name string, reg *subseqRegion) {
	if !s.seen[name] {
		s.seen[name] = true
		s.order = append(s.order, name)
	}
	if reg != nil {
		s.regions[name] = append(s.regions[name], *reg)
	}
}

// looksLikeBED reports whether a non-comment line, split on whitespace, has at
// least three fields and the second and third fields are integers — matching how
// upstream seqtk auto-detects a BED file versus a plain name list.
func looksLikeBED(fields []string) bool {
	if len(fields) < 3 {
		return false
	}
	if _, err := strconv.Atoi(fields[1]); err != nil {
		return false
	}
	if _, err := strconv.Atoi(fields[2]); err != nil {
		return false
	}
	return true
}

// parseSubseqSpec reads a name list or a BED file from r. It auto-detects the
// format from the first non-comment line: if that line looks like BED (>= 3
// whitespace/tab fields with integer second and third fields) the whole file is
// parsed as BED; otherwise it is treated as a name list (one name per line,
// anything after the first whitespace-delimited token is ignored).
func parseSubseqSpec(r io.Reader) (*subseqSpec, error) {
	spec := newSubseqSpec()
	scanner := bufio.NewScanner(r)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024)

	decided := false
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") ||
			strings.HasPrefix(line, "track") || strings.HasPrefix(line, "browser") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) == 0 {
			continue
		}
		if !decided {
			spec.isBED = looksLikeBED(fields)
			decided = true
		}
		if spec.isBED {
			if !looksLikeBED(fields) {
				// Tolerate the occasional malformed line in a BED file by skipping it.
				continue
			}
			start, _ := strconv.Atoi(fields[1])
			end, _ := strconv.Atoi(fields[2])
			spec.add(fields[0], &subseqRegion{Start: start, End: end})
		} else {
			spec.add(fields[0], nil)
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return spec, nil
}

// peekIsFastq peeks at the first non-whitespace byte of r (via a *bufio.Reader
// it returns) to decide whether the stream is FASTQ ('@') or FASTA ('>').
func peekIsFastq(r io.Reader) (*bufio.Reader, bool) {
	br := bufio.NewReaderSize(r, 64*1024)
	for {
		b, err := br.Peek(1)
		if err != nil {
			return br, false
		}
		if b[0] == ' ' || b[0] == '\t' || b[0] == '\n' || b[0] == '\r' {
			if _, err := br.Discard(1); err != nil {
				return br, false
			}
			continue
		}
		return br, b[0] == '@'
	}
}

// writeFastaRecord writes a single FASTA record (header + sequence) to w,
// wrapping sequence lines at lineLen characters (lineLen <= 0 means no wrapping).
func writeFastaRecord(w *bufio.Writer, header string, seq []byte, lineLen int) error {
	if _, err := fmt.Fprintf(w, ">%s\n", header); err != nil {
		return err
	}
	if lineLen <= 0 {
		if _, err := w.Write(seq); err != nil {
			return err
		}
		return w.WriteByte('\n')
	}
	for len(seq) > 0 {
		end := lineLen
		if end > len(seq) {
			end = len(seq)
		}
		if _, err := w.Write(seq[:end]); err != nil {
			return err
		}
		if err := w.WriteByte('\n'); err != nil {
			return err
		}
		seq = seq[end:]
	}
	return nil
}

// emitSubseqRecord writes the output for one input sequence (identified by name,
// full header and sequence bytes) according to the parsed spec. In name-list
// mode it emits the whole record with its original header; in BED mode it emits
// one FASTA record per region, named "name:start+1-end", clamping end to the
// sequence length and skipping regions whose start is at or past the end. It
// returns the number of records written.
func emitSubseqRecord(w *bufio.Writer, spec *subseqSpec, name, header string, seq []byte, lineLen int) (int, error) {
	regions := spec.regions[name]
	if len(regions) == 0 {
		// Name-list mode (or BED entry with no regions): emit the whole record.
		if err := writeFastaRecord(w, header, seq, lineLen); err != nil {
			return 0, err
		}
		return 1, nil
	}
	n := 0
	for _, reg := range regions {
		start := reg.Start
		end := reg.End
		if end < 0 || end > len(seq) {
			end = len(seq)
		}
		if start < 0 {
			start = 0
		}
		if start >= len(seq) || start >= end {
			fmt.Fprintf(os.Stderr, "[seqtk subseq] warning: region %s:%d-%d is out of range; skipped\n", name, reg.Start+1, reg.End)
			continue
		}
		h := fmt.Sprintf("%s:%d-%d", name, start+1, end)
		if err := writeFastaRecord(w, h, seq[start:end], lineLen); err != nil {
			return n, err
		}
		n++
	}
	return n, nil
}

// Subseq extracts subsequences from a FASTA/FASTQ stream (in) given a region
// specification (regionSpec), which is either a list of sequence names (one per
// line) or a BED file of regions; the format is auto-detected. Output is always
// FASTA, written to w with sequence lines wrapped at lineLen characters (0 = no
// wrapping). Names present in the spec but absent from the input, and BED
// regions that fall outside their sequence, produce a warning on stderr and are
// skipped rather than causing an error.
func Subseq(in io.Reader, regionSpec io.Reader, w io.Writer, lineLen int) error {
	spec, err := parseSubseqSpec(regionSpec)
	if err != nil {
		return err
	}

	br, isFastq := peekIsFastq(in)
	bw := bufio.NewWriter(w)

	emitted := make(map[string]bool)

	if isFastq {
		reader := fastq.NewReader(br, fastq.Phred33)
		for {
			rec, err := reader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if !spec.seen[rec.ID] {
				continue
			}
			if _, err := emitSubseqRecord(bw, spec, rec.ID, rec.Description, rec.Sequence, lineLen); err != nil {
				return err
			}
			emitted[rec.ID] = true
		}
	} else {
		reader := fasta.NewReader(br)
		for {
			rec, err := reader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if !spec.seen[rec.ID] {
				continue
			}
			if _, err := emitSubseqRecord(bw, spec, rec.ID, rec.Description, rec.Sequence, lineLen); err != nil {
				return err
			}
			emitted[rec.ID] = true
		}
	}

	for _, name := range spec.order {
		if !emitted[name] {
			fmt.Fprintf(os.Stderr, "[seqtk subseq] warning: sequence %q not found in input; skipped\n", name)
		}
	}

	return bw.Flush()
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
