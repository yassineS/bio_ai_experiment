// Package seqtk provides core functionality for sequence processing.
// This is a Go reimplementation of seqtk, a fast FASTA/Q processor.
package seqtk

import (
	"bufio"
	"compress/bzip2"
	"compress/gzip"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
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

// MergePE interleaves two FASTA/FASTQ streams (in1, in2) into a single stream
// written to w, producing read1[0], read2[0], read1[1], read2[1], ... The two
// inputs must have the same record count and the same format (auto-detected
// from the first non-whitespace byte of each: '>' => FASTA, '@' => FASTQ); if
// the two streams disagree on format, an error is returned. If one stream runs
// short of records before the other, an error identifying the shorter stream
// (by 1-based name "in1"/"in2") and the pair index where the mismatch was
// detected is returned. The output preserves the input format (FASTA in =>
// FASTA out, FASTQ in => FASTQ out).
func MergePE(in1, in2 io.Reader, w io.Writer) error {
	br1, isFastq1 := peekIsFastq(in1)
	br2, isFastq2 := peekIsFastq(in2)
	if isFastq1 != isFastq2 {
		return fmt.Errorf("mergepe: input formats differ (in1 is %s, in2 is %s)",
			fmtName(isFastq1), fmtName(isFastq2))
	}

	if isFastq1 {
		return mergePEFastq(br1, br2, w)
	}
	return mergePEFasta(br1, br2, w)
}

// fmtName returns "FASTQ" or "FASTA" for use in error messages.
func fmtName(isFastq bool) string {
	if isFastq {
		return "FASTQ"
	}
	return "FASTA"
}

func mergePEFastq(in1, in2 io.Reader, w io.Writer) error {
	r1 := fastq.NewReader(in1, fastq.Phred33)
	r2 := fastq.NewReader(in2, fastq.Phred33)
	wr := fastq.NewWriter(w, fastq.Phred33)

	pair := 0
	for {
		rec1, err1 := r1.Read()
		rec2, err2 := r2.Read()

		if err1 == io.EOF && err2 == io.EOF {
			break
		}
		if err1 == io.EOF && err2 == nil {
			return fmt.Errorf("mergepe: in1 is shorter than in2 (mismatch at pair %d, in1 has %d records)", pair+1, pair)
		}
		if err2 == io.EOF && err1 == nil {
			return fmt.Errorf("mergepe: in2 is shorter than in1 (mismatch at pair %d, in2 has %d records)", pair+1, pair)
		}
		if err1 != nil && err1 != io.EOF {
			return fmt.Errorf("mergepe: error reading in1 at pair %d: %w", pair+1, err1)
		}
		if err2 != nil && err2 != io.EOF {
			return fmt.Errorf("mergepe: error reading in2 at pair %d: %w", pair+1, err2)
		}

		if err := wr.Write(rec1); err != nil {
			return err
		}
		if err := wr.Write(rec2); err != nil {
			return err
		}
		pair++
	}
	return wr.Flush()
}

func mergePEFasta(in1, in2 io.Reader, w io.Writer) error {
	r1 := fasta.NewReader(in1)
	r2 := fasta.NewReader(in2)
	wr := fasta.NewWriter(w, 0)

	pair := 0
	for {
		rec1, err1 := r1.Read()
		rec2, err2 := r2.Read()

		if err1 == io.EOF && err2 == io.EOF {
			break
		}
		if err1 == io.EOF && err2 == nil {
			return fmt.Errorf("mergepe: in1 is shorter than in2 (mismatch at pair %d, in1 has %d records)", pair+1, pair)
		}
		if err2 == io.EOF && err1 == nil {
			return fmt.Errorf("mergepe: in2 is shorter than in1 (mismatch at pair %d, in2 has %d records)", pair+1, pair)
		}
		if err1 != nil && err1 != io.EOF {
			return fmt.Errorf("mergepe: error reading in1 at pair %d: %w", pair+1, err1)
		}
		if err2 != nil && err2 != io.EOF {
			return fmt.Errorf("mergepe: error reading in2 at pair %d: %w", pair+1, err2)
		}

		if err := wr.Write(rec1); err != nil {
			return err
		}
		if err := wr.Write(rec2); err != nil {
			return err
		}
		pair++
	}
	return wr.Flush()
}

// CutNOptions holds parameters for CutN.
type CutNOptions struct {
	// MinN is the minimum length of a run of Ns required to trigger a cut.
	// Must be >= 1.
	MinN int
	// EmitGaps, if true, causes the cut-out N-runs to be written to GapsW
	// as BED-like records: "<chrom>\t<start0>\t<end>\tN\n" (0-based half-open).
	EmitGaps bool
	// GapsW is the writer for BED-format gap records when EmitGaps is true.
	// If nil and EmitGaps is true, gap records are silently dropped.
	GapsW io.Writer
}

// CutN reads a FASTA or FASTQ stream from in (auto-detected via the first
// non-whitespace byte: '>' => FASTA, '@' => FASTQ) and writes FASTA fragments
// to w, splitting each input sequence at runs of 'N' or 'n' of length >=
// opts.MinN. Each retained fragment is emitted as a new FASTA record named
// "<orig-name>:<start>-<end>", where coordinates are 1-based inclusive (start
// is the position of the first retained base, end the position of the last).
// If a record has no qualifying N-run it is emitted unchanged with its
// original header and no coordinate suffix. If a sequence is entirely Ns or
// the only fragments would be empty (e.g. leading/trailing N-runs only), no
// record is emitted for that input.
//
// When opts.EmitGaps is true and opts.GapsW is non-nil, the cut-out N-runs
// are written to opts.GapsW as BED-like records: "<chrom>\t<start0>\t<end>\tN"
// where coordinates are 0-based half-open.
//
// Returns an error if opts.MinN < 1, or on I/O errors.
func CutN(in io.Reader, w io.Writer, opts CutNOptions) error {
	if opts.MinN < 1 {
		return fmt.Errorf("cutN: -n/--min-n must be >= 1 (got %d)", opts.MinN)
	}

	br, isFastq := peekIsFastq(in)
	bw := bufio.NewWriter(w)

	emit := func(name string, seq []byte) error {
		runs := findNRuns(seq, opts.MinN)
		if opts.EmitGaps && opts.GapsW != nil {
			for _, r := range runs {
				if _, err := fmt.Fprintf(opts.GapsW, "%s\t%d\t%d\tN\n", name, r[0], r[1]); err != nil {
					return err
				}
			}
		}
		if len(runs) == 0 {
			// No qualifying N-run: still emit the record with the
			// upstream "name:1-len" suffix so downstream tooling sees a
			// uniform header layout. Skip empty sequences entirely.
			if len(seq) == 0 {
				return nil
			}
			header := fmt.Sprintf("%s:1-%d", name, len(seq))
			return writeFastaRecord(bw, header, seq, 0)
		}
		// Build fragment intervals between runs and emit non-empty ones.
		prev := 0
		for _, r := range runs {
			if r[0] > prev {
				header := fmt.Sprintf("%s:%d-%d", name, prev+1, r[0])
				if err := writeFastaRecord(bw, header, seq[prev:r[0]], 0); err != nil {
					return err
				}
			}
			prev = r[1]
		}
		if prev < len(seq) {
			header := fmt.Sprintf("%s:%d-%d", name, prev+1, len(seq))
			if err := writeFastaRecord(bw, header, seq[prev:], 0); err != nil {
				return err
			}
		}
		return nil
	}

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
			if err := emit(rec.ID, rec.Sequence); err != nil {
				return err
			}
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
			if err := emit(rec.ID, rec.Sequence); err != nil {
				return err
			}
		}
	}

	return bw.Flush()
}

// findNRuns returns the [start, end) (0-based, half-open) positions of every
// maximal run of 'N' or 'n' in seq whose length is >= minN. Runs are returned
// in left-to-right order.
func findNRuns(seq []byte, minN int) [][2]int {
	var runs [][2]int
	i := 0
	for i < len(seq) {
		if seq[i] == 'N' || seq[i] == 'n' {
			j := i + 1
			for j < len(seq) && (seq[j] == 'N' || seq[j] == 'n') {
				j++
			}
			if j-i >= minN {
				runs = append(runs, [2]int{i, j})
			}
			i = j
		} else {
			i++
		}
	}
	return runs
}

// SampleSeed is the default RNG seed used by `seqtk sample` when no
// `-s SEED` is given. Matches upstream's `if (kr == 0) kr = kr_srand(11);`
// at reference_code/seqtk/seqtk.c:1254.
const SampleSeed = 11

// Sample randomly samples a fraction of sequences. Uses the default
// upstream seed (SampleSeed = 11) so that the streaming `r < fraction`
// decision is reproducible and byte-for-byte identical to
// `seqtk sample`.
//
// Mirrors the streaming (no -2) path of `stk_sample` in
// reference_code/seqtk/seqtk.c: one drand call per record, keep iff the
// draw is strictly below `fraction`.
func Sample(input io.Reader, output io.Writer, fraction float64, isFastq bool, encoding fastq.QualityEncoding) error {
	return SampleSeeded(input, output, fraction, isFastq, encoding, SampleSeed)
}

// SampleSeeded is the explicit-seed form of Sample, equivalent to
// `seqtk sample -s SEED <in> FRACTION`. The seed feeds the in-tree
// MT19937-64 port (`krand`), which produces the same `kr_drand`
// sequence as upstream so output is byte-for-byte identical.
func SampleSeeded(input io.Reader, output io.Writer, fraction float64, isFastq bool, encoding fastq.QualityEncoding, seed uint64) error {
	if fraction <= 0 || fraction > 1 {
		return fmt.Errorf("fraction must be between 0 and 1")
	}
	kr := newKrand(seed)
	if isFastq {
		return sampleFastq(input, output, fraction, encoding, kr)
	}
	return sampleFasta(input, output, fraction, kr)
}

func sampleFasta(input io.Reader, output io.Writer, fraction float64, kr *krand) error {
	reader := fasta.NewReader(input)
	bw := bufio.NewWriter(output)

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		// Upstream draws one drand per record (the increment that bumps
		// n_seqs is also unconditional) so the per-record decision is
		// deterministic given seed + record index.
		r := kr.drand()
		if r < fraction {
			// Upstream's stk_printseq with UINT_MAX line_len emits the
			// sequence on a single un-wrapped line — preserve that here
			// rather than fasta.NewWriter's default 80-char wrap.
			if _, err := fmt.Fprintf(bw, ">%s\n", record.Description); err != nil {
				return err
			}
			if _, err := bw.Write(record.Sequence); err != nil {
				return err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return err
			}
		}
	}

	return bw.Flush()
}

func sampleFastq(input io.Reader, output io.Writer, fraction float64, encoding fastq.QualityEncoding, kr *krand) error {
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

		r := kr.drand()
		if r < fraction {
			if err := writer.Write(record); err != nil {
				return err
			}
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

// TrimfqMott implements upstream `seqtk trimfq`'s modified Mott
// algorithm (reference_code/seqtk/seqtk.c:361-441). Each base's quality
// `q` (clamped to [36, 127]) contributes `param - 10^(-(q-33)/10)` to a
// running sum `s`; if `s` exceeds the running maximum the trim window
// is extended; if `s` drops below zero the window is reset. The kept
// window is `[beg, end)`. A `minLen` floor (default 30) triggers a
// fallback window-based scan picking the highest-mean-quality `minLen`
// window. Phred+33 is assumed (upstream hardcodes the q_int2real table
// at offset 33).
//
// param is the error-rate threshold (upstream default 0.05). minLen is
// the minimum window length floor (upstream default 30). Both must be
// > 0 for non-trivial behaviour.
func TrimfqMott(input io.Reader, output io.Writer, param float64, minLen int) error {
	reader := fastq.NewReader(input, fastq.Phred33)
	bw := bufio.NewWriter(output)
	q2real := mottPhredErrTable()
	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		beg, end := mottTrimWindow(rec.Quality, param, minLen, q2real)
		// Match upstream's record emitter (seqtk.c:428-436): use '@'
		// for FASTQ records (rec.Quality non-empty) and '>' otherwise;
		// emit name + optional comment + trimmed seq + optional
		// "+\n<qual>".
		hasQual := len(rec.Quality) > 0
		if hasQual {
			if err := bw.WriteByte('@'); err != nil {
				return err
			}
		} else {
			if err := bw.WriteByte('>'); err != nil {
				return err
			}
		}
		if _, err := bw.WriteString(rec.Description); err != nil {
			return err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
		if _, err := bw.Write(rec.Sequence[beg:end]); err != nil {
			return err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
		if hasQual {
			if _, err := bw.WriteString("+\n"); err != nil {
				return err
			}
			if _, err := bw.Write(rec.Quality[beg:end]); err != nil {
				return err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return err
			}
		}
	}
	return bw.Flush()
}

// mottPhredErrTable precomputes `pow(10, -(q-33)/10)` for q in 0..127,
// matching upstream's `for (i = 0; i < 128; ++i) q_int2real[i] = pow(10., -(i - 33) / 10.);`
// at seqtk.c:395-396.
func mottPhredErrTable() [128]float64 {
	var t [128]float64
	for i := 0; i < 128; i++ {
		t[i] = math.Pow(10, -float64(i-33)/10.0)
	}
	return t
}

// mottTrimWindow runs upstream's per-record Mott loop over `qual`
// (Phred+33 ASCII) and returns the [beg, end) trim window. minLen is
// the floor below which a sliding-window fallback runs (seqtk.c:417-426).
func mottTrimWindow(qual []byte, param float64, minLen int, q2real [128]float64) (beg, end int) {
	n := len(qual)
	if n <= minLen {
		return 0, n
	}
	// Main Mott pass.
	var tmp, sBeg, sEnd int
	end = n
	s, max := 0.0, 0.0
	for i := 0; i < n; i++ {
		q := int(qual[i])
		if q < 36 {
			q = 36
		}
		if q > 127 {
			q = 127
		}
		s += param - q2real[q]
		if s > max {
			max = s
			sBeg = tmp
			sEnd = i + 1
		}
		if s < 0 {
			s = 0
			tmp = i + 1
		}
	}
	beg, end = sBeg, sEnd
	if max == 0 {
		// Upstream comment: "max never set; all low qual, just give
		// first min_len bp" (seqtk.c:414-415).
		beg, end = 0, minLen
	}
	if end-beg < minLen {
		// Fallback: pick the minLen window with the highest summed
		// raw quality. Mirrors seqtk.c:417-426.
		is := 0
		for i := 0; i < minLen; i++ {
			is += int(qual[i]) - 33
		}
		imax := is
		beg = 0
		for i := minLen; i < n; i++ {
			is += int(qual[i]) - int(qual[i-minLen])
			if imax < is {
				imax = is
				beg = i - minLen + 1
			}
		}
		end = beg + minLen
	}
	return beg, end
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
