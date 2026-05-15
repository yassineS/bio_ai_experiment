package prinseq

import (
	"bufio"
	"fmt"
	"io"
	"math"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// Stats holds sequence statistics
type Stats struct {
	NumReads            int            `json:"num_reads"`
	TotalBases          int            `json:"total_bases"`
	MinLength           int            `json:"min_length"`
	MaxLength           int            `json:"max_length"`
	AvgLength           float64        `json:"avg_length"`
	GCContent           float64        `json:"gc_content"`
	AvgQuality          float64        `json:"avg_quality,omitempty"` // Only for FASTQ
	NumNs               int            `json:"num_ns"`
	LengthDistribution  map[int]int    `json:"length_distribution,omitempty"`
	QualityDistribution map[int]int    `json:"quality_distribution,omitempty"`
	BaseComposition     map[string]int `json:"base_composition,omitempty"`
	Dinucleotides       map[string]int `json:"dinucleotides,omitempty"`
	PositionalQuality   []float64      `json:"positional_quality,omitempty"`
}

// CalculateStats computes statistics for FASTA or FASTQ files
func CalculateStats(reader io.Reader, isFastq bool) (*Stats, error) {
	return CalculateStatsWithEncoding(reader, isFastq, "sanger")
}

// CalculateStatsWithEncoding computes statistics with a specific quality encoding
func CalculateStatsWithEncoding(reader io.Reader, isFastq bool, qualType string) (*Stats, error) {
	stats := &Stats{
		MinLength: -1,
	}

	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 1024*1024), 1024*1024) // 1MB buffer

	if isFastq {
		offset := 33
		if qualType == "illumina" {
			offset = 64
		}
		return calculateFastqStatsWithOffset(scanner, stats, offset)
	}
	return calculateFastaStats(scanner, stats)
}

func calculateFastaStats(scanner *bufio.Scanner, stats *Stats) (*Stats, error) {
	var currentSeq string
	gcCount := 0

	for scanner.Scan() {
		line := scanner.Text()
		if len(line) == 0 {
			continue
		}

		if line[0] == '>' {
			// Process previous sequence
			if currentSeq != "" {
				processSequence(currentSeq, &gcCount, stats)
				currentSeq = ""
			}
		} else {
			currentSeq += line
		}
	}

	// Process last sequence
	if currentSeq != "" {
		processSequence(currentSeq, &gcCount, stats)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	// Calculate averages
	if stats.NumReads > 0 {
		stats.AvgLength = float64(stats.TotalBases) / float64(stats.NumReads)
		stats.GCContent = float64(gcCount) / float64(stats.TotalBases) * 100.0
	}

	return stats, nil
}

func calculateFastqStats(scanner *bufio.Scanner, stats *Stats) (*Stats, error) {
	return calculateFastqStatsWithOffset(scanner, stats, 33)
}

func calculateFastqStatsWithOffset(scanner *bufio.Scanner, stats *Stats, offset int) (*Stats, error) {
	gcCount := 0
	totalQuality := 0.0
	lineNum := 0

	var seq string
	var qual string

	for scanner.Scan() {
		line := scanner.Text()
		mod := lineNum % 4

		switch mod {
		case 0: // Header line
			if len(line) == 0 || line[0] != '@' {
				return nil, fmt.Errorf("invalid FASTQ format at line %d", lineNum+1)
			}
		case 1: // Sequence line
			seq = line
		case 2: // Plus line
			if len(line) == 0 || line[0] != '+' {
				return nil, fmt.Errorf("invalid FASTQ format at line %d", lineNum+1)
			}
		case 3: // Quality line
			qual = line
			if len(qual) != len(seq) {
				return nil, fmt.Errorf("quality length (%d) doesn't match sequence length (%d) at line %d",
					len(qual), len(seq), lineNum+1)
			}
			// Process the complete FASTQ record
			processSequence(seq, &gcCount, stats)
			totalQuality += calculateAvgQualityScoreWithOffset(qual, offset)
			seq = ""
			qual = ""
		}

		lineNum++
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	// Check if we have incomplete records
	if lineNum%4 != 0 {
		return nil, fmt.Errorf("incomplete FASTQ record")
	}

	// Calculate averages
	if stats.NumReads > 0 {
		stats.AvgLength = float64(stats.TotalBases) / float64(stats.NumReads)
		stats.GCContent = float64(gcCount) / float64(stats.TotalBases) * 100.0
		stats.AvgQuality = totalQuality / float64(stats.NumReads)
	}

	return stats, nil
}

func processSequence(seq string, gcCount *int, stats *Stats) {
	seqLen := len(seq)
	stats.NumReads++
	stats.TotalBases += seqLen

	if stats.MinLength == -1 || seqLen < stats.MinLength {
		stats.MinLength = seqLen
	}
	if seqLen > stats.MaxLength {
		stats.MaxLength = seqLen
	}

	// Count GC and Ns
	for _, base := range seq {
		switch base {
		case 'G', 'C', 'g', 'c':
			*gcCount++
		case 'N', 'n':
			stats.NumNs++
		}
	}
}

func calculateAvgQualityScore(qual string) float64 {
	return calculateAvgQualityScoreWithOffset(qual, 33)
}

func calculateAvgQualityScoreWithOffset(qual string, offset int) float64 {
	if len(qual) == 0 {
		return 0.0
	}

	total := 0
	for _, q := range qual {
		// Use specified offset (33 for Phred+33, 64 for Phred+64)
		total += int(q) - offset
	}
	return float64(total) / float64(len(qual))
}

// trimSequence applies various trimming operations to a sequence
func trimSequence(seq, qual string, opts FilterOptions) (string, string) {
	seqBytes := []byte(seq)
	var qualBytes []byte
	if qual != "" {
		qualBytes = []byte(qual)
	}

	// Apply trim_left (fixed position from left)
	if opts.TrimLeft > 0 && len(seqBytes) > opts.TrimLeft {
		seqBytes = seqBytes[opts.TrimLeft:]
		if len(qualBytes) > 0 {
			qualBytes = qualBytes[opts.TrimLeft:]
		}
	}

	// Apply trim_right (fixed position from right)
	if opts.TrimRight > 0 && len(seqBytes) > opts.TrimRight {
		seqBytes = seqBytes[:len(seqBytes)-opts.TrimRight]
		if len(qualBytes) > 0 {
			qualBytes = qualBytes[:len(qualBytes)-opts.TrimRight]
		}
	}

	// Apply trim_left_p (percentage from left)
	if opts.TrimLeftP > 0 && len(seqBytes) > 0 {
		trimPos := (len(seqBytes) * opts.TrimLeftP) / 100
		if trimPos > 0 && trimPos < len(seqBytes) {
			seqBytes = seqBytes[trimPos:]
			if len(qualBytes) > 0 {
				qualBytes = qualBytes[trimPos:]
			}
		}
	}

	// Apply trim_right_p (percentage from right)
	if opts.TrimRightP > 0 && len(seqBytes) > 0 {
		trimPos := (len(seqBytes) * opts.TrimRightP) / 100
		if trimPos > 0 && trimPos < len(seqBytes) {
			seqBytes = seqBytes[:len(seqBytes)-trimPos]
			if len(qualBytes) > 0 {
				qualBytes = qualBytes[:len(qualBytes)-trimPos]
			}
		}
	}

	// Apply trim_ns_left (trim poly-N from left)
	if opts.TrimNsLeft > 0 {
		seqBytes, qualBytes = trimPolyNLeft(seqBytes, qualBytes, opts.TrimNsLeft)
	}

	// Apply trim_ns_right (trim poly-N from right)
	if opts.TrimNsRight > 0 {
		seqBytes, qualBytes = trimPolyNRight(seqBytes, qualBytes, opts.TrimNsRight)
	}

	// Apply trim_tail_left (trim poly-A/T from left)
	if opts.TrimTailLeft > 0 {
		seqBytes, qualBytes = trimPolyATLeft(seqBytes, qualBytes, opts.TrimTailLeft)
	}

	// Apply trim_tail_right (trim poly-A/T from right)
	if opts.TrimTailRight > 0 {
		seqBytes, qualBytes = trimPolyATRight(seqBytes, qualBytes, opts.TrimTailRight)
	}

	// Apply quality-based trimming from left
	if opts.TrimQualL > 0 && len(qualBytes) > 0 {
		seqBytes, qualBytes = trimQualityLeft(seqBytes, qualBytes, opts.TrimQualL, phredOffset(opts.QualType))
	}

	// Apply quality-based trimming from right
	if opts.TrimQualR > 0 && len(qualBytes) > 0 {
		seqBytes, qualBytes = trimQualityRight(seqBytes, qualBytes, opts.TrimQualR, phredOffset(opts.QualType))
	}

	return string(seqBytes), string(qualBytes)
}

func trimPolyNLeft(seq, qual []byte, minLen int) ([]byte, []byte) {
	nCount := 0
	for i := 0; i < len(seq); i++ {
		if seq[i] == 'N' || seq[i] == 'n' {
			nCount++
		} else {
			break
		}
	}

	if nCount >= minLen {
		seq = seq[nCount:]
		if len(qual) > 0 {
			qual = qual[nCount:]
		}
	}
	return seq, qual
}

func trimPolyNRight(seq, qual []byte, minLen int) ([]byte, []byte) {
	nCount := 0
	for i := len(seq) - 1; i >= 0; i-- {
		if seq[i] == 'N' || seq[i] == 'n' {
			nCount++
		} else {
			break
		}
	}

	if nCount >= minLen {
		seq = seq[:len(seq)-nCount]
		if len(qual) > 0 {
			qual = qual[:len(qual)-nCount]
		}
	}
	return seq, qual
}

// trimPolyATLeft trims a poly-A or poly-T tail from the 5'-end of seq when
// the run is at least minLen bases long. The logic mirrors upstream prinseq:
// (1) match either `^A{minLen}` or `^T{minLen}`, picking whichever exists;
// (2) extend the run with more of the SAME base or with Ns; (3) if the run
// covers the entire read, trim everything (the caller will length-filter it
// out separately, matching upstream's `good = 0` branch).
//
// Before this fix we treated A and T as interchangeable, so a run that
// alternated A and T (or that started with A and continued with T) would
// over-trim by one or more bases relative to upstream.
func trimPolyATLeft(seq, qual []byte, minLen int) ([]byte, []byte) {
	if len(seq) < minLen {
		return seq, qual
	}
	// Detect a homogeneous A-run or T-run of length minLen at the head.
	var anchor byte
	switch {
	case allEqualCase(seq[:minLen], 'A'):
		anchor = 'A'
	case allEqualCase(seq[:minLen], 'T'):
		anchor = 'T'
	default:
		return seq, qual
	}
	// Extend with the same base or N.
	end := minLen
	for end < len(seq) && (matchesCase(seq[end], anchor) || seq[end] == 'N' || seq[end] == 'n') {
		end++
	}
	seq = seq[end:]
	if len(qual) > 0 {
		qual = qual[end:]
	}
	return seq, qual
}

// trimPolyATRight is the 3'-end counterpart of trimPolyATLeft.
func trimPolyATRight(seq, qual []byte, minLen int) ([]byte, []byte) {
	if len(seq) < minLen {
		return seq, qual
	}
	tail := seq[len(seq)-minLen:]
	var anchor byte
	switch {
	case allEqualCase(tail, 'A'):
		anchor = 'A'
	case allEqualCase(tail, 'T'):
		anchor = 'T'
	default:
		return seq, qual
	}
	start := len(seq) - minLen
	for start > 0 && (matchesCase(seq[start-1], anchor) || seq[start-1] == 'N' || seq[start-1] == 'n') {
		start--
	}
	seq = seq[:start]
	if len(qual) > 0 {
		qual = qual[:start]
	}
	return seq, qual
}

// matchesCase reports whether b equals base or its lowercase form.
func matchesCase(b, base byte) bool {
	return b == base || b == (base|0x20)
}

// allEqualCase reports whether every byte in s equals base or its
// lowercase form.
func allEqualCase(s []byte, base byte) bool {
	for _, b := range s {
		if !matchesCase(b, base) {
			return false
		}
	}
	return true
}

// phredOffset returns the ASCII offset used to decode quality characters for
// the given prinseq quality-type string. It returns 64 for Illumina 1.3-1.7
// style encodings ("illumina") and 33 for Sanger/Phred+33 (the default).
func phredOffset(qualType string) int {
	if qualType == "illumina" {
		return 64
	}
	return 33
}

func trimQualityLeft(seq, qual []byte, threshold, offset int) ([]byte, []byte) {
	trimPos := 0
	for i := 0; i < len(qual); i++ {
		if int(qual[i])-offset < threshold {
			trimPos = i + 1
		} else {
			break
		}
	}

	if trimPos > 0 && trimPos < len(seq) {
		seq = seq[trimPos:]
		qual = qual[trimPos:]
	}
	return seq, qual
}

func trimQualityRight(seq, qual []byte, threshold, offset int) ([]byte, []byte) {
	trimPos := len(qual)
	for i := len(qual) - 1; i >= 0; i-- {
		if int(qual[i])-offset < threshold {
			trimPos = i
		} else {
			break
		}
	}

	if trimPos < len(seq) {
		seq = seq[:trimPos]
		qual = qual[:trimPos]
	}
	return seq, qual
}

// FilterOptions holds filtering parameters
type FilterOptions struct {
	MinLen        int
	MaxLen        int
	MinGC         float64
	MaxGC         float64
	MinQual       float64
	MaxNsP        float64 // Max percentage of Ns
	MaxNsN        int     // Max number of Ns
	TrimLeft      int
	TrimRight     int
	TrimLeftP     int // Trim percentage from left
	TrimRightP    int // Trim percentage from right
	TrimQualL     int // Quality threshold for left trimming
	TrimQualR     int // Quality threshold for right trimming
	TrimNsLeft    int // Trim poly-N from left
	TrimNsRight   int // Trim poly-N from right
	TrimTailLeft  int // Trim poly-A/T from left
	TrimTailRight int // Trim poly-A/T from right
	MinQualMean   float64
	MaxQualMean   float64
	Derep         int       // Duplicate removal mode (1=exact, 4=reverse complement)
	DerepMin      int       // Minimum occurrences to keep
	QualType      string    // Quality encoding type: "sanger" (Phred+33) or "illumina" (Phred+64)
	OutBad        io.Writer // Writer for rejected sequences (optional)
	LcMethod      string    // Low complexity method: "dust" or "entropy"
	LcThreshold   float64   // Low complexity threshold (7 for dust, 70 for entropy)
}

// Filter filters a FASTA/FASTQ file based on the given options
func Filter(reader io.Reader, writer io.Writer, isFastq bool, opts FilterOptions) error {
	if isFastq {
		return filterFastq(reader, writer, opts)
	}
	return filterFasta(reader, writer, opts)
}

// getQualityEncoding returns the appropriate quality encoding based on options
func getQualityEncoding(qualType string) fastq.QualityEncoding {
	if qualType == "illumina" {
		return fastq.Phred64
	}
	return fastq.Phred33 // Default to sanger
}

func filterFasta(reader io.Reader, writer io.Writer, opts FilterOptions) error {
	fastaReader := fasta.NewReader(reader)
	fastaWriter := fasta.NewWriter(writer, 80)

	var badWriter *fasta.Writer
	if opts.OutBad != nil {
		badWriter = fasta.NewWriter(opts.OutBad, 80)
	}

	seenSeqs := make(map[string]int) // For duplicate tracking

	for {
		record, err := fastaReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		seq := string(record.Sequence)

		// Apply trimming
		seq, _ = trimSequence(seq, "", opts)

		// Check for duplicates if derep is enabled
		if opts.Derep > 0 {
			if shouldFilterDuplicate(seq, seenSeqs, opts) {
				if badWriter != nil {
					record.Sequence = []byte(seq)
					if err := badWriter.Write(record); err != nil {
						return err
					}
				}
				continue
			}
		}

		// Apply filters
		if shouldFilterSequence(seq, "", opts) {
			// Write to bad output if enabled
			if badWriter != nil {
				record.Sequence = []byte(seq)
				if err := badWriter.Write(record); err != nil {
					return err
				}
			}
			continue
		}

		// Update record with trimmed sequence
		record.Sequence = []byte(seq)

		// Write filtered record
		if err := fastaWriter.Write(record); err != nil {
			return err
		}
	}

	if err := fastaWriter.Flush(); err != nil {
		return err
	}
	if badWriter != nil {
		return badWriter.Flush()
	}
	return nil
}

// writePrinseqFastq writes one FASTQ record in the upstream PRINSEQ-lite
// layout, which repeats the sequence header on the "+" separator line
// (e.g. "+read1\n" rather than the bare "+\n" emitted by sickle / fastp).
// Upstream switches to a bare "+" only when -no_qual_header is supplied.
// Centralising the write here lets us swap encodings later without
// touching the filter loop.
func writePrinseqFastq(w *bufio.Writer, rec *fastq.Record) error {
	if _, err := fmt.Fprintf(w, "@%s\n", rec.Description); err != nil {
		return err
	}
	if _, err := w.Write(rec.Sequence); err != nil {
		return err
	}
	if err := w.WriteByte('\n'); err != nil {
		return err
	}
	if _, err := fmt.Fprintf(w, "+%s\n", rec.Description); err != nil {
		return err
	}
	if _, err := w.Write(rec.Quality); err != nil {
		return err
	}
	return w.WriteByte('\n')
}

func filterFastq(reader io.Reader, writer io.Writer, opts FilterOptions) error {
	encoding := getQualityEncoding(opts.QualType)
	fastqReader := fastq.NewReader(reader, encoding)
	bw := bufio.NewWriter(writer)

	var bbw *bufio.Writer
	if opts.OutBad != nil {
		bbw = bufio.NewWriter(opts.OutBad)
	}

	seenSeqs := make(map[string]int) // For duplicate tracking

	for {
		record, err := fastqReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		seq := string(record.Sequence)
		qual := string(record.Quality)

		// Apply trimming
		seq, qual = trimSequence(seq, qual, opts)

		// Check for duplicates if derep is enabled
		if opts.Derep > 0 {
			if shouldFilterDuplicate(seq, seenSeqs, opts) {
				if bbw != nil {
					record.Sequence = []byte(seq)
					record.Quality = []byte(qual)
					if err := writePrinseqFastq(bbw, record); err != nil {
						return err
					}
				}
				continue
			}
		}

		// Apply filters
		if shouldFilterSequence(seq, qual, opts) {
			// Write to bad output if enabled
			if bbw != nil {
				record.Sequence = []byte(seq)
				record.Quality = []byte(qual)
				if err := writePrinseqFastq(bbw, record); err != nil {
					return err
				}
			}
			continue
		}

		// Update record with trimmed sequence and quality
		record.Sequence = []byte(seq)
		record.Quality = []byte(qual)

		// Write filtered record in upstream-compatible format.
		if err := writePrinseqFastq(bw, record); err != nil {
			return err
		}
	}

	if err := bw.Flush(); err != nil {
		return err
	}
	if bbw != nil {
		return bbw.Flush()
	}
	return nil
}

func shouldFilterSequence(seq, qual string, opts FilterOptions) bool {
	seqLen := len(seq)

	// Length filters
	if opts.MinLen > 0 && seqLen < opts.MinLen {
		return true
	}
	if opts.MaxLen > 0 && seqLen > opts.MaxLen {
		return true
	}

	// GC content filter
	gcCount := 0
	nCount := 0
	for _, base := range seq {
		switch base {
		case 'G', 'C', 'g', 'c':
			gcCount++
		case 'N', 'n':
			nCount++
		}
	}

	if seqLen > 0 {
		gcContent := float64(gcCount) / float64(seqLen) * 100.0
		if opts.MinGC > 0 && gcContent < opts.MinGC {
			return true
		}
		if opts.MaxGC > 0 && gcContent > opts.MaxGC {
			return true
		}
	}

	// N content filters
	if opts.MaxNsN > 0 && nCount > opts.MaxNsN {
		return true
	}
	if opts.MaxNsP > 0 && seqLen > 0 {
		nPercent := float64(nCount) / float64(seqLen) * 100.0
		if nPercent > opts.MaxNsP {
			return true
		}
	}

	// Quality filters (only for FASTQ). Honour the per-options Phred
	// offset so that `-phred64` inputs are decoded against ASCII 64.
	// Before this fix the mean was always computed under Phred+33,
	// silently mis-classifying Phred+64 reads as much higher quality
	// than they actually were.
	if qual != "" {
		avgQual := calculateAvgQualityScoreWithOffset(qual, phredOffset(opts.QualType))
		if opts.MinQualMean > 0 && avgQual < opts.MinQualMean {
			return true
		}
		if opts.MaxQualMean > 0 && avgQual > opts.MaxQualMean {
			return true
		}
	}

	// Complexity filter
	if opts.LcMethod != "" && opts.LcThreshold > 0 {
		if opts.LcMethod == "dust" {
			score := calculateDustScore(seq)
			if score > opts.LcThreshold {
				return true
			}
		} else if opts.LcMethod == "entropy" {
			score := calculateEntropy(seq)
			if score < opts.LcThreshold {
				return true
			}
		}
	}

	return false
}

// calculateDustScore computes DUST score for low-complexity filtering
// Lower scores indicate higher complexity; higher scores indicate low complexity
func calculateDustScore(seq string) float64 {
	if len(seq) < 3 {
		return 0.0
	}

	// Count triplet frequencies
	triplets := make(map[string]int)
	for i := 0; i <= len(seq)-3; i++ {
		triplet := seq[i : i+3]
		triplets[triplet]++
	}

	// Calculate DUST score as sum of (count * (count-1) / 2) for each triplet
	score := 0.0
	for _, count := range triplets {
		if count > 1 {
			score += float64(count * (count - 1) / 2)
		}
	}

	// Normalize by sequence length
	if len(seq) > 0 {
		score = score / float64(len(seq)-2) * 10.0
	}

	return score
}

// calculateEntropy computes Shannon entropy for complexity filtering
// Higher scores indicate higher complexity; lower scores indicate low complexity
func calculateEntropy(seq string) float64 {
	if len(seq) == 0 {
		return 0.0
	}

	// Count base frequencies
	counts := make(map[rune]int)
	for _, base := range seq {
		counts[base]++
	}

	// Calculate Shannon entropy
	entropy := 0.0
	length := float64(len(seq))
	for _, count := range counts {
		if count > 0 {
			p := float64(count) / length
			entropy -= p * math.Log2(p)
		}
	}

	// Normalize to percentage (0-100)
	// Maximum entropy for DNA is log2(4) = 2.0, so normalize to 100
	return (entropy / 2.0) * 100.0
}

func shouldFilterDuplicate(seq string, seenSeqs map[string]int, opts FilterOptions) bool {
	if opts.Derep == 0 {
		return false
	}

	// Mode 1: exact duplicate
	if opts.Derep&1 != 0 {
		seenSeqs[seq]++
		if seenSeqs[seq] >= opts.DerepMin {
			return true
		}
	}

	// Mode 4: reverse complement exact duplicate
	if opts.Derep&4 != 0 {
		revComp := reverseComplement(seq)
		seenSeqs[revComp]++
		if seenSeqs[revComp] >= opts.DerepMin {
			return true
		}
	}

	return false
}

func reverseComplement(seq string) string {
	complement := map[rune]rune{
		'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
		'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
		'N': 'N', 'n': 'n',
	}

	result := make([]rune, len(seq))
	for i, base := range seq {
		if comp, ok := complement[base]; ok {
			result[len(seq)-1-i] = comp
		} else {
			result[len(seq)-1-i] = base
		}
	}
	return string(result)
}

// FilterPaired filters paired-end FASTA/FASTQ files
func FilterPaired(reader1, reader2 io.Reader, writer1, writer2 io.Writer, isFastq bool, opts FilterOptions) error {
	if isFastq {
		return filterPairedFastq(reader1, reader2, writer1, writer2, opts)
	}
	return filterPairedFasta(reader1, reader2, writer1, writer2, opts)
}

func filterPairedFasta(reader1, reader2 io.Reader, writer1, writer2 io.Writer, opts FilterOptions) error {
	fastaReader1 := fasta.NewReader(reader1)
	fastaReader2 := fasta.NewReader(reader2)
	fastaWriter1 := fasta.NewWriter(writer1, 80)
	fastaWriter2 := fasta.NewWriter(writer2, 80)

	seenSeqs := make(map[string]int)

	for {
		record1, err1 := fastaReader1.Read()
		record2, err2 := fastaReader2.Read()

		if err1 == io.EOF && err2 == io.EOF {
			break
		}
		if err1 == io.EOF || err2 == io.EOF {
			return fmt.Errorf("paired files have different number of sequences")
		}
		if err1 != nil {
			return err1
		}
		if err2 != nil {
			return err2
		}

		seq1 := string(record1.Sequence)
		seq2 := string(record2.Sequence)

		// Apply trimming
		seq1, _ = trimSequence(seq1, "", opts)
		seq2, _ = trimSequence(seq2, "", opts)

		// Check duplicates (consider both reads together)
		if opts.Derep > 0 {
			combinedSeq := seq1 + "|" + seq2
			if shouldFilterDuplicate(combinedSeq, seenSeqs, opts) {
				continue
			}
		}

		// Filter: if either read fails, both are filtered
		if shouldFilterSequence(seq1, "", opts) || shouldFilterSequence(seq2, "", opts) {
			continue
		}

		// Update records
		record1.Sequence = []byte(seq1)
		record2.Sequence = []byte(seq2)

		// Write both records
		if err := fastaWriter1.Write(record1); err != nil {
			return err
		}
		if err := fastaWriter2.Write(record2); err != nil {
			return err
		}
	}

	if err := fastaWriter1.Flush(); err != nil {
		return err
	}
	return fastaWriter2.Flush()
}

func filterPairedFastq(reader1, reader2 io.Reader, writer1, writer2 io.Writer, opts FilterOptions) error {
	encoding := getQualityEncoding(opts.QualType)
	fastqReader1 := fastq.NewReader(reader1, encoding)
	fastqReader2 := fastq.NewReader(reader2, encoding)
	bw1 := bufio.NewWriter(writer1)
	bw2 := bufio.NewWriter(writer2)

	seenSeqs := make(map[string]int)

	for {
		record1, err1 := fastqReader1.Read()
		record2, err2 := fastqReader2.Read()

		if err1 == io.EOF && err2 == io.EOF {
			break
		}
		if err1 == io.EOF || err2 == io.EOF {
			return fmt.Errorf("paired files have different number of sequences")
		}
		if err1 != nil {
			return err1
		}
		if err2 != nil {
			return err2
		}

		seq1 := string(record1.Sequence)
		qual1 := string(record1.Quality)
		seq2 := string(record2.Sequence)
		qual2 := string(record2.Quality)

		// Apply trimming
		seq1, qual1 = trimSequence(seq1, qual1, opts)
		seq2, qual2 = trimSequence(seq2, qual2, opts)

		// Check duplicates (consider both reads together)
		if opts.Derep > 0 {
			combinedSeq := seq1 + "|" + seq2
			if shouldFilterDuplicate(combinedSeq, seenSeqs, opts) {
				continue
			}
		}

		// Filter: if either read fails, both are filtered
		if shouldFilterSequence(seq1, qual1, opts) || shouldFilterSequence(seq2, qual2, opts) {
			continue
		}

		// Update records
		record1.Sequence = []byte(seq1)
		record1.Quality = []byte(qual1)
		record2.Sequence = []byte(seq2)
		record2.Quality = []byte(qual2)

		// Write both records in upstream-compatible FASTQ format.
		if err := writePrinseqFastq(bw1, record1); err != nil {
			return err
		}
		if err := writePrinseqFastq(bw2, record2); err != nil {
			return err
		}
	}

	if err := bw1.Flush(); err != nil {
		return err
	}
	return bw2.Flush()
}

// CalculateEnhancedStats computes enhanced statistics including distributions and dinucleotides
func CalculateEnhancedStats(reader io.Reader, isFastq bool) (*Stats, error) {
	return CalculateEnhancedStatsWithEncoding(reader, isFastq, "sanger")
}

// CalculateEnhancedStatsWithEncoding computes enhanced statistics with specific quality encoding
func CalculateEnhancedStatsWithEncoding(reader io.Reader, isFastq bool, qualType string) (*Stats, error) {
	stats := &Stats{
		MinLength:           -1,
		LengthDistribution:  make(map[int]int),
		QualityDistribution: make(map[int]int),
		BaseComposition:     make(map[string]int),
		Dinucleotides:       make(map[string]int),
	}

	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 1024*1024), 1024*1024) // 1MB buffer

	if isFastq {
		offset := 33
		if qualType == "illumina" {
			offset = 64
		}
		return calculateEnhancedFastqStats(scanner, stats, offset)
	}
	return calculateEnhancedFastaStats(scanner, stats)
}

func calculateEnhancedFastaStats(scanner *bufio.Scanner, stats *Stats) (*Stats, error) {
	var currentSeq string
	gcCount := 0

	for scanner.Scan() {
		line := scanner.Text()
		if len(line) == 0 {
			continue
		}

		if line[0] == '>' {
			// Process previous sequence
			if currentSeq != "" {
				processEnhancedSequence(currentSeq, "", &gcCount, stats)
				currentSeq = ""
			}
		} else {
			currentSeq += line
		}
	}

	// Process last sequence
	if currentSeq != "" {
		processEnhancedSequence(currentSeq, "", &gcCount, stats)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	// Calculate averages
	if stats.NumReads > 0 {
		stats.AvgLength = float64(stats.TotalBases) / float64(stats.NumReads)
		stats.GCContent = float64(gcCount) / float64(stats.TotalBases) * 100.0
	}

	return stats, nil
}

func calculateEnhancedFastqStats(scanner *bufio.Scanner, stats *Stats, offset int) (*Stats, error) {
	gcCount := 0
	totalQuality := 0.0
	lineNum := 0
	maxLen := 0

	var seq string
	var qual string

	for scanner.Scan() {
		line := scanner.Text()
		mod := lineNum % 4

		switch mod {
		case 0: // Header line
			if len(line) == 0 || line[0] != '@' {
				return nil, fmt.Errorf("invalid FASTQ format at line %d", lineNum+1)
			}
		case 1: // Sequence line
			seq = line
		case 2: // Plus line
			if len(line) == 0 || line[0] != '+' {
				return nil, fmt.Errorf("invalid FASTQ format at line %d", lineNum+1)
			}
		case 3: // Quality line
			qual = line
			if len(qual) != len(seq) {
				return nil, fmt.Errorf("quality length (%d) doesn't match sequence length (%d) at line %d",
					len(qual), len(seq), lineNum+1)
			}
			// Process the complete FASTQ record
			processEnhancedSequence(seq, qual, &gcCount, stats)
			avgQual := calculateAvgQualityScoreWithOffset(qual, offset)
			totalQuality += avgQual

			// Update quality distribution
			qualInt := int(avgQual)
			stats.QualityDistribution[qualInt]++

			// Update positional quality
			if len(seq) > maxLen {
				maxLen = len(seq)
				// Extend positional quality array
				for len(stats.PositionalQuality) < maxLen {
					stats.PositionalQuality = append(stats.PositionalQuality, 0)
				}
			}
			for i, q := range qual {
				qualVal := float64(int(q) - offset)
				if i < len(stats.PositionalQuality) {
					stats.PositionalQuality[i] += qualVal
				}
			}

			seq = ""
			qual = ""
		}

		lineNum++
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	// Check if we have incomplete records
	if lineNum%4 != 0 {
		return nil, fmt.Errorf("incomplete FASTQ record")
	}

	// Calculate averages
	if stats.NumReads > 0 {
		stats.AvgLength = float64(stats.TotalBases) / float64(stats.NumReads)
		stats.GCContent = float64(gcCount) / float64(stats.TotalBases) * 100.0
		stats.AvgQuality = totalQuality / float64(stats.NumReads)

		// Average positional quality
		for i := range stats.PositionalQuality {
			stats.PositionalQuality[i] /= float64(stats.NumReads)
		}
	}

	return stats, nil
}

func processEnhancedSequence(seq, qual string, gcCount *int, stats *Stats) {
	seqLen := len(seq)
	stats.NumReads++
	stats.TotalBases += seqLen

	if stats.MinLength == -1 || seqLen < stats.MinLength {
		stats.MinLength = seqLen
	}
	if seqLen > stats.MaxLength {
		stats.MaxLength = seqLen
	}

	// Update length distribution
	stats.LengthDistribution[seqLen]++

	// Count bases, GC, Ns, and dinucleotides
	var prevBase byte
	for i, base := range seq {
		// Base composition
		baseStr := string(base)
		stats.BaseComposition[baseStr]++

		switch base {
		case 'G', 'C', 'g', 'c':
			*gcCount++
		case 'N', 'n':
			stats.NumNs++
		}

		// Dinucleotides
		if i > 0 {
			dinuc := string([]byte{prevBase, byte(base)})
			stats.Dinucleotides[dinuc]++
		}
		prevBase = byte(base)
	}
}
